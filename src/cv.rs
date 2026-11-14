//! Валідація: leave-one-year-out (LOYO), вкладена LOYO, метрики на агрегатах і в підгрупах,
//! підбір гіперпараметрів.
//!
//! Спостереження — **всі** точки S-NPP відкладеного року (зокрема там, де MODIS не бачив нічого),
//! прогноз — лише з MODIS. Тобто перевіряється саме те, що модель робить для 2000–2011 років.

use crate::data::{FxMap, LC_GROUP_NAMES, NLC};
use crate::model::{self, ConstModel, Ctx, Hyper, Predictor, Sensor, EVAL_LEVEL_1DEG, EVAL_LEVEL_2DEG};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

pub const MAPE_MIN: f64 = 20.0;
/// z для двостороннього 90 % інтервалу.
pub const Z90: f64 = 1.644_853_6;

/// Квантиль 0,95 розподілу Стьюдента (двосторонній 90 % інтервал) для `df` ступенів свободи.
pub fn t95(df: usize) -> f64 {
    const T: [f64; 30] = [
        6.314, 2.920, 2.353, 2.132, 2.015, 1.943, 1.895, 1.860, 1.833, 1.812, 1.796, 1.782, 1.771, 1.761, 1.753, 1.746, 1.740, 1.734, 1.729,
        1.725, 1.721, 1.717, 1.714, 1.711, 1.708, 1.706, 1.703, 1.701, 1.699, 1.697,
    ];
    match df {
        0 => f64::INFINITY,
        1..=30 => T[df - 1],
        _ => Z90 + 2.4 / df as f64,
    }
}

/// Множник для прогнозного 90 % інтервалу з σ, оціненої за `n` залишками:
/// квантиль Стьюдента × √(1 + 1/n) (враховує невизначеність самої σ на малій вибірці).
pub fn pred_mult(n: usize) -> f64 {
    if n < 2 {
        return f64::INFINITY;
    }
    t95(n - 1) * (1.0 + 1.0 / n as f64).sqrt()
}

pub fn par_map<R: Send>(n: usize, threads: usize, f: impl Fn(usize) -> R + Sync) -> Vec<R> {
    let next = AtomicUsize::new(0);
    let results: Mutex<Vec<Option<R>>> = Mutex::new((0..n).map(|_| None).collect());
    let th = threads.max(1).min(n.max(1));
    std::thread::scope(|s| {
        for _ in 0..th {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::SeqCst);
                if i >= n {
                    break;
                }
                let r = f(i);
                results.lock().unwrap()[i] = Some(r);
            });
        }
    });
    results.into_inner().unwrap().into_iter().map(|x| x.expect("усі задачі виконані")).collect()
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Metrics {
    pub n: usize,
    pub sum_obs: f64,
    pub sum_pred: f64,
    /// Σpred / Σobs − 1
    pub bias: f64,
    /// Σ|pred − obs| / Σobs
    pub wape: f64,
    /// Σ|pred − obs| / n — у точках VIIRS на одиницю агрегації
    pub mae: f64,
    /// середнє |pred − obs| / obs для одиниць з obs ≥ MAPE_MIN
    pub mape: f64,
    /// R² на ln(1 + x)
    pub r2_log: f64,
}

pub fn metrics(pairs: &[(f64, f64)]) -> Metrics {
    let n = pairs.len();
    if n == 0 {
        return Metrics { bias: f64::NAN, wape: f64::NAN, mae: f64::NAN, mape: f64::NAN, r2_log: f64::NAN, ..Default::default() };
    }
    let so: f64 = pairs.iter().map(|p| p.0).sum();
    let sp: f64 = pairs.iter().map(|p| p.1).sum();
    let sae: f64 = pairs.iter().map(|p| (p.1 - p.0).abs()).sum();
    let ape: Vec<f64> = pairs.iter().filter(|p| p.0 >= MAPE_MIN).map(|p| (p.1 - p.0).abs() / p.0).collect();
    let lo: Vec<f64> = pairs.iter().map(|p| p.0.ln_1p()).collect();
    let mean = lo.iter().sum::<f64>() / n as f64;
    let ss_tot: f64 = lo.iter().map(|x| (x - mean).powi(2)).sum();
    let ss_res: f64 = pairs.iter().map(|p| (p.1.max(0.0).ln_1p() - p.0.ln_1p()).powi(2)).sum();
    Metrics {
        n,
        sum_obs: so,
        sum_pred: sp,
        bias: if so > 0.0 { sp / so - 1.0 } else { f64::NAN },
        wape: if so > 0.0 { sae / so } else { f64::NAN },
        mae: sae / n as f64,
        mape: if ape.is_empty() { f64::NAN } else { ape.iter().sum::<f64>() / ape.len() as f64 },
        r2_log: if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { f64::NAN },
    }
}

#[derive(Clone, Copy, Default)]
pub struct Cell {
    pub o: f64,
    pub p: f64,
    pub fo: f64,
    pub fp: f64,
}

/// Агрегати одного відкладеного року.
#[derive(Default, Clone)]
pub struct YearAgg {
    pub year: i32,
    pub ua: [(f64, f64); 12],
    /// ключ: (регіон 2°) << 16 | місяць << 8 | день/ніч << 4 | група покриву
    pub fine: FxMap<u64, Cell>,
    pub r1: FxMap<u32, (f64, f64)>,
    pub tot: (f64, f64),
    pub frp_tot: (f64, f64),
}

pub fn collect(ctx: &Ctx, m: &dyn Predictor, year: i32) -> YearAgg {
    let mut ya = YearAgg { year, ..Default::default() };
    let Some(yi) = ctx.yi(year) else { return ya };
    let ds = ctx.ds;
    for &i in &ctx.ov_by_year[yi] {
        let i = i as usize;
        let a = &ds.agg[i];
        let (p, pf) = m.predict(ctx, i);
        let (o, of) = (a.v as f64, a.fv);
        let mo = ctx.month[i] as u64;
        ya.ua[mo as usize].0 += o;
        ya.ua[mo as usize].1 += p;
        let r2 = (ctx.lv.idx[EVAL_LEVEL_2DEG][i] / 24) as u64;
        let key = (r2 << 16) | (mo << 8) | ((ds.dn[i] as u64) << 4) | ds.lc_of(ds.cell[i]) as u64;
        let e = ya.fine.entry(key).or_default();
        e.o += o;
        e.p += p;
        e.fo += of;
        e.fp += pf;
        let s1 = ctx.lv.idx[EVAL_LEVEL_1DEG][i];
        let e = ya.r1.entry((s1 / 24) * 12 + s1 % 12).or_default();
        e.0 += o;
        e.1 += p;
        ya.tot.0 += o;
        ya.tot.1 += p;
        ya.frp_tot.0 += of;
        ya.frp_tot.1 += pf;
    }
    ya
}

#[derive(Clone, Default)]
pub struct Eval {
    pub ua: Metrics,
    pub r2deg: Metrics,
    pub r1deg: Metrics,
    pub frp_r2deg: Metrics,
    /// середнє |Σpred/Σobs − 1| річних сум по Україні
    pub annual_err: f64,
    /// середня відстань між сезонними профілями (частки місяців у році), 0…1
    pub profile_err: f64,
    /// покриття 90 % інтервалу для річних сум; σ — з інших років (leave-one-out)
    pub cov_annual: (usize, usize),
    /// те саме для місячних сум по Україні (obs ≥ 20)
    pub cov_month: (usize, usize),
    /// (тип підгрупи, назва, метрики на рівні 2° × місяць)
    pub subgroups: Vec<(String, String, Metrics)>,
    pub years: Vec<YearAgg>,
}

const SEASONS: [(&str, [usize; 3]); 4] = [("зима (гру–лют)", [11, 0, 1]), ("весна (бер–тра)", [2, 3, 4]), ("літо (чер–сер)", [5, 6, 7]), ("осінь (вер–лис)", [8, 9, 10])];

pub fn summarize(years: Vec<YearAgg>) -> Eval {
    let keep = |p: &(f64, f64)| p.0 > 0.0 || p.1 > 0.0;
    let ua: Vec<(f64, f64)> = years.iter().flat_map(|y| y.ua.iter().copied()).filter(keep).collect();
    let r1: Vec<(f64, f64)> = years.iter().flat_map(|y| y.r1.values().copied()).filter(keep).collect();

    // одиниці 2° × місяць (сума по день/ніч і покриву) + розрізи
    let mut unit: FxMap<(i32, u64), Cell> = FxMap::default();
    let mut unit_dn: FxMap<(i32, u64, u64), (f64, f64)> = FxMap::default();
    let mut unit_lc: FxMap<(i32, u64, u64), (f64, f64)> = FxMap::default();
    for y in &years {
        for (&k, c) in &y.fine {
            let u = k >> 8 << 8; // регіон + місяць
            let e = unit.entry((y.year, u)).or_default();
            e.o += c.o;
            e.p += c.p;
            e.fo += c.fo;
            e.fp += c.fp;
            let e = unit_dn.entry((y.year, u, (k >> 4) & 15)).or_default();
            e.0 += c.o;
            e.1 += c.p;
            let e = unit_lc.entry((y.year, u, k & 15)).or_default();
            e.0 += c.o;
            e.1 += c.p;
        }
    }
    let r2: Vec<(f64, f64)> = unit.values().map(|c| (c.o, c.p)).filter(keep).collect();
    let fr: Vec<(f64, f64)> = unit.values().map(|c| (c.fo, c.fp)).filter(keep).collect();

    let mut subgroups = Vec::new();
    for (name, months) in SEASONS {
        let v: Vec<(f64, f64)> =
            unit.iter().filter(|((_, u), _)| months.contains(&(((u >> 8) & 255) as usize))).map(|(_, c)| (c.o, c.p)).filter(keep).collect();
        subgroups.push(("сезон".to_string(), name.to_string(), metrics(&v)));
    }
    for (dn, name) in [(0u64, "день"), (1, "ніч")] {
        let v: Vec<(f64, f64)> = unit_dn.iter().filter(|((_, _, d), _)| *d == dn).map(|(_, p)| *p).filter(keep).collect();
        subgroups.push(("день/ніч".to_string(), name.to_string(), metrics(&v)));
    }
    for g in 0..NLC as u64 {
        let v: Vec<(f64, f64)> = unit_lc.iter().filter(|((_, _, l), _)| *l == g).map(|(_, p)| *p).filter(keep).collect();
        if !v.is_empty() {
            subgroups.push(("земний покрив".to_string(), LC_GROUP_NAMES[g as usize].to_string(), metrics(&v)));
        }
    }
    // Групувати за фактом некоректно для оцінки зміщення: одиниці з випадково малим фактом
    // завжди виглядають «завищеними» (регресія до середнього). Правильно — за прогнозом.
    for (name, lo, hi) in [("низька (<20 точок)", 0.0, 20.0), ("середня (20–200)", 20.0, 200.0), ("висока (>200)", 200.0, f64::INFINITY)] {
        let v: Vec<(f64, f64)> = r2.iter().copied().filter(|p| p.1 >= lo && p.1 < hi).collect();
        subgroups.push(("активність за прогнозом".to_string(), name.to_string(), metrics(&v)));
    }
    for (name, lo, hi) in [("низька (<20 точок)", 0.0, 20.0), ("середня (20–200)", 20.0, 200.0), ("висока (>200)", 200.0, f64::INFINITY)] {
        let v: Vec<(f64, f64)> = r2.iter().copied().filter(|p| p.0 >= lo && p.0 < hi).collect();
        subgroups.push(("активність за фактом".to_string(), name.to_string(), metrics(&v)));
    }

    // річні суми, профіль, покриття
    let ann: Vec<(f64, f64)> = years.iter().map(|y| y.tot).filter(|t| t.0 > 0.0 && t.1 > 0.0).collect();
    let annual_err = if ann.is_empty() { f64::NAN } else { ann.iter().map(|t| (t.1 / t.0 - 1.0).abs()).sum::<f64>() / ann.len() as f64 };
    let mut prof = Vec::new();
    for y in &years {
        let (so, sp) = y.ua.iter().fold((0.0, 0.0), |a, x| (a.0 + x.0, a.1 + x.1));
        if so > 0.0 && sp > 0.0 {
            prof.push(0.5 * y.ua.iter().map(|x| (x.0 / so - x.1 / sp).abs()).sum::<f64>());
        }
    }
    let profile_err = if prof.is_empty() { f64::NAN } else { prof.iter().sum::<f64>() / prof.len() as f64 };
    let loo_cov = |res: &[(i32, f64)]| -> (usize, usize) {
        let mut inside = 0;
        for &(yr, r) in res {
            let other: Vec<f64> = res.iter().filter(|x| x.0 != yr).map(|x| x.1).collect();
            if other.is_empty() {
                continue;
            }
            let s = (other.iter().map(|x| x * x).sum::<f64>() / other.len() as f64).sqrt();
            if r.abs() <= pred_mult(other.len()) * s {
                inside += 1;
            }
        }
        (inside, res.len())
    };
    let res_a: Vec<(i32, f64)> = years.iter().filter(|y| y.tot.0 > 0.0 && y.tot.1 > 0.0).map(|y| (y.year, (y.tot.1 / y.tot.0).ln())).collect();
    let res_m: Vec<(i32, f64)> = years
        .iter()
        .flat_map(|y| y.ua.iter().filter(|x| x.0 >= MAPE_MIN && x.1 > 0.0).map(move |x| (y.year, (x.1 / x.0).ln())))
        .collect();

    Eval {
        ua: metrics(&ua),
        r2deg: metrics(&r2),
        r1deg: metrics(&r1),
        frp_r2deg: metrics(&fr),
        annual_err,
        profile_err,
        cov_annual: loo_cov(&res_a),
        cov_month: loo_cov(&res_m),
        subgroups,
        years,
    }
}

pub fn loyo(ctx: &Ctx, sensor: Sensor, h: Hyper, fit_years: &[i32], threads: usize) -> Eval {
    let ys = par_map(fit_years.len(), threads, |k| {
        let y = fit_years[k];
        let yw = ctx.year_weights(fit_years, Some(y));
        let f = model::fit(ctx, sensor, h, &yw);
        collect(ctx, &f, y)
    });
    summarize(ys)
}

pub fn loyo_const(ctx: &Ctx, sensor: Sensor, fit_years: &[i32], threads: usize) -> Eval {
    let ys = par_map(fit_years.len(), threads, |k| {
        let y = fit_years[k];
        let yw = ctx.year_weights(fit_years, Some(y));
        let f = ConstModel::fit(ctx, sensor, &yw);
        collect(ctx, &f, y)
    });
    summarize(ys)
}

pub fn eval_years(ctx: &Ctx, m: &dyn Predictor, years: &[i32], threads: usize) -> Eval {
    let ys = par_map(years.len(), threads, |k| collect(ctx, m, years[k]));
    summarize(ys)
}

pub fn score(e: &Eval) -> f64 {
    let s = e.r2deg.wape;
    if s.is_finite() {
        s
    } else {
        f64::INFINITY
    }
}

pub const LAMBDAS: [f64; 10] = [3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0, 30000.0, 100000.0];
pub const SMOOTHS: [f64; 3] = [0.0, 0.25, 0.5];
pub const DEPTHS: [usize; 5] = [1, 2, 3, 4, 5];

pub struct TuneEntry {
    pub sensor: Sensor,
    pub hyper: Hyper,
    pub ua: Metrics,
    pub r2deg: Metrics,
    pub r1deg: Metrics,
}

/// Підбір гіперпараметрів покоординатним спуском (до 3 проходів) від [`Hyper::default`];
/// критерій — LOYO WAPE на рівні 2°×місяць. Потім окремо — спосіб FRP за LOYO WAPE FRP.
pub fn tune(ctx: &Ctx, sensor: Sensor, fit_years: &[i32], threads: usize, verbose: bool, log: &mut Vec<TuneEntry>) -> (Hyper, Eval) {
    let mut cache: Vec<(Hyper, f64)> = Vec::new();
    let mut eval = |h: Hyper, log: &mut Vec<TuneEntry>| -> f64 {
        if let Some((_, s)) = cache.iter().find(|(c, _)| *c == h) {
            return *s;
        }
        let e = loyo(ctx, sensor, h, fit_years, threads);
        let s = score(&e);
        if verbose {
            eprintln!(
                "    {:<10} λ={:<6} sm={:<4} розмір={:<5} low={:<5} глиб={} патерн={:<5} інтенс={:<5} lc={:<5} пропущ={} EM={} → WAPE 2°×міс={:.4} bias={:+.4}",
                sensor.name(),
                h.lambda,
                h.smooth,
                h.bins,
                h.modis_low,
                h.depth,
                h.pattern,
                h.intensity,
                h.lc,
                h.miss_alloc,
                h.miss_model,
                s,
                e.r2deg.bias
            );
        }
        log.push(TuneEntry { sensor, hyper: h, ua: e.ua, r2deg: e.r2deg, r1deg: e.r1deg });
        cache.push((h, s));
        s
    };
    let mut best = Hyper::default();
    let mut best_s = eval(best, log);
    let pattern_ok = sensor == Sensor::TerraAqua;
    let lc_ok = ctx.lv_lc.is_some();
    for _sweep in 0..3 {
        let mut changed = false;
        // розподіл пропущеної маси і її модель — найвпливовіші параметри, тому перевіряються першими
        for &param in &[8usize, 9, 0, 1, 2, 3, 4, 5, 6, 7] {
            let b = best;
            let cands: Vec<Hyper> = match param {
                0 => LAMBDAS.iter().map(|&x| Hyper { lambda: x, ..b }).collect(),
                1 => SMOOTHS.iter().map(|&x| Hyper { smooth: x, ..b }).collect(),
                2 => [false, true].iter().map(|&x| Hyper { bins: x, ..b }).collect(),
                3 => [true, false].iter().map(|&x| Hyper { modis_low: x, ..b }).collect(),
                4 => DEPTHS.iter().map(|&x| Hyper { depth: x, ..b }).collect(),
                5 if pattern_ok => [false, true].iter().map(|&x| Hyper { pattern: x, ..b }).collect(),
                6 => [false, true].iter().map(|&x| Hyper { intensity: x, ..b }).collect(),
                7 if lc_ok => [false, true].iter().map(|&x| Hyper { lc: x, ..b }).collect(),
                8 => (0..4u8).map(|x| Hyper { miss_alloc: x, ..b }).collect(),
                9 => (0..3u8).map(|x| Hyper { miss_model: x, ..b }).collect(),
                _ => Vec::new(),
            };
            for c in cands {
                let s = eval(c, log);
                if s < best_s - 1e-9 {
                    best = c;
                    best_s = s;
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let mut best_e = loyo(ctx, sensor, best, fit_years, threads);
    for mode in 0..model::FRP_MODES {
        if mode == best.frp_mode {
            continue;
        }
        let h = Hyper { frp_mode: mode, ..best };
        let e = loyo(ctx, sensor, h, fit_years, threads);
        log.push(TuneEntry { sensor, hyper: h, ua: e.ua, r2deg: e.r2deg, r1deg: e.r1deg });
        if e.frp_r2deg.wape < best_e.frp_r2deg.wape {
            best = h;
            best_e = e;
        }
    }
    (best, best_e)
}

/// Вкладена LOYO: для кожного зовнішнього року гіперпараметри підбираються лише на інших роках.
/// Дає чесну (не оптимістичну) оцінку якості всієї процедури разом із підбором.
pub fn nested(ctx: &Ctx, sensor: Sensor, fit_years: &[i32], threads: usize) -> (Eval, Vec<(i32, Hyper)>) {
    let mut ys = Vec::new();
    let mut chosen = Vec::new();
    for &y in fit_years {
        let inner: Vec<i32> = fit_years.iter().copied().filter(|&x| x != y).collect();
        let mut dummy = Vec::new();
        let (h, _) = tune(ctx, sensor, &inner, threads, false, &mut dummy);
        let f = model::fit(ctx, sensor, h, &ctx.year_weights(&inner, None));
        ys.push(collect(ctx, &f, y));
        chosen.push((y, h));
        eprintln!("      вкладена CV {}: рік {y} → {:?}", sensor.name(), h);
    }
    (summarize(ys), chosen)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn metrics_basic() {
        let m = metrics(&[(100.0, 110.0), (50.0, 40.0), (0.0, 1.0)]);
        assert!((m.bias - (151.0 / 150.0 - 1.0)).abs() < 1e-12);
        assert!((m.wape - 21.0 / 150.0).abs() < 1e-12);
        assert!((m.mae - 7.0).abs() < 1e-12);
        assert!((m.mape - (0.1 + 0.2) / 2.0).abs() < 1e-12);
        let p = par_map(10, 3, |i| i * i);
        assert_eq!(p[9], 81);
    }
}
