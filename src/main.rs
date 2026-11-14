//! `harmonize` — hotspots.csv.gz (контракт 1) → cell_day.csv.gz + model_meta.json (контракт 2).

use fire_harmonizer::cv::{self, Eval, Metrics, TuneEntry};
use fire_harmonizer::data::{self, AQUA, LC_GROUP_NAMES, NOAA20, NOAA21, SAT_NAMES, SNPP, TERRA};
use fire_harmonizer::date::{self, fmt_date};
use fire_harmonizer::gz;
use fire_harmonizer::json::J;
use fire_harmonizer::model::{self, Ctx, Fitted, Hyper, Predictor, Sensor, BIN_NAMES, INTENSITY_NAMES, LEVEL_NAMES, PATTERN_NAMES};
use fire_harmonizer::rng::Rng;
use std::io::Write;
use std::time::Instant;

const HELP: &str = "\
harmonize — гармонізація MODIS ↔ VIIRS S-NPP (Fire Calendar)

ВИКОРИСТАННЯ:
  harmonize --input hotspots.csv.gz --out out/ [опції]

ОПЦІЇ:
  --input PATH          вхідний файл (контракт 1); .gz, звичайний CSV або '-' (stdin)  [hotspots.csv.gz]
  --out DIR             тека для результатів                                          [out]
  --model-version STR   версія моделі (model_meta.model_version)                      [v3]
  --persistent-mask M:Y маскувати мікроклітинки ~0,005° з детекціями в ≥M місяцях і ≥Y роках
                        (заводи без позначки type 2); 0 = вимкнути                    [30:8]
  --nested              додатково вкладена LOYO-крос-валідація (чесні метрики, ~×8 довше)
  --fit-years A-B       роки навчання (до дрейфу орбіт Terra/Aqua)                    [2012-2019]
  --viirs-drop-low      не брати VIIRS confidence=l в опорний ряд
  --no-tune             не підбирати гіперпараметри, а взяти вибрані на даних України
                        2012–2019 (у ~10 разів швидше; для інших регіонів краще з підбором)
  --bootstrap N         кількість bootstrap-реплік для довірчих інтервалів            [200]
  --gap-threshold X     мін. прогноз (точок/добу по Україні), щоб нуль S-NPP вважати збоєм  [15]
  --threads N           кількість потоків                                             [усі ядра]
  --seed N              зерно ГПВЧ для bootstrap                                      [42]
  -h, --help            ця довідка
";

struct Args {
    input: String,
    out: String,
    version: String,
    fit: (i32, i32),
    viirs_drop_low: bool,
    tune: bool,
    bootstrap: usize,
    gap_threshold: f64,
    threads: usize,
    seed: u64,
    nested: bool,
    persist: Option<(usize, usize)>,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        input: "hotspots.csv.gz".into(),
        out: "out".into(),
        version: "v3".into(),
        fit: (2012, 2019),
        viirs_drop_low: false,
        tune: true,
        bootstrap: 200,
        gap_threshold: 15.0,
        threads: std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4),
        seed: 42,
        nested: false,
        persist: Some((30, 8)),
    };
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        let mut val = |name: &str| it.next().ok_or_else(|| format!("{name}: потрібне значення"));
        match k.as_str() {
            "-h" | "--help" => {
                print!("{HELP}");
                std::process::exit(0);
            }
            "--input" => a.input = val(&k)?,
            "--out" => a.out = val(&k)?,
            "--model-version" => a.version = val(&k)?,
            "--fit-years" => {
                let v = val(&k)?;
                let (x, y) = v.split_once('-').ok_or("--fit-years: формат A-B, напр. 2012-2019")?;
                a.fit = (x.parse().map_err(|_| "--fit-years: рік")?, y.parse().map_err(|_| "--fit-years: рік")?);
            }
            "--viirs-drop-low" => a.viirs_drop_low = true,
            "--no-tune" => a.tune = false,
            "--bootstrap" => a.bootstrap = val(&k)?.parse().map_err(|_| "--bootstrap: число")?,
            "--gap-threshold" => a.gap_threshold = val(&k)?.parse().map_err(|_| "--gap-threshold: число")?,
            "--threads" => a.threads = val(&k)?.parse::<usize>().map_err(|_| "--threads: число")?.max(1),
            "--seed" => a.seed = val(&k)?.parse().map_err(|_| "--seed: число")?,
            "--nested" => a.nested = true,
            "--persistent-mask" => {
                let v = val(&k)?;
                a.persist = if v == "0" {
                    None
                } else {
                    let (m, y) = v.split_once(':').ok_or("--persistent-mask: формат M:Y, напр. 30:8, або 0")?;
                    Some((m.parse().map_err(|_| "--persistent-mask: число")?, y.parse().map_err(|_| "--persistent-mask: число")?))
                };
            }
            other => return Err(format!("невідома опція `{other}` (див. --help)")),
        }
    }
    Ok(a)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("ПОМИЛКА: {e}");
        std::process::exit(1);
    }
}

fn mj(m: &Metrics) -> J {
    J::obj()
        .set("n_units", m.n)
        .set("sum_obs", m.sum_obs)
        .set("sum_pred", m.sum_pred)
        .set("bias", m.bias)
        .set("wape", m.wape)
        .set("mae", m.mae)
        .set("mape", m.mape)
        .set("r2_log", m.r2_log)
}

fn ej(e: &Eval) -> J {
    J::obj()
        .set("ukraine_x_month", mj(&e.ua))
        .set("block2deg_x_month", mj(&e.r2deg))
        .set("block1deg_x_month", mj(&e.r1deg))
        .set("frp_block2deg_x_month", mj(&e.frp_r2deg))
        .set("annual_total_abs_error", e.annual_err)
        .set("seasonal_profile_error", e.profile_err)
        .set("coverage90_annual", J::Arr(vec![J::Int(e.cov_annual.0 as i64), J::Int(e.cov_annual.1 as i64)]))
        .set("coverage90_monthly", J::Arr(vec![J::Int(e.cov_month.0 as i64), J::Int(e.cov_month.1 as i64)]))
        .set(
            "subgroups_block2deg_x_month",
            J::Arr(e.subgroups.iter().map(|(g, n, m)| J::obj().set("group", g.as_str()).set("name", n.as_str()).set("wape", m.wape).set("bias", m.bias).set("mae", m.mae).set("n_units", m.n)).collect()),
        )
}

fn hj(h: &Hyper) -> J {
    J::obj()
        .set("lambda", h.lambda)
        .set("month_smoothing", h.smooth)
        .set("size_class_correction", h.bins)
        .set("modis_low_confidence_used", h.modis_low)
        .set("depth", h.depth)
        .set("leaf_level", LEVEL_NAMES[h.depth])
        .set("frp_mode", model::FRP_MODE_NAMES[h.frp_mode as usize])
        .set("terra_aqua_pattern", h.pattern)
        .set("frp_intensity", h.intensity)
        .set("landcover", h.lc)
        .set("missed_mass_allocation", model::MISS_ALLOC_NAMES[h.miss_alloc as usize])
        .set("missed_mass_em_model", model::MISS_MODEL_NAMES[h.miss_model as usize])
}

fn ranges(days: &[i32]) -> Vec<(i32, i32)> {
    let mut out: Vec<(i32, i32)> = Vec::new();
    for &d in days {
        match out.last_mut() {
            Some(r) if r.1 + 1 == d => r.1 = d,
            _ => out.push((d, d)),
        }
    }
    out
}

fn rj(r: &[(i32, i32)]) -> J {
    J::Arr(r.iter().map(|&(a, b)| J::Arr(vec![fmt_date(a).into(), fmt_date(b).into()])).collect())
}

fn io<T>(r: std::io::Result<T>, what: &str) -> Result<T, String> {
    r.map_err(|e| format!("{what}: {e}"))
}

fn run() -> Result<(), String> {
    let a = parse_args()?;
    let t_all = Instant::now();
    let diag = format!("{}/diagnostics", a.out);
    io(std::fs::create_dir_all(&diag), &diag)?;

    // ---------- 1. Дані ----------
    eprintln!("[1/7] Читаю {} …", a.input);
    let t = Instant::now();
    let rd = io(gz::open_read(&a.input), &a.input)?;
    let ds = data::load(rd, a.viirs_drop_low, a.persist)?;
    let st = &ds.stats;
    eprintln!(
        "      рядків {}, некоректних {}, type 2/3 виключено {}, NRT без type {} (з них замасковано {}), записів клітинка×доба×Д/Н {} — {:.1} с",
        st.lines, st.bad_lines, st.type_excluded, st.nrt_untyped, st.nrt_masked, ds.cell.len(), t.elapsed().as_secs_f64()
    );
    if let Some((m, y)) = a.persist {
        eprintln!(
            "      маска стабільних джерел (≥{m} місяців, ≥{y} років): {} мікроклітинок, прибрано {} точок",
            st.persistent_cells, st.persistent_masked
        );
    }
    for s in 0..5 {
        if let (Some(f), Some(l)) = (ds.first[s], ds.last[s]) {
            eprintln!("      {:<7} {} … {}  ({} точок)", SAT_NAMES[s], fmt_date(f), fmt_date(l), st.by_sat[s]);
        }
    }
    let t0 = ds.first[TERRA].ok_or("у даних немає Terra")?;
    let t1 = ds.last[TERRA].unwrap();
    let a0 = ds.first[AQUA].ok_or("у даних немає Aqua")?;
    let a1 = ds.last[AQUA].unwrap();
    let v0 = ds.first[SNPP].ok_or("у даних немає S-NPP — нема на чому навчати модель")?;
    let v1 = ds.last[SNPP].unwrap();
    let ov_start = v0.max(a0);
    let ov_end = v1.min(t1).min(a1);
    if ov_start > ov_end {
        return Err("немає періоду, коли одночасно працюють Terra, Aqua і S-NPP".into());
    }

    let mut ctx = Ctx::new(&ds);
    let n = ds.cell.len();
    let dmin = *ds.day.iter().min().unwrap();
    let dmax = *ds.day.iter().max().unwrap();
    let nd = (dmax - dmin + 1) as usize;
    let di = |d: i32| (d - dmin) as usize;
    let mut daily_v = vec![0.0f64; nd];
    let mut daily_m = vec![0.0f64; nd];
    for i in 0..n {
        daily_v[di(ds.day[i])] += ds.agg[i].v as f64;
        daily_m[di(ds.day[i])] += ds.agg[i].modis_raw() as f64;
    }

    // ---------- 2. Період перекриття і збої ----------
    eprintln!("[2/7] Період перекриття {} … {}; шукаю дні збою сенсорів …", fmt_date(ov_start), fmt_date(ov_end));
    let (oy0, _, _) = date::civil_from_days(ov_start);
    let (oy1, _, _) = date::civil_from_days(ov_end);
    let fit_years: Vec<i32> = (a.fit.0..=a.fit.1).filter(|&y| y >= oy0 && y <= oy1).collect();
    if fit_years.len() < 3 {
        return Err(format!("замало років навчання в періоді перекриття: {fit_years:?} (потрібно ≥ 3)"));
    }
    let stress_years: Vec<i32> = (oy0..=oy1).filter(|y| !fit_years.contains(y)).collect();

    ctx.set_overlap(|d| d >= ov_start && d <= ov_end);
    let yw_fit = ctx.year_weights(&fit_years, None);
    let pre = model::fit(&ctx, Sensor::TerraAqua, Hyper::default(), &yw_fit);
    let mut daily_pred = vec![0.0f64; nd];
    for i in 0..n {
        let d = ds.day[i];
        if d >= v0 && d <= v1 {
            daily_pred[di(d)] += pre.predict(&ctx, i).0;
        }
    }
    let mut gap_v = vec![false; nd];
    let mut gap_m = vec![false; nd];
    for d in v0..=v1 {
        if daily_v[di(d)] == 0.0 && daily_pred[di(d)] >= a.gap_threshold {
            gap_v[di(d)] = true;
        }
    }
    for d in ov_start..=ov_end {
        if daily_m[di(d)] == 0.0 && daily_v[di(d)] >= 3.0 * a.gap_threshold {
            gap_m[di(d)] = true;
        }
    }
    let gap_v_days: Vec<i32> = (v0..=v1).filter(|&d| gap_v[di(d)]).collect();
    let gap_m_days: Vec<i32> = (ov_start..=ov_end).filter(|&d| gap_m[di(d)]).collect();
    // Підозра на відсутність даних до появи S-NPP: ≥3 доби поспіль у квітні–вересні без жодної точки MODIS.
    let mut zero_pre: Vec<i32> = Vec::new();
    for d in t0..v0.min(dmax + 1) {
        let (_, m, _) = date::civil_from_days(d);
        if (4..=9).contains(&m) && daily_m[di(d)] == 0.0 {
            zero_pre.push(d);
        }
    }
    let suspected: Vec<(i32, i32)> = ranges(&zero_pre).into_iter().filter(|r| r.1 - r.0 >= 2).collect();
    eprintln!(
        "      днів збою S-NPP: {}, днів збою MODIS у перекритті: {}, підозрілих прогалин до {}: {}",
        gap_v_days.len(),
        gap_m_days.len(),
        fmt_date(v0),
        suspected.len()
    );
    ctx.set_overlap(|d| d >= ov_start && d <= ov_end && !gap_v[di(d)] && !gap_m[di(d)]);

    // ---------- 3. Підбір гіперпараметрів ----------
    eprintln!(
        "[3/7] LOYO-крос-валідація на {:?}{} (земний покрив у даних: {}) …",
        fit_years,
        if a.tune { " і підбір гіперпараметрів" } else { ", гіперпараметри фіксовані (--no-tune)" },
        if ds.has_lc { "так" } else { "ні" }
    );
    let t = Instant::now();
    let mut tlog: Vec<TuneEntry> = Vec::new();
    let select = |sensor: Sensor, tlog: &mut Vec<TuneEntry>| -> (Hyper, Eval) {
        if a.tune {
            cv::tune(&ctx, sensor, &fit_years, a.threads, true, tlog)
        } else {
            let h = Hyper::v3(sensor);
            let e = cv::loyo(&ctx, sensor, h, &fit_years, a.threads);
            tlog.push(TuneEntry { sensor, hyper: h, ua: e.ua, r2deg: e.r2deg, r1deg: e.r1deg });
            (h, e)
        }
    };
    let (h_ta, ev_ta) = select(Sensor::TerraAqua, &mut tlog);
    let (h_t, ev_t) = select(Sensor::Terra, &mut tlog);
    let base_ta = cv::loyo_const(&ctx, Sensor::TerraAqua, &fit_years, a.threads);
    let base_t = cv::loyo_const(&ctx, Sensor::Terra, &fit_years, a.threads);
    eprintln!("      вибрано Terra+Aqua: {:?}", h_ta);
    eprintln!("      вибрано Terra:      {:?}", h_t);
    eprintln!(
        "      WAPE 2°×міс: константа {:.4} → модель {:.4} (Terra+Aqua);  константа {:.4} → модель {:.4} (Terra) — {:.1} с",
        base_ta.r2deg.wape,
        ev_ta.r2deg.wape,
        base_t.r2deg.wape,
        ev_t.r2deg.wape,
        t.elapsed().as_secs_f64()
    );

    let nested = if a.nested {
        eprintln!("[3b] Вкладена LOYO-крос-валідація (підбір усередині кожного фолду) …");
        let t = Instant::now();
        let (e_ta, c_ta) = cv::nested(&ctx, Sensor::TerraAqua, &fit_years, a.threads);
        let (e_t, c_t) = cv::nested(&ctx, Sensor::Terra, &fit_years, a.threads);
        eprintln!(
            "      вкладена: WAPE 2°×міс {:.4} (Terra+Aqua), {:.4} (Terra) — {:.1} с",
            e_ta.r2deg.wape,
            e_t.r2deg.wape,
            t.elapsed().as_secs_f64()
        );
        Some((e_ta, c_ta, e_t, c_t))
    } else {
        None
    };

    // ---------- 4. Фінальна модель і стрес-тест ----------
    eprintln!("[4/7] Фінальне навчання і стрес-тест на {:?} …", stress_years);
    let fit_ta = model::fit(&ctx, Sensor::TerraAqua, h_ta, &yw_fit);
    let fit_t = model::fit(&ctx, Sensor::Terra, h_t, &yw_fit);
    let stress = cv::eval_years(&ctx, &fit_ta, &stress_years, a.threads);

    let model_day = |d: i32| d < v0 || d > v1 || gap_v[di(d)];
    let pick = |d: i32| -> &Fitted {
        if d < a0 {
            &fit_t
        } else {
            &fit_ta
        }
    };

    // ---------- 5. Bootstrap ----------
    // Інтервал = невизначеність коефіцієнтів (bootstrap по роках) × залишкова мінливість
    // окремого року (σ з LOYO: RMS ln(pred/obs) річних сум). Без другого множника інтервал
    // завузький (90 %-інтервал накривав лише ~половину років).
    let sigma_of = |e: &Eval| -> f64 {
        let v: Vec<f64> = e.years.iter().filter(|y| y.tot.0 > 0.0 && y.tot.1 > 0.0).map(|y| (y.tot.1 / y.tot.0).ln()).collect();
        if v.len() < 2 {
            0.0
        } else {
            (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt()
        }
    };
    // Поправка на малу вибірку: σ з ~8 років → квантиль Стьюдента і √(1+1/n) замість z = 1,645.
    let nres = |e: &Eval| e.years.iter().filter(|y| y.tot.0 > 0.0 && y.tot.1 > 0.0).count();
    let (sig_ta, sig_t) = (
        sigma_of(&ev_ta) * cv::pred_mult(nres(&ev_ta)) / cv::Z90,
        sigma_of(&ev_t) * cv::pred_mult(nres(&ev_t)) / cv::Z90,
    );
    eprintln!("[5/7] Bootstrap по роках: {} реплік; σ року: Terra+Aqua {:.4}, Terra {:.4} …", a.bootstrap, sig_ta, sig_t);
    let t = Instant::now();
    let ny = ctx.ny;
    let boots: Vec<Vec<f64>> = cv::par_map(a.bootstrap, a.threads, |b| {
        let mut rng = Rng::new(a.seed.wrapping_add(7919 * b as u64));
        let mut yw = vec![0.0; ny];
        for _ in 0..fit_years.len() {
            let y = fit_years[rng.below(fit_years.len())];
            yw[(y - ctx.y0) as usize] += 1.0;
        }
        let fta = model::fit(&ctx, Sensor::TerraAqua, h_ta, &yw);
        let ft = model::fit(&ctx, Sensor::Terra, h_t, &yw);
        let mut tot_ta = vec![0.0; ny];
        let mut tot_t = vec![0.0; ny];
        for i in 0..n {
            let d = ds.day[i];
            if model_day(d) {
                let yk = (ctx.year[i] - ctx.y0) as usize;
                if d < a0 {
                    tot_t[yk] += ft.predict(&ctx, i).0;
                } else {
                    tot_ta[yk] += fta.predict(&ctx, i).0;
                }
            }
        }
        (0..ny).map(|k| tot_ta[k] * (sig_ta * rng.normal()).exp() + tot_t[k] * (sig_t * rng.normal()).exp()).collect()
    });
    eprintln!("      {:.1} с", t.elapsed().as_secs_f64());

    // ---------- 6. cell_day.csv.gz ----------
    let cell_day_path = format!("{}/cell_day.csv.gz", a.out);
    eprintln!("[6/7] Пишу {cell_day_path} …");
    let t = Instant::now();
    let mut w = io(gz::create(&cell_day_path), &cell_day_path)?;
    io(writeln!(w, "cell_id,date,modis_count,viirs_count,count_h,frp_h"), &cell_day_path)?;
    let mut rows = 0u64;
    // річні агрегати: [count_h, frp_h, model_part, observed_part, modis_based, snpp_obs]
    let mut annual = vec![[0.0f64; 6]; ny];
    let mut i = 0;
    while i < n {
        let (c, d) = (ds.cell[i], ds.day[i]);
        let md = model_day(d);
        let f = pick(d);
        let (mut mc, mut vc, mut ch, mut fh, mut mb) = (0u32, 0u32, 0.0f64, 0.0f64, 0.0f64);
        let mut j = i;
        while j < n && ds.cell[j] == c && ds.day[j] == d {
            let ag = &ds.agg[j];
            mc += ag.modis_raw();
            vc += ag.v;
            let (p, pf) = f.predict(&ctx, j);
            mb += p;
            if md {
                ch += p;
                fh += pf;
            } else {
                ch += ag.v as f64;
                fh += ag.fv;
            }
            j += 1;
        }
        let yk = (ctx.year[i] - ctx.y0) as usize;
        annual[yk][0] += ch;
        annual[yk][1] += fh;
        annual[yk][if md { 2 } else { 3 }] += ch;
        annual[yk][4] += mb;
        annual[yk][5] += vc as f64;
        if mc > 0 || vc > 0 {
            io(writeln!(w, "{c},{},{mc},{vc},{ch:.3},{fh:.2}", fmt_date(d)), &cell_day_path)?;
            rows += 1;
        }
        i = j;
    }
    io(w.finish(), &cell_day_path)?;
    eprintln!("      рядків {rows} — {:.1} с", t.elapsed().as_secs_f64());

    // ---------- 7. Діагностика і model_meta.json ----------
    eprintln!("[7/7] Діагностика і model_meta.json …");
    let q = |v: &mut Vec<f64>, p: f64| -> f64 {
        if v.is_empty() {
            return f64::NAN;
        }
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let x = p * (v.len() - 1) as f64;
        let lo = x.floor() as usize;
        let hi = x.ceil() as usize;
        v[lo] + (v[hi] - v[lo]) * (x - lo as f64)
    };
    let mut csv = String::from(
        "year,terra,aqua,snpp,noaa20,noaa21,count_h,count_h_p05,count_h_p95,frp_h,count_from_modis,snpp_observed,model_share\n",
    );
    let mut ci_json: Vec<J> = Vec::new();
    for (yk, an) in annual.iter().enumerate() {
        let y = ctx.y0 + yk as i32;
        let raw = ds.raw_year.get(&y).copied().unwrap_or([0; 5]);
        let mut bs: Vec<f64> = boots.iter().map(|b| b[yk] + an[3]).collect();
        let (p05, p95) = if an[2] > 0.0 && !bs.is_empty() { (q(&mut bs, 0.05), q(&mut bs, 0.95)) } else { (an[0], an[0]) };
        let share = if an[0] > 0.0 { an[2] / an[0] } else { 0.0 };
        csv.push_str(&format!(
            "{y},{},{},{},{},{},{:.1},{:.1},{:.1},{:.1},{:.1},{:.0},{:.3}\n",
            raw[TERRA], raw[AQUA], raw[SNPP], raw[NOAA20], raw[NOAA21], an[0], p05, p95, an[1], an[4], an[5], share
        ));
        if an[2] > 0.0 {
            ci_json.push(
                J::obj().set("year", y).set("count_h", an[0]).set("p05", p05).set("p95", p95).set("model_share", share),
            );
        }
    }
    let p = format!("{diag}/annual.csv");
    io(std::fs::write(&p, csv), &p)?;

    // відношення по роках (дрейф орбіт): по днях перекриття без збоїв
    let mut ry = vec![[0.0f64; 5]; ny]; // V, T, A, FV, FM
    for i in 0..n {
        let d = ds.day[i];
        if d >= ov_start && d <= ov_end && !gap_v[di(d)] && !gap_m[di(d)] {
            let ag = &ds.agg[i];
            let r = &mut ry[(ctx.year[i] - ctx.y0) as usize];
            r[0] += ag.v as f64;
            r[1] += (ag.t_hi + ag.t_lo) as f64;
            r[2] += (ag.a_hi + ag.a_lo) as f64;
            r[3] += ag.fv;
            r[4] += ag.ft_hi + ag.ft_lo + ag.fa_hi + ag.fa_lo;
        }
    }
    let pred_by_year: Vec<(i32, f64, f64, &str)> = ev_ta
        .years
        .iter()
        .map(|y| (y.year, y.tot.0, y.tot.1, "loyo"))
        .chain(stress.years.iter().map(|y| (y.year, y.tot.0, y.tot.1, "stress")))
        .collect();
    let mut csv = String::from("year,snpp,terra,aqua,k_terra_aqua,k_terra,k_aqua,k_frp,model_pred,pred_over_obs,eval\n");
    let mut ratio_json: Vec<J> = Vec::new();
    for (yk, r) in ry.iter().enumerate() {
        if r[0] == 0.0 {
            continue;
        }
        let y = ctx.y0 + yk as i32;
        let div = |a: f64, b: f64| if b > 0.0 { a / b } else { f64::NAN };
        let pr = pred_by_year.iter().find(|p| p.0 == y);
        let (mp, ratio, ev) = match pr {
            Some(p) => (p.2, div(p.2, p.1), p.3),
            None => (f64::NAN, f64::NAN, ""),
        };
        csv.push_str(&format!(
            "{y},{:.0},{:.0},{:.0},{:.4},{:.4},{:.4},{:.4},{:.1},{:.4},{ev}\n",
            r[0],
            r[1],
            r[2],
            div(r[0], r[1] + r[2]),
            div(r[0], r[1]),
            div(r[0], r[2]),
            div(r[3], r[4]),
            mp,
            ratio
        ));
        ratio_json.push(
            J::obj()
                .set("year", y)
                .set("k_terra_aqua", div(r[0], r[1] + r[2]))
                .set("k_terra", div(r[0], r[1]))
                .set("k_aqua", div(r[0], r[2]))
                .set("pred_over_obs", ratio)
                .set("eval", ev),
        );
    }
    let p = format!("{diag}/ratio_by_year.csv");
    io(std::fs::write(&p, csv), &p)?;

    let mut csv = String::from("year,month,eval,snpp_observed,predicted_from_modis\n");
    for (set, e) in [("loyo", &ev_ta), ("stress", &stress)] {
        for y in &e.years {
            for (m, (o, pp)) in y.ua.iter().enumerate() {
                csv.push_str(&format!("{},{},{set},{o:.0},{pp:.1}\n", y.year, m + 1));
            }
        }
    }
    let p = format!("{diag}/validation_monthly.csv");
    io(std::fs::write(&p, csv), &p)?;

    let mut csv = String::from("model,lambda,smooth,bins,modis_low,depth,pattern,intensity,landcover,miss_alloc,miss_model,frp_mode,wape_ua_month,wape_2deg_month,wape_1deg_month,bias_2deg_month,r2log_2deg_month\n");
    for e in &tlog {
        let h = e.hyper;
        csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{},{},{:.5},{:.5},{:.5},{:.5},{:.5}\n",
            e.sensor.name(),
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
            h.frp_mode,
            e.ua.wape,
            e.r2deg.wape,
            e.r1deg.wape,
            e.r2deg.bias,
            e.r2deg.r2_log
        ));
    }
    let p = format!("{diag}/cv_grid.csv");
    io(std::fs::write(&p, csv), &p)?;

    let mut csv = String::from("model,eval,group,name,n_units,wape,bias,mae,r2_log\n");
    let mut sets: Vec<(&str, &str, &Eval)> = vec![("terra_aqua", "loyo", &ev_ta), ("terra_only", "loyo", &ev_t), ("terra_aqua", "stress", &stress)];
    if let Some((e_ta, _, e_t, _)) = &nested {
        sets.push(("terra_aqua", "nested", e_ta));
        sets.push(("terra_only", "nested", e_t));
    }
    for (m, ev, e) in sets {
        for (g, nme, mt) in &e.subgroups {
            csv.push_str(&format!("{m},{ev},{g},{nme},{},{:.5},{:.5},{:.3},{:.5}\n", mt.n, mt.wape, mt.bias, mt.mae, mt.r2_log));
        }
    }
    let p = format!("{diag}/subgroups.csv");
    io(std::fs::write(&p, csv), &p)?;

    let mut csv = String::from("model,leaf_level,landcover,block_lon_min,block_lat_min,block_deg,daynight,month,k,q,rho,q_frp,phi\n");
    for f in [&fit_ta, &fit_t] {
        let lv = ctx.levels(f.hyper.lc);
        let dep = model::leaf_level(&f.hyper);
        for s in 0..lv.n[dep] {
            if lv.parent[dep][s] == u32::MAX {
                continue;
            }
            let nm = lv.nm[dep];
            let (region, dn, mo) = (s / (2 * nm), (s / nm) % 2, s % nm);
            let (bx, by, lcg) = lv.regions[dep].get(region).copied().unwrap_or((u32::MAX, u32::MAX, 255));
            let sz = lv.block[dep];
            let (lon, lat, deg) = if sz > 0 {
                ((bx * sz) as f64 / 10.0 - 180.0, (by * sz) as f64 / 10.0 - 90.0, sz as f64 / 10.0)
            } else {
                (f64::NAN, f64::NAN, f64::NAN)
            };
            let lcn = if lcg < 6 { LC_GROUP_NAMES[lcg as usize] } else { "усі" };
            csv.push_str(&format!(
                "{},{},{},{lon:.1},{lat:.1},{deg:.1},{},{},{:.5},{:.5},{:.5},{:.5},{:.5}\n",
                f.sensor.name(),
                LEVEL_NAMES[f.hyper.depth],
                lcn,
                if dn == 0 { "D" } else { "N" },
                if nm == 12 { (mo + 1).to_string() } else { "all".into() },
                f.k[s],
                f.q[s],
                f.rho[s],
                f.qf[s],
                f.phi[s]
            ));
        }
    }
    let p = format!("{diag}/coefficients.csv");
    io(std::fs::write(&p, csv), &p)?;

    let beta_j = |f: &Fitted| -> J {
        let names = |n: &[&str]| J::Arr(n.iter().map(|&s| s.into()).collect());
        J::obj()
            .set(
                "size",
                J::obj()
                    .set("bins_modis_points_per_cell_day", names(&BIN_NAMES))
                    .set("day", J::Arr(f.beta_s[0].iter().map(|&x| J::Num(x)).collect()))
                    .set("night", J::Arr(f.beta_s[1].iter().map(|&x| J::Num(x)).collect())),
            )
            .set(
                "terra_aqua_pattern",
                J::obj()
                    .set("levels", names(&PATTERN_NAMES))
                    .set("day", J::Arr(f.beta_p[0].iter().map(|&x| J::Num(x)).collect()))
                    .set("night", J::Arr(f.beta_p[1].iter().map(|&x| J::Num(x)).collect())),
            )
            .set(
                "frp_intensity",
                J::obj()
                    .set("levels", names(&INTENSITY_NAMES))
                    .set("thresholds_mw_per_detection", J::Arr(f.thr.iter().map(|&x| J::Num(x)).collect()))
                    .set("day", J::Arr(f.beta_i[0].iter().map(|&x| J::Num(x)).collect()))
                    .set("night", J::Arr(f.beta_i[1].iter().map(|&x| J::Num(x)).collect())),
            )
    };
    let fit_j = |f: &Fitted| -> J {
        J::obj()
            .set("hyperparameters", hj(&f.hyper))
            .set("k_pooled", f.k_pooled)
            .set("missed_share_pooled", f.q_pooled)
            .set("frp_ratio_pooled", f.rho_pooled)
            .set("multipliers_beta", beta_j(f))
            .set("calibration_n_strata", f.cal.len())
            .set(
                "missed_mass_gamma",
                J::obj()
                    .set("intensity", J::obj().set("levels", J::Arr(INTENSITY_NAMES.iter().map(|&s| s.into()).collect())).set("day", J::Arr(f.gm_i[0].iter().map(|&x| J::Num(x)).collect())).set("night", J::Arr(f.gm_i[1].iter().map(|&x| J::Num(x)).collect())))
                    .set("neighbourhood", J::obj().set("levels", J::Arr(model::NB_NAMES.iter().map(|&s| s.into()).collect())).set("day", J::Arr(f.gm_n[0].iter().map(|&x| J::Num(x)).collect())).set("night", J::Arr(f.gm_n[1].iter().map(|&x| J::Num(x)).collect())))
                    .set("size", J::obj().set("levels", J::Arr(BIN_NAMES.iter().map(|&s| s.into()).collect())).set("day", J::Arr(f.gm_s[0].iter().map(|&x| J::Num(x)).collect())).set("night", J::Arr(f.gm_s[1].iter().map(|&x| J::Num(x)).collect()))),
            )
            .set("frp_intensity_gamma", J::obj().set("day", J::Arr(f.gf_i[0].iter().map(|&x| J::Num(x)).collect())).set("night", J::Arr(f.gf_i[1].iter().map(|&x| J::Num(x)).collect())))
    };

    let total_model_days = (dmin..=dmax).filter(|&d| model_day(d) && d >= t0).count();
    let notes = format!(
        "count_h і frp_h виражені в шкалі VIIRS S-NPP (375 м). \
З {v0s} по {v1s} count_h = спостережений S-NPP; виняток — {ngv} діб збою S-NPP, для них підставлено модель. \
З {a0s} до {v0s} значення пораховано моделлю з Terra+Aqua, з {t0s} по {a0e} — окремою моделлю лише з Terra (у цей час Aqua ще не працювала). \
Модель v3 — ієрархічне відношення VIIRS/MODIS зі стягуванням (день/ніч → місяць → [земний покрив] → блоки 4°/2°/1°/0,5°), \
множниками за розміром події, патерном Terra/Aqua та інтенсивністю FRP (які саме увімкнені — див. model.*.hyperparameters), \
часткою «пропущених» дрібних пожеж, яких не бачить MODIS, і калібруванням сум на рівні день/ніч × місяць. \
Ці пропущені пожежі розподілено по клітинках із MODIS у тій самій страті, тому суми по області й місяцю незміщені, \
а окрема клітинко-доба до 2012 року — лише статистична оцінка. Через це до 2012 року менше діб «з вогнем» на малих полігонах, \
і для малих областей календар краще агрегувати по тижнях. \
Модель навчено на {fy0}–{fy1} — до дрейфу орбіт: Terra дрейфує з 2020, Aqua з 2021; \
у роки стрес-тесту {sy} співвідношення змінюється через дрейф і війну. \
Метрики в `metrics` — LOYO-крос-валідація на роках навчання, агрегат «блок 2° × місяць». \
type 2/3 виключено; в NRT, де type не заповнено, промислові джерела прибрано маскою з SP ({mask} мікроклітинок, {nm} точок). \
Маска стабільних джерел за даними (мікроклітинки з детекціями в багатьох місяцях і роках): {pc} мікроклітинок, {pm} точок. \
NOAA-20/21 у ряд не входять. Підозрілі прогалини в даних до {v0s} див. `suspected_no_data_ranges`: нуль у ці дні може означати «немає даних».",
        v0s = fmt_date(v0),
        v1s = fmt_date(v1),
        ngv = gap_v_days.len(),
        a0s = fmt_date(a0),
        t0s = fmt_date(t0),
        a0e = fmt_date(a0 - 1),
        fy0 = fit_years[0],
        fy1 = fit_years[fit_years.len() - 1],
        sy = if stress_years.is_empty() {
            "—".to_string()
        } else {
            format!("{}–{}", stress_years[0], stress_years[stress_years.len() - 1])
        },
        mask = st.static_mask_size,
        pc = st.persistent_cells,
        pm = st.persistent_masked,
        nm = st.nrt_masked,
    );

    let meta = J::obj()
        .set("model_version", a.version.as_str())
        .set("reference_scale", "VIIRS_SNPP")
        .set("fit_years", J::Arr(fit_years.iter().map(|&y| J::Int(y as i64)).collect()))
        .set("holdout_years", J::Arr(stress_years.iter().map(|&y| J::Int(y as i64)).collect()))
        .set("metrics", J::obj().set("mape", ev_ta.r2deg.mape).set("r2", ev_ta.r2deg.r2_log).set("bias", ev_ta.r2deg.bias))
        .set("notes", notes)
        .set(
            "metrics_definition",
            "LOYO: для кожного року навчання модель вчиться на решті й прогнозує його S-NPP лише з MODIS. \
Одиниця — блок 2°×місяць (сума всіх точок S-NPP, зокрема там, де MODIS не бачив нічого). \
bias = Σpred/Σobs − 1; r2 = R² на ln(1+x); mape — по одиницях зі спостереженням ≥ 20; wape = Σ|pred−obs|/Σobs.",
        )
        .set(
            "periods",
            J::obj()
                .set("terra_start", fmt_date(t0))
                .set("aqua_start", fmt_date(a0))
                .set("snpp_start", fmt_date(v0))
                .set("snpp_end", fmt_date(v1))
                .set("overlap_start", fmt_date(ov_start))
                .set("overlap_end", fmt_date(ov_end))
                .set("terra_only", J::Arr(vec![fmt_date(t0).into(), fmt_date(a0 - 1).into()]))
                .set("model_days_total", total_model_days),
        )
        .set(
            "filters",
            J::obj()
                .set("excluded_types", J::Arr(vec![J::Int(2), J::Int(3)]))
                .set("nrt_static_mask", "мікроклітинки 0,005° з SP-точок type 2/3, сусідство 3×3")
                .set("viirs_satellite", "SNPP")
                .set("viirs_low_confidence_dropped", a.viirs_drop_low)
                .set("modis_count_column", "Terra + Aqua, усі confidence (сирі)"),
        )
        .set("model", J::obj().set("terra_aqua", fit_j(&fit_ta)).set("terra_only", fit_j(&fit_t)))
        .set(
            "cv",
            J::obj()
                .set("terra_aqua", ej(&ev_ta))
                .set("terra_only", ej(&ev_t))
                .set("baseline_v0_terra_aqua", ej(&base_ta))
                .set("baseline_v0_terra_only", ej(&base_t))
                .set("selection_criterion", "мінімум WAPE на рівні блок 2° × місяць")
                .set("configs_evaluated", tlog.len()),
        )
        .set(
            "cv_nested",
            match &nested {
                None => J::Null,
                Some((e_ta, c_ta, e_t, c_t)) => {
                    let ch = |c: &Vec<(i32, Hyper)>| J::Arr(c.iter().map(|(y, h)| J::obj().set("held_out_year", *y).set("hyperparameters", hj(h))).collect());
                    J::obj()
                        .set("description", "вкладена LOYO: гіперпараметри підбираються без відкладеного року — чесна оцінка")
                        .set("terra_aqua", ej(e_ta))
                        .set("terra_only", ej(e_t))
                        .set("chosen_terra_aqua", ch(c_ta))
                        .set("chosen_terra_only", ch(c_t))
                }
            },
        )
        .set("hyperparameters_tuned", a.tune)
        .set("stress_test", J::obj().set("years", J::Arr(stress_years.iter().map(|&y| J::Int(y as i64)).collect())).set("terra_aqua", ej(&stress)))
        .set("ratio_by_year", J::Arr(ratio_json))
        .set("uncertainty_bootstrap", J::obj().set("replicates", a.bootstrap).set("residual_sigma_log_terra_aqua", sig_ta).set("residual_sigma_log_terra_only", sig_t).set("method", "bootstrap по роках навчання × lognormal залишок року (σ з LOYO, з поправкою Стьюдента на малу вибірку); 5–95 перцентилі річних сум").set("annual", J::Arr(ci_json)))
        .set("viirs_gap_days", J::Arr(gap_v_days.iter().map(|&d| fmt_date(d).into()).collect()))
        .set("modis_gap_days", J::Arr(gap_m_days.iter().map(|&d| fmt_date(d).into()).collect()))
        .set("suspected_no_data_ranges", rj(&suspected))
        .set(
            "load_stats",
            J::obj()
                .set("lines", st.lines)
                .set("bad_lines", st.bad_lines)
                .set("type_2_3_excluded", st.type_excluded)
                .set("nrt_untyped", st.nrt_untyped)
                .set("nrt_masked", st.nrt_masked)
                .set("viirs_low_dropped", st.viirs_low_dropped)
                .set("static_mask_cells", st.static_mask_size)
                .set("persistent_mask", match a.persist { Some((m, y)) => J::Arr(vec![J::Int(m as i64), J::Int(y as i64)]), None => J::Null })
                .set("persistent_mask_cells", st.persistent_cells)
                .set("persistent_mask_points_removed", st.persistent_masked)
                .set("cell_day_rows", rows),
        )
        .set("runtime_seconds", t_all.elapsed().as_secs_f64());
    let p = format!("{}/model_meta.json", a.out);
    io(std::fs::write(&p, meta.render()), &p)?;

    eprintln!(
        "Готово за {:.1} с. LOYO 2°×міс: WAPE {:.3}, bias {:+.3}, R²(log) {:.3}, MAPE {:.3}",
        t_all.elapsed().as_secs_f64(),
        ev_ta.r2deg.wape,
        ev_ta.r2deg.bias,
        ev_ta.r2deg.r2_log,
        ev_ta.r2deg.mape
    );
    Ok(())
}
