//! Модель гармонізації MODIS → шкала VIIRS S-NPP (версія v3).
//!
//! Для запису i (клітинка 0,1° × доба UTC × день/ніч) з `m` точками MODIS:
//!
//! ```text
//! count_h = ( m · k[s] · βₛ · βₚ · βᵢ   +   a(m) · g[s] · γᵢ · γₙ · γₛ ) · c[день/ніч, місяць]
//!             └─ детектована частина ─┘     └──── пропущена маса ────┘
//! ```
//!
//! - `s` — листова страта ієрархії: день/ніч × місяць × [земний покрив] × просторовий блок;
//! - `k[s]` — відношення VIIRS/MODIS там, де MODIS щось бачив;
//! - `βₛ` — розмір події (1, 2, 3–4, 5–8, 9–16, 17+ точок MODIS); `βₚ` — патерн Terra/Aqua;
//!   `βᵢ` — інтенсивність (середній FRP однієї детекції MODIS, терцилі);
//! - `g[s]` — «пропущена маса»: точки VIIRS у клітинко-добах, де MODIS не бачив нічого
//!   (~59 % сигналу), на одиницю ваги розподілу `a(m)`;
//! - `γ` — множники пропущеної маси за інтенсивністю, сусідством (інші клітинко-доби з MODIS
//!   у вікні 3×3 × ±1 доба) і розміром, оцінені EM на одиницях «страта × рік»;
//! - `c` — калібрування сум на навчальних роках.
//!
//! Коефіцієнти страт оцінюються ієрархічно зі стягуванням до батьківської страти
//! (гамма-пуассонівська оцінка): `r_L[s] = (N + λ·r_{L−1}[parent]) / (D + λ)`.
//! Множники β — ітеративним пропорційним підганянням (IPF), що еквівалентно пуассонівській
//! GLM з лог-зв'язком `log μ = log m + α_s + βₛ + βₚ + βᵢ`.
//!
//! Які саме складові вмикати, вирішує крос-валідація (див. [`crate::cv::tune`]); набір,
//! вибраний на даних України 2012–2019, — у [`Hyper::v3_terra_aqua`] і [`Hyper::v3_terra_only`].

use crate::data::{Agg, Dataset, FxMap};
use crate::date;

pub const NB: usize = 6;
pub const NP: usize = 3;
pub const NI: usize = 3;
/// Класи сусідства: скільки інших клітинко-діб з MODIS у вікні 3×3 клітинки × ±1 доба.
pub const NNB: usize = 4;
pub const NB_NAMES: [&str; NNB] = ["0", "1–2", "3–6", "7+"];
pub const BIN_NAMES: [&str; NB] = ["1", "2", "3-4", "5-8", "9-16", "17+"];
pub const PATTERN_NAMES: [&str; NP] = ["лише Aqua", "лише Terra", "Terra і Aqua"];
pub const INTENSITY_NAMES: [&str; NI] = ["низька", "середня", "висока"];
/// Назви «просторової глибини» (однакові з земним покривом і без нього).
pub const LEVEL_NAMES: [&str; 6] =
    ["день/ніч", "день/ніч×місяць", "4°×місяць", "2°×місяць", "1°×місяць", "0,5°×місяць"];
pub const MAX_DEPTH: usize = 5;
/// Рівні ієрархії **без** земного покриву, на яких рахуються регіональні метрики.
pub const EVAL_LEVEL_2DEG: usize = 3;
pub const EVAL_LEVEL_1DEG: usize = 4;
pub const FRP_MODES: u8 = 5;
pub const FRP_MODE_NAMES: [&str; 5] = [
    "ratio_to_modis_frp",
    "count_h_x_mean_viirs_frp",
    "constant_ratio",
    "hierarchical_total_ratio",
    "separate_detected_and_missed_x_intensity",
];
pub const MISS_ALLOC_NAMES: [&str; 4] = ["proportional_to_detected", "proportional_to_m", "proportional_to_sqrt_m", "presence"];
pub const MISS_MODEL_NAMES: [&str; 3] = ["off", "unit_stratum_x_day", "unit_stratum_x_year"];

#[inline]
pub fn nb_bin(n: u16) -> usize {
    match n {
        0 => 0,
        1..=2 => 1,
        3..=6 => 2,
        _ => 3,
    }
}

#[inline]
pub fn bin_of(m: u32) -> usize {
    match m {
        0 | 1 => 0,
        2 => 1,
        3..=4 => 2,
        5..=8 => 3,
        9..=16 => 4,
        _ => 5,
    }
}

#[inline]
fn alloc_w(mode: u8, m: f64) -> f64 {
    match mode {
        1 => m,
        2 => m.sqrt(),
        _ => 1.0,
    }
}

#[inline]
fn ibin(r: f64, thr: &[f64; 2]) -> usize {
    if r < thr[0] {
        0
    } else if r < thr[1] {
        1
    } else {
        2
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Sensor {
    /// Terra + Aqua (з 04.07.2002).
    TerraAqua,
    /// Лише Terra (01.11.2000 – 03.07.2002, коли Aqua ще не працювала).
    Terra,
}
impl Sensor {
    pub fn name(&self) -> &'static str {
        match self {
            Sensor::TerraAqua => "terra_aqua",
            Sensor::Terra => "terra_only",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hyper {
    /// Сила стягування до батьківської страти, у «псевдоточках MODIS».
    pub lambda: f64,
    /// Вага сусідніх місяців при оцінці страти (0 = без згладжування).
    pub smooth: f64,
    /// Множник βₛ за розміром події.
    pub bins: bool,
    /// Чи брати MODIS low-confidence (conf < 30) як предиктор.
    pub modis_low: bool,
    /// Просторова глибина листа: 1 = без простору, 2 = 4°, 3 = 2°, 4 = 1°, 5 = 0,5°.
    pub depth: usize,
    /// Спосіб FRP (див. [`FRP_MODE_NAMES`]).
    pub frp_mode: u8,
    /// Множник βₚ за патерном Terra/Aqua.
    pub pattern: bool,
    /// Множник βᵢ за інтенсивністю (FRP на детекцію MODIS).
    pub intensity: bool,
    /// Рівень земного покриву в ієрархії страт (потрібна колонка `lc`).
    pub lc: bool,
    /// Розподіл пропущеної маси по клітинко-добах з MODIS: 0 = пропорційно до детектованої
    /// частини; 1 = ∝ m; 2 = ∝ √m; 3 = порівну на кожну клітинко-добу («присутність»).
    pub miss_alloc: u8,
    /// Модель пропущеної маси (множники γ, EM; лише якщо `miss_alloc ≠ 0`):
    /// 0 = вимкнено; 1 = одиниця EM «страта × доба»; 2 = одиниця «страта × рік».
    pub miss_model: u8,
}

impl Default for Hyper {
    /// Стартова точка підбору гіперпараметрів.
    fn default() -> Self {
        Hyper {
            lambda: 30.0,
            smooth: 0.25,
            bins: true,
            modis_low: true,
            depth: 4,
            frp_mode: 1,
            pattern: false,
            intensity: false,
            lc: false,
            miss_alloc: 0,
            miss_model: 0,
        }
    }
}

impl Hyper {
    /// Параметри v3, вибрані крос-валідацією на даних України 2012–2019 (Terra+Aqua).
    pub fn v3_terra_aqua() -> Hyper {
        Hyper {
            lambda: 100.0,
            smooth: 0.5,
            bins: true,
            modis_low: true,
            depth: 4,
            frp_mode: 4,
            pattern: false,
            intensity: true,
            lc: false,
            miss_alloc: 3,
            miss_model: 2,
        }
    }
    /// Параметри v3, вибрані крос-валідацією на даних України 2012–2019 (лише Terra).
    pub fn v3_terra_only() -> Hyper {
        Hyper {
            lambda: 100.0,
            smooth: 0.25,
            bins: true,
            modis_low: true,
            depth: 3,
            frp_mode: 4,
            pattern: false,
            intensity: true,
            lc: true,
            miss_alloc: 3,
            miss_model: 2,
        }
    }
    pub fn v3(sensor: Sensor) -> Hyper {
        match sensor {
            Sensor::TerraAqua => Hyper::v3_terra_aqua(),
            Sensor::Terra => Hyper::v3_terra_only(),
        }
    }
}

/// Рівень ієрархії, що відповідає листовим стратам для цих гіперпараметрів.
#[inline]
pub fn leaf_level(h: &Hyper) -> usize {
    if h.lc {
        h.depth + 1
    } else {
        h.depth
    }
}

/// (кількість MODIS, FRP MODIS, патерн: 0 = лише Aqua / 1 = лише Terra / 2 = обидва)
#[inline]
pub fn predictor(a: &Agg, s: Sensor, low: bool) -> (u32, f64, usize) {
    let (mut am, mut af, mut tm, mut tf) = (a.a_hi, a.fa_hi, a.t_hi, a.ft_hi);
    if low {
        am += a.a_lo;
        af += a.fa_lo;
        tm += a.t_lo;
        tf += a.ft_lo;
    }
    match s {
        Sensor::TerraAqua => {
            let pat = if am > 0 && tm > 0 {
                2
            } else if tm > 0 {
                1
            } else {
                0
            };
            (am + tm, af + tf, pat)
        }
        Sensor::Terra => (tm, tf, 1),
    }
}

// ---------------------------------------------------------------- ієрархія страт

struct LevelSpec {
    month: bool,
    lc: bool,
    block: u32,
}

/// Ієрархія страт. Індекс страти на рівні L: `(region * 2 + dn) * nm + month`,
/// де region — щільний id пари (просторовий блок, група земного покриву).
pub struct Levels {
    pub n: Vec<usize>,
    pub idx: Vec<Vec<u32>>,
    pub parent: Vec<Vec<u32>>,
    pub nm: Vec<usize>,
    /// (bx, by, lc) для кожного регіону рівня; u32::MAX / 255 — «немає».
    pub regions: Vec<Vec<(u32, u32, u8)>>,
    pub block: Vec<u32>,
}

impl Levels {
    pub fn build(ds: &Dataset, month: &[u8], use_lc: bool) -> Levels {
        let mut specs = vec![LevelSpec { month: false, lc: false, block: 0 }, LevelSpec { month: true, lc: false, block: 0 }];
        if use_lc {
            specs.push(LevelSpec { month: true, lc: true, block: 0 });
        }
        for block in [40, 20, 10, 5] {
            specs.push(LevelSpec { month: true, lc: use_lc, block });
        }
        let nl = specs.len();
        let nrec = ds.cell.len();
        let nm: Vec<usize> = specs.iter().map(|s| if s.month { 12 } else { 1 }).collect();
        let mut idx: Vec<Vec<u32>> = (0..nl).map(|_| Vec::with_capacity(nrec)).collect();
        let mut regions: Vec<Vec<(u32, u32, u8)>> = vec![Vec::new(); nl];
        let mut maps: Vec<FxMap<(u32, u32, u8), u32>> = (0..nl).map(|_| FxMap::default()).collect();
        for i in 0..nrec {
            let c = ds.cell[i];
            let (cx, cy) = (c % 3600, c / 3600);
            let lcg = ds.lc_of(c);
            let dn = ds.dn[i] as usize;
            let mo = month[i] as usize;
            for (l, sp) in specs.iter().enumerate() {
                let region = if !sp.lc && sp.block == 0 {
                    0
                } else {
                    let key = (
                        if sp.block > 0 { cx / sp.block } else { u32::MAX },
                        if sp.block > 0 { cy / sp.block } else { u32::MAX },
                        if sp.lc { lcg } else { 255 },
                    );
                    let next = maps[l].len() as u32;
                    let regs = &mut regions[l];
                    *maps[l].entry(key).or_insert_with(|| {
                        regs.push(key);
                        next
                    })
                };
                let m = if nm[l] == 1 { 0 } else { mo };
                idx[l].push(((region as usize * 2 + dn) * nm[l] + m) as u32);
            }
        }
        let n: Vec<usize> = (0..nl)
            .map(|l| {
                let r = if !specs[l].lc && specs[l].block == 0 { 1 } else { regions[l].len().max(1) };
                r * 2 * nm[l]
            })
            .collect();
        let mut parent: Vec<Vec<u32>> = vec![Vec::new(); nl];
        for l in 1..nl {
            parent[l] = vec![u32::MAX; n[l]];
            for i in 0..nrec {
                parent[l][idx[l][i] as usize] = idx[l - 1][i];
            }
        }
        Levels { n, idx, parent, nm, regions, block: specs.iter().map(|s| s.block).collect() }
    }
}

// ---------------------------------------------------------------- контекст

/// Дані + похідні масиви, спільні для всіх підгонок.
pub struct Ctx<'a> {
    pub ds: &'a Dataset,
    pub year: Vec<i32>,
    pub month: Vec<u8>,
    pub y0: i32,
    pub ny: usize,
    /// Записи періоду перекриття (обидва сенсори працюють, без днів збою) — лише на них вчиться модель.
    pub overlap: Vec<bool>,
    pub ov_by_year: Vec<Vec<u32>>,
    /// Кількість інших клітинко-діб з MODIS у вікні 3×3 клітинки × ±1 доба: Terra+Aqua і лише Terra.
    pub nb_ta: Vec<u16>,
    pub nb_t: Vec<u16>,
    /// Ієрархія без земного покриву (також для метрик і калібрування).
    pub lv: Levels,
    /// Ієрархія з земним покривом (якщо у даних є `lc`).
    pub lv_lc: Option<Levels>,
}

impl<'a> Ctx<'a> {
    pub fn new(ds: &'a Dataset) -> Self {
        let n = ds.cell.len();
        let mut year = Vec::with_capacity(n);
        let mut month = Vec::with_capacity(n);
        for &d in &ds.day {
            let (y, m, _) = date::civil_from_days(d);
            year.push(y);
            month.push((m - 1) as u8);
        }
        let y0 = *year.iter().min().unwrap();
        let y1 = *year.iter().max().unwrap();
        let ny = (y1 - y0 + 1) as usize;
        let lv = Levels::build(ds, &month, false);
        let lv_lc = ds.has_lc.then(|| Levels::build(ds, &month, true));
        let (nb_ta, nb_t) = neighbourhood(ds);
        Ctx { ds, year, month, y0, ny, overlap: vec![false; n], ov_by_year: vec![Vec::new(); ny], nb_ta, nb_t, lv, lv_lc }
    }

    /// Сусідство запису i: модель «лише Terra» рахує його лише з Terra (Aqua тоді ще не було).
    #[inline]
    pub fn nb(&self, sensor: Sensor, i: usize) -> u16 {
        match sensor {
            Sensor::TerraAqua => self.nb_ta[i],
            Sensor::Terra => self.nb_t[i],
        }
    }

    #[inline]
    pub fn levels(&self, lc: bool) -> &Levels {
        if lc {
            self.lv_lc.as_ref().unwrap_or(&self.lv)
        } else {
            &self.lv
        }
    }

    pub fn set_overlap(&mut self, f: impl Fn(i32) -> bool) {
        let ds = self.ds;
        for v in &mut self.ov_by_year {
            v.clear();
        }
        for i in 0..ds.cell.len() {
            let ok = f(ds.day[i]);
            self.overlap[i] = ok;
            if ok {
                self.ov_by_year[(self.year[i] - self.y0) as usize].push(i as u32);
            }
        }
    }

    pub fn yi(&self, y: i32) -> Option<usize> {
        if y < self.y0 {
            return None;
        }
        let k = (y - self.y0) as usize;
        (k < self.ny).then_some(k)
    }

    pub fn year_weights(&self, years: &[i32], exclude: Option<i32>) -> Vec<f64> {
        let mut w = vec![0.0; self.ny];
        for &y in years {
            if Some(y) == exclude {
                continue;
            }
            if let Some(k) = self.yi(y) {
                w[k] += 1.0;
            }
        }
        w
    }
}

/// Для кожного запису: скільки **інших** клітинко-діб з MODIS у вікні 3×3 клітинки × ±1 доба.
fn neighbourhood(ds: &Dataset) -> (Vec<u16>, Vec<u16>) {
    let n = ds.cell.len();
    let (mut minx, mut maxx, mut miny, mut maxy) = (u32::MAX, 0u32, u32::MAX, 0u32);
    for &c in &ds.cell {
        minx = minx.min(c % 3600);
        maxx = maxx.max(c % 3600);
        miny = miny.min(c / 3600);
        maxy = maxy.max(c / 3600);
    }
    let pad = 1u32;
    let (x0, y0) = (minx.saturating_sub(pad), miny.saturating_sub(pad));
    let (w, h) = ((maxx + pad - x0 + 1) as usize, (maxy + pad - y0 + 1) as usize);
    let gi = |c: u32| -> usize { ((c / 3600 - y0) as usize) * w + (c % 3600 - x0) as usize };
    // активні дні по клітинках: [0] Terra+Aqua, [1] лише Terra
    let mut days: [Vec<Vec<i32>>; 2] = [vec![Vec::new(); w * h], vec![Vec::new(); w * h]];
    for i in 0..n {
        let a = &ds.agg[i];
        let g = gi(ds.cell[i]);
        if a.modis_raw() > 0 {
            days[0][g].push(ds.day[i]);
        }
        if a.t_hi + a.t_lo > 0 {
            days[1][g].push(ds.day[i]);
        }
    }
    for src in days.iter_mut() {
        for v in src.iter_mut() {
            v.sort_unstable();
            v.dedup();
        }
    }
    let count = |dsrc: &Vec<Vec<i32>>, i: usize| -> u16 {
        let (c, d) = (ds.cell[i], ds.day[i]);
        let (cx, cy) = ((c % 3600 - x0) as i32, (c / 3600 - y0) as i32);
        let mut cnt = 0u32;
        for dy in -1..=1 {
            for dx in -1..=1 {
                let (nx, ny) = (cx + dx, cy + dy);
                if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                    continue;
                }
                let v = &dsrc[ny as usize * w + nx as usize];
                let lo = v.partition_point(|&x| x < d - 1);
                let hi = v.partition_point(|&x| x <= d + 1);
                cnt += (hi - lo) as u32;
            }
        }
        let own = dsrc[gi(c)].binary_search(&d).is_ok() as u32;
        (cnt - own).min(u16::MAX as u32) as u16
    };
    let nb_ta = (0..n).map(|i| count(&days[0], i)).collect();
    let nb_t = (0..n).map(|i| count(&days[1], i)).collect();
    (nb_ta, nb_t)
}

// ---------------------------------------------------------------- прогноз

pub trait Predictor: Sync {
    /// (гармонізована кількість, гармонізований FRP) для запису i.
    fn predict(&self, ctx: &Ctx, i: usize) -> (f64, f64);
}

/// Навчена модель однієї підмоделі (Terra+Aqua або лише Terra).
#[derive(Clone)]
pub struct Fitted {
    pub sensor: Sensor,
    pub hyper: Hyper,
    /// Детектована частина: k по листових стратах, множники β, межі терцилів інтенсивності (МВт).
    pub k: Vec<f64>,
    pub beta_s: [[f64; NB]; 2],
    pub beta_p: [[f64; NP]; 2],
    pub beta_i: [[f64; NI]; 2],
    pub thr: [f64; 2],
    /// Пропущена маса: q (для `miss_alloc = 0`), g і множники γ.
    pub q: Vec<f64>,
    pub g: Vec<f64>,
    pub gm_i: [[f64; NI]; 2],
    pub gm_n: [[f64; NNB]; 2],
    pub gm_s: [[f64; NB]; 2],
    /// FRP (різні способи, див. [`FRP_MODE_NAMES`]).
    pub rho: Vec<f64>,
    pub qf: Vec<f64>,
    pub phi: Vec<f64>,
    pub rho_tot: Vec<f64>,
    pub rho_all: f64,
    pub phi_det: Vec<f64>,
    pub phi_miss: Vec<f64>,
    pub gf_i: [[f64; NI]; 2],
    /// Калібрування сум по стратах «день/ніч × місяць».
    pub cal: Vec<f64>,
    pub cal_frp: Vec<f64>,
    pub k_pooled: f64,
    pub q_pooled: f64,
    pub rho_pooled: f64,
    pub frp_per_count: f64,
}

impl Predictor for Fitted {
    #[inline]
    fn predict(&self, ctx: &Ctx, i: usize) -> (f64, f64) {
        let a = &ctx.ds.agg[i];
        let h = &self.hyper;
        let (m, fm, pat) = predictor(a, self.sensor, h.modis_low);
        if m == 0 {
            return (0.0, 0.0);
        }
        let l = ctx.levels(h.lc).idx[leaf_level(h)][i] as usize;
        let l1 = ctx.lv.idx[1][i] as usize;
        let dn = ctx.ds.dn[i] as usize;
        let ib = if fm > 0.0 { ibin(fm / m as f64, &self.thr) } else { 1 };
        let sb = bin_of(m);
        let det = m as f64 * self.k[l] * (self.beta_s[dn][sb] * self.beta_p[dn][pat] * self.beta_i[dn][ib]);
        let miss = if h.miss_alloc == 0 {
            det * self.q[l]
        } else {
            let nbb = nb_bin(ctx.nb(self.sensor, i));
            alloc_w(h.miss_alloc, m as f64) * self.g[l] * self.gm_i[dn][ib] * self.gm_n[dn][nbb] * self.gm_s[dn][sb]
        };
        let c = det + miss;
        let f = match h.frp_mode {
            4 => det * self.phi_det[l] * self.gf_i[dn][ib] + miss * self.phi_miss[l],
            0 if fm > 0.0 => fm * self.rho[l] * (1.0 + self.qf[l]),
            2 if fm > 0.0 => fm * self.rho_all,
            3 if fm > 0.0 => fm * self.rho_tot[l],
            1 => c * self.phi[l],
            _ => c * self.frp_per_count,
        };
        (c * self.cal[l1], f * self.cal_frp[l1])
    }
}

/// Базова лінія для порівняння: одна константа k = ΣVIIRS / ΣMODIS.
pub struct ConstModel {
    pub sensor: Sensor,
    pub k: f64,
    pub kf: f64,
}
impl ConstModel {
    pub fn fit(ctx: &Ctx, sensor: Sensor, yw: &[f64]) -> Self {
        let (mut sm, mut sv, mut sfm, mut sfv) = (0.0, 0.0, 0.0, 0.0);
        for (y, &w) in yw.iter().enumerate() {
            if w <= 0.0 {
                continue;
            }
            for &i in &ctx.ov_by_year[y] {
                let a = &ctx.ds.agg[i as usize];
                let (m, fm, _) = predictor(a, sensor, true);
                sm += w * m as f64;
                sfm += w * fm;
                sv += w * a.v as f64;
                sfv += w * a.fv;
            }
        }
        ConstModel { sensor, k: if sm > 0.0 { sv / sm } else { 1.0 }, kf: if sfm > 0.0 { sfv / sfm } else { 1.0 } }
    }
}
impl Predictor for ConstModel {
    fn predict(&self, ctx: &Ctx, i: usize) -> (f64, f64) {
        let (m, fm, _) = predictor(&ctx.ds.agg[i], self.sensor, true);
        (m as f64 * self.k, fm * self.kf)
    }
}

// ---------------------------------------------------------------- навчання

/// Ієрархічне стягування відношення num/den згори донизу; повертає відношення на кожному рівні.
pub fn hier(lv: &Levels, depth: usize, num_leaf: &[f64], den_leaf: &[f64], lam: f64, smooth: f64) -> Vec<Vec<f64>> {
    let mut num: Vec<Vec<f64>> = (0..=depth).map(|l| vec![0.0; lv.n[l]]).collect();
    let mut den: Vec<Vec<f64>> = (0..=depth).map(|l| vec![0.0; lv.n[l]]).collect();
    num[depth].copy_from_slice(num_leaf);
    den[depth].copy_from_slice(den_leaf);
    for l in (1..=depth).rev() {
        let (nlo, nhi) = num.split_at_mut(l);
        let (dlo, dhi) = den.split_at_mut(l);
        for s in 0..lv.n[l] {
            let p = lv.parent[l][s];
            if p == u32::MAX {
                continue;
            }
            nlo[l - 1][p as usize] += nhi[0][s];
            dlo[l - 1][p as usize] += dhi[0][s];
        }
    }
    let tn: f64 = num[0].iter().sum();
    let td: f64 = den[0].iter().sum();
    let pooled = if td > 0.0 { tn / td } else { 0.0 };
    let ratio = |a: f64, b: f64, prior: f64| if b + lam > 0.0 { (a + lam * prior) / (b + lam) } else { prior };
    let mut r: Vec<Vec<f64>> = Vec::with_capacity(depth + 1);
    r.push(num[0].iter().zip(&den[0]).map(|(&a, &b)| ratio(a, b, pooled)).collect());
    for l in 1..=depth {
        let nm = lv.nm[l];
        let mut rl = vec![pooled; lv.n[l]];
        for (s, out) in rl.iter_mut().enumerate() {
            let p = lv.parent[l][s];
            if p == u32::MAX {
                continue;
            }
            let (mut a, mut b) = (num[l][s], den[l][s]);
            if smooth > 0.0 && nm == 12 {
                let base = s - s % 12;
                let mo = s % 12;
                let prev = base + (mo + 11) % 12;
                let next = base + (mo + 1) % 12;
                a += smooth * (num[l][prev] + num[l][next]);
                b += smooth * (den[l][prev] + den[l][next]);
            }
            *out = ratio(a, b, r[l - 1][p as usize]);
        }
        r.push(rl);
    }
    r
}

/// Навчальний запис (копія потрібних полів у компактному вигляді).
struct Tr {
    i: u32,
    day: i32,
    year: i32,
    l: u32,
    dn: u8,
    sb: u8,
    pb: u8,
    ib: u8,
    nbb: u8,
    w: f64,
    m: f64,
    fm: f64,
    v: f64,
    fv: f64,
}

/// Зважені терцилі (значення, вага).
fn tertiles(vals: &mut [(f64, f64)]) -> [f64; 2] {
    if vals.is_empty() {
        return [f64::INFINITY, f64::INFINITY];
    }
    vals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let tot: f64 = vals.iter().map(|x| x.1).sum();
    let mut out = [vals[vals.len() - 1].0; 2];
    let mut acc = 0.0;
    let mut k = 0;
    for &(v, w) in vals.iter() {
        acc += w;
        while k < 2 && acc >= tot * (k as f64 + 1.0) / 3.0 {
            out[k] = v;
            k += 1;
        }
        if k == 2 {
            break;
        }
    }
    out
}

/// Оновлення множників з нормуванням до зваженого середнього 1:
/// `new[b] = (num[b] + lam·prior) / (den[b] + lam)`, потім ділення на зважене середнє.
fn normalized(num: &[f64], den: &[f64], lam: f64, prior: f64) -> Vec<f64> {
    let new: Vec<f64> = num.iter().zip(den).map(|(&a, &b)| (a + lam * prior) / (b + lam)).collect();
    let tot: f64 = den.iter().sum();
    let totu: f64 = den.iter().zip(&new).map(|(&d, &x)| d * x).sum();
    let s = if tot > 0.0 && totu > 0.0 { totu / tot } else { 1.0 };
    new.into_iter().map(|x| x / s).collect()
}

/// Навчання на записах періоду перекриття з вагами років `yw`
/// (0 — рік не бере участі; > 1 — повтор року в bootstrap).
pub fn fit(ctx: &Ctx, sensor: Sensor, h: Hyper, yw: &[f64]) -> Fitted {
    let ds = ctx.ds;
    let h = Hyper { depth: h.depth.clamp(1, MAX_DEPTH), lc: h.lc && ctx.lv_lc.is_some(), pattern: h.pattern && sensor == Sensor::TerraAqua, ..h };
    let lv = ctx.levels(h.lc);
    let d = leaf_level(&h);
    let nl = lv.n[d];
    let leaf = &lv.idx[d];

    // ---- навчальні записи
    let mut tr: Vec<Tr> = Vec::new();
    let mut ints: Vec<(f64, f64)> = Vec::new();
    for (y, &w) in yw.iter().enumerate() {
        if w <= 0.0 {
            continue;
        }
        for &i in &ctx.ov_by_year[y] {
            let i = i as usize;
            let a = &ds.agg[i];
            let (m, fm, pat) = predictor(a, sensor, h.modis_low);
            if m > 0 && fm > 0.0 {
                ints.push((fm / m as f64, w * m as f64));
            }
            tr.push(Tr {
                i: i as u32,
                day: ds.day[i],
                year: ctx.year[i],
                l: leaf[i],
                dn: ds.dn[i],
                sb: bin_of(m) as u8,
                pb: pat as u8,
                ib: 1,
                nbb: nb_bin(ctx.nb(sensor, i)) as u8,
                w,
                m: m as f64,
                fm,
                v: a.v as f64,
                fv: a.fv,
            });
        }
    }
    let thr = if h.intensity { tertiles(&mut ints) } else { [f64::INFINITY, f64::INFINITY] };
    for t in &mut tr {
        t.ib = if t.m > 0.0 && t.fm > 0.0 { ibin(t.fm / t.m, &thr) as u8 } else { 1 };
    }

    let (mut sm, mut sv, mut sfm, mut sfv, mut sv_all, mut sfv_all) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for t in &tr {
        if t.m > 0.0 {
            sm += t.w * t.m;
            sv += t.w * t.v;
            sfm += t.w * t.fm;
            sfv += t.w * t.fv;
        }
        sv_all += t.w * t.v;
        sfv_all += t.w * t.fv;
    }
    let k_pool = if sm > 0.0 { sv / sm } else { 1.0 };
    let frp_m_pool = if sm > 0.0 { (sfm / sm).max(1e-3) } else { 10.0 };
    let rho_pool = if sfm > 0.0 { sfv / sfm } else { 1.0 };
    let frp_per_count = if sv_all > 0.0 { sfv_all / sv_all } else { 5.0 };

    // ---- детектована частина: k і множники β (IPF)
    let mut bs = [[1.0f64; NB]; 2];
    let mut bp = [[1.0f64; NP]; 2];
    let mut bi = [[1.0f64; NI]; 2];
    let fac = |t: &Tr, bs: &[[f64; NB]; 2], bp: &[[f64; NP]; 2], bi: &[[f64; NI]; 2]| {
        bs[t.dn as usize][t.sb as usize] * bp[t.dn as usize][t.pb as usize] * bi[t.dn as usize][t.ib as usize]
    };
    let any = h.bins || h.pattern || h.intensity;
    let iters = if any { 10 } else { 1 };
    let lamb = h.lambda * k_pool;
    let mut k_leaf: Vec<f64> = Vec::new();
    for it in 0..iters {
        let mut num = vec![0.0; nl];
        let mut den = vec![0.0; nl];
        for t in &tr {
            if t.m > 0.0 {
                num[t.l as usize] += t.w * t.v;
                den[t.l as usize] += t.w * t.m * fac(t, &bs, &bp, &bi);
            }
        }
        k_leaf = hier(lv, d, &num, &den, h.lambda, h.smooth).pop().unwrap();
        if !any || it + 1 == iters {
            continue;
        }
        for which in 0..3 {
            let (on, nb) = match which {
                0 => (h.bins, NB),
                1 => (h.pattern, NP),
                _ => (h.intensity, NI),
            };
            if !on {
                continue;
            }
            let mut bn = [[0.0f64; NB]; 2];
            let mut bd = [[0.0f64; NB]; 2];
            for t in &tr {
                if t.m <= 0.0 {
                    continue;
                }
                let (dn, b, own) = match which {
                    0 => (t.dn as usize, t.sb as usize, bs[t.dn as usize][t.sb as usize]),
                    1 => (t.dn as usize, t.pb as usize, bp[t.dn as usize][t.pb as usize]),
                    _ => (t.dn as usize, t.ib as usize, bi[t.dn as usize][t.ib as usize]),
                };
                bn[dn][b] += t.w * t.v;
                bd[dn][b] += t.w * t.m * k_leaf[t.l as usize] * fac(t, &bs, &bp, &bi) / own;
            }
            for dn in 0..2 {
                let upd = normalized(&bn[dn][..nb], &bd[dn][..nb], lamb, 1.0);
                for (b, &v) in upd.iter().enumerate() {
                    match which {
                        0 => bs[dn][b] = v,
                        1 => bp[dn][b] = v,
                        _ => bi[dn][b] = v,
                    }
                }
            }
        }
    }

    // ---- FRP: відношення до FRP MODIS
    let mut fnum = vec![0.0; nl];
    let mut fden = vec![0.0; nl];
    for t in &tr {
        if t.m > 0.0 {
            fnum[t.l as usize] += t.w * t.fv;
            fden[t.l as usize] += t.w * t.fm;
        }
    }
    let rho = hier(lv, d, &fnum, &fden, h.lambda * frp_m_pool, h.smooth).pop().unwrap();

    // ---- пропущена маса: q (пропорційний розподіл)
    let mut vm = vec![0.0; nl];
    let mut p = vec![0.0; nl];
    let mut fvm = vec![0.0; nl];
    let mut pf = vec![0.0; nl];
    for t in &tr {
        let l = t.l as usize;
        if t.m > 0.0 {
            p[l] += t.w * t.m * k_leaf[l] * fac(t, &bs, &bp, &bi);
            pf[l] += t.w * t.fm * rho[l];
        } else {
            vm[l] += t.w * t.v;
            fvm[l] += t.w * t.fv;
        }
    }
    let q_pooled = {
        let (a, b): (f64, f64) = (vm.iter().sum(), p.iter().sum());
        if b > 0.0 {
            a / b
        } else {
            0.0
        }
    };
    let q = hier(lv, d, &vm, &p, h.lambda * k_pool, h.smooth).pop().unwrap();

    // ---- пропущена маса: g[s] і множники γ (EM)
    let mut gm_i = [[1.0f64; NI]; 2];
    let mut gm_n = [[1.0f64; NNB]; 2];
    let mut gm_s = [[1.0f64; NB]; 2];
    let gam = |t: &Tr, gi: &[[f64; NI]; 2], gn: &[[f64; NNB]; 2], gs: &[[f64; NB]; 2]| {
        let dn = t.dn as usize;
        gi[dn][t.ib as usize] * gn[dn][t.nbb as usize] * gs[dn][t.sb as usize]
    };
    let mut g = vec![0.0; nl];
    if h.miss_alloc != 0 {
        let tot_alloc: f64 = tr.iter().filter(|t| t.m > 0.0).map(|t| t.w * alloc_w(h.miss_alloc, t.m)).sum();
        let tot_m: f64 = tr.iter().filter(|t| t.m > 0.0).map(|t| t.w * t.m).sum();
        let scale = if tot_m > 0.0 { tot_alloc / tot_m } else { 1.0 };
        let em_on = h.miss_model != 0;
        let by_year = h.miss_model == 2;
        // одиниця EM: листова страта × рік (або × доба); щільні індекси, щоб EM працював на масивах
        let unit = |t: &Tr| ((t.l as u64) << 32) | if by_year { t.year as u32 as u64 } else { t.day as u32 as u64 };
        let mut uid: Vec<u32> = vec![u32::MAX; tr.len()];
        let mut target: Vec<f64> = Vec::new();
        if em_on {
            let mut map: FxMap<u64, u32> = FxMap::default();
            for (j, t) in tr.iter().enumerate() {
                if t.m > 0.0 || t.v > 0.0 {
                    let n = map.len() as u32;
                    let id = *map.entry(unit(t)).or_insert(n);
                    if id as usize >= target.len() {
                        target.push(0.0);
                    }
                    uid[j] = id;
                    if t.m <= 0.0 {
                        target[id as usize] += t.w * t.v;
                    }
                }
            }
        }
        let iters = if em_on { 10 } else { 1 };
        for it in 0..iters {
            let mut aw = vec![0.0; nl];
            for t in &tr {
                if t.m > 0.0 {
                    aw[t.l as usize] += t.w * alloc_w(h.miss_alloc, t.m) * gam(t, &gm_i, &gm_n, &gm_s);
                }
            }
            g = hier(lv, d, &vm, &aw, h.lambda * scale, h.smooth).pop().unwrap();
            if !em_on || it + 1 == iters {
                continue;
            }
            let pred = |t: &Tr, gi: &[[f64; NI]; 2], gn: &[[f64; NNB]; 2], gs: &[[f64; NB]; 2]| {
                t.w * alloc_w(h.miss_alloc, t.m) * g[t.l as usize] * gam(t, gi, gn, gs)
            };
            for which in 0..3 {
                // E-крок: ціль одиниці ділиться між її активними записами пропорційно до прогнозу
                let mut psum = vec![0.0f64; target.len()];
                for (j, t) in tr.iter().enumerate() {
                    if t.m > 0.0 {
                        psum[uid[j] as usize] += pred(t, &gm_i, &gm_n, &gm_s);
                    }
                }
                let nb = [NI, NNB, NB][which];
                let mut num = [[0.0f64; NB]; 2];
                let mut den = [[0.0f64; NB]; 2];
                for (j, t) in tr.iter().enumerate() {
                    if t.m <= 0.0 {
                        continue;
                    }
                    let p = pred(t, &gm_i, &gm_n, &gm_s);
                    let u = uid[j] as usize;
                    let attributed = if psum[u] > 0.0 { target[u] * p / psum[u] } else { 0.0 };
                    let b = [t.ib as usize, t.nbb as usize, t.sb as usize][which];
                    num[t.dn as usize][b] += attributed;
                    den[t.dn as usize][b] += p;
                }
                // M-крок: мультиплікативне оновлення зі стягуванням і нормуванням до середнього 1
                for dn in 0..2 {
                    let r_all = {
                        let (a, b): (f64, f64) = ((0..nb).map(|b| num[dn][b]).sum(), (0..nb).map(|b| den[dn][b]).sum());
                        if b > 0.0 {
                            a / b
                        } else {
                            1.0
                        }
                    };
                    let upd = normalized(&num[dn][..nb], &den[dn][..nb], lamb.max(1e-9), r_all);
                    for (b, &v) in upd.iter().enumerate() {
                        match which {
                            0 => gm_i[dn][b] *= v,
                            1 => gm_n[dn][b] *= v,
                            _ => gm_s[dn][b] *= v,
                        }
                    }
                }
            }
        }
    }

    let qf = hier(lv, d, &fvm, &pf, h.lambda * frp_m_pool * rho_pool, h.smooth).pop().unwrap();

    // ---- FRP: φ і повне відношення
    let mut pn = vec![0.0; nl];
    let mut pd = vec![0.0; nl];
    let mut rd = vec![0.0; nl];
    let mut sfm_all = 0.0;
    for t in &tr {
        pn[t.l as usize] += t.w * t.fv;
        pd[t.l as usize] += t.w * t.v;
        rd[t.l as usize] += t.w * t.fm;
        sfm_all += t.w * t.fm;
    }
    let phi = hier(lv, d, &pn, &pd, h.lambda * k_pool, h.smooth).pop().unwrap();
    let rho_tot = hier(lv, d, &pn, &rd, h.lambda * frp_m_pool, h.smooth).pop().unwrap();
    let rho_all = if sfm_all > 0.0 { sfv_all / sfm_all } else { 1.0 };

    // ---- FRP: окремо детектована і пропущена частини + множник інтенсивності
    let (mut dnum, mut dden, mut mnum, mut mden) = (vec![0.0; nl], vec![0.0; nl], vec![0.0; nl], vec![0.0; nl]);
    for t in &tr {
        let l = t.l as usize;
        if t.m > 0.0 {
            dnum[l] += t.w * t.fv;
            dden[l] += t.w * t.v;
        } else {
            mnum[l] += t.w * t.fv;
            mden[l] += t.w * t.v;
        }
    }
    let phi_det = hier(lv, d, &dnum, &dden, h.lambda * k_pool, h.smooth).pop().unwrap();
    let phi_miss = hier(lv, d, &mnum, &mden, h.lambda * k_pool, h.smooth).pop().unwrap();
    let mut gf_i = [[1.0f64; NI]; 2];
    if h.intensity {
        let mut num = [[0.0f64; NI]; 2];
        let mut den = [[0.0f64; NI]; 2];
        for t in &tr {
            if t.m > 0.0 {
                let det = t.w * t.m * k_leaf[t.l as usize] * fac(t, &bs, &bp, &bi);
                num[t.dn as usize][t.ib as usize] += t.w * t.fv;
                den[t.dn as usize][t.ib as usize] += det * phi_det[t.l as usize];
            }
        }
        for dn in 0..2 {
            let upd = normalized(&num[dn], &den[dn], 20.0 * frp_per_count, 1.0);
            gf_i[dn].copy_from_slice(&upd);
        }
    }

    let n1 = ctx.lv.n[1];
    let mut f = Fitted {
        sensor,
        hyper: h,
        k: k_leaf,
        beta_s: bs,
        beta_p: bp,
        beta_i: bi,
        thr,
        q,
        g,
        gm_i,
        gm_n,
        gm_s,
        rho,
        qf,
        phi,
        rho_tot,
        rho_all,
        phi_det,
        phi_miss,
        gf_i,
        cal: vec![1.0; n1],
        cal_frp: vec![1.0; n1],
        k_pooled: k_pool,
        q_pooled,
        rho_pooled: rho_pool,
        frp_per_count,
    };

    // ---- калібрування сум: на навчальних роках сума прогнозу = сума S-NPP у кожній страті «день/ніч × місяць»
    let (mut v1, mut p1, mut fv1, mut fp1) = (vec![0.0; n1], vec![0.0; n1], vec![0.0; n1], vec![0.0; n1]);
    for t in &tr {
        let (c, fr) = f.predict(ctx, t.i as usize);
        let j = ctx.lv.idx[1][t.i as usize] as usize;
        v1[j] += t.w * t.v;
        fv1[j] += t.w * t.fv;
        p1[j] += t.w * c;
        fp1[j] += t.w * fr;
    }
    let lc = 20.0;
    let lf = 20.0 * frp_per_count;
    for j in 0..n1 {
        f.cal[j] = (v1[j] + lc) / (p1[j] + lc);
        f.cal_frp[j] = (fv1[j] + lf) / (fp1[j] + lf);
    }
    f
}
