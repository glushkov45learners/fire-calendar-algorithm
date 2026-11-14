//! Читання `hotspots.csv.gz` (контракт 1) і агрегація до записів «клітинка × доба × день/ніч».
//!
//! Фільтри:
//! - `type` 2 і 3 (заводи, факели, офшор) виключаються; їхні координати з SP-продукту
//!   формують «статичну маску» (мікроклітинки ~0,005°);
//! - точки без `type` (NRT, 2026 рік) перевіряються маскою 3×3 мікроклітинки: промислові
//!   джерела в NRT не мають позначки, тому без маски вони б потрапили в ряд;
//! - NOAA-20/21 рахуються лише в діагностиці (опорний ряд — тільки S-NPP);
//! - MODIS low/high confidence зберігаються окремо: модель сама вирішує (через CV),
//!   чи включати low-confidence MODIS як предиктор.

use crate::date;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::io::BufRead;

#[derive(Default, Clone, Copy)]
pub struct Fx(u64);
impl Fx {
    #[inline]
    fn add(&mut self, i: u64) {
        self.0 = (self.0.rotate_left(5) ^ i).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}
impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.add(b as u64);
        }
    }
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }
    fn write_i32(&mut self, i: i32) {
        self.add(i as u32 as u64);
    }
}
pub type FxMap<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;
pub type FxSet<K> = HashSet<K, BuildHasherDefault<Fx>>;

pub const TERRA: usize = 0;
pub const AQUA: usize = 1;
pub const SNPP: usize = 2;
pub const NOAA20: usize = 3;
pub const NOAA21: usize = 4;
pub const SAT_NAMES: [&str; 5] = ["Terra", "Aqua", "SNPP", "NOAA20", "NOAA21"];

/// Групи земного покриву (IGBP, колонка `lc`). Клітинці 0,1° присвоюється група,
/// найчастіша серед її точок.
pub const NLC: usize = 6;
pub const LC_GROUP_NAMES: [&str; NLC] =
    ["ліс (IGBP 1–5)", "савани/чагарники (6–9)", "луки (10)", "ріллі (12)", "мозаїка рілля/природне (14)", "інше (11, 13, 15–17)"];
pub fn lc_group(lc: u8) -> u8 {
    match lc {
        1..=5 => 0,
        6..=9 => 1,
        10 => 2,
        12 => 3,
        14 => 4,
        _ => 5,
    }
}

/// Сума по одній клітинці-добі-(день|ніч). `_lo` — точки MODIS з conf = `l`.
#[derive(Clone, Copy, Default, Debug)]
pub struct Agg {
    pub t_hi: u32,
    pub t_lo: u32,
    pub a_hi: u32,
    pub a_lo: u32,
    pub v: u32,
    pub ft_hi: f64,
    pub ft_lo: f64,
    pub fa_hi: f64,
    pub fa_lo: f64,
    pub fv: f64,
}
impl Agg {
    /// Сира кількість MODIS (Terra + Aqua, усі confidence) — колонка `modis_count`.
    pub fn modis_raw(&self) -> u32 {
        self.t_hi + self.t_lo + self.a_hi + self.a_lo
    }
}

#[derive(Default, Debug, Clone)]
pub struct LoadStats {
    pub lines: u64,
    pub bad_lines: u64,
    pub type_excluded: u64,
    pub nrt_untyped: u64,
    pub nrt_masked: u64,
    pub viirs_low_dropped: u64,
    pub static_mask_size: usize,
    /// Маска стабільних джерел за даними: мікроклітинки з детекціями в багатьох місяцях і роках.
    pub persistent_cells: usize,
    pub persistent_masked: u64,
    pub by_sat: [u64; 5],
}

/// Записи, відсортовані за (cell, day, dn). Struct-of-arrays заради швидкості.
pub struct Dataset {
    pub cell: Vec<u32>,
    pub day: Vec<i32>,
    pub dn: Vec<u8>,
    pub agg: Vec<Agg>,
    pub first: [Option<i32>; 5],
    pub last: [Option<i32>; 5],
    /// Кількість точок по роках і супутниках (після фільтра type/маски, усі confidence).
    pub raw_year: BTreeMap<i32, [u64; 5]>,
    pub stats: LoadStats,
    /// Група земного покриву для клітинки (якщо у вхідному файлі є колонка `lc`).
    pub cell_lc: FxMap<u32, u8>,
    pub has_lc: bool,
}

impl Dataset {
    #[inline]
    pub fn lc_of(&self, cell: u32) -> u8 {
        *self.cell_lc.get(&cell).unwrap_or(&5)
    }
}

struct Pt {
    day: i32,
    sat: u8,
    low: bool,
    frp: f64,
    dn: u8,
    cell: u32,
    micro: (i32, i32),
    lc: u8,
}

fn micro(lat: f64, lon: f64) -> (i32, i32) {
    ((lat * 200.0).round() as i32, (lon * 200.0).round() as i32)
}

struct Builder {
    map: FxMap<u64, Agg>,
    first: [Option<i32>; 5],
    last: [Option<i32>; 5],
    raw_year: BTreeMap<i32, [u64; 5]>,
    viirs_drop_low: bool,
    stats: LoadStats,
    lc_counts: FxMap<u32, [u32; 18]>,
}

impl Builder {
    fn accept(&mut self, p: &Pt) {
        let s = p.sat as usize;
        let (y, _, _) = date::civil_from_days(p.day);
        self.raw_year.entry(y).or_insert([0; 5])[s] += 1;
        self.stats.by_sat[s] += 1;
        self.first[s] = Some(self.first[s].map_or(p.day, |f| f.min(p.day)));
        self.last[s] = Some(self.last[s].map_or(p.day, |f| f.max(p.day)));
        if s >= NOAA20 {
            return;
        }
        if s == SNPP && self.viirs_drop_low && p.low {
            self.stats.viirs_low_dropped += 1;
            return;
        }
        if (p.lc as usize) < 18 {
            self.lc_counts.entry(p.cell).or_insert([0; 18])[p.lc as usize] += 1;
        }
        let key = ((p.cell as u64) << 24) | ((p.day as u64) << 1) | p.dn as u64;
        let a = self.map.entry(key).or_default();
        match (s, p.low) {
            (TERRA, false) => {
                a.t_hi += 1;
                a.ft_hi += p.frp
            }
            (TERRA, true) => {
                a.t_lo += 1;
                a.ft_lo += p.frp
            }
            (AQUA, false) => {
                a.a_hi += 1;
                a.fa_hi += p.frp
            }
            (AQUA, true) => {
                a.a_lo += 1;
                a.fa_lo += p.frp
            }
            _ => {
                a.v += 1;
                a.fv += p.frp
            }
        }
    }
}

/// `persist`: Some((мін. різних місяців, мін. різних років)) — маскувати мікроклітинки ~0,005°,
/// де детекції трапляються стабільно (заводи, факели, які FIRMS не позначив як type 2).
pub fn load<R: BufRead>(mut r: R, viirs_drop_low: bool, persist: Option<(usize, usize)>) -> Result<Dataset, String> {
    let mut line = String::new();
    if r.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
        return Err("порожній вхідний файл".into());
    }
    let header: Vec<String> = line
        .trim_end()
        .trim_start_matches('\u{feff}')
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let col = |n: &str| header.iter().position(|h| h == n).ok_or_else(|| format!("у заголовку немає колонки `{n}`"));
    let i_date = col("acq_date")?;
    let i_lat = col("lat")?;
    let i_lon = col("lon")?;
    let i_sat = col("satellite")?;
    let i_conf = col("conf")?;
    let i_frp = col("frp")?;
    let i_dn = col("daynight")?;
    let i_type = col("type")?;
    let i_prod = col("product")?;
    let i_cell = col("cell_id")?;
    let i_lc = header.iter().position(|h| h == "lc");
    let need = [i_date, i_lat, i_lon, i_sat, i_conf, i_frp, i_dn, i_type, i_prod, i_cell, i_lc.unwrap_or(0)].into_iter().max().unwrap() + 1;
    if need > 32 {
        return Err("забагато колонок у заголовку".into());
    }

    let mut b = Builder {
        map: FxMap::default(),
        first: [None; 5],
        last: [None; 5],
        raw_year: BTreeMap::new(),
        viirs_drop_low,
        stats: LoadStats::default(),
        lc_counts: FxMap::default(),
    };
    let mut static_mask: FxSet<(i32, i32)> = FxSet::default();
    let mut pending: Vec<Pt> = Vec::new();
    let mut typed: Vec<Pt> = Vec::new();

    loop {
        line.clear();
        let n = r.read_line(&mut line).map_err(|e| format!("помилка читання: {e}"))?;
        if n == 0 {
            break;
        }
        let l = line.trim_end_matches(&['\n', '\r'][..]);
        if l.is_empty() {
            continue;
        }
        b.stats.lines += 1;
        let mut f: [&str; 32] = [""; 32];
        let mut nf = 0;
        for s in l.split(',') {
            if nf < 32 {
                f[nf] = s;
            }
            nf += 1;
        }
        if nf < need {
            b.stats.bad_lines += 1;
            continue;
        }
        let parsed = (|| -> Option<(Pt, Option<u8>, bool)> {
            let day = date::parse_date(f[i_date])?;
            let sat = match f[i_sat] {
                "Terra" => 0u8,
                "Aqua" => 1,
                "SNPP" => 2,
                "NOAA20" => 3,
                "NOAA21" => 4,
                _ => return None,
            };
            let lat: f64 = f[i_lat].parse().ok()?;
            let lon: f64 = f[i_lon].parse().ok()?;
            let cell: u32 = f[i_cell].parse().ok()?;
            let frp: f64 = f[i_frp].parse().unwrap_or(0.0);
            let frp = if frp.is_finite() && frp > 0.0 { frp } else { 0.0 };
            let low = f[i_conf] == "l";
            let dn = if f[i_dn] == "N" { 1 } else { 0 };
            let typ: Option<u8> = f[i_type].trim().parse().ok();
            let sp = f[i_prod] == "SP";
            if day < 0 || day >= (1 << 23) || cell >= (1 << 23) {
                return None;
            }
            let lc = match i_lc {
                Some(k) => f[k].trim().parse::<u8>().ok().filter(|&x| x < 18).unwrap_or(255),
                None => 255,
            };
            Some((Pt { day, sat, low, frp, dn, cell, micro: micro(lat, lon), lc }, typ, sp))
        })();
        let Some((p, typ, sp)) = parsed else {
            b.stats.bad_lines += 1;
            continue;
        };
        match typ {
            Some(2) | Some(3) => {
                b.stats.type_excluded += 1;
                if sp {
                    static_mask.insert(p.micro);
                }
            }
            Some(_) => typed.push(p),
            None => {
                b.stats.nrt_untyped += 1;
                pending.push(p);
            }
        }
    }

    b.stats.static_mask_size = static_mask.len();
    let near = |set: &FxSet<(i32, i32)>, (a, o): (i32, i32)| -> bool {
        for da in -1..=1 {
            for db in -1..=1 {
                if set.contains(&(a + da, o + db)) {
                    return true;
                }
            }
        }
        false
    };
    let mut cand: Vec<Pt> = typed;
    for p in pending {
        if near(&static_mask, p.micro) {
            b.stats.nrt_masked += 1;
        } else {
            cand.push(p);
        }
    }
    // --- маска стабільних джерел за даними ---
    let mut pmask: FxSet<(i32, i32)> = FxSet::default();
    if let Some((min_months, min_years)) = persist {
        let mut pairs: Vec<(u64, u32)> = cand
            .iter()
            .map(|p| {
                let (y, m, _) = date::civil_from_days(p.day);
                (((p.micro.0 as u32 as u64) << 32) | (p.micro.1 as u32 as u64), (y * 12 + m as i32 - 1) as u32)
            })
            .collect();
        pairs.sort_unstable();
        pairs.dedup();
        let mut i = 0;
        while i < pairs.len() {
            let key = pairs[i].0;
            let (mut months, mut years, mut last_y) = (0usize, 0usize, u32::MAX);
            while i < pairs.len() && pairs[i].0 == key {
                months += 1;
                let y = pairs[i].1 / 12;
                if y != last_y {
                    years += 1;
                    last_y = y;
                }
                i += 1;
            }
            if months >= min_months && years >= min_years {
                pmask.insert(((key >> 32) as u32 as i32, key as u32 as i32));
            }
        }
    }
    b.stats.persistent_cells = pmask.len();
    for p in &cand {
        if !pmask.is_empty() && near(&pmask, p.micro) {
            b.stats.persistent_masked += 1;
        } else {
            b.accept(p);
        }
    }

    let mut v: Vec<(u64, Agg)> = b.map.into_iter().collect();
    v.sort_unstable_by_key(|x| x.0);
    let n = v.len();
    let (mut cell, mut day, mut dn, mut agg) =
        (Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n));
    for (k, a) in v {
        cell.push((k >> 24) as u32);
        day.push(((k >> 1) & ((1 << 23) - 1)) as i32);
        dn.push((k & 1) as u8);
        agg.push(a);
    }
    if cell.is_empty() {
        return Err("після фільтрів не лишилось жодної точки".into());
    }
    let mut cell_lc: FxMap<u32, u8> = FxMap::default();
    for (c, cnt) in &b.lc_counts {
        let (best, _) = cnt.iter().enumerate().max_by_key(|(_, &v)| v).unwrap();
        cell_lc.insert(*c, lc_group(best as u8));
    }
    let has_lc = i_lc.is_some() && !cell_lc.is_empty();
    Ok(Dataset { cell, day, dn, agg, first: b.first, last: b.last, raw_year: b.raw_year, stats: b.stats, cell_lc, has_lc })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parse_and_mask() {
        let csv = "acq_date,acq_time,lat,lon,instrument,satellite,conf,conf_raw,frp,daynight,type,product,cell_id\n\
2015-01-01,0931,47.08959,37.612,VIIRS,SNPP,n,n,4.45,D,2,SP,4934176\n\
2015-01-01,0931,47.08959,37.612,VIIRS,SNPP,n,n,4.45,D,0,SP,4934176\n\
2015-01-01,0931,47.08,37.61,MODIS,Aqua,l,20,9.0,D,0,SP,4934176\n\
2026-09-01,0219,47.0896,37.6121,VIIRS,SNPP,n,n,8.36,N,,NRT,4934176\n\
2026-09-01,0219,45.25808,31.68034,MODIS,Aqua,n,53,8.36,N,,NRT,4869316\n\
bad,line\n";
        let ds = load(csv.as_bytes(), false, None).unwrap();
        assert_eq!(ds.stats.type_excluded, 1);
        assert_eq!(ds.stats.nrt_masked, 1); // NRT точка біля заводу
        assert_eq!(ds.stats.bad_lines, 1);
        assert_eq!(ds.cell.len(), 2); // (4934176, 2015-01-01, D) і (4869316, 2026-09-01, N)
        // сортування за cell_id: спершу 4869316, потім 4934176
        assert_eq!(ds.cell[0], 4869316);
        assert_eq!(ds.agg[0].a_hi, 1);
        assert_eq!(ds.agg[1].v, 1);
        assert_eq!(ds.agg[1].a_lo, 1);
    }
}
