//! `check` — перевірка результатів `harmonize`.
//!
//! 1. `cell_day.csv.gz` відповідає контракту 2: заголовок, унікальний ключ (cell_id, date),
//!    невід'ємні значення, жодного порожнього рядка, дата у форматі YYYY-MM-DD.
//! 2. `model_meta.json` містить обов'язкові поля контракту.
//!
//! Використання: check --out out/

use fire_harmonizer::data::FxSet;
use fire_harmonizer::date::parse_date;
use fire_harmonizer::gz;
use std::collections::BTreeMap;
use std::io::BufRead;

fn main() {
    let mut out = "out".to_string();
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        match k.as_str() {
            "--out" => out = it.next().expect("--out DIR"),
            "-h" | "--help" => {
                println!("check --out DIR");
                return;
            }
            o => {
                eprintln!("невідома опція {o}");
                std::process::exit(2);
            }
        }
    }
    let mut fails = 0;
    let mut report = |ok: bool, msg: String| {
        println!("{} {msg}", if ok { "[OK]  " } else { "[FAIL]" });
        if !ok {
            fails += 1;
        }
    };

    // ---- 1. cell_day.csv.gz
    let path = format!("{out}/cell_day.csv.gz");
    match gz::open_read(&path) {
        Err(e) => report(false, format!("{path}: {e}")),
        Ok(r) => {
            let mut lines = r.lines();
            let header = lines.next().and_then(|l| l.ok()).unwrap_or_default();
            report(
                header == "cell_id,date,modis_count,viirs_count,count_h,frp_h",
                format!("заголовок cell_day: `{header}`"),
            );
            let mut keys: FxSet<(u32, i32)> = FxSet::default();
            let (mut n, mut dup, mut neg, mut empty, mut bad) = (0u64, 0u64, 0u64, 0u64, 0u64);
            let mut by_year: BTreeMap<i32, (f64, f64, f64)> = BTreeMap::new();
            for l in lines {
                let Ok(l) = l else {
                    bad += 1;
                    continue;
                };
                n += 1;
                let f: Vec<&str> = l.split(',').collect();
                if f.len() != 6 {
                    bad += 1;
                    continue;
                }
                let (Ok(c), Some(d), Ok(m), Ok(v), Ok(ch), Ok(fh)) = (
                    f[0].parse::<u32>(),
                    parse_date(f[1]),
                    f[2].parse::<u32>(),
                    f[3].parse::<u32>(),
                    f[4].parse::<f64>(),
                    f[5].parse::<f64>(),
                ) else {
                    bad += 1;
                    continue;
                };
                if !keys.insert((c, d)) {
                    dup += 1;
                }
                if ch < 0.0 || fh < 0.0 || !ch.is_finite() || !fh.is_finite() {
                    neg += 1;
                }
                if m == 0 && v == 0 {
                    empty += 1;
                }
                let y: i32 = f[1][..4].parse().unwrap_or(0);
                let e = by_year.entry(y).or_default();
                e.0 += m as f64;
                e.1 += v as f64;
                e.2 += ch;
            }
            report(n > 0, format!("рядків даних: {n}"));
            report(bad == 0, format!("некоректних рядків: {bad}"));
            report(dup == 0, format!("дублікатів ключа (cell_id, date): {dup}"));
            report(neg == 0, format!("від'ємних або нескінченних count_h/frp_h: {neg}"));
            report(empty == 0, format!("рядків без жодної точки (modis_count = viirs_count = 0): {empty}"));
            println!("\n  рік    modis_count  viirs_count      count_h");
            for (y, (m, v, c)) in &by_year {
                println!("  {y}  {m:>11.0}  {v:>11.0}  {c:>11.0}");
            }
            println!();
        }
    }

    // ---- 2. model_meta.json
    let mp = format!("{out}/model_meta.json");
    match std::fs::read_to_string(&mp) {
        Err(e) => report(false, format!("{mp}: {e}")),
        Ok(s) => {
            for k in ["model_version", "reference_scale", "fit_years", "holdout_years", "metrics", "notes", "mape", "r2", "bias"] {
                report(s.contains(&format!("\"{k}\":")), format!("model_meta.json має поле `{k}`"));
            }
            report(s.contains("\"VIIRS_SNPP\"") || s.contains("\"MODIS\""), "reference_scale — VIIRS_SNPP або MODIS".into());
        }
    }

    println!("\n{}", if fails == 0 { "УСІ ПЕРЕВІРКИ ПРОЙДЕНО" } else { "Є ПОМИЛКИ — див. [FAIL] вище" });
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
