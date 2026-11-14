//! Дати як кількість діб від 1970-01-01 (UTC). Алгоритм Howard Hinnant.

pub fn days_from_civil(y: i32, m: u32, d: u32) -> i32 {
    let y = if m <= 2 { y - 1 } else { y };
    let era = (if y >= 0 { y } else { y - 399 }) / 400;
    let yoe = (y - era * 400) as u32;
    let mp = if m > 2 { m - 3 } else { m + 9 };
    let doy = (153 * mp + 2) / 5 + d - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe as i32 - 719_468
}

pub fn civil_from_days(z: i32) -> (i32, u32, u32) {
    let z = z + 719_468;
    let era = (if z >= 0 { z } else { z - 146_096 }) / 146_097;
    let doe = (z - era * 146_097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i32 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// `YYYY-MM-DD` → доба. Повертає `None` для некоректного рядка.
pub fn parse_date(s: &str) -> Option<i32> {
    let b = s.as_bytes();
    if b.len() < 10 || b[4] != b'-' || b[7] != b'-' {
        return None;
    }
    let num = |r: &[u8]| -> Option<u32> {
        let mut v = 0u32;
        for &c in r {
            if !c.is_ascii_digit() {
                return None;
            }
            v = v * 10 + (c - b'0') as u32;
        }
        Some(v)
    };
    let y = num(&b[0..4])? as i32;
    let m = num(&b[5..7])?;
    let d = num(&b[8..10])?;
    if !(1..=12).contains(&m) || !(1..=31).contains(&d) {
        return None;
    }
    Some(days_from_civil(y, m, d))
}

pub fn fmt_date(z: i32) -> String {
    let (y, m, d) = civil_from_days(z);
    format!("{y:04}-{m:02}-{d:02}")
}

/// День року, 1..=366.
pub fn doy(z: i32) -> u32 {
    let (y, _, _) = civil_from_days(z);
    (z - days_from_civil(y, 1, 1)) as u32 + 1
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip() {
        for z in -1000..40_000 {
            let (y, m, d) = civil_from_days(z);
            assert_eq!(days_from_civil(y, m, d), z);
        }
        assert_eq!(parse_date("1970-01-01"), Some(0));
        assert_eq!(parse_date("2012-01-20"), Some(15_359));
        assert_eq!(fmt_date(15_359), "2012-01-20");
        assert_eq!(doy(parse_date("2020-12-31").unwrap()), 366);
        assert_eq!(parse_date("2020-13-01"), None);
    }
}
