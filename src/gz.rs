//! gzip на чистому Rust, без зовнішніх утиліт і crate'ів — працює однаково на Windows, Linux, macOS.
//!
//! - Читання: повний inflate (stored / fixed / dynamic Huffman), кілька gzip-членів підряд,
//!   перевірка CRC32 і довжини. Сумісний з файлами від `gzip`, Go `compress/gzip`, Python тощо.
//! - Запис: deflate з LZ77 (хеш-ланцюжки, вікно 32 КБ) і фіксованими кодами Хаффмана.

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};

fn bad(msg: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("gzip: {msg}"))
}

// ---------------------------------------------------------------- CRC32

fn crc_table() -> [u32; 256] {
    let mut t = [0u32; 256];
    for (n, e) in t.iter_mut().enumerate() {
        let mut c = n as u32;
        for _ in 0..8 {
            c = if c & 1 != 0 { 0xEDB8_8320 ^ (c >> 1) } else { c >> 1 };
        }
        *e = c;
    }
    t
}

pub fn crc32(data: &[u8]) -> u32 {
    let t = crc_table();
    let mut c = 0xFFFF_FFFFu32;
    for &b in data {
        c = t[((c ^ b as u32) & 0xFF) as usize] ^ (c >> 8);
    }
    c ^ 0xFFFF_FFFF
}

// ---------------------------------------------------------------- спільні таблиці deflate

const LBASE: [u16; 29] = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258];
const LEXT: [u8; 29] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
const DBASE: [u16; 30] = [
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289,
    16385, 24577,
];
const DEXT: [u8; 30] = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13];

fn fixed_litlen_lengths() -> [u8; 288] {
    let mut l = [0u8; 288];
    for (i, x) in l.iter_mut().enumerate() {
        *x = match i {
            0..=143 => 8,
            144..=255 => 9,
            256..=279 => 7,
            _ => 8,
        };
    }
    l
}

fn rev(mut code: u32, len: u32) -> u32 {
    let mut r = 0;
    for _ in 0..len {
        r = (r << 1) | (code & 1);
        code >>= 1;
    }
    r
}

/// Канонічні коди Хаффмана з довжин (RFC 1951, 3.2.2).
fn canonical(lengths: &[u8]) -> Vec<u32> {
    let mut bl_count = [0u32; 16];
    for &l in lengths {
        bl_count[l as usize] += 1;
    }
    bl_count[0] = 0;
    let mut next = [0u32; 16];
    let mut code = 0;
    for b in 1..16 {
        code = (code + bl_count[b - 1]) << 1;
        next[b] = code;
    }
    lengths
        .iter()
        .map(|&l| {
            if l == 0 {
                0
            } else {
                let c = next[l as usize];
                next[l as usize] += 1;
                c
            }
        })
        .collect()
}

// ---------------------------------------------------------------- inflate

struct Bits<'a> {
    data: &'a [u8],
    pos: usize,
    buf: u64,
    cnt: u32,
}

impl<'a> Bits<'a> {
    #[inline]
    fn fill(&mut self) {
        while self.cnt <= 56 && self.pos < self.data.len() {
            self.buf |= (self.data[self.pos] as u64) << self.cnt;
            self.pos += 1;
            self.cnt += 8;
        }
    }
    #[inline]
    fn bits(&mut self, n: u32) -> io::Result<u32> {
        if n == 0 {
            return Ok(0);
        }
        if self.cnt < n {
            self.fill();
            if self.cnt < n {
                return Err(bad("неочікуваний кінець даних"));
            }
        }
        let v = (self.buf & ((1u64 << n) - 1)) as u32;
        self.buf >>= n;
        self.cnt -= n;
        Ok(v)
    }
    fn align_byte(&mut self) {
        let drop = self.cnt % 8;
        self.buf >>= drop;
        self.cnt -= drop;
    }
    /// Кількість байтів вхідних даних, реально спожитих (після вирівнювання на байт).
    fn consumed(&self) -> usize {
        self.pos - (self.cnt / 8) as usize
    }
}

struct Huff {
    table: Vec<u32>, // (sym << 4) | len; 0 = недійсний код
    maxlen: u32,
}

impl Huff {
    fn new(lengths: &[u8]) -> io::Result<Huff> {
        let maxlen = lengths.iter().copied().max().unwrap_or(0) as u32;
        if maxlen == 0 {
            return Ok(Huff { table: vec![0; 2], maxlen: 1 });
        }
        if maxlen > 15 {
            return Err(bad("довжина коду > 15"));
        }
        let codes = canonical(lengths);
        let mut table = vec![0u32; 1 << maxlen];
        for (sym, (&l, &c)) in lengths.iter().zip(&codes).enumerate() {
            if l == 0 {
                continue;
            }
            let l = l as u32;
            let r = rev(c, l);
            let step = 1usize << l;
            let mut k = r as usize;
            while k < table.len() {
                table[k] = ((sym as u32) << 4) | l;
                k += step;
            }
        }
        Ok(Huff { table, maxlen })
    }
    #[inline]
    fn decode(&self, b: &mut Bits) -> io::Result<u32> {
        if b.cnt < self.maxlen {
            b.fill();
        }
        let e = self.table[(b.buf & ((1u64 << self.maxlen) - 1)) as usize];
        let l = e & 15;
        if l == 0 || l > b.cnt {
            return Err(bad("пошкоджений код Хаффмана"));
        }
        b.buf >>= l;
        b.cnt -= l;
        Ok(e >> 4)
    }
}

/// Розпаковує один deflate-потік; повертає кількість спожитих байтів.
fn inflate(data: &[u8], out: &mut Vec<u8>) -> io::Result<usize> {
    let mut b = Bits { data, pos: 0, buf: 0, cnt: 0 };
    let fixed_lit = Huff::new(&fixed_litlen_lengths())?;
    let fixed_dist = Huff::new(&[5u8; 30])?;
    loop {
        let last = b.bits(1)?;
        let btype = b.bits(2)?;
        match btype {
            0 => {
                b.align_byte();
                let len = b.bits(16)?;
                let nlen = b.bits(16)?;
                if len != (!nlen & 0xFFFF) {
                    return Err(bad("stored-блок: LEN/NLEN не збігаються"));
                }
                for _ in 0..len {
                    out.push(b.bits(8)? as u8);
                }
            }
            1 => codes(&mut b, out, &fixed_lit, &fixed_dist)?,
            2 => {
                let hlit = b.bits(5)? as usize + 257;
                let hdist = b.bits(5)? as usize + 1;
                let hclen = b.bits(4)? as usize + 4;
                const ORD: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
                let mut cl = [0u8; 19];
                for &o in ORD.iter().take(hclen) {
                    cl[o] = b.bits(3)? as u8;
                }
                let clh = Huff::new(&cl)?;
                let mut lens = vec![0u8; hlit + hdist];
                let mut i = 0;
                while i < hlit + hdist {
                    let sym = clh.decode(&mut b)?;
                    match sym {
                        0..=15 => {
                            lens[i] = sym as u8;
                            i += 1;
                        }
                        16 => {
                            if i == 0 {
                                return Err(bad("повтор без попередньої довжини"));
                            }
                            let prev = lens[i - 1];
                            for _ in 0..3 + b.bits(2)? {
                                if i >= lens.len() {
                                    return Err(bad("забагато довжин"));
                                }
                                lens[i] = prev;
                                i += 1;
                            }
                        }
                        17 | 18 => {
                            let n = if sym == 17 { 3 + b.bits(3)? } else { 11 + b.bits(7)? };
                            for _ in 0..n {
                                if i >= lens.len() {
                                    return Err(bad("забагато довжин"));
                                }
                                lens[i] = 0;
                                i += 1;
                            }
                        }
                        _ => return Err(bad("недійсний символ довжини коду")),
                    }
                }
                let lit = Huff::new(&lens[..hlit])?;
                let dist = Huff::new(&lens[hlit..])?;
                codes(&mut b, out, &lit, &dist)?;
            }
            _ => return Err(bad("недійсний тип блоку")),
        }
        if last == 1 {
            b.align_byte();
            return Ok(b.consumed());
        }
    }
}

fn codes(b: &mut Bits, out: &mut Vec<u8>, lit: &Huff, dist: &Huff) -> io::Result<()> {
    loop {
        let sym = lit.decode(b)?;
        if sym < 256 {
            out.push(sym as u8);
        } else if sym == 256 {
            return Ok(());
        } else {
            let s = (sym - 257) as usize;
            if s >= 29 {
                return Err(bad("недійсний символ довжини"));
            }
            let len = LBASE[s] as usize + b.bits(LEXT[s] as u32)? as usize;
            let ds = dist.decode(b)? as usize;
            if ds >= 30 {
                return Err(bad("недійсний символ відстані"));
            }
            let d = DBASE[ds] as usize + b.bits(DEXT[ds] as u32)? as usize;
            if d > out.len() {
                return Err(bad("відстань виходить за межі"));
            }
            let start = out.len() - d;
            if d >= len {
                out.extend_from_within(start..start + len);
            } else {
                for k in 0..len {
                    let c = out[start + k];
                    out.push(c);
                }
            }
        }
    }
}

/// Розпаковує gzip (один або кілька членів підряд) з перевіркою CRC32.
pub fn gunzip(data: &[u8]) -> io::Result<Vec<u8>> {
    let mut out = Vec::with_capacity(data.len() * 4);
    let mut p = 0;
    while p < data.len() {
        let d = &data[p..];
        if d.len() < 18 || d[0] != 0x1f || d[1] != 0x8b {
            if p > 0 && d.iter().all(|&x| x == 0) {
                break; // нульове доповнення в кінці файлу
            }
            return Err(bad("це не gzip-файл"));
        }
        if d[2] != 8 {
            return Err(bad("непідтримуваний метод стиснення"));
        }
        let flg = d[3];
        let mut h = 10;
        if flg & 4 != 0 {
            let xlen = d[h] as usize | (d[h + 1] as usize) << 8;
            h += 2 + xlen;
        }
        for bit in [8u8, 16] {
            if flg & bit != 0 {
                while h < d.len() && d[h] != 0 {
                    h += 1;
                }
                h += 1;
            }
        }
        if flg & 2 != 0 {
            h += 2;
        }
        if h >= d.len() {
            return Err(bad("пошкоджений заголовок"));
        }
        let start = out.len();
        let used = inflate(&d[h..], &mut out)?;
        let t = h + used;
        if t + 8 > d.len() {
            return Err(bad("немає трейлера CRC32/ISIZE"));
        }
        let crc = u32::from_le_bytes([d[t], d[t + 1], d[t + 2], d[t + 3]]);
        let isize = u32::from_le_bytes([d[t + 4], d[t + 5], d[t + 6], d[t + 7]]);
        if crc32(&out[start..]) != crc {
            return Err(bad("CRC32 не збігається — файл пошкоджений"));
        }
        if (out.len() - start) as u32 != isize {
            return Err(bad("довжина не збігається"));
        }
        p += t + 8;
    }
    Ok(out)
}

// ---------------------------------------------------------------- deflate (LZ77 + фіксований Хаффман)

struct BitW {
    out: Vec<u8>,
    buf: u64,
    cnt: u32,
}
impl BitW {
    #[inline]
    fn put(&mut self, v: u32, n: u32) {
        self.buf |= (v as u64) << self.cnt;
        self.cnt += n;
        while self.cnt >= 8 {
            self.out.push(self.buf as u8);
            self.buf >>= 8;
            self.cnt -= 8;
        }
    }
    fn flush(&mut self) {
        if self.cnt > 0 {
            self.out.push(self.buf as u8);
            self.buf = 0;
            self.cnt = 0;
        }
    }
}

pub fn deflate_fixed(data: &[u8]) -> Vec<u8> {
    const WSIZE: usize = 1 << 15;
    const HBITS: u32 = 16;
    const MAX_CHAIN: usize = 48;
    let lens = fixed_litlen_lengths();
    let codes = canonical(&lens);
    let lit: Vec<(u32, u32)> = (0..288).map(|i| (rev(codes[i], lens[i] as u32), lens[i] as u32)).collect();
    let dcode: Vec<u32> = (0..30).map(|i| rev(i, 5)).collect();
    let mut len_sym = [0u8; 259];
    for l in 3..=258usize {
        let mut s = 0;
        while s + 1 < 29 && LBASE[s + 1] as usize <= l {
            s += 1;
        }
        len_sym[l] = s as u8;
    }
    let dist_sym = |d: usize| -> usize {
        let mut s = 0;
        while s + 1 < 30 && DBASE[s + 1] as usize <= d {
            s += 1;
        }
        s
    };

    let mut w = BitW { out: Vec::with_capacity(data.len() / 3), buf: 0, cnt: 0 };
    w.put(1, 1); // BFINAL
    w.put(1, 2); // BTYPE = 01 (фіксований Хаффман)
    let n = data.len();
    let mut head = vec![u32::MAX; 1 << HBITS];
    let mut prev = vec![u32::MAX; WSIZE];
    let hash = |i: usize| -> usize {
        let v = (data[i] as u32) | (data[i + 1] as u32) << 8 | (data[i + 2] as u32) << 16;
        (v.wrapping_mul(0x9E37_79B1) >> (32 - HBITS)) as usize
    };
    let insert = |i: usize, head: &mut Vec<u32>, prev: &mut Vec<u32>| {
        if i + 2 < n {
            let h = hash(i);
            prev[i & (WSIZE - 1)] = head[h];
            head[h] = i as u32;
        }
    };
    let mut i = 0;
    while i < n {
        let (mut best_len, mut best_d) = (0usize, 0usize);
        if i + 2 < n {
            let h = hash(i);
            let mut cand = head[h];
            let mut chain = 0;
            let max = (n - i).min(258);
            while cand != u32::MAX && chain < MAX_CHAIN {
                let c = cand as usize;
                let d = i - c;
                if d > WSIZE - 1 || d == 0 {
                    break;
                }
                if data[c + best_len.min(max - 1)] == data[i + best_len.min(max - 1)] {
                    let mut l = 0;
                    while l < max && data[c + l] == data[i + l] {
                        l += 1;
                    }
                    if l > best_len {
                        best_len = l;
                        best_d = d;
                        if l == max {
                            break;
                        }
                    }
                }
                let nxt = prev[c & (WSIZE - 1)];
                if nxt != u32::MAX && (nxt as usize) >= c {
                    break;
                }
                cand = nxt;
                chain += 1;
            }
        }
        if best_len >= 3 {
            let s = len_sym[best_len] as usize;
            let (c, l) = lit[257 + s];
            w.put(c, l);
            w.put((best_len - LBASE[s] as usize) as u32, LEXT[s] as u32);
            let ds = dist_sym(best_d);
            w.put(dcode[ds], 5);
            w.put((best_d - DBASE[ds] as usize) as u32, DEXT[ds] as u32);
            for k in i..i + best_len {
                insert(k, &mut head, &mut prev);
            }
            i += best_len;
        } else {
            let (c, l) = lit[data[i] as usize];
            w.put(c, l);
            insert(i, &mut head, &mut prev);
            i += 1;
        }
    }
    let (c, l) = lit[256];
    w.put(c, l);
    w.flush();
    w.out
}

pub fn gzip(data: &[u8]) -> Vec<u8> {
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 255];
    out.extend_from_slice(&deflate_fixed(data));
    out.extend_from_slice(&crc32(data).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

// ---------------------------------------------------------------- API для решти коду

/// Відкриває файл на читання: `-` = stdin, `*.gz` = розпакування в пам'яті, інакше звичайний файл.
pub fn open_read(path: &str) -> io::Result<Box<dyn BufRead>> {
    const CAP: usize = 1 << 20;
    if path == "-" {
        return Ok(Box::new(BufReader::with_capacity(CAP, io::stdin())));
    }
    if path.ends_with(".gz") {
        let mut raw = Vec::new();
        File::open(path)?.read_to_end(&mut raw)?;
        let data = gunzip(&raw)?;
        return Ok(Box::new(Cursor::new(data)));
    }
    Ok(Box::new(BufReader::with_capacity(CAP, File::open(path)?)))
}

pub enum Out {
    Plain(BufWriter<File>),
    Gz { path: String, buf: Vec<u8> },
}

/// Створює файл на запис; для `*.gz` дані стискаються під час `finish()`.
pub fn create(path: &str) -> io::Result<Out> {
    if path.ends_with(".gz") {
        File::create(path)?; // одразу перевіряємо, що шлях доступний для запису
        Ok(Out::Gz { path: path.to_string(), buf: Vec::with_capacity(1 << 24) })
    } else {
        Ok(Out::Plain(BufWriter::with_capacity(1 << 20, File::create(path)?)))
    }
}

impl Out {
    pub fn finish(self) -> io::Result<()> {
        match self {
            Out::Plain(mut w) => w.flush(),
            Out::Gz { path, buf } => std::fs::write(path, gzip(&buf)),
        }
    }
}

impl Write for Out {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        match self {
            Out::Plain(w) => w.write(data),
            Out::Gz { buf, .. } => {
                buf.extend_from_slice(data);
                Ok(data.len())
            }
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self {
            Out::Plain(w) => w.flush(),
            Out::Gz { .. } => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip() {
        let mut s = String::new();
        for i in 0..20_000 {
            s.push_str(&format!("{},2020-04-{:02},{},{},{:.3}\n", 4_934_176 + i % 37, i % 28 + 1, i % 5, i % 11, i as f64 * 0.37));
        }
        for data in [s.as_bytes(), b"" as &[u8], b"a", b"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"] {
            let z = gzip(data);
            assert_eq!(gunzip(&z).unwrap(), data);
        }
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
    }
}
