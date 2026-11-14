//! Мінімальний JSON-писач (для model_meta.json).

pub enum J {
    Null,
    Bool(bool),
    Num(f64),
    Int(i64),
    Str(String),
    Arr(Vec<J>),
    Obj(Vec<(String, J)>),
}

impl J {
    pub fn obj() -> J {
        J::Obj(Vec::new())
    }
    pub fn set(mut self, k: &str, v: impl Into<J>) -> J {
        if let J::Obj(ref mut o) = self {
            o.push((k.to_string(), v.into()));
        }
        self
    }
    pub fn render(&self) -> String {
        let mut s = String::new();
        self.w(&mut s, 0);
        s.push('\n');
        s
    }
    fn w(&self, s: &mut String, ind: usize) {
        match self {
            J::Null => s.push_str("null"),
            J::Bool(b) => s.push_str(if *b { "true" } else { "false" }),
            J::Int(i) => s.push_str(&i.to_string()),
            J::Num(x) => {
                if x.is_finite() {
                    s.push_str(&fmt_num(*x))
                } else {
                    s.push_str("null")
                }
            }
            J::Str(t) => esc(s, t),
            J::Arr(a) => {
                if a.is_empty() {
                    s.push_str("[]");
                    return;
                }
                let simple = a.iter().all(|x| matches!(x, J::Num(_) | J::Int(_) | J::Str(_) | J::Bool(_) | J::Null));
                if simple {
                    s.push('[');
                    for (i, x) in a.iter().enumerate() {
                        if i > 0 {
                            s.push_str(", ");
                        }
                        x.w(s, ind);
                    }
                    s.push(']');
                } else {
                    s.push_str("[\n");
                    for (i, x) in a.iter().enumerate() {
                        pad(s, ind + 2);
                        x.w(s, ind + 2);
                        if i + 1 < a.len() {
                            s.push(',');
                        }
                        s.push('\n');
                    }
                    pad(s, ind);
                    s.push(']');
                }
            }
            J::Obj(o) => {
                if o.is_empty() {
                    s.push_str("{}");
                    return;
                }
                s.push_str("{\n");
                for (i, (k, v)) in o.iter().enumerate() {
                    pad(s, ind + 2);
                    esc(s, k);
                    s.push_str(": ");
                    v.w(s, ind + 2);
                    if i + 1 < o.len() {
                        s.push(',');
                    }
                    s.push('\n');
                }
                pad(s, ind);
                s.push('}');
            }
        }
    }
}

fn pad(s: &mut String, n: usize) {
    for _ in 0..n {
        s.push(' ');
    }
}

fn esc(s: &mut String, t: &str) {
    s.push('"');
    for c in t.chars() {
        match c {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c if (c as u32) < 0x20 => s.push_str(&format!("\\u{:04x}", c as u32)),
            c => s.push(c),
        }
    }
    s.push('"');
}

fn fmt_num(x: f64) -> String {
    if x == x.trunc() && x.abs() < 1e15 {
        format!("{x:.1}")
    } else {
        let t = format!("{x:.6}");
        let t = t.trim_end_matches('0');
        if t.ends_with('.') {
            format!("{t}0")
        } else {
            t.to_string()
        }
    }
}

impl From<f64> for J {
    fn from(x: f64) -> J {
        J::Num(x)
    }
}
impl From<i64> for J {
    fn from(x: i64) -> J {
        J::Int(x)
    }
}
impl From<i32> for J {
    fn from(x: i32) -> J {
        J::Int(x as i64)
    }
}
impl From<u64> for J {
    fn from(x: u64) -> J {
        J::Int(x as i64)
    }
}
impl From<usize> for J {
    fn from(x: usize) -> J {
        J::Int(x as i64)
    }
}
impl From<bool> for J {
    fn from(x: bool) -> J {
        J::Bool(x)
    }
}
impl From<&str> for J {
    fn from(x: &str) -> J {
        J::Str(x.to_string())
    }
}
impl From<String> for J {
    fn from(x: String) -> J {
        J::Str(x)
    }
}
impl From<Vec<J>> for J {
    fn from(x: Vec<J>) -> J {
        J::Arr(x)
    }
}
