//! Fire Calendar — гармонізація гарячих точок MODIS (Terra/Aqua) і VIIRS S-NPP.
//!
//! Модулі:
//! - [`data`]  — читання `hotspots.csv.gz` (контракт 1), фільтри, агрегація до клітинка×доба×день/ніч;
//! - [`model`] — ієрархічна модель перерахунку MODIS → шкала VIIRS S-NPP;
//! - [`cv`]    — leave-one-year-out крос-валідація, метрики, підбір гіперпараметрів;
//! - [`date`], [`gz`], [`json`], [`rng`] — допоміжне (без зовнішніх залежностей).

pub mod cv;
pub mod data;
pub mod date;
pub mod gz;
pub mod json;
pub mod model;
pub mod rng;
