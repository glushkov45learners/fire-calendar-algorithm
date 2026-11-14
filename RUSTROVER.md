# Як запустити fire-harmonizer у RustRover

Покрокова інструкція: від встановлення до готових файлів `cell_day.csv.gz` і `model_meta.json` для передачі Міші (loader → БД).

**Що вийде в результаті.** Модель (версія v3) прочитає `hotspots.csv` (2,39 млн точок) і приблизно за 2–5 хвилин (залежно від кількості ядер процесора) створить теку `out/`:

- `cell_day.csv.gz` — гармонізований ряд за контрактом 2 (~737 тис. рядків);
- `model_meta.json` — опис моделі, метрики, параметри;
- `diagnostics/` — таблиці для графіків у презентації.

**Системні вимоги:** Windows 10/11, macOS або Linux; ~2 ГБ вільного місця (з них ~1,5 ГБ займе Rust); ~300 МБ оперативної пам'яті під час роботи. Інтернет потрібен лише для встановлення Rust; сам проєкт не має зовнішніх залежностей і збирається офлайн.

---

## Крок 1. Встановити Rust (компілятор і cargo)

RustRover — це лише редактор. Сам компілятор Rust ставиться окремо, через **rustup**.

### Windows

1. Відкрий **https://rustup.rs** і завантаж `rustup-init.exe` (64-bit).
2. Запусти його. Якщо побачиш повідомлення, що бракує **Visual Studio C++ Build Tools**:
   - вибери варіант **1) Quick install via the Visual Studio Community installer**;
   - у вікні Visual Studio Installer має бути позначено **«Desktop development with C++»** (Розробка класичних застосунків на C++). Встанови, це ~2–6 ГБ і 5–15 хвилин;
   - після цього rustup продовжить сам.
3. На питання про тип встановлення натисни **Enter** (варіант 1, *default*).
4. Дочекайся `Rust is installed now. Great!` і закрий вікно.

### macOS / Linux

У терміналі:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Натисни **Enter** на питання про тип встановлення. На macOS система може попросити встановити Command Line Tools — погодься.

### Перевірка

Відкрий **новий** термінал (PowerShell на Windows) і виконай:

```bash
rustc --version
cargo --version
```

Має з'явитися щось на кшталт `rustc 1.9x.x` і `cargo 1.9x.x`. Підійде будь-яка версія від 1.75.

> Якщо RustRover уже був відкритий під час встановлення Rust — **перезапусти його**, інакше він не побачить cargo.

---

## Крок 2. Встановити RustRover

1. Завантаж з **https://www.jetbrains.com/rust/** і встанови.
2. Під час першого запуску вибери ліцензію **Non-commercial use** (безкоштовно для навчання й некомерційних проєктів) і увійди в акаунт JetBrains.

> Якщо на кроці 4 RustRover скаже, що не знайшов тулчейн, і запропонує **Install Rustup** — це те саме, що крок 1; можна погодитись.

---

## Крок 3. Підготувати теку проєкту

1. **Розпакуй архів** `fire-harmonizer-v3.tar.gz`:
   - Windows 10/11: у PowerShell у теці з архівом виконай `tar -xzf fire-harmonizer-v3.tar.gz` (tar вбудований у Windows), або розпакуй через 7-Zip;
   - macOS / Linux: `tar -xzf fire-harmonizer-v3.tar.gz`.
2. Отримаєш теку `fire-harmonizer` з такою структурою:

   ```
   fire-harmonizer/
   ├── Cargo.toml          ← за цим файлом RustRover розпізнає проєкт
   ├── README.md           ← опис моделі й опцій
   ├── RUSTROVER.md        ← ця інструкція
   ├── .run/               ← 5 готових конфігурацій запуску
   ├── data/               ← сюди покласти hotspots.csv
   ├── examples/           ← результати на даних України для порівняння
   ├── docs/history/       ← як модель розвивалась (v1 → v4)
   └── src/
       ├── main.rs         ← програма harmonize
       ├── model.rs, cv.rs, data.rs, gz.rs, …
       └── bin/
           └── check.rs    ← перевірка результатів
   ```

3. **Скопіюй `hotspots.csv` у теку `data/`**, щоб вийшло `fire-harmonizer/data/hotspots.csv`.
   - Підійде і стиснений `hotspots.csv.gz`, але тоді в конфігураціях 2, 4 і 5 заміни `data/hotspots.csv` на `data/hotspots.csv.gz` (див. крок 7).

> Шлях до теки краще тримати **без кирилиці та пробілів** (наприклад `C:\dev\fire-harmonizer`). Програма з кирилицею впорається, але деякі інструменти збірки на Windows — не завжди.

---

## Крок 4. Відкрити проєкт у RustRover

1. **File → Open…** (або **Open** на стартовому екрані).
2. Вибери **теку** `fire-harmonizer` (ту, де лежить `Cargo.toml`) і натисни **OK**.
3. На питання **Trust Project?** натисни **Trust Project**.
4. Дочекайся, поки в нижньому рядку стану зникне індикатор індексації. Зазвичай це 10–60 секунд.
5. Перевір, що RustRover бачить Rust: **File → Settings** (macOS: **RustRover → Settings**) → **Rust**. У полі **Toolchain location** має стояти шлях на кшталт:
   - Windows: `C:\Users\<ім'я>\.cargo\bin`
   - macOS / Linux: `/Users/<ім'я>/.cargo/bin` або `/home/<ім'я>/.cargo/bin`

   і поруч версія `rustc`. Якщо поле порожнє — вкажи цей шлях вручну і натисни **Apply**.

---

## Крок 5. Знайти готові конфігурації запуску

У правому верхньому куті вікна, ліворуч від зеленої кнопки ▶, є випадний список конфігурацій. У ньому має бути п'ять пунктів:

| Конфігурація | Команда cargo | Що робить |
|---|---|---|
| **1. Тести** | `test --release` | юніт-тести (дати, парсинг, маска промислових джерел, метрики, gzip) |
| **2. Гармонізувати (повний підбір)** | `run --release --bin harmonize -- --input data/hotspots.csv --out out` | **основний запуск**: підбір гіперпараметрів + результат (2–5 хв) |
| **3. Перевірити результат** | `run --release --bin check -- --out out` | перевіряє `out/` на відповідність контракту 2 |
| **4. Швидкий прогін (без підбору)** | `run --release --bin harmonize -- --input data/hotspots.csv --out out_fast --no-tune --bootstrap 50` | готові параметри v3, ~15 с. Для даних України `cell_day.csv.gz` такий самий, як у конфігурації 2 |
| **5. Чесна оцінка (вкладена CV)** | `run --release --bin harmonize -- --input data/hotspots.csv --out out_nested --nested` | те саме, що 2, плюс вкладена крос-валідація для звіту (~15 хв на 2 ядрах) |

Запуск: вибери конфігурацію у списку → натисни зелений **▶ (Run)** або **Shift+F10**. Вивід з'явиться у вікні **Run** унизу.

### Якщо конфігурацій у списку немає

Створи їх вручну (одна хвилина на кожну):

1. **Run → Edit Configurations…**
2. **+** (Add New Configuration) → **Cargo**.
3. **Name:** будь-яка назва, наприклад `2. Гармонізувати`.
4. **Command:** скопіюй з колонки «Команда cargo» таблиці вище, наприклад:
   ```
   run --release --bin harmonize -- --input data/hotspots.csv --out out
   ```
5. **Working directory:** тека проєкту `fire-harmonizer` (там, де `Cargo.toml`). Зазвичай RustRover підставляє її сам.
6. **OK**.

> **Важливо:** завжди залишай `--release`. Без нього програма збирається в режимі налагодження і працює в 10–30 разів повільніше.

---

## Крок 6. Запустити по порядку

### 6.1. Тести — «1. Тести»

Перший запуск ще й компілює проєкт (~1–2 хвилини). Наступні — секунди. Очікуваний кінець виводу:

```
test result: ok. 4 passed; 0 failed; ...
```

### 6.2. Гармонізація — «2. Гармонізувати (повний підбір)»

Очікуваний вивід (час на твоєму комп'ютері буде іншим):

```
[1/7] Читаю data/hotspots.csv …
      рядків 2391912, некоректних 0, type 2/3 виключено 276373, NRT без type 147252 (з них замасковано 11011), записів клітинка×доба×Д/Н 773601 — 3.1 с
      маска стабільних джерел (≥30 місяців, ≥8 років): 112 мікроклітинок, прибрано 14233 точок
      Terra   2000-11-01 … 2026-11-10  (…)
      …
[2/7] Період перекриття 2012-01-25 … 2026-11-01; шукаю дні збою сенсорів …
      днів збою S-NPP: 51, днів збою MODIS у перекритті: 5, підозрілих прогалин до 2012-01-25: 12
[3/7] LOYO-крос-валідація на [2012, …, 2019] і підбір гіперпараметрів (земний покрив у даних: так) …
    terra_aqua λ=30 sm=0.25 … пропущ=0 EM=0 → WAPE 2°×міс=…        ← ~100 рядків підбору, це нормально
      вибрано Terra+Aqua: Hyper { lambda: 100.0, smooth: 0.5, bins: true, modis_low: true, depth: 4, frp_mode: 4, pattern: false, intensity: true, lc: false, miss_alloc: 3, miss_model: 2 }
      вибрано Terra:      Hyper { lambda: 100.0, smooth: 0.25, bins: true, modis_low: true, depth: 3, frp_mode: 4, pattern: false, intensity: true, lc: true, miss_alloc: 3, miss_model: 2 }
      WAPE 2°×міс: константа 0.2371 → модель 0.1420 (Terra+Aqua);  константа 0.2363 → модель 0.1780 (Terra)
[4/7] Фінальне навчання і стрес-тест на [2020, 2021, 2022, 2023, 2024, 2025, 2026] …
[5/7] Bootstrap по роках: 200 реплік; σ року: Terra+Aqua 0.0507, Terra 0.0504 …
[6/7] Пишу out/cell_day.csv.gz …
      рядків 736739
[7/7] Діагностика і model_meta.json …
Готово за 135.0 с. LOYO 2°×міс: WAPE 0.142, bias +0.002, R²(log) 0.848, MAPE 0.267
```

> Програма пише хід роботи в **stderr**, тому у вікні Run ці рядки можуть бути **червоними**. Це не помилка. Справжня помилка виглядає як рядок `ПОМИЛКА: …` і код завершення `exit code 1`.

### 6.3. Перевірка — «3. Перевірити результат»

```
[OK]   заголовок cell_day: `cell_id,date,modis_count,viirs_count,count_h,frp_h`
[OK]   рядків даних: 736739
[OK]   некоректних рядків: 0
[OK]   дублікатів ключа (cell_id, date): 0
[OK]   від'ємних або нескінченних count_h/frp_h: 0
[OK]   рядків без жодної точки (modis_count = viirs_count = 0): 0
  … таблиця сум по роках …
[OK]   model_meta.json має поле `model_version`
  …
УСІ ПЕРЕВІРКИ ПРОЙДЕНО
```

---

## Крок 7. Змінити параметри запуску (за потреби)

**Run → Edit Configurations…** → вибери конфігурацію 2 → допиши опції в поле **Command** після `--`. Приклади:

| Що потрібно | Command |
|---|---|
| Інша назва версії в `model_meta.json` | `run --release --bin harmonize -- --input data/hotspots.csv --out out --model-version v3.1` |
| Вхід стиснений | `… -- --input data/hotspots.csv.gz --out out` |
| Інші роки навчання | `… -- --input data/hotspots.csv --out out --fit-years 2013-2019` |
| Без підбору параметрів (готові параметри v3) | `… -- --input data/hotspots.csv --out out --no-tune` |
| Без маски стабільних джерел | `… -- --input data/hotspots.csv --out out --persistent-mask 0` |
| Чесні метрики з вкладеною крос-валідацією (~×8 довше) | `… -- --input data/hotspots.csv --out out --nested` |
| Не брати VIIRS з low confidence | `… -- --input data/hotspots.csv --out out --viirs-drop-low` |

Повна довідка: конфігурація з командою `run --release --bin harmonize -- --help`.

---

## Крок 8. Що лежить в `out/` і як це читати

Теку `out/` видно в дереві проєкту ліворуч; якщо її немає — натисни правою кнопкою на корінь проєкту → **Reload from Disk**.

| Файл | Кому | Зміст |
|---|---|---|
| `cell_day.csv.gz` | Міша (loader) | контракт 2: `cell_id, date, modis_count, viirs_count, count_h, frp_h` |
| `model_meta.json` | Міша (loader) | версія, шкала, роки, метрики, примітки (`notes`), коефіцієнти, дні збою |
| `diagnostics/annual.csv` | презентація | сирі ряди по супутниках і гармонізований `count_h` з інтервалом 5–95 % |
| `diagnostics/ratio_by_year.csv` | презентація | відношення VIIRS/MODIS по роках — видно дрейф орбіт |
| `diagnostics/validation_monthly.csv` | презентація | місяць × рік: факт S-NPP проти прогнозу з MODIS |
| `diagnostics/cv_grid.csv` | судді | усі перевірені варіанти моделі з метриками |
| `diagnostics/subgroups.csv` | аналіз | якість у підгрупах: сезони, день/ніч, земний покрив, рівень активності |
| `diagnostics/coefficients.csv` | відтворюваність | коефіцієнти по стратах |

`.csv` і `.json` відкриваються прямо в RustRover подвійним кліком. `.gz` — у будь-якому архіваторі або через `check`.

Чесні метрики (конфігурація 5) — у `out_nested/model_meta.json` → розділ `cv_nested`. Для порівняння: у `examples/` лежать результати того самого прогону на даних України; історія версій v1 → v4 — у `docs/history/`.

### На що дивитися, щоб переконатися, що все гаразд

| Перевірка | Де | Норма на наших даних |
|---|---|---|
| Модель краща за константу | лог, рядок `WAPE 2°×міс: константа … → модель …` | 0,237 → **0,142** |
| Немає систематичного зміщення | `model_meta.json` → `metrics.bias` | **+0,002** (норма: \|bias\| < 0,05) |
| Немає стрибка у 2012 році | `annual.csv`: `count_from_modis` проти `snpp_observed` на 2012–2019 | відношення 0,90–1,10 по роках, у середньому ≈ 1,0 |
| Дрейф орбіт видно | `ratio_by_year.csv`: `k_terra` | ~6,2–7,5 до 2019 року → 9,8–11,4 у 2024–2026 |
| Збій S-NPP знайдено | `model_meta.json` → `viirs_gap_days` | є діапазон 2022-07-27 … 2022-08-10 |

---

## Крок 9. Передати результат

Файли `out/cell_day.csv.gz` і `out/model_meta.json` кладуться в бакет за контрактом:

```bash
gsutil cp out/cell_day.csv.gz out/model_meta.json gs://<бакет>/v1/model/v3/
```

(`v3` у кінці шляху має збігатися з `--model-version`.)

---

## Альтернатива: запуск через термінал RustRover

Унизу вікна є вкладка **Terminal** (або **Alt+F12**). Ті самі дії командами:

```bash
cargo test --release
cargo run --release --bin harmonize -- --input data/hotspots.csv --out out
cargo run --release --bin check -- --out out
```

---

## Типові проблеми

| Симптом | Причина | Що робити |
|---|---|---|
| `error: linker 'link.exe' not found` (Windows) | Не встановлені C++ Build Tools | Встанови **Visual Studio Build Tools** з компонентом «Desktop development with C++», перезапусти RustRover |
| `cargo: command not found` / RustRover не бачить тулчейн | Rust встановлено після запуску IDE або PATH не оновився | Перезапусти RustRover; у **Settings → Rust → Toolchain location** вкажи `…\.cargo\bin` |
| `ПОМИЛКА: data/hotspots.csv: … (os error 2)` | Файлу немає за цим шляхом | Перевір, що файл лежить саме в `fire-harmonizer/data/` і **Working directory** конфігурації — тека проєкту |
| `ПОМИЛКА: у заголовку немає колонки …` | Файл не за контрактом 1 | Перевір перший рядок CSV; зайві колонки (як `lc`) не заважають |
| Конфігурацій 1–5 немає в списку | IDE не підхопила теку `.run/` | Створи вручну (крок 5) або запускай з **Terminal** |
| Кирилиця у вікні Run показується «кракозябрами» | Кодування консолі | **Settings → Editor → General → Console → Default Encoding: UTF-8**. На самі файли це не впливає — вони завжди в UTF-8 |
| Дуже повільно (хвилини замість секунд) | Запуск без `--release` | Перевір, що в Command є `--release` |
| Антивірус блокує `.exe` у `target\` (Windows) | Новий незнайомий файл | Додай теку проєкту у винятки антивірусу |

Якщо з'являється інша помилка — скопіюй повний текст із вікна Run, за ним легко знайти причину.
