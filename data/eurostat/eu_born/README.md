# Eurostat extracts: EU-born immigrants (robustness calibration)

Extracted 2026-10-16 from the Eurostat dissemination API (JSON-stat 2.0). Each JSON
file is the unmodified API response.

| file | dataset (release date in the file) | query |
|---|---|---|
| `migr_pop3ctb_EU27_2024.json` | Population on 1 January by age group, sex and country of birth (updated 2026-09-25) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/migr_pop3ctb?geo=EU27_2020&time=2024&sex=T&unit=NR&lang=EN` |
| `edat_lfs_9915_EU27_2024.json` | Population in private households by educational attainment level, country of birth and degree of urbanisation, % (EU-LFS; updated 2026-09-10) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/edat_lfs_9915?geo=EU27_2020&time=2024&sex=T&age=Y25-64&lang=EN` |
| `edat_lfs_9912_EU27_2024.json` | Population in private households by educational attainment level and country of birth, % (EU-LFS; cross-check, identical to 9915 at deg_urb = TOTAL) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/edat_lfs_9912?geo=EU27_2020&time=2024&sex=T&age=Y25-64&lang=EN` |
| `edat_lfs_9912_EU27_2020-2025_ED0-2_ED5-8.json` | same as above, 2020–2025, levels 0–2 and 5–8 (used only to check which year the baseline shares come from) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/edat_lfs_9912?geo=EU27_2020&sex=T&age=Y25-64&isced11=ED5-8&isced11=ED0-2&sinceTimePeriod=2020&lang=EN` |

Country-of-birth codes: `NAT` reporting country (natives); `FOR` foreign country;
`EU27_2020_FOR` EU27 countries (from 2020) except the reporting country (EU-born
immigrants); `NEU27_2020_FOR` non-EU27 countries nor the reporting country.

## Population aged 25–64, EU27_2020, 1 January 2024 (migr_pop3ctb, sex = T, unit = NR)

Sum of the eight five-year groups Y25-29, …, Y60-64.

| country of birth | persons 25–64 |
|---|---|
| NAT | 196 210 940 |
| FOR | 42 027 781 |
| EU27_2020_FOR | 10 919 080 |
| NEU27_2020_FOR | 31 108 701 |

EU27_2020_FOR + NEU27_2020_FOR = FOR exactly. NAT and FOR reproduce the baseline
counts in `config.py` (POP_25_64 = NAT + FOR = 238 238 721,
FOREIGN_25_64 = 42 027 781), so this release of migr_pop3ctb is the one the baseline used.

## Tertiary share (ISCED 5–8), ages 25–64, EU27_2020, 2024 (edat_lfs_9915, deg_urb = TOTAL, %)

| country of birth | levels 0–2 | levels 3–4 | levels 5–8 |
|---|---|---|---|
| NAT | 16.9 | 46.5 | 36.6 |
| FOR | 33.3 | 34.1 | 32.5 |
| EU27_2020_FOR | 25.1 | 40.0 | 34.8 |
| NEU27_2020_FOR | 36.3 | 32.1 | 31.6 |
| TOTAL | 19.6 | 44.4 | 35.9 |

Release difference: the baseline (`config.py`, appendix table `tab:education`) uses
36.8 (NAT) and 32.9 (FOR) for 2024, with 33.1 (FOR, levels 0–2) and 36.1 (TOTAL, levels
5–8). The current release gives 36.6, 32.5, 33.3 and 35.9 for the same cells. No year
2020–2025 in the current release has 36.8/32.9 (file `..._2020-2025_...`), so the
baseline shares come from an earlier release of the 2024 data. The earlier release is
not available from the API and was not recovered. The only flagged cells (`u`, low
reliability) are the "no response" country-of-birth category.

## EU-born targets (used: `MIG_EU_BORN=1`, Appendix N)

Model population = all residents aged 25–64; residents born outside the EU are counted
among the natives:

- EU-born immigrant share, 25–64: 10 919 080 / 238 238 721 = 0.045833 (4.58%).
- Tertiary share among EU-born immigrants: 0.348.
- High-skill EU-born immigrants: 0.045833 × 0.348 = 0.015950 (1.60%).
- Tertiary share among natives: unchanged (0.368 in the baseline release; 0.366 in the
  current release).

Not used: natives + EU-born only as the population (non-EU-born excluded from the
economy): 10 919 080 / 207 130 020 = 0.052716.
