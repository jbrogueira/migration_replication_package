# Eurostat extracts (dissemination API, JSON-stat 2.0)

Extracted 2026-10-07. Each file is the unmodified API response.

| file | dataset | query |
|---|---|---|
| `earn_ses22_30.json` | Mean annual earnings by sex, economic activity and educational attainment (2022), SES; DOI 10.2908/EARN_SES22_30 | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/earn_ses22_30?geo=EU27_2020&sex=T&lang=EN` |
| `earn_ses22_04.json` | Employees by sex, economic activity and educational attainment (2022), SES; DOI 10.2908/EARN_SES22_04 | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/earn_ses22_04?geo=EU27_2020&sex=T&lang=EN` |
| `earn_ses22_16.json` | Mean hourly earnings by sex, economic activity and educational attainment (2022), SES; DOI 10.2908/EARN_SES22_16 (cross-check only) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/earn_ses22_16?geo=EU27_2020&sex=T&lang=EN` |
| `ilc_di08.json` | Mean and median income by educational attainment level, EU-SILC (previous target; equivalised disposable income, not used) | `https://ec.europa.eu/eurostat/api/dissemination/statistics/1.0/data/ilc_di08?geo=EU27_2020&time=2024&lang=EN` |

Skill-premium target (`config.py`): SES 2022, EU27_2020, indic_se = ERN (gross
earnings), sizeclas = GE10 (enterprises with 10 employees or more), nace_r2 =
B-S_X_O (industry, construction and services except public administration),
unit = EUR, sex = T. Mean annual gross earnings: levels 0–2 = 27 165, levels 3–4 =
34 282, levels 5–8 = 50 948. Employees: 21 567 667, 60 849 551, 43 914 560.
Non-college = employee-weighted mean of levels 0–2 and 3–4 = 32 419.56.
Premium = 50 948 / 32 419.56 = 1.5715.
