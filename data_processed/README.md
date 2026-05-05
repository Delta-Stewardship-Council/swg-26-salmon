# American & Stanislaus River RST Data Processing

## Overview

Raw rotary screw trap (RST) data for the **American River** and **Stanislaus River** were read into [Positron](https://positron.posit.co/) and processed using **R**.

---

## Processing Date

`2026-05-05`

---

## Processing Steps

1. **Date formatting** — Date fields were converted to `YYYY-MM-DD` format.
2. **Calendar and water year columns** — Calendar year and water year columns were added to the dataset.
3. **Column header standardization** — All column header names were made uniform using lowercase letters with periods between words (e.g., `flow.cfs`, `water.temp`).
4. **Data join** — Clean catch data and environmental data were combined using a `left_join`.

---
