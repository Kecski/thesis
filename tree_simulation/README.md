# BDMM-Prime template + generator

| file | purpose |
|---|---|
| `template.xml` | Your analysis XML with `{{PLACEHOLDERS}}` for data and parameters |
| `make_beast_xml.py` | Fills the template from a FASTA file and a TOML config |
| `config_example.toml` | Example config (reproduces your original analysis) |
| `example.fasta` | Your 359 sequences, as an example input |

## Use
```bash
python make_beast_xml.py config_example.toml            # writes analysis.xml
python make_beast_xml.py config.toml --fasta my.fasta -o my_run.xml
```
Needs Python >= 3.11, no extra packages.

## Input requirements
* **FASTA**: aligned (equal lengths), nucleotides only, unique tip names. The tip name is the
  header up to the first space. Tip names must not contain `,` `=` spaces or quotes.
* **Tip dates** (decimal years): parsed from the tip name with `date_regex` (default: name ends in
  `_2002.2012`), or from a CSV (`dates_csv`, header row, columns `name,date`).

## Re skyline and priors
Everything is set in the `[re_prior]` table (vectors, one entry per interval) and `[model]`:

```toml
[model]
re_change_times = [0.1, 0.2, 0.3]   # n-1 change times (ages, years back from last sample)
# re_interval_length = 0.1          # alternative: equal-length intervals

[re_prior]
distribution       = "LogNormal"    # "LogNormal" or "Normal"; single value or one per interval
mean_in_real_space = true           # LogNormal only; single value or one per interval
mean    = [1.0, 1.1, 1.2, 1.3]
sigma   = [0.3, 0.3, 0.3, 0.3]      # or a single value for all intervals
initial = [1.04, 1.07, 0.94, 1.27]  # optional starting values
```

* **Number of intervals** = length of the vectors. Change the vectors and
  `re_change_times` together; with `re_interval_length` you only change the vectors.
* **Constant Re** (one interval): use single values (`mean = 1.2`, `sigma = 0.4`) or length-1 vectors
  and `re_change_times = []` (or omit it). The `changeTimes` element is then left out of the XML.
* **LogNormal, `mean_in_real_space = true`**: `mean` is the mean of Re, `sigma` the SD of log(Re).
  **`false`**: `mean` is M = mean of log(Re) (median of Re = exp(M)), `sigma` the SD of log(Re).
* **Normal**: `mean` and `sigma` are on the Re scale.
* Rates (e.g. `become_uninfectious_rate`) are per **year**.
* Interval order: entry 1 = `index 0` of the `ReEpi` vector. Check once that this is the interval
  you intend (e.g. give one interval a very tight, extreme prior and see which part of the logged
  `ReSPEpi` moves), then reorder if needed.

## Requirements inside BEAST
Packages: BEAST.base, BDMM-Prime, SA, **feast** (provides `feast.function.Slice`, used for the
per-interval priors). Set `feast_version` in the config to add feast to the `required` attribute.

## Changing the template
Add a `{{NEW_NAME}}` placeholder to `template.xml` and a matching entry to the `values` dict in
`make_beast_xml.py`. The script stops if the template contains a placeholder it does not know.
