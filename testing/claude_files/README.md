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

## Re skyline
* `n` `[[re_interval]]` blocks need `n-1` entries in `re_change_times`.
* Change times are **ages** (years back from the most recent sample), as in your original
  (`timesAreAges="true"`). Rates (e.g. `become_uninfectious_rate`) are per **year**.
* Interval order: block 1 = `index 0` of the `ReEpi` vector. Check once that this is the
  interval you intend (e.g. give one interval a very tight, extreme prior and see which part of
  the logged `ReSPEpi` moves), then reorder the blocks if needed.
* Priors: `LogNormal` (`mean` = mean of Re in real space, `sigma` = SD on log scale) or
  `Normal` (both on the Re scale).

## Requirements inside BEAST
Packages: BEAST.base, BDMM-Prime, SA, **feast** (provides `feast.function.Slice`, used for the
per-interval priors). Set `feast_version` in the config to add feast to the `required` attribute.

## Changing the template
Add a `{{NEW_NAME}}` placeholder to `template.xml` and a matching entry to the `values` dict in
`make_beast_xml.py`. The script stops if the template contains a placeholder it does not know.
