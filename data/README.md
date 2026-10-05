# Data: `psilodep.csv`

A copy of `data.csv` from the **`data-depression-psiloctr`** dataset (Metapsy / Sypres
Collaboration). It holds psilocybin-assisted therapy vs. control comparisons for adults with
depressive symptoms. The last search update was on 2025-10-31.

- Repository: <https://github.com/metapsy-project/data-depression-psiloctr>
- DOI: [10.5281/zenodo.15714852](https://doi.org/10.5281/zenodo.15714852)
- Analysis pipeline: <https://sypres.io/docs/datasets/psilodep-meta-analysis/>
- Data standard (variable definitions): <https://docs.metapsy.org/data-preparation/format/>
- License: Open Data Commons Attribution License (ODC-By 1.0),
  <http://opendatacommons.org/licenses/by/1.0/>

**Format note:** the file is separated by semicolons and uses a decimal comma (`27,670`), so read it with
`read.csv2()` and not `read.csv()`.

Each row is one comparison × one instrument × one time point. There are 221 rows from 15 studies.
