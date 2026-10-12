# Notebooks

This directory contains the notebooks used to calculate paper values and make figures.

Read the whole-genome notebooks at [Phasing T2T notebooks](https://josephlalli.github.io/phasing_T2T/):

- [Manuscript values](https://josephlalli.github.io/phasing_T2T/calc_figures_for_paper.html)
- [Per-variant calculations](https://josephlalli.github.io/phasing_T2T/calc_per_variant_figures_for_paper.html)
- [Publication figures](https://josephlalli.github.io/phasing_T2T/make_plots.html)

The HTML files in `notebooks_whole_genome/` can also be opened locally. They show saved notebook outputs; `make_plots.html` shows the Docker figure run.

`calc_figures_for_paper.ipynb` calculates all values reported in the text of the paper. It only uses smaller summary parquet files and is designed to run on a laptop in Docker.

`calc_per_variant_figures_for_paper.ipynb` calculates values that require loading the full variant dataframes into memory. It requires at least 400GB of RAM.

The notebooks in this directory use the test-region data. `notebooks_whole_genome/` contains the whole-genome notebooks and figure code.

To produce figures from the supplied whole-genome inputs, run this from the repository root on your host machine:

```bash
bash scripts/run_make_plots.sh
```

The runner checks the 188 supplied inputs and executes only `notebooks_whole_genome/make_plots.ipynb` in Docker. Generate Figure 5 separately with `scripts/figure5/Figure_5_script.R`, then combine its panels with `scripts/figure5/stitch_svgs.py`.

See the [main README](../README.md#produce-figures-from-summary-statistics) for input files, Docker settings, and figure filenames.

Input data:

- `intermediate_data/` and `imputation_statistics/` contain test-region summary and imputation results.
- `intermediate_data_whole_genome/` and `imputation_statistics_whole_genome/` contain whole-genome results.
- `SHAPEIT5_switch_output/` contains switch-error calculations.

The test run writes figures and tables to `figures/` and `tables/`. The whole-genome run writes them to `figures_whole_genome/` and `tables_whole_genome/`.
