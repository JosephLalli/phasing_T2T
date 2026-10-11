# Notebooks

This directory contains the notebooks used to calculate paper values and make figures.

`calc_figures_for_paper.ipynb` calculates all values reported in the text of the paper. It only uses smaller summary parquet files and is designed to run on a laptop in Docker.

`calc_per_variant_figures_for_paper.ipynb` calculates values that require loading the full variant dataframes into memory. Run it on a server in Docker.

The notebooks in this directory use the test-region data. `notebooks_whole_genome/` contains the whole-genome notebooks and figure code.

`notebooks_whole_genome/make_plots.py` is the Jupytext `py:percent` source for `make_plots.ipynb`. Edit the Python file, then sync and validate the notebook. Set `PHASING_T2T_RUN_PROFILE=whole_genome` for publication inputs or `test` for test-region inputs.

Run the host `scripts/run_make_plots.sh` wrapper after leaving the interactive container; it checks the 188 supplied input files and runs the whole-genome notebook in Docker. Run Figure 5 separately with `scripts/figure5/Figure_5_script.R` and `scripts/figure5/stitch_svgs.py`. `figure_composition.py` combines figure panels and `utilities.py` contains plotting and data helpers. Run notebook, Python, and R commands inside the Docker container described in the [main README](../README.md#run-downloads-and-manual-analysis-in-docker).

See the [main README](../README.md#publication-figure-code) for input files, Docker instructions, and figure filenames.

Input data:

- `intermediate_data/` and `imputation_statistics/` contain test-region summary and imputation results.
- `intermediate_data_whole_genome/` and `imputation_statistics_whole_genome/` contain whole-genome results.
- `SHAPEIT5_switch_output/` contains switch-error calculations.

The test run writes figures and tables to `figures/` and `tables/`. The whole-genome run writes them to `figures_whole_genome/` and `tables_whole_genome/`.
