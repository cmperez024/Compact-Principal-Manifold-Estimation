# Reproduce `interval_S1_full.png`

Run the plotting script directly from any working directory, supplying its path:

```powershell
Rscript path/to/code/plot_1d_figures.R
```

The script includes its own data generation, persistence calculation, and
CompactPME fitting code. Missing inputs are generated automatically; it never
executes another simulation script. It writes `interval_S1_full.png` beside
`REPRODUCE.md` and stores its own inputs in `results_1d/` there.
`generate_1d_results.R` remains an optional way to precompute those same inputs.
Full calculations are expensive. `PME_SMOKE_TEST=1` runs a reduced validation
experiment; set `PME_OUTPUT_DIR` to a temporary directory when using this mode.
