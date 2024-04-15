## Version 0.5.5: Now that I can handle them, I want more lineage definitions!

- `usher_barcodes()` can now download `usher_barcodes_with_gisaid.csv`.
- `covariants()` runs a python script to get lineage definitions from [Hodcroft Lab](https://github.com/hodcroftlab/covariants), the data behind [covariants.org](covariants.org).
    - Requires the `reticulate` package and an installation of python (the lineage definitions are stored in a python file, and I'm not a good enough programmer to parse it).
- `plot_lineage_defs2()` can handle cases where mutations/lineages are completely missing from one of the defs. 
- `pkgdown()` site complete.

## Version 0.5.4: What will I do with all these lineage definitions?

- Visualize which mutations were used the most with `plot_actual_defs(res, type = "used")`!
- Visualize the coverage for each mutation used with `plot_actual_defs(res, type = "coverage")`!
- Colours are always going to be a pain, but `plot_lineage_defs2()` has been overhauled to be a little easier.
- Bugfixes:
    - error in plot and autoplot when `by_col` was NULL.

## Version 0.5.3: What does that lineage definition look like to you?

- Added a plotting function `plot_lineage_defs()` for a single lineage definition matrix. Made it pretty.
- Added a plotting function `plot_lineage_defs2()` that compares two lineage definition matrices. Made it pretty.
    - Allows for user-specified colours, which was hard and probably won't be used much but I'm happy with it.
    - Accepts output of `provoc()`, from which it will extract the actual lineage defintions used.

## Version 0.5.2: Call it what you want.

- Documentation has been significantly updated.
- BREAKING CHANGE: If you use the argument `by = "sra"`, the result will have a column labelled `sra` (not `group`, as before).

## Version 0.5.1: "Lineages of Concern" doesn't have the same ring to it.

- BREAKING CHANGE: All names changed from "variant" to "lineage" (except in "variant of concern").
    - `varmat` has been changed to `lineage_defs` because that's much clearer.
- Minor documentation updates, mainly failed inheritance.
