
## 0.5.3: What does that lineage definition look like to you?

- Added a plotting function `plot_lineage_defs()` for a single lineage definition matrix. Made it pretty.
- Added a plotting function `plot_lineage_defs2()` that compares two lineage definition matrices. Made it pretty.
    - Allows for user-specified colours, which was hard and probably won't be used much but I'm happy with it.

## 0.5.2: Call it what you want.

- Documentation has been significantly updated.
- BREAKING CHANGE: If you use the argument `by = "sra"`, the result will have a column labelled `sra` (not `group`, as before).

## 0.5.1: "Lineages of Concern" doesn't have the same ring to it.

- BREAKING CHANGE: All names changed from "variant" to "lineage" (except in "variant of concern").
    - `varmat` has been changed to `lineage_defs` because that's much clearer.
- Minor documentation updates, mainly failed inheritance.
