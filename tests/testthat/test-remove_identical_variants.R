test_that("ensuring only unique lineages are kept", {
  lineage_defs <- simulate_lineage_defs()
  coco <- simulate_coco(lineage_defs)
  lin_df <- as.data.frame(t(lineage_defs))
  names(lin_df) <- paste0("lin_", names(lin_df))
  lin_df$mutation <- rownames(lin_df)
  fused_df <- dplyr::left_join(coco, lin_df, by = "mutation")
  fused_df$lin_BA.3 <- as.vector(fused_df$lin_BA.1)
  fused_df$lin_BA.4 <- as.vector(fused_df$lin_BA.2)
  fused_df$lin_BA.5 <- as.vector(fused_df$lin_BA.2)
  updated_fused_df <- provoc:::remove_identical_lineages(fused_df, annihilate = TRUE)
  expect_equal(ncol(updated_fused_df), 6)
})

test_that("ensuring when annilihate is false, no columns are removed", {
  lineage_defs <- simulate_lineage_defs()
  coco <- simulate_coco(lineage_defs)
  lin_df <- as.data.frame(t(lineage_defs))
  names(lin_df) <- paste0("lin_", names(lin_df))
  lin_df$mutation <- rownames(lin_df)
  fused_df <- dplyr::left_join(coco, lin_df, by = "mutation")
  fused_df$lin_BA.3 <- as.vector(fused_df$lin_BA.1)
  fused_df$lin_BA.4 <- as.vector(fused_df$lin_BA.2)
  fused_df$lin_BA.5 <- as.vector(fused_df$lin_BA.2)
  updated_fused_df <- provoc:::remove_identical_lineages(fused_df, annihilate = FALSE)
  expect_equal(ncol(updated_fused_df), 9)
})