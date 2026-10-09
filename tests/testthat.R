library(testthat)
library(BGVAR)

# Runs every test-*.R file in tests/testthat, loading helper-*.R first.
# Current coverage: BGVAR inputs, Minnesota, SSVS, Normal-Gamma, Horseshoe,
# impulse responses (irf), FEVD, historical decomposition (hd), and prediction.
# New test files are included automatically; no individual source() calls needed.
test_check("BGVAR")
