context("impute")

test_that("Impute mean functions properly", {
    mx <- cbind(a = c(NA, 2, 3, 4), b = c(1, 6, 5, NA), c = c(3, 6, 7, 8))
    mx_imputed_cols <- imputeBetasMatrixByMean(mx, axis = 2)
    mx_imputed_rows <- imputeBetasMatrixByMean(mx, axis = 1)
    expect_true(mx_imputed_cols[1,1] == 3)
    expect_true(mx_imputed_cols[4,2] == 4)
    expect_true(mx_imputed_rows[1,1] == 2)
    expect_true((mx_imputed_rows[4,2]) == 6)
})
