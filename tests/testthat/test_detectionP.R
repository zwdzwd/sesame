context("detectionP")
test_that("test='detectionP' gives correct errors", {
    sdf <- sesameDataGet("EPIC.1.SigDF")
    expect_is(pOOBAH(sdf), "SigDF")
})
