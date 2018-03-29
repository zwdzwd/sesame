context("background")
test_that("test=background subtraction gives correct warning", {
    set.seed(1)
    sset <- makeExampleTinyEPICDataSet()
    expect_warning(noobsb(sset))
})
