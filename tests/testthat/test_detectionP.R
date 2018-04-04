context("detectionP")
test_that("test='detectionP' gives correct errors", {
    sset <- makeExampleSeSAMeDataSet()

    ## normal output is numeric
    ssetc <- sset
    expect_is(detectionPnegEcdf(ssetc), "SigSet")

    ## missing negative control, issue error
    ssetc@ctl <- ssetc@ctl[1:3,]
    expect_error(detectionPnegEcdf(ssetc))
})
