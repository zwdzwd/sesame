context('sesamize')

test_that("RGChannelSetToSigSet gives correct results", {
    library(FlowSorted.Blood.450k)
    rgSet <- FlowSorted.Blood.450k[,1:3]
    ssets <- RGChannelSetToSigSet(rgSet, BPPARAM=MulticoreParam(3))
    expect_is(rgSet, "RGChannelSet")
    expect_equal(length(ssets), 3)

    ## now convert back to RGChannelSet to check
    rgSet2 <- SigSetToRGChannelSet(ssets)
    expect_equal(ncol(rgSet), ncol(rgSet2))
    expect_equal(annotation(rgSet)['array'], annotation(rgSet2)['array'])

    ## check old and new mSet is the same
    mSet <- preprocessRaw(rgSet)
    mSet2 <- preprocessRaw(rgSet2)
    expect_equal(colnames(mSet), colnames(mSet2))
    expect_equal(getMeth(mSet), getMeth(mSet2))
    expect_equal(getUnmeth(mSet), getUnmeth(mSet2))
    expect_equal(getBeta(mSet), getBeta(mSet2))
})

test_that("SigSetToRGChannelSet gives correct results", {
    library(minfi)
    sset <- makeExampleSeSAMeDataSet()
    rgSet <- SigSetToRGChannelSet(sset)
    expect_is(rgSet, "RGChannelSet")
    expect_equal(ncol(rgSet), 1)

    ## when represented as sset, the two are the same, given
    ## the rwo orders are the same
    sset2 <- RGChannelSetToSigSet(rgSet)[[1]]
    expect_equal(sset@IG, sset2@IG[rownames(sset@IG),])
    expect_equal(sset@IR, sset2@IR[rownames(sset@IR),])
    expect_equal(sset@II, sset2@II[rownames(sset@II),])
    expect_equal(sset@oobG, sset2@oobG[rownames(sset@oobG),])
    expect_equal(sset@oobR, sset2@oobR[rownames(sset@oobR),])
    expect_equal(sset@ctl, sset2@ctl[rownames(sset@ctl),])
})
