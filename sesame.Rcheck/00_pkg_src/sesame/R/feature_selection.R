

getSignatureU <- function(
    betas, grouping, u_max = 0.2, m_min = 0.7,
    max_na_in = 0, max_na_out = 0) {
    
    groups <- unique(grouping)
    is_na <- is.na(betas)
    sigs <- lapply(groups, function(g) {
        m1 <- rowMeans(betas[,grouping==g], na.rm=TRUE) < u_max
        m2 <- rowMeans(betas[,grouping!=g], na.rm=TRUE) > m_min
        ps1 <- rowSums(is_na[,grouping==g]) <= max_na_in
        ps2 <- rowSums(is_na[,grouping!=g]) <= max_na_out
        names(which(m1 & m2 & ps1 & ps2)) })
    names(sigs) <- groups
    sigs
}

getSignatureUTop <- function(
    betas, grouping, n=100,
    max_na_in = 0, max_na_out = 0) {
    
    groups <- unique(grouping)
    is_na <- is.na(betas)
    sigs <- lapply(groups, function(g) {
        mean1 <- rowMeans(betas[,grouping == g], na.rm=TRUE)
        mean0 <- rowMeans(betas[,grouping != g], na.rm=TRUE)
        ps1 <- rowSums(is_na[,grouping == g]) <= max_na_in
        ps2 <- rowSums(is_na[,grouping != g] <= max_na_out)
        head(names(sort((mean1 - mean0)[ps1 & ps2])), n=n)
    })
    names(sigs) <- groups
    sigs
}

clusterWithSignature <- function(betas, grouping, sigs) {
    pbs <- do.call(c, lapply(names(sigs), function(g) {
        if (length(sigs[[g]]) > 5)
            rownames(row.cluster(betas[intersect(
                rownames(betas), sigs[[g]]),])$mat)
        else
            NULL
    }))
    spl <- do.call(c, lapply(names(sigs), function(g) {
        colnames(column.cluster(betas[,grouping == g])$mat)
    }))
    betas[pbs, spl]
}

clusterWithSampleGrouping <- function(
    betas, grouping, groups=unique(grouping)) {
    
    do.call(cbind, lapply(groups, function(g) {
        column.cluster(betas[,grouping == g])$mat
    }))
}

clusterWithinRowGroups <- function(betas, sigs) {
    do.call(rbind, lapply(sigs, function(x) {
        row.cluster(betas[x,])$mat
    }))
}
