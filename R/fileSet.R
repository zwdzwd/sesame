
#' openSesame pipeline with file-backed storage
#'
#' @param map_path path of file to be mapped (beta values file)
#' @param idat_dir source IDAT directory
#' @param BPPARAM get parallel with MulticoreParam(2)
#' @param inc bytes per item data storage. increase to 8 if precision
#' is important. Most cases 32-bit representation is enough.
#' @return a sesame::fileSet
#' @import BiocParallel
#' @examples
#'
#' openSesameToFile('mybetas',
#'     system.file('extdata',package='sesameData'))
#' 
#' @export
openSesameToFile <- function(
    map_path, idat_dir, BPPARAM=SerialParam(), inc = 4) {
    samples <- basename(searchIDATprefixes(idat_dir))
    sdf <- readIDATpair(file.path(idat_dir, samples[1]))

    fset <- initFileSet(map_path, sdfPlatform(sdf), samples, inc = inc)
    
    message(
        'Mapping ', length(samples), ' ', sdfPlatform(sdf),
        ' samples to ', map_path, '.')
    
    returned <- bplapply(samples, function(sample) {
        try({
            ## message('.', appendLF=FALSE)
            betas <- openSesame(file.path(idat_dir, sample))
            mapFileSet(fset, sample, betas)
            TRUE
        })
    }, BPPARAM = BPPARAM)

    succeeded <- !vapply(
        returned, function(x) inherits(x, "try-error"), logical(1))

    message(
        'Successfully processed ', sum(succeeded),
        ' IDATs (', sum(!succeeded),' failed).')
    fset
}

#' initialize a fileSet class by allocating appropriate storage
#'
#' @param map_path path of file to map
#' @param platform EPIC, HM450 or HM27, consistent with sdfPlatform(sdf)
#' @param samples sample names
#' @param probes probe names
#' @param inc bytes per unit data storage
#' @return a sesame::fileSet object
#' @examples
#'
#' fset <- initFileSet('mybetas2', 'HM27', c('s1','s2'))
#' @export
initFileSet <- function(map_path, platform, samples,
    probes = NULL, inc = 4) {

    if (is.null(probes)) {
        addr <- sesameDataGet(paste0(platform,'.address'))
        probes <- addr$ordering$Probe_ID
    }
    
    fset <- structure(list(
        map_path = map_path,
        platform = platform,
        probes = probes,
        samples = sort(samples),
        inc = inc # bytes per storage
    ), class='fileSet')
    fset$n <- length(fset$probes)
    fset$m <- length(fset$samples)

    message(
        'Allocating space for ', fset$m, ' ',
        platform, ' samples at ', map_path, '.')
    
    ## initialize to NA, this can be time-consuming
    ## but this should make parallel writing safer
    fh <- file(map_path, 'wb')
    for (i in seq_along(samples)) {
        writeBin(
            as.numeric(rep(NA, times = fset$n)),
            fh, size = inc)
    }
    close(fh)
    saveRDS(fset, paste0(map_path,'_idx.rds'))
    fset
}

#' Read an existing fileSet from storage
#'
#' This function only reads the meta-data.
#' @param map_path path of file to map (should contain valid _idx.rds index)
#' @return a sesame::fileSet object
#' @examples
#'
#' ## create two samples
#' fset <- initFileSet('mybetas2', 'HM27', c('s1','s2'))
#'
#' ## a hypothetical numeric array (can be beta values, intensities etc)
#' hypothetical <- setNames(runif(fset$n), fset$probes)
#'
#' ## map the numeric to file
#' mapFileSet(fset, 's1', hypothetical)
#'
#' ## read it from file
#' fset <- readFileSet('mybetas2')
#'
#' ## get data
#' sliceFileSet(fset, 's1', 'cg00000292')
#' 
#' @export
readFileSet <- function(map_path) {
    fset <- readRDS(paste0(map_path,'_idx.rds'))
    fset$map_path <- map_path # just in case
    fset
}

#' Deposit data of one sample to a fileSet (and hence to file)
#'
#' @param fset a sesame::fileSet, as obtained via readFileSet
#' @param sample sample name as a string
#' @param named_values value vector named by probes
#' @return a sesame::fileSet
#' @examples
#'
#' ## create two samples
#' fset <- initFileSet('mybetas2', 'HM27', c('s1','s2'))
#'
#' ## a hypothetical numeric array (can be beta values, intensities etc)
#' hypothetical <- setNames(runif(fset$n), fset$probes)
#'
#' ## map the numeric to file
#' mapFileSet(fset, 's1', hypothetical)
#'
#' ## get data
#' sliceFileSet(fset, 's1', 'cg00000292')
#' 
#' @export
mapFileSet <- function(fset, sample, named_values) {

    con <- file(fset$map_path, 'r+b')
    sample_index <- which(fset$samples == sample)
    if (length(sample_index) != 1) {
        stop('There are no or more than one sample ', sample, ' in fileSet\n')
    }
    binWriteNumeric(
        con,
        named_values[match(fset$probes, names(named_values))],
        beg = (sample_index-1) * fset$n,
        inc = fset$inc)
    close(con)
    fset
}

#' Slice a fileSet with samples and probes
#'
#' @param fset a sesame::fileSet, as obtained via readFileSet
#' @param samples samples to query (default to all samples)
#' @param probes probes to query (default to all probes)
#' @param memmax maximum items to read from file to memory, to protect from
#' accidental memory congestion.
#' @return a numeric matrix of length(samples) columns and length(probes) rows
#' @examples
#'
#' ## create two samples
#' fset <- initFileSet('mybetas2', 'HM27', c('s1','s2'))
#'
#' ## a hypothetical numeric array (can be beta values, intensities etc)
#' hypothetical <- setNames(runif(fset$n), fset$probes)
#'
#' ## map the numeric to file
#' mapFileSet(fset, 's1', hypothetical)
#'
#' ## get data
#' sliceFileSet(fset, 's1', 'cg00000292')
#' 
#' @export
sliceFileSet <- function(
    fset, samples = fset$samples,
    probes = fset$probes, memmax = 10^5) {

    sample_indices <- match(samples, fset$samples)
    if (any(is.na(sample_indices))) {
        message(
            sum(is.na(sample_indices)),
            ' sample(s) are nonexistent')
        samples <- samples[!is.na(sample_indices)]
        sample_indices <- sample_indices[!is.na(sample_indices)]
    }

    probe_indices <- match(probes, fset$probes)
    if (any(is.na(probe_indices))) {
        message(
            sum(is.na(probe_indices)),
            ' probe(s) are nonexistent')
        probes <- probes[!is.na(probe_indices)]
        probe_indices <- probe_indices[!is.na(probe_indices)]
    }

    if (length(probes) * length(samples) > memmax) {
        stop('Too many items retrieved (memmax = ', memmax, '\n')
    }

    con <- file(fset$map_path, 'rb')
    res <- do.call(cbind, lapply(setNames(
        sample_indices, samples), function(s_ind) {
        vapply(probe_indices, function(p_ind) {
            binReadNumeric(
                con, s_ind, p_ind, fset$n, inc=fset$inc)
        }, numeric(1))
    }))
    close(con)
    rownames(res) <- probes
    res
}



#' Print a fileSet
#'
#' @param x a sesame::fileSet
#' @param ... stuff for print
#' @return string representation
#' @examples
#'
#' fset <- initFileSet('mybetas2', 'HM27', c('s1','s2'))
#' fset
#' 
#' @export
print.fileSet <- function(x, ...) {
    message('File Set for', x$n,
        'probes and', x$m, 'samples.\n')
}

binWriteNumeric <- function(con, num_array, inc=4, beg=0) {
    seek(con, beg * inc, origin='start', rw = 'write')
    writeBin(as.numeric(num_array), con, size=inc)
}

binReadNumeric <- function(con, s_ind, p_ind, p_len, n=1, inc=4) {
    seek(con, ((s_ind - 1) * p_len + p_ind - 1) * inc, origin='start')
    readBin(con, 'numeric', n, size=inc)
}
