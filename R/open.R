#' List supported prepSesame functions
#'
#' @return a data frame with code, func, description
#' @examples
#' prepSesameList()
#' @export
prepSesameList <- function() {
    x <- data.frame(rbind(
        c("Q", "qualityMask", "Mask probes of poor design"),
        c("G", "prefixMaskButCG", "Mask all but cg- probes"),
        c("H", "prefixMaskButC", "Mask all but cg- and ch-probes"),
        c("C", "inferInfiniumIChannel", "Infer channel for Infinium-I probes"),
        c("D", "dyeBiasNL", "Dye bias correction (non-linear)"),
        c("P", "pOOBAH", "Detection p-value masking using oob"),
        c("B", "noob", "Background subtraction using oob"),
        c("S", "inferSpecies", "Set species-specific mask"),
        c("T", "inferStrain", "Set strain-specific mask"),
        c("M", "matchDesign", "Match Inf-I/II in beta distribution")))
    colnames(x) <- c("code", "func", "description")
    x
}

#' Apply a chain of sesame preprocessing functions in an arbitrary order
#' 
#' Notes on the order of operation:
#' 1. qualityMask and inferSpecies should go before noob and pOOBAH,
#' otherwise the background is too high because of Multi,
#' uk and other probes
#' 2. dyeBias correction needs to happen early
#' 3. channel inference before dyebias
#' 4. noob should happen last, pOOBAH before noob because noob modifies oob
#' 
#' @param sdf SigDF
#' @param code code that indicates preprocessing functions and their
#' execution order (functions on the left is executed first).
#' @return SigDF
#' @examples
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' sdf1 <- prepSesame(sdf, "QCDPB")
#' @export
prepSesame <- function(sdf, code = "QCDPB") {
    cfuns <- prepSesameList()
    
    codes <- str_split(code,"")[[1]]
    stopifnot(all(codes %in% cfuns$code))
    x <- sdf
    for(c1 in codes) {
        x <- get(cfuns[cfuns$code == c1, "func"])(x) }
    x
}

#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#' 
#' If the input is an IDAT prefix or a \code{SigDF}, the output is
#' the beta value numerics.
#' 
#' @param x SigDF(s), IDAT prefix(es)
#' @param platform optional platform string
#' @param prep preprocessing code, see ?prepSesame
#' @param manifest optional dynamic manifest
#' @param func either getBetas or getAFs
#' @param ... parameters to getBetas
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @return a numeric vector for processed beta values
#' @import BiocParallel
#' @examples
#' sdf <- sesameDataGet('HM450.1.TCGA.PAAD')$sdf
#' IDATprefixes <- searchIDATprefixes(
#'     system.file("extdata", "", package = "sesameData"))
#' betas <- openSesame(IDATprefixes)
#' @export
openSesame <- function(
    x, platform = '', prep = "QCDPB", manifest = NULL,
    func = getBetas, BPPARAM=SerialParam(), ...) {

    ## expand if a directory
    if (length(x) == 1 && is(x, 'character') && dir.exists(x)) {
        x <- searchIDATprefixes(x)
    }

    if (is(x, "SigDF")) {
        func(prepSesame(x, prep), ...)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            func(prepSesame(readIDATpair(
                x, platform = platform, manifest = manifest), prep), ...)
        } else { # multiple IDAT prefixes / SigDFs
            do.call(
                cbind, bplapply(x, openSesame,
                    platform = platform, fun = func,
                    manifest = manifest, BPPARAM=BPPARAM, ...))
        }
    } else if (is(x, "list") && is(x[[1]], "SigDF")) {
        do.call(
            cbind, bplapply(x, openSesame,
                platform = platform, fun = func,
                manifest = manifest, BPPARAM=BPPARAM, ...))
    } else {
        stop("Unsupported input")
    }
}

