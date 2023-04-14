#' List supported prepSesame functions
#'
#' @return a data frame with code, func, description
#' @examples
#' prepSesameList()
#' @export
prepSesameList <- function() {
    x <- data.frame(rbind(
        c("0", "resetMask", "Reset mask to all FALSE"),
        c("Q", "qualityMask", "Mask probes of poor design"),
        c("G", "prefixMaskButCG", "Mask all but cg- probes"),
        c("H", "prefixMaskButC", "Mask all but cg- and ch-probes"),
        c("C", "inferInfiniumIChannel", "Infer channel for Infinium-I probes"),
        c("D", "dyeBiasNL", "Dye bias correction (non-linear)"),
        c("E", "dyeBiasL", "Dye bias correction (linear)"),
        c("P", "pOOBAH", "Detection p-value masking using oob"),
        c("I", "ELBAR", "Mask background-dominated readings"),
        c("B", "noob", "Background subtraction using oob"),
        c("S", "inferSpecies", "Set species-specific mask"),
        c("T", "inferStrain", "Set strain-specific mask (mouse)"),
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
#' @param prep code that indicates preprocessing functions and their
#' execution order (functions on the left is executed first).
#' @param prep_args optional argument list to individual functions, e.g.,
#' prepSesame(sdf, prep_args=list(Q=list(mask_names = "design_issue")))
#' sets qualityMask(sdf, mask_names = "design_issue")
#' @return SigDF
#' @examples
#' sdf <- sesameDataGet("MM285.1.SigDF")
#' sdf1 <- prepSesame(sdf, "QCDPB")
#' @export
prepSesame <- function(sdf, prep = "QCDPB", prep_args = NULL) {
    cfuns <- prepSesameList()
    
    codes <- str_split(prep,"")[[1]]
    stopifnot(all(codes %in% cfuns$code))
    x <- sdf
    for(c1 in codes) {
        x <- do.call(get(cfuns[cfuns$code == c1, "func"]),
            c(list(x), prep_args[[c1]]))
    }
    x
}

wrap_openSesame1 <- function(func, ret, ...) {
    if (is.null(func)) {
        ret
    } else {
        func(ret, ...)
    }
}

wrap_openSesame <- function(x, ret) {
    if (all(vapply(ret, is.numeric, logical(1))) &&
            length(unique(vapply(ret, length, integer(1)))) == 1) {
        ## getBetas, getAFs, ...
        ret <- do.call(cbind, ret)
        if (is.null(colnames(ret)) &&
                is.character(x) && length(x) == ncol(ret)) {
            colnames(ret) <- basename(x)
        }
        ret
    } else { # others
        if (is.null(names(ret)) &&
                is.character(x) && length(x) == length(ret)) {
            names(ret) <- basename(x)
        }
        ret

    }
}

#' The openSesame pipeline
#'
#' This function is a simple wrapper of noob + nonlinear dye bias 
#' correction + pOOBAH masking.
#'
#' Please use mask=FALSE to turn off masking.
#' 
#' If the input is an IDAT prefix or a \code{SigDF}, the output is
#' the beta value numerics.
#' 
#' @param x SigDF(s), IDAT prefix(es)
#' @param prep preprocessing code, see ?prepSesame
#' @param prep_args optional preprocessing argument list, see ?prepSesame
#' @param manifest optional dynamic manifest
#' @param func either getBetas or getAFs, if NULL, then return SigDF list
#' @param platform optional platform string
#' @param ... parameters to getBetas
#' @param BPPARAM get parallel with MulticoreParam(n)
#' @return a numeric vector for processed beta values
#' @import BiocParallel
#' @examples
#'
#' in_dir <- system.file("extdata", "", package = "sesameData")
#' betas <- openSesame(in_dir)
#' ## or
#' IDATprefixes <- searchIDATprefixes(in_dir)
#' betas <- openSesame(IDATprefixes)
#'
#' @export
openSesame <- function(
    x, prep = "QCDPB", prep_args = NULL, manifest = NULL,
    func = getBetas, BPPARAM=SerialParam(), platform = "", ...) {

    ## expand if a directory
    if (length(x) == 1 && is(x, 'character') && dir.exists(x)) {
        x <- searchIDATprefixes(x)
    }

    if (is(x, "SigDF")) {
        wrap_openSesame1(func, prepSesame(x, prep, prep_args), ...)
    } else if (is(x, 'character')) {
        if (length(x) == 1) {
            wrap_openSesame1(func, prepSesame(readIDATpair(
                x, platform = platform, manifest = manifest),
                prep, prep_args), ...)
        } else { # multiple IDAT prefixes / SigDFs
            wrap_openSesame(x, bplapply(x, openSesame,
                platform = platform, prep = prep, prep_args = prep_args,
                func = func, manifest = manifest, BPPARAM=BPPARAM, ...))
        }
    } else if (is(x, "list") && is(x[[1]], "SigDF")) {
        wrap_openSesame(x, bplapply(x, openSesame,
            platform = platform, prep = prep, prep_args = prep_args,
            fun = func, manifest = manifest, BPPARAM=BPPARAM, ...))
    } else {
        stop("Unsupported input")
    }
}

