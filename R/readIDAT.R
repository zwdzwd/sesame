## the following copied from illuminaio package but with additional fixes to support de-identified files.
readIDAT <- function(file) {
    stopifnot(is.character(file) || length(file) != 0)
    file <- path.expand(file)
    if(!file.exists(file)) {
        stop("File not found: ", file)
    }

    if(grepl("\\.gz", file))
        con <- gzfile(file, "rb")
    else
        con <- file(file, "rb")
    on.exit({ close(con) })

    ## Assert file format
    magic <- readChar(con, nchars=4)
    if (magic != "IDAT") {
        stop("Invalid IDAT magic number: ", magic, " (expected 'IDAT')")
    }

    ## Read IDAT file format version
    version <- readBin(con, what="integer", size=4, n=1, signed = TRUE, endian = "little")
    if (version == 3) {
        res <- readIDAT_nonenc(file)
    } else {
        stop("Unsupported IDAT version: ", version)
    }
    res
}

readIDAT_nonenc <- function(file, what = c("all", "IlluminaID", "nSNPsRead")) {
    readByte <- function(con, n=1, ...) {
        readBin(con, what="integer", n=n, size=1, endian="little", signed=FALSE)
    }

    readShort <- function(con, n=1, ...) {
        readBin(con, what="integer", n=n, size=2, endian="little", signed=FALSE)
    }

    readInt <- function(con, n=1, ...) {
        readBin(con, what="integer", n=n, size=4, endian="little", signed=TRUE)
    }

    readLong <- function(con, n=1, ...) {
        readBin(con, what="integer", n=n, size=8, endian="little", signed=TRUE)
    }

    readString <- function(con, ...) {
        ## From [1] https://code.google.com/p/glu-genetics/source/browse/glu/lib/illumina.py#86:
        ## String data are encoded as a sequence of one or more length
        ## bytes followed by the specified number of data bytes.
        ##
        ## The lower 7 bits of each length byte encodes the bits that
        ## comprise the length of the following byte string.  When the
        ## most significant bit it set, then an additional length byte
        ## follows with 7 additional high bits to be added to the current
        ## length.  The following string lengths are accommodated by
        ## increasing sequences of length bytes:
        ##
        ## length  maximum
        ## bytes   length
        ## ------  --------
        ##   1       127 B
        ##   2        16 KB
        ##   3         2 MB
        ##   4       256 MB
        ##   5        32 GB
        ##
        ## While this seems like a sensible progression, there is some
        ## uncertainty about this interpretation, since the longest of
        ## string observed in the wild has been of length 6,264 with
        ## two length bytes.
        ##
        ## EXAMPLES: by HB (2015-09-11)
        ## Length   Len bytes  Iterations (n, m, k, shift)
        ## 1        (1)
        ## 127      (127)       -> n=128
        ## 128      (128,1)     -> n=0  ,m=1  ,shift=7, k=128 -> n=128
        ## 255      (255,1)     -> n=128,m=1  ,shift=7, k=128 -> n=255
        ## 256      (128,2)     -> n=0  ,m=2  ,shift=7, k=256 -> n=256
        ## 257      (129,2)     -> n=1  ,m=2  ,shift=7, k=256 -> n=257
        ## 512      (128,4)     -> n=0  ,m=4  ,shift=7, k=512 -> n=512
        ## 81921    (129,128,5) -> n=1  ,m=128,shift=7, k=0 -> n=1+0=1
        ##                      -> n=0  ,m=5  ,shift=14,k=81920
        ##                                              -> n=1+81920=81921

        ## Parse the number of characters to read
        m <- readByte(con, n=1)
        n <- m %% 128

        shift <- 0L
        while (m %/% 128 == 1) {
            ## Read next length byte ...
            m <- readByte(con, n=1)

            ## ... which represents the next 7 hi-bits
            shift <- shift + 7L
            k <- (m %% 128) * 2^shift

            ## Total number of bytes to read
            n <- n + k
        }

        ## Now read all bytes/characters
        readChar(con, nchars=n, useBytes=TRUE)
    }

    readField <- function(con, field) {
        switch(field,
            "IlluminaID" = readInt(con = con, n=nSNPsRead),
            "SD" = readShort(con = con, n=nSNPsRead),
            "Mean"= readShort(con = con, n=nSNPsRead),
            "NBeads" = readByte(con = con, n=nSNPsRead),
            "MidBlock" = {
                nMidBlockEntries <- readInt(con = con, n=1)
                MidBlock <- readInt(con = con, n=nMidBlockEntries)
            },
            "RedGreen" = readInt(con = con, n=1),
            "MostlyNull" = readString(con = con),
            "Barcode" = readString(con = con),
            "ChipType" = readString(con = con),
            "MostlyA" = readString(con = con),
            "Unknown.1" = readString(con = con),
            "Unknown.2" = readString(con = con),
            "Unknown.3" = readString(con = con),
            "Unknown.4" = readString(con = con),
            "Unknown.5" = readString(con = con),
            "Unknown.6" = readString(con = con),
            "Unknown.7" = readString(con = con),
            "RunInfo" = {
                ## make sure it's >= 0.
                nRunInfoBlocks <- max(0,readInt(con = con, n=1))
                naValue <- as.character(NA)
                RunInfo <- matrix(naValue, nrow=nRunInfoBlocks, ncol=5)
                colnames(RunInfo) <- c("RunTime", "BlockType", "BlockPars",
                    "BlockCode", "CodeVersion")
                for (ii in seq_len(nRunInfoBlocks)) {
                    for (jj in seq_len(5)) {
                        RunInfo[ii,jj] <- readString(con = con)
                    }
                }
                RunInfo
            },
            stop("readIDAT_nonenc: unknown field"))
    }

    if(! (is.character(file) || try(isOpen(file))))
        stop("'file' must be a character or an open, seekable connection")
    what <- match.arg(what)

    if(is.character(file)) {
        stopifnot(length(file) == 1)
        file <- path.expand(file)
        stopifnot(file.exists(file))
        fileSize <- file.info(file)$size
        if(grepl("\\.gz$", file))
            con <- gzfile(file, "rb")
        else
            con <- file(file, "rb")
        on.exit({
            close(con)
        })
    } else {
        con <- file
        fileSize <- 0
    }

    if(!isSeekable(con))
        stop("File connection must be seekable")

    ## Assert file format
    magic <- readChar(con, nchars=4)
    if (magic != "IDAT") {
        stop("Invalid IDAT magic number: ", magic, " (expected 'IDAT')")
    }

    ## Read IDAT file format version
    version <- readLong(con, n=1)
    if (version < 3) {
        stop("Unsupported IDAT version: ", version)
    }

    ## Number of fields
    nFields <- readInt(con, n=1)

    fields <- matrix(0, nrow=nFields, ncol=3)
    colnames(fields) <- c("fieldCode", "byteOffset", "Bytes")
    for (ii in seq_len(nFields)) {
        fields[ii,"fieldCode"] <- readShort(con, n=1)
        fields[ii,"byteOffset"] <- readLong(con, n=1)
    }

    knownCodes <- c(
        "nSNPsRead"  = 1000,
        "IlluminaID" =  102,
        "SD"         =  103,
        "Mean"       =  104,
        "NBeads"     =  107,
        "MidBlock"   =  200,
        "RunInfo"    =  300,
        "RedGreen"   =  400,
        "MostlyNull" =  401, # 'Manifest', cf [1].
        "Barcode"    =  402,
        "ChipType"   =  403,
        "MostlyA"    =  404, # 'Stripe', cf [1].
        "Unknown.1"  =  405,
        "Unknown.2"  =  406, # 'Sample ID', cf [1].
        "Unknown.3"  =  407,
        "Unknown.4"  =  408, # 'Plate', cf [1].
        "Unknown.5"  =  409, # 'Well', cf [1].
        "Unknown.6"  =  410,
        "Unknown.7"  =  510
    )

    nNewFields <- 1
    rownames(fields) <- paste("Null", seq_len(nFields))
    for (ii in seq_len(nFields)) {
        temp <- match(fields[ii,"fieldCode"], knownCodes)
        if (!is.na(temp)) {
            rownames(fields)[ii] <- names(knownCodes)[temp]
        } else {
            rownames(fields)[ii] <- paste("newField", nNewFields, sep=".")
            nNewFields <- nNewFields + 1
        }
    }

    stopifnot(min(fields[, "byteOffset"]) == fields["nSNPsRead", "byteOffset"])

    seek(con, where=fields["nSNPsRead", "byteOffset"], origin="start")
    nSNPsRead <- readInt(con, n=1)
    if(what == "nSNPsRead")
        return(nSNPsRead)
    if(what == "IlluminaID") {
        where <- fields["IlluminaID", "byteOffset"]
        seek(con, where = where, origin = "start")
        res <- readField(con = con, field = "IlluminaID")
        return(as.character(res))
    }

    res <- rownames(fields)
    names(res) <- res
    res <- res[order(fields[res, "byteOffset"])]
    res <- res[names(res) != "nSNPsRead"]
    ## only these fields will be searched.
    ## gzipped IDAT has some issue finding other fields.
    res <- res[c("IlluminaID", "SD", "Mean", "NBeads")]
    res <- lapply(res, function(xx) {
        where <- fields[xx, "byteOffset"]
        seek(con, where = where, origin = "start")
        readField(con = con, field = xx)
    })

    Unknowns <-
        list(MostlyNull=res$MostlyNull,
            MostlyA=res$MostlyA,
            Unknown.1=res$Unknown.1,
            Unknown.2=res$Unknown.2,
            Unknown.3=res$Unknown.3,
            Unknown.4=res$Unknown.4,
            Unknown.5=res$Unknown.5
        )

    Quants <- cbind(res$Mean, res$SD, res$NBeads)
    colnames(Quants) <- c("Mean", "SD", "NBeads")
    rownames(Quants) <- as.character(res$IlluminaID)

    res <- list(
        fileSize=fileSize,
        versionNumber=version,
        nFields=nFields,
        fields=fields,
        nSNPsRead=nSNPsRead,
        Quants=Quants,
        MidBlock=res$MidBlock,
        RedGreen=res$RedGreen,
        Barcode=res$Barcode,
        ChipType=res$ChipType,
        RunInfo=res$RunInfo,
        Unknowns=Unknowns
    )
    res
}
