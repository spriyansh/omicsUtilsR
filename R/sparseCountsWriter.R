#' Convert 10x Counts to Desired Format
#'
#' This function converts 10x genomics sparse count matrices to either
#' H5 or MM format and saves to the specified path.
#'
#' @param sparseCounts A sparse matrix of 10x genomics counts.
#' @param format A character string specifying the desired output format. Can be either "h5" or "mm".
#' @param file.name The name of the output file.
#' @param inputPath The directory path where the output file should be saved.
#' @return This function does not return a value. It saves the converted counts to the specified path.
#' @importFrom rhdf5 h5createFile h5createGroup h5write
#' @importFrom DropletUtils write10xCounts
#' @export
convert10xCounts <- function(sparseCounts, format = "h5", file.name, inputPath) {

    file.name <- paste0(inputPath, file.name)

    if (format == "h5") {
        if (file.exists(file.name)) {
            unlink(file.name)
            h5createFile(file.name)
        } else {
            h5createFile(file.name)
        }

        h5createGroup(file.name, group = "matrix")
        h5createGroup(file.name, group = "matrix/features")

        suppressMessages(h5write(sparseCounts@Dimnames[2][[1]],
                                 file.name,
                                 name = "matrix/barcodes",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))

        suppressMessages(h5write(sparseCounts@x,
                                 file.name,
                                 name = "matrix/data",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        suppressMessages(h5write(sparseCounts@i,
                                 file.name,
                                 name = "matrix/indices",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        suppressMessages(h5write(sparseCounts@p,
                                 file.name,
                                 name = "matrix/indptr",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        suppressMessages(h5write(c(nrow(sparseCounts), ncol(sparseCounts)),
                                 file.name,
                                 name = "matrix/shape",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        suppressMessages(h5write(sparseCounts@Dimnames[1][[1]],
                                 file.name,
                                 name = "matrix/features/id",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        suppressMessages(h5write("genome",
                                 file.name,
                                 name = "matrix/features/_all_tag_keys",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))

        suppressMessages(h5write(rep("Gene Expression", length(sparseCounts@Dimnames[1][[1]])),
                                 file.name,
                                 name = "matrix/features/feature_type",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))

        suppressMessages(h5write(rep("cr_index_hs", length(sparseCounts@Dimnames[1][[1]])),
                                 file.name,
                                 name = "matrix/features/genome",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))

        suppressMessages(h5write(sparseCounts@Dimnames[1][[1]],
                                 file.name,
                                 name = "matrix/features/name",
                                 write.attributes = T,
                                 createnewfile = F, native = TRUE
        ))
        message(paste("Written to", file.name, "Successfully"))
    } else if (format == "mm") {

        write10xCounts(x = sparseCounts, path = file.name)
        message(paste("Written to", file.name, "Successfully"))

    } else {
        stop("Supported formats are H5 and MM")
    }
}
