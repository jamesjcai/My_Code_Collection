#' Basic methods for class demand
#' This document lists a series of basic methods for the class demand

#' @param object Object of class demand
#' @return returns summary information about the demand object
#' @rdname printDeMAND
#' @examples
#' data(inputExample)
#' dobj <- demandClass(exp=bcellExp, anno=bcellAnno, network=bcellNetwork)
#' printDeMAND(dobj)
#' @export
printDeMAND <- function(x) {
        message("An object of class demand")
        message("Slot exp:")
        message("\tExpression data of ", nrow(x@exp), " features by ", ncol(x@exp), " samples")
        message("Slot anno:")
        message("\tAnnotation of ", nrow(x@anno), " probes and ", length(unique(x@anno[,2])), " genes")
        message("Slot network:")
        message("\tInteractome with ", length(unique(x@network)), " nodes, ", nrow(x@anno), " edges")
        message("Slot moa (head):")
        if (length(x@moa)>0) print(head(x@moa), row.names=FALSE)
        else message("\tEmpty")
        message("Slot KLD (head):")
        if (dim(x@KLD)[1]>0) print(head(x@KLD), row.names=FALSE) 
        else message("\tEmpty")
}

