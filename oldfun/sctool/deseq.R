# library( DESeq )

estimateSizeFactorsForMatrixx <- function(counts, locfunc=stats::median,
                                         geoMeans, controlGenes) {
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    loggeomeans <- rowMeans(log(counts))
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}


dtable<-read.csv('tmp.txt',header=FALSE)
sfHeLa <- estimateSizeFactorsForMatrixx( dtable )
nCountsHela <-t(t(dtable)/sfHeLa)

View(nCountsHela)
