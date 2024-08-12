#' @title Extracting properties from genomic regions
#'
#' @rdname extractRegionProperty-methods
#' @docType methods
#'
#' @description Methods for annotating the query \link{GenomicRanges} with various properties on the subject regions,
#' and the subject regions could be either \link{GenomicRanges} or \link{GRangesList} objects.
#'
#' The \code{extractRegionProperty} function extract the annotations on each query ranges from the properties of their overlapping/nearest subject regions.
#'
#' Other functions are specific implementations of \code{extractRegionProperty},
#' and their purpose is to summarize the genomic features frequently involved in functional genomics analysis.
#'
#' \code{extractRegionOverlap} computes the overlapping of each query ranges with the subject regions.
#' Likewise, the \code{extractRegionLength}, \code{extractRegionLetterFrequency}, and \code{extractRegionScores} functions
#' compute the Lengths, Sequence Contents, and \link{GenomicScores} of the overlapped subjects.
#'
#' \code{extractRegionYCount} and \code{extractRegionNearestDistToY} quantify the density/distance of an annotation (Y)
#' over the subject regions, the resulting properties can capture the clustering effect between the subject region and the annotation object.
#'
#' The functions \code{extractDistToRegion5end} & \code{extractDistToRegion3end} retrieve the distance between the query and ends of the subject, and
#' the \code{extractRegionRelativePosition} function extract the relative position of the query in the matched regions.
#' Specifically, the relative position is calculated by dividing the distance of the query to the 5'-end of the region by the length of the region.
#'
#' @param x A \link{GRanges} object for the query.
#' @param region A \link{GRanges} or \link{GRangesList} object for the subject.
#' If not provided (\code{region=NULL}), the returned properties will be calculated directly from \code{x}.
#' @param property A vector specifying the properties/attributes of the \code{region}.
#'
#' @param ambiguityMethod By default, for \code{logical} \code{property} input, if x overlaps multiple regions,
#' as long as any of the properties is TRUE, the returned value will be TRUE.
#' For other types of \code{property} input, the returned property will be the average of the properties over multiple overlapped regions.
#' If \code{ambiguityMethod} is \code{"mean"}, \code{"sum"}, \code{"min"}, or \code{"max"},
#' then the mean, sum, minimum, and maximum values of the >1 mapping will be returned in the output value.
#'
#' @param maxgap,minoverlap,type See \code{?\link{findOverlaps}} in the \bold{IRanges} package for a description of these arguments.
#'
#' @param nomapValue When \code{nomapValue} is \code{"NA"}, \code{"zero"}, or \code{"FALSE"},
#' the \code{x} that do not match the \code{region} will return \code{NA}, \code{0}, and \code{FALSE} respectively.
#' If \code{nomapValue} is \code{"nearest"}, the not matched \code{x} will be set to be the property value on its nearest \code{region}.
#'
#' @param sequence A \link{BSgenome} or \link{XStringSet} object.
#' See the \link{available.genomes} function for how to install a genome.
#'
#' @param gscores A \link{GScores} object.
#' See the vignette for more details: \code{vignette("GenomicScores", package = "GenomicScores")}.
#'
#' @param as.prob,letters See \code{?\link{letterFrequency}} in the \bold{Biostrings} package for a description of these arguments.
#'
#' @param ignore.strand When set to \code{TRUE}, the strand information is ignored in the overlap calculations.
#'
#' @param output.logical If \code{TRUE} then the returned value is \code{logical}, otherwise \code{numeric}.
#'
#' @param efficient \code{TRUE} if only internally extract the properties on \code{regions} overlap with \code{x},
#' the option makes the computation more efficient if the number of gnomic regions is much larger than the query.
#'
#' @param missingScores Specifying the imputation methods on missing scores in \code{extractRegionScores},
#' default is returning zero for unknown gscores under the region.
#'
#' @param y A \link{GRanges} or \link{GRangesList} object for the clustering annotation.
#'
#' @param normalize If \code{TRUE}, then returned count in \code{extractRegionYCount} will be divided by the length of the region.
#' The normalized count can be understood as the density of the annotation on the area.
#'
#' @param maxDist A value of the maximum nucleotide distance returned in \code{extractRegionNearestDistToY}, default: \code{3e+06}.
#'
#' @param ... Additional arguments, it passes to \link{letterFrequency} in method \code{extractRegionLetterFrequency}, and it passes to \link{score} in method \code{extractRegionScores}.
#'
#' @details
#' For specific extractor functions, when only x is provided with \code{region = NULL},
#' the functions directly use the query ranges to calculate the corresponding properties.
#'
#' In case the region ranges are provided,
#' the extractor functions first find the index where x overlaps with the region and then matches the property of the region to the query GRanges according to the overlapping index.
#'
#' For \code{nomapValue = "nearest"},
#' when ranges in \code{x} do not overlap with ranges in \code{regions}, the extractor functions will additionally match the unmapped x with its nearest region.
#'
#' For \code{extractRegionNearestDistToY}, if ranges in \code{y} overlap the ranges in \code{x},
#' the overlapped ranges in \code{y} are skipped and the ranges precede or follow \code{x} are used to calculate the distances.
#'
#' @return
#' A logical or a numeric vector with the same length as \code{x}.
#'
#' For \code{extractRegionProperty}, the vector type depends on the type of \code{property}.
#'
#' For \code{extractRegionOverlap} either logical when \code{output.logical = TRUE} or a numeric otherwise.
#'
#' @author
#' Zhen Wei
#'
#' @seealso
#' \itemize{
#' \item{}{The \link{genomeDerivedFeatures} function for the extraction of genome-derived features.}
#' \item{}{The \link{sequenceDerivedFeatures} function for the extraction of sequence-derived features.}
#' \item{}{The \link{topologyOnTranscripts} function for generating the topology of annotations on transcripts of genes.}
#' }
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## BASIC EXAMPLES
#' ## ---------------------------------------------------------------------
#' library(GenomicRanges)
#'
#' ## Build the query GRanges object:
#' x_gr <- GRanges(rep(c("chr1", "chr2"), c(5, 15)),
#'                 IRanges(c(sample(11874:12127, 5), sample(38814:41527,15)), width=100),
#'                 strand=Rle(c("+", "-"), c(5, 15)))
#' x_gr
#'
#' ## The region GRanges or GRangesList object:
#' exons_gr <- GRanges(c("chr1","chr2","chr2"),
#'                     IRanges(start=c(11874,38814,45440),end=c(12227,41627,46588)),
#'                      strand=c("+","-","-"))
#' genes_grl <- GRangesList(gene1=exons_gr[1],gene2=exons_gr[c(2,3)])
#'
#' ## Extract lengths of the query:
#' extractRegionLength(x_gr)
#'
#' ## Extract length of the exon overlapped by the query:
#' extractRegionLength(x_gr, exons_gr)
#'
#' ## Extract exonic length of the genes overlapped by the query:
#' extractRegionLength(x_gr, genes_grl)
#'
#' ## Exract self defined property on exons overlapped by the query:
#' exons_property <- c(1,6,8)
#' extractRegionProperty(x_gr, exons_gr, exons_property)
#'
#' ## ---------------------------------------------------------------------
#' ## MORE FEATURE EXTRACTORS
#' ## ---------------------------------------------------------------------
#'
#' ## Quantifying clustering effect with annotation y:
#'
#' ## Self clustering
#' extractRegionYCount(x_gr, x_gr)
#'
#' ## Self clustering on overlapping exons/genes
#' extractRegionYCount(x_gr, x_gr, exons_gr)
#' extractRegionYCount(x_gr, x_gr, genes_grl)
#'
#' ## Clustering with another annotation:
#' y_gr <- GRanges(rep(c("chr1", "chr2"), c(50, 50)),
#'                 IRanges(c(sample(11874:12127, 50),
#'                           sample(38814:41527,50)), width=1),
#'                 strand=Rle(c("+", "-"), c(50, 50)))
#'
#' extractRegionYCount(x_gr, y_gr)
#' extractRegionYCount(x_gr, y_gr, exons_gr)
#' extractRegionYCount(x_gr, y_gr, genes_grl)
#' extractRegionNearestDistToY(x_gr, y_gr)
#'
#' ## Relative position on exons/genes:
#' extractRegionRelativePosition(x_gr, exons_gr)
#' extractRegionRelativePosition(x_gr, genes_grl)
#'
#' ## Distance to 5'end/3'end of exons/genes:
#' extractDistToRegion5end(x_gr, exons_gr)
#' extractDistToRegion5end(x_gr, genes_grl)
#' extractDistToRegion3end(x_gr, exons_gr)
#' extractDistToRegion3end(x_gr, genes_grl)
#'
#' ## Extract features that depends on sequence/annotation packages:
#'
#' \donttest{
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg19
#'
#' ## GC content of the query:
#' extractRegionLetterFrequency(x_gr, bsgenome, letters="GC")
#'
#' ## GC contents of exons/genes containing the query:
#' extractRegionLetterFrequency(x_gr, bsgenome, exons_gr, letters="GC")
#' extractRegionLetterFrequency(x_gr, bsgenome, genes_grl, letters="GC")
#'
#' library(phastCons100way.UCSC.hg19)
#' gscores <- phastCons100way.UCSC.hg19
#'
#' ## PhastCons scores of the query:
#' extractRegionScores(x_gr, gscores)
#'
#' ## PhastCons scores of exons/genes containing the query:
#' extractRegionScores(x_gr, gscores, exons_gr)
#' extractRegionScores(x_gr, gscores, genes_grl)
#'
#' }
#'

#' ## ---------------------------------------------------------------------
#' ## VISUALIZE FEATURE IN PREDICTION MODEL
#' ## ---------------------------------------------------------------------
#'
#' ## Load the GRanges of the m6A miCLIP dataset prepared for the classification model:
#' GSE63753_abcam <- readRDS(system.file("extdata", "GSE63753_abcam.rds", package = "predictiveFeatures"))
#'
#' ## The metadata column of the GRanges is a 0/1 vector, 1 for the positive m6A site, 0 is the negative DRACH:
#' GSE63753_abcam
#' table(GSE63753_abcam$target)
#'
#' ## Extract exon length overlapped by the DRACH sites:
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' GSE63753_abcam$exon_length <- extractRegionLength(GSE63753_abcam, exons(TxDb.Hsapiens.UCSC.hg19.knownGene))
#'
#' ## Plot the logistic regression fit with cubic splines:
#' library(ggplot2)
#'
#' ggplot(na.omit(as.data.frame(mcols(GSE63753_abcam))), aes(log(exon_length), target)) +
#'   geom_smooth(formula = y ~ splines::ns(x, 3), method = "glm", method.args = list(family = "binomial")) +
#'   geom_hline(yintercept = 0.5, linetype = 3) +
#'   scale_x_continuous(limits = c(4,9)) +
#'   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#'   theme_classic() + labs(x = "log(length of exon)", y = "prob of m6A = 1", title = "LR fit with cubic splines")
#'
#' @importFrom BSgenome getSeq
#' @importFrom GenomicRanges findOverlaps countOverlaps distanceToNearest split resize nearest precede follow distance trim
#' @importFrom BiocGenerics width start end unlist unique strand
#' @importFrom Biostrings letterFrequency
#' @importFrom methods is
#' @importFrom GenomicScores score
#' @importFrom S4Vectors elementNROWS runValue queryHits subjectHits
#' @importFrom GenomicFeatures mapToTranscripts
#' @importFrom stats na.omit
#' @aliases extractRegionProperty
#' @export
setMethod("extractRegionProperty",
          "GRanges",
          function(x,
                   region,
                   property,
                   ambiguityMethod = c("auto", "mean", "sum", "min", "max"),
                   maxgap = 0L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "FALSE", "nearest"),
                   ignore.strand = FALSE) {
            stopifnot(is(x, "GRanges"))
            stopifnot(is(region, "GRanges") | is(region, "GRangesList"))
            stopifnot(length(property) == length(region))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            return_property <- rep(NA, length(x))
            
            fol <-
              findOverlaps(x,
                           region,
                           maxgap = maxgap,
                           minoverlap = minoverlap,
                           type = type)
            property_mapped <- property[subjectHits(fol)]
            
            if (!is.null(property_mapped)) {
              if (ambiguityMethod == "auto") {
                if (is.logical(property)) {
                  weighted_properties <- tapply(property_mapped, queryHits(fol), any)
                } else{
                  weighted_properties <- tapply(property_mapped, queryHits(fol), mean)
                }
              } else{
                weighted_properties <-
                  tapply(property_mapped, queryHits(fol), eval(parse(text = ambiguityMethod)))
              }
              return_property[as.numeric(names(weighted_properties))] <-
                weighted_properties
            }
            
            if (anyNA(return_property)) {
              if (nomapValue == "nearest") {
                if (!is.null(property_mapped)) {
                  if (is(region, "GRanges")) {
                    nomap_indx <- which(is.na(return_property))
                    return_property[nomap_indx] <-
                      property[nearest(x[nomap_indx], region)]
                    return_property[is.na(return_property)] <- 0
                  } else{
                    region_ranges <- range_unique_grl(region)
                    nomap_indx <- which(is.na(return_property))
                    return_property[nomap_indx] <-
                      property[nearest(x[nomap_indx], unlist(region_ranges))]
                    rm(region_ranges, nomap_indx)
                    return_property[is.na(return_property)] <- 0
                  }
                }
              } else if (nomapValue == "zero") {
                return_property[is.na(return_property)] <- 0
              } else if (nomapValue == "FALSE") {
                return_property[is.na(return_property)] <- FALSE
              }
            }
            return(return_property)
          })

range_unique_grl <- function(x) {
  names(x) <- seq_along(x)
  x_ranges <- range(x)
  indx_multi <- elementNROWS(x_ranges) > 1
  range_multi <- x_ranges[indx_multi]
  range_multi_unlist <- unlist(range_multi)
  rm(range_multi)
  range_dedup <-
    range_multi_unlist[!duplicated(names(range_multi_unlist))]
  rm(range_multi_unlist)
  x_ranges[indx_multi] <-
    split(range_dedup, names(range_dedup))[names(range_dedup)]
  rm(indx_multi, range_dedup)
  return(x_ranges)
}

#' @aliases extractRegionOverlap
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
setMethod("extractRegionOverlap",
          "GRanges",
          function(x,
                   region = NULL,
                   ambiguityMethod = c("auto", "mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   ignore.strand = FALSE,
                   output.logical = TRUE) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            
            if (is.null(region)) {
              return(rep(TRUE, length(x)))
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = rep(TRUE, length(region)),
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = "FALSE",
                ignore.strand = ignore.strand
              )
            }
            
            if (!output.logical) {
              region_property <- as.numeric(region_property)
            }
            
            return(region_property)
          })

#' @aliases extractRegionLength
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
setMethod("extractRegionLength",
          "GRanges",
          function(x,
                   region = NULL,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest"),
                   ignore.strand = FALSE) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            nomapValue <- match.arg(nomapValue)
            type <- match.arg(type)
            
            if (is.null(region)) {
              length_property <- width(x)
            } else if (is(region, "GRanges")) {
              length_property <- width(region)
            } else if (is(region, "GRangesList")) {
              region <- region[elementNROWS(region) != 0]
              length_property <- sum(width(region))
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            
            if (is.null(region)) {
              return(length_property)
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = length_property,
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = nomapValue,
                ignore.strand = ignore.strand
              )
              return(region_property)
            }
          })

#' @aliases extractRegionLetterFrequency
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
#sequence can be a BSgenome or a XStringSet object
setMethod("extractRegionLetterFrequency",
          "GRanges",
          function(x,
                   sequence,
                   region = NULL,
                   letters = "GC",
                   as.prob = TRUE,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest"),
                   ignore.strand = FALSE,
                   efficient = TRUE,
                   ...) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            
            if (efficient & !is.null(region)) {
              region <- Efficient_subRegion(region, x, nomapValue = nomapValue)
            }
            
            if (is.null(region)) {
              letter_freq_property <- MemorySaving_LetterFreq(
                sequence = sequence,
                gr = x,
                letters = letters,
                as.prob = as.prob,
                ...
              )
            } else if (is(region, "GRanges")) {
              letter_freq_property <- MemorySaving_LetterFreq(
                sequence = sequence,
                gr = region,
                letters = letters,
                as.prob = as.prob,
                ...
              )
            } else if (is(region, "GRangesList")) {
              region <- region[elementNROWS(region) != 0]
              names(region) <- seq_along(region)
              region_ulst <- unlist(region)
              letter_freq_property <-
                MemorySaving_LetterFreq(
                  sequence = sequence,
                  gr = region_ulst,
                  letters = letters,
                  as.prob = FALSE,
                  ...
                )
              if (!is.null(letter_freq_property)) {
                letter_freq_property <-
                  tapply(letter_freq_property, names(region_ulst), sum)
              }
              rm(region_ulst)
              
              region_widths <- sum(width(region))
              if (as.prob) {
                letter_freq_property <-
                  letter_freq_property[names(region_widths)] / region_widths
              } else{
                letter_freq_property <- letter_freq_property[names(region_widths)]
              }
              rm(region_widths)
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            
            if (is.null(region)) {
              return(as.vector(letter_freq_property))
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = as.vector(letter_freq_property),
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = nomapValue,
                ignore.strand = ignore.strand
              )
            }
            return(region_property)
          })

Efficient_subRegion <- function(region, x, nomapValue) {
  region_sub <- subsetByOverlaps(region, x, maxgap = 0L)
  if (nomapValue == "nearest") {
    x_nomap <- subsetByOverlaps(x, region, maxgap = 0L, invert = TRUE)
    if (is(region, "GRanges")) {
      region_sub <-
        c(region[na.omit(nearest(x_nomap, region))], region_sub)
    } else{
      region <- region[elementNROWS(region) != 0]
      region_ranges <- range_unique_grl(region)
      region_sub <-
        c(region[unique(na.omit(nearest(x_nomap, unlist(
          region_ranges
        ))))], region_sub)
      rm(region_ranges)
    }
    rm(x_nomap)
  }
  names(region_sub) <- NULL
  region <- unique(region_sub)
  rm(region_sub)
  region <- sort(region)
  return(region)
}

MemorySaving_LetterFreq <-
  function(sequence,
           gr,
           sequence_chunk_size = 1e9,
           letters = "GC",
           as.prob = TRUE,
           ...) {
    stopifnot(is(gr, "GRanges"))
    cumSum_indx <- cumsum(as.numeric(width(gr)))
    i = 0
    letterFreq_property <- c()
    while (!all(cumSum_indx <= sequence_chunk_size * i)) {
      i = i + 1
      gr_i <-
        gr[cumSum_indx <= sequence_chunk_size * i &
             cumSum_indx > sequence_chunk_size * (i - 1)]
      gr_seq_i <- getSeq(
        sequence,
        seqnames(gr_i),
        start = start(gr_i),
        end = end(gr_i),
        strand = strand(gr_i),
        as.character = FALSE
      )
      rm(gr_i)
      letterFreq_property <-
        c(
          letterFreq_property,
          letterFrequency(gr_seq_i, as.prob = as.prob, letters = letters, ...)
        )
      rm(gr_seq_i)
    }
    return(letterFreq_property)
  }

#' @aliases extractRegionScores
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
#gscores is a GScores object
setMethod("extractRegionScores",
          "GRanges",
          function(x,
                   gscores,
                   region = NULL,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest"),
                   missingScores = c("zero", "mean", "none"),
                   ignore.strand = FALSE,
                   efficient = TRUE,
                   ...) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            
            if (efficient & !is.null(region)) {
              region <- Efficient_subRegion(region, x, nomapValue = nomapValue)
            }
            
            missingScores <- match.arg(missingScores)
            #Define a helper function for missing scores imputation
            helper <- function(property, method) {
              if (method == "mean") {
                property[is.na(property)] <- mean(property)
              } else if (method == "zero") {
                property[is.na(property)] <- 0
              } else{
              }
              return(property)
            }
            
            if (is.null(region)) {
              gscores_property <- score(gscores, x, ...)
              gscores_property <- helper(gscores_property, missingScores)
            } else if (is(region, "GRanges")) {
              if (length(region) != 0) {
                gscores_property <- score(gscores, region, ...)
                gscores_property <- helper(gscores_property, missingScores)
              } else{
                gscores_property <-
                  rep(ifelse(missingScores == "zero", 0, NA), length(region))
              }
            } else if (is(region, "GRangesList")) {
              if (length(region) != 0) {
                region <- region[elementNROWS(region) != 0]
                names(region) <- seq_along(region)
                region_ulst <- unlist(region)
                gscores_property <- score(gscores, region_ulst, ...)
                gscores_property <- helper(gscores_property, missingScores)
                gscores_property <-
                  tapply(gscores_property, names(region_ulst), mean)
              } else{
                gscores_property <-
                  rep(ifelse(missingScores == "zero", 0, NA), length(region))
              }
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            
            if (is.null(region)) {
              return(gscores_property)
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = gscores_property,
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = nomapValue,
                ignore.strand = ignore.strand
              )
            }
            return(region_property)
          })

#' @aliases extractRegionYCount
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
#x is the GRanges to be extractated
#y is the GRanges that represent units of count
#if normalize == TRUE, then the count will be divided by the region length
setMethod("extractRegionYCount",
          "GRanges",
          function(x,
                   y,
                   region = NULL,
                   normalize = FALSE,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest"),
                   ignore.strand = FALSE,
                   efficient = TRUE) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            
            if (efficient & !is.null(region)) {
              region <- Efficient_subRegion(region, x, nomapValue = nomapValue)
            }
            
            if (is.null(region)) {
              cc_property <- countOverlaps(x,
                                           y,
                                           maxgap = maxgap,
                                           minoverlap = minoverlap,
                                           type = type)
              if (normalize)
                cc_property <- cc_property / width(x)
            } else if (is(region, "GRanges") | is(region, "GRangesList")) {
              if (is(region, "GRangesList"))
                region <- region[elementNROWS(region) != 0]
              
              cc_property <- countOverlaps(region,
                                           y,
                                           maxgap = maxgap,
                                           minoverlap = minoverlap,
                                           type = type)
              if (normalize) {
                if (is(region, "GRanges")) {
                  cc_property <- cc_property / width(region)
                } else{
                  cc_property <- cc_property / sum(width(region))
                }
              }
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            
            if (is.null(region)) {
              return(cc_property)
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = cc_property,
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = nomapValue,
                ignore.strand = ignore.strand
              )
            }
            return(region_property)
          })

#' @aliases extractRegionNearestDistToY
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
#In practice, consider exclude y within x at the beginning
setMethod("extractRegionNearestDistToY",
          "GRanges",
          function(x,
                   y,
                   region = NULL,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest"),
                   maxDist = 3e6,
                   ignore.strand = FALSE) {
            stopifnot(is(x, "GRanges"))
            stopifnot(is(y, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            if (is.null(region)) {
              dtn_property <- rep(maxDist, length(x))
              dtn_property <- distanceToNearesti(x, y, maxDist = 3e6)
            } else if (is(region, "GRanges")) {
              dtn_property <- rep(maxDist, length(region))
              dtn_property <-
                distanceToNearesti(region, y, maxDist = 3e6, ignore.strand = ignore.strand)
            } else if (is(region, "GRangesList")) {
              region <- region[elementNROWS(region) != 0]
              names(region) <- seq_along(region)
              region_ulst <- unlist(region)
              dtn_property <- rep(maxDist, length(region_ulst))
              dtn_property <-
                distanceToNearesti(region_ulst,
                                   y,
                                   maxDist = 3e6,
                                   ignore.strand = ignore.strand)
              dtn_property <- tapply(dtn_property, names(region_ulst), min)
              rm(region_ulst)
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            if (is.null(region)) {
              return(dtn_property)
            } else{
              region_property <- extractRegionProperty(
                x = x,
                property = dtn_property,
                region = region,
                ambiguityMethod = ambiguityMethod,
                maxgap = maxgap,
                minoverlap = minoverlap,
                type = type,
                nomapValue = nomapValue,
                ignore.strand = FALSE
              )
            }
            return(region_property)
          })

distanceToNearesti <- function(x,
                               y,
                               maxDist = 3e6,
                               ignore.strand = FALSE) {
  precede_indx <- precede(x, y)
  follow_indx <- follow(x, y)
  nna_prec <- !is.na(precede_indx)
  nna_follow <- !is.na(follow_indx)
  dist_precede <- rep(maxDist, length(x))
  dist_follow <- rep(maxDist, length(x))
  dist_precede[nna_prec] <-
    distance(x[nna_prec], y[precede_indx[nna_prec]], ignore.strand = ignore.strand)
  dist_follow[nna_follow] <-
    distance(x[nna_follow], y[follow_indx[nna_follow]], ignore.strand = ignore.strand)
  rm(precede_indx, follow_indx, nna_prec, nna_follow)
  logidx <- dist_follow > dist_precede
  dist_precede[logidx] <- dist_follow[logidx]
  return(dist_precede)
}

#' @aliases extractRegionRelativePosition
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
setMethod("extractRegionRelativePosition",
          "GRanges",
          function(x,
                   region = NULL,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   nomapValue = c("NA", "zero"),
                   ignore.strand = FALSE) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            nomapValue <- match.arg(nomapValue)
            if(nomapValue == "NA"){
              nomapValue <- NA
            }else{
              nomapValue <- 0
            }
            
            if (is.null(region)) {
              rrp_property <- rep(nomapValue, length(x))
            } else if (is(region, "GRanges") | is(region, "GRangesList")) {
              if (is(region, "GRanges")) {
                region_grl <- split(region, seq_along(region))
              } else{
                region <- region[elementNROWS(region) != 0]
                more_strand_region <-
                  which(elementNROWS(runValue(strand(region))) == 1)
                region <- grl_resolve_multi_strand(region)
                names(region) <- seq_along(region)
                region_grl <- region
              }
              rrp_property <- rep(nomapValue, length(x))
              map2tx <-
                mapToTranscripts(x, region_grl, ignore.strand = ignore.strand)
              relpos <-
                start(map2tx) / sum(width(region_grl))[map2tx$transcriptsHits]
              weighted_relpos <-
                tapply(relpos, map2tx$xHits, eval(parse(text = ambiguityMethod[1])))
              rm(relpos)
              rrp_property[as.numeric(names(weighted_relpos))] <-
                weighted_relpos
              rm(map2tx, weighted_relpos)
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            return(rrp_property)
          })

grl_resolve_multi_strand <- function(grl) {
  indx_multi <- elementNROWS(runValue(strand(grl))) > 1
  if (any(indx_multi)) {
    gr_ambiguous <- unlist(grl[indx_multi])
    names(gr_ambiguous) <-
      paste0(names(gr_ambiguous), strand(gr_ambiguous))
    gr_resolved <- split(gr_ambiguous, names(gr_ambiguous))
    rm(gr_ambiguous)
    return(c(grl[!indx_multi], gr_resolved))
  } else{
    return(grl)
  }
}

#' @aliases extractDistToRegion5end
#' @rdname extractRegionProperty-methods
#' @docType methods
#' @export
setMethod("extractDistToRegion5end",
          "GRanges",
          function(x,
                   region = NULL,
                   ignore.strand = FALSE,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest")) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            
            if (is.null(region)) {
              d2five_property <- width(x)
            } else if (is(region, "GRanges") | is(region, "GRangesList")) {
              if (is(region, "GRanges")) {
                region_grl <- split(region, seq_along(region))
              } else{
                region <- region[elementNROWS(region) != 0]
                region <- grl_resolve_multi_strand(region)
                names(region) <- seq_along(region)
                region_grl <- region
              }
              
              d2five_property <- rep(NA, length(x))
              map2tx <-
                mapToTranscripts(x, region_grl, ignore.strand = ignore.strand)
              d2five <- start(map2tx) - 1
              weighted_d2five <-
                tapply(d2five, map2tx$xHits, eval(parse(text = ambiguityMethod[1])))
              rm(d2five)
              d2five_property[as.numeric(names(weighted_d2five))] <-
                weighted_d2five
              rm(map2tx, weighted_d2five)
              if (nomapValue == "nearest") {
                indx_nomap <- which(is.na(d2five_property))
                d2n <- distanceToNearest(x[indx_nomap],
                                         resize(
                                           unlist(range_unique_grl(region_grl)),
                                           width = 1,
                                           fix = "start"
                                         ))
                d2five_property[indx_nomap][queryHits(d2n)] <-
                  mcols(d2n)[["distance"]]
                rm(d2n, indx_nomap)
                if (anyNA(d2five_property))
                  d2five_property[is.na(d2five_property)] <-
                  max(d2five_property, na.rm = TRUE)
              } else if (nomapValue == "zero") {
                d2five_property[is.na(d2five_property)] <- 0
              }
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            return(d2five_property)
          })

#' @rdname extractRegionProperty-methods
#' @aliases extractDistToRegion3end
#' @docType methods
#' @export
setMethod("extractDistToRegion3end",
          "GRanges",
          function(x,
                   region = NULL,
                   ignore.strand = FALSE,
                   ambiguityMethod = c("mean", "sum", "min", "max"),
                   maxgap = -1L,
                   minoverlap = 0L,
                   type = c("any", "start", "end", "within", "equal"),
                   nomapValue = c("NA", "zero", "nearest")) {
            stopifnot(is(x, "GRanges"))
            ambiguityMethod <- match.arg(ambiguityMethod)
            type <- match.arg(type)
            nomapValue <- match.arg(nomapValue)
            
            if (is.null(region)) {
              d2three_property <- rep(nomapValue, length(x))
            } else if (is(region, "GRanges") | is(region, "GRangesList")) {
              if (is(region, "GRanges")) {
                region_grl <- split(region, seq_along(region))
              } else{
                region <- region[elementNROWS(region) != 0]
                region <- grl_resolve_multi_strand(region)
                names(region) <- seq_along(region)
                region_grl <- region
              }
              d2three_property <- rep(NA, length(x))
              map2tx <-
                mapToTranscripts(x, region_grl, ignore.strand = ignore.strand)
              d2three <-
                sum(width(region_grl))[map2tx$transcriptsHits] - start(map2tx)
              weighted_d2three <-
                tapply(d2three, map2tx$xHits, eval(parse(text = ambiguityMethod[1])))
              rm(d2three)
              d2three_property[as.numeric(names(weighted_d2three))] <-
                weighted_d2three
              rm(map2tx, weighted_d2three)
              if (nomapValue == "nearest") {
                indx_nomap <- which(is.na(d2three_property))
                d2n <- distanceToNearest(x[indx_nomap],
                                         resize(
                                           unlist(range_unique_grl(region_grl)),
                                           width = 1,
                                           fix = "end"
                                         ))
                d2three_property[indx_nomap][queryHits(d2n)] <-
                  mcols(d2n)[["distance"]]
                rm(d2n, indx_nomap)
                if (anyNA(d2three_property))
                  d2three_property[is.na(d2three_property)] <-
                  max(d2three_property, na.rm = TRUE)
              } else if (nomapValue == "zero") {
                d2three_property[is.na(d2three_property)] <- 0
              }
            } else{
              stop("`region` should be either `GRanges` or `GRangesList`")
            }
            return(d2three_property)
          })
