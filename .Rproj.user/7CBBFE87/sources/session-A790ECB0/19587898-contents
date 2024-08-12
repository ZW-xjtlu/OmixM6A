#' @title Extract genome-derived predictive features from interval-based data
#'
#' @description A function to extract genome-derived features from the input \code{\link{GRanges}} object and the annotations specified by TxDb-like objects.
#'
#' @param x A \link{GRanges} object for the genomic ranges to be extracted.
#' @param transcriptdb A \link{TxDb} or \link{EnsDb} object for the transcript annotation.
#' @param sequence A \link{BSgenome} or \link{XStringSet} object for the genome sequence.
#' See the \link{available.genomes} function for how to install a genome.
#' @param gscores A \link{GScores} object specifying the genomic scores to be extracted, or a \code{list} of \code{GScores} for creating properties using more than one type of genomic scores.
#' See the vignette for more details: \code{vignette("GenomicScores", package = "GenomicScores")}.
#' @param clusteringY A \link{GRanges} or \link{GRangesList} specifying the object for clustering annotation, or a \code{list} of \code{GRanges}/\code{GRangesList} for creating clustering properties on more than one annotations.
#' @param extraRegions A \link{GRanges} or \link{GRangesList} providing the additional region to extract properties, can be a \code{list} of \code{GRanges} or \code{GRangesList} to supplement more than one regions.
#' @param flankSizes An \code{integer} vector indimessageing the number of bases flanked by \code{x} when calculating the \code{x} centered properties; default\code{=(25*2^(0:6))}.
#' @param ambiguityMethod By default, except the overlapping features, if \code{x[i]} overlaps multiple regions, the returned property will be the average of the properties over all overlapped regions by \code{x[i]}.
#' For the overlapping features, as long as x overlaps any range of the region, the returned value is 1 .
#'
#' If \code{ambiguityMethod} is \code{"mean"}, \code{"sum"}, \code{"min"}, or \code{"max"}, then the mean, sum, minimum, and maximum values of the >1 mapping will be returned.
#' @param nomapValue When \code{nomapValue} is \code{"NA"} or \code{"zero"}, the \code{x} that do not match the region will return \code{NA} and \code{0} respectively.
#'
#' If \code{nomapValue} is \code{"nearest"}, the not matched \code{x} will be set to be the properties on its nearest region.
#'
#' @param annotSeqnames \code{TRUE} or \code{FALSE}. If \code{TRUE}, the seqnames in \code{x} are output as features, with each level represented by a dummy variable column.
#' @param annotBiotype \code{TRUE} or \code{FALSE}. If \code{TRUE}, the transcript biotypes defined in the \code{EnsDb} are output as features, with each biotype represented by a dummy variable column.
#‘
#' @return A \code{data.frame} object whose number of rows is the length of \code{x}.
#' The column types in the \code{data.frame} are all numeric.
#'
#' @details
#' The function first extract multiple genome-derived properties on x and 13 default region types defined by the TxDb-like object,
#' and it will then assign the region's properties to \code{x} if \code{x} overlapped with the region.
#' The specific behaviors of the assignments can be modified through the parameters of \code{ambiguityMethod} and \code{nomapValue}.
#'
#' Particularly, the 13 default genome region types include \emph{exons}, \emph{introns}, \emph{5'UTR} (exons only), \emph{full 5'UTR}, \emph{CDS} (exons only), \emph{full CDS},
#' \emph{3'UTR} (exons only), \emph{full 3'UTR}, \emph{mature transcripts}, \emph{full transcripts}, \emph{genes} (exons only), \emph{full genes}, and \emph{promoters}.
#'
#' For each region type, the function can generate 5 types of basic properties accordingly:
#' \emph{overlap with region},
#' \emph{region length},
#' \emph{relative position of \code{x} on the region},
#' \emph{distance to region 5'end}, and
#' \emph{distance to region 3'end.}
#'
#' In addition, when the corresponding annotation objects are provided, the function can add more types of properties.
#' Specifically, GC contents of regions are extracted if \code{sequence} is provided;
#' Genomic Scores of the regions are extracted if \code{gscores} is provided; and
#' Clustering effects of annotations on regions, including Count, Densities, and the Distance to the nearest annotation
#' are calculated if \code{clusteringY} is provided.
#'
#' Genomic regions other the 13 basic types can be defined by supplying the \code{list} of \code{GRanges} objects at \code{extraRegions}.
#' Similarly, the input for \code{gscores}, and \code{clusteringY} can be a list so more than one properties can be added for these genomic metrics.
#' Please note that the \code{names} of the \code{list} elements will be used to generate the feature names in the output.
#'
#' Finally, the function will generate 3 extra features to describe the properties uniquely defined at the gene's level:
#' \emph{gene's exon number},
#' \emph{gene's transcript isoform number} and
#' \emph{meta-tx topology}.
#'
#' Please see the package vignettes for the detailed information of those features.
#'
#' @examples
#' ## ---------------------------------------------------------------------
#' ## SIMPLE EXAMPLE
#' ## ---------------------------------------------------------------------
#'
#' ## Load the hg19 TxDb object for human transcript annotation (UCSC hg19):
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' ## Define the Granges to be extracted:
#' set.seed(01)
#'
#' X <- GRanges(rep(c("chr1", "chr2"), c(15, 15)),
#'              IRanges(c(sample(11874:12127, 15), sample(38814:41527,15)), width=5),
#'              strand=Rle(c("+", "-"), c(15, 15)))
#'
#' ## Extract the basic set of properties using the genomic regions defined in TxDb:
#' gfeatures <- genomeDerivedFeatures(X,
#'                                    transcriptdb=TxDb.Hsapiens.UCSC.hg19.knownGene)
#' str(gfeatures)
#'
#' ## ---------------------------------------------------------------------
#' ## ADD MORE PROPERTIES AND REGIONS
#' ## ---------------------------------------------------------------------
#' \donttest{
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(phastCons100way.UCSC.hg19)
#'
#' ## Extract more genomic properties derived from genome sequence and phastCons score:
#' gfeatures_expanded <- genomeDerivedFeatures(X,
#'                                             transcriptdb=TxDb.Hsapiens.UCSC.hg19.knownGene,
#'                                             sequence=BSgenome.Hsapiens.UCSC.hg19,
#'                                             gscores=phastCons100way.UCSC.hg19,
#'                                             clusteringY=X) #quantify clustering properties of x itself
#' str(gfeatures_expanded)
#'
#' ## Extract region properties features from Ensemble transcript annotation (EnsDb):
#' library(EnsDb.Hsapiens.v75)
#'
#' gfeatures_ensdb <- genomeDerivedFeatures(X,
#'                                          transcriptdb=EnsDb.Hsapiens.v75,
#'                                          sequence=BSgenome.Hsapiens.UCSC.hg19,
#'                                          gscores=phastCons100way.UCSC.hg19,
#'                                          clusteringY=X,
#'                                          annotBiotype=TRUE) #extract transcript biotypes in EnsDb
#' str(gfeatures_ensdb)
#' }
#'
#' @seealso
#'
#' \itemize{
#' \item{}{The \link{extractRegionProperty} for crafting of various properties on genomic regions.}
#' \item{}{The \link{sequenceDerivedFeatures} for extraction of sequence-derived features under different encoding schema.}
#' \item{}{The \link{topologyOnTranscripts} for calculation of the meta-tx topology on transcripts of genes.}
#' }
#'
#' @author Zhen Wei
#'
#' @importFrom GenomicFeatures exonsBy intronsByTranscript transcripts threeUTRsByTranscript fiveUTRsByTranscript cdsBy exons genes promoters transcriptsBy
#' @importFrom IRanges subsetByOverlaps %over%
#' @importFrom ensembldb exonsBy intronsByTranscript transcripts threeUTRsByTranscript fiveUTRsByTranscript cdsBy exons genes promoters
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<- genome<- keepStandardChromosomes seqlengths<- seqnames isCircular isCircular<- seqlevels seqlevels<-
#' @importFrom S4Vectors mcols<- mcols
#' @importFrom gtools mixedsort
#' @importFrom utils packageVersion
#' @export
genomeDerivedFeatures <- function(x,
                                  transcriptdb,
                                  sequence = NULL,
                                  gscores = NULL,
                                  clusteringY = NULL,
                                  extraRegions = NULL,
                                  flankSizes = (25 * 2 ^ (0:6)),
                                  ambiguityMethod = c("auto", " mean", "sum", "min", "max"),
                                  nomapValue = c("zero", "NA", "nearest"),
                                  annotSeqnames = FALSE,
                                  annotBiotype = FALSE) {
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  genome(x) <- NA
  isCircular(x) <- NA
  
  if (!is.null(clusteringY)) {
    if (is.list(clusteringY)) {
      if (is.null(names(clusteringY))) {
        names(clusteringY) <- paste0("y_", seq_along(clusteringY))
      }
      for (i in seq_along(clusteringY)) {
        stopifnot(is(clusteringY[[i]], "GRanges"))
        genome(clusteringY[[i]]) <- NA
        isCircular(clusteringY[[i]]) <- NA
      }
    } else{
      stopifnot(is(clusteringY, "GRanges"))
      genome(clusteringY) <- NA
      isCircular(clusteringY) <- NA
    }
  }
  
  if (!is.null(extraRegions)) {
    if (is.list(extraRegions)) {
      if (is.null(names(extraRegions))) {
        names(extraRegions) <-
          paste0("extraRegion_", seq_along(extraRegions))
      }
      for (i in seq_along(extraRegions)) {
        stopifnot(is(extraRegions[[i]], "GRanges") |
                    is(extraRegions[[i]], "GRangesList"))
        genome(extraRegions[[i]]) <- NA
        isCircular(extraRegions[[i]]) <- NA
      }
    } else{
      stopifnot(is(extraRegions, "GRanges") |
                  is(extraRegions, "GRangesList"))
      genome(extraRegions) <- NA
      isCircular(extraRegions) <- NA
      extraRegions <- list(extraRegion = extraRegions)
    }
  }
  
  if (!is.null(gscores)) {
    if (is.list(gscores)) {
      if (is.null(names(gscores))) {
        names(gscores) <- paste0("gscores_", seq_along(gscores))
      }
      for (i in seq_along(gscores)) {
        stopifnot(is(gscores[[i]], "GScores"))
      }
    } else{
      stopifnot(is(gscores, "GScores"))
    }
  }
  
  if (annotSeqnames) {
    xSeqlevels <- as.character(unique(seqnames(x)))
    xSeqlevels <- mixedsort(xSeqlevels)
  }
  
  ##Fill unkown attributes of x with the attributes of BSgenome
  if (anyNA(isCircular(x)) & !is.null(sequence)) {
    isCircularGenome <- isCircular(sequence)
    isCircular(x) <- isCircularGenome[names(isCircular(x))]
    rm(isCircularGenome)
  }
  if (anyNA(seqlengths(x)) & !is.null(sequence)) {
    seqLenGenome <- seqlengths(sequence)
    seqlengths(x) <- seqLenGenome[names(seqlengths(x))]
    rm(seqLenGenome)
  }
  
  message_env <- new.env()
  message_env$count <- 1
  message(
    "#############################################################################################################################################\n"
  )
  message(
    "##                                                     predictiveFeatures Package Version: ",
    paste0(
      packageVersion("predictiveFeatures"),
      paste(rep(" ", 53 - nchar(
        packageVersion("predictiveFeatures")
      )), collapse = ""),
      "##\n"
    )
  )
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message_env$total <- calculate_total(
    transcriptdb = transcriptdb,
    bsgenome = sequence,
    gscores = gscores,
    clusteringY = clusteringY,
    flankSizes = flankSizes,
    annotBiotype = annotBiotype,
    extraRegions = extraRegions
  ) + ifelse(annotSeqnames, length(xSeqlevels), 0)
  message(
    "##                                   The total number of genome-derived features to be extracted: ",
    paste0(message_env$total, paste(rep(
      " ", 41 - nchar(as.character(message_env$total))
    ), collapse = ""), "##\n")
  )
  
  names(x) <- seq_along(x)
  x <- sort(x)
  
  #Features for x
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties over flanking regions of x                                             ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <- EnumerateRegionFeatures(
    x,
    bsgenome = sequence,
    gscores = gscores,
    clusteringY = clusteringY,
    flankSizes = flankSizes,
    ambiguityMethod = ambiguityMethod,
    nomapValue = nomapValue,
    message_env = message_env
  )
  
  if (annotSeqnames) {
    for (i in xSeqlevels) {
      Message_i(paste0("seqname_", i),
                paste0("seqname of x is ", i),
                message_env)
      
      X[[paste0("seqname_x_", i)]] <- as.numeric(seqnames(x) == i)
      
      message("Done")
    }
    
    rm(xSeqlevels)
    
  }
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exons                                                               ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(exons(transcriptdb), x),
        region_name = "exons",
        region_info = "exons",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of introns                                                             ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(unlist(intronsByTranscript(transcriptdb)), x),
        region_name = "introns",
        region_info = "introns",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exonic 5'UTR                                                        ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(fiveUTRsByTranscript(transcriptdb), x),
        region_name = "exonicFivePrimeUTR",
        region_info = "exonic 5'UTR",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of full 5'UTR                                                          ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(unlist(range(
          fiveUTRsByTranscript(transcriptdb)
        )), x),
        region_name = "fullFivePrimeUTR",
        region_info = "full 5'UTR",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exonic CDS                                                          ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(cdsBy(transcriptdb, by = "tx"), x),
        region_name = "exonicCDS",
        region_info = "exonic CDS",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of full CDS                                                            ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(unlist(range(
          cdsBy(transcriptdb, by = "tx")
        )), x),
        region_name = "fullCDS",
        region_info = "full CDS",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exonic 3'UTR                                                        ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(threeUTRsByTranscript(transcriptdb), x),
        region_name = "exonicThreePrimeUTR",
        region_info = "exonic 3'UTR",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of full 3'UTR                                                          ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(unlist(range(
          threeUTRsByTranscript(transcriptdb)
        )), x),
        region_name = "fullThreePrimeUTR",
        region_info = "full 3'UTR",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exonic transcripts                                                  ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(exonsBy(transcriptdb, "tx"), x),
        region_name = "exonicTranscripts",
        region_info = "exonic tx",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of full transcripts                                                    ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(transcripts(transcriptdb), x),
        region_name = "fullTranscripts",
        region_info = "full tx",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of intronic transcripts                                                ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(intronsByTranscript(transcriptdb), x),
        region_name = "intronicTranscripts",
        region_info = "intronic tx",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of exonic genes                                                        ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(exonsBy(transcriptdb, "gene"), x),
        region_name = "exonicGenes",
        region_info = "exonic genes",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  
  if (is(transcriptdb, "TxDb")) {
    fullGenes <-
      unlist(genes(transcriptdb, single.strand.genes.only = FALSE))
  } else{
    fullGenes <- genes(transcriptdb)
  }
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties on full genes                                                       ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <-
    cbind(
      X,
      EnumerateRegionFeatures(
        x,
        region = EnsureUCSC(fullGenes, x),
        region_name = "fullGenes",
        region_info = "full genes",
        bsgenome = sequence,
        gscores = gscores,
        clusteringY = clusteringY,
        flankSizes = flankSizes,
        ambiguityMethod = ambiguityMethod,
        nomapValue = nomapValue,
        message_env = message_env
      )
    )
  rm(fullGenes)
  
  promo <- promoters(transcriptdb)
  
  if (any(isCircular(promo), na.rm = TRUE)) {
    promo <- promo[seqnames(promo) != names(which(isCircular(promo)))]
  }
  
  
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  message("##                                               Extract properties of promoters                                                           ##\n")                                              ##\n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
  )
  X <- cbind(
    X,
    EnumerateRegionFeatures(
      x,
      region = EnsureUCSC(promo, x),
      region_name = "promoters",
      region_info = "promoters",
      bsgenome = sequence,
      gscores = gscores,
      clusteringY = clusteringY,
      flankSizes = flankSizes,
      ambiguityMethod = ambiguityMethod,
      nomapValue = nomapValue,
      message_env = message_env
    )
  )
  rm(promo)
  
  if (!is.null(extraRegions)) {
    for (i in names(extraRegions)) {
      message(
        "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
      )
      message("##                                               Extract properties of ", i , paste0(paste0(rep(
        " ", 68 - nchar(i)
      ), collapse = "")), "##\n")
      message(
        "## --------------------------------------------------------------------------------------------------------------------------------------- ##\n"
      )
      
    }
    
    
    X <-
      cbind(
        X,
        EnumerateRegionFeatures(
          x,
          region = EnsureUCSC(extraRegions[[i]], x),
          region_name = i,
          region_info = i,
          bsgenome = sequence,
          gscores = gscores,
          clusteringY = clusteringY,
          flankSizes = flankSizes,
          ambiguityMethod = ambiguityMethod,
          nomapValue = nomapValue,
          message_env = message_env
        )
      )
    
  }
  
  
  #EnsDb Biotype
  if (annotBiotype & (!is(transcriptdb, "EnsDb"))) {
    warning(
      "The transcriptdb is not an EnsDb object, hence the argument setting `annotBiotype=TRUE` is ignored."
    )
  }
  
  if (is(transcriptdb, "EnsDb") & annotBiotype) {
    message(
      "## --------------------------------------------------------------------------------------------------------------------------------------- ## \n"
    )
    message("##                                              Extracting transcript biotype                                                              ## \n")                                              ##\n")
    message(
      "## --------------------------------------------------------------------------------------------------------------------------------------- ## \n"
    )
    
    message("Extract features of transcript biotype ...\n")
    TxEnsDb <-
      EnsureUCSC(transcripts(transcriptdb), x, clean_columns = FALSE)
    for (i in unique(TxEnsDb$tx_biotype)) {
      Message_i(i, paste0("EnsDb biotype ", i), message_env)
      X[[i]] <-
        suppressWarnings(x %over% TxEnsDb[TxEnsDb$tx_biotype == i])
      message("Done")
    }
  }
  
  #Some additional calculations
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ## \n"
  )
  message("##                                                Extracting other genomic metrics                                                         ## \n")                                              ##\n")
  message(
    "## --------------------------------------------------------------------------------------------------------------------------------------- ## \n"
  )
  
  introns <- EnsureUCSC(unlist(intronsByTranscript(transcriptdb)), x)
  Message_i("log2_NearestDistToJunction",
            "log2(nearest distance to splicing junctions + 1)",
            message_env)
  NearestDistToJunction <- rep(NA, length(x))
  introns <- unlist(intronsByTranscript(transcriptdb))
  dnn <- distanceToNearest(x, c(resize(introns, 1, "start"), resize(introns, 1, "end")))
  NearestDistToJunction[queryHits(dnn)] <- mcols(dnn)$distance
  if(nomapValue == "zero") { 
    NearestDistToJunction[is.na(NearestDistToJunction)] <- 0 
    }
  X[["log2_NearestDistToJunction"]] <- log2(
  NearestDistToJunction + 1
  )
  message("Done")
  rm(introns, dnn, NearestDistToJunction)
  
  exbg <- EnsureUCSC(exonsBy(transcriptdb, "gene"), x)
  Message_i("log2_GeneExonNumber",
            "log2(exon number of genes + 1)",
            message_env)
  X[["log2_GeneExonNumber"]] <- log2(
    extractRegionProperty(
      x,
      property = elementNROWS(exbg),
      region = range(exbg),
      nomapValue = nomapValue
    ) + 1
  )
  message("Done")
  rm(exbg)
  
  Message_i("log2_TxIsoformNumber",
            "log2(transcript isoform number of genes + 1)",
            message_env)
  txbg <- EnsureUCSC(transcriptsBy(transcriptdb, "gene"), x)
  X[["log2_TxIsoformNumber"]] <- log2(
    extractRegionProperty(
      x,
      property = elementNROWS(txbg),
      region = range(txbg),
      nomapValue = nomapValue
    ) + 1
  )
  message("Done")
  rm(txbg)
  
  Message_i("MetaTxTopology", "meta-transcript topology", message_env)
  
  X[["MetaTxTopology"]] <-
    suppressWarnings(topologyOnTranscripts(x, transcriptdb))
  
  X[["MetaTxTopology"]][is.na(X[["MetaTxTopology"]])] <- 0
  
  message("Done")
  
  message(
    "#############################################################################################################################################\n"
  )
  message("##                                      The generation of genomic predictive features has been completed.                                 ##\n")
  message(
    "#############################################################################################################################################\n"
  )
  
  indx_features <- match(as.character(seq_along(x)), names(x))
  
  X <- X[indx_features,]
  
  rm(indx_features)
  
  rownames(X) <- NULL
  
  return(X)
}

EnumerateRegionFeatures <- function(x,
                                    region = NULL,
                                    region_name = NULL,
                                    region_info = NULL,
                                    bsgenome = NULL,
                                    gscores = NULL,
                                    clusteringY = NULL,
                                    flankSizes = (25 * 2 ^ (0:6)),
                                    ambiguityMethod = c("auto", "mean", "sum", "min", "max"),
                                    nomapValue = c("zero", "NA", "nearest"),
                                    message_env) {
  ambiguityMethod <- match.arg(ambiguityMethod)
  nomapValue <- match.arg(nomapValue)
  #Length of x
  
  if (is.null(region)) {
    Message_i("log2_length_x", "log2(length of x + 1)", message_env)
    X <-
      data.frame(log2_length_x = log2(extractRegionLength(x) + 1))
    message("Done")
    
    #GC content  of x and its flanks
    isCircular(x)[isCircular(x)] <- FALSE
    if (!is.null(bsgenome)) {
      Message_i("GC_x", "GC content of x", message_env)
      X[["GC_x"]] <-
        quiet(extractRegionLetterFrequency(trim(x), bsgenome))
      message("Done")
      for (i in flankSizes) {
        Message_i(
          paste0("GC_flank_", i, "_x"),
          paste0("GC content of x flanked by ", i, ""),
          message_env
        )
        X[[paste0("GC_flank_", i, "_x")]] <-
          quiet(extractRegionLetterFrequency(trim(x + i), bsgenome))
        message("Done")
      }
    }
    
    #Gscores of x and its flanks
    if (!is.null(gscores)) {
      if (is.list(gscores)) {
        for (i in seq_along(gscores)) {
          Message_i(paste0(names(gscores)[i], "_x"),
                    paste0(names(gscores)[i], " scores of x"),
                    message_env)
          X[[paste0(names(gscores)[i], "_x")]] <-
            quiet(extractRegionScores(trim(x), gscores[[i]], missingScores = "zero"))
          message("Done")
          for (j in flankSizes) {
            Message_i(
              paste0(names(gscores)[i], "_flank_", j, "_x"),
              paste0(names(gscores)[i], " scores of x flanked by ", j),
              message_env
            )
            X[[paste0(names(gscores)[i], "_flank_", j, "_x")]] <-
              quiet(extractRegionScores(trim(x + j), gscores[[i]], missingScores = "zero"))
            message("Done")
          }
        }
      } else{
        Message_i("gscores_x", "gscores of x", message_env)
        X[["gscores_x"]] <-
          quiet(extractRegionScores(trim(x), gscores, missingScores = "zero"))
        message("Done")
        for (j in flankSizes) {
          Message_i(
            paste0("gscores_flank_", j, "_x"),
            paste0("gscores of x flanked by ", j),
            message_env
          )
          X[[paste0("gscores_flank_", j, "_x")]] <-
            quiet(extractRegionScores(trim(x + j), gscores, missingScores = "zero"))
          message("Done")
        }
      }
    }
    
    #Count y on x and its flanks
    if (!is.null(clusteringY)) {
      if (is.list(clusteringY)) {
        for (i in seq_along(clusteringY)) {
          Message_i(
            paste0("log2_", names(clusteringY)[i], "Count_x"),
            paste0("log2(count of ", names(clusteringY)[i] , " on x + 1)"),
            message_env
          )
          X[[paste0("log2_", names(clusteringY)[i], "Count_x")]] <-
            quiet(log2(
              extractRegionYCount(x, clusteringY[[i]], normalize = FALSE) + 1
            ))
          message("Done")
          for (j in flankSizes) {
            Message_i(
              paste0(
                "log2_",
                names(clusteringY)[i],
                "Count_flank_",
                j,
                "_x"
              ),
              paste0(
                "log2(count of ",
                names(clusteringY)[i] ,
                " on x flanked by ",
                j,
                " + 1)"
              ),
              message_env
            )
            X[[paste0("log2_",
                      names(clusteringY)[i],
                      "Count_flank_",
                      j,
                      "_x")]] <-
              quiet(log2(
                extractRegionYCount(trim(x + j), clusteringY[[i]], normalize = FALSE) + 1
              ))
            message("Done")
          }
          #Quantify nearest distance of y to x
          Message_i(
            paste0("log2_nearestDist2", names(clusteringY)[i], "_x"),
            paste0(
              "log2(x nearest distance to ",
              names(clusteringY)[i],
              " + 1)"
            ),
            message_env
          )
          X[[paste0("log2_nearestDist2", names(clusteringY)[i], "_x")]] <-
            quiet(log2(extractRegionNearestDistToY(x, y = clusteringY[[i]]) + 1))
          message("Done")
        }
      } else{
        Message_i("log2_yCount_x", "log2(count of y on x + 1)", message_env)
        X[["log2_yCount_x"]] <-
          quiet(log2(
            extractRegionYCount(x, clusteringY, normalize = FALSE) + 1
          ))
        message("Done")
        for (i in flankSizes) {
          Message_i(
            paste0("log2_yCount_flank_", i, "_x"),
            paste0("log2(Count of y on x flanked by ", i, " + 1)"),
            message_env
          )
          X[[paste0("log2_yCount_flank_", i, "_x")]] <-
            quiet(log2(
              extractRegionYCount(trim(x + i), clusteringY, normalize = FALSE) + 1
            ))
          message("Done")
        }
        #Quantify nearest distance of y to x
        Message_i("log2_nearestDist2y_x",
                  "log2(nearest distance of y to x + 1)",
                  message_env)
        X[["log2_nearestDist2y_x"]] <-
          quiet(log2(extractRegionNearestDistToY(x, y = clusteringY) + 1))
        message("Done")
      }
    }
    return(X)
  } else {
    #Overlap region
    Message_i(
      paste0("overlap_", region_name),
      paste0("Overlap with ", region_info),
      message_env
    )
    X <- data.frame(
      overlap_region = extractRegionOverlap(
        x,
        region,
        output.logical = FALSE,
        ambiguityMethod = ambiguityMethod
      )
    )
    names(X) = paste0("overlap_", region_name)
    message("Done")
    
    #Reset argument
    if (ambiguityMethod == "auto")
      ambiguityMethod <- "mean"
    
    #Length of region
    Message_i(
      paste0("log2_length_", region_name),
      paste0("log2(length of ", region_info , " + 1)"),
      message_env
    )
    
    if (length(X != 0)) {
      X[[paste0("log2_length_", region_name)]] <-
        log2(
          extractRegionLength(
            x,
            region,
            nomapValue = nomapValue,
            ambiguityMethod = ambiguityMethod
          ) + 1
        )
    }
    message("Done")
    
    #GC content of region
    if (!is.null(bsgenome)) {
      Message_i(paste0("GC_", region_name),
                paste0("GC content of ", region_info),
                message_env)
      X[[paste0("GC_", region_name)]] <-
        quiet(
          extractRegionLetterFrequency(
            x,
            bsgenome,
            region,
            nomapValue = nomapValue,
            ambiguityMethod = ambiguityMethod
          )
        )
      message("Done")
    }
    
    #GScores on region
    if (!is.null(gscores)) {
      if (is.list(gscores)) {
        for (i in seq_along(gscores)) {
          Message_i(
            paste0(names(gscores)[i], "_", region_name),
            paste0(names(gscores)[i], " scores of ", region_info),
            message_env
          )
          X[[paste0(names(gscores)[i], "_", region_name)]] <-
            quiet(
              extractRegionScores(
                trim(x),
                gscores[[i]],
                region,
                nomapValue = nomapValue,
                ambiguityMethod = ambiguityMethod
              )
            )
          message("Done")
        }
      } else{
        Message_i(
          paste0("gscores_", region_name),
          paste0("gscores of ", region_info),
          message_env
        )
        X[[paste0("gscores_", region_name)]] <-
          quiet(
            extractRegionScores(
              trim(x),
              gscores,
              region,
              nomapValue = nomapValue,
              ambiguityMethod = ambiguityMethod
            )
          )
        message("Done")
      }
    }
    
    #Count y on region
    if (!is.null(clusteringY)) {
      if (is.list(clusteringY)) {
        for (i in seq_along(clusteringY)) {
          Message_i(
            paste0(names(clusteringY)[i], "Density_", region_name),
            paste0("Density of ", names(clusteringY)[i], " on ", region_info),
            message_env
          )
          X[[paste0(names(clusteringY)[i], "Density_", region_name)]] <-
            quiet(
              extractRegionYCount(
                x,
                clusteringY[[i]],
                region,
                normalize = TRUE,
                nomapValue = nomapValue,
                ambiguityMethod = ambiguityMethod
              )
            )
          message("Done")
          
          Message_i(
            paste0("log2_", names(clusteringY)[i], "Count_", region_name),
            paste0(
              "log2(count of ",
              names(clusteringY)[i],
              " on ",
              region_info,
              " + 1)"
            ),
            message_env
          )
          X[[paste0("log2_", names(clusteringY)[i], "Count_", region_name)]] <-
            quiet(log2(
              extractRegionYCount(
                x,
                clusteringY[[i]],
                region,
                normalize = FALSE,
                nomapValue = nomapValue,
                ambiguityMethod = ambiguityMethod
              ) + 1
            ))
          message("Done")
          
          Message_i(
            paste0(
              "log2_nearestDist2",
              names(clusteringY)[i],
              "_",
              region_name
            ),
            paste0(
              "log2(",
              region_info,
              " nearest distance to ",
              names(clusteringY)[i],
              " + 1)"
            ),
            message_env
          )
          X[[paste0("log2_nearestDist2",
                    names(clusteringY)[i],
                    "_",
                    region_name)]] <- quiet(log2(
                      extractRegionNearestDistToY(
                        x,
                        clusteringY[[i]],
                        nomapValue = nomapValue,
                        region,
                        ambiguityMethod = ambiguityMethod
                      ) + 1
                    ))
          message("Done")
        }
      } else{
        Message_i(
          paste0("yDensity_", region_name),
          paste0("density of y on ", region_info),
          message_env
        )
        X[[paste0("yDensity_", region_name)]] <-
          quiet(
            extractRegionYCount(
              x,
              clusteringY,
              region,
              normalize = TRUE,
              nomapValue = nomapValue,
              ambiguityMethod = ambiguityMethod
            )
          )
        message("Done")
        
        Message_i(
          paste0("log2_yCount_", region_name),
          paste0("log2(count of y on ", region_info , " + 1)"),
          message_env
        )
        X[[paste0("log2_yCount_", region_name)]] <-
          quiet(log2(
            extractRegionYCount(
              x,
              clusteringY,
              region,
              normalize = FALSE,
              nomapValue = nomapValue,
              ambiguityMethod = ambiguityMethod
            ) + 1
          ))
        message("Done")
        
        Message_i(
          paste0("log2_", region_name, "_nearestDist2y"),
          paste0("log2(", region_info, " nearest distance to y + 1)"),
          message_env
        )
        X[[paste0("log2_", region_name, "_nearestDist2y")]] <-
          quiet(log2(
            extractRegionNearestDistToY(
              x,
              clusteringY,
              nomapValue = nomapValue,
              region,
              ambiguityMethod = ambiguityMethod
            ) + 1
          ))
        message("Done")
      }
    }
    #Relative Position
    Message_i(
      paste0("relativePOS_", region_name),
      paste0("relative position of x on ", region_info),
      message_env
    )
    X[[paste0("relativePOS_", region_name)]] <-
      quiet(extractRegionRelativePosition(x, region,
                                          nomapValue = nomapValue))
    message("Done")
    
    #Distance to region five prime end
    Message_i(
      paste0("log2_dist5prime_", region_name),
      paste0("log2(distance of x to 5' end of ", region_info, " + 1)"),
      message_env
    )
    X[[paste0("log2_dist5prime_", region_name)]] <-
      quiet(log2(
        extractDistToRegion5end(
          x,
          region,
          nomapValue = nomapValue,
          ambiguityMethod = ambiguityMethod
        ) + 1
      ))
    message("Done")
    
    #Distance to region three prime end
    Message_i(
      paste0("log2_dist3prime_", region_name),
      paste0("log2(distance of x to 3' end of ", region_info, " + 1)"),
      message_env
    )
    X[[paste0("log2_dist3prime_", region_name)]] <-
      quiet(log2(
        extractDistToRegion3end(
          x,
          region,
          nomapValue = nomapValue,
          ambiguityMethod = ambiguityMethod
        ) + 1
      ))
    message("Done")
    
    return(X)
  }
}

EnsureUCSC <- function(range_object, x, clean_columns = TRUE) {
  if (all(seqlevelsStyle(range_object) %in% "UCSC")) {
    
  } else{
    seqlevelsStyle(range_object) <- "UCSC"
    range_object <-
      keepStandardChromosomes(range_object, pruning.mode = c("coarse"))
    seqlevels(range_object) <-
      c(seqlevels(range_object), setdiff(seqlevels(x), seqlevels(range_object)))
    new_seqlengths <- seqlengths(x)[names(seqlengths(range_object))]
    names(new_seqlengths) <- names(seqlengths(range_object))
    seqlengths(range_object) <- new_seqlengths
    new_IsCircular <- isCircular(x)[names(isCircular(range_object))]
    names(new_IsCircular) <- names(isCircular(range_object))
    isCircular(range_object) <- new_IsCircular
  }
  if (clean_columns) {
    mcols(range_object) <- NULL
    names(range_object) <- NULL
  }
  return(range_object)
}


calculate_total <- function(transcriptdb,
                            bsgenome,
                            gscores,
                            clusteringY,
                            flankSizes,
                            annotBiotype,
                            extraRegions) {
  f_num <- length(flankSizes)
  if (is.null(gscores)) {
    s_num <- 0
  } else if (is.list(gscores)) {
    s_num <- length(gscores)
  } else{
    s_num <- 1
  }
  
  if (is.null(bsgenome)) {
    I_genome <- 0
  } else{
    I_genome <- 1
  }
  
  if (is.null(clusteringY)) {
    c_num <- 0
  } else if (is.list(clusteringY)) {
    c_num <- length(clusteringY)
  } else{
    c_num <- 1
  }
  self_feature_num <-
    1 + I_genome * (1 + f_num) + s_num * (1 + f_num) + c_num * (1 + f_num +
                                                                  1)
  region_feature_num <-
    (14 + length(extraRegions)) * (2 + I_genome + s_num + 3 * c_num + 3)
  if (annotBiotype & is(transcriptdb, "EnsDb")) {
    additional_feature_num <-
      4 + length(unique(transcripts(transcriptdb)$tx_biotype))
  } else{
    additional_feature_num <- 4
  }
  return(self_feature_num + region_feature_num + additional_feature_num)
}

Message_i <- function(Feature_name,
                      Feature_info,
                      message_env) {
  fcount <- get("count", envir = message_env)
  message(
    "Feature ",
    fcount,
    "/",
    get("total", envir = message_env),
    "; Name: '",
    Feature_name,
    "'; Info: '",
    Feature_info,
    "'; Extracting ... ",
    appendLF = FALSE
  )
  assign("count", fcount + 1, envir = message_env)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(suppressMessages(suppressWarnings(x))))
}

##To Do: support the selection of properties
