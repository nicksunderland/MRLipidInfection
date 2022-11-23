#' GenomeRegion class
#'
#' This class define an area of the genome and handles various parameters passed to the LD clumping
#' algorithms. This is useful when working with SNPs surrounding multiple gene targets of interest.
#' To perform normal genome wide SNP clumping simple create a GenomeRegion object with
#' 'genome_wide=TRUE'.
#'
#' @slot name character. The name of the gene or area of the genome referred to.
#' @slot assembly character. The assembly information, e.g. GRCh37
#' @slot chromosome numeric. The chromosome the gene is on.
#' @slot start numeric. The start (number of bases).
#' @slot end numeric. The end (number of bases).
#' @slot cis_tol when defining SNPs on a gene, set a tolerance of how many bases either side of the
#'   gene
#' @slot genome_wide this is usually FALSE, if TRUE this is treated as a object placeholder and will
#'   refer to the 'whole genome'
#' @slot clump_kb Clumping window, default is 10000. (vector of length = length(exposure_data); or
#'   1, in which case it will be recycled)
#' @slot clump_r2 Clumping r2 cutoff.  (vector of length = length(exposure_data); or 1, in which
#'   case it will be recycled)
#' @slot clump_p1 Clumping sig level for index SNPs, default is 1.  (vector of length =
#'   length(exposure_data); or 1, in which case it will be recycled)
#' @slot clump_p2 Clumping sig level for secondary SNPs, default is 1.  (vector of length =
#'   length(exposure_data); or 1, in which case it will be recycled)
#' @slot pop Super-population to use as reference panel. Default = "EUR". Options are EUR, SAS, EAS,
#'   AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel with a
#'   slightly different set of markers
#'
#' @importFrom methods new
#'
#' @return a GenomeRegion object
#' @export
#'
GenomeRegion <- setClass(

  # The class name
  Class = "GenomeRegion",

  # Class variables
  representation(
    genome_wide = "logical",
    name        = "character",
    assembly    = "character",
    chromosome  = "numeric",
    start       = "numeric",
    end         = "numeric",
    cis_tol     = "numeric",
    clump_r2    = "numeric",
    clump_kb    = "numeric",
    clump_p1    = "numeric",
    clump_p2    = "numeric",
    pop         = "character"),

  # Constructor
  prototype(
    genome_wide = TRUE,
    name        = NA_character_,
    assembly    = "GRCh37",
    chromosome  = NA_real_,
    start       = NA_real_,
    end         = NA_real_,
    cis_tol     = 3e5,
    clump_r2    = 0.001,
    clump_kb    = 10000,
    clump_p1    = 1,
    clump_p2    = 1,
    pop         = "EUR")
)

setMethod("show",
          "GenomeRegion",
          function(object) {
            cat("Object of class 'GenomeRegion':\n")
            cat("-------------------------------\n")
            t <- tibble::tibble("Parameter"=character(), "Value"=character())
            if(object@genome_wide)
            {
                t <- tibble::add_row(t, Parameter="Genome-wide", Value="TRUE")
            }
            else
            {
                t <- tibble::add_row(t, Parameter="Genome-wide", Value="FALSE")
            }
            t <- tibble::add_row(t, Parameter="Name", Value=object@name)
            t <- tibble::add_row(t, Parameter="Assembly", Value=object@assembly)

            if(!object@genome_wide)
            {
                t <- tibble::add_row(t, Parameter="Gene region", Value="")
                t <- tibble::add_row(t, Parameter="Chromosome", Value=as.character(object@chromosome))
                t <- tibble::add_row(t, Parameter="Start", Value=as.character(object@start))
                t <- tibble::add_row(t, Parameter="End", Value=as.character(object@end))
                t <- tibble::add_row(t, Parameter="Cis-tolerance", Value=as.character(object@cis_tol))
            }
            else
            {
                t <- tibble::add_row(t, Parameter="Gene region", Value="")
                t <- tibble::add_row(t, Parameter="Chromosome", Value=NA_character_)
                t <- tibble::add_row(t, Parameter="Start", Value=NA_character_)
                t <- tibble::add_row(t, Parameter="End", Value=NA_character_)
                t <- tibble::add_row(t, Parameter="Cis-tolerance", Value=NA_character_)
            }
            t <- tibble::add_row(t, Parameter="Clumping params", Value="")
            t <- tibble::add_row(t, Parameter="R2", Value=as.character(object@clump_r2))
            t <- tibble::add_row(t, Parameter="kb", Value=as.character(object@clump_kb))
            t <- tibble::add_row(t, Parameter="p1", Value=as.character(object@clump_p1))
            t <- tibble::add_row(t, Parameter="p1", Value=as.character(object@clump_p2))
            t <- tibble::add_row(t, Parameter="p2", Value=as.character(object@pop))

            print(t, row.names = FALSE, dims = FALSE, classes = FALSE)
          }
)
