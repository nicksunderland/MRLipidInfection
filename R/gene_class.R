#' Gene class
#'
#' @slot name character. The name of the gene.
#' @slot assembly character. The assembly information, e.g. GRCh37
#' @slot chromosome numeric. The chromosome the gene is on.
#' @slot start numeric. The start (number of bases).
#' @slot end numeric. The end (number of bases).
#' @slot cis_tol when defining SNPs on a gene, set a tolerance of how many bases either side of the gene
#' @slot genome_wide this is usually FALSE, if TRUE this is treated as a object placeholder and will refer to the 'whole genome'
#'
#' @importFrom methods new
#'
#' @return a Gene object
#' @export
#'
Gene <- setClass(

  # The class name
  Class = "Gene",

  # Class variables
  representation(
    genome_wide = "logical",
    name  = "character",
    assembly = "character",
    chromosome = "numeric",
    start = "numeric",
    end = "numeric",
    cis_tol = "numeric"),

  # Constructor
  prototype(
    genome_wide = FALSE,
    name  = NA_character_,
    assembly = "GRCh37",
    chromosome = NA_real_,
    start = NA_real_,
    end = NA_real_,
    cis_tol = 3e5)
)
