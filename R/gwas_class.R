utils::globalVariables(".")

#' GWAS
#'
#' This class handles GWAS data making it easy to import external GWASs, manipulate them and quickly
#' retrieve the results by loaded cached .csv files. It uses TwoSampleMR::format_data s
#'
#' @slot type "exposure" or "outcome".
#' @slot file_path path to a gwas file that can be read by data.table::fread()
#' @slot config path to a yaml file with the column config information
#' @slot name a name for this GWAS
#' @slot sig_pval significance level at which to filter SNPs
#' @slot col_map list of .
#' @slot overwrite whether or not to overwrite the cached data files
#' @slot data a data.table object with the parse GWAS SNP data
#' @slot use_flag the SNPs to use in the analysis
#' @slot gen_reg a GenomeRegion object defining how to filter the SNP data
#'
#' @importFrom magrittr %>%
#' @importFrom yaml read_yaml
#' @importFrom RCurl url.exists
#' @importFrom urltools suffix_extract domain
#' @importFrom dplyr mutate across filter distinct
#' @importFrom rlang exec
#' @importFrom TwoSampleMR format_data
#' @importFrom data.table as.data.table fwrite fread
#' @importFrom methods callNextMethod validObject
#'
#' @return a GWAS object
#' @export
#'
GWAS <- setClass(

  # The class name
  Class = "GWAS",

  # Class variables
  representation(
    type      = "character",
    name      = "character",
    file_path = "character",
    config    = "character",
    overwrite = "logical",
    sig_pval  = "numeric",
    col_map   = "list",
    data      = "data.table",
    use_flag  = "logical",
    gen_reg   = "GenomeRegion"),

  # Constructor
  prototype(
    type      = NA_character_,
    file_path = NA_character_,
    config    = NA_character_,
    col_map   = list(),
    name      = NA_character_,
    sig_pval  = 5e-8,
    overwrite = FALSE,
    data      = data.table::data.table(),
    gen_reg   = GenomeRegion(),
    use_flag  = logical())
)

#' initialize
#'
#' @param .Object d
#' @param ... d
#'
#' @return a GWAS object
#'
setMethod(
  f          = "initialize",
  signature  = "GWAS",
  definition = function(.Object, ...)
  {
      # Calling a custom initialize so need to first call the standard one first to fill slots
      .Object <- callNextMethod(.Object, ...)

      # Read in the config file
      .Object@col_map <- tryCatch(yaml::read_yaml(.Object@config), error=function(e) stop("A valid config file must be provided"))

      # Extract the significant SNPs
      .Object <- extract_snps(.Object)

      #TODO: Other initialise functions

      # Make sure a valid object and return
      validObject(.Object)
      return(.Object)
  }
)


setValidity("GWAS", function(object)
{
    if(any(is.na(c(object@type, object@file_path, object@config, object@name))))
    {
        "Arguments 'type', 'file_path', 'config' and 'name' must be provided"
    }
    else
    {
        TRUE
    }

  #TODO: Other validation checks
  # colname checks after loading
  # length of use_flag
  # config file parsing
  # type either exposure or outcome

}
)

#' extract_snps
#'
#' @param object a GWAS objects
#'
#' @return GWAS object with SNPs extracted to data field
#' @export
#'
setGeneric("extract_snps", function(object) standardGeneric("extract_snps"))


#' Title
#'
#' @param object d
#'
#' @importFrom methods validObject
#' @importFrom dplyr mutate across filter distinct
#' @importFrom rlang :=
#'
#' @return d
#' @export
#'
setMethod(
  f          = "extract_snps",
  signature  = "GWAS",
  definition = function(object)
  {
      message(paste0("Extracting significant SNPs from: ", object@file_path))

      tryCatch(
      {
          # Try to download if path is a URL
          if(RCurl::url.exists(object@file_path))
          {
              # Split the URL into its parts
              url_parts <- urltools::suffix_extract(urltools::domain(object@file_path))

              # Download the file and set the filepath to it
              object@file_path <- get_web_file(object@file_path, directory=url_parts$domain)
          }

          # The path to write to (place it alongside the main data file)
          cached_snps_fp <- paste0(dirname(object@file_path), "/sig_snps_", object@name, ".csv")

          # Add the phenotype and id colnames to the mapping if it is not in the data
          phenotype_col_added <- FALSE
          id_col_added <- FALSE
          if(is.null(object@col_map$phenotype_col))
          {
              object@col_map <- c(object@col_map, list("phenotype_col" = "phenotype_col"))
              phenotype_col_added <- TRUE
          }

          if(is.null(object@col_map$id_col))
          {
              object@col_map <- c(object@col_map, list("id_col" = "id_col"))
              id_col_added <- TRUE
          }

          # If extracted SNPs file already exists return that
          if(file.exists(cached_snps_fp) & !object@overwrite)
          {
              warning("Extracted GWAS SNPs file already exists, change 'overwrite' parameter to TRUE if you wish to re-extract.\n")

              # Get the data
              tryCatch(
                expr = {
                    object@data <- data.table::fread(cached_snps_fp)
                },
                error = function(e){
                    stop(paste("Failed to extract from cached GWAS file:", cached_snps_fp))
                }
              )
          }
          else
          {
              # Read the data from the file path
              object@data <- data.table::fread(object@file_path) %>%

                # Make sure p value numeric
                dplyr::mutate(!!rlang::sym(object@col_map$pval_col) :=  as.numeric(!!rlang::sym(object@col_map$pval_col))) %>%

                # Filter out the SNPs with high pvalues
                dplyr::filter(!!rlang::sym(object@col_map$pval_col) < object@sig_pval) %>%

                # If no phenotype name specified set to object name
                {if(phenotype_col_added)
                 {
                      dplyr::mutate(., phenotype_col = object@name)
                 }
                 else
                 {
                      if(!object@col_map$phenotype_col %in% names(.))
                      {
                          stop("Phenotype column name provided in config file but not found in data")
                      }
                      else
                      {.}
                 }
                } %>%

                # If no id specified set to object name
                {if(id_col_added)
                 {
                      dplyr::mutate(., id_col = object@name)
                 }
                 else
                 {
                      if(!object@col_map$id_col %in% names(.))
                      {
                          stop("Phenotype column name provided in config file but not found in data")
                      }
                      else
                      {.}
                 }
                } %>%

                # Ensure unique SNPs
                dplyr::distinct(!!rlang::sym(object@col_map$snp_col), .keep_all=TRUE) %>%

                # Check if now empty - if so return empty df with colnames changed
                {if(nrow(.)==0)
                 {
                      dplyr::rename(., !!!object@col_map)
                 }
                 else
                 {
                      # Run the format_data(), do after filtering to make it quicker
                      rlang::exec(TwoSampleMR::format_data, dat=., !!!object@col_map)
                 }
                } %>%

                # Return as data.table
                data.table::as.data.table()

              # Write out the significant SNPs
              data.table::fwrite(object@data, cached_snps_fp)
          }
      },
      error = function(e) stop(paste("Failed to extract from GWAS file:", object@file_path, "\n", e))
      )

      # Set the GenomeRegion to itself (which has the side effect of creating the use_flag vector)
      genome_region(object) <- object@gen_reg

      validObject(object)
      return(object)
  }
)



#' set genome_region
#'
#' @param object a valid GWAS object
#' @param value a valid GenomeRegion object
#'
#' @importFrom methods callNextMethod validObject
#'
#' @return a valid GWAS object with GenomeRegion set
#' @export
#'
setGeneric("genome_region<-", function(object, value) standardGeneric("genome_region<-"))

#' set genome_region
#'
#' @param object a valid GWAS object
#' @param value a valid GenomeRegion object
#'
#' @importFrom methods callNextMethod validObject
#'
#' @return a valid GWAS object with GenomeRegion set
#' @export
#'
setMethod(
  f          = "genome_region<-",
  signature  = c("GWAS", "GenomeRegion"),
  definition = function(object, value)
  {
      # Set the GenomeRegion object and check valid
      object@gen_reg <- value
      validObject(object)

      # Set gene_col and id.exposure to "genome_wide"
      if(object@type == "exposure")
      {
          object@data$id.exposure <- object@gen_reg@name
          object@data$gene.exposure <- object@gen_reg@name
      }
      #TODO: do we need one for outcome?, if I expand the class to outcome handling

      # If not a genome_wide GenomeRegion object placeholder, do some SNP filtering
      if(!object@gen_reg@genome_wide)
      {
          # Logical flag as to which SNPs to use in the analysis
          object@use_flag <- object@data$chr.exposure == object@gen_reg@chromosome &
                             object@data$pos.exposure  > (object@gen_reg@start - object@gen_reg@cis_tol) &
                             object@data$pos.exposure  < (object@gen_reg@end   + object@gen_reg@cis_tol)
      }
      else
      {
          # Set use_flag to all TRUE
          object@use_flag <- rep(TRUE, nrow(object@data))
      }

      # Recheck valid and return
      validObject(object)
      return(object)
  }
)


#' clump_snps
#'
#' @param object a GWAS object
#' @param do_local run on local machine, if false does on IEU servers
#'
#' @importFrom methods validObject
#'
#' @return a GWAS object with clumped data
#' @export
#'
setGeneric("clump_snps", function(object, do_local=FALSE) standardGeneric("clump_snps"))



#' Title
#'
#' @param object d
#' @param do_local d
#'
#' @return d
#' @export
#'
setMethod(
  f          = "clump_snps",
  signature  = "GWAS",
  definition = function(object, do_local=FALSE)
    {
      # If no significant SNPs (either from pvalue filtering or GenomeRegion filtering)
      if(sum(object@use_flag)==0)
      {
        validObject(object)
        return(object)
      }

      # Where to do the processing (local or IEU servers)
      if(do_local)
      {
          tryCatch(
            expr = {
                # Clump the SNPs (locally)
                d <- ieugwasr::ld_clump(dat = dplyr::tibble(rsid=object@data$SNP[object@use_flag],
                                                            pval=object@data$pval.exposure[object@use_flag]),
                                        plink_bin = genetics.binaRies::get_plink_binary(),
                                        bfile     = get_b_file(),
                                        clump_kb  = object@gen_reg@clump_kb,
                                        clump_r2  = object@gen_reg@clump_r2,
                                        clump_p1  = object@gen_reg@clump_p1,
                                        clump_p2  = object@gen_reg@clump_p2,
                                        pop       = object@gen_reg@pop)

                # Set the flag to only those that are left after clumping
                object@use_flag <- object@data$SNP %in% d$SNP
            },
            error = function(e){
                stop(paste("Error clumping SNPs locally. Message:", e))
            }
          )
      }
      else
      {
          tryCatch(
            expr = {
              # Clump the SNPs (on IEU servers)
              d <- TwoSampleMR::clump_data(dat       = object@data[object@use_flag, ],
                                           clump_kb  = object@gen_reg@clump_kb,
                                           clump_r2  = object@gen_reg@clump_r2,
                                           clump_p1  = object@gen_reg@clump_p1,
                                           clump_p2  = object@gen_reg@clump_p2,
                                           pop       = object@gen_reg@pop)

              # Set the flag to only those that are left after clumping
              object@use_flag <- object@data$SNP %in% d$SNP
            },
            error = function(e){
              stop(paste("Error clumping SNPs remotely Message:", e))
            }
          )
      }

      # Check valid and return
      validObject(object)
      return(object)
    }
)


#' get_default_gwas_config
#'
#' opens the default config file and returns the path
#'
#' @param open_file whether to open file in R
#'
#' @importFrom utils file.edit
#'
#' @return file path
#' @export
#'
get_default_gwas_config <- function(open_file=FALSE)
{
    # Path to the default config file
    fp <- system.file("extdata/config/default_gwas.yml", package="MRLipidInfection")

    if(open_file)
    {
        # Open the file for editing
        utils::file.edit(fp)
    }

    # Return the path
    return(fp)
}
