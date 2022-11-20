#' clump_snps
#'
#' @param exposure_data a data.table consistent with output from TwoSampleMR::format_data; must be a
#'   named list
#' @param do_local whether to do the processing on the local computer or IEU servers
#' @param clump_kb Clumping window, default is 10000. (vector of length = length(exposure_data); or
#'   1, in which case it will be recycled)
#' @param clump_r2 Clumping r2 cutoff.  (vector of length = length(exposure_data); or 1, in which
#'   case it will be recycled)
#' @param clump_p1 Clumping sig level for index SNPs, default is 1.  (vector of length =
#'   length(exposure_data); or 1, in which case it will be recycled)
#' @param clump_p2 Clumping sig level for secondary SNPs, default is 1.  (vector of length =
#'   length(exposure_data); or 1, in which case it will be recycled)
#' @param pop Super-population to use as reference panel. Default = "EUR". Options are EUR, SAS,
#'   EAS, AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel
#'   with a slightly different set of markers
#'
#' @return clumped data.table, or list of them
#' @export
#'
clump_snps <- function(exposure_data,
                       do_local=FALSE,
                       clump_kb = 10000,
                       clump_r2 = 0.001,
                       clump_p1 = 1,
                       clump_p2 = 1,
                       pop = "EUR"){

    # The length of all the clumping parameters should be either 1 or length(exposure_data)
    if(is.data.frame(exposure_data)){
        n <- 1
    }else if(is.list(exposure_data)){
        n <- length(exposure_data)
    }else{
        stop("exposure_data input must be a data.frame or list of data.frames")
    }

    # Make a list of the clumping arguments
    clump_args <- list("clump_kb"=clump_kb,
                       "clump_r2"=clump_r2,
                       "clump_p1"=clump_p1,
                       "clump_p2"=clump_p2,
                       "pop"     =pop)

    # Make sure there are either 1 or n replications of each argument
    clump_args <- purrr::map2(.x = clump_args,
                              .y = purrr::map(clump_args, ~length(.x)),
                              .f = function(arg, len, n){
                                  if(len==n)
                                  {
                                      return(arg)
                                  }
                                  else if(len==1)
                                  {
                                      return(rep(arg, n))
                                  }
                                  else
                                  {
                                      stop("argumnents must be length 1 or length(exposure_data)")
                                  }
                                }, n
                               )

    # Where to do the processing (local or IEU servers)
    if(do_local)
    {
        d1 <- purrr::pmap(.l = list("x"         = exposure_data,
                                    "clump_kb_" = clump_args$clump_kb,
                                    "clump_r2_" = clump_args$clump_r2,
                                    "clump_p1_" = clump_args$clump_p1,
                                    "clump_p2_" = clump_args$clump_p2,
                                    "pop_"      = clump_args$pop),
                          .f = function(x, clump_kb_, clump_r2_, clump_p1_, clump_p2_, pop_){

                            # If there are no SNPs just return the empty data.table
                            if(nrow(x)==0) return(x)

                            # If do_local is TRUE then perform analysis on local machine
                            d2 <- ieugwasr::ld_clump(dat = dplyr::tibble(rsid=x$SNP,
                                                                         pval=x$pval.exposure),
                                                    plink_bin = genetics.binaRies::get_plink_binary(),
                                                    bfile     = get_b_file(),
                                                    clump_kb  = clump_kb_,
                                                    clump_r2  = clump_r2_,
                                                    clump_p1  = clump_p1_,
                                                    clump_p2  = clump_p2_,
                                                    pop       = pop_)
                            return(d2)
                          })

    }
    else
    {
        d1 <- purrr::pmap(.l = list("x"         = exposure_data,
                                    "clump_kb_" = clump_args$clump_kb,
                                    "clump_r2_" = clump_args$clump_r2,
                                    "clump_p1_" = clump_args$clump_p1,
                                    "clump_p2_" = clump_args$clump_p2,
                                    "pop_"      = clump_args$pop),
                          .f = function(x, clump_kb_, clump_r2_, clump_p1_, clump_p2_, pop_){

                            # If there are no SNPs just return the empty data.table
                            if(nrow(x)==0) return(x)

                            # If do_local is FALSE then perform analysis on IEU servers
                            d2 <- TwoSampleMR::clump_data(dat = x,
                                                          clump_kb  = clump_kb_,
                                                          clump_r2  = clump_r2_,
                                                          clump_p1  = clump_p1_,
                                                          clump_p2  = clump_p2_,
                                                          pop       = pop_)
                            return(d2)
                          })
    }

    # Return the list of data.tables
    return(d1)
}
