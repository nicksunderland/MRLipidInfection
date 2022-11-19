#' clump_snps
#'
#' @param exposure_data a data.table consistent with output from TwoSampleMR::format_data (or named list of them)
#' @param do_local whether to do the processing on the local computer or IEU servers
#' @param clump_kb Clumping window, default is 10000.
#' @param clump_r2 Clumping r2 cutoff.
#' @param clump_p1 Clumping sig level for index SNPs, default is 1
#' @param clump_p2 Clumping sig level for secondary SNPs, default is 1.
#' @param popSuper-population to use as reference panel. Default = "EUR". Options are EUR, SAS, EAS,
#'   AFR, AMR. 'legacy' also available - which is a previously used verison of the EUR panel with a
#'   slightly different set of markers
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

    if(do_local)
    {
        d1 <- purrr::map(.x = exposure_data,
                         .f = function(x, ...){

                            # If there are no SNPs just return the empty data.table
                            if(nrow(x)==0) return(x)

                            # If do_local is TRUE then perform analysis on local machine
                            d2 <- ieugwasr::ld_clump(dat = dplyr::tibble(rsid=x$SNP,
                                                                        pval=x$pval.exposure),
                                                    plink_bin = genetics.binaRies::get_plink_binary(),
                                                    bfile     = get_b_file(),
                                                    ...)
                            return(d2)
                          }, clump_kb, clump_r2, clump_p1, clump_p2, pop)

    }
    else
    {
        d1 <- purrr::map(.x = exposure_data,
                         .f = function(x, ...){

                           # If there are no SNPs just return the empty data.table
                           if(nrow(x)==0) return(x)

                           # If do_local is FALSE then perform analysis on IEU servers
                           TwoSampleMR::clump_data(dat = x, ...)

                         }, clump_kb, clump_r2, clump_p1, clump_p2, pop)
    }

    # Return the list of data.tables
    return(d1)
}
