#' get_proxy_snps
#'
#' @param snps a character list or vector of SNPs
#' @param ancestry the ancestry to use default "EUR"
#' @param ... additional paramters to pass to gwasvcf::get_ld_proxies
#'
#' @importFrom gwasvcf check_plink set_plink get_ld_proxies
#' @importFrom genetics.binaRies get_plink_binary
#'
#' @return a dataframe
#' @export
get_proxy_snps <- function(snps, ancestry="EUR", ...){

    ## Warning as they occur
    options(warn=1)

    ## Ensure that we know where the PLINK binary is
    if(!suppressMessages(gwasvcf::check_plink()))
    {
        plink_path <- genetics.binaRies::get_plink_binary()

        gwasvcf::set_plink(plink_path)

        if(!gwasvcf::check_plink()) stop("Error setting path to PLINK binary")
    }

    ## Use only unqiue SNPs
    snps_len <- length(snps)
    snps     <- unique(snps)
    if(snps_len != length(snps))
    {
        warning(paste0(as.character(snps_len-length(snps)), " duplicate SNPs were removed."))
    }

    ## Get the path to the reference population (only 1000 genomes supported at present)
    ref_path <- get_b_file()

    ## Run the LD proxy search
    proxies <- gwasvcf::get_ld_proxies(snps, ref_path, ...)

    return(proxies)
}
