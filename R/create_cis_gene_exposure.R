#' create_cis_gene_exposure
#'
#' @param exposure_data exposure data.table, or list of, as returned by TwoSampleMR::format_data
#' @param gene a Gene object or list of Gene objects
#'
#' @return a list of data.tables, one for each exposure:gene combination
#' @export
#'
create_cis_gene_exposure <- function(exposure_data, gene){

    # The exposure name
    if(!is.null(names(exposure_data)))
    {
        exp_names <- names(exposure_data)
    }
    else if(!any(duplicated(purrr::map_chr(exposure_data, ~ unique(.x$exposure)))))
    {
        exp_names <- purrr::map_chr(exposure_data, ~ unique(.x$exposure))
    }
    else
    {
        exp_names <- paste0("exposure_", 1:length(exposure_data))
    }

    # The gene name(s)
    gene_names <- purrr::map_chr(gene, ~ .x@name)

    # The possible data combinations
    name_vec <- paste0(rep(gene_names, each=length(exposure_data)), "_", rep(exp_names, times=length(gene)))
    exp <- rep(exposure_data, times=length(gene))
    gen <- rep(gene, each=length(exposure_data))

    # Filter the exposures on the gene locations
    d <- purrr::pmap(.l= list(x=exp, y=gen, z=name_vec),
                     .f = function(x,y,z){

                       # Set the exposure name
                       x <- x |> dplyr::mutate(exposure = z)

                       # If not a genome_wide Gene object placeholder, do some SNP filtering
                       if(!y@genome_wide)
                       {
                          x <- x |>
                           dplyr::filter(.data$chr.exposure == y@chromosome,
                                         .data$pos.exposure > (y@start - y@cis_tol) &
                                         .data$pos.exposure < (y@end   + y@cis_tol))
                       }

                       return(x)
                      })

    # Name the combinations
    names(d) <- name_vec

    # Return the list
    return(d)
}

