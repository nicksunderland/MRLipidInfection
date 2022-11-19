#' extract_significant_snps
#'
#' @param filepath a file path to extract from
#' @param pval p-value filter (values less than this are kept)
#' @param pval_col the name of the column holding the p-value
#' @param snp_col the name of the column holding the snp id
#' @param overwrite whether or not to overwrite the cache .csv file
#' @param name the name of the exposure / outcome
#'
#' @importFrom dplyr mutate filter distinct .data
#' @importFrom data.table fread
#' @importFrom readr write_csv
#'
#' @return the data table; and attribute file path to the filtered .csv file
#' @export
#'
extract_significant_snps <- function(filepath, name=NULL, pval=5e-8, pval_col="pvalue", snp_col="SNP", overwrite=FALSE)
{
    message(paste0("Extracting significant SNPs from: ", basename(filepath)))

    # If csv_out name was not set, set to filename
    if(is.null(name))
    {
        name <- basename(filepath)
    }

    # The path to write to (place it alongside the main data file)
    fp <- paste0(dirname(filepath), "/sig_snps_", name, ".csv")

    # If file already exists
    if(file.exists(fp) & !overwrite)
    {
        warning("File already exists, change 'overwrite' parameter to TRUE if you wish to re-extract.\n")

        # Get the data
        d <- readr::read_csv(fp, show_col_types=F)
    }
    else
    {
        # Read the data from the file
        d <- data.table::fread(filepath) |>

            # Ensure that the pvalue is numeric
            dplyr::mutate(pvalue = as.numeric(.data[[pval_col]])) |>

            # Filter for SNPs with a pvalue <5e-8
            dplyr::filter(.data[[pval_col]] < pval) |>

            # Ensure unique SNPs
            dplyr::distinct(.data[[snp_col]], .keep_all=TRUE)

        # Write out the significant SNPs
        readr::write_csv(d, fp)
    }

    # Set the exposure name
    d$exposure  <- name
    d$phenotype <- name
    d$id.exposure <- name

    # Attached the filepath and name attribute
    attr(d, "filepath") <- fp
    attr(d, "name") <- name

    # Return the data
    return(d)
}
