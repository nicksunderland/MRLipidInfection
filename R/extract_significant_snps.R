#' extract_significant_snps
#'
#' @param filepath a file path to extract from
#' @param pval p-value filter (values less than this are kept)
#' @param pval_col the name of the column holding the p-value
#' @param snp_col the name of the column holding the snp id
#' @param csv_out a name to label the output .csv with, default is to use the filepath filename
#' @param overwrite whether or not to overwrite the cache .csv file
#'
#' @importFrom dplyr mutate filter distinct .data
#' @importFrom data.table fread
#' @importFrom readr write_csv
#'
#' @return the data table; and attribute file path to the filtered .csv file
#' @export
#'
extract_significant_snps <- function(filepath, csv_out, pval, pval_col, snp_col, overwrite=FALSE)
{
    message(paste0("Extracting significant SNPs from: ", basename(filepath)))

    # If csv_out name was not set, set to filename
    if(is.null(csv_out))
    {
        csv_out <- basename(filepath)
    }

    # The path to write to (place it alongside the main data file)
    fp <- paste0(dirname(filepath), "/sig_snps_", csv_out, ".csv")

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

    # Attached the filepath attribute
    attr(d, "filepath") <- fp

    # Return the data
    return(d)
}
