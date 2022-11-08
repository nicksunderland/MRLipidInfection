#' Gets the FinnGen data - downloads into inst/extdata
#'
#' @param filename FinnGen file name, a string
#' @param type FinnGen data type
#' @param overwrite whether or not to overwrite the file
#' @return a tibble containing the data
#' @importFrom vroom vroom
#' @importFrom utils download.file
#' @export
#'
getFinnGen <- function(filename, type="summary_stats", overwrite=FALSE)
{
    ## Parameter checks
    stopifnot("Type must be one of 'summary_stats', 'finemapping', 'annotations', or 'covid'" =
                type %in% c('summary_stats', 'finemapping', 'annotations', 'covid'),
              "Filename must end in .gz" = grepl(".*.gz$", filename),
              "Overwrite parameter must be logical" = is.logical(overwrite))

    ## FINNGEN website base URL
    url <- "https://storage.googleapis.com/finngen-public-data-r7/"

    ## FinnGen data category
    url <- paste0(url, type)

    ## Filepath to data
    fp <- paste0(system.file("extdata", "finngen_data", package = "MRLipidInfection"), "/", filename)

    ## Download the data
    if(!file.exists(fp) | overwrite)
    {
        options(timeout=600)
        f <- download.file(url = paste0(url, "/", filename), destfile = fp)
    }
    else
    {
        f <- 0
        message("File was already downloaded, if you want to overwrite it change the overwirte parameter to 'TRUE'")
    }

    ## Download successful or data already downloaded
    if(f==0)
    {
        t <- vroom::vroom(fp, show_col_types=FALSE)
        return(t)
    }
}
