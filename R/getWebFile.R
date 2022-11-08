#' Gets a web based file - downloads into inst/extdata
#'
#' @param url url to file to download, a string
#' @param overwrite whether or not to overwrite the file
#' @return a tibble containing the data
#' @importFrom vroom vroom
#' @importFrom fs path_sanitize
#' @importFrom httr GET
#' @export
#'
getWebFile <- function(url, overwrite=FALSE)
{
    ## Parameter checks
    stopifnot("Overwrite parameter must be logical" = is.logical(overwrite))

    ## Filepath to data
    fp <- paste0(system.file("extdata", "other_files", package = "MRLipidInfection"), "/", sub("/", "", fs::path_sanitize(url)))

    ## Download the data
    if(!file.exists(fp) | overwrite)
    {
        try({
              httr::GET(url, httr::user_agent("Chrome/102.0.5005.61"), httr::write_disk(fp, overwrite=overwrite))
              f <- 0
            })
    }
    else
    {
        f <- 0
        message("File was already downloaded, if you want to overwrite it change the overwirte parameter to 'TRUE'")
    }

    ## Download successful or data already downloaded
    if(f==0)
    {
        return(fp)
    }
}
