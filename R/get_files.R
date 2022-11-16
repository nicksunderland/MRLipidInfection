#' get_web_file()
#'
#' Gets a web based file and downloads it into the packages downloads directory (inst/extdata)
#'
#' @param url url to file to download, a string
#' @param directory either NA (will download to main downloads directory), a full directory path
#'   (will download here), or a string (will be used to create a subdirectory within the main
#'   downloads directory)
#' @param overwrite whether or not to overwrite the file
#' @return a tibble containing the data
#' @importFrom httr GET
#' @importFrom utils tail download.file
#' @export
#'
get_web_file <- function(url, directory=NA, overwrite=FALSE)
{
    ## Parameter checks
    stopifnot("Overwrite parameter must be logical" = is.logical(overwrite))

    ## Get the directory path to download to
    dir_path <- NULL
    if(is.na(directory))
    {
        ## The base downloads directory
        dir_path <- get_downloads_dir()
    }
    else if(dir.exists(directory))
    {
        ## Some other valid directory that is given
        dir_path <- directory
    }
    else
    {
        ## A subd-irectory of the base downloads directory
        dir_path <- get_downloads_dir(subdir = directory)
    }

    ## Create the filename from the url
    fname <- utils::tail(strsplit(url, "/")[[1]],1)

    ## The path to the file
    fp <- paste0(dir_path, "/", fname)

    ## Download the data
    if(!file.exists(fp) | overwrite)
    {
        tryCatch(
          {
              utils::download.file(url=url, destfile=fp, method='curl')
              return(fp)
          },
          error=function(e){

              tryCatch(
                {
                    httr::GET(url, httr::user_agent("Chrome/102.0.5005.61"), httr::write_disk(fp, overwrite=overwrite))
                    return(fp)
                },
                error=function(e){
                    stop("Tried to download file with 'download.file()' and 'httr::GET()' methods, both failed.")
                    return(NULL)
                }
              )
          }
        )
    }
    else
    {
        message("File was already downloaded, if you want to overwrite it change the overwirte parameter to 'TRUE'")
        return(fp)
    }
}


#' get_finngen()
#'
#' Gets the FinnGen data - downloads into package downloads directory
#'
#' @param filename FinnGen file name, a string
#' @param type FinnGen data type
#' @param overwrite whether or not to overwrite the file
#' @return a tibble containing the data
#' @export
#'
get_finngen <- function(filename, type="summary_stats", overwrite=FALSE)
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

  ## Download the data
  fp <- get_web_file(url, directory="finngen", overwrite = overwrite)

  return(fp)
}

#' get_downloads_dir()
#'
#' @param subdir if provided with return the full file path to this subdirectory, within the main
#'   downloads directory. If it does not exist, it will be created.
#'
#' @return the absoltue file path to the downloads folder
#' @export
#'
get_downloads_dir <- function(subdir=NA){

  ## If asking for a sub-directory, create here
  sub <- ""
  if(!is.na(subdir))
  {
      sub <- paste0("/", subdir)
  }

  ## Check the extdata directory exists
  dl_dir <- paste0(system.file(package="MRLipidInfection"), "/extdata")
  if(!dir.exists(dl_dir))
  {
    message(paste0("Creating downloads directory at: ", dl_dir))
    dir.create(dl_dir)
  }

  ## The full directory
  dl_dir <- paste0(system.file(package="MRLipidInfection"), "/extdata", sub)
  if(!dir.exists(dl_dir))
  {
      message(paste0("Creating downloads directory at: ", dl_dir))
      dir.create(dl_dir)
  }

  ## Return the file path
  if(dir.exists(dl_dir))
  {
      return(dl_dir)
  }
  else
  {
      stop("Downloads directory does not exist and creating one failed")
  }
}

#' get_b_file()
#'
#' @param ancestry default is European 'EUR'
#'
#' @return the file path to the ld reference panel 'B' files
#' @export
#'
get_b_file <- function(ancestry="EUR"){

    fp  <- paste0(get_downloads_dir(), "/ld_reference_panel/", ancestry)
    fps <- paste0(fp, c(".bed", ".bim", ".fam"))

    if(!all(file.exists(fps)))
    {
        stop(paste0("Could not find b file on path: ", fp,
                    " Have you run 'setupMRLipidInfection()' at least once?"))
    }
    else
    {
        return(fp)
    }
}
