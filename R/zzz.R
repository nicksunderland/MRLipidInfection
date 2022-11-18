#' setupMRLipidInfection
#'
#' @return side effect setup functions (see below)
#' @export
#'
setupMRLipidInfection <- function(){

    message("Checking MRLipidInfection package data...")

    ## Force create
    get_downloads_dir()

    ## Check and download the LD reference panel
    ld_ref_panel_data_check()

    ## Check PLINK installed
    check_local_plink()

    message("...checking complete")

}

ld_ref_panel_data_check <- function(){

    message("\tLD reference panel data:\t", appendLF=FALSE)

    ## Path to the downloads directory ld_reference_panel folder
    dl_dir <- get_downloads_dir(subdir = "ld_reference_panel")

    ## Number of files within the ld_reference_panel directory if it exists
    num_files <- length(list.files(dl_dir))

    ## If there are no ld_reference_panel files, download them
    ref <- tryCatch(
        {
            ## If no files, try to download them
            if(num_files == 0)
            {
                message("1000 Genomes LD Reference Panel required.")
                message("Starting download...")
                url <- "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
                fp  <- paste0(dl_dir, "/ld_reference_panel.tgz")
                f   <- download.file(url=url, destfile=fp, method='curl')

                ## download.file successfull
                if(f==0)
                {
                    fs = utils::untar(fp, exdir=dl_dir)

                    ## extraction successfull
                    if(f==0)
                    {
                        file.remove(fp)
                        message("LD Reference Panel download complete")
                    }
                }
            }
        },
        error=function(e)
        {
            stop("There was an error downloading and/or extracting the 1000 Genomes LD Reference Panel files")
        }
    )

    message("Done")
}

check_local_plink <- function(){

    message("\tLocal PLINK binary:\t\t", appendLF=FALSE)

    genetics.binaRies_installed <- requireNamespace("genetics.binaRies", quietly = TRUE)
    if(!genetics.binaRies_installed){
        devtools::install_github("explodecomputer/genetics.binaRies")
    }

    plink_bin <- genetics.binaRies::get_plink_binary()
    if(!file.exists(plink_bin))
    {
        warning(paste0("Failed to find PLINK binary at:\n", plink_bin))
    }
    else
    {
        message("Done")
    }


}

