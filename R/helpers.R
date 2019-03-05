#' Show all available scripts
#'
#' @rdname show_scripts
#' @param category Name of type of analysis. (singlecell for now)
#' @export
#' @examples
#' show_scripts() # show categories
#' show_scripts("singlecell")
show_scripts = function(category=NULL){
    if (is.null(category))
        return(list.files(file.path(system.file("rmarkdown", package="hbcABC"), "Rscripts")))
    list.files(file.path(system.file("rmarkdown", package="hbcABC"),
              "Rscripts",
              category))
}

#' Copy a given script to other place
#'
#' @rdname copy_script
#' @param category Name of type of analysis.
#' @param fn Name of the file to copy.
#' @param dest Final folder to copy the file to.
#' @param dry Whether to run or not the copy action.
#' @export
#' @examples
#' copy_script("singlecell", "markers_seurat.R", dry=TRUE)
copy_script = function(category, fn, dest=NULL, dry=FALSE, ...){
    if (is.null(dest))
        dest = getwd()
    if (!file.exists(dest))
        dir.create(dest, recursive = TRUE)
    dest = file.path(dest, fn)
    origin =  file.path(system.file("rmarkdown", package="hbcABC"),
                        "Rscripts",
                        category,
                        fn)
    if (!file.exists(origin)){
        message("Use show_scripts to find accesible files")
        stop("This file doesn't exist ", origin)
    }
    if (!dry)
        file.copy(origin, dest, ...)
    message("copy ", fn, " to ", dest)
}
