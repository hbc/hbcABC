#' Show all available scripts
#'
#' @rdname show_scripts
#' @param category Name of type of analysis. (singlecell for now)
#' @export
show_scripts = function(category="singlecell"){
    list.files(file.path(system.file("rmarkdown", package="hbcABC"),
              "Rscripts",
              category))
}

#' Copy a given script to other place
#'
#' @rdname copy_script
#' @param category Name of type of analysis.
#' @param fn Name of the file to copy.
#' @export
copy_script = function(category, fn, dest=NULL, ...){
    if (is.null(dest))
        dest = getwd()
    dest = file.path(dest, fn)
    file.copy(file.path(system.file("rmarkdown", package="hbcABC"),
                        "Rscripts",
                        category,
                        fn),
              dest, ...)
}
