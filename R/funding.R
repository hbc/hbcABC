#' Return acknowledgement section for analysis done by HBC
#'
#' @rdname funding
#' @param author Name of the author/s to be added to the text.
#' @param service Name of the analysis done by the core.
#' @export
#' @examples
#' funding_general("Lorena Pantnao", "RNASeq...")
#' funding_catalyst()
funding_general = function(author, service){
    paste("The authors would like to thank", author , "of the Harvard Chan Bioinformatics Core, Harvard T.H. Chan School of Public Health, Boston, MA for assistance with", service,".")
}
#' @rdname funding
#' @export
funding_catalyst = function(){
    paste("The project described was conducted with the support of Harvard Catalyst | The Harvard Clinical and Translational Science Center (NIH award #UL1 RR 025758 and financial contributions from participating institutions). The content is solely the responsibility of the authors and does not
necessarily represent the official views of the National Center for Research Resources or theNational Institutes of Health.")
}
#' @rdname funding
#' @export
funding_cfar = function(author){
    paste("Work contributed by ", author, "to this abstract/publication/presentation/grant proposal was made possible with help from the Harvard University Center for AIDS Research (CFAR), an NIH funded program (P30 AI060354), which is supported by the following NIH Co-Funding and Participating Institutes and Centers: NIAID, NCI, NICHD, NHLBI, NIDA, NIMH, NIA, NIDDK, NIGMS, FIC, and OAR. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.")
}
#' @rdname funding
#' @export
funding_hsci = function(author){
    paste("Work by", author, "at the Harvard Chan
Bioinformatics Core was funded by the \"HSCI Center for Stem Cell Bioinformatics\"")
}
#' @rdname funding
#' @export
funding_tnt = function(author){
    paste("Work by", author, "at the Harvard Chan Bioinformatics Core was funded by the Harvard Medical School Tools and Technology Committee")
}
#' @rdname funding
#' @export
funding_niehs = function(author){
   paste("Work by", author, "at the Harvard Chan Bioinformatics Core was supported by Award Number P30 ES00002 from the National Institute of Environmental Health Sciences. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institute of Environmental Health Sciences or the National Institutes of Health.")
}
