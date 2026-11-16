#' @keywords internal
"_PACKAGE"

# mgcv is in Depends (attached, not just imported) so its formula symbols
# (s, te) and family constructors are visible to the user when writing
# structural terms and inspecting the returned model. The explicit importFrom
# below is what the package's own code binds against (and keeps R CMD check
# quiet about the Depends-but-not-imported namespace).
#' @importFrom mgcv gam bam nb tw betar scat concurvity gam.check
#' @importFrom stats AIC Gamma binomial gaussian poisson
#' @importFrom stats formula predict residuals
#' @importFrom graphics par plot
NULL
