#' Random starts with m variables
#'
#' @param m order of interaction to include.
#' @param nvars number of varialbes that could be included in the interaction.
#' @param quad whether to include quadratic terms or not.
#' @param numrs number of random starts to compute.
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#'
#' @return A matrix of random starts by row
#' @export
#'
#' @examples
#' rstart(2,5,quad=FALSE,numrs=3) #doesn't include quadratic terms
#' rstart(2,5,quad=TRUE,numrs=3) #does include quadratic terms
rstart = function(m,nvars,quad = FALSE,numrs = 1) {
  t(replicate(numrs,sample(
    1:nvars,size = m,replace = quad
  )))
}
