
#' @title Distribution of the Difference between Two F Beta Scores
#' @description A function to show a distribution of the difference between two F beta scores.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param beta beta.

#' @examples # fbscore <- F1.cond.two.b(N=50, s=20, pp1=9/11, ps1=0.75, beta=1)
#' @examples # plotF1(fbscore, type=c("Density", "CDF")[1])

#' @export
F1.cond.two.b <- function(N, s,
                          pp1=9/11, ps1=0.75,
                          pp2=NULL, ps2=NULL,
                          beta=1) {

  if(is.null(pp2) | is.null(ps2)) {
    F1dist0 <- F1dist1 <- F1.cond.b(N, s, pp1, ps1, beta=beta)
  } else {
    F1dist0 <- F1.cond.b(N, s, pp1, ps1, beta)
    F1dist1 <- F1.cond.b(N, s, pp2, ps2, beta)
  }

  d <- 1:nrow(F1dist0)
  b <- 1:nrow(F1dist1)
  xf10 <- F1dist0$f1s
  yf10 <- F1dist0$pf1

  xf11 <- F1dist1$f1s
  yf11 <- F1dist1$pf1


  idx.d <- (yf10 > 1E-6)
  idx.b <- (yf11 > 1E-6)

  d <- d[idx.d]
  b <- b[idx.b]

  val1 <- outer(d, b, function(d, b) xf10[d] - xf11[b])
  val2 <- outer(d, b, function(d, b) yf10[d] * yf11[b])

  pf1 <- tapply(as.vector(val2), as.vector(val1), sum)
  f1s <- as.numeric(names(pf1))
  return(data.frame(f1s=f1s, pf1=as.numeric(pf1)))

}
