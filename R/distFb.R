
#' @title Distribution of F Beta Score
#' @description A function to show a distribution of F beta score.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp Precision.
#' @param ps Sensitivity.
#' @param beta beta.

#' @examples # fbscore <- F1.cond.b(N=50, s=20, pp=9/11, ps=0.75, beta=1)
#' @examples # plotF1(fbscore, type=c("Density", "CDF")[1])

#' @export
F1.cond.b <- function(N, s, pp=9/11, ps=0.75, beta=1) {

  #pa <- paCalculation_v1(N, s, pp, ps)$pa
  pa <- paCalculation(N, s, pp, ps, beta)$pa

  d <- 0:s
  b <- 0:(N-s)
  d_probs <- dbinom(d, s, prob=ps)
  b_probs <- dbinom(b, N-s, prob=pa)

  idx.d <- (d_probs > 1E-10)
  idx.b <- (b_probs > 1E-10)

  d <- d[idx.d]
  b <- b[idx.b]

  one.plus.beta2 <- (1+beta^2)
  s.beta2 <- s*beta^2

  val1 <- outer(d, b, function(d, b) ifelse(d == 0, 0, one.plus.beta2 * d / (d + b + s.beta2)))
  val2 <- outer(d, b, function(d, b) d_probs[d + 1] * b_probs[b + 1])

  pf1 <- tapply(as.vector(val2), as.vector(val1), sum)
  f1s <- as.numeric(names(pf1))
  return(data.frame(f1s=f1s, pf1=as.numeric(pf1) ))

}

