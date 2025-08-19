
#' @title Confidence Interval Estimate for F.beta of One Classifier
#' @description A function to calculate the CI of the F.beta of one classifier.
#' @importFrom utils tail
#' @importFrom stats qnorm

#' @param d Number of true positives.
#' @param b Number of false positives.
#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param beta beta.
#' @param alpha Significance level.
#' @param normal.approx normal.approx = TRUE or FALSE.

#' @examples # ciFbOne(d=9, b=2, N=50, s=20)

#' @export
ciFbOne <- function(d=9, b=2, N=50, s=20,
                    beta=1,
                    alpha=0.05,
                    normal.approx=TRUE) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if (b==0) b <- 0.5

  pp.hat <- (d+0.5)/(d+b+1) # Jeffrey’s non-informative prior
  ps.hat <- (d+0.5)/(s+1)

  #pp.hat <- d/(d+b)
  #ps.hat <- d/s

  fb.raw <- one.plus.beta2*d/(d + b + s*beta2)
  fb <- one.plus.beta2*pp.hat*ps.hat/(ps.hat + pp.hat*beta2)

  F1.score <- F1.cond.b(N, s, pp.hat, ps.hat, beta)
  f1s <- round(F1.score$f1s, 10)
  pf1 <- F1.score$pf1


  if( length(f1s) > 1E+5 & normal.approx ) {

    mean1 <- sum(f1s * pf1)
    var1 <- sum(f1s^2 * pf1) - mean1^2

    Fb.L <- max(mean1 - qnorm(1-alpha/2)*sqrt(var1), 0)

    Fb.U <- min(mean1 + qnorm(1-alpha/2)*sqrt(var1), 1)

    Fb.vec <- c(fb, Fb.L, Fb.U, fb.raw)

  } else {

    C.obs <- alpha/2

    cumsum.pf1 <- cumsum(pf1)
    cumsum.rev.pf1 <- cumsum(rev(pf1))

    Fb.L <- ifelse(any(cumsum.pf1<=C.obs), tail(f1s[cumsum.pf1<=C.obs],1), 0)

    Fb.U <- ifelse(any(cumsum.rev.pf1<=C.obs), tail(rev(f1s)[cumsum.rev.pf1<=C.obs],1), 1)

    Fb.vec <- c(fb, Fb.L, Fb.U, fb.raw)

  }

  names(Fb.vec) <- c("Fb.est", paste0((1-alpha)*100,"%CI", c(".L", ".U")), "Fb.est.raw")

  return(Fb.vec)

}

