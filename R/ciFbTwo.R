
#' @title Confidence Interval Estimate for F.beta Difference
#' @description A function to calculate the CI of the F.beta difference between
#' two classifiers.
#' @importFrom utils tail
#' @importFrom stats qnorm

#' @param d1 Number of true positives for 1st classifier.
#' @param b1 Number of false positives for 1st classifier.
#' @param d2 Number of true positives for 2nd classifier.
#' @param b2 Number of false positives for 2nd classifier.
#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param beta beta.
#' @param alpha Significance level.
#' @param normal.approx normal.approx = TRUE or FALSE.

#' @examples # ciFbTwo(d1=9, b1=2, d2=10, b2=3, N=50, s=20)

#' @export
ciFbTwo <- function(d1=9, b1=2,
                    d2=9, b2=2,
                    N=50, s=20,
                    beta=1,
                    alpha=0.05,
                    normal.approx=TRUE) {


  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if (b1==0) b1 <- 0.5
  if (b2==0) b2 <- 0.5

  pp1.hat <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
  ps1.hat <- (d1+0.5)/(s+1)

  pp2.hat <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
  ps2.hat <- (d2+0.5)/(s+1)

  #pp1.hat <- d1/(d1+b1)
  #ps1.hat <- d1/s

  #pp2.hat <- d2/(d2+b2)
  #ps2.hat <- d2/s

  fb1.raw <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2.raw <- one.plus.beta2*d2/(d2 + b2 + s*beta2)
  fb1 <- one.plus.beta2*pp1.hat*ps1.hat/(ps1.hat + pp1.hat*beta2)
  fb2 <- one.plus.beta2*pp2.hat*ps2.hat/(ps2.hat + pp2.hat*beta2)

  dfb.raw <- fb1.raw - fb2.raw
  dfb <- fb1 - fb2

  F1dist1 <- F1.cond.b(N=N, s=s, pp=pp1.hat, ps=ps1.hat, beta)
  F1dist2 <- F1.cond.b(N=N, s=s, pp=pp2.hat, ps=ps2.hat, beta)

  x1 <- 1:nrow(F1dist1)
  x2 <- 1:nrow(F1dist2)
  xf11 <- F1dist1$f1s
  xf12 <- F1dist2$f1s
  yf11 <- F1dist1$pf1
  yf12 <- F1dist2$pf1

  idx.d <- (yf11 > 1E-6)
  idx.b <- (yf12 > 1E-6)
  x1 <- x1[idx.d]
  x2 <- x2[idx.b]

  if( (length(xf11) > 1E+3 | length(xf12) > 1E+3) & normal.approx ) {

    mean1 <- sum(xf11 * yf11)
    var1 <- sum(xf11^2 * yf11) - mean1^2

    mean2 <- sum(xf12 * yf12)
    var2 <- sum(xf12^2 * yf12) - mean2^2

    dmean <- mean1 - mean2
    dvar <- (var1 + var2)

    dFb.L <- dmean - qnorm(1-alpha/2)*sqrt(dvar)

    dFb.U <- dmean + qnorm(1-alpha/2)*sqrt(dvar)

    dFb.vec <- c(dfb, dFb.L, dFb.U, dfb.raw)

  } else {

    val1 <- outer(x1, x2, function(d, b) xf11[d] - xf12[b])
    val2 <- outer(x1, x2, function(d, b) yf11[d] * yf12[b])

    pf1 <- tapply(as.vector(val2), as.vector(val1), sum)
    pf1 <- pf1[pf1>1E-10]
    pf1 <- pf1/sum(pf1)
    f1s <- round(as.numeric(names(pf1)), 10)
    pf1 <- as.numeric(pf1)

    C.obs <- alpha/2

    cumsum.pf1 <- cumsum(pf1)
    cumsum.rev.pf1 <- cumsum(rev(pf1))

    dFb.L <- ifelse(any(cumsum.pf1<=C.obs), tail(f1s[cumsum.pf1<=C.obs],1), 0)

    dFb.U <- ifelse(any(cumsum.rev.pf1<=C.obs), tail(rev(f1s)[cumsum.rev.pf1<=C.obs],1), 1)

    dFb.vec <- c(dfb, dFb.L, dFb.U, dfb.raw)

  }

  names(dFb.vec) <- c("dFb.est", paste0((1-alpha)*100,"%CI", c(".L", ".U")), "dFb.est.raw")

  return(dFb.vec)

}
