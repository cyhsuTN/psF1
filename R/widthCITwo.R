

#' @title Computation of CI Widths for Comparative F.beta Scores in Precision Analysis
#' @description A function to compute the expected width of a confidence
#' interval for differences in F.beta scores during the study planning stage.
#' @importFrom utils tail
#' @importFrom stats rbinom

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param N2 Number of classification instances in 2nd independent dataset
#' used for 2nd classifier. N2 = NULL denotes a common dataset is used to
#' compare both classifiers.
#' @param s2 Number of positives in 2nd independent dataset.
#' @param beta beta.
#' @param alpha One minus confidence level. Default of alpha = 0.05, i.e., 95 confidence interval.
#' @param method c("analytic", "simulation"). Default = "analytic".
#' @param normal.approx normal.approx = c("auto", "TRUE", "FALSE").
#' @param replicates Simulation replicates when method = "simulation".
#' @param seed Simulation seed when method = "simulation".

#' @examples # widthCI.Two(N=50, s=20, pp1=0.6, ps1=0.8, pp2=0.7, ps2=0.9)

#' @export
widthCI.Two <- function(N, s,
                        pp1, ps1,
                        pp2, ps2,
                        N2 = NULL, s2 = NULL,
                        beta = 1,
                        alpha = 0.05,
                        method = c("analytic", "simulation"),
                        normal.approx=c("auto", "TRUE", "FALSE"),
                        replicates = 1000,
                        seed = 2026) {

  method <- match.arg(method)
  normal.approx <- match.arg(normal.approx)

  if(is.null(N2)) {N2 <- N; s2 <- s}
  if(is.null(s2)) s2 <- s

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb1 <- one.plus.beta2*pp1*ps1/(ps1 + pp1*beta2)
  fb2 <- one.plus.beta2*pp2*ps2/(ps2 + pp2*beta2)

  dfb <- fb1 - fb2

  F1dist1 <- F1.cond.b(N=N, s=s, pp=pp1, ps=ps1, beta)
  F1dist2 <- F1.cond.b(N=N2, s=s2, pp=pp2, ps=ps2, beta)

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

  if (method=="analytic") {

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

    dFb.vec <- c(dfb, dFb.L, dFb.U, dFb.U - dFb.L)

  } else {

    f <- function(N, s, pp1, ps1,
                  N2, s2, pp2, ps2,
                  beta, alpha, replicates, seed) {

      pa1 <- paCalculation(N=N, s=s, pp=pp1, ps=ps1, beta=beta)$pa
      pa2 <- paCalculation(N=N2, s=s2, pp=pp2, ps=ps2, beta=beta)$pa

      test.Hsu <- function() {
        b1 <- rbinom(1, size=N-s, prob=pa1)
        d1 <- rbinom(1, size=s, prob=ps1)
        b2 <- rbinom(1, size=N2-s2, prob=pa2)
        d2 <- rbinom(1, size=s2, prob=ps2)

        out <- ciFbTwo(d1=d1, b1=b1,
                       d2=d2, b2=b2,
                       N=N, s=s,
                       N2=N2, s2=s2,
                       alpha=0.05, normal.approx="auto")
        out[2:3]
      }

      set.seed(seed)
      ci <- replicate(n=replicates, test.Hsu())
      dFb.L <- mean( apply(ci, 2, function(x) x[1] ) ) # average length of lower CI
      dFb.U <- mean( apply(ci, 2, function(x) x[2] ) ) # average length of upper CI
      width <- dFb.U - dFb.L #mean( apply(ci, 2, function(x) x[2]-x[1] ) ) # average length of CI
      c(dFb.L, dFb.U, width)
    }

    dFb.vec <- c(dfb, f(N, s, pp1, ps1,
                      N2, s2, pp2, ps2,
                      beta, alpha, replicates, seed))

  }

  dFb.vec <- c(dFb.vec, tail(dFb.vec,1)/2)
  names(dFb.vec) <- c("dFb", paste0((1-alpha)*100,"%CI", c(".L", ".U")),
                      "Width", "Half.width")

  return(dFb.vec)

}
