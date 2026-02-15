

#' @title CI Width Estimation for Precision Analysis
#' @description A function to compute the expected widths of confidence
#' intervals for F.beta scores during the planning phase of a study.
#' @importFrom utils tail
#' @importFrom stats rbinom

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp Precision/PPV.
#' @param ps Sensitivity.
#' @param beta beta.
#' @param alpha One minus confidence level. Default of alpha = 0.05, i.e., 95 confidence interval.
#' @param method c("simulation", "fixed"). Suggest "simulation" for small sample sizes.
#' Default = "simulation".
#' @param normal.approx normal.approx = c("auto", "TRUE", "FALSE").
#' @param repliactes Simulation replicates when method = "simulation".
#' @param seed Simulation seed when method = "simulation".

#' @examples # widthCI(N=50, s=20, pp=9/11, ps=0.75)

#' @export
widthCI <- function(N, s, pp, ps,
                    beta=1,
                    alpha=0.05,
                    method=c("simulation", "fixed"),
                    normal.approx=c("auto", "TRUE", "FALSE"),
                    repliactes=5000,
                    seed=2026) {

  method <- match.arg(method)
  normal.approx <- match.arg(normal.approx)

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb <- one.plus.beta2*pp*ps/(ps + pp*beta2)

  F1.score <- F1.cond.b(N, s, pp, ps, beta)
  f1s <- round(F1.score$f1s, 10)
  pf1 <- F1.score$pf1

  if (method=="fixed") {

    C.obs <- alpha/2

    cumsum.pf1 <- cumsum(pf1)
    cumsum.rev.pf1 <- cumsum(rev(pf1))

    Fb.L <- ifelse(any(cumsum.pf1<=C.obs), tail(f1s[cumsum.pf1<=C.obs],1), 0)

    Fb.U <- ifelse(any(cumsum.rev.pf1<=C.obs), tail(rev(f1s)[cumsum.rev.pf1<=C.obs],1), 1)

    Fb.vec <- c(fb, Fb.L, Fb.U, Fb.U - Fb.L)

  } else {

    f <- function(N, s, pp, ps, beta, alpha, repliactes, seed) {
      pa <- paCalculation(N=N, s=s, pp=pp, ps=ps, beta=beta)$pa
      test.Hsu <- function() {
        b <- rbinom(1, size=N-s, prob=pa)
        d <- rbinom(1, size=s, prob=ps)
        out <- ciFbOne(d=d, b=b, N=N, s=s, alpha=alpha, normal.approx=normal.approx)
        out[2:3]
      }

      set.seed(seed)
      ci <- replicate(n=repliactes, test.Hsu())
      Fb.L <- mean( apply(ci, 2, function(x) x[1] ) ) # average length of lower CI
      Fb.U <- mean( apply(ci, 2, function(x) x[2] ) ) # average length of upper CI
      width <- Fb.U - Fb.L #mean( apply(ci, 2, function(x) x[2]-x[1] ) ) # average length of CI
      c(Fb.L, Fb.U, width)
    }

    Fb.vec <- c(fb, f(N, s, pp, ps, beta, alpha, repliactes, seed))

  }

  names(Fb.vec) <- c("Fb", paste0((1-alpha)*100,"%CI", c(".L", ".U")), "Width")

  return(Fb.vec)

}
