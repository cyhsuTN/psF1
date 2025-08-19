
#' @title Power Calculation for F.beta Comparison Between Two Classifiers
#' @description A function to calculate power for testing the difference in F.beta between two classifiers.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot pnorm qnorm qt


#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param Fb1 F.beta score of 1st classifier.
#' @param pp1 Precision of 1st classifier.
#' @param ps1 Sensitivity of 1st classifier.
#' @param Fb2 F.beta score of 2nd classifier.
#' @param pp2 Precision of 2nd classifier.
#' @param ps2 Sensitivity of 2nd classifier.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = TRUE or FALSE.
#' @param alpha.level alpha level.

#' @examples # powerFbTwo(N=50, s=35,
#' @examples #            Fb1=0.783, pp1=NULL, ps1=NULL,
#' @examples #            Fb2=0.95, pp2=NULL, ps2=NULL,
#' @examples #            beta=1,
#' @examples #            alternative=c("two.sided", "less", "greater")[1],
#' @examples #            normal.approx=TRUE,
#' @examples #            alpha.level=0.05)

#' @export
powerFbTwo <- function(N=50, s=20,
                       Fb1=0.783, pp1=NULL, ps1=NULL,
                       Fb2=0.9, pp2=NULL, ps2=NULL,
                       beta=1,
                       alternative = c("two.sided", "less", "greater")[1],
                       normal.approx=TRUE,
                       alpha.level=0.05) {


  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb1)) Fb1 <- one.plus.beta2*pp1*ps1/(beta2*pp1+ps1)
  if(is.null(Fb2)) Fb2 <- one.plus.beta2*pp2*ps2/(beta2*pp2+ps2)

  Fb1.beta2 <- Fb1*beta2
  Fb2.beta2 <- Fb2*beta2


  z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
  z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))

  if(is.null(pp2) & is.null(ps2)) {

    if(!is.null(pp1)) {
      xx0 <- pp1
    } else if(is.null(pp1) & !is.null(ps1)) {
      xx0 <- Fb1.beta2*ps1/(one.plus.beta2*ps1 - Fb1)
    } else {
      pp.low <- min(round(Fb1/(one.plus.beta2-Fb1.beta2), 1) + 0.05, 1)
      if(Fb1 <  0.5) {
        xx0 <- seq(pp.low, 0.5, 0.05)
      } else {
        xx0 <- seq(pp.low, 1, 0.05)
      }
    }

    allf1 <- lapply(xx0, function(x) {
      pp1 <- x
      ps1 <- Fb1.beta2*x/(one.plus.beta2*x - Fb1)

      res <- try({
        F1dist0 <- F1.cond.b(N, s, pp1, ps1, beta)
        x0 <- 1:nrow(F1dist0)
        xf10 <- F1dist0$f1s
        yf10 <- F1dist0$pf1

        idx.d <- (yf10 > 1E-6)
        x0 <- x0[idx.d]

        list(pp1=pp1, ps1=ps1, x0=x0, xf10=xf10, yf10=yf10)

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(NULL)
      } else {
        return(res)
      }

    })
    allf1 <- Filter(Negate(is.null), allf1)


    pp.low <- min(round(Fb2/(one.plus.beta2-Fb2.beta2), 1) + 0.05)
    if(Fb2 <  0.5) {
      xx <- seq(pp.low, 0.5, 0.05)
    } else {
      xx <- seq(pp.low, 1, 0.05)
    }

    yy <- sapply(xx, function(x) {
      pp2 <- x
      ps2 <- Fb2.beta2*x/(one.plus.beta2*x - Fb2)

      res <- try({

        pw.s <- sapply(allf1, function(u) {

          pp1 <- u$pp1
          ps1 <- u$ps1
          x0 <- u$x0
          xf10 <- u$xf10
          yf10 <- u$yf10

          F1dist1 <- F1.cond.b(N, s, pp2, ps2, beta)
          x1 <- 1:nrow(F1dist1)
          xf11 <- F1dist1$f1s
          yf11 <- F1dist1$pf1

          idx.b <- (yf11 > 1E-6)
          x1 <- x1[idx.b]

          ## Null distribution
          pp.bar <- (pp1+pp2)/2
          ps.bar <- (ps1+ps2)/2

          F1dist12 <- F1.cond.b(N, s, pp.bar, ps.bar, beta)
          x12 <- 1:nrow(F1dist12)
          xf112 <- F1dist12$f1s
          yf112 <- F1dist12$pf1

          idx.12 <- (yf112 > 1E-6)
          x12 <- x12[idx.12]

          if( length(xf112) > 500 & normal.approx ) {

            mean12 <- sum(xf112 * yf112)
            var12 <- sum(xf112^2 * yf112) - mean12^2

            mean0 <- sum(xf10 * yf10)
            var0 <- sum(xf10^2 * yf10) - mean0^2

            mean1 <- sum(xf11 * yf11)
            var1 <- sum(xf11^2 * yf11) - mean1^2

            dmean.bar.00 <- 0
            dvar.bar.00 <- (var12 + var12)

            dmean.bar.01 <- mean0 - mean1
            dvar.bar.01 <- (var0 + var1)

            if(alternative=="greater") {
              cr <-   z.alpha * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
            } else if(alternative=="less") {
              cr <- - z.alpha * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
            } else {
              cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
              cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
              pw.v <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
                pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

            }

          } else {

            ## Null distribution
            val1.null <- outer(x12, x12, function(d, b) xf112[d] - xf112[b])
            val2.null <- outer(x12, x12, function(d, b) yf112[d] * yf112[b])

            pf1.null <- tapply(as.vector(val2.null), as.vector(val1.null), sum)
            pf1.null <- pf1.null[pf1.null>1E-10]
            f1s.null <- as.numeric(names(pf1.null))
            pf1.null <- as.numeric(pf1.null)

            cumsum.pf1 <- cumsum(pf1.null)
            cumsum.rev.pf1 <- cumsum(rev(pf1.null))

            val1.alt <- outer(x0, x1, function(d, b) xf10[d] - xf11[b])
            val2.alt <- outer(x0, x1, function(d, b) yf10[d] * yf11[b])

            ## Alternative distribution
            pf1.alt <- tapply(as.vector(val2.alt), as.vector(val1.alt), sum)
            pf1.alt <- pf1.alt[pf1.alt>1E-10]
            pf1.alt <- pf1.alt/sum(pf1.alt)
            f1s.alt <- as.numeric(names(pf1.alt))
            pf1.alt <- as.numeric(pf1.alt)

            id.left <- findInterval(f1s.alt, f1s.null)
            id.right <- findInterval(-f1s.alt, -rev(f1s.null))

            if(alternative=="greater") {
              pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
            } else if(alternative=="less") {
              pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
            } else {
              prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
              prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

              C.obs.vec <- pmin(prob.left, prob.right)

              id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
              id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

              pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
              pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

              pvalue <- pv.left + pv.right
            }

            reject1 <- pvalue<alpha.level
            pw.v <- sum(pf1.alt[reject1])
          }

          c(pw.v, u$pp1, u$ps1)

        })

        min.pw.idx <- order(pw.s[1,])
        c(pw.s[,min.pw.idx[1]], Fb1, pp2, ps2)

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(rep(NA, 6))
      } else {
        return(res)
      }

    })

    min.pw.idx <- order(yy[1,])
    out <- yy[,min.pw.idx[1]]

  } else {
    if(is.null(pp2) & !is.null(ps2)) {
      pp2 <- Fb2.beta2*ps2/(one.plus.beta2*ps2 - Fb2)
    } else if(!is.null(pp2) & is.null(ps2)) {
      ps2 <- Fb2.beta2*pp2/(one.plus.beta2*pp2 - Fb2)
    }
    pw <- Power.Two.Classifier.Cond.Test1.beta(N=N, s=s, pp2=pp2, ps2=ps2,
                                               Fb1=Fb1, pp1=pp1, ps1=ps1,
                                               beta=beta,
                                               alternative=alternative,
                                               normal.approx=normal.approx,
                                               alpha.level=alpha.level)
    out <- pw
  }

  out <- c(out, Fb2)
  names(out) <- c("Power", "pp1", "ps1", "Fb1", "pp2", "ps2", "Fb2")
  return(round(out, 5))

}





### Not export
Power.Two.Classifier.Cond.Test1.beta <- function(N=50, s=20,
                                                 pp2=0.9, ps2=0.9,
                                                 Fb1=0.783, pp1=NULL, ps1=NULL,
                                                 beta=1,
                                                 alternative = c("two.sided", "less", "greater")[1],
                                                 normal.approx=TRUE,
                                                 alpha.level=0.05) {


  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb1)) Fb1 <- one.plus.beta2*pp1*ps1/(beta2*pp1+ps1)

  Fb1.beta2 <- Fb1*beta2

  z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
  z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))


  F1dist1 <- F1.cond.b(N, s, pp2, ps2, beta)
  x1 <- 1:nrow(F1dist1)
  xf11 <- F1dist1$f1s
  yf11 <- F1dist1$pf1

  idx.b <- (yf11 > 1E-6)
  x1 <- x1[idx.b]

  if(is.null(pp1) & is.null(ps1)) {

    pp.low <- min(round(Fb1/(one.plus.beta2-Fb1.beta2), 1) + 0.05, 1)
    if(Fb1 <  0.5) {
      xx <- seq(pp.low, 0.5, 0.05)
    } else {
      xx <- seq(pp.low, 1, 0.05)
    }

    yy <- sapply(xx, function(x) {

      res <- try({
        pp1 <- x
        ps1 <- Fb1.beta2*x/(one.plus.beta2*x - Fb1)

        F1dist0 <- F1.cond.b(N, s, pp1, ps1, beta)
        x0 <- 1:nrow(F1dist0)
        xf10 <- F1dist0$f1s
        yf10 <- F1dist0$pf1

        idx.d <- (yf10 > 1E-6)
        x0 <- x0[idx.d]

        ## Null distribution
        pp.bar <- (pp1+pp2)/2
        ps.bar <- (ps1+ps2)/2

        F1dist12 <- F1.cond.b(N, s, pp.bar, ps.bar, beta)
        x12 <- 1:nrow(F1dist12)
        xf112 <- F1dist12$f1s
        yf112 <- F1dist12$pf1

        idx.12 <- (yf112 > 1E-6)
        x12 <- x12[idx.12]

        if( length(xf112) > 500 & normal.approx ) {

          mean12 <- sum(xf112 * yf112)
          var12 <- sum(xf112^2 * yf112) - mean12^2

          mean0 <- sum(xf10 * yf10)
          var0 <- sum(xf10^2 * yf10) - mean0^2

          mean1 <- sum(xf11 * yf11)
          var1 <- sum(xf11^2 * yf11) - mean1^2

          dmean.bar.00 <- 0
          dvar.bar.00 <- (var12 + var12)

          dmean.bar.01 <- mean0 - mean1
          dvar.bar.01 <- (var0 + var1)

          if(alternative=="greater") {
            cr <-   z.alpha * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
          } else if(alternative=="less") {
            cr <- - z.alpha * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
          } else {
            cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
            cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
            pw.v <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
              pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

          }

        } else {

          val1.null <- outer(x12, x12, function(d, b) xf112[d] - xf112[b])
          val2.null <- outer(x12, x12, function(d, b) yf112[d] * yf112[b])

          val1.alt <- outer(x0, x1, function(d, b) xf10[d] - xf11[b])
          val2.alt <- outer(x0, x1, function(d, b) yf10[d] * yf11[b])

          ## Null distribution
          pf1.null <- tapply(as.vector(val2.null), as.vector(val1.null), sum)
          pf1.null <- pf1.null[pf1.null>1E-10]
          pf1.null <- pf1.null/sum(pf1.null)
          f1s.null <- as.numeric(names(pf1.null))
          pf1.null <- as.numeric(pf1.null)

          ## Alternative distribution
          pf1.alt <- tapply(as.vector(val2.alt), as.vector(val1.alt), sum)
          pf1.alt <- pf1.alt[pf1.alt>1E-10]
          pf1.alt <- pf1.alt/sum(pf1.alt)
          f1s.alt <- as.numeric(names(pf1.alt))
          pf1.alt <- as.numeric(pf1.alt)

          cumsum.pf1 <- cumsum(pf1.null)
          cumsum.rev.pf1 <- cumsum(rev(pf1.null))

          id.left <- findInterval(f1s.alt, f1s.null)
          id.right <- findInterval(-f1s.alt, -rev(f1s.null))

          if(alternative=="greater") {
            pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
          } else if(alternative=="less") {
            pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
          } else {
            prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
            prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

            C.obs.vec <- pmin(prob.left, prob.right)

            id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
            id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

            pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
            pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

            pvalue <- pv.left + pv.right
          }

          reject1 <- pvalue<alpha.level
          pw.v <- sum(pf1.alt[reject1])
        }

        c(pw.v, pp1, ps1)

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(rep(NA, 3))
      } else {
        return(res)
      }

    })
    min.pw.idx <- order(yy[1,])
    out <- yy[,min.pw.idx[1]]

  } else {

    if(is.null(pp1) & !is.null(ps1)) {
      pp1 <- Fb1.beta2*ps1/(one.plus.beta2*ps1 - Fb1)
    } else if(!is.null(pp1) & is.null(ps1)) {
      ps1 <- Fb1.beta2*pp1/(one.plus.beta2*pp1 - Fb1)
    }
    pw <- Power.Two.Classifier.Cond.Test.beta(N=N, s=s, pp2=pp2, ps2=ps2, pp1=pp1, ps1=ps1,
                                              beta=beta,
                                              alternative=alternative,
                                              normal.approx=normal.approx,
                                              alpha.level=alpha.level)
    out <- pw[1:3]
  }

  out <- c(out, Fb1, pp2, ps2)
  names(out) <- c("Power", "pp1", "ps1", "Fb1", "pp2", "ps2")
  return(out)

}



### Not export
Power.Two.Classifier.Cond.Test.beta <- function(N=50, s=20,
                                                pp2=0.9, ps2=0.9,
                                                pp1=9/11, ps1=0.75,
                                                beta=1,
                                                alternative = c("two.sided", "less", "greater")[1],
                                                normal.approx=TRUE,
                                                alpha.level=0.05) {

  ## Null distribution
  pp.bar <- (pp1+pp2)/2
  ps.bar <- (ps1+ps2)/2

  F1dist12 <- F1.cond.b(N, s, pp.bar, ps.bar, beta)
  x12 <- 1:nrow(F1dist12)
  xf112 <- F1dist12$f1s
  yf112 <- F1dist12$pf1

  idx.12 <- (yf112 > 1E-6)
  x12 <- x12[idx.12]


  ## Alternative distribution
  F1dist0 <- F1.cond.b(N, s, pp1, ps1, beta)
  F1dist1 <- F1.cond.b(N, s, pp2, ps2, beta)

  x0 <- 1:nrow(F1dist0)
  x1 <- 1:nrow(F1dist1)
  xf10 <- F1dist0$f1s
  yf10 <- F1dist0$pf1

  xf11 <- F1dist1$f1s
  yf11 <- F1dist1$pf1

  idx.d <- (yf10 > 1E-6)
  idx.b <- (yf11 > 1E-6)
  x0 <- x0[idx.d]
  x1 <- x1[idx.b]

  if( length(xf112) > 500 & normal.approx ) {

    mean12 <- sum(xf112 * yf112)
    var12 <- sum(xf112^2 * yf112) - mean12^2

    mean0 <- sum(xf10 * yf10)
    var0 <- sum(xf10^2 * yf10) - mean0^2

    mean1 <- sum(xf11 * yf11)
    var1 <- sum(xf11^2 * yf11) - mean1^2

    dmean.bar.00 <- 0
    dvar.bar.00 <- (var12 + var12)

    dmean.bar.01 <- mean0 - mean1
    dvar.bar.01 <- (var0 + var1)

    z.alpha      <- (0.4*qnorm(1-alpha.level, 0, 1) + 0.6*qt(1-alpha.level,   df=N))
    z.alpha.half <- (0.4*qnorm(1-alpha.level/2, 0, 1) + 0.6*qt(1-alpha.level/2, df=N))

    if(alternative=="greater") {
      cr <-   z.alpha * sqrt(dvar.bar.00)
      pw <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)
    } else if(alternative=="less") {
      cr <- - z.alpha * sqrt(dvar.bar.00)
      pw <- pnorm(cr, mean=dmean.bar.01, sd=sqrt(dvar.bar.01))
    } else {
      cr.L <- - z.alpha.half * sqrt(dvar.bar.00)
      cr.R <-   z.alpha.half * sqrt(dvar.bar.00)
      pw <- pnorm(cr.L, mean=dmean.bar.01, sd=sqrt(dvar.bar.01)) +
        pnorm(cr.R, mean=dmean.bar.01, sd=sqrt(dvar.bar.01), lower.tail=FALSE)

    }

  } else {

    val1.null <- outer(x12, x12, function(d, b) xf112[d] - xf112[b])
    val2.null <- outer(x12, x12, function(d, b) yf112[d] * yf112[b])

    val1.alt <- outer(x0, x1, function(d, b) xf10[d] - xf11[b])
    val2.alt <- outer(x0, x1, function(d, b) yf10[d] * yf11[b])

    ## null distribution
    pf1.null <- tapply(as.vector(val2.null), as.vector(val1.null), sum)
    pf1.null <- pf1.null[pf1.null>1E-10]
    pf1.null <- pf1.null/sum(pf1.null)
    f1s.null <- as.numeric(names(pf1.null))
    pf1.null <- as.numeric(pf1.null)

    ## alternative distribution
    pf1.alt <- tapply(as.vector(val2.alt), as.vector(val1.alt), sum)
    pf1.alt <- pf1.alt[pf1.alt>1E-10]
    pf1.alt <- pf1.alt/sum(pf1.alt)
    f1s.alt <- as.numeric(names(pf1.alt))
    pf1.alt <- as.numeric(pf1.alt)

    cumsum.pf1 <- cumsum(pf1.null)
    cumsum.rev.pf1 <- cumsum(rev(pf1.null))

    id.left <- findInterval(f1s.alt, f1s.null)
    id.right <- findInterval(-f1s.alt, -rev(f1s.null))

    if(alternative=="greater") {
      pvalue <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])
    } else if(alternative=="less") {
      pvalue <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])
    } else {
      prob.left  <- ifelse(id.left==0, 0, cumsum.pf1[pmax(1, id.left)])  #### NOTE: y[0] is excluded
      prob.right <- ifelse(id.right==0, 0, cumsum.rev.pf1[pmax(1, id.right)])

      C.obs.vec <- pmin(prob.left, prob.right)

      id.pv.left <- findInterval(C.obs.vec, cumsum.pf1)
      id.pv.right <- findInterval(C.obs.vec, cumsum.rev.pf1)

      pv.left  <- ifelse(id.pv.left==0, 0, cumsum.pf1[pmax(1, id.pv.left)])
      pv.right <- ifelse(id.pv.right==0, 0, cumsum.rev.pf1[pmax(1, id.pv.right)])

      pvalue <- pv.left + pv.right
    }

    reject1 <- pvalue<alpha.level
    pw <- sum(pf1.alt[reject1])

  }

  out <- c(pw, pp1, ps1, pp2, ps2)
  names(out) <- c("Power", "pp1", "ps1", "pp2", "ps2")
  return(out)


}

