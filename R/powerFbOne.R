
#' @title Power Calculation for F.beta of One Classifier
#' @description A function to calculate power for testing the F.beta of one classifier.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot qnorm pnorm

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param Fb1 F.beta score of the new classifier.
#' @param pp1 Precision of the new classifier.
#' @param ps1 Sensitivity of the new classifier.
#' @param Fb Target F.beta score.
#' @param pp Target precision.
#' @param ps Target sensitivity.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = TRUE or FALSE.
#' @param alpha.level alpha level.

#' @examples # powerFbOne(N=50, s=35,
#' @examples #            Fb1=0.9, pp1=NULL, ps1=NULL,
#' @examples #            Fb=0.783, pp=NULL, ps=NULL,
#' @examples #            beta=1,
#' @examples #            alternative=c("two.sided", "less", "greater")[1],
#' @examples #            normal.approx=TRUE,
#' @examples #            alpha.level=0.05)

#' @export
powerFbOne <- function(N=50, s=20,
                       Fb1=0.9, pp1=NULL, ps1=NULL,
                       Fb=0.783, pp=NULL, ps=NULL,
                       beta=1,
                       alternative=c("two.sided", "less", "greater")[1],
                       normal.approx=TRUE,
                       alpha.level=0.05) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb)) Fb <- one.plus.beta2*pp*ps/(beta2*pp+ps)
  if(is.null(Fb1)) Fb1 <- one.plus.beta2*pp1*ps1/(beta2*pp1+ps1)

  Fb.beta2 <- Fb*beta2
  Fb1.beta2 <- Fb1*beta2


  z.alpha      <- qnorm(1-alpha.level, 0, 1)
  z.alpha.half <- qnorm(1-alpha.level/2, 0, 1)

  if(is.null(pp1) & is.null(ps1)) {

    if(!is.null(pp)) {
      xx <- pp
    } else if(is.null(pp) & !is.null(ps)) {
      xx <- Fb.beta2*ps/(one.plus.beta2*ps - Fb)
    } else {
      xx <- seq(round(Fb/(one.plus.beta2-Fb.beta2), 2)+0.01, 1, 0.01)
    }

    allf1 <- lapply(xx, function(x) {
      pp <- x
      ps <- Fb.beta2*x/(one.plus.beta2*x - Fb)

      res <- try({
        F1.score <- F1.cond.b(N, s, pp, ps, beta)
        f1s <- F1.score$f1s
        pf1 <- F1.score$pf1

        if( length(f1s) > 1E+4 & normal.approx ) {

          list(pp=pp, ps=ps, f1s=f1s, pf1=pf1)

        } else {

          cumsum.pf1 <- cumsum(pf1)
          cumsum.rev.pf1 <- cumsum(rev(pf1))
          list(pp=pp, ps=ps,
               f1s=f1s, pf1=pf1,
               cumsum.pf1=cumsum.pf1, cumsum.rev.pf1=cumsum.rev.pf1)

        }

      }, silent = T)

      if (inherits(res, "try-error")) {
        return(NULL)
      } else {
        return(res)
      }

    })
    allf1 <- Filter(Negate(is.null), allf1)

    xx <- seq(round(Fb1/(one.plus.beta2-Fb1.beta2), 2)+0.01, 1, 0.01)
    yy <- sapply(xx, function(x) {
      pp1 <- x
      ps1 <- Fb1.beta2*x/(one.plus.beta2*x - Fb1)

      res <- try({
        F1.score1 <- F1.cond.b(N, s, pp1, ps1, beta)
        f1s1 <- F1.score1$f1s
        pf11 <- F1.score1$pf1

        pw.s <- sapply(allf1, function(u) {
          f1s <- u$f1s
          pf1 <- u$pf1

          if( (length(f1s) > 1E+4 | length(f1s1) > 1E+4) & normal.approx ) {

            mean0 <- sum(f1s * pf1)
            var0 <- sum(f1s^2 * pf1) - mean0^2

            mean1 <- sum(f1s1 * pf11)
            var1 <- sum(f1s1^2 * pf11) - mean1^2

            if(alternative=="greater") {
              cr <- mean0 + z.alpha * sqrt(var0)
              pw.v <- pnorm(cr, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)
            } else if(alternative=="less") {
              cr <- mean0 - z.alpha * sqrt(var0)
              pw.v <- pnorm(cr, mean=mean1, sd=sqrt(var1))
            } else {
              cr.L <- mean0 - z.alpha.half * sqrt(var0)
              cr.R <- mean0 + z.alpha.half * sqrt(var0)
              pw.v <- pnorm(cr.L, mean=mean1, sd=sqrt(var1)) +
                pnorm(cr.R, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)

            }

          } else {

            cumsum.pf1 <- u$cumsum.pf1
            cumsum.rev.pf1 <- u$cumsum.rev.pf1

            id.left <- findInterval(f1s1, f1s)
            id.right <- findInterval(-f1s1, -rev(f1s))

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
            pw.v <- sum(pf11[reject1])

          }

          c(pw.v, u$pp, u$ps)

        })

        min.pw.idx <- order(pw.s[1,])
        c(pw.s[,min.pw.idx[1]], Fb, pp1, ps1)

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
    if(is.null(pp1) & !is.null(ps1)) {
      pp1 <- Fb1.beta2*ps1/(one.plus.beta2*ps1 - Fb1)
    } else if(!is.null(pp1) & is.null(ps1)) {
      ps1 <- Fb1.beta2*pp1/(one.plus.beta2*pp1 - Fb1)
    }
    pw <- Power.One.Classifier.Cond.Test1.beta(N=N, s=s, pp1=pp1, ps1=ps1,
                                               Fb=Fb, pp=pp, ps=ps,
                                               beta=beta,
                                               alternative=alternative,
                                               normal.approx=normal.approx,
                                               alpha.level=alpha.level)
    out <- pw
  }

  out <- c(out, Fb1)
  names(out) <- c("Power", "pp", "ps", "Fb", "pp1", "ps1", "Fb1")
  return(round(out, 5))
}


### Not export
Power.One.Classifier.Cond.Test1.beta <- function(N=50, s=20,
                                                 pp1=0.9, ps1=0.9,
                                                 Fb=0.783, pp=NULL, ps=NULL,
                                                 beta=1,
                                                 alternative = c("two.sided", "less", "greater")[1],
                                                 normal.approx=TRUE,
                                                 alpha.level=0.05) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  if(is.null(Fb)) Fb <- one.plus.beta2*pp*ps/(beta2*pp+ps)

  Fb.beta2 <- Fb*beta2


  z.alpha      <- qnorm(1-alpha.level, 0, 1)
  z.alpha.half <- qnorm(1-alpha.level/2, 0, 1)

  if(is.null(pp) & is.null(ps)) {

    F1.score1 <- F1.cond.b(N, s, pp1, ps1, beta)
    f1s1 <- F1.score1$f1s
    pf11 <- F1.score1$pf1

    xx <- seq(round(Fb/(one.plus.beta2-Fb.beta2), 2)+0.01, 1, 0.01)
    yy <- sapply(xx, function(x) {
      pp <- x
      ps <- Fb.beta2*x/(one.plus.beta2*x - Fb)

      res <- try({
        F1.score <- F1.cond.b(N, s, pp, ps, beta)
        f1s <- F1.score$f1s
        pf1 <- F1.score$pf1

        if( (length(f1s) > 1E+5 | length(f1s1) > 1E+5) & normal.approx ) {

          mean0 <- sum(f1s * pf1)
          var0 <- sum(f1s^2 * pf1) - mean0^2

          mean1 <- sum(f1s1 * pf11)
          var1 <- sum(f1s1^2 * pf11) - mean1^2

          if(alternative=="greater") {
            cr <- mean0 + z.alpha * sqrt(var0)
            pw.v <- pnorm(cr, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)
          } else if(alternative=="less") {
            cr <- mean0 - z.alpha * sqrt(var0)
            pw.v <- pnorm(cr, mean=mean1, sd=sqrt(var1))
          } else {
            cr.L <- mean0 - z.alpha.half * sqrt(var0)
            cr.R <- mean0 + z.alpha.half * sqrt(var0)
            pw.v <- pnorm(cr.L, mean=mean1, sd=sqrt(var1)) +
              pnorm(cr.R, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)

          }

        } else {

          cumsum.pf1 <- cumsum(pf1)
          cumsum.rev.pf1 <- cumsum(rev(pf1))

          id.left <- findInterval(f1s1, f1s)
          id.right <- findInterval(-f1s1, -rev(f1s))

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
          pw.v <- sum(pf11[reject1])

        }

        c(pw.v, pp, ps)

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
    if(is.null(pp) & !is.null(ps)) {
      pp <- Fb.beta2*ps/(one.plus.beta2*ps - Fb)
    } else if(!is.null(pp) & is.null(ps)) {
      ps <- Fb.beta2*pp/(one.plus.beta2*pp - Fb)
    }
    pw <- Power.One.Classifier.Cond.Test.beta(N=N, s=s, pp1=pp1, ps1=ps1, pp=pp, ps=ps,
                                              beta=beta,
                                              alternative=alternative,
                                              normal.approx=normal.approx,
                                              alpha.level=alpha.level)
    out <- pw[1:3]
  }


  out <- c(out, Fb, pp1, ps1)
  names(out) <- c("Power", "pp", "ps", "Fb", "pp1", "ps1")
  return(out)

}



### Not export
Power.One.Classifier.Cond.Test.beta <- function(N=50, s=20,
                                                pp1=0.9, ps1=0.9,
                                                pp=9/11, ps=0.75,
                                                beta=1,
                                                alternative = c("two.sided", "less", "greater")[1],
                                                normal.approx=TRUE,
                                                alpha.level=0.05) {


  F1.score <- F1.cond.b(N, s, pp, ps, beta)
  f1s <- F1.score$f1s
  pf1 <- F1.score$pf1

  F1.score1 <- F1.cond.b(N, s, pp1, ps1, beta)
  f1s1 <- F1.score1$f1s
  pf11 <- F1.score1$pf1

  if( (length(f1s) > 1E+5 | length(f1s1) > 1E+5) & normal.approx ) {

    mean0 <- sum(f1s * pf1)
    var0 <- sum(f1s^2 * pf1) - mean0^2

    mean1 <- sum(f1s1 * pf11)
    var1 <- sum(f1s1^2 * pf11) - mean1^2

    z.alpha      <- qnorm(1-alpha.level, 0, 1)
    z.alpha.half <- qnorm(1-alpha.level/2, 0, 1)

    if(alternative=="greater") {
      cr <- mean0 + z.alpha * sqrt(var0)
      pw <- pnorm(cr, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)
    } else if(alternative=="less") {
      cr <- mean0 - z.alpha * sqrt(var0)
      pw <- pnorm(cr, mean=mean1, sd=sqrt(var1))
    } else {
      cr.L <- mean0 - z.alpha.half * sqrt(var0)
      cr.R <- mean0 + z.alpha.half * sqrt(var0)
      pw <- pnorm(cr.L, mean=mean1, sd=sqrt(var1)) +
        pnorm(cr.R, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)

    }


  } else {

    cumsum.pf1 <- cumsum(pf1)
    cumsum.rev.pf1 <- cumsum(rev(pf1))

    id.left <- findInterval(f1s1, f1s)
    id.right <- findInterval(-f1s1, -rev(f1s))

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
    pw <- sum(pf11[reject1])

  }



  out <- c(pw, pp, ps, pp1, ps1)
  names(out) <- c("Power", "pp", "ps", "pp1", "ps1")
  return(out)

}


