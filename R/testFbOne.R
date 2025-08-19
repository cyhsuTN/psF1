
#' @title Test for F.beta of One Classifier
#' @description A function to test the F.beta of one classifier.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot pnorm

#' @param d Number of true positives.
#' @param b Number of false positives.
#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param Fb Target F.beta score.
#' @param pp Target precision.
#' @param ps Target sensitivity.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = TRUE or FALSE.

#' @examples # testFbOne(d=9, b=2, N=50, s=20,
#' @examples #           Fb=0.783, pp=NULL, ps=NULL,
#' @examples #           beta=1,
#' @examples #           alternative=c("two.sided", "less", "greater")[1],
#' @examples #           normal.approx=TRUE)

#' @export
testFbOne <- function(d=9, b=2, N=50, s=20,
                      Fb=0.783, pp=NULL, ps=NULL,
                      beta=1,
                      alternative=c("two.sided", "less", "greater")[1],
                      normal.approx=TRUE) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb <- one.plus.beta2*d/(d + b + s*beta2)

  Fb.beta2 <- Fb*beta2


  if(is.null(pp) & is.null(ps)) {

    xx <- seq(round(Fb/(one.plus.beta2-Fb.beta2), 2)+0.01, 1, 0.01)
    yy <- sapply(xx, function(x) {
      ps <- Fb.beta2*x/(one.plus.beta2*x - Fb)
      res <- try({
        One.Classifier.Cond.Test.beta(d=d, b=b, N=N, s=s, pp=x, ps=ps,
                                      beta=beta,
                                      alternative=alternative,
                                      normal.approx=normal.approx)
      }, silent = T)

      if (inherits(res, "try-error")) {
        return(rep(NA, 5))
      } else {
        return(res)
      }
    })
    sup.idx <- order(yy[2,], decreasing = T)[1]
    out <- yy[,sup.idx]
    names(out) <- c("fb", "P-value", "Fb", "pp.argsup", "ps.argsup")

  } else {

    if(is.null(pp) & !is.null(ps)) {
      pp <- Fb.beta2*ps/(one.plus.beta2*ps - Fb)
    } else if(!is.null(pp) & is.null(ps)) {
      ps <- Fb.beta2*pp/(one.plus.beta2*pp - Fb)
    }
    out <- One.Classifier.Cond.Test.beta(d=d, b=b, N=N, s=s, pp=pp, ps=ps,
                                         beta=beta,
                                         alternative=alternative,
                                         normal.approx=normal.approx)
  }

  round(out, 5)

}


### Not export
One.Classifier.Cond.Test.beta <- function(d=9, b=2, N=50, s=20,
                                          pp=9/11, ps=0.75,
                                          beta=1,
                                          alternative=c("two.sided", "less", "greater")[1],
                                          normal.approx=TRUE) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb <- one.plus.beta2*d/(d + b + s*beta2)

  Fb <- one.plus.beta2*pp*ps/(beta2*pp + ps)

  F1.score <- F1.cond.b(N, s, pp, ps, beta)
  f1s <- round(F1.score$f1s, 10); fb <- round(fb, 10)
  pf1 <- F1.score$pf1

  if( length(f1s) > 1E+5 & normal.approx ) {

    mean1 <- sum(f1s * pf1)
    var1 <- sum(f1s^2 * pf1) - mean1^2

    if(alternative=="greater") {
      pvalue <- pnorm(fb, mean=mean1, sd=sqrt(var1), lower.tail=FALSE)
    } else if(alternative=="less") {
      pvalue <- pnorm(fb, mean=mean1, sd=sqrt(var1))
    } else {
      ppp <- pnorm(fb, mean=mean1, sd=sqrt(var1), lower.tail=TRUE)
      pvalue <- 2 * ifelse(ppp<0.5, ppp,
                           pnorm(fb, mean=mean1, sd=sqrt(var1), lower.tail=FALSE))
    }

  } else {

    if(alternative=="greater") {
      pvalue <- sum(pf1[f1s>fb]) + 0.5*sum(pf1[f1s=fb])
    } else if(alternative=="less") {
      pvalue <- sum(pf1[f1s<fb]) + 0.5*sum(pf1[f1s=fb])
    } else {
      C.obs <- min(sum(pf1[f1s<=fb]), sum(pf1[f1s>=fb]))

      cumsum.pf1 <- cumsum(pf1)
      cumsum.rev.pf1 <- cumsum(rev(pf1))

      pvalue <- as.numeric(
        ifelse(!any(cumsum.rev.pf1<=C.obs), 0,
               tail(cumsum.rev.pf1[cumsum.rev.pf1<=C.obs],1)) +
          ifelse(!any(cumsum.pf1<=C.obs), 0,
                 tail(cumsum.pf1[cumsum.pf1<=C.obs],1)))
    }

  }

  c(fb=fb, 'P-value'=pvalue, Fb=Fb, pp=pp, ps=ps)

}

