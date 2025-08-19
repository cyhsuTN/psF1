
#' @title Test for Difference in F.beta Between Two Classifiers
#' @description A function to test the difference in F.beta between two classifiers.
#' @importFrom graphics abline
#' @importFrom utils tail
#' @importFrom stats dbinom uniroot

#' @param d1 Number of true positives for 1st classifier.
#' @param b1 Number of false positives for 1st classifier.
#' @param d2 Number of true positives for 2nd classifier.
#' @param b2 Number of false positives for 2nd classifier.
#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param Fb0 Null distribution is generated using the sample F.beta if Fb0 = "est";
#' Null distribution is generated using all possible F.beta if Fb0 = NULL.
#' @param beta beta.
#' @param alternative A character string specifying the alternative hypothesis:
#' c("two.sided", "less", "greater").
#' @param normal.approx normal.approx = TRUE or FALSE.

#' @examples # testFbTwo(d1=13, b1=5,
#' @examples #          d2=12, b2=3,
#' @examples #          Fb0="est",
#' @examples #          N=50, s=20,
#' @examples #          beta=1,
#' @examples #          alternative=c("two.sided", "less", "greater")[1],
#' @examples #          normal.approx=TRUE)

#' @export
testFbTwo <- function(d1=9, b1=2,
                      d2=12, b2=3,
                      Fb0="est",
                      N=50, s=20,
                      beta=1,
                      alternative=c("two.sided", "less", "greater")[1],
                      normal.approx=TRUE) {


  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb1 <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2 <- one.plus.beta2*d2/(d2 + b2 + s*beta2)

  if(!is.null(Fb0) & !is.numeric(Fb0)) {

    if(is.null(d1) | is.null(b1) | is.null(d2) | is.null(b2)) {
      stop("pp.bar and ps.bar cannot be estimated")
    }

    if (b1==0) b1 <- 0.5
    if (b2==0) b2 <- 0.5

    #pp1 <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
    #ps1 <- (d1+0.5)/(s+1)
    #pp2 <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
    #ps2 <- (d2+0.5)/(s+1)

    pp1 <- d1/(d1+b1); pp2 <- d2/(d2+b2)
    ps1 <- d1/s; ps2 <- d2/s
    pp.ave <- (pp1+pp2)/2; ps.ave <- (ps1+ps2)/2
    #pp.ave <- (d1+d2)/(d1+d2+b1+b2); ps.ave <- (d1+d2)/(s+s)

    out <- Two.Classifier.Cond.Test.beta(d1, b1,
                                         d2, b2,
                                         pp=pp.ave, ps=ps.ave,
                                         N=N, s=s,
                                         beta=beta,
                                         alternative=alternative,
                                         normal.approx=normal.approx)

    names(out) <- c("fb1", "fb2", "dfb", "P-value", "Fb0.est", "pp.bar", "ps.bar")

  } else {

    if(is.null(Fb0)) {

      if(is.null(d1) | is.null(b1) | is.null(d2) | is.null(b2)) {
        print("Search for (0, 1) instead of 95% intervals of hat{fb} due to d1, b1, d2, and b2 missing")

        fb.hat <- round(0.5 * (fb1 + fb2), 2)
        Fb0all <- seq( min(0.3, round(fb.hat,1)), 0.7, 0.1)

      } else {

        if (b1==0) b1 <- 0.5
        if (b2==0) b2 <- 0.5

        #pp1 <- (d1+0.5)/(d1+b1+1) # Jeffrey’s non-informative prior
        #ps1 <- (d1+0.5)/(s+1)
        #pp2 <- (d2+0.5)/(d2+b2+1) # Jeffrey’s non-informative prior
        #ps2 <- (d2+0.5)/(s+1)

        pp1 <- d1/(d1+b1); pp2 <- d2/(d2+b2)
        ps1 <- d1/s; ps2 <- d2/s
        pp.ave <- (pp1+pp2)/2; ps.ave <- (ps1+ps2)/2
        #pp.ave <- (d1+d2)/(d1+d2+b1+b2); ps.ave <- (d1+d2)/(s+s)

        F1dist0 <- F1.cond.b(N=N, s=s, pp=pp.ave, ps=ps.ave, beta)

        cumsum.pf1 <- cumsum(F1dist0$pf1)
        cumsum.rev.pf1 <- cumsum(rev(F1dist0$pf1))

        id.left <- findInterval(0.025, cumsum.pf1)
        id.right <- findInterval(0.025, cumsum.rev.pf1)

        Fb0all <- seq( max(round(F1dist0$f1[id.left]-0.05, 1), 0),
                       min(round(rev(F1dist0$f1)[id.right]+0.05, 1), 1), 0.1)

      }

    } else {

      Fb0all <- Fb0

    }

    outall <- matrix(NA, nrow=length(Fb0all), ncol=7)
    for(i in 1:length(Fb0all)) {
      Fb0 <- Fb0all[i]

      Fb0.beta2 <- Fb0*beta2

      pp.low <- min(round(Fb0/(one.plus.beta2-Fb0.beta2), 1) + 0.05, 1)
      if(Fb0 <  0.5) {
        xx <- seq(pp.low, 0.5, 0.05)
      } else {
        xx <- seq(pp.low, 1, 0.05)
      }

      yy <- sapply(xx, function(x) {
        res <- try({
          ps <- Fb0.beta2*x/(one.plus.beta2*x - Fb0)

          Two.Classifier.Cond.Test.beta(d1, b1,
                                        d2, b2,
                                        pp=x, ps=ps,
                                        N=N, s=s,
                                        beta=beta,
                                        alternative=alternative,
                                        normal.approx=normal.approx)
        }, silent = T)

        if (inherits(res, "try-error")) {
          return(rep(NA, 7))
        } else {
          return(res)
        }
      })
      sup.idx <- order(yy[4,], decreasing = T)[1]
      out1 <- yy[,sup.idx]
      outall[i,] <- out1
    }

    out <- outall[which.max(outall[,4]),]
    names(out) <- c("fb1", "fb2", "dfb", "P-value", "Fb0.argsup", "pp.argsup", "ps.argsup")


  }

  round(out, 5)

}


### Not export
Two.Classifier.Cond.Test.beta <- function(d1=9, b1=2,
                                          d2=9, b2=2,
                                          pp=9/11, ps=0.75,
                                          N=50, s=20,
                                          beta=1,
                                          alternative=c("two.sided", "less", "greater")[1],
                                          normal.approx=TRUE) {


  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)

  fb1 <- one.plus.beta2*d1/(d1 + b1 + s*beta2)
  fb2 <- one.plus.beta2*d2/(d2 + b2 + s*beta2)

  Fb0 <- one.plus.beta2*pp*ps/(beta2*pp + ps)

  fb <- fb1 - fb2

  F1dist0 <- F1.cond.b(N=N, s=s, pp=pp, ps=ps, beta)

  x0 <- x1 <- 1:nrow(F1dist0)
  xf10 <- xf11 <- F1dist0$f1s
  yf10 <- yf11 <- F1dist0$pf1

  idx.d <- idx.b <- (yf10 > 1E-6)
  x0 <- x0[idx.d]
  x1 <- x1[idx.b]

  if( (length(xf10) > 1E+3 | length(xf11) > 1E+3) & normal.approx ) {

    mean0 <- sum(xf10 * yf10)
    var0 <- sum(xf10^2 * yf10) - mean0^2

    mean1 <- sum(xf11 * yf11)
    var1 <- sum(xf11^2 * yf11) - mean1^2

    dmean <- mean0 - mean1
    dvar <- (var0 + var1)

    if(alternative=="greater") {
      pvalue <- pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=FALSE)
    } else if(alternative=="less") {
      pvalue <- pnorm(fb, mean=dmean, sd=sqrt(dvar))
    } else {
      ppp <- pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=TRUE)
      pvalue <- 2 * ifelse(ppp<0.5, ppp,
                           pnorm(fb, mean=dmean, sd=sqrt(dvar), lower.tail=FALSE))
    }

  } else {

    val1 <- outer(x0, x1, function(d, b) xf10[d] - xf11[b])
    val2 <- outer(x0, x1, function(d, b) yf10[d] * yf11[b])

    pf1 <- tapply(as.vector(val2), as.vector(val1), sum)
    pf1 <- pf1[pf1>1E-10]
    pf1 <- pf1/sum(pf1)
    f1s <- round(as.numeric(names(pf1)), 10); fb <- round(fb, 10)
    pf1 <- as.numeric(pf1)

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
      pvalue <- min(pvalue, 1)
    }

  }

  c(fb1=fb1, fb2=fb2, dfb=fb, 'P-value'=pvalue, Fb0=Fb0, pp=pp, ps=ps)

}






