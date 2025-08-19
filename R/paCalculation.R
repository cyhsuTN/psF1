
#' @title Solving for pa
#' @description A function to solve pa.
#' @importFrom stats dbinom uniroot

#' @param N Number of classification instances.
#' @param s Number of positives.
#' @param pp Precision.
#' @param ps Sensitivity.
#' @param beta beta.

#' @examples # paCalculation(N=50, s=20, pp=9/11, ps=0.75)


#' @export
paCalculation <- function(N, s, pp=9/11, ps=0.75, beta=1) {

  beta2 <- beta^2
  one.plus.beta2 <- (1+beta2)
  s.beta2 <- s*beta2

  Fb <- one.plus.beta2*pp*ps/(pp*beta2+ps)

  if(ps>1 | pp>1) {
    stop(paste0("pp > 1 or ps > 1"))
  }

  if(ps<0 | pp<0) {
    stop(paste0("pp < 0 or ps < 0"))
  }

  d <- 0:s
  d_probs <- dbinom(d, s, prob=ps)
  idx.d <- d_probs > 1E-10
  d <- d[idx.d]

  test1 <- function(pa) {
    b <- 0:(N-s)
    b_probs <- dbinom(b, N-s, prob=pa)

    idx.b <- b_probs > 1E-10
    b <- b[idx.b]

    cond.pre <- sum(outer(d, b, function(d, b) {
      ifelse(d == 0, 0, one.plus.beta2 * d / (d + b + s.beta2) * d_probs[d + 1] * b_probs[b + 1])
    }))

    cond.pre
  }

  if(pp == 1) {
    pa <- 0
    details <- NULL
  } else if(test1(0) < Fb) {
    stop(paste0("No pa satisfies E(fb|S,N-S,ps,pa)=Fb. Fb needs to be less than ", round(test1(0), 2),
                ". Please increase s."))
  } else if(test1(1) > Fb) {
    stop(paste0("No pa satisfies E(fb|S,N-S,ps,pa)=Fb. Fb needs to be larger than ", round(test1(1), 2),
                ". Please decrease s."))
  } else {
    f.pa <- function(b) test1(b) - Fb
    solve.pa <- uniroot(f.pa, interval = c(0, 1))
    pa <- solve.pa$root
    details <- solve.pa
  }

  list(pa = pa, details = details)
}


### Not export
paCalculation_v1 <- function(N, s, pp=9/11, ps=0.75) {

  if(ps>1 | pp>1) {
    stop(paste0("pp > 1 or ps > 1"))
  }

  if(ps<0 | pp<0) {
    stop(paste0("pp < 0 or ps < 0"))
  }

  d <- 0:s
  d_probs <- dbinom(d, s, prob=ps)
  idx.d <- d_probs > 1E-10
  d <- d[idx.d]

  test1 <- function(pa) {
    b <- 0:(N-s)
    b_probs <- dbinom(b, N-s, prob=pa)

    idx.b <- b_probs > 1E-10
    b <- b[idx.b]

    cond.pre <- sum(outer(d, b, function(d, b) {
      ifelse(d == 0, 0, d / (d + b) * d_probs[d + 1] * b_probs[b + 1])
    }))

    cond.pre
  }

  if(pp == 1) {
    pa <- 0
    details <- NULL
  } else if(test1(0) < pp) {
    stop(paste0("No pa satisfies E(ppv|S,N-S,ps,pa)=pp. pp needs to be less than ", round(test1(0), 2),
                ". Please increase s."))
  } else if(test1(1) > pp) {
    stop(paste0("No pa satisfies E(ppv|S,N-S,ps,pa)=pp. pp needs to be larger than ", round(test1(1), 2),
                ". Please decrease s."))
  } else {
    f.pa <- function(b) test1(b) - pp
    solve.pa <- uniroot(f.pa, interval = c(0, 1))
    pa <- solve.pa$root
    details <- solve.pa
  }

  list(pa = pa, details = details)
}



