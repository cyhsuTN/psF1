
#' @title Plot Distribution of F Beta Score
#' @description A function to plot a distribution of F beta score.

#' @param F1.score Objective from F1.cond.b or F1.cond.two.b.
#' @param type Density or CDF to show.
#' @param xlab Name to show in x-lab.
#' @param ylab Name to show in y-lab.
#' @param main Main title in figure.
#' @param ... Arguments in graphics::plot function.

#' @examples # fbscore <- F1.cond.b(N=50, s=20, pp=9/11, ps=0.75, beta=1)
#' @examples # plotF1(fbscore)

#' @export
plotF1 <- function(F1.score,
                   type=c("Density", "CDF")[1],
                   xlab="F_beta score",
                   ylab=NULL,
                   main="Distribution of F_beta score",
                   ...) {

  f1s <- F1.score$f1s

  if(type=="Density") {
    pf1 <- F1.score$pf1
    type1 <- "h"
    if(is.null(ylab)) ylab <- "Density"
  } else if(type=="CDF") {
    pf1 <- cumsum(F1.score$pf1)
    type1 <- "l"
    if(is.null(ylab)) ylab <- "CDF"
  } else {
    stop("What type is used: Density or CDF?")
  }

  plot(f1s, pf1, type=type1, xlab=xlab, ylab=ylab, main=main, ...)

}
