#############################################################
#
#	quantile.simp.circular function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: September, 30, 2008
#	Version: 0.2-1
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

quantile.simp.circular <- function(x, probs, all=FALSE) {
  nx <- length(x)
  diameters <- rep(0, choose(nx, 2))
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  if (nx < 2) stop('x must have at least length 2')
  diameters <- .C("circdiam", x = as.double(x), nx = as.integer(nx),
               diameters = as.double(diameters),
               DUP = FALSE, NAOK = FALSE, PACKAGE = "localdepth")$diameters
  res <- quantile.default(diameters, probs)
  if (all) {
     res <- list(quantile=res, stats=diameters, call=match.call())
  }
  return(res)
}

