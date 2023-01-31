


##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--  or do  help(data=index)  for the standard data sets.

## The function is currently defined as
corvif<- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  #vif part
  form <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz <- data.frame(fooy=1,dataz)
  lm_mod <- lm(form,dataz)
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

mygamplot2<-function (x, residuals = FALSE, rug = TRUE, se = TRUE, pages = 0,
                      select = NULL, scale = -1, n = 100, n2 = 40, pers = FALSE,
                      theta = 30, phi = 30, jit = FALSE, xlab = NULL, ylab = NULL,
                      main = NULL, ylim = NULL, xlim = NULL, too.far = 0.1, all.terms = FALSE,
                      shade = FALSE, shade.col = "gray80", shift = 0, trans = I,
                      ...)
{
  
  OUT<-c(0,0,0,0,0)
  
  
  
  
  w.resid <- NULL
  if (length(residuals) > 1) {
    if (length(residuals) == length(x$residuals))
      w.resid <- residuals
    else warning("residuals argument to plot.gam is wrong length: ignored")
    partial.resids <- TRUE
  }
  else partial.resids <- residuals
  m <- length(x$smooth)
  order <- attr(x$pterms, "order")
  if (all.terms)
    n.para <- sum(order == 1)
  else n.para <- 0
  if (m + n.para == 0)
    stop("No terms to plot - nothing for plot.gam() to do.")
  if (se) {
    if (is.numeric(se))
      se2.mult <- se1.mult <- se
    else {
      se1.mult <- 2
      se2.mult <- 1
    }
    if (se1.mult < 0)
      se1.mult <- 0
    if (se2.mult < 0)
      se2.mult <- 0
  }
  
  
  
  else se1.mult <- se2.mult <- 1
  if (se && x$Vp[1, 1] <= 0) {
    se <- FALSE
    warning("No variance estimates available")
  }
  n.plots <- m + n.para
  if (pages > n.plots)
    pages <- n.plots
  if (pages < 0)
    pages <- 0
  if (pages != 0) {
    ppp <- n.plots%/%pages
    if (n.plots%%pages != 0) {
      ppp <- ppp + 1
      while (ppp * (pages - 1) >= n.plots) pages <- pages -
        1
      if (n.plots%%pages)
        last.pages <- 0
      else last.ppp <- n.plots - ppp * pages
    }
    else last.ppp <- 0
    c <- trunc(sqrt(ppp))
    if (c < 1)
      c <- 1
    r <- ppp%/%c
    if (r < 1)
      r <- 1
    while (r * c < ppp) r <- r + 1
    while (r * c - ppp > c && r > 1) r <- r - 1
    while (r * c - ppp > r && c > 1) c <- c - 1
    oldpar <- par(mfrow = c(r, c))
  }
  else {
    ppp <- 1
    oldpar <- par()
  }
  if (partial.resids) {
    fv.terms <- predict(x, type = "terms")
    if (is.null(w.resid))
      w.resid <- x$residuals * sqrt(x$weights)
  }
  pd <- list()
  
  i <- 1
  if (m > 0)
    for (i in 1:m) {
      if (x$smooth[[i]]$dim == 1) {
        raw <- x$model[x$smooth[[i]]$term]
        xx <- seq(min(raw), max(raw), length = n)
        if (x$smooth[[i]]$by != "NA") {
          by <- rep(1, n)
          dat <- data.frame(x = xx, by = by)
          names(dat) <- c(x$smooth[[i]]$term, x$smooth[[i]]$by)
        }
        else {
          dat <- data.frame(x = xx)
          names(dat) <- x$smooth[[i]]$term
        }
        X <- PredictMat(x$smooth[[i]], dat)
        first <- x$smooth[[i]]$first.para
        last <- x$smooth[[i]]$last.para
        p <- x$coefficients[first:last]
        fit <- X %*% p
        if (se)
          se.fit <- sqrt(rowSums((X %*% x$Vp[first:last,
                                             first:last]) * X))
        edf <- sum(x$edf[first:last])
        xterm <- x$smooth[[i]]$term
        if (is.null(xlab))
          xlabel <- xterm
        else xlabel <- xlab
        if (is.null(ylab))
          ylabel <- paste("s(", xterm, ",", as.character(round(edf,
                                                               2)), ")", sep = "")
        else ylabel <- ylab
        pd.item <- list(fit = fit, dim = 1, x = xx, ylab = ylabel,
                        xlab = xlabel, raw = raw[[1]])
        if (partial.resids) {
          pd.item$p.resid <- fv.terms[, length(order) +
                                        i] + w.resid
        }
        if (se)
          pd.item$se = se.fit * se1.mult
        pd[[i]] <- pd.item
        rm(pd.item)
      }
      else if (x$smooth[[i]]$dim == 2) {
        xterm <- x$smooth[[i]]$term[1]
        if (is.null(xlab))
          xlabel <- xterm
        else xlabel <- xlab
        yterm <- x$smooth[[i]]$term[2]
        if (is.null(ylab))
          ylabel <- yterm
        else ylabel <- ylab
        raw <- data.frame(x = x$model[xterm][[1]], y = x$model[yterm][[1]])
        n2 <- max(10, n2)
        xm <- seq(min(raw$x), max(raw$x), length = n2)
        ym <- seq(min(raw$y), max(raw$y), length = n2)
        xx <- rep(xm, n2)
        yy <- rep(ym, rep(n2, n2))
        if (too.far > 0)
          exclude <- exclude.too.far(xx, yy, raw$x, raw$y,
                                     dist = too.far)
        else exclude <- rep(FALSE, n2 * n2)
        if (x$smooth[[i]]$by != "NA") {
          by <- rep(1, n)
          dat <- data.frame(x = xx, y = yy, by = by)
          names(dat) <- c(xterm, yterm, x$smooth[[i]]$by)
        }
        else {
          dat <- data.frame(x = xx, y = yy)
          names(dat) <- c(xterm, yterm)
        }
        X <- PredictMat(x$smooth[[i]], dat)
        first <- x$smooth[[i]]$first.para
        last <- x$smooth[[i]]$last.para
        p <- x$coefficients[first:last]
        fit <- X %*% p
        fit[exclude] <- NA
        if (se) {
          se.fit <- sqrt(rowSums((X %*% x$Vp[first:last,
                                             first:last]) * X))
          se.fit[exclude] <- NA
        }
        edf <- sum(x$edf[first:last])
        if (is.null(main)) {
          if (is.null(x$smooth[[i]]$margin))
            title <- paste("s(", xterm, ",", yterm, ",",
                           as.character(round(edf, 2)), ")", sep = "")
          else title <- paste("te(", xterm, ",", yterm,
                              ",", as.character(round(edf, 2)), ")", sep = "")
        }
        else title <- main
        pd.item <- list(fit = fit, dim = 2, xm = xm,
                        ym = ym, ylab = ylabel, xlab = xlabel, title = title,
                        raw = raw)
        if (is.null(ylim))
          pd.item$ylim <- range(ym)
        else pd.item$ylim <- ylim
        if (is.null(xlim))
          pd.item$xlim <- range(xm)
        else pd.item$xlim <- xlim
        if (se)
          pd.item$se = se.fit * se2.mult
        pd[[i]] <- pd.item
        rm(pd.item)
      }
      else {
        pd[[i]] <- list(dim = x$smooth[[i]]$dim)
      }
    }
  
  
  
  
  if (se) {
    k <- 0
    if (scale == -1 && is.null(ylim))
      if (m > 0)
        for (i in 1:m) {
          if (pd[[i]]$dim == 1) {
            ul <- pd[[i]]$fit + pd[[i]]$se
            ll <- pd[[i]]$fit - pd[[i]]$se
            if (k == 0) {
              ylim <- c(min(ll), max(ul))
              k <- 1
            }
            else {
              if (min(ll) < ylim[1])
                ylim[1] <- min(ll)
              if (max(ul) > ylim[2])
                ylim[2] <- max(ul)
            }
            if (partial.resids) {
              ul <- max(pd[[i]]$p.resid, na.rm = TRUE)
              if (ul > ylim[2])
                ylim[2] <- ul
              ll <- min(pd[[i]]$p.resid, na.rm = TRUE)
              if (ll < ylim[1])
                ylim[1] <- ll
            }
          }
        }
    j <- 1
    if (m > 0)
      for (i in 1:m) {
        if (is.null(select) || i == select) {
          if (interactive() && is.null(select) && pd[[i]]$dim <
                3 && i > 1 && (i - 1)%%ppp == 0)
            readline("Press return for next page....")
          if (pd[[i]]$dim == 1) {
            ul <- pd[[i]]$fit + pd[[i]]$se
            ll <- pd[[i]]$fit - pd[[i]]$se
            if (scale == 0 && is.null(ylim)) {
              ylimit <- c(min(ll), max(ul))
              if (partial.resids) {
                max.r <- max(pd[[i]]$p.resid, na.rm = TRUE)
                if (max.r > ylimit[2])
                  ylimit[2] <- max.r
                min.r <- min(pd[[i]]$p.resid, na.rm = TRUE)
                if (min.r < ylimit[1])
                  ylimit[1] <- min.r
              }
            }
            if (!is.null(ylim))
              ylimit <- ylim
            if (shade) {
              plot(pd[[i]]$x, (pd[[i]]$fit + shift),
                   type = "n", xlab = pd[[i]]$xlab, ylim =  (ylimit +
                                                               shift), xlim = xlim, ylab = pd[[i]]$ylab,
                   main = main, ...)
              polygon(c(pd[[i]]$x, pd[[i]]$x[n:1], pd[[i]]$x[1]),
                      (c(ul, ll[n:1], ul[1]) + shift),
                      col = shade.col, border = NA)
              lines(pd[[i]]$x,  (pd[[i]]$fit + shift))
            }
            else {
              plot(pd[[i]]$x,  (pd[[i]]$fit + shift),
                   type = "l", xlab = pd[[i]]$xlab, ylim =  (ylimit +
                                                               shift), xlim = xlim, ylab = pd[[i]]$ylab,
                   main = main, ...)
              
              print(shift)
              ###--->
              OUT1<-cbind(pd[[i]]$x,pd[[i]]$fit,ul,ll,rep(i,length(ll)))
              OUT<-rbind(OUT,OUT1)
              
              if (is.null(list(...)[["lty"]])) {
                lines(pd[[i]]$x,  (ul + shift), lty = 2,
                      ...)
                lines(pd[[i]]$x,  (ll + shift), lty = 2,
                      ...)
              }
              else {
                lines(pd[[i]]$x,  (ul + shift), ...)
                lines(pd[[i]]$x,  (ll + shift), ...)
              }
            }
            if (partial.resids) {
              if (is.null(list(...)[["pch"]]))
                points(pd[[i]]$raw,  (pd[[i]]$p.resid +
                                        shift), pch = ".", ...)
              else points(pd[[i]]$raw,  (pd[[i]]$p.resid +
                                           shift), ...)
            }
            if (rug) {
              if (jit)
                rug(jitter(as.numeric(pd[[i]]$raw)),
                    ...)
              else rug(as.numeric(pd[[i]]$raw), ...)
            }
          }
          else if (pd[[i]]$dim == 2) {
            if (pers) {
              if (!is.null(main))
                pd[[i]]$title <- main
              persp(pd[[i]]$xm, pd[[i]]$ym, matrix( (pd[[i]]$fit +
                                                       shift), n2, n2), xlab = pd[[i]]$xlab,
                    ylab = pd[[i]]$ylab, zlab = pd[[i]]$title,
                    ylim = pd[[i]]$ylim, xlim = pd[[i]]$xlim,
                    theta = theta, phi = phi, ...)
            }
            else {
              if (rug) {
                if (is.null(list(...)[["pch"]]))
                  points(pd[[i]]$raw$x, pd[[i]]$raw$y,
                         pch = ".", ...)
                else points(pd[[i]]$raw$x, pd[[i]]$raw$y,
                            ...)
              }
            }
          }
          else {
            warning("no automatic plotting for smooths of more than two variables")
          }
        }
        j <- j + pd[[i]]$dim
      }
  }
  else {
    k <- 0
    if (scale == -1 && is.null(ylim))
      if (m > 0)
        for (i in 1:m) {
          if (pd[[i]]$dim == 1) {
            if (k == 0) {
              if (partial.resids)
                ylim <- range(pd[[i]]$p.resid, na.rm = TRUE)
              else ylim <- range(pd[[i]]$fit)
              k <- 1
            }
            else {
              if (partial.resids) {
                if (min(pd[[i]]$p.resid) < ylim[1])
                  ylim[1] <- min(pd[[i]]$p.resid, na.rm = TRUE)
                if (max(pd[[i]]$p.resid) > ylim[2])
                  ylim[2] <- max(pd[[i]]$p.resid, na.rm = TRUE)
              }
              else {
                if (min(pd[[i]]$fit) < ylim[1])
                  ylim[1] <- min(pd[[i]]$fit)
                if (max(pd[[i]]$fit) > ylim[2])
                  ylim[2] <- max(pd[[i]]$fit)
              }
            }
          }
        }
    j <- 1
    if (m > 0)
      for (i in 1:m) {
        if (is.null(select) || i == select) {
          if (interactive() && pd[[i]]$dim < 3 && i >
                1 && (i - 1)%%ppp == 0)
            readline("Press return for next page....")
          if (pd[[i]]$dim == 1) {
            if (scale == 0 && is.null(ylim)) {
              if (partial.resids)
                ylimit <- range(pd[[i]]$p.resid, na.rm = TRUE)
              else ylimit <- range(pd[[i]]$fit)
            }
            if (!is.null(ylim))
              ylimit <- ylim
            plot(pd[[i]]$x,  (pd[[i]]$fit + shift),
                 type = "l", , xlab = pd[[i]]$xlab, ylab = pd[[i]]$ylab,
                 ylim =  (ylimit + shift), xlim = xlim,
                 main = main, ...)
            if (rug) {
              if (jit)
                rug(jitter(as.numeric(pd[[i]]$raw)),
                    ...)
              else rug(as.numeric(pd[[i]]$raw), ...)
            }
            if (partial.resids) {
              if (is.null(list(...)[["pch"]]))
                points(pd[[i]]$raw,  (pd[[i]]$p.resid +
                                        shift), pch = ".", ...)
              else points(pd[[i]]$raw,  (pd[[i]]$p.resid +
                                           shift), ...)
            }
          }
          else if (pd[[i]]$dim == 2) {
            if (!is.null(main))
              pd[[i]]$title <- main
            if (pers) {
              persp(pd[[i]]$xm, pd[[i]]$ym, matrix( (pd[[i]]$fit +
                                                       shift), n2, n2), xlab = pd[[i]]$xlab,
                    ylab = pd[[i]]$ylab, zlab = pd[[i]]$title,
                    theta = theta, phi = phi, xlim = pd[[i]]$xlim,
                    ylim = pd[[i]]$ylim, ...)
            }
            else {
              contour(pd[[i]]$xm, pd[[i]]$ym, matrix( (pd[[i]]$fit +
                                                         shift), n2, n2), xlab = pd[[i]]$xlab,
                      ylab = pd[[i]]$ylab, main = pd[[i]]$title,
                      xlim = pd[[i]]$xlim, ylim = pd[[i]]$ylim,
                      ...)
              if (rug) {
                if (is.null(list(...)[["pch"]]))
                  points(pd[[i]]$raw$x, pd[[i]]$raw$y,
                         pch = ".", ...)
                else points(pd[[i]]$raw$x, pd[[i]]$raw$y,
                            ...)
              }
            }
          }
          else {
            warning("no automatic plotting for smooths of more than one variable")
          }
        }
        j <- j + pd[[i]]$dim
      }
  }
  if (n.para > 0) {
    class(x) <- c("gam", "glm", "lm")
    if (is.null(select)) {
      attr(x, "para.only") <- TRUE
      if (interactive() && m && i%%ppp == 0)
        #readline("Press return for next page....")
        termplot(x, se = se, rug = rug, col.se = 1, col.term = 1)
    }
    else {
      if (select > m) {
        select <- select - m
        term.labels <- attr(x$pterms, "term.labels")
        term.labels <- term.labels[order == 1]
        if (select <= length(term.labels)) {
          if (interactive() && m && i%%ppp == 0)
            #readline("Press return for next page....")
            termplot(x, terms = term.labels[select], se = se,
                     rug = rug, col.se = 1, col.term = 1)
        }
      }
    }
  }
  if (pages > 0)
    par(oldpar)
  
  n<-dim(OUT)[1]
  OUT[2:n,]
}





myvif<-function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}


panel.cor <-function(x, y, digits=1, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.9/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex * r * 2)
}
# JRS - cex I belive changes font size, here its multiplied by the correlation, making strong correlations obvious, but
# its still too small, so I added the * 2


panel.hist<-function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}


panel.lines2<-function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                        cex = 1, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)){
    tmp=lm(y[ok]~x[ok])
    abline(tmp)}
}




panel.smooth2 <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),
                           cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
          col = 1, ...)
}




