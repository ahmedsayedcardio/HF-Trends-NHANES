multinom_wtd_mod <- function (K = 1) 
{
  if (K < 1) 
    stop("number of categories must be at least 2")
  stats <- list()
  for (i in 1:K) {
    stats[[i]] <- make.link("identity")
    fam <- structure(list(link = "identity", canonical = "none", 
                          linkfun = stats[[i]]$linkfun, mu.eta = stats[[i]]$mu.eta), 
                     class = "family")
    fam <- fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }
  residuals <- function(object, type = c("deviance")) {
    type <- match.arg(type)
    p <- object$family$predict(object$family, eta = object$fitted.values)[[1]]
    pc <- apply(p, 1, function(x) which(max(x) == x)[1]) - 
      1
    n <- length(pc)
    sgn <- rep(-1, n)
    sgn[pc == object$y] <- 1
    sgn * sqrt(-2 * log(pmax(.Machine$double.eps, p[1:n + 
                                                      object$y * n])))
  }
  predict <- function(family, se = FALSE, eta = NULL, y = NULL, 
                      X = NULL, beta = NULL, off = NULL, Vb = NULL) {
    if (is.null(eta)) {
      discrete <- is.list(X)
      lpi <- attr(X, "lpi")
      if (is.null(lpi)) {
        lpi <- list(1:ncol(X))
      }
      K <- length(lpi)
      nobs <- if (discrete) 
        nrow(X$kd)
      else nrow(X)
      eta <- matrix(0, nobs, K)
      if (se) {
        ve <- matrix(0, nobs, K)
        ce <- matrix(0, nobs, K * (K - 1)/2)
      }
      ii <- 0
      for (i in 1:K) {
        if (discrete) {
          eta[, i] <- Xbd(X$Xd, beta, k = X$kd, ks = X$ks, 
                          ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, 
                          drop = X$drop, lt = X$lpid[[i]])
        }
        else {
          Xi <- X[, lpi[[i]], drop = FALSE]
          eta[, i] <- Xi %*% beta[lpi[[i]]]
        }
        if (!is.null(off[[i]])) 
          eta[, i] <- eta[, i] + off[[i]]
        if (se) {
          ve[, i] <- if (discrete) 
            diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, 
                     ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, 
                     drop = X$drop, nthreads = 1, lt = X$lpid[[i]], 
                     rt = X$lpid[[i]])
          else drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], 
                                               lpi[[i]]]) * Xi)))
          if (i < K) 
            for (j in (i + 1):K) {
              ii <- ii + 1
              ce[, ii] <- if (discrete) 
                diagXVXd(X$Xd, Vb, k = X$kd, ks = X$ks, 
                         ts = X$ts, dt = X$dt, v = X$v, qc = X$qc, 
                         drop = X$drop, nthreads = 1, lt = X$lpid[[i]], 
                         rt = X$lpid[[j]])
              else drop(pmax(0, rowSums((Xi %*% Vb[lpi[[i]], 
                                                   lpi[[j]]]) * X[, lpi[[j]]])))
            }
        }
      }
    }
    else {
      se <- FALSE
    }
    gamma <- cbind(1, exp(eta))
    beta <- rowSums(gamma)
    gamma <- gamma/beta
    vp <- gamma * 0
    if (se) {
      for (j in 1:(K + 1)) {
        if (j == 1) 
          dp <- -gamma[, -1, drop = FALSE]/beta
        else {
          dp <- -gamma[, j] * gamma[, -1, drop = FALSE]
          dp[, j - 1] <- gamma[, j] * (1 - gamma[, j])
        }
        vp[, j] <- rowSums(dp^2 * ve)
        ii <- 0
        for (i in 1:K) if (i < K) 
          for (k in (i + 1):K) {
            ii <- ii + 1
            vp[, j] <- vp[, j] + 2 * dp[, i] * dp[, 
                                                  k] * ce[, ii]
          }
        vp[, j] <- sqrt(pmax(0, vp[, j]))
      }
      return(list(fit = gamma, se.fit = vp))
    }
    list(fit = gamma)
  }
  postproc <- expression({
    multinom <- list()
    object$y <- round(object$y)
    multinom$nj <- tabulate(object$y + 1)
    multinom$n <- sum(multinom$nj)
    multinom$K <- length(multinom$nj) - 1
    multinom$gamma <- c(1, solve(diag(multinom$n/multinom$nj[-1], 
                                      multinom$K) - matrix(1, multinom$K, multinom$K), 
                                 rep(1, multinom$K)))
    multinom$gamma <- log(multinom$gamma/sum(multinom$gamma))
    object$null.deviance <- -2 * sum(multinom$gamma[object$y + 
                                                      1])
  })
  ncv <- function(X, y, wt, nei, beta, family, llf, H = NULL, 
                  Hi = NULL, R = NULL, offset = NULL, dH = NULL, db = NULL, 
                  deriv = FALSE, nt = 1) {
    gamlss.ncv(X, y, wt, nei, beta, family, llf, H = H, 
               Hi = Hi, R = R, offset = offset, dH = dH, db = db, 
               deriv = deriv, nt = nt)
  }
  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0, 
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL, D = NULL, 
                 eta = NULL, ncv = FALSE, sandwich = FALSE) {
    n <- length(y)
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      discrete <- is.list(X)
      K <- length(jj)
      eta <- matrix(1, n, K + 1)
      if (is.null(offset)) 
        offset <- list()
      offset[[K + 1]] <- 0
      for (i in 1:K) if (is.null(offset[[i]])) 
        offset[[i]] <- 0
      for (i in 1:K) eta[, i + 1] <- offset[[i]] + if (discrete) 
        Xbd(X$Xd, coef, k = X$kd, ks = X$ks, ts = X$ts, 
            dt = X$dt, v = X$v, qc = X$qc, drop = X$drop, 
            lt = X$lpid[[i]])
      else X[, jj[[i]], drop = FALSE] %*% coef[jj[[i]]]
    }
    else {
      l2 <- 0
      K <- ncol(eta)
      eta <- cbind(1, eta)
    }
    if (K != family$nlp) 
      stop("number of linear predictors doesn't match")
    y <- round(y)
    if (min(y) < 0 || max(y) > K) 
      stop("response not in 0 to number of predictors + 1")
    ee <- exp(eta[, -1, drop = FALSE])
    beta <- 1 + rowSums(ee)
    alpha <- log(beta)
    l0 <- wt * (eta[1:n + y * n] - alpha)
    l <- sum(l0)
    l1 <- matrix(0, n, K)
    if (deriv > 0) {
      for (i in 1:K) l1[, i] <- wt * ee[, i]/beta
      l2 <- matrix(0, n, K * (K + 1)/2)
      ii <- 0
      b2 <- beta^2
      for (i in 1:K) for (j in i:K) {
        ii <- ii + 1
        l2[, ii] <- if (i == j) 
          -l1[, i] + (wt * (ee[, i]^2/b2))
        else wt * ((ee[, i] * ee[, j])/b2)
      }
      for (i in 1:K) l1[, i] <- (wt * as.numeric(y == i)) - l1[, 
                                                        i]
    }
    l3 <- l4 <- 0
    tri <- family$tri
    if (deriv > 1) {
      l3 <- matrix(0, n, (K * (K + 3) + 2) * K/6)
      ii <- 0
      b3 <- b2 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) {
        ii <- ii + 1
        if (i == j && j == k) {
          l3[, ii] <- l2[, tri$i2[i, i]] + (wt * (2 * ee[, 
                                                  i]^2/b2 - 2 * ee[, i]^3/b3))
        }
        else if (i != j && j != k & i != k) {
          l3[, ii] <- wt * (-2 * (ee[, i] * ee[, j] * ee[, 
                                                   k])/b3)
        }
        else {
          kk <- if (i == j) 
            k
          else j
          l3[, ii] <- l2[, tri$i2[i, kk]] - (wt * (2 * (ee[, 
                                                    i] * ee[, j] * ee[, k])/b3))
        }
      }
    }
    if (deriv > 3) {
      l4 <- matrix(0, n, (6 + K * 11 + K^2 * 6 + K^3) * 
                     K/24)
      ii <- 0
      b4 <- b3 * beta
      for (i in 1:K) for (j in i:K) for (k in j:K) for (l in k:K) {
        ii <- ii + 1
        uni <- unique(c(i, j, k, l))
        nun <- length(uni)
        if (nun == 1) {
          l4[, ii] <- l3[, tri$i3[i, i, i]] + 4 * ee[, 
                                                     i]^2/b2 - 10 * ee[, i]^3/b3 + 6 * ee[, i]^4/b4
        }
        else if (nun == 4) {
          l4[, ii] <- 6 * ee[, i] * ee[, j] * ee[, k] * 
            ee[, l]/b4
        }
        else if (nun == 3) {
          l4[, ii] <- l3[, tri$i3[uni[1], uni[2], uni[3]]] + 
            6 * ee[, i] * ee[, j] * ee[, k] * ee[, l]/b4
        }
        else if (sum(uni[1] == c(i, j, k, l)) == 2) {
          l4[, ii] <- l3[, tri$i3[uni[1], uni[2], uni[2]]] - 
            2 * ee[, uni[1]]^2 * ee[, uni[2]]/b3 + 6 * 
            ee[, i] * ee[, j] * ee[, k] * ee[, l]/b4
        }
        else {
          if (sum(uni[1] == c(i, j, k, l)) == 1) 
            uni <- uni[2:1]
          l4[, ii] <- l3[, tri$i3[uni[1], uni[1], uni[2]]] - 
            4 * ee[, uni[1]]^2 * ee[, uni[2]]/b3 + 6 * 
            ee[, i] * ee[, j] * ee[, k] * ee[, l]/b4
        }
      }
    }
    if (deriv) {
      ret <- gamlss.gH(X, jj, l1, l2, tri$i2, l3 = l3, 
                       i3 = tri$i3, l4 = l4, i4 = tri$i4, d1b = d1b, 
                       d2b = d2b, deriv = deriv - 1, fh = fh, D = D, 
                       sandwich = sandwich)
      if (ncv) {
        ret$l1 = l1
        ret$l2 = l2
        ret$l3 = l3
      }
    }
    else ret <- list()
    ret$l <- l
    ret
  }
  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1, 
       sandwich = TRUE)$lbb
  }
  rd <- function(mu, wt, scale) {
    p <- exp(cbind(0, mu))
    p <- p/rowSums(p)
    cp <- t(apply(p, 1, cumsum))
    apply(cp, 1, function(x) min(which(x > runif(1)))) - 
      1
  }
  initialize <- expression({
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      jj <- attr(x, "lpi")
      if (is.list(x)) {
        start <- rep(0, max(unlist(jj)))
        for (k in 1:length(jj)) {
          yt1 <- 6 * as.numeric(y == k) - 3
          R <- suppressWarnings(mchol(XWXd(x$Xd, w = rep(1, 
                                                         length(y)), k = x$kd, ks = x$ks, ts = x$ts, 
                                           dt = x$dt, v = x$v, qc = x$qc, nthreads = 1, 
                                           drop = x$drop, lt = x$lpid[[k]]) + crossprod(E[, 
                                                                                          jj[[k]]])))
          Xty <- XWyd(x$Xd, rep(1, length(y)), yt1, 
                      x$kd, x$ks, x$ts, x$dt, x$v, x$qc, x$drop, 
                      lt = x$lpid[[k]])
          piv <- attr(R, "pivot")
          rrank <- attr(R, "rank")
          startji <- rep(0, ncol(R))
          if (rrank < ncol(R)) {
            R <- R[1:rrank, 1:rrank]
            piv <- piv[1:rrank]
          }
          startji[piv] <- backsolve(R, forwardsolve(t(R), 
                                                    Xty[piv]))
          startji[!is.finite(startji)] <- 0
          start[jj[[k]]] <- startji
        }
      } else {
        start <- rep(0, ncol(x))
        for (k in 1:length(jj)) {
          yt1 <- 6 * as.numeric(y == k) - 3
          x1 <- x[, jj[[k]], drop = FALSE]
          e1 <- E[, jj[[k]], drop = FALSE]
          if (use.unscaled) {
            qrx <- qr(rbind(x1, e1))
            x1 <- rbind(x1, e1)
            startji <- qr.coef(qr(x1), c(yt1, rep(0, 
                                                  nrow(E))))
            startji[!is.finite(startji)] <- 0
          } else startji <- pen.reg(x1, e1, yt1)
          start[jj[[k]]] <- startji
        }
      }
    }
  })
  structure(list(family = "multinom", ll = ll, link = NULL, 
                 nlp = round(K), rd = rd, ncv = ncv, sandwich = sandwich, 
                 tri = trind.generator(K), initialize = initialize, postproc = postproc, 
                 residuals = residuals, predict = predict, linfo = stats, 
                 d2link = 1, d3link = 1, d4link = 1, ls = 1, available.derivs = 2, 
                 discrete.ok = TRUE), class = c("general.family", "extended.family", 
                                                "family"))
}
