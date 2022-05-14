#' @title Posterior plots for \code{lawbl} object
#'
#' @description Provide posterior plots based on the factorial eigenvalues of a \code{lawbl} object.
#' For \code{PEFA} or \code{FEFA}, only true factors will be plotted.
#'
#' @name plot_lawbl
#'
#' @param object A \code{lawbl} object
#'
#' @param what A list of options for what to plot.
#'
#' \itemize{
#'      \item \code{trace}: The trace of each factor's eigenvalue.
#'     \item \code{density}: The trace of each factor's eigenvalue.
#'     \item \code{EPSR}: Estimated Potential Scale Reduction (Gelman-Rubin diagnostics) for each factor.
#'  }
#'
#' @param istart Starting point of the Markov chain for plotting.
#'
#' @param iend Ending point of the Markov chain for plotting; -1 for the actual final point.
#'
#' @return NULL
#'
#' @importFrom graphics plot
#' @importFrom graphics par
#'
#' @export
#'
#' @examples
#' \donttest{
#' dat <- sim18cfa0$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[7:8,2]<-Q[13:14,3]<-1
#'
#'
#' m0 <- pcfa(dat = dat, Q = Q, LD = FALSE,burn = 1000, iter = 1000)
#' plot_lawbl(m0) # trace
#' plot_lawbl(m0, what='density')
#' plot_lawbl(m0, what='EPSR')
#' }
plot_lawbl <- function(object, what = "trace",istart =1, iend = -1) {
    # if (class(obj) != "lawbl")
    #     stop("It must be a lawbl object.", call. = F)

    oldmar <- par("mar")
    par(mar = rep(2, 4))

    Q <- object$Q
    iter <- object$iter
    K <- ncol(Q)

    if (iend == -1 | iend > iter)
        iend <- iter
    iter <- iend - istart + 1

    # poq <- which(Q != 0, arr.ind = T)
    # eig_arr <- array(0, dim = c(iter, K))
    # for (k in 1:K) {
    #     ind1 <- (poq[, 2] == k)
    #     eig_arr[, k] <- rowSums(obj$LA[, ind1]^2)
    # }
    # colnames(eig_arr) <- paste0("F", c(1:K))
    # mobj <- mcmc(eig_arr)

    eigen <-object$Eigen[istart:iend,]
    colnames(eigen) <- paste0("F", c(1:K))
    eig_eps <- object$eig_eps
    if(is.null(eig_eps)) eig_eps <- .1
    TF_ind<-(colMeans(eigen)>eig_eps)
    eig_arr<- as.matrix(eigen[,TF_ind])
    mobj <- mcmc(eig_arr)

    x1 <- mcmc(mobj[1:(iter/2), ])
    x2 <- mcmc(mobj[(iter/2 + 1):iter, ])
    # x1<-mcmc(mobj[1:(iter/2)]);x2<-mcmc(mobj[(iter/2+1):iter])
    xx <- mcmc.list(x1, x2)
    # gelman.diag(xx)

    switch(what,
           EPSR = gelman.plot(xx, autoburnin = F,xlab = "", ylab = "PSRF"),
           trace = plot(mcmc(eig_arr), density = F,xlab = "",cex.main = 1),
           density = plot(mcmc(eig_arr), trace = F, xlab = "",cex.main = 1),
        stop(sprintf("Can not plot element '%s'", what), call. = FALSE))

    par(mar = oldmar) #reset to old mar
}
