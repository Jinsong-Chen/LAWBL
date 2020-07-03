#' @title Factorial eigenvalue plot
#'
#' @description Provide plots based on the factorial eigenvalues of a \code{pcfa} objects.
#'
#'
#' @name plot_eigen
#'
#' @param obj A \code{lawbl} object
#'
#' @param what A list of options for what to plot.
#'
#' \itemize{
#'      \item \code{trace}: The trace of each factor's eigenvalue.
#'     \item \code{density}: The trace of each factor's eigenvalue.
#'     \item \code{PGR}: The pseudo Gelman-Rubin diagnostics of each factor's eigenvalue.
#'  }
#'
#'
#' @return NULL
#'
#' @importFrom graphics plot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- sim18cfa0$dat
#' J <- ncol(dat)
#' K <- 3
#' Q<-matrix(-1,J,K);
#' Q[1:2,1]<-Q[9:10,2]<-Q[13:14,3]<-1
#'
#'
#' mod0 <- pcfa(dat = dat, Q = Q, LD = F,burn = 2000, iter = 2000)
#' plot_eigen(mod0) # trace
#' plot_eigen(mod0, what='density')
#' plot_eigen(mod0, what='PGR')
#' }
plot_eigen <- function(obj, what = "trace") {
    if (class(obj) != "lawbl")
        stop("It must be a lawbl object.", call. = F)
    Q <- obj$Q
    poq <- which(Q != 0, arr.ind = T)
    iter <- obj$iter
    K <- ncol(Q)
    eig_arr <- array(0, dim = c(iter, K))
    for (k in 1:K) {
        ind1 <- (poq[, 2] == k)
        eig_arr[, k] <- rowSums(obj$LA[, ind1]^2)
    }
    colnames(eig_arr) <- paste0("F", c(1:K))
    mobj <- mcmc(eig_arr)
    x1 <- mcmc(mobj[1:(iter/2), ])
    x2 <- mcmc(mobj[(iter/2 + 1):iter, ])
    # x1<-mcmc(mobj[1:(iter/2)]);x2<-mcmc(mobj[(iter/2+1):iter])
    xx <- mcmc.list(x1, x2)
    # gelman.diag(xx)

    switch(what,
           APSR = gelman.plot(xx, autoburnin = F,xlab = "", ylab = "PSRF"),
           trace = plot(mcmc(eig_arr), density = F,xlab = "",cex.main = 1),
           density = plot(mcmc(eig_arr), trace = F, xlab = "",cex.main = 1),
        stop(sprintf("Can not plot element '%s'", what), call. = FALSE))

}
