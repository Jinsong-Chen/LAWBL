test_that("checking pcfa with LI data", {
  dat <- sim18cfa0$dat
  J <- ncol(dat) # no. of items
  K <- 3 # no. of factors
  Q<-matrix(-1,J,K); # -1 for unspecified items
  Q[1:2,1]<-Q[7:8,2]<-Q[13:14,3]<-1 # 1 for specified items

  m0 <- pcfa(dat = dat, Q = Q, LD = FALSE, burn = 1000, iter = 1000)

  expect_equal(TRUE, is.list(m0))
})
