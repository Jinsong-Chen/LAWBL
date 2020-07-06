test_that("checking pcfa with LD data", {
  dat <- sim18cfa1$dat
  J <- ncol(dat) # no. of items
  K <- 3 # no. of factors
  Q<-matrix(-1,J,K); # -1 for unspecified items
  Q[1:6,1]<-Q[7:12,2]<-Q[13:18,3]<-1 # 1 for specified items

  m1 <- pcfa(dat = dat, Q = Q,burn = 1000, iter = 1000)

  expect_equal(TRUE, is.list(m1))
})
