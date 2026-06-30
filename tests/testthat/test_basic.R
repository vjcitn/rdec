context("basic checks")

test_that("m2ll computation succeeds with bare column names", {
 library(rdec)
 data(rdtest)
 f1 = rdec(rdec.y~rdec.x, id=rdec.id, S=rdec.s,
   omega.init = c(.5,.5), omega.low=c(.01,.01),
   omega.high=c(.95, .95), data=rdtest)
 expect_true(abs(f1$m2ll - 352.8368)<.001)
})

test_that("m2ll matches when id is passed as explicit vector", {
 library(rdec)
 data(rdtest)
 f1 = rdec(rdec.y~rdec.x, id=rdtest$rdec.id, S=rdtest$rdec.s,
   omega.init = c(.5,.5), omega.low=c(.01,.01),
   omega.high=c(.95, .95), data=rdtest)
 expect_true(abs(f1$m2ll - 352.8368)<.001)
})

test_that("character id is tolerated", {
 library(rdec)
 data(rdtest)
 rdtest2 <- rdtest
 rdtest2$rdec.id <- paste0("subj", rdtest$rdec.id)
 f1 = rdec(rdec.y~rdec.x, id=rdec.id, S=rdec.s,
   omega.init = c(.5,.5), omega.low=c(.01,.01),
   omega.high=c(.95, .95), data=rdtest2)
 expect_true(abs(f1$m2ll - 352.8368)<.001)
})
