# Unit tests for gREM functions in supplementaryRscript.R

# To run, go to an R console and do
# > library(testthat) 
# > test_file('/home/tim/Dropbox/liz-paper/lucasMoorcroftManuscript/supplementary-material/testSupplementaryRscript.R')


context("Test all functions in supplementaryRscript.R")

source('/home/tim/Dropbox/liz-paper/lucasMoorcroftManuscript/supplementary-material/supplementaryRscript.R')


# calcProfileWidth() tests

test_that('Profile width is within 0 and 2*r', {
        alpha <- seq(0, 2*pi, length.out=100)
        theta <- seq(0, 2*pi, length.out=100)
        paras <- expand.grid(alpha, theta)
        paras <- cbind(paras, 1)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(sum(p>2*pi) == 0, is_true())
        expect_that(sum(p<0) == 0, is_true())
})


test_that('Profile width is 0 when alpha is 0', {
        theta <- seq(0, 2*pi, length.out=100)
        paras <- cbind(theta, 0, 1)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(all(p==0), is_true())
})


test_that('Profile width is 0 when r is 0', {
        alpha <- seq(0, 2*pi, length.out=100)
        theta <- seq(0, 2*pi, length.out=100)
        paras <- expand.grid(alpha, theta)
        paras <- cbind(paras, 0)
        colnames(paras) <- c('theta','alpha','r')
        p <- sapply(1:nrow(paras), function(x) do.call(calcProfileWidth,as.list(paras[x,])))
        expect_that(all(p==0), is_true())
})


# calcDensity() tests

test_that('Density returns error if p is 0', {
        expect_that(calcDensity(2, 0.1, 0.1, 0, 2, 2), throws_error())
        expect_that(calcDensity(2, 0.1, 0.1, 10, 0, 2), throws_error())
        expect_that(calcDensity(2, 0.1, 0.1, 10, 2, 0), throws_error())
        expect_that(calcDensity(2, 0, 0.1, 10, 2, 2), throws_error())
})

test_that('z of zero gives error', {
        expect_that(calcDensity(0, 0.1, 0.1, 2, 2, 2), throws_error())
})


# calcAbundance() tests
