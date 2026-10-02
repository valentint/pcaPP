## this will render the output independent from the version of the package
suppressPackageStartupMessages(library(pcaPP))

##  Sparse loadings
set.seed (0)
##x <- data.Zou()
data(iris)
x <- iris[, 1:4]

l1median_NM (x)$par
l1median_CG (x)$par
l1median_BFGS (x)$par
l1median_NLM (x)$par
l1median_HoCr (x)$par
l1median_VaZh (x)$par

# compare with coordinate-wise median:
apply(x,2,median)

pc <- PCAgrid(x)
pc
summary(pc)
pc$loadings
pc$scores

##  Test "covPC() errors with k = 1 for a princomp object with multiple variables"
library(pcaPP)

X <- rbind(
  c(3, 0),
  c(-1, 0),
  c(-1, 0),
  c(-1, 0),
  c(0, 1),
  c(0, -1)
)
pc <- princomp(X)
covPC(pc, k=2)      # this works
covPC(pc, k=1)      # this should work also
 
Y <- rbind(
  c(3, 0, 1),
  c(-1, 0, 1),
  c(-1, 0, 0),
  c(-1, 0, 1),
  c(0, 1, 0),
  c(0, -1, 1)
)

pc <- princomp(Y)

cc <- covPC(pc, k=3) 
cc0 <- matrix(c(2.000000e+00, -1.301043e-17,  0.1666667,
                -1.604619e-17,  3.333333e-01, -0.1666667,
                1.666667e-01, -1.666667e-01,  0.2222222), nrow=3)
all.equal(cc$cov, cc0, tolerance=1e-6)

cc <- covPC(pc, k=2) 
cc0 <- matrix(c(1.999525491,  0.003756655,  0.1720978,
                0.003756655,  0.303592164, -0.2096643,
                0.172097769, -0.209664319,  0.1600593), nrow=3)
all.equal(cc$cov, cc0, tolerance=1e-6)

covPC(pc, k=1) 
