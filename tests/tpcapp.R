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

X <- rbind(
  c(3, 0),
  c(-1, 0),
  c(-1, 0),
  c(-1, 0),
  c(0, 1),
  c(0, -1)
)
pc <- princomp(X)
cc_k2 <- covPC(pc, k=2)      # this works
cc_k1 <- covPC(pc, k=1)      # this should work also
 

##  Test "opt.BIC errors with k.max = 2 on a 6 * 2 matrix"
axis_data <- rbind(
  c(3, 0),
  c(-1, 0),
  c(-1, 0),
  c(-1, 0),
  c(0, 1),
  c(0, -1)
)

## This works
oo_k1 <- opt.BIC(axis_data, k.max=1, n.lambda=5, method="sd", maxiter=5, center=colMeans, scale=NULL)


## This should work also
oo_k2 <- opt.BIC(axis_data, k.max=2, n.lambda=5, method="sd", maxiter=5, center=colMeans, scale=NULL)


