library(htetree)

set.seed(1)
n <- 200
dat <- data.frame(
  Y = rnorm(n),
  W = rbinom(n, 1, 0.5),
  X1 = rnorm(n),
  X2 = rnorm(n)
)

fit <- htetree::causalTree(
  Y ~ X1 + X2,
  data = dat,
  treatment = dat$W,
  split.Rule = "med",
  split.Bucket = FALSE,   
  minsize = 20
)

print(fit)

print(dat)


fit_ct <- htetree::causalTree(
  Y ~ X1 + X2,
  data = dat,
  treatment = dat$W,
  split.Rule = "CT",
  split.Bucket = FALSE,
  minsize = 20
)

print(fit_ct)