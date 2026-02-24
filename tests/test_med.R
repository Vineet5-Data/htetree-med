library(htetree)

set.seed(1)
n <- 200
dat <- data.frame(
  Y = rnorm(n),
  W = rbinom(n, 1, 0.5),
  X1 = rnorm(n),
  X2 = rnorm(n)
)

fit1 <- htetree::causalTree(
  Y ~ X1 + X2,
  data = dat,
  treatment = dat$W,
  split.Rule = "medD",
  split.Bucket = FALSE,   
  minsize = 20
)

print(fit1)

print(dat)


fit_ct <- htetree::causalTree(
  Y ~ X1 + X2,
  data = dat,
  treatment = dat$W,
  split.Rule = "fit",
  split.Bucket = FALSE,
  minsize = 20
)

print(fit_ct)