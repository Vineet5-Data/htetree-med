#' Robust Causal Effect Forests
#'
#' Construct a causal forest that can dispatch to either the standard
#' ``causalTree``/``honest.causalTree`` functions or a user‑implemented
#' ``robustcausalTree`` based on a new set of splitting rules (e.g.
#' ``"MAD"``, ``"LMS"`` and ``"MSD"``).  This function mirrors the
#' behaviour of the original ``causalForest`` from the htetree package
#' while adding logic to recognise the new robust rules and route
#' accordingly.  For standard split rules the existing htetree tree
#' builders are used unchanged.  The function accepts the same
#' arguments as ``causalForest``, with two important differences:
#'
#' * The ``split.Rule`` argument may include the new robust rules
#'   ``"MAD"``, ``"LMS"`` or ``"MSD"`` (case insensitive).  These
#'   rules are mapped to integer codes and passed to the C layer via
#'   the tree builder.  Bucketised variants of the robust rules are
#'   automatically handled by ``robustcausalTree``—they do not have
#'   "D" appended to the rule name.
#' * When a robust rule is selected, this function calls
#'   ``robustcausalTree`` (or ``robustcausalTree`` in honesty mode)
#'   instead of the standard ``causalTree``/``honest.causalTree``.
#'
#' The sampling logic, data preparation, and forest bookkeeping are
#' inherited from ``causalForest``.  See ``?causalForest`` for a full
#' description of the arguments and return value.
#'
#' @inheritParams causalForest
#' @param split.Rule Character string giving the splitting rule.  In
#'   addition to the standard rules recognised by ``causalTree``, the
#'   values ``"MAD"``, ``"LMS"`` and ``"MSD"`` are supported.  These
#'   robust rules are handled by calling ``robustcausalTree``.  If
#'   ``split.Bucket`` is ``TRUE``, the robust rules will still be
#'   passed unchanged to ``robustcausalTree``; bucketised behaviour is
#'   handled internally by that function.
#' @param ... Additional arguments are passed on to the underlying
#'   tree builder (``causalTree``, ``honest.causalTree`` or
#'   ``robustcausalTree`` as appropriate).  Arguments that are not
#'   recognised by ``rpart.control`` will trigger an error.
#'
#' @return An object of class ``causalForest`` with the usual
#'   components ``trees``, ``fsample``, ``cov_sample`` etc.  When
#'   robust rules are used the individual tree objects will be of
#'   class returned by ``robustcausalTree``.  Predictions work
#'   transparently via ``predict.causalForest``.
#' @export
robustCausalForest <- function(formula, data, treatment, na.action = na.causalTree,
                               split.Rule = "CT", double.Sample = TRUE, split.Honest = TRUE, split.Bucket = FALSE,
                               bucketNum = 5, bucketMax = 100, cv.option = "CT", cv.Honest = TRUE,
                               minsize = 2L, propensity, control, split.alpha = 0.5, cv.alpha = 0.5,
                               sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = .5,
                               mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees = nrow(data),
                               cost = FALSE, weights = FALSE, ncolx, ncov_sample, ...)

  {
  ##--------------------------------------------------------------------
  ## Determine which split rule is requested and whether it is robust.
  ## We include the new rules "MAD", "LMS" and "MSD" alongside the
  ## existing htetree rules.  Matching is case insensitive and uses
  ## pmatch for partial matching, consistent with causalTree.


  standard_rule_for_index <- c("TOT","CT","fit","tstats","TOTD","CTD","fitD","tstatsD",
                               "user","userD","policy","policyD")

  .robust_rules <- c("MAD","LMS","MSD")

  .map_robust_split_rule_int <- function(split.Rule) {
      if (!is.character(split.Rule) || length(split.Rule) != 1L)
          stop("`split.Rule` must be a single character value.")
      if (!split.Rule %in% .robust_rules)
          stop("robustcausalTree only supports split.Rule in {MAD, LMS, MSD}. For standard rules use htetree::causalTree().")
      length(.standard_rules_for_index) + match(split.Rule, .robust_rules)  # 13/14/15
  }



  ##--------------------------------------------------------------------
  ## Handle bucketised splitting.  For standard rules we append a
  ## trailing "D" to request the discrete variant.  Robust rules do
  ## not have discrete variants at the R level: bucketised behaviour
  ## should be implemented internally by robustcausalTree.
  # robust rules: either disable bucketing for now (recommended),
  # or implement robust bucketing inside the MAD/LMS/MSD C split code.
  if (isTRUE(split.Bucket)) {
      stop("split.Bucket=TRUE is not supported for MAD/LMS/MSD yet. Use htetree rules for bucketized splitting.")
  }

  split.Rule.int <- .map_robust_split_rule_int(split.Rule)

  ##--------------------------------------------------------------------
  ## Validate split.Honest and related parameters.  This closely
  ## follows the logic in causalTree.  We re‑compute the numeric
  ## index after possibly appending "D" so that the indices match
  ## those expected by the C code.  When new robust rules are added
  ## they should be appended to this vector in the same order as
  ## valid.rules above so that pmatch gives consistent indices.
  ## example below is converting the split rule into integer list like (1,2,3,...,15)
  all.rules.for.index <- c("TOT", "CT", "fit", "tstats",
                           "TOTD", "CTD", "fitD", "tstatsD",
                           "user", "userD", "policy", "policyD",
                           "MAD", "LMS", "MSD")
  split.Rule.int <- pmatch(split.Rule, all.rules.for.index)
  if (is.na(split.Rule.int)) {
    stop("Could not match split.Rule to index list")
  }

 # ----------------------------------------------------------------------------------
  ## Determine default honesty behaviour.  TOT/TOTD never uses
  ## honesty, so split.Honest is forced to FALSE for those rules.
  if (split.Rule.int %in% c(1, 5)) {
    if (!missing(split.Honest)) {
      warning("split.Honest is not used for TOT/TOTD splitting; setting to FALSE")
    }
    if (!missing(split.alpha)) {
      warning("split.alpha is not used for TOT/TOTD; forcing split.Honest to FALSE")
    }
    split.Honest <- FALSE
  } else {
    if (missing(split.Honest)) {
      split.Honest <- TRUE
      warning("The default split.Honest = TRUE for your chosen splitting rule.")
    }
  }
  ## Validate split.Honest logical
  split.Honest.num <- pmatch(split.Honest, c(TRUE, FALSE))
  if (is.na(split.Honest.num)) stop("Invalid split.Honest input; must be TRUE or FALSE")

  ## Adjust split.alpha and split.gamma based on honesty and rule index
  if (split.Honest && split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) {
    if (missing(split.alpha)) {
      split.alpha <- 0.5
    } else if (split.alpha < 0 || split.alpha > 1) {
      stop("split.alpha must be between 0 and 1")
    }
    if (missing(split.gamma)) {
      split.gamma <- 0.5
    } else if (split.gamma < 0 || split.gamma > 1) {
      stop("split.gamma must be between 0 and 1")
    }
  } else if (split.Rule.int %in% c(2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)) {
    ## dishonest (adaptive) splitting: force alpha to 1 and validate gamma
    if (split.alpha != 1) warning("For dishonest splitting, split.alpha is forced to 1")
    split.alpha <- 1
    if (missing(split.gamma)) {
      split.gamma <- 0.5
    } else if (split.gamma < 0 || split.gamma > 1) {
      stop("split.gamma must be between 0 and 1")
    }
  }

  ## Ensure cv.gamma and split.gamma exist.  If the user does not
  ## supply these parameters via '...', default to 0.5.  causalTree
  ## requires these when splitting or cross validating policies.  We
  ## intentionally set cv.gamma only after checking honesty so that
  ## user‐supplied values are respected.
  if (!exists("split.gamma", inherits = FALSE)) split.gamma <- 0.5
  if (!exists("cv.gamma", inherits = FALSE)) cv.gamma <- 0.5

  ##--------------------------------------------------------------------
  ## Propensity score sanity check.  Only applies to TOT/TOTD and fit
  ## rules.  Robust rules do not use propensity directly at this level.
  ## Add the robust rules in the below to get the propensity score.
  if (split.Rule %in% c("TOT", "TOTD", "fit", "fitD")) {
    if (!missing(propensity) && (propensity < 0 || propensity > 1)) {
      stop("Propensity must lie in [0, 1]")
    }
  }

  ##--------------------------------------------------------------------
  ## Pre‑process the formula and data.  This is largely copied from
  ## causalForest.  We construct a data frame containing only the
  ## covariates and response specified in the formula and append
  ## treatment.  The variable ``ncolx`` (the number of covariates) and
  ## ``ncov_sample`` (the number of covariates sampled for each tree)
  ## must be supplied by the user; otherwise an error is thrown.
  vars <- all.vars(formula)
  if (length(vars) < 2) {
    stop("Formula must have a response and at least one predictor")
  }
  yname <- vars[1]
  xnames <- vars[-1]
  treatmentdf <- data.frame(treatment)
  ## Subset data to only x and y.  Do not rely on causalTree's model.frame
  if (inherits(data, "data.table", TRUE) == 1) {
    ## convert to data.table and align treatment
    datax <- data[, ..xnames]
    datay <- data[, yname, with = FALSE]
    treatmentdt <- data.table(treatment)
    data <- cbind(datax, datay, treatmentdt)
  } else if (inherits(data, "data.frame")) {
    data <- data[, c(xnames, yname)]
    data <- cbind(data, treatmentdf)
  } else {
    stop("Data must be a data.frame or data.table")
  }
  num.obs <- nrow(data)

  ## The number of columns in the covariate matrix.  Must be
  ## specified by the caller because robustCausalForest does not know
  ## how many variables will be sampled otherwise.
  if (missing(ncolx)) {
    stop("Argument 'ncolx' (number of covariates) must be provided")
  }
  if (missing(ncov_sample)) {
    stop("Argument 'ncov_sample' (number of covariates sampled per tree) must be provided")
  }

  ## Limit the sample size used to construct each tree
  sample.size <- min(sample.size.total, num.obs)
  if (double.Sample) {
    train.size <- round(sample.size.train.frac * sample.size)
    est.size   <- sample.size - train.size
  }

  ## Initialise the causal forest object.  init.causalForest comes from
  ## htetree and sets up the storage for trees and inbag indicators.
  cf.obj <- init.causalForest(formula = formula, data = data,
                              treatment = treatment, weights = weights,
                              cost = cost, num.trees = num.trees,
                              ncov_sample = ncov_sample)

  ##--------------------------------------------------------------------
  ## Build each tree in turn.  We replicate the sampling logic of
  ## causalForest: draw a random subset of observations, split into
  ## training and estimation sets if honesty is requested, sample a
  ## subset of covariates, rename columns accordingly, and fit the
  ## appropriate tree.  When the robust rule is selected we call
  ## robustcausalTree; otherwise we call honest.causalTree or
  ## causalTree as in causalForest.
  for (tree.index in seq_len(num.trees)) {
    ## draw observations
    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)
    if (double.Sample) {
      train.idx       <- full.idx[1:train.size]
      reestimation.idx <- full.idx[(train.size + 1):sample.size]
    }
    ## randomise covariates to pick ncov_sample for this tree
    cov_sample <- sample.int(ncolx)
    cov_sample <- cov_sample[1:ncov_sample]
    ## construct formula part and record variable names
    fsample  <- ""
    name_all <- character(0)
    for (ii in seq_len(ncov_sample)) {
      nextx <- xnames[[cov_sample[ii]]]
      if (ii == 1L) {
        fsample <- nextx
        name_all <- c(name_all, nextx)
      } else {
        fsample <- paste0(fsample, "+", nextx)
        name_all <- c(name_all, nextx)
      }
    }

    ## name_all holds covariate names, append response and weight
    name_all <- c(name_all, yname, "w")
    ## Store sampling information into the forest object
    cf.obj$cov_sample[tree.index, ] <- cov_sample
    cf.obj$nameall_sample[tree.index, ] <- name_all
    cf.obj$fsample[[tree.index]] <- fsample
    ## Prepare data for tree fitting
    if (inherits(data, "data.table", TRUE) == 1) {
      ## training and estimation sets as data.tables
      if (double.Sample) {
        dataTree   <- data.table(data[train.idx, ])
        dataEstim  <- data.table(data[reestimation.idx, ])
      } else {
        dataTree   <- data.table(data[full.idx, ])
      }
      ## pick only sampled covariates plus response and treatment
      treeRange  <- c(cov_sample, (ncolx + 1):ncol(dataTree))
      dataTree  <- dataTree[, ..treeRange]
      if (double.Sample) {
        estimRange <- c(cov_sample, (ncolx + 1):ncol(dataEstim))
        dataEstim  <- dataEstim[, ..estimRange]
      }
    } else {
      ## data.frame case
      if (double.Sample) {
        dataTree   <- data.frame(data[train.idx, ])
        dataEstim  <- data.frame(data[reestimation.idx, ])
      } else {
        dataTree   <- data.frame(data[full.idx, ])
      }
      dataTree <- dataTree[, c(cov_sample, (ncolx + 1):ncol(dataTree))]
      if (double.Sample) {
        dataEstim <- dataEstim[, c(cov_sample, (ncolx + 1):ncol(dataEstim))]
      }
    }
    ## rename columns: sampled covariates, response, then weight
    names(dataTree) <- name_all
    if (double.Sample) names(dataEstim) <- name_all
    ## Compose the formula for the tree using the sampled covariates
    formula.tmp <- as.formula(paste0(yname, "~", fsample))
    ## Dispatch to the appropriate tree builder
    if (double.Sample) {
      if (is.robust) {
        ## robust tree with honesty
        tree.obj <- robustcausalTree(formula.tmp,
                                     data = dataTree,
                                     treatment = treatmentdf[train.idx, ],
                                     est_data = dataEstim,
                                     est_treatment = treatmentdf[reestimation.idx, ],
                                     split.Rule = split.Rule,
                                     split.Honest = split.Honest,
                                     split.Bucket = split.Bucket,
                                     bucketNum = bucketNum, bucketMax = bucketMax,
                                     cv.option = cv.option, cv.Honest = cv.Honest,
                                     minsize = nodesize, split.alpha = split.alpha,
                                     cv.alpha = cv.alpha, xval = 0,
                                     HonestSampleSize = est.size, cp = 0, ...)
      } else {
        ## standard honest tree
        tree.obj <- honest.causalTree(formula.tmp,
                                      data = dataTree,
                                      treatment = treatmentdf[train.idx, ],
                                      est_data = dataEstim,
                                      est_treatment = treatmentdf[reestimation.idx, ],
                                      split.Rule = split.Rule,
                                      split.Honest = split.Honest,
                                      split.Bucket = split.Bucket,
                                      bucketNum = bucketNum, bucketMax = bucketMax,
                                      cv.option = cv.option, cv.Honest = cv.Honest,
                                      minsize = nodesize, split.alpha = split.alpha,
                                      cv.alpha = cv.alpha, xval = 0,
                                      HonestSampleSize = est.size, cp = 0, ...)
      }
    } else {
      ## No double sampling: build a single tree on full.idx
      if (is.robust) {
        tree.obj <- robustcausalTree(formula.tmp,
                                     data = dataTree,
                                     treatment = treatmentdf[full.idx, ],
                                     split.Rule = split.Rule,
                                     split.Honest = split.Honest,
                                     split.Bucket = split.Bucket,
                                     bucketNum = bucketNum, bucketMax = bucketMax,
                                     cv.option = cv.option, cv.Honest = cv.Honest,
                                     x = FALSE, y = TRUE,
                                     split.alpha = split.alpha, cv.alpha = cv.alpha,
                                     split.gamma = split.gamma, cv.gamma = cv.gamma, ...)
      } else {
        tree.obj <- causalTree(formula.tmp,
                               data = dataTree,
                               treatment = treatmentdf[full.idx, ],
                               na.action = na.action,
                               split.Rule = split.Rule,
                               split.Honest = split.Honest,
                               split.Bucket = split.Bucket,
                               bucketNum = bucketNum, bucketMax = bucketMax,
                               cv.option = cv.option, cv.Honest = cv.Honest,
                               x = FALSE, y = TRUE,
                               split.alpha = split.alpha, cv.alpha = cv.alpha,
                               split.gamma = split.gamma, cv.gamma = cv.gamma, ...)
      }
    }
    ## Save the fitted tree and update inbag indicators
    cf.obj$trees[[tree.index]] <- tree.obj
    cf.obj$inbag[full.idx, tree.index] <- 1
    if (double.Sample) {
      cf.obj$inbag.Est[reestimation.idx, tree.index] <- 1
    }
  }
  cf.obj
}
