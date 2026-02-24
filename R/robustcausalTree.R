#' Robust Causal Effect Regression and Estimation Trees
#'
#' Fit a \code{robustcausalTree} model to get an \code{rpart} object using
#' robust splitting rules like MAD.
#'
#' @inheritParams causalTree
#' @export
robustcausalTree <- function(formula, data, weights, treatment, subset,
                             na.action = na.causalTree,
                             split.Rule, split.Honest, HonestSampleSize, split.Bucket, bucketNum = 5,
                             bucketMax = 100, cv.option, cv.Honest, minsize = 2L,
                             x = FALSE, y = TRUE, propensity, control, split.alpha = 0.5,
                             cv.alpha = 0.5, cv.gamma = 0.5, split.gamma = 0.5,
                             cost, ...) {

    Call <- match.call()

    indx <- match(c("formula", "data", "weights", "subset"), names(Call), nomatch = 0L)
    if (indx[1] == 0L) stop("a 'formula' argument is required")

    temp <- Call[c(1L, indx)]
    temp$na.action <- na.action
    temp[[1L]] <- quote(stats::model.frame)
    names(treatment) <- rownames(data)
    m <- eval.parent(temp)
    treatment <- treatment[(rownames(m))]

    Terms <- attr(m, "terms")
    if (any(attr(Terms, "order") > 1L))
        stop("Trees cannot handle interaction terms")

    Y <- model.response(m)
    wt <- model.weights(m)
    if (any(wt < 0)) stop("negative weights not allowed")
    if (!length(wt)) wt <- rep(1, nrow(m))

    offset <- model.offset(m)
    X <- causalTree.matrix(m)
    nobs <- nrow(X)
    nvar <- ncol(X)

    if (missing(treatment)) stop("You should input the treatment status vector.")
    if (sum(treatment %in% c(0,1)) != nobs) stop("The treatment status should be 1 or 0 only.")
    if (sum(treatment) == 0 || sum(treatment) == nobs) stop("The data only contains treated or controlled cases.")
    if (missing(propensity)) propensity <- sum(treatment) / nobs

    if (missing(split.Rule)) {
        split.Rule <- "TOT"
        warning("The default split rule is 'TOT'.")
    }

    if (missing(split.Bucket)) {
        split.Bucket <- FALSE
        bucketNum <- 0
        bucketMax <- 0
    }

    split.Bucket.num <- pmatch(split.Bucket, c(TRUE, FALSE))
    if (is.na(split.Bucket.num)) stop("Invalid split.Bucket input.")

    if (!split.Bucket) {
        bucketNum <- 0
        bucketMax <- 0
    } else {
        if (missing(bucketMax)) bucketMax <- 100
        if (missing(bucketNum)) bucketNum <- 5
        split.Rule <- paste(split.Rule, 'D', sep = '')
    }

    # Add robust rules here to generate integers 13
    all.rules <- c("TOT", "CT", "fit", "tstats", "TOTD", "CTD", "fitD", "tstatsD",
                   "user", "userD", "policy", "policyD", "MAD")
    split.Rule.int <- pmatch(split.Rule, all.rules)

    if (is.na(split.Rule.int)) stop("Invalid splitting rule.")
    split.Rule <- all.rules[split.Rule.int]

    if (split.Rule.int %in% c(1, 5)) {
        if (!missing(split.Honest)) warning("split.Honest is not used in your chosen splitting rule.")
        if (!missing(split.alpha)) warning("split.alpha is not used. split.Honest set to FALSE")
        split.Honest <- FALSE
    } else {
        if (missing(split.Honest)) {
            split.Honest <- TRUE
            warning("The default split.Honest = TRUE for your chosen splitting rule.")
        }
    }

    split.Honest.num <- pmatch(split.Honest, c(T, F))
    if (is.na(split.Honest.num)) stop("Invalid split.Honest input.")

    # Include robust rules (13) in the honesty check
    honest.rules <- c(2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13)
    if (split.Honest == TRUE && split.Rule.int %in% honest.rules) {
        if (missing(split.alpha)) split.alpha <- 0.5 else if (split.alpha > 1 || split.alpha < 0) stop("Invalid split.alpha.")
        if (missing(split.gamma)) split.gamma <- 0.5 else if (split.gamma > 1 || split.gamma < 0) stop("Invalid split.gamma.")
    } else if (split.Rule.int %in% honest.rules) {
        if (split.alpha != 1) warning("For dishonest(adaptive) splitting, split.alpha = 1.")
        split.alpha <- 1
        if (missing(split.gamma)) split.gamma <- 0.5 else if (split.gamma > 1 || split.gamma < 0) stop("Invalid split.gamma.")
    }

    if (split.Rule %in% c("TOT", "TOTD", "fit", "fitD")) {
        if (propensity > 1 || propensity < 0) stop("Propensity score should be between 0 and 1.")
    }

    xvar <- apply(X, 2, var)
    method <- "anova"
    method.int <- 1

    if (missing(cv.option)) {
        warning("Miss 'cv.option', choose not to do cross validations.")
        cv.option <- "none"
        xval <- 0
    }

    if (missing(cv.Honest)) cv.Honest <- TRUE
    cv.Honest.num <- pmatch(cv.Honest, c(T, F))
    if (is.na(cv.Honest.num)) stop("Invalid cv.Honest.")

    if (cv.option %in% c("CT", "fit", "user", "policy")) {
        cv.option <- paste0(cv.option, ifelse(cv.Honest, 'H', 'A'))
    }

    cv.option.num <- pmatch(cv.option, c("TOT", "matching", "fitH", "fitA", "CTH", "CTA", "userH", "userA", "policyH", "policyA", "none"))
    if (is.na(cv.option.num)) stop("Invalid cv option.")

    if (cv.option.num %in% c(1, 2, 4, 6, 8)) {
        if (!missing(cv.alpha)) warning("cv.alpha is not used in your chosen cross validation method.")
    }

    if (missing(cv.alpha)) cv.alpha <- 0.5
    if (missing(cv.gamma)) cv.gamma <- 0.5
    if (missing(HonestSampleSize)) HonestSampleSize <- nobs

    HonestSampleSize <- as.integer(HonestSampleSize)

    init <- get(paste("htetree", method, sep = "."), envir = environment())(Y, offset, wt)
    ns <- asNamespace("htetree")
    if (!is.null(init$print)) environment(init$print) <- ns
    if (!is.null(init$summary)) environment(init$summary) <- ns
    if (!is.null(init$text)) environment(init$text) <- ns

    Y <- init$y
    xlevels <- .getXlevels(Terms, m)
    cats <- rep(0L, ncol(X))
    if (!is.null(xlevels)) cats[match(names(xlevels), colnames(X))] <- unlist(lapply(xlevels, length))

    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(rpart.control))
        indx <- match(names(extraArgs), controlargs, nomatch = 0L)
        if (any(indx == 0L)) stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 0L]), domain = NA)
    }

    controls <- causalTree.control(...)
    if (!missing(control)) controls[names(control)] <- control

    xval <- controls$xval
    if (is.null(xval) || (length(xval) == 1L && xval == 0L) || method == "user") {
        xgroups <- 0L
        xval <- 0L
    } else if (length(xval) == 1L) {
        control_idx <- which(treatment == 0)
        treat_idx <- which(treatment == 1)
        xgroups <- rep(0, nobs)
        xgroups[control_idx] <- sample(rep(1L:xval, length = length(control_idx)), length(control_idx), replace = F)
        xgroups[treat_idx] <- sample(rep(1L:xval, length = length(treat_idx)), length(treat_idx), replace = F)
    } else if (length(xval) == nobs) {
        xgroups <- xval
        xval <- length(unique(xgroups))
    } else {
        if (!is.null(attr(m, "na.action"))) {
            temp <- as.integer(attr(m, "na.action"))
            xval <- xval[-temp]
            if (length(xval) == nobs) {
                xgroups <- xval
                xval <- length(unique(xgroups))
            } else stop("Wrong length for 'xval'")
        } else stop("Wrong length for 'xval'")
    }

    if (missing(cost)) cost <- rep(1, nvar) else {
        if (length(cost) != nvar) stop("Cost vector is the wrong length")
        if (any(cost <= 0)) stop("Cost vector must be positive")
    }

    tfun <- function(x) if (is.matrix(x)) rep(is.ordered(x), ncol(x)) else is.ordered(x)
    labs <- sub("^`(.*)`$", "\\1", attr(Terms, "term.labels"))
    isord <- unlist(lapply(m[labs], tfun))

    storage.mode(X) <- "double"
    storage.mode(wt) <- "double"
    storage.mode(treatment) <- "double"
    minsize <- as.integer(minsize)

    # [Inference] The .Call targets the main C entry point for causalTree.
    # If you registered a new name in init.c, change "causalTree" to "robustcausalTree".
    ctfit <- .Call("causalTree",
                   ncat = as.integer(cats * !isord),
                   split_Rule = as.integer(split.Rule.int),
                   bucketNum = as.integer(bucketNum),
                   bucketMax = as.integer(bucketMax),
                   method = as.integer(method.int),
                   crossmeth = as.integer(cv.option.num),
                   crossHonest = as.integer(cv.Honest.num),
                   as.double(unlist(controls)),
                   minsize,
                   as.double(propensity),
                   as.integer(xval),
                   as.integer(xgroups),
                   as.double(t(init$y)),
                   X,
                   wt,
                   treatment,
                   as.integer(init$numy),
                   as.double(cost),
                   as.double(xvar),
                   as.double(split.alpha),
                   as.double(cv.alpha),
                   as.integer(HonestSampleSize),
                   as.double(cv.gamma))

    # The rest of the object formatting remains identical to causalTree
    nsplit <- nrow(ctfit$isplit)
    ncat <- if (!is.null(ctfit$csplit)) nrow(ctfit$csplit) else 0L

    if (nsplit == 0L) xval <- 0L

    numcp <- ncol(ctfit$cptable)
    temp <- if (nrow(ctfit$cptable) == 3L) c("CP", "nsplit", "rel error") else c("CP", "nsplit", "rel error", "xerror", "xstd")
    dimnames(ctfit$cptable) <- list(temp, 1L:numcp)

    tname <- c("<leaf>", colnames(X))
    splits <- matrix(c(ctfit$isplit[, 2:3], ctfit$dsplit), ncol = 5L,
                     dimnames = list(tname[ctfit$isplit[, 1L] + 1L],
                                     c("count", "ncat", "improve", "index", "adj")))
    index <- ctfit$inode[, 2L]

    nadd <- sum(isord[ctfit$isplit[, 1L]])
    if (nadd > 0L) {
        newc <- matrix(0L, nadd, max(cats))
        cvar <- ctfit$isplit[, 1L]
        indx <- isord[cvar]
        cdir <- splits[indx, 2L]
        ccut <- floor(splits[indx, 4L])
        splits[indx, 2L] <- cats[cvar[indx]]
        splits[indx, 4L] <- ncat + 1L:nadd

        for (i in 1L:nadd) {
            newc[i, 1L:(cats[(cvar[indx])[i]])] <- -as.integer(cdir[i])
            newc[i, 1L:ccut[i]] <- as.integer(cdir[i])
        }
        catmat <- if (ncat == 0L) newc else {
            cs <- ctfit$csplit
            ncs <- ncol(cs); ncc <- ncol(newc)
            if (ncs < ncc) cs <- cbind(cs, matrix(0L, nrow(cs), ncc - ncs))
            rbind(cs, newc)
        }
        ncat <- ncat + nadd
    } else catmat <- ctfit$csplit

    if (nsplit == 0L) {
        frame <- data.frame(row.names = 1L, var = "<leaf>", n = ctfit$inode[, 5L],
                            wt = ctfit$dnode[, 3L], dev = ctfit$dnode[, 1L],
                            yval = ctfit$dnode[, 4L], complexity = ctfit$dnode[, 2L],
                            ncompete = 0L, nsurrogate = 0L)
    } else {
        temp <- ifelse(index == 0L, 1L, index)
        svar <- ifelse(index == 0L, 0L, ctfit$isplit[temp, 1L])
        frame <- data.frame(row.names = ctfit$inode[, 1L], var = tname[svar + 1L],
                            n = ctfit$inode[, 5L], wt = ctfit$dnode[, 3L],
                            dev = ctfit$dnode[, 1L], yval = ctfit$dnode[, 4L],
                            complexity = ctfit$dnode[, 2L],
                            ncompete = pmax(0L, ctfit$inode[, 3L] - 1L),
                            nsurrogate = ctfit$inode[, 4L])
    }

    if (is.null(init$summary)) stop("Initialization routine is missing the 'summary' function")
    functions <- if (is.null(init$print)) list(summary = init$summary) else list(summary = init$summary, print = init$print)
    if (!is.null(init$text)) functions <- c(functions, list(text = init$text))
    if (method == "user") functions <- c(functions, mlist)

    where <- ctfit$which
    names(where) <- row.names(X)

    ans <- list(frame = frame, where = where, call = Call, terms = Terms,
                cptable = t(ctfit$cptable), method = method, control = controls,
                functions = functions, numresp = init$numresp)
    if (nsplit) ans$splits = splits
    if (ncat > 0L) ans$csplit <- catmat + 2L
    if (nsplit) ans$variable.importance <- importance(ans)

    if (y) ans$y <- Y
    if (x) { ans$x <- X; ans$wt <- wt }

    ans$ordered <- isord
    if (!is.null(attr(m, "na.action"))) ans$na.action <- attr(m, "na.action")
    if (!is.null(xlevels)) attr(ans, "xlevels") <- xlevels
    if (method == "class") attr(ans, "ylevels") <- init$ylevels
    class(ans) <- "rpart"

    if (ncol(ans$cptable) >= 4) ans$cptable[,4] <- ans$cptable[,4] / ans$cptable[1, 4]

    ans
}
