# Vendored from josephwb/chronos (GPL-3) for downstream PCR use.
# Source: https://github.com/josephwb/chronos
# Included here so chronos bootstrap CI methods are available locally.
# Not used by the benchmark performance runner.

#########################
### Quantile Function ###
#########################

##' @title Quantile-based age estimation
##' @description Quantile function that computes the age corresponding to a
##' particular probability for the upper bound of a distribution of ages
##' (Strauss & Sadler 1989, Gingerich & Uhen 1998, Solow 2003). The method of
##' Strauss & Sadler (1989) assumes that the distribution of fossil ages is
##' uniform and their formula depends on the fossil ages range and the number of
##' fossil ages. The method of Solow (2003) is a general method for non-uniform
##' distributions and depends on the temporal gap between the oldest and the
##' second oldest fossil ages. Both methods assume that fossil ages are
##' independent samples from the same distribution (only relevant for the two
##' oldest ages for Solow's method), therefore, fossils should be as independent
##' as possible (ideally from different geological formations and different
##' regions).
##' 
##' In the particular case where there are only two fossil ages, Strauss &
##' Sadler's and Solow's methods converge to the same result; the quantile
##' functions are simply Xn/(1-P), and the likelihood function is 1/X.
##' @param p The desired probability level (0 < p < 1). A vector of
##' probabilities can be provided. Default p = 0.5
##' @param ages Either a vector of fossil ages or a matrix with two columns: the
##' first with the minimum age bounds (upper stratigraphic bounds) and the
##' second with the maximum age bounds (lower stratigraphic bounds) of each
##' fossil (in rows). A minimum of 2 are required
##' @param method The method for modelling age uncertainty. A number of options
##' is available:
##' \itemize{
##' \item{\code{"StraussSadler"}} {The (default) method of Strauss & Sadler 
##' (1989) assumes that the sampled fossil ages are uniformly distributed in 
##' time, and a warning is returned if a Kolmogorov-Smirnov test rejects the 
##' uniformity hypothesis.}
##' \item{\code{“Beta”}} {The method is a different parameterization of the
##' Strauss & Sadler method (Wang et al., 2009) that uses the qbeta function,
##' since the ratio between the observed maximum fossil age and the clade age
##' for the Strauss & Sadler model is distributed according to a Beta
##' distribution with parameters N and 1 (Wang et al. 2009); this should give
##' the same result as "StraussSadler".}
##' \item{\code{Solow}} {The method of Solow (2003) does not assume uniformity 
##' in sampled fossil ages, but is based on the two oldest ages only.}
##' \item{\code{NorrisPenG} or \code{NorrisGLin}} {The method of Norris et al.
##' (2015) based on the two oldest fossils only and the log-logistic
##' distribution. "NorrisPenG" is used when the precise phylogenetic placement
##' of fossils is not known, whereas "NorrisGLin" is used when one fossil from
##' each daughter lineage is used.}
##' \item{\code{"RobertsSolow"}} {The method of Roberts & Solow (2003) does not
##' assume uniformity in sampling and is based on the fact that the joint
##' distribution of the k oldest fossil ages can be modeled with the same
##' Weibull distribution.}
##' }
##' @param k The number of fossil ages to use in when
##' \code{method="RobertsSolow"}, in which only the k oldest fossils are used.
##' Default k = 5.
##' @return A numeric value (or vector of numeric values, if multiple p values
##' are provided) representing the age estimate of the clade origin given the
##' method a p value provided
##' @importFrom stats qbeta qlogis
##' @examples
##' \dontrun{
##'   # The following demonstrates how inferences depend on p and method
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5))
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Beta")
##'   qage(p=c(0.1, 0.5, 0.9), ages=c(54, 30, 25, 14, 5), method="Solow")
##'   }
##' @importFrom Rdpack reprompt
##' @references
##' \insertRef{Gingerich1998}{chronos}
##' 
##' \insertRef{Norris2015}{chronos}
##' 
##' \insertRef{Solow2003}{chronos}
##' 
##' \insertRef{Strauss1989}{chronos}
##' 
##' \insertRef{Wang2007}{chronos}
##' 
##' \insertRef{Wang2009}{chronos}
##' 
##' \insertRef{Wang2010}{chronos}
##' @keywords models
##' @export
qage <- function(p=0.5, ages,  method="StraussSadler", k=5) {
 	
	if (length(ages) < 2) {
		stop("More than one fossil age is needed for inference.", call.=FALSE)
	}
  if (p < 0 || p > 1.0) {
    stop("Valid values for p are: 0 < p < 1.")
  }
 	
 	ages <- sort(ages, decreasing=TRUE)
	
	if (method=="StraussSadler") {
		n <- length(ages)
		range <- max(ages) - min(ages)
		AGE <- range*(1-p)^-(1/n) + min(ages)
	} else if (method=="Beta") {
		n <- length(ages)
		range <- max(ages) - min(ages)
		X <- range/stats::qbeta(p=p, shape1=n, shape2=1, lower.tail=FALSE)
		AGE <- X + min(ages)
	} else if (method=="Solow") {
		AGE <- ages[1] + (ages[1]-ages[2])*p/(1-p)
	} else if (method=="NorrisPenG") {
		PenG <- ages[1] - ages[2]
		UltG <- exp(stats::qlogis(p=p, location=log(PenG)))
		AGE <- UltG + ages[1]
	} else if (method=="NorrisGLin") {
		GLin <- ages[1] - ages[2]
		UltG <- exp(stats::qlogis(p=p, location=log(GLin/2)))
		AGE <- UltG + ages[1]
	} else if (method=="RobertsSolow") {
		# Select the k oldest observations
	  if (k > length(ages)) {
	    print(paste0("NOTE: ", k, " oldest fossil specified, but only ",
        length(ages), " fossils available. Using all fossils in calculations."))
	    k <- length(ages)
	  }
	  ages <- ages[1:k]
		
		# Estimate the shape paremter of the joint Weibull distribution
		t <- numeric()
		
		for (i in 1:(k-2)) {
			t[i] <- (ages[1] - ages[k]) / (ages[1] - ages[i+1])
		}
		
		v <- 1/(k-1) * sum(log(t))
		
		# Calculate Su (Robert & Solow 2003, = c(alpha) Solow 2005 eq. 18)
		Su <- (-log(1-p)/k)^-v
		
		# But Su must be >= 1 (otherwise the estimate may be younger than older fossil) so:
		Su <- pmax(Su, 1); # use pmax instead of max to allow for vectors of p and Su
		
		# Calculate the quantile (Solow 2005 eq. 17) (Something seems wrong)	
		#AGE <- ( ages[1] - Su*ages[k]) / ( 1 - Su)
		
		# Calculate the quantile (Roberts & Solow 2003)
		AGE <- ages[1] + ( ages[1] - ages[k]) / ( Su - 1 )
	} else {
		stop("Method '", method, "' is not recognized.", call.=FALSE)
	}
	
	return(AGE)
}



##################################
### Random Clade Age Generator ###
##################################

# This function uses the quantile function above to generate random ages from the
# corresponding models.

##' @title Random clade age generator
##' @description This function uses the [qage()] function above to generate
##' random clade ages from the corresponding models.
##' @param n Number of samples to be drawn. Default n = 1000.
##' @param ages Either a vector of fossil ages or a matrix with two columns: the
##' first with the minimum age bounds (upper stratigraphic bounds) and the
##' second with the maximum age bounds (lower stratigraphic bounds) of each
##' fossil (in rows). A minimum of 2 are required.
##' @param ... Other options passed to [qage()].
##' @details If ages are known exactly, only minimum ages are used. If some or
##' all ages have uncertainty, typically upper and lower bounds defined by
##' bracketing geological strata, age maximums and minimums are set for each
##' fossil.
##' @return A numeric vector of length \code{n} representing simulated clade
##' ages.
##' @importFrom stats runif
##' @examples
##' \dontrun{
##'   # Example using the default method of Strauss & Sadler (1989)
##'   hist(rage(n=10000, ages=c(50, 30, 25, 14, 3.5)), breaks=500, freq=FALSE,
##'     xlim=c(50,120), xlab="Clade Age",
##'     main="Strauss & Sadler Example\nn=1000, ages=c(50, 30, 25, 14, 3.5)")
##'   graphics::curve(dexp(x-50, rate=0.07), col='blue', lwd=2, add=TRUE)
##'   graphics::legend("right", legend="dexp", col="blue", lty=1, lwd=2,
##'     box.lty=0)
##'   
##'   # Example using the method of Solow (2003)
##'   hist(rage(n=10000, ages=c(50, 30, 25, 14, 3.5), method="Solow"),
##'     breaks=500, freq=FALSE, xlim=c(50,120), xlab="Clade Age",
##'     main="Solow Example\nn=1000, min.ages=c(50, 30, 25, 14, 3.5)")
##'     # heavy on intermediate values
##'   graphics::curve(dexp(x-50, rate=0.04), col='blue', lwd=2, add=TRUE)
##'   graphics::curve(dlnorm(x-50, meanlog=log(50-30), sdlog=pi/sqrt(3)),
##'     col='orange', lwd=2, add=TRUE)
##'   graphics::legend("right", legend=c("dexp", "dlnorm"), col=c("blue",
##'     "orange"), lty=1, lwd=2, box.lty=0)
##'   
##'   # Strauss & Sadler (1989) method incorporating fossil age uncertainty
##'   hist(rage(n=10000, ages=cbind(c(50, 30, 25, 14, 3.5),
##'     c(56, 35, 25, 14, 6))), breaks=500, freq=FALSE, xlim=c(50,120), xlab="Age",
##'     main="Solow Example (with fossil age uncertainty); n=1000,
##'       ages=cbind(c(50, 30, 25, 14, 3.5), c(56, 35, 25, 14, 6))")
##'   graphics::curve(dlnorm(x-50, meanlog=2.6, sdlog=0.9), col='blue', lwd=2, add=TRUE)
##'   graphics::legend("right", legend="dlnorm", col="blue", lty=1, lwd=2, box.lty=0)
##'   }
##' @seealso [qage()]
##' @export
rage <- function(ages, n=1000, ...) {
	ages <- as.matrix(ages)

	if (ncol(ages) > 2) {
		stop("'ages' must be a vector or a two-column matrix.")
	}
	if (n < 0) {
	  stop("n must be a positive integer.")
	}

	# When fossil ages are known exactly
	if (ncol(ages) == 1) {
		rA <- qage(p=stats::runif(n), ages=ages[,1], ...)
	} else if (ncol(ages) == 2) {
		# When fossil ages are known from time intervals (max-min)

		# Check that bounds are correct for all fossils
		if (any(apply(ages, 1, diff) < 0)) {
			stop("Detected at least one fossil with older bound < younger bound.")
		}
	
		rA <- numeric(length=n)
		for(i in 1:n) {
			A <- stats::runif(nrow(ages), min=ages[,1], max=ages[,2])
			rA[i] <- qage(p=stats::runif(1), ages=A, ...)
		}
	}
	return(rA)
}



## chronosCI.R (2021-10-28)
## Copyright 2021 Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep

##' @title Confidence Intervals of Chronograms
##' @description Computes the confidence intervals of a chronogram using
##' different methods.
##' @param chronogram a chronogram output by [ape::chronos()].
##' @param pml.output an unrooted tree output by [phangorn::optim.pml()].
##' @param B the number of replications.
##' @param type the bootstrap method to be used. Can be either a character
##' string among "nonparametric", "semiparametric", or "parametric", or an
##' integer between 1 and 3.
##' @param calibration a data frame with the calibration points of the tree.
##' @param method a character string specifying the method used to assess
##' calibration uncertainty. See [qage()] for available methods.
##' @param control see [ape::chronos()]. By default, these are the same than
##' used to estimate the chronogram with [chronos()].
##' @param quiet a logical value. By default, the progress of the
##' computations is printed.
##' @param trees a logical value specifying whether to return the bootstrap
##' trees produced by phangorn (FALSE, by default).
##' @details The details of the methods are presented in the manuscript below.
##'
##' The labels (or taxa names) of the first argument (\code{chronogram})
##' must be present in the second argument (\code{pml.output}), and this
##' second one must also include the outgroup.
##' @return a matrix with the 50\% and 95\% lower and upper bounds of the
##' confidence intervals.
##' @references Paradis, E., Claramunt, S., Brown, J., and Schliep, K. Confidence
##' intervals in molecular dating by penalized likelihood. (in preparation)
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom ape chronos.control Nedge.phylo drop.tip root makeNodeLabel Ntip branching.times chronos
##' @importFrom stats rpois
##' @importFrom phangorn bootstrap.pml pml.control
##' @importFrom stats quantile
##' @examples
##' \dontrun{
##' ##---- Should be DIRECTLY executable !! ----
##' ##-- ==>  Define data, use random,
##' ##--	or do  help(data=index)for the standard data sets.
##' }
##' @seealso [qage()] [drawChronosCI()] [ape::chronos()]
##' [phangorn::bootstrap.pml()]
##' @keywords models
##' @export
chronosCI <-
    function(chronogram, pml.output, B = 100, type = "semiparametric",
             calibration = NULL, method = "StraussSadler",
             control = NULL, quiet = FALSE, trees = FALSE)
{
    getArgument <- function(name, default) {
        i <- match(name, args)
        if (is.na(i)) return(default)
        res <- thecall[[i]]
        if (!is.character(res)) eval(res) else res
    }
    thecall <- attr(chronogram, "call")
    args <- names(thecall) # the arguments used in the call
    model <- getArgument("model", "correlated")
    lambda <- getArgument("lambda", 1)
    if (is.null(control)) control <- getArgument("control", chronos.control())

    BOOTSTRAPMETHODS <- c("nonparametric", "semiparametric", "parametric")
    if (is.character(type))
        type <- pmatch(type, BOOTSTRAPMETHODS)
    if (is.na(type) || type > length(BOOTSTRAPMETHODS) || !is.numeric(type))
        stop("bootstrap method cannot be found")
    switch(type, {
        smoothRates <- FALSE
        Poisson <- FALSE
    }, {
        smoothRates <- TRUE
        Poisson <- FALSE
    }, {
        smoothRates <- FALSE
        Poisson <- TRUE
    })

    treeML <- pml.output$tree
    taxa <- treeML$tip.label
    ingroup <- chronogram$tip.label
    outgroup <- taxa[! taxa %in% ingroup]
    same.tip.set <- length(taxa) == length(ingroup) && setequal(taxa, ingroup)
    if ((length(outgroup) < 1 || all(is.na(outgroup))) && !same.tip.set)
        stop("no outgroup identified from the arguments")

    if (Poisson) {
        if (!quiet) cat("Generating Poisson branch lengths...")
        S <- length(attr(pml.output$data, "index"))
        tr0 <- treeML
        tr0$tip.label <- NULL
        TR <- vector("list", B)
        N <- ape::Nedge.phylo(treeML)
        for (i in 1:B) {
            tr0$edge.length <- stats::rpois(N, S * treeML$edge.length) / S
            TR[[i]] <- tr0
        }
        class(TR) <- "multiPhylo"
        attr(TR, "TipLabel") <- taxa
    } else {
        if (!quiet) cat("Running phangorn bootstrap...")
        TR <- phangorn::bootstrap.pml(pml.output, B, optNni = FALSE,
                            control = phangorn::pml.control(trace = 0L))
    }
    if (!quiet) cat("\n")

    ## root and drop the outgroup from the bootstrap trees when an explicit
    ## outgroup is present. For same-tip rooted trees used in PCR's parametric
    ## wrapper, keep the bootstrap trees as-is.
    if (same.tip.set) {
        if (Poisson) {
            for (i in seq_along(TR)) TR[[i]]$tip.label <- taxa
        }
        rTR <- TR
    } else {
        rTR <- lapply(TR, function(x) ape::drop.tip(ape::root(x, outgroup), outgroup))
    }

    ## make labels associated to the edges of the chronogram:
    tmp <- c(ingroup, ape::makeNodeLabel(chronogram, "md5sum")$node.label)
    original.edge.labels <- tmp[chronogram$edge[, 2]]
    original.ints <- chronogram$edge[, 2] > ape::Ntip(chronogram)
    ## id. for the 1st bootstrap tree:
    tr1 <- rTR[[1]]
    tmp <- c(tr1$tip.label, ape::makeNodeLabel(tr1, "md5sum")$node.label)
    boot.edge.labels <- tmp[tr1$edge[, 2]]
    boot.ints <- tr1$edge[, 2] > ape::Ntip(tr1)
    ## to reorder the results below
    needReorder <- !identical(original.edge.labels, boot.edge.labels)
    if (needReorder)
        o <- match(original.edge.labels[original.ints], boot.edge.labels[boot.ints])
    ## a finir....

    ## get bootstrap branch lengths:
    if (smoothRates && !Poisson) {
        EL <- sapply(rTR, "[[", "edge.length")
        NE <- nrow(EL)
        bootBL <- matrix(NA_real_, NE, B)
        for (i in 1:NE) {
            f <- getDensel(EL[i, ])
            F <- Cdensel(f)
            bootBL[i, ] <- rbl(B, F)
        }
        for (i in 1:B) rTR[[i]]$edge.length <- bootBL[, i]
    }

    ##
    CALS <- list()
    if (!is.null(calibration)) {
        calNodes <- unique(calibration$node)
        for (nd in calNodes) {
            i <- which(calibration$node == nd)
            if (length(i) == 1) next
            ages <- as.matrix(calibration[i, c("age.min", "age.max")])
            cals <- rage(ages, n = B, method = method)
            CALS <- c(CALS, list(cals))
            names(CALS)[length(CALS)] <- nd
        }
    } else {
        cal <- data.frame(node = length(rTR[[1]]$tip.label) + 1L,
                          age.min = 1, age.max = 1, soft.bounds = FALSE)
    }

    CHR <- vector("list", B)
    OK <- rep(FALSE, B)

    if (length(CALS)) {
        nodes2 <- as.integer(names(CALS))
        DUP <- duplicated.default(calibration$node)
    } else {
        cal <- calibration
    }

    boot_cal_base <- if (length(CALS)) calibration[!DUP, ] else cal

    for (i in 1:B) {
        if (!quiet) cat("\rRunning chronos bootstrap:", i, "/", B)
        if (length(CALS)) {
            cal <- boot_cal_base
            ages2 <- sapply(CALS, "[", i)
            cal[match(names(ages2), cal$node), 2:3] <- ages2
        } else {
            cal <- boot_cal_base
        }
        cal <- .chronos_remap_calibration_nodes(chronogram, rTR[[i]], cal)
        cal <- .chronos_prepare_calibration_for_ape(cal)
        if (is.null(cal) || !nrow(cal)) next

        chr_try <- try(
            ape::chronos(rTR[[i]], quiet = TRUE, model = model,
                         calibration = cal, control = control),
            silent = TRUE
        )
        if (!inherits(chr_try, "try-error")) {
            CHR[[i]] <- chr_try
            OK[i] <- TRUE
        }
    }
    if (!any(OK)) stop("chronosCI bootstrap produced no successful chronos replicates")
    ok_idx <- which(OK)
    target_sig <- .chronos_signature_index(chronogram)
    target_nodes <- as.integer(names(target_sig))
    BT <- matrix(
        NA_real_,
        nrow = length(target_nodes),
        ncol = length(ok_idx),
        dimnames = list(as.character(target_nodes), NULL)
    )
    for (j in seq_along(ok_idx)) {
        chr_j <- CHR[[ok_idx[[j]]]]
        bt_j <- ape::branching.times(chr_j)
        sig_j <- .chronos_signature_index(chr_j)
        node_by_sig_j <- stats::setNames(as.integer(names(sig_j)), sig_j)
        matched_nodes_j <- unname(node_by_sig_j[target_sig])
        keep_j <- is.finite(matched_nodes_j)
        if (any(keep_j)) {
            BT[keep_j, j] <- bt_j[as.character(matched_nodes_j[keep_j])]
        }
    }
    q_probs <- c(0.025, 0.25, 0.75, 0.975)
    CI <- sapply(
        seq_len(nrow(BT)),
        function(r) {
            vals <- BT[r, is.finite(BT[r, ])]
            if (!length(vals)) return(rep(NA_real_, length(q_probs)))
            as.numeric(stats::quantile(vals, probs = q_probs, names = FALSE))
        }
    )
    colnames(CI) <- rownames(BT)
    if (!quiet) cat("\nDone.\n")
    if (!trees) return(CI)
    out_trees <- CHR[OK]
    class(out_trees) <- "multiPhylo"
    list(CI = CI, trees = out_trees)
}

isUnimodal <- function(density)
{
    v <- rle(sign(diff(density$y)))$values
    if (length(v) > 3) return(FALSE)
    if (length(v) == 3 && all(v == c(1, 0, -1))) return(TRUE)
    if (length(v) == 2 && all(v == c(1, -1))) return(TRUE)
    FALSE
}

getDensel <- function(el)
{
    start <- stats::bw.nrd0(el)
    inc <- start / 10
    b <- start
    repeat {
        densel <- stats::density(el, bw = b)
        if (isUnimodal(densel)) {
            xsmall <- densel$x < 1e-8
            if (any(xsmall)) {
                w <- which(xsmall)
                densel$x <- densel$x[-w]
                total <- sum(densel$y)
                extra <- sum(densel$y[w])
                densel$y <- densel$y[-w]
                densel$y <- densel$y * total / (total - extra)
            }
            break
        }
        b <- b + inc
    }
    densel
}

Cdensel <- function(f)
{
    n <- length(f$y)
    list(x = (f$x[-1] + f$x[-n]) / 2,
         y = cumsum(diff(f$x) * (f$y[-1] + f$y[-n]) / 2))
}

## random branch lengths
rbl <- function(n, F) {
    z <- numeric(n)
    p <- stats::runif(n)
    N <- length(F$y)
    if (any(s <- p < F$y[1])) z[s] <- F$x[1]
    if (any(l <- p > F$y[N])) z[l] <- F$x[N]
    m <- which(!s & !l)
    if (!length(m)) return(z)
    for (k in m) {
        j <- which(F$y > p[k])[1]
        i <- j - 1L
        b <- (F$y[j] - F$y[i]) / (F$x[j] - F$x[i])
        a <- F$y[i] - b * F$x[i]
        z[k] <- (p[k] - a) / b
    }
    z
}



##' @title Draw Confidences Intervals on Phylogenies
##' @description This is a low-level plotting command to draw the confidence
##' intervals on the node of a tree as rectangles with coloured backgrounds.
##' @param CI output from [chronosCI()] or a similar matrix.
##' @param col95 colour used for the 95% intervals; by default: transparent
##' red.
##' @param col50 colour used for the 50% intervals; by default: transparent
##' blue.
##' @param height the height of the boxes.
##' @param legend a logical value.
##' @param \dots arguments passed to [graphics::legend()]
##' @details The matrix \code{CI} must have four rows and as many columns as the
##' number of nodes of the tree. The first and fourth rows give the lower and
##' upper bounds of the 95% confidence intervals. The second and third rows
##' give the lower and upper bounds of the 50% confidence intervals.
##' @return NULL
##' @author Emmanuel Paradis, Santiago Claramunt, Joseph Brown, Klaus Schliep
##' @importFrom graphics legend rect yinch
##' @examples
##' \dontrun{
##' ##---- Should be DIRECTLY executable !! ----
##' ##-- ==>  Define data, use random,
##' ##--	or do  help(data=index)for the standard data sets.
##' }
##' @seealso [chronosCI()]
##' @keywords aplot
##' @export
drawChronosCI <- function(CI, col95 = "#FF00004D", col50 = "#0000FF4D",
                          height = NULL, legend = TRUE, ...)
{
  lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  if (is.null(height)) height <- graphics::yinch(0.2)
  ## the nodes are renumbered => need to use labels
  direction <- lastPP$direction
  left_right <- FALSE
  if(direction=="rightwards" || direction=="leftwards"){
  	left_right <- TRUE
  	if(direction=="rightwards") CI <- max(lastPP$xx) - CI
  	if(direction=="leftwards") CI <- min(lastPP$xx) + CI
  	Y <- lastPP$yy - height / 2
  }
  else {
  	if(direction=="downwards") CI <- min(lastPP$yy) + CI
  	if(direction=="upwards") CI <- max(lastPP$yy) - CI
  	Y <- lastPP$xx - height / 2
  }
  L <- CI[1, ]
  R <- CI[4, ]
  B <- Y[ - seq_len(lastPP$Ntip)]
  T <- B + height
  if(left_right) graphics::rect(L, B, R, T, col = col95, border = NULL)
  else graphics::rect(B, L, T, R, col = col95, border = NULL)
  L <- CI[2, ]
  R <- CI[3, ]
  if(left_right) graphics::rect(L, B, R, T, col = col50, border = NULL)
  else graphics::rect(B, L, T, R, col = col50, border = NULL)
  if (!identical(legend, FALSE)) {
  	loc <- if (is.logical(legend)) "topleft" else legend
  	graphics::legend(loc, legend = c("95% CI", "50% CI"), pch = 22,
  					 pt.bg = c(col95, col50), col = c(col95, col50), ...)
  }
}


# Compatibility wrappers for PCR scripts ---------------------------------
# These wrappers preserve the previous helper entry points while routing the
# chronos side to the actual Paradis/Brown/Claramunt/Schliep bootstrap code.

.paradis_ci_matrix_to_df <- function(dated_tree, CI) {
  n_tip <- ape::Ntip(dated_tree)
  internal_nodes <- seq.int(n_tip + 1L, n_tip + dated_tree$Nnode)
  depth <- ape::node.depth.edgelength(dated_tree)
  root_depth <- max(depth[seq_len(n_tip)], na.rm = TRUE)
  node_age <- pmax(root_depth - depth, 0)
  if (!is.matrix(CI) || nrow(CI) != 4L) {
    stop('chronosCI() must return a 4 x N CI matrix.')
  }
  ci_nodes <- suppressWarnings(as.integer(colnames(CI)))
  if (length(ci_nodes) == ncol(CI) && all(is.finite(ci_nodes)) && all(ci_nodes %in% internal_nodes)) {
    out <- data.frame(
      node = ci_nodes,
      age = as.numeric(node_age[ci_nodes]),
      ci_lo = as.numeric(CI[1, ]),
      ci_hi = as.numeric(CI[4, ]),
      q50_lo = as.numeric(CI[2, ]),
      q50_hi = as.numeric(CI[3, ]),
      se = NA_real_,
      stringsAsFactors = FALSE
    )
    out <- out[match(internal_nodes, out$node), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  bt_nodes <- suppressWarnings(as.integer(names(ape::branching.times(dated_tree))))
  if (length(bt_nodes) != length(internal_nodes) || any(!is.finite(bt_nodes))) {
    bt_nodes <- internal_nodes
  }
  if (ncol(CI) == length(bt_nodes)) {
    return(data.frame(
      node = bt_nodes,
      age = as.numeric(node_age[bt_nodes]),
      ci_lo = as.numeric(CI[1, ]),
      ci_hi = as.numeric(CI[4, ]),
      q50_lo = as.numeric(CI[2, ]),
      q50_hi = as.numeric(CI[3, ]),
      se = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  if (!is.null(colnames(CI)) && all(nzchar(colnames(CI)))) {
    ci_nodes <- suppressWarnings(as.integer(colnames(CI)))
    keep <- is.finite(ci_nodes) & ci_nodes %in% internal_nodes
    if (!any(keep)) {
      stop('CI matrix column names do not map to internal nodes.')
    }
    CI <- CI[, keep, drop = FALSE]
    ci_nodes <- ci_nodes[keep]
    out <- data.frame(
      node = ci_nodes,
      age = as.numeric(node_age[ci_nodes]),
      ci_lo = as.numeric(CI[1, ]),
      ci_hi = as.numeric(CI[4, ]),
      q50_lo = as.numeric(CI[2, ]),
      q50_hi = as.numeric(CI[3, ]),
      se = NA_real_,
      stringsAsFactors = FALSE
    )
    missing_nodes <- setdiff(internal_nodes, out$node)
    if (length(missing_nodes)) {
      fill <- data.frame(
        node = missing_nodes,
        age = as.numeric(node_age[missing_nodes]),
        ci_lo = as.numeric(node_age[missing_nodes]),
        ci_hi = as.numeric(node_age[missing_nodes]),
        q50_lo = as.numeric(node_age[missing_nodes]),
        q50_hi = as.numeric(node_age[missing_nodes]),
        se = 0,
        stringsAsFactors = FALSE
      )
      out <- rbind(out, fill)
    }
    out <- out[match(internal_nodes, out$node), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  if (ncol(CI) != length(internal_nodes)) {
    stop('CI matrix column count does not match internal-node count.')
  }
  data.frame(
    node = internal_nodes,
    age = as.numeric(node_age[internal_nodes]),
    ci_lo = as.numeric(CI[1, ]),
    ci_hi = as.numeric(CI[4, ]),
    q50_lo = as.numeric(CI[2, ]),
    q50_hi = as.numeric(CI[3, ]),
    se = NA_real_,
    stringsAsFactors = FALSE
  )
}

.legacy_calibration_to_df <- function(calibration = NULL, cal_min = NULL, cal_max = NULL) {
  if (!is.null(calibration)) {
    calibration <- as.data.frame(calibration, stringsAsFactors = FALSE)
    if ("age_min" %in% names(calibration) && !("age.min" %in% names(calibration))) {
      names(calibration)[names(calibration) == "age_min"] <- "age.min"
    }
    if ("age_max" %in% names(calibration) && !("age.max" %in% names(calibration))) {
      names(calibration)[names(calibration) == "age_max"] <- "age.max"
    }
    if (!("soft.bounds" %in% names(calibration))) calibration$soft.bounds <- FALSE
    if ("node" %in% names(calibration)) calibration$node <- suppressWarnings(as.integer(calibration$node))
    if ("age.min" %in% names(calibration)) calibration$age.min <- as.numeric(calibration$age.min)
    if ("age.max" %in% names(calibration)) calibration$age.max <- as.numeric(calibration$age.max)
    return(calibration)
  }
  cal_min_names <- if (is.null(cal_min)) character(0) else names(cal_min)
  cal_max_names <- if (is.null(cal_max)) character(0) else names(cal_max)
  keys <- union(cal_min_names, cal_max_names)
  if (!length(keys)) return(NULL)
  nodes <- suppressWarnings(as.integer(keys))
  age_min <- vapply(keys, function(k) {
    val <- cal_min[[k]]
    if (length(val)) as.numeric(val[[1L]]) else NA_real_
  }, numeric(1))
  age_max <- vapply(keys, function(k) {
    val <- cal_max[[k]]
    if (length(val)) as.numeric(val[[1L]]) else age_min[[match(k, keys)]]
  }, numeric(1))
  data.frame(
    node = nodes,
    age.min = age_min,
    age.max = age_max,
    soft.bounds = FALSE,
    stringsAsFactors = FALSE
  )
}

.chronos_prepare_calibration_for_ape <- function(calibration) {
  if (is.null(calibration) || !nrow(calibration)) return(calibration)
  keep <- intersect(c("node", "age.min", "age.max", "soft.bounds"), names(calibration))
  out <- calibration[, keep, drop = FALSE]
  if (!("soft.bounds" %in% names(out))) out$soft.bounds <- FALSE
  out$node <- suppressWarnings(as.integer(out$node))
  out$age.min <- as.numeric(out$age.min)
  out$age.max <- as.numeric(out$age.max)
  out$soft.bounds <- as.logical(out$soft.bounds)
  out
}

.chronos_children_lookup <- function(phy) {
  total <- ape::Ntip(phy) + phy$Nnode
  kids <- vector("list", total)
  for (k in seq_len(nrow(phy$edge))) {
    kids[[phy$edge[k, 1L]]] <- c(kids[[phy$edge[k, 1L]]], phy$edge[k, 2L])
  }
  kids
}

.chronos_descendant_tip_labels <- function(phy, node, children = NULL) {
  if (is.null(children)) children <- .chronos_children_lookup(phy)
  n_tip <- ape::Ntip(phy)
  stack <- node
  out <- integer(0)
  while (length(stack)) {
    cur <- stack[[1L]]
    stack <- stack[-1L]
    if (cur <= n_tip) {
      out <- c(out, cur)
    } else {
      stack <- c(children[[cur]], stack)
    }
  }
  sort(phy$tip.label[unique(out)])
}

.chronos_signature_index <- function(phy) {
  n_tip <- ape::Ntip(phy)
  int_nodes <- seq.int(n_tip + 1L, n_tip + phy$Nnode)
  children <- .chronos_children_lookup(phy)
  sig <- vapply(
    int_nodes,
    function(nd) paste(.chronos_descendant_tip_labels(phy, nd, children), collapse = "\r"),
    FUN.VALUE = character(1)
  )
  names(sig) <- as.character(int_nodes)
  sig
}

.chronos_remap_calibration_nodes <- function(source_phy, target_phy, calibration) {
  if (is.null(calibration) || !nrow(calibration)) return(calibration)
  if (!inherits(source_phy, "phylo") || !inherits(target_phy, "phylo")) return(calibration)
  if (!setequal(source_phy$tip.label, target_phy$tip.label)) return(calibration)
  if (!("node" %in% names(calibration))) return(calibration)

  if (all(c("taxonA", "taxonB") %in% names(calibration))) {
    pair_nodes <- mapply(
      function(a, b) {
        if (!is.character(a) || !is.character(b)) return(NA_integer_)
        if (!(a %in% source_phy$tip.label) || !(b %in% source_phy$tip.label)) return(NA_integer_)
        nd <- ape::getMRCA(source_phy, c(a, b))
        if (is.null(nd) || !is.finite(nd)) return(NA_integer_)
        as.integer(nd)
      },
      calibration$taxonA,
      calibration$taxonB
    )
    use_pair_nodes <- is.finite(pair_nodes)
    if (any(use_pair_nodes)) {
      calibration$node[use_pair_nodes] <- as.integer(pair_nodes[use_pair_nodes])
    }
  }

  src_sig <- .chronos_signature_index(source_phy)
  tgt_sig <- .chronos_signature_index(target_phy)
  tgt_node_by_sig <- stats::setNames(as.integer(names(tgt_sig)), tgt_sig)
  src_nodes <- as.integer(calibration$node)
  src_key <- as.character(src_nodes)
  src_match <- src_sig[src_key]
  mapped <- unname(tgt_node_by_sig[src_match])
  keep <- is.finite(mapped)
  if (!all(keep)) {
    calibration <- calibration[keep, , drop = FALSE]
    mapped <- mapped[keep]
  }
  calibration$node <- as.integer(mapped)
  rownames(calibration) <- NULL
  calibration
}

.chronos_call_arg <- function(chronogram, name, default = NULL) {
  thecall <- attr(chronogram, "call")
  if (!is.language(thecall)) return(default)
  args <- names(thecall)
  i <- match(name, args)
  if (is.na(i)) return(default)
  res <- thecall[[i]]
  if (!is.character(res)) eval(res) else res
}

.annotate_chronogram_call <- function(dated_tree, model = NULL, lambda = NULL,
                                      nb_rate_cat = NULL, control = NULL) {
  if (!inherits(dated_tree, "phylo")) stop("dated_tree must be a phylo object.")
  if (is.null(model)) model <- .chronos_call_arg(dated_tree, "model", "correlated")
  if (is.null(lambda)) lambda <- .chronos_call_arg(dated_tree, "lambda", 1)
  if (is.null(control)) control <- .chronos_call_arg(dated_tree, "control", ape::chronos.control())
  if (!is.list(control)) control <- ape::chronos.control()
  if (identical(model, "discrete") && is.finite(as.numeric(nb_rate_cat))) {
    control$nb.rate.cat <- as.integer(nb_rate_cat)
  }
  attr(dated_tree, "call") <- substitute(
    chronos(
      phy = phy,
      lambda = LAMBDA,
      model = MODEL,
      quiet = TRUE,
      calibration = calib,
      control = CTRL
    ),
    list(LAMBDA = lambda, MODEL = model, CTRL = control)
  )
  dated_tree
}

.minimal_pml_output_proxy <- function(phy, n_sites = 1000L) {
  if (!inherits(phy, 'phylo')) stop('phy must be a phylo object.')
  n_sites <- as.integer(n_sites)
  if (!is.finite(n_sites) || n_sites < 1L) n_sites <- 1000L
  data_proxy <- structure(list(), index = seq_len(n_sites))
  list(tree = phy, data = data_proxy)
}

chronos_ci <- function(phy, dated_tree, pml_output = NULL, calibration = NULL,
                       B = 100L, type = 'parametric', method = 'StraussSadler',
                       control = NULL, quiet = FALSE, trees = FALSE,
                       n_sites = 1000L, cal_min = NULL, cal_max = NULL,
                       model = NULL, lambda = NULL, nb_rate_cat = NULL, ...) {
  calibration <- .legacy_calibration_to_df(calibration, cal_min, cal_max)
  dated_tree <- .annotate_chronogram_call(
    dated_tree,
    model = model,
    lambda = lambda,
    nb_rate_cat = nb_rate_cat,
    control = control
  )
  calibration <- .chronos_remap_calibration_nodes(phy, dated_tree, calibration)
  if (is.null(pml_output)) {
    if (!(identical(type, 'parametric') || identical(type, 3L))) {
      stop(
        'chronos_ci() wraps the josephwb/chronos bootstrap implementation. ',
        'Nonparametric and semiparametric modes require a phangorn pml.output object.'
      )
    }
    pml_output <- .minimal_pml_output_proxy(phy, n_sites = n_sites)
  }
  CI <- chronosCI(
    chronogram = dated_tree,
    pml.output = pml_output,
    B = B,
    type = type,
    calibration = calibration,
    method = method,
    control = control,
    quiet = quiet,
    trees = trees
  )
  if (is.list(CI) && !is.null(CI$CI)) {
    if (isTRUE(trees)) return(CI)
    CI <- CI$CI
  }
  .paradis_ci_matrix_to_df(dated_tree, CI)
}

treepl_ci <- function(...) {
  stop(
    'treepl_ci() has been removed from the local helper. ',
    'Use an external bootstrap workflow for treePL uncertainty and feed the resulting CI table into PCR explicitly.'
  )
}

ci_width_summary <- function(candidate, ci_df, source = 'bootstrap_ci') {
  need <- c('ci_lo', 'ci_hi')
  if (!all(need %in% names(ci_df))) {
    stop('CI data frame must contain ci_lo and ci_hi columns.')
  }
  widths <- as.numeric(ci_df$ci_hi) - as.numeric(ci_df$ci_lo)
  widths <- widths[is.finite(widths) & widths >= 0]
  if (!length(widths)) {
    return(data.frame(
      candidate = candidate,
      uncertainty_source = source,
      uncertainty_n_bars = 0L,
      uncertainty_mean_width_ma = NA_real_,
      uncertainty_median_width_ma = NA_real_,
      uncertainty_sd_width_ma = NA_real_,
      uncertainty_min_width_ma = NA_real_,
      uncertainty_max_width_ma = NA_real_,
      uncertainty_q25_width_ma = NA_real_,
      uncertainty_q75_width_ma = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    candidate = candidate,
    uncertainty_source = source,
    uncertainty_n_bars = length(widths),
    uncertainty_mean_width_ma = mean(widths),
    uncertainty_median_width_ma = stats::median(widths),
    uncertainty_sd_width_ma = if (length(widths) > 1L) stats::sd(widths) else 0,
    uncertainty_min_width_ma = min(widths),
    uncertainty_max_width_ma = max(widths),
    uncertainty_q25_width_ma = as.numeric(stats::quantile(widths, 0.25, names = FALSE, type = 8)),
    uncertainty_q75_width_ma = as.numeric(stats::quantile(widths, 0.75, names = FALSE, type = 8)),
    stringsAsFactors = FALSE
  )
}
