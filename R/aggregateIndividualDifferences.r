
### aggregateIndividualDifferences aggregates the frequencies of variables to
### measure the distances between lects, but it does so in a transparent way
### by using individual differences scaling.

aggregateIndividualDifferences <- function(initialization, k=2){

	# calculate the city block distances (normalized between 0 and 1)
	cb.distances <- list()
	var.names <- c()
	i = 1
	j = 1
	for (relfreqs in initialization$rel.frequencies){
		mask <- apply(relfreqs, 2, is.nan)
		relfreqs[mask] <- 0
		d <- 0.5 * dist(relfreqs, method="manhattan")
		w <- initialization$llr.pvalues[[j]]
		# weighting
		w[w < 0.05] <- 0
		w[w >= 0.05] <- 1
		dw <- as.matrix(d) * (1-w)
		rownames(dw) <- rownames(w)
		colnames(dw) <- colnames(w)
		if (sum(dw) > 0) {
          cb.distances[[i]] <- as.dist(dw)
          var.names <- rbind(var.names, names(initialization$variable.names)[j])
	      i = i + 1
		}
		j = j + 1
	}
	
	# load smacof
	library(smacof)
	# find a three dimensional, non metric solution, with diagonal constraint
	fit <- smacofIndDiff(cb.distances, constraint="diagonal", ndim=k, 
      metric=F, verbose=TRUE, eps=0.0001, itmax=1000)

	# Group Stimulus Space
	gss <- fit$gspace

	# Configuration Weights
	cws <- c()
	for (k in fit$cweights){
		w <- diag(k)
		cws <- rbind(cws, w)
	}
	rownames(cws) <- var.names

	# prepare output
	out = list("distances"=cb.distances,
			   "fit"=fit,
			   "groupStimulusSpace"=gss,
			   "configurationWeights"=cws)
}
