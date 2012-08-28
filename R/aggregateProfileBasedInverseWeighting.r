### aggregateProfileBasedInverseWeighting is a function that performs 
### aggregation according to the profile based approach in Speelman et al 2003
### but with inverse weighting, boosting low frequency profiles as in GIW

aggregateProfileBasedInverseWeighting <- function(initialization, dstat="llr", 
	cweight="w"){
	print("Aggregate in the profile based way")
	
	# calculate the city block distances (normalized between 0 and 1)
	cb.distances <- list()
	i = 1
	for (relfreqs in initialization$rel.frequencies){
		mask <- apply(relfreqs, 2, is.nan)
		relfreqs[mask] <- 0
		cb.distances[[i]] <- 0.5 * dist(relfreqs, method="manhattan")
		i = i + 1
	}
	
	# calculate conceptual weight (speelman et al 2003)
	w.cweight <- list()
	i <- 1
	cweightmat <- matrix(nrow=length(initialization$lect.names), 
					ncol=length(initialization$lect.names), 
					dimnames=list(initialization$lect.names,
								  initialization$lect.names))				  
	while (i <= length(initialization$raw.frequencies)){
		for (var1 in initialization$lect.names){
			for (var2 in initialization$lect.names){
				sumx <- initialization$lect.frequencies[var1] + initialization$lect.frequencies[var2] # size of lects
				x <- sum( as.vector(initialization$raw.frequencies[[i]][var1,]) ) + sum( as.vector(initialization$raw.frequencies[[i]][var2,]) )
				amountvariants <- length(initialization$raw.frequencies)
				cweight <- (sumx - x) / (sumx * (amountvariants - 1))
				cweightmat[var1,var2] <- cweight
			}
		}
		w.cweight[[i]] <- cweightmat
		i <- i + 1
	}
	
	# user input: which statistical measure
	if (dstat == "llr"){
		pvalues <- initialization$llr.pvalues
	}
	if (dstat == "fe"){
		pvalues <- initialization$fe.pvalues
	}

	cweightmats <- w.cweight

	variable.names <- initialization$variable.names
	lect.names <- initialization$lect.names

	# initialize output
	out = matrix(ncol=length(initialization$lect.names), 
				 nrow=length(initialization$lect.names))
	out[is.na(out)] <- 0

	# calculate aggregate distance
	print("Aggregating the distances")
	i = 1
	while (i <= length(cb.distances)){
		# turn the pvalues into a filter
		pvaluesFilter <- initialization$llr.pvalues[[i]]
		pvaluesFilter[pvaluesFilter < 0.05 ] <- 1
		pvaluesFilter[pvaluesFilter >= 0.05 & pvaluesFilter < 1] <- 0

		# distances after weighing and filtering
		distance <- as.matrix(cb.distances[[i]])*cweightmats[[i]]*pvaluesFilter
		distance[is.na(distance) ] <- 0 
		out = out + distance		
		i = i + 1 
	}

	return(as.dist(out))

}

