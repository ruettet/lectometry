### aggregateGIW aggregates the frequencies to calculate the distances between 
### the lects with the Gewichteter Identitätswert algorithm of Goebl 1983. The
### idea is that this is a binary check, naive (no semantic ties) but frequency
### sensitive and with an inverse frequency weight (high frequency = low weight,
### low frequency = high weight), based on speelman 2009

aggregateGIW <- function(initialization){
	print("Aggregate with the Gewichteter Identitätswert")
	i = 1
	outlist <- list()
	for (rawfreqs in initialization$raw.frequencies){
		ccov = sum(rowSums(rawfreqs) > 0)
		out <- matrix(nrow=length(initialization$lect.names), ncol=length(initialization$lect.names), 
			dimnames=list(initialization$lect.names, initialization$lect.names))
		for (var1 in rownames(rawfreqs)){
			for (var2 in rownames(rawfreqs)){
				d = c()
				Tvar1 <- rawfreqs[var1,] >= max(rawfreqs[var1,])
				Tvar2 <- rawfreqs[var2,] >= max(rawfreqs[var2,])			
				for (w1 in names(Tvar1[Tvar1 == TRUE])){
					for (w2 in names(Tvar1[Tvar2 == TRUE])){
						dc = 0
						if (var1 == var2){
							dc = 0
						}
						if (var1 != var2 & w1 != w2){
							dc = 1
						}
						if (var1 != var2 & w1 == w2){
							wccov = sum(rawfreqs[,w1] > 0)
							dc = (wccov - 1) / ccov
						}
						d = cbind(d, dc)
					}
				}	
				if (length(d) > 1){
					ad = mean(d)
				}
				else{
					ad = d[1]
				}
				out[var1,var2] <- ad
			}
		}
		outlist[[i]] <- out
		i = i + 1
	}
	aggregated <- Reduce("+", outlist) / length(outlist)
	return(as.dist(aggregated))
}

