# lectometry.r
# functions to perform a lectometric analysis with a profile-based background

### initialize
### this function takes the input data t, the column in which the linguistic 
### functions in variables, the column in which the realizations reside of these 
### functions reside in variants, and the formula that is in the form of 
### variant ~ lectalfeature1 + lectalfeature2 + ... in formula
### output of this initialization fase is all the raw and relative frequencies,
### sizes of the lectal interactions, etc., everything except the aggregation

initialize <- function(t, variables, variants, formula){

	# extract the profile information, ie which variants for which variables
	# this may take a while
	print("creating profile information")
	variablesIndex <- which(colnames(t) == variables)
	variantsIndex <- which(colnames(t) == variants)

	fun = list()
	i <- 1
	for (variable in levels( t[,variablesIndex] )){
		fun[[i]] <- levels(as.factor(as.vector(
					t[as.vector(t[,variablesIndex] == variable), 
					variantsIndex])))
		i = i + 1
	} 

	# make the contingency table and its data frame according to the formula
	print("Creating data frame from formula")
	t.cont <- ftable(formula, data=t)
	t.cont.df <- as.data.frame(t.cont)

	# overview of the lectal features
	x.cat <- all.vars(formula,functions=FALSE)
	amountOfLin <- length(strsplit(as.character(formula)[2], "+", 
		fixed=TRUE)[[1]])
	x.ext <- x.cat[-c(1,amountOfLin)]
	x.rev <- lapply(t.cont.df[,x.ext],levels)
	x.fac <- rev(expand.grid(rev(x.rev)))

	# initialize the outputs
	print("Initialize output")
	raw.frequencies <- list() # contains the raw frequencies per variable
	rel.frequencies <- list() # contains the relative frequencies per variable
	# contains the size of the lects (sum of variable frequencies)
	lect.frequencies <- vector(length=nrow(x.fac))
	fe.pvalues <- list() # contains the fisher exact p values per variable
	llr.pvalues <- list() # contains the llr p values per variable
	# contains the variables and their variants
	variable.names <- fun 
	names(variable.names) <- levels(t.cont.df[,colnames(t.cont.df)==variables])
	lect.names <- c() # contains the names of the lectal combinations

	# go through the variables and gather the raw frequenties per variable 
	# across the lects
	print("Going through the variables, this will take a while")

	i <- 1
	for (variablename in names(variable.names)){ # variable per variable

		# grab the raw frequencies
		rawfreqs = c()
		for (variant in variable.names[[i]]){ # variant per variant
			# subset of the current variant of the variable
			ss <- t.cont.df[
				t.cont.df[, colnames(t.cont.df) == variables] == variablename & 
				t.cont.df[, colnames(t.cont.df) == variants] == variant, ]
			# keep the frequencies of this variant across the lects
			rawfreqs <- cbind(rawfreqs, ss$Freq)
		}

		# make sure to get the right lect names
		amountOfLin <- length(strsplit(as.character(formula)[2], "+", 
			fixed=TRUE)[[1]])
		x.ext <- all.vars(formula,functions=FALSE)[-c(1,amountOfLin)]
		nms <- ss[,x.ext]
		lect.names <- c()
		k <- 1
		while (k <= nrow(nms)){
			lect.name <- ""
			l <- 1
			while (l <= ncol(nms[k,])){
				lect.name <- paste(lect.name, nms[k,l])
				l = l + 1
			}
			lect.names <- cbind(lect.names, as.vector(lect.name))
			k = k + 1
		}
		lect.names <- as.vector(lect.names)
	
		# calculate the relative frequencies
		relfreqs <- prop.table(rawfreqs, 1)

		# give some names
		colnames(rawfreqs) <- variable.names[[i]]
		rownames(rawfreqs) <- lect.names
		colnames(relfreqs) <- variable.names[[i]]
		rownames(relfreqs) <- lect.names

		# calculate lect sizes
		lect.frequencies <- lect.frequencies + rowSums(rawfreqs)

		# calculate the fisher exact p value
		pmat <- matrix(nrow=length(lect.names), ncol=length(lect.names), 
			dimnames=list(lect.names, lect.names))
		for (var1 in rownames(rawfreqs)){
			for (var2 in rownames(rawfreqs)){
				p <- tryCatch(fisher.test(cbind(rawfreqs[var1,], 
						rawfreqs[var2,]))$p.value, 
						error=function(dealWithIt){1})
				pmat[var1, var2] <- p
			}
		}
		fe.pvalues[[i]] <- pmat
		
		# calculate the log likelihood ratio p value
		# LLR g test implementation from 
		# http://www.psych.ualberta.ca/~phurd/cruft/g.test.r
		source("g.test.r")
		pmat <- matrix(nrow=length(lect.names), ncol=length(lect.names), 
			dimnames=list(lect.names, lect.names))
		for (var1 in rownames(rawfreqs)){
			for (var2 in rownames(rawfreqs)){
				p <- tryCatch(g.test(cbind(rawfreqs[var1,], 
						rawfreqs[var2,]))$p.value, 
						error=function(dealWithIt){1})
				pmat[var1, var2] <- p
			}
		}
		llr.pvalues[[i]] <- pmat
		
		# fill the outputs
		raw.frequencies[[i]] <- rawfreqs
		rel.frequencies[[i]] <- relfreqs
	
		# increase the counter
		i <- i + 1
	}

	print("Finalizing output")
	out = list(data=t,
			   formula=formula,
			   all.frequencies=t.cont,
			   all.dataframe=t.cont.df,
			   raw.frequencies=raw.frequencies,
			   rel.frequencies=rel.frequencies,
			   lect.frequencies=lect.frequencies,
			   fe.pvalues=fe.pvalues,
			   llr.pvalues=llr.pvalues,
			   variable.names=variable.names,
			   lect.names=lect.names)

	return(out)	   

}

