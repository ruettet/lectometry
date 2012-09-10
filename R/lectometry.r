# lectometry.r
# functions to perform a lectometric analysis with a profile-based background

### initialize
### this function takes the input data t and the formula that is in the form of 
### variable + variant ~ lectalfeature1 + lectalfeature2 + ...
### Output of this initialization fase is all the raw and relative frequencies,
### sizes of the lectal interactions, etc., everything except the aggregation

initialize <- function(t, formula){

	# make sure that there are no empty levels
	t <- droplevels(t)

	# make the contingency table and its data frame according to the formula
	cat("Creating data frame from formula\n")
	t.cont <- ftable(formula, data=t)
	t.cont.df <- as.data.frame(t.cont)

	# overview of the lectal features
	x.cat <- all.vars(formula,functions=FALSE)
	amountOfLin <- length(strsplit(as.character(formula)[2], "+", 
		fixed=TRUE)[[1]])
	x.ext <- x.cat[-c(1,amountOfLin)]
	x.rev <- lapply(t.cont.df[,x.ext],levels)
	x.fac <- rev(expand.grid(rev(x.rev)))

	# in case of a profile-based approach, there are two linguistic factors
	if (amountOfLin == 2){

		# extract the profile information, ie which variants for which variables
		# this may take a while
		cat("creating profile information\n")
		x.lin <- strsplit(as.character(formula)[2], "+", fixed=TRUE)[[1]]
		variables <- gsub(" ","", x.lin[1] , fixed=TRUE)
		variants <- gsub(" ","", x.lin[2] , fixed=TRUE)
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

		# initialize the outputs
		cat("Initialize output\n")
		raw.frequencies <- list() # contains the raw frequencies per variable
		rel.frequencies <- list() # contains the relative frequencies per variable
		# contains the size of the lects (sum of variable frequencies)
		lect.frequencies <- vector(length=nrow(x.fac))
		fe.pvalues <- list() # contains the fisher exact p values per variable
		llr.pvalues <- list() # contains the llr p values per variable
		# contains the variables and their variants
		variable.names <- fun # let us put the fun in variable.names!
		names(variable.names) <- levels(t.cont.df[,colnames(t.cont.df)==variables])
		lect.names <- c() # contains the names of the lectal combinations

		# go through the variables and gather the raw frequenties per variable 
		# across the lects
		cat("Going through the variables, this will take a while\n")

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

		cat("Finalizing output\n")
		out = list(type="variable-based",
				   data=t,
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
	}

	if (amountOfLin == 1){

		# initialize the outputs
		cat("Initialize output\n")
		raw.frequencies <- t.cont
		rel.frequencies <- prop.table(t.init$all.frequencies,2)
		# contains the size of the lects (sum of variable frequencies)
		lect.frequencies <- rowSums(t.cont)
		variable.names <- c() 		# contains names of the variables
		lect.names <- rownames(t.cont) # contains the names of the lectal combinations	

		cat("Making the names\n")
	
		# make sure to get the right variable names
		x.lin <- strsplit(as.character(formula)[2], "+", fixed=TRUE)[[1]]
		variables <- gsub(" ","", x.lin[1] , fixed=TRUE)
		variablesIndex <- which(colnames(t) == variables)
		variable.names = levels(t[,variablesIndex])
		
		# make sure to get the right lect names
		nms <- unique(t[,x.ext])
		nmsout <- c()
		k <- 1 # amount of rows
		while (k <= nrow(nms)){
			l <- 1 # amount of columns
			nm <- ""
			while (l <= ncol(nms)){
				nm <- paste(nm,nms[k,l], sep=" ")
				l = l + 1
			}
			nmsout <- cbind(nmsout, as.vector(nm))
			k = k + 1
		}
		lect.names <- as.vector(nmsout)

		cat("Making the frequency tables\n")
		# make a nice table of the raw frequencies
		raw.frequencies <- t.cont[c(1:nrow(t.cont)),]
		colnames(raw.frequencies) <- variable.names
		rownames(raw.frequencies) <- lect.names
		
		# and derive a nice relative frequencies table
		rel.frequencies <- prop.table(raw.frequencies, 2)
	
		# prepare output
		cat("Finalizing output\n")
		out = list(type="word-based",
				   data=t,
				   formula=formula,
				   all.frequencies=t.cont,
				   all.dataframe=t.cont.df,
				   raw.frequencies=raw.frequencies,
				   rel.frequencies=rel.frequencies,
				   lect.frequencies=lect.frequencies,
				   variable.names=variable.names,
				   lect.names=lect.names)		
	}

	return(out)	   

}

