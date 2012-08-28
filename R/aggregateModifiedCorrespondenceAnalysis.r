

### aggregateModifiedCorrespondenceAnalysis implements the efforts of Plevoets
### 2008 to modify Correspondence Analysis to comply with the ideas of profile-
### based lectometry.

aggregateModifiedCorrespondenceAnalysis <- function(init, way=2, nf=2, b=3000, 
												cl=0.95, rel=TRUE, std=FALSE, 
												phi=FALSE) {

	# get the variables from the initialization
	formula=init$formula
	frame=init$data
	fun=init$variable.names

	if (nf!=2&&b>0) {	# Security check before the script is executed
		stop("Confidence regions necessarily assume nf=2 !")
		}
	if (b>0) {
		library(ellipse)
		# Confidence ellipses require the library "ellipse"
		x.rep <- list(LIN=list(),EXT=list())
		# List "x.rep" will contain the bootstrap coordinates
		}
	x.cat <- all.vars(formula,functions=FALSE)
	# Vector "x.cat" contains all the factors in the analysis
	x.lin <- x.cat[1]	# Object "x.lin" is the factor of the linguistic vars
	x.ext <- x.cat[-1]	# Vector "x.ext" contains all external factors
	x.rev <- lapply(frame[,x.ext],levels)
	# List "x.rev" contains the levels of all external factors
	x.fac <- rev(expand.grid(rev(x.rev)))
	# Object "x.fac" contains all combinations of the levels of the external 
	# factors
	x.way <- list()
	# List "x.way" will contain the logical indices for each main category and 
	# interaction
	for (n in 1:way) {
		if (n==0) {
			break
			}
		x.com <- combn(x.ext,n,simplify=TRUE)
		for (m in 1:ncol(x.com)) {
			x.int <- interaction(x.fac[,x.com[,m]],sep=".")
			for (lev in levels(x.int)) {
				x.way[[lev]] <- x.int==lev
				}
			}
		}
	# NOTE: function "multicorrespart" is still under development
	# Objects "x.fac" and "x.way" are metadata on the columns of the frequency 
	# table (below) They serve to identify the appropriate columns for each main 
	# category and interaction. In later versions of the script, they will 
	# probably be replaced with more elegant code
	x.for <- formula(paste(c("~",paste(x.cat,collapse="+")),collapse=""))
	# Formula "x.for" is a redefinition for the computations
	i <- 0			# Zeroth iteration for the coordinates in the actual sample
	while (i<=b) {		# Iteration from zero until b
		if (i>0) {	# Bootstrap resample
			x.dat <- frame[sample(1:nrow(frame),replace=TRUE),]
			}
		else {		# Actual sample
			x.dat <- frame
			}
		x.tab <- ftable(xtabs(formula=x.for,data=x.dat,drop.unused=FALSE),
						row.vars=x.lin)
		dimnames(x.tab) <- list(attr(x.tab,"row.vars")[[x.lin]],1:ncol(x.tab))
		# Frequency table "x.tab" with the linguistic variants as the rows
		# The columns are the combinations of all levels of the external factors
		x.dev <- matrix(nrow=nrow(x.tab),ncol=ncol(x.tab),
						dimnames=dimnames(x.tab))
		# Matrix "x.dev" will contain the numerator of the within-residuals
		for (sub in fun) {	# Iteration over the profiles
			s.tab <- x.tab[sub,]	# Sub-table
			s.exp <- outer(apply(s.tab,1,sum),apply(s.tab,2,sum))/sum(s.tab)
			# Expected values
			x.dev[sub,] <- s.tab - s.exp
			}
		t.row <- apply(x.tab,1,sum)	# Overall row totals
		t.col <- apply(x.tab,2,sum)	# Overall column totals
		if (phi) {	# Phi-square coordinates
			x.tot <- 1
			}
		else {		# Chi-square coordinates
			x.tot <- sum(x.tab)
			}
		x.chi <- sqrt(x.tot) * x.dev / outer(sqrt(t.row),sqrt(t.col))
		# Matrix "x.chi" contains the within-residuals
		x.chi[is.nan(x.chi)] <- 0
		# Security check: Within-residuals of empty rows or columns are set to 0
		if (rel) {	# Row weights
			w.row <- sqrt(1/t.row)
			}
		else {		# No weights
			w.row <- rep(1,nrow(x.tab))
			}
		w.row[is.infinite(1/t.row)] <- NaN
		# Security check: Empty rows are not to be weighed
		if (i==0) {	# Actual coordinates
			x.svd <- svd(x.chi)	# Singular Value Decomposition
			x.val <- x.svd$d^2
			names(x.val) <- 1:length(x.val)
			# Vector "x.val" contains the eigenvalues, used in the output
			x.var <- sqrt(x.val[1:nf])
			if (length(x.var)>1) {
				x.var <- diag(x.var)
				dimnames(x.var) <- list(1:nf,1:nf)
				}
			# Matrix "x.var" contains the singular values, used in the comp
			x.uuu <- x.svd$u[,1:nf]	# Matrix U
			x.vvv <- x.svd$v[,1:nf]	# Matrix V
			}
		x.row <- (w.row * x.chi) %*% x.vvv %*% solve(x.var)	# Matrix R
		dimnames(x.row) <- list(rownames(x.tab),1:nf)
		x.col <- matrix(nrow=length(x.way),ncol=nf,
						dimnames=list(names(x.way),1:nf))
		# Matrix "x.col" will be matrix C
		# Its rows are all external categories from "x.way"
		for (lab in rownames(x.col)) {	# Iteration over external categories
			c.sum <- apply(cbind(x.tab[,x.way[[lab]]],0),1,sum)
			# Vector "c.sum" summates the appropriate columns for category "lab"
			c.dev <- numeric(nrow(x.tab))
			names(c.dev) <- rownames(x.tab)
			# Vector "c.dev" will contain the numerator of the within-residuals
			c.wgt <- numeric(nrow(x.tab))
			names(c.wgt) <- rownames(x.tab)
			# Vector "c.wgt" will contain the profile-based column weights
			for (cub in fun) {	# Iteration over the profiles
				c.sub <- c.sum[cub]
				c.row <- t.row[cub]
				c.tot <- sum(c.sub)
				c.dev[cub] <- c.sub - c.tot*c.row/sum(c.row)
				c.wgt[cub] <- rep(1/c.tot,length(cub))
				}
			c.wgt[is.infinite(c.wgt)] <- 0
			# Security check: Weights of empty profiles are set to zero
			c.chi <- sqrt(x.tot) * c.dev / sqrt(t.row*sum(c.sum))
			# Vector "c.chi" contains the within-residuals for category "lab"
			c.chi[is.infinite(1/t.row)] <- 0
			# Security check: Within-residuals of empty rows are set to zero
			if (rel) {	# Column weights
				w.col <- sqrt(c.wgt)
				}
			else {		# No weights
				w.col <- rep(1,nrow(x.tab))
				}
			x.col[lab,] <- t(w.col * c.chi) %*% x.uuu %*% solve(x.var)
			# Coordinate for category "lab"
			}
		if (!std) {	# Principal (or standard) coordinates
			x.row <- x.row %*% x.var
			x.col <- x.col %*% x.var
			}
		if (i>0) {	# Bootstrap coordinates
			x.rep$LIN[[i]] <- as.matrix(x.row)
			x.rep$EXT[[i]] <- as.matrix(x.col)
			# LIN contains the coordinates of the linguistic variants
			# EXT contains the coordinates of the external categories
			}
		else {	# Actual coordinates
			x.out <- list(LIN=x.row,EXT=x.col,VAL=x.val)
			# List "x.out" will be the output
			}
		i <- i+1	# Update iteration
		}
	if (b>0) {	# Computation of the confidence regions
		# ACKNOWLEDGEMENT: the code below is from Greenacre (2007: 250-252)
		x.lab <- list(LIN=levels(frame[,x.lin]),EXT=names(x.way))
		# List "x.lab" contains all levels of all factors
		x.ell <- list(LIN=list(),EXT=list())
		x.hul <- list(LIN=list(),EXT=list())
		# List "x.ell" will contain the ellipses
		# List "x.hul" will contain the convex hulls
		for (j in names(x.rep)) {	# Iteration over "LIN" and "EXT"
			j.ell <- list()
			j.hul <- list()
			j.rep <- x.rep[[j]]
			for (k in x.lab[[j]]) {	# Iteration over the levels
				k.mat <- matrix(nrow=b,ncol=nf)	# Dummy matrix
				for (l in 1:b) {	# Iteration over the bootstrap replicates
					l.vec <- j.rep[[l]]
					k.mat[l,] <- l.vec[k,]
					}
				if (all(is.nan(k.mat))) {
					# Security check: Levels with no coordinates are skipped
					j.ell[[k]] <- NA
					j.hul[[k]] <- NA
					next
					}
				k.mat <- k.mat[!is.nan(apply(k.mat,1,sum)),]
				# Security check: Empty coordinates are removed
				k.ell <- ellipse(cov(k.mat),centre=apply(k.mat,2,mean),level=cl)
				# REMINDER: Function "ellipse" requires library "ellipse"
				repeat {
					k.pts <- chull(k.mat)
					if ((nrow(k.mat[-k.pts,])/b)<cl) {
						break
						}
					k.mat <- k.mat[-k.pts,]
					}
				k.hul <- k.mat[c(k.pts,k.pts[1]),]	# Convex hull
				dimnames(k.ell) <- list(1:nrow(k.ell),1:ncol(k.ell))
				dimnames(k.hul) <- list(1:nrow(k.hul),1:ncol(k.hul))
				j.ell[[k]] <- k.ell
				j.hul[[k]] <- k.hul
				}
			x.ell[[j]] <- j.ell
			x.hul[[j]] <- j.hul
			}
		x.out$ELLIPS <- x.ell	# ELLIPS contains the ellipses in the output
		x.out$HULL <- x.hul	# HULL contains the convex hulls in the output
		}
	x.out
	}

### aggregateEuclidean aggregates the frequencies to calculate the distances
### between the lects with the euclidean distance. This is a naive approach that
### ignores the semantic ties between the variants.

aggregateEuclidean <- function(t.init){
	dmat <- dist(t.init$all.frequencies, method="euclidean")
	rownames(dmat) <- t.init$lect.names
	colnames(dmat) <- t.init$lect.names
	return(dmat)
}


