\name{initialize}
\alias{initialize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initializes the workspace by taking a long format dataset and turning it into a contingency table
}
\description{
This function is the starting point of the lectometry package. The other functions rely on its output.
}
\usage{
initialize(t, variables, variants, formula)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{variables}{
%%     ~~Describe \code{variables} here~~
}
  \item{variants}{
%%     ~~Describe \code{variants} here~~
}
  \item{formula}{
%%     ~~Describe \code{formula} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (t, variables, variants, formula) 
{
    print("creating profile information")
    variablesIndex <- which(colnames(t) == variables)
    variantsIndex <- which(colnames(t) == variants)
    fun = list()
    i <- 1
    for (variable in levels(t[, variablesIndex])) {
        fun[[i]] <- levels(as.factor(as.vector(t[as.vector(t[, 
            variablesIndex] == variable), variantsIndex])))
        i = i + 1
    }
    print("Creating data frame from formula")
    t.cont <- ftable(formula, data = t)
    t.cont.df <- as.data.frame(t.cont)
    x.cat <- all.vars(formula, functions = FALSE)
    amountOfLin <- length(strsplit(as.character(formula)[2], 
        "+", fixed = TRUE)[[1]])
    x.ext <- x.cat[-c(1, amountOfLin)]
    x.rev <- lapply(t.cont.df[, x.ext], levels)
    x.fac <- rev(expand.grid(rev(x.rev)))
    print("Initialize output")
    raw.frequencies <- list()
    rel.frequencies <- list()
    lect.frequencies <- vector(length = nrow(x.fac))
    fe.pvalues <- list()
    llr.pvalues <- list()
    variable.names <- fun
    names(variable.names) <- levels(t.cont.df[, colnames(t.cont.df) == 
        variables])
    lect.names <- c()
    print("Going through the variables, this will take a while")
    i <- 1
    for (variablename in names(variable.names)) {
        rawfreqs = c()
        for (variant in variable.names[[i]]) {
            ss <- t.cont.df[t.cont.df[, colnames(t.cont.df) == 
                variables] == variablename & t.cont.df[, colnames(t.cont.df) == 
                variants] == variant, ]
            rawfreqs <- cbind(rawfreqs, ss$Freq)
        }
        amountOfLin <- length(strsplit(as.character(formula)[2], 
            "+", fixed = TRUE)[[1]])
        x.ext <- all.vars(formula, functions = FALSE)[-c(1, amountOfLin)]
        nms <- ss[, x.ext]
        lect.names <- c()
        k <- 1
        while (k <= nrow(nms)) {
            lect.name <- ""
            l <- 1
            while (l <= ncol(nms[k, ])) {
                lect.name <- paste(lect.name, nms[k, l])
                l = l + 1
            }
            lect.names <- cbind(lect.names, as.vector(lect.name))
            k = k + 1
        }
        lect.names <- as.vector(lect.names)
        relfreqs <- prop.table(rawfreqs, 1)
        colnames(rawfreqs) <- variable.names[[i]]
        rownames(rawfreqs) <- lect.names
        colnames(relfreqs) <- variable.names[[i]]
        rownames(relfreqs) <- lect.names
        lect.frequencies <- lect.frequencies + rowSums(rawfreqs)
        pmat <- matrix(nrow = length(lect.names), ncol = length(lect.names), 
            dimnames = list(lect.names, lect.names))
        for (var1 in rownames(rawfreqs)) {
            for (var2 in rownames(rawfreqs)) {
                p <- tryCatch(fisher.test(cbind(rawfreqs[var1, 
                  ], rawfreqs[var2, ]))$p.value, error = function(dealWithIt) {
                  1
                })
                pmat[var1, var2] <- p
            }
        }
        fe.pvalues[[i]] <- pmat
        source("g.test.r")
        pmat <- matrix(nrow = length(lect.names), ncol = length(lect.names), 
            dimnames = list(lect.names, lect.names))
        for (var1 in rownames(rawfreqs)) {
            for (var2 in rownames(rawfreqs)) {
                p <- tryCatch(g.test(cbind(rawfreqs[var1, ], 
                  rawfreqs[var2, ]))$p.value, error = function(dealWithIt) {
                  1
                })
                pmat[var1, var2] <- p
            }
        }
        llr.pvalues[[i]] <- pmat
        raw.frequencies[[i]] <- rawfreqs
        rel.frequencies[[i]] <- relfreqs
        i <- i + 1
    }
    print("Finalizing output")
    out = list(data = t, formula = formula, all.frequencies = t.cont, 
        all.dataframe = t.cont.df, raw.frequencies = raw.frequencies, 
        rel.frequencies = rel.frequencies, lect.frequencies = lect.frequencies, 
        fe.pvalues = fe.pvalues, llr.pvalues = llr.pvalues, variable.names = variable.names, 
        lect.names = lect.names)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
