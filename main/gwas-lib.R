#-------------------------------------------------------------
#-------------------------------------------------------------
getOperatingSystem <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin')
    os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(as.character (os))
}
#-------------------------------------------------------------
# Calculate the inflation factor from -log10 values
# It can fire warning, here they are hidign
#-------------------------------------------------------------
calculateInflationFactor <- function (scores)
{
	oldw <- getOption("warn")
	options(warn = -1)

	remove <- which(is.na(scores))
	if (length(remove)>0) 
		x <- sort(scores[-remove],decreasing=TRUE)
	else 
		x <- sort(scores,decreasing=TRUE)

	pvalues = 10^-x
	chisq <- na.omit (qchisq(1-pvalues,1))
	delta  = round (median(chisq)/qchisq(0.5,1), 3)

	options (warn = oldw)

	return (list(delta=delta, scores=x))
}



#-------------------------------------------------------------
# Adjust both pValues and threshold (for Bonferroni)
# Calculate threshold to decide SNPs significance
#-------------------------------------------------------------
adjustPValues <- function (level, pValues, method="FDR") 
{
	m <- length(pValues)
	if (method=="Bonferroni") {
		threshold = -log10(level/m)
	} else if (method=="FDR") {
		tmp <- cbind(pValues,calculateQValue (pValues))
		tmp <- tmp[order(tmp[,2]),]
		if (tmp[1,2] > level) {
			threshold <- -log10(tmp[1,1])*1.2
		} else {
			k <- max(which(tmp[,2] < level))
			threshold <- -log10(mean(tmp[k:(k+1),1]))
		}
	}

	return (list (threshold=threshold, pValues=pValues))
}


calculateQValue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
            print("ERROR: p-values not in valid range.")
            return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }

        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        if (pi0 <= 0) {
            print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
            return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
            idx <- sort.list(x)
            fc <- factor(x)
            nl <- length(levels(fc))
            bin <- as.integer(fc)
            tbl <- tabulate(bin)
            cs <- cumsum(tbl)
            tbl <- rep(cs, tbl)
            tbl[idx] <- tbl
            return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
            qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 
                1)
        }
        return(qvalue)
}

#-------------------------------------------------------------
# Print a log message with the parameter
#-------------------------------------------------------------
msg <- function (...) 
{
  messages = unlist (list (...))
  cat (">>>>", messages, "\n")
}

msgmsg <- function (...) 
{
  messages = unlist (list (...))
  cat ("\t>>", messages, "\n")
}

msgmsgmsg <- function (...)
{
  messages = unlist (list (...))
  message ("\t\t>", messages)
}

msgError <- function (...) {
		messages = unlist (list (...))
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
		cat (messages, "\n")
		cat (strrep("-", sum(sapply(messages, nchar))),"\n")
}

#-------------------------------------------------------------
# Add label to filename and new extension (optional)
#-------------------------------------------------------------
addLabel <- function (filename, label, newExt=NULL)  {
	nameext = strsplit (filename, split="[.]")
	name    = nameext [[1]][1] 
	if (is.null (newExt))
		ext     = nameext [[1]][2] 
	else
		ext     = newExt
	newName = paste0 (nameext [[1]][1], "-", label, ".", ext )
	return (newName)
}
#----------------------------------------------------------
# Util to print head of data
# Fast head, for debug only
#----------------------------------------------------------
view <- function (data, n=5,m=6) {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 || nrow (data) < 5) n = nrow(data)
		if (m==0 || ncol (data) < 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)
}
viewx <- function (data, n=5,m=6) {
	filename = deparse (substitute (data))
	name = paste (deparse (substitute (data)),":  ")
	if (is.null (dim (data))) {
		dimensions = paste (length (data))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (length (data) < 6) n = length(data)
		print (data[1:n])
	}else {
		dimensions = paste0 (unlist (dim (data)),sep=c(" x ",""))
		message (name, class(data), " : (", paste0 (dimensions),")")
		if (n==0 || nrow (data) <= 5) n = nrow(data)
		if (m==0 || ncol (data) <= 6) m = ncol(data)
		print (data[1:n,1:m])
	}
	write.csv (data, paste0("x-", filename, ".csv"), quote=F, row.names=F)

	quit()
}



