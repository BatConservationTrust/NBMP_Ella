# Code to carry out generalised addidtive models (GAMs)

# NB: This code has been adapted from Rachel Fewsters's GAM code which is available online at
# https://www.stat.auckland.ac.nz/~fewster/gams/R/
# and was developed to analyse breeding bird survey trends
# Fewster et al. 2000. Ecology

# NBMP code

# load mgcv package for GAM
library(mgcv)

## example for duabenton's hibernation counts
## Covariates = month, temperature. These were determined to have a significant effect previously

# read in file
hdaub<-read.csv("hdaub.csv", header=T, sep=",")
# change name of column 4 to count
colnames(hdaub)[4]<-"count"


indsp.func <- function(sp, dfvec = c(4, 7, 10, 15, 20, 33))
{
	if(length(sp$count[is.na(sp$count)]) > 0) stop(
			"no missing data allowed")


	fit.func <- function(dfval)
	{
                gam.df <- gam(count~as.factor(site) + s(year, fx=TRUE, k=dfval+1) +tpc1grp + as.factor(month), 
			family = poisson(link = log), data = sp)                  
                pred.df <- predict.gam(gam.df,type="terms")
                srow <- length(pred.df[1,])   
                
               Nentries <- length(sp$count)
                reqd.entries <- (1:Nentries)[!duplicated(sp$year)] #
                years.entries <- sp$year[!duplicated(sp$year)]
                years.pred <- pred.df[reqd.entries, srow]

                years.ord <- sort(years.entries)
                pred.ord <- years.pred[order(years.entries)]                
                ind.df <- exp(pred.ord)/exp(pred.ord[2])   
 
                
                
		cat("df ", dfval, " complete \n")
		ind.df

	}
	ind.all <- sapply(dfvec, fit.func)
	dimnames(ind.all) <- list(character(0), paste("df", as.character(dfvec),
		sep = ""))
	ind.all
}
inddaub<-indsp.func(daub,3,4,5,6,7,8)
inddaub
# shows the population indexes for each year using the different degrees of freedom

# here we chose 4 df = approx 0.3* the number of survey years

bootstrap.func <- function(sp, defree = 4)
{
  uniq.site <- unique(sp$site)
  Nentries <- length(sp$count)
  Nyears <- length(unique(sp$year))
  Nsites <- length(uniq.site)	#
  # take sample of sites to be included in the resample 
  # (bootstrap across sites):
  sam <- sample(uniq.site, replace = T)	#
  
  # elements.func lists the rows of the data frame sp that need to 
  # be included in the resample, for an individual site x in the sample sam:
  elements.func <- function(x)
  {
    (1:Nentries)[sp$site == x]
  }
  elements <- lapply(sam, elements.func)	#
  # elements is the set of rows of the data frame to be included in the
  # resample (with repetitions listed separately).  It is a ragged list. 
  dr.site <- rep(1:Nsites, as.vector(unlist(lapply(elements, length))))	#
  # dr.site is the vector of site levels for the data resample: these should
  # go from 1 to Nsites even though some sites have been included in the 
  # resample more than once, because we want to fit the model to the
  # repetitions as if they were different sites.  However, year and count
  # for the resample should be extracted directly from the original data
  # frame, sp. 
  # 
  elements <- as.vector(unlist(elements))
  data.resample <- data.frame(site = dr.site, year = sp$year[elements],
                              count = sp$count[elements], temp = sp$tpc1grp[elements], month = as.factor(month)[elements])	#
  
  # fit GAM on this df using MGCV library:
  dr.gam <- gam(count~ s(year, fx=TRUE, k=defree+1) + as.factor(site) + cloud + temp + month, family = 
                  poisson(link = log), data = data.resample)	#
  # dr.gam is the GAM on the data resample.
  #
  # Now make out a new data frame consisting of all years (in case some were
  # omitted from data.resample) and all sites in data.resample.
  # The predict feature doesn't seem to work unless all sites are included.
  
  year.base <- min(sp$year)-1
  dr.new <- expand.grid(year=year.base+1:Nyears, site=sort(unique(dr.site)), temp = unique(temp$data.resample), month = unique(month$data.resample))
  result <- predict.gam(dr.gam, newdata=dr.new, type = "terms")                
  srow <- length(result[1,])   
  #
  # srow is the row containing the smooth term in the predict object
  #
  pred.allyears <-  result[1:Nyears, srow]
  # *** For old MGCV versions, change line above to
  #     pred.allyears <-  result[srow, 1:Nyears]
  indices <- exp(pred.allyears)/exp(pred.allyears[2])        #
  indices
  
  
}


outer.boot.func <- function(sp, defree = 2, Nrep = 399, indexfile)
{
  # outer.boot.func: 
  # outer function for bootstrapping.  
  # Calls the function bootstrap.func that will compute a 
  # bootstrap resample, fit the GAM and find the indices.
  # 
  # sp is the dataframe eg daub
  # Nrep is number of bootstrap replicates; defree is degrees of freedom
  # for the GAM fit.  Indexfile is the name of a file in the working
  # directory (see below): it should be in quotes (" ").  
  # 
  # 
  # NOTE:
  # The file indexfile is a safeguard against premature termination of the
  # bootstrap routine: all results are written to this file as the 
  # function proceeds.  If the program crashes, the early results can be
  # salvaged by:
  #   sp.bootind.crash <- matrix(scan("indexfile"), nrow = k, byrow = T)
  # where k is the number of replicates completed before the crash.
  # The remaining replicates can then be made up by restarting the
  # oter.boot.func function.
  #
  # ** IMPORTANT NOTE **
  # When restarting any Splus simulation after a crash, be sure to reset the
  # random seed before resuming: otherwise the previous results (before the 
  # crash) will just be repeated.  The easiest way to reset the seed is
  # simply to generate some new random numbers: e.g. runif(100).
  
  #      
  # 
  # START OF FUNCTION
  # some housekeeping jobs first:
  
  if(!is.character(indexfile)) stop(
    "filenames have to be character strings")	#
  # the dataset sp is assumed to have NO MISSING ENTRIES:
  if(length(sp$count[is.na(sp$count)]) > 0) stop(
    "no missing data allowed")	#
  # start bootstrap loop:
  replicate.func <- function(i)
  {
    print(i)
    rep.indices <- bootstrap.func(sp, defree = defree)
    write(rep.indices, file = indexfile, ncolumns = 4, append = T)
    rep.indices
  }
  index.matrix <- sapply(1:Nrep, replicate.func)
  matrix(unlist(index.matrix), byrow = T, nrow = Nrep)
}
daub.bootind.399.2 <- outer.boot.func(daub, 2, 399, "daub.bootind.399.2")
#conf should be chosen so that (Nrep+1)*(1-conf)/2 is an integer eg. 399 + 1 = 400 * (1-0.95 = -0.5) = 200 / 2 = 100


diffop.func <- function(x, h, d, interval)
{
  # 
  # EXAMPLE: 
  # this function is usually called from within sp.plot via the
  # function indmat.derivmat.func
  # first some housekeeping
  if(is.na(match(d, c(2, 4, 6)))) stop("d must be 2, 4 or 6.")	#
  N <- length(x)	# x is recorded at N equally spaced points
  hvec <- 1:(N/2)
  hvec <- (hvec - 1) %/% (d/2)
  hvec[hvec > h] <- h
  if(N %% 2)
    hvec <- c(hvec, max(hvec), rev(hvec))	# (N odd)
  else hvec <- c(hvec, rev(hvec))
  hvec[hvec == 0] <- 1	#
  # similarly, dvec stores the actual values (2, 4, or 6) for the 
  # difference operator d along the series:
  dvec <- rep(d, N)
  if(d == 6) {
    dvec[c(2, N - 1)] <- 2
    dvec[c(3, N - 2)] <- 4
  }
  else if(d == 4) {
    dvec[c(2, N - 1)] <- 2
  }
  dvec[c(1, N)] <- 0	#
  # Now code the various difference operators approximating 2nd derivative, 
  # corresponding to d=2, 4, or 6. 
  # xinds is the set of indices of x for which the derivative is to be
  # calculated (usually all those with d taking a certain value).
  # hvals are the h-values corresponding to those indices in xinds. 
  diff.op2 <- function(xinds, hvals)
  {
    (x[xinds + hvals] - 2 * x[xinds] + x[xinds - hvals])/hvals^2
  }
  diff.op4 <- function(xinds, hvals)
  {
    ( - x[xinds + 2 * hvals] + 16 * x[xinds + hvals] - 30 * x[xinds
                                                              ] + 16 * x[xinds - hvals] - x[xinds - 2 * hvals])/(12 * 
                                                                                                                   hvals^2)
  }
  diff.op6 <- function(xinds, hvals)
  {
    (2 * x[xinds + 3 * hvals] - 27 * x[xinds + 2 * hvals] + 270 * x[
      xinds + hvals] - 490 * x[xinds] + 270 * x[xinds - hvals
                                                ] - 27 * x[xinds - 2 * hvals] + 2 * x[xinds - 3 * hvals
                                                                                      ])/(180 * hvals^2)
  }
  ##
  secderiv <- numeric(N)
  secderiv[0] <- 0
  secderiv[N] <- 0	#
  # apply the relevant difference operators throughout the series to get the
  # estimate secderiv of the second derivative vector:
  secderiv[dvec == 2] <- diff.op2((1:N)[dvec == 2], hvec[dvec == 2])
  secderiv[dvec == 4] <- diff.op4((1:N)[dvec == 4], hvec[dvec == 4])
  secderiv[dvec == 6] <- diff.op6((1:N)[dvec == 6], hvec[dvec == 6])	#
  # correct for interval^2 as described above:
  secderiv <- secderiv/(interval^2)
  secderiv
}

indmat.derivmat.func <- function(sp.bootind, h, d, interval)
{
  # indmat.derivmat.func:
  # Calculates the matrix of 2nd derivative estimates from the
  # matrix sp.bootind of bootstrapped abundance indices.
  # 
  # sp.bootind is a matrix of the form sp.bootind.Nrep.df 
  
  #
  # The output from this function would normally be named 
  # sp.bootderiv.Nrep.df, although this function is usually only called
  # from within sp.plot.  
  
  t(apply(sp.bootind, 1, diffop.func, h = h, d = d, interval = interval))
}


index.ci.func <- function(index.matrix, conf = 0.95)
{
  # index.ci.func: 
  # Given a matrix of bootstrapped index values of the form 
  # index.matrix = sp.bootind.Nrep.df, this function extracts 
  # the lower and upper (100*conf)% confidence limits for the abundance
  # indices. 
  #
  # conf should be chosen so that (Nrep+1)*(1-conf)/2 is an integer:
  # otherwise it will not be possible to find integers l and u such that the
  # l'th and u'th ordered values give us exactly a conf% confidence interval. 
  # 
  # syntax:  sp.ind.ci <- index.ci.func(sp.bootind.Nrep.df, 0.95)
  #
  # Usually called from within the function sp.plot. 
  # 
  Nyears <- ncol(index.matrix)
  Nrep <- nrow(index.matrix)	#
  # Nrep is the number of bootstrap replicates in the index matrix. 
  alp <- (1 - conf)/2	# 
  # Interpolation is not ideal, so discard some rows of the bootstrap
  # matrix if (Nrep+1)*alp is not an integer:
  if(abs((Nrep + 1) * alp - round((Nrep + 1) * alp)) > 1e-05)
    stop("need to discard rows of index.matrix, or change conf, so that (Nrep+1)*(1-conf)/2 is an integer."
    )
  
  lowerpt <- (Nrep + 1) * alp
  upperpt <- (Nrep + 1) * (1 - alp)
  
  # The confidence interval goes from the (Nrep+1)*alp 'th ordered
  # bootstrap value (low) to the (Nrep+1)*(1-alp) 'th ordered bootstrap
  # value (high). 
  inner.func <- function(yr)
  {
    sort(index.matrix[, yr])[c(lowerpt, upperpt)]
  }
  index.ci <- sapply(seq(1, Nyears), inner.func)
  dimnames(index.ci) <- list(c("lower", "upper"), NULL)
  index.ci
}



deriv.ci.func <- function(deriv.matrix, conf = 0.95)
{
  # deriv.ci.func: 
  # Given a matrix of bootstrapped 2nd derivative values (of the form 
  # deriv.matrix = sp.bootderiv.Nrep.df), this function extracts 
  # the lower and upper (100*conf)% confidence limits for the second
  # derivatives.
  #
  # conf should be chosen so that (Nrep+1)*(1-conf)/2 is an integer:
  # otherwise it will not be possible to find integers l and u such that the
  # l'th and u'th ordered values give us exactly a conf% confidence interval. 
  # 
  # syntax:  sp.deriv.ci <- deriv.ci.func(sp.bootderiv.Nrep.df, 0.95)
  #
  # However, this function is rarely called directly: almost always called
  # from within the function sp.plot. 
  # 
  Nyears <- ncol(deriv.matrix)
  Nrep <- nrow(deriv.matrix)	#
  # Nrep is the number of bootstrap replicates in the 2nd derivative 
  # matrix. 
  
  alp <- (1 - conf)/2	# 
  # Interpolation is not ideal, so discard some rows of the bootstrap
  # matrix if (Nrep+1)*alp is not an integer:
  if(abs((Nrep + 1) * alp - round((Nrep + 1) * alp)) > 1e-05)
    stop("need to discard rows of deriv.matrix, or change\n\tconf, so that (Nrep+1)*(1-conf)/2 is an integer."
    )
  lowerpt <- (Nrep + 1) * alp
  upperpt <- (Nrep + 1) * (1 - alp)
  
  # The confidence interval goes from the (Nrep+1)*alp 'th ordered
  # bootstrap value (low) to the (Nrep+1)*(1-alp) 'th ordered bootstrap
  # value (high). 
  inner.func <- function(yr)
  {
    sort(deriv.matrix[, yr])[c(lowerpt, upperpt)]
  }
  deriv.ci <- sapply(seq(1, Nyears), inner.func)
  dimnames(deriv.ci) <- list(c("lower", "upper"), NULL)
  deriv.ci
}


sp.plot <- function(sp, defree, sp.bootind, h, d, interval, conf = 0.95,
                    add = F, col = 2, cex = 2)
{
 
  # "sp" is the name of the species data frame, and must be a character
  #    string: e.g. "daub". 
  #    
  # "defree" is the df of the required GAM fit 
  # "sp.bootind" is a matrix of bootstrapped indices
  #   
  # See function diffop.func for explanation of "h", "d", and "interval". 
  # "conf" gives confidence level for the ABUNDANCE INDICES confidence
  #    intervals: *not* the confidence level for the changepoints.
  #    The changepoints are automatically calculated at the 5% significance
  #    level. This can be changed if required by manually editing the value
  #    "deriv.conf" in the function body. 
  # "add" determines whether to add the plot to a current plot, or
  #    whether to start a new plot.  The default is to start a new plot. 
  # "col": determines what colour the new plot 
  # "cex" determines the size of the text and changepoints on the plot:
  #    increase it for larger text, decrease it for smaller text. 
  #
  deriv.conf <- 0.95
  par(font = 3)	#
  # housekeeping:
  if(!is.character(sp))
    stop("first argument must be a character string.")
  dfname <- paste("df", defree, sep = "")
  sp.ind.ci <- index.ci.func(sp.bootind, conf = conf)
  Nyears <- length(sp.ind.ci[1,  ])
  yearvec <- seq(2011, 2015)	#
  #change yearvec to the actual years of the survey 
  #
  lower <- min(sp.ind.ci[1,  ])
  upper <- max(sp.ind.ci[2,  ])
  indsp <- eval(parse(text = paste("ind", sp, sep = "")))    
  if(!add) {
    plot(yearvec, indsp[, dfname], type = "l", ylim = c(lower, 
                                                        upper), xlab = "Year", ylab = "Index of Abundance (base year = 2012)", cex = cex, lwd = 2, bty
         = "l", font = 3)
    textint <- 0.2/1.7 * (upper - lower)	#
    # textint is for beautification and nothing else
    text(min(yearvec) - 1, upper + textint, "index", cex = cex, 
         font = 3)
    text(max(yearvec), lower - textint, "year", cex = cex, font = 3
    )
  }
  if(add) {
    par(col = col)
    lines(yearvec, indsp[, dfname], lwd = 2)
  }
  lines(yearvec, sp.ind.ci[1,  ], lty = 4, lwd = 2)
  lines(yearvec, sp.ind.ci[2,  ], lty = 4, lwd = 2)	#
  # Now calculate matrix of bootstrapped second derivatives:
  sp.bootderiv <- indmat.derivmat.func(sp.bootind = sp.bootind, h = h, d
                                       = d, interval = interval)
  sp.deriv.ci <- deriv.ci.func(sp.bootderiv, conf = deriv.conf)	#
  # Calculate the changepoints:
  # the lower changepoints are those where the *upturn* of the curve
  # is statistically significant (i.e. the lower confidence limit is 
  # greater than 0);
  # the upper changepoints are those where the *downturn* of the curve
  # is statistically significant (i.e. the upper confidence limit is 
  # less than 0);
   changepts.lower <- seq(1, Nyears)[sp.deriv.ci["lower",  ] > 0]
  changepts.upper <- seq(1, Nyears)[sp.deriv.ci["upper",  ] < 0]
  if((length(changepts.lower > 0) & (!is.na(changepts.lower[1])))) {
   points(changepts.lower + min(yearvec) - 1, indsp[   
   changepts.lower, dfname], pch = 1, cex = cex)
    changepts.lower <- changepts.lower + min(yearvec) - 1
    }
   else changepts.lower <- "none"
  if((length(changepts.upper > 0) & (!is.na(changepts.upper[1])))) {
   points(changepts.upper + min(yearvec) - 1, indsp[
    changepts.upper, dfname], pch = 16, cex = cex)
    changepts.upper <- changepts.upper + min(yearvec) - 1
    }
  else changepts.upper <- "none"
  par(col = 1)
  list(upturns = changepts.lower, downturns = changepts.upper)
}

sp.plot("daub", 2, daub.bootind.399.4, 1, 6, 1, conf = 0.95,
        add = F, col = 2, cex = 2)
abline(1,0)
#gives the population indicies for each year
inddaub

## to get the width of the confidence intervals in each year
daub.ind.ci <-index.ci.func(daub.bootind.399.4, 0.95)
daub.ind.ci

daub.bootderiv.399.4 <- indmat.derivmat.func(daub.bootind.399.4, 1, 6, 1)
daub.deriv.ci <- deriv.ci.func(daub.bootderiv.399.4, 0.95)
