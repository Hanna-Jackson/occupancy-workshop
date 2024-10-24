
expit <- function(x) 1/(1+exp(-x))
logit <- function(x) log(x/(1-x))

## Make a colour transparent 
makeTransparent <- function(..., alpha=0.5) {
   if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
   alpha = floor(255*alpha)  
   newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
   .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
   }
   apply(newColor, 2, .makeTransparent, alpha=alpha)
}


## Nice little pdf function
pdf.f <- function(f, file, ...) {
  cat(sprintf('Writing %s\n', file))
  pdf(file, ...)
  on.exit(dev.off())
  f()
}


## Set Colour Scheme
generic.col     <-"black" 

fire.col       <- "darkgoldenrod1"
dark.fire.col  <- "darkgoldenrod3"

forest.col        <- "darkslategray4"
dark.forest.col   <-"cadetblue4"
dark.forest.col.2 <-"darkslategray"


## ~~~ Making a flat JAGS.visit for plotting later ~~~ 
JAGS.visit.flat <- as.data.frame(rbind(JAGS.visit[,1,],JAGS.visit[,2,]))

JAGS.visit.flat$siteID <- rep(site.info$siteID, times=2)

visit.flat <- merge(JAGS.visit.flat,
                    JAGS.site[,c('siteID', 'sitetype','canopyopenness')],
                    by='siteID')



## Make a column of JAGS.site and visit.flat that will store colours for fire and forest 
JAGS.site$sitetypecol <- NA
JAGS.site[which(JAGS.site$sitetype==0),"sitetypecol"] <- fire.col
JAGS.site[which(JAGS.site$sitetype==1),"sitetypecol"] <- forest.col

JAGS.site$darksitetypecol <- NA
JAGS.site[which(JAGS.site$sitetype==0),"darksitetypecol"] <- dark.fire.col
JAGS.site[which(JAGS.site$sitetype==1),"darksitetypecol"] <- dark.forest.col

visit.flat$sitetypecol <- NA
visit.flat[which(visit.flat$sitetype==0),"sitetypecol"] <- fire.col
visit.flat[which(visit.flat$sitetype==1),"sitetypecol"] <- forest.col

## Subsetting JAGS.site to fire and forest for ease of use in plotting later 
fire   <- JAGS.site[JAGS.site$sitetype==0,]
forest <- JAGS.site[JAGS.site$sitetype==1,]

## Summarize trap hours out to the site level 
hoursout.site   <- apply(JAGS.visit[,,'HoursOut'],1,sum)
hoursout.fire   <- hoursout.site[JAGS.site$sitetype==0] 
hoursout.forest <- hoursout.site[JAGS.site$sitetype==1]

## For convenience, save the data we fed into JAGS in it's unscaled version too for plotting later: 
unscaled <- list(X                 = JAGS.obs,
                 nsite             = dim(JAGS.obs)[1],
                 nvisit            = dim(JAGS.obs)[2],
                 nspecies          = dim(JAGS.obs)[3],
                 ## Site
                 sitetype            = JAGS.site$sitetype,           
                 canopyopenness      = JAGS.site$canopyopenness,     
                 varcanopyopenness   = JAGS.site$varcanopyopenness,  
                 fireyear            = as.numeric(as.character(2019-(JAGS.site$fireyear))), 
                 openflowerabundsite = JAGS.site$openflowerabundsite,
                 flowerSRsite        = JAGS.site$flowerSRsite,
                 ## Visit
                 hoursout             = JAGS.visit[,,"HoursOut"],     
                 trapjulian           = JAGS.visit[,,"TrapJulian"],    
                 openflowerabundvisit = JAGS.visit[,,"openflowerabundvisit"],
                 flowerSRvisit        = JAGS.visit[,,"flowerSRvisit"]
                 )
