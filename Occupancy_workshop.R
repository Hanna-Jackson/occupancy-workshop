
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                   Occupancy Model Workshop
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##                                                 \\    //
##             Written by:                         (o>  <o)
##                  Hanna Jackson               \\_//)  (\\_//
##                                               \_/_)  (_\_/ 
##                                                _|_    _|_   

## Set your working directory, or start a new project if that's your preference
setwd("~/Dropbox/Bumble Bee Habitat/workshop")

## Prepare the data
## If you'd like to, you can take a look at what's in here but its not necessary
source('source/Prep_data.R')

## Load in the prepared data that was just created by sourcing
load(file='saved_objects/JAGS.data.RData', verbose=TRUE)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                   What data do we need?
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~ These are data frames that contained our data ~~~~~  

## Each bumble bee observation 
head(bumble.obs)

## Each of the sites where we surveyed for bumble bees 
head(site.info)

## Each of the visits to each of those sites (2 visits per site) 
head(visit.info)


## JAGS cannot read lists, and data frames are just fancy lists so
## JAGS cannot read data frames. We can only give it vectors,
## matricies, and arrays. So that's why the previous file took those
## data frames and made the following objects containing the same data
## just in a different format: 


## ~~~~ And here is that info stored in matricies & arrays instead ~~~~~

## 3 dimensional array of site x visit x species
JAGS.obs

## matrix of site x variable
JAGS.site

## 3 dimensional array of site x visit x variable 
JAGS.visit 


## ~~~~~ To use these arrays, all we'll really need to know is how to index them ~~~~~
JAGS.obs['2010.F4',,]
JAGS.obs[,'visit.2',]
JAGS.obs[,,'melanopygus']

JAGS.obs['2010.F4','visit.2','melanopygus']

JAGS.obs[20,,]
JAGS.obs[20,2,5]


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                   Prep for model run 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Download JAGS: https://sourceforge.net/projects/mcmc-jags/ 

## Lets make sure you have the right packages installed!
for (ii in  c("coda","rjags","R2jags")){
  if( isFALSE(ii %in% installed.packages()[,"Package"]) ){
    install.packages(ii)
  }
}

## Load those required packages 
library(coda)
library(rjags)
library(R2jags)

## My initial values - JAGS needs this 
my.inits <- function(){
  list(Z = array(1,
                 dim = c(length(dimnames(JAGS.obs)[["site"]]),
                         length(dimnames(JAGS.obs)[["species"]]))))
}

my.scale <- function(x){
   return( (x-mean(x))/sd(x) )
}

##Taking the log of open flower abundance 
JAGS.visit[,,"openflowerabundvisit"] <- log(JAGS.visit[,,"openflowerabundvisit"] + 1)
JAGS.site [,"openflowerabundsite"]   <- log(JAGS.site[,"openflowerabundsite"] + 1)

## Package everything into a list for use in JAGS! 
my.data <- list(X                 = JAGS.obs,
                nsite             = dim(JAGS.obs)[1],
                nvisit            = dim(JAGS.obs)[2],
                nspecies          = dim(JAGS.obs)[3],
                ## ~~~ Site Variables
                sitetype            = JAGS.site$sitetype,
                region              = JAGS.site$region,
                canopyopenness      = my.scale(JAGS.site$canopyopenness),
                latitude            = my.scale(JAGS.site$latitude),
                openflowerabundsite = my.scale(JAGS.site$openflowerabundsite),
                flowerSRsite        = my.scale(JAGS.site$flowerSRsite),
                fireyear            = as.numeric(as.character(2019-(JAGS.site$fireyear))), 
                ## ~~~ Visit Variables 
                hoursout             = my.scale(JAGS.visit[,,"HoursOut"]),     
                trapjulian           = my.scale(JAGS.visit[,,"TrapJulian"]),    
                openflowerabundvisit = my.scale(JAGS.visit[,,"openflowerabundvisit"]),
                flowerSRvisit        = my.scale(JAGS.visit[,,"flowerSRvisit"])
                )



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                   Run the first model 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Source the model file - loads my.model and my.params objects 
model <- 'Model_no_detection_predictors'

source(paste0('models/', model, ".R"))

## Run the model! 
bugs <- jags(data       = my.data,  
             inits      = my.inits,
             parameters.to.save = my.params, 
             model.file = my.model,
             n.iter     = 11000, ##110000, 
             n.burnin   = 1000, ##10000,
             n.thin     = 10, ##100,
             n.chains   = 3,
             working.directory = NULL)


## ~~~~~~~~~~~~~~ Save the run ~~~~~~~~~~~~~~

save(my.data, bugs, args, my.params,
     file = paste0('model_output/', model, '.RData'))

## Code to load what we just saved in case you overwrite/remove it
load(file = paste0('model_output/', model, '.RData'), verbose=TRUE)


## ~~~~~~~~~~~~~~ Look at the results ~~~~~~~~~~~


## JAGS gives us this summary of parameters out that we can look at:
summ    <- bugs$BUGSoutput$summary
columns <- c('mean','2.5%','97.5%', 'Rhat')
summ[,columns]

## It also gives us the non-summarized raw chains
sims.arr <- bugs$BUGSoutput$sims.array
sims.mat <- bugs$BUGSoutput$sims.matrix  


## Let's understand what those raw chains look like before we go any further: 

## From the output array 
canopy.chain <- sims.arr[,,'psi.canopyopenness']
plot(canopy.chain[,1],
     col='red', type='l',
     ylab='psi.canopyopenness',xlab='Iteration')
lines(canopy.chain[,2], col='blue')
lines(canopy.chain[,3], col='green')
abline(h=0)
abline(h=mean(canopy.chain), col='black', lty=2,lwd=3)


## From the output matrix 
head(sims.mat)

canopy.chain <- sims.mat[,'psi.canopyopenness']
plot(canopy.chain, type='l')
abline(h=mean(canopy.chain), col='red')
mean(canopy.chain)

## And we can see thsat the mean(canopy.chain) is the same as the
## estimate from the summary object! 
summ[,c("mean","2.5%",'97.5%')]


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                           Plotting 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Next we need to plot our model run!

## Prepare for plotting
##    Setting colours, loading pdf.f, cleaning up some objects to make
##    plotting easier 
source('source/Prepare_for_plotting.R')

## ~~~ And now to actually plotting! ~~~

## Just like we sourced code earlier in the data prep phase, we're
## going to source files that contain the actual plotting code.
## The only difference is that this time instead of the source command
## processing a lot of data, it's just going to create and load into
## our workshapce a huge plotting function. We can then call that
## function to plot our model. 

## The only catch to that is that we won't call these functions
## directly, we're going to call them via the pdf.f function that is
## in our workspace from when we sourced
## Prepare_for_plotting.R.

## All this pdf.f() function does is let us choose the file path and
## name, as well as the width, height, and any other plotting
## parameters that are in the function itself. When you run it it will
## save a pdf of the figure to. 


## ~~~ Figure 3: Data Boxplots ~~~  
source("source/Fig3_Data_Boxplots.R")
pdf.f(f = plot.data.boxplots,      ## The function loaded in when we sourced
      file = paste0('plots/distofdata_', model, '_.pdf'),     ## The file path you want your figure to save to 
      width = 7,                   ## PDF width
      height = 7)                  ## PDF height 

## ~~~ Figure 4: Effect sizes ~~~
source("source/Fig4_Effect_Sizes.R")
pdf.f(f=make.effect.size.fig,       
      file = paste0("plots/effectsize_", model, ".pdf"),  
      width =7,                         
      height = 5)

## ~~~ Figure 5: Occupancy vs Covariates ~~~
source("source/Fig5_Occupancy.R")
pdf.f(f = make.occ.figure,       
      file = paste0('plots/occ_',model,".pdf"),  
      width = 8,                         
      height = 8)  






## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##       And all of that again for the other model
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

model <- 'Model_all_detection_predictors'

source(paste0('models/', model, ".R"))

## Run the model! 
bugs <- jags(data       = my.data,  
             inits      = my.inits,
             parameters.to.save = my.params, 
             model.file = my.model,
             n.iter     = 11000, ##110000, 
             n.burnin   = 1000, ##10000,
             n.thin     = 10, ##100,
             n.chains   = 3,
             working.directory = NULL)


## ~~~~~~~~~~~~~~ Save the run ~~~~~~~~~~~~~~

save(my.data, bugs, args, my.params,
     file = paste0('model_output/', model, '.RData'))

## Code to load what we just saved in case you overwrite/remove it
load(file = paste0('model_output/', model, '.RData'), verbose=TRUE)

## ~~~~~~~~~~~~~~ Look at the results ~~~~~~~~~~~

summ    <- bugs$BUGSoutput$summary
columns <- c('mean','2.5%','97.5%', 'Rhat')
summ[,columns]

## Making objects 
sims.arr <- bugs$BUGSoutput$sims.array
sims.mat <- bugs$BUGSoutput$sims.matrix  


## ~~~~~~~~~~~~~~ Plotting ~~~~~~~~~~~
## Prepare for plotting
source('source/Prepare_for_plotting.R')

## ~~~ Figure 4: Effect sizes ~~~
source("source/Fig4_Effect_Sizes.R")
pdf.f(f=make.effect.size.fig,       
      file = paste0("plots/effectsize_", model, ".pdf"),  
      width =7,                         
      height = 5)

## ~~~ Figure 5: Occupancy vs Covariates ~~~
source("source/Fig5_Occupancy.R")
pdf.f(f = make.occ.figure,       
      file = paste0('plots/occ_',model,".pdf"),  
      width = 8,                         
      height = 8)  

## ~~~ Figure 6: Detection vs Covariates ~~~
source("source/Fig6_Detection.R")
pdf.f(f = make.det.fig,       
      file = paste0('plots/detection_',model,".pdf"),  
      width = 11.3,                         
      height = 7)


