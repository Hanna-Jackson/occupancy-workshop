

## Load pre-cleaned and summarized data 
load(file='data/prepped_data.RData',verbose=TRUE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##             Prep bumble bee observations
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bumble.obs <- na.omit(bumble.obs)

## Use that object to make a 3D array (site x visit x species)
JAGS.arr  <- array(0,
                   dim = c(nsite    = nrow(site.info),
                           nvisit   = max(bumble.obs$Round),
                           nspecies = length(species.names)),
                   dimnames = list(site = site.info$siteID,
                                   visit = paste0('visit.',1:2),
                                   species = species.names[order(species.names)])
                   )

for(ii in 1:nrow(bumble.obs)){
  site <- bumble.obs$siteID[ii]
  visit <- bumble.obs$Round[ii]
  species <- bumble.obs$Species[ii]
  abund <- bumble.obs$Abundance[ii]

  JAGS.arr[site,visit,species] <- JAGS.arr[site, visit,species] + abund
  
}

## Turn all non-zero values to 1 (1 = site is occupied) 
JAGS.arr  <- (JAGS.arr > 0) * 1



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                    Prep site variables
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

JAGS.site <- site.info[,c('siteID','sitetype',
                          'canopyopenness','openflowerabundsite',
                          'flowerSRsite','latitude','region')]



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                   Prep visit variables
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
visit.vars <- c("TrapJulian", "HoursOut", "openflowerabundvisit", "flowerSRvisit")

## Turn it into a 3 dimensional array
JAGS.visit <- array(NA,
                    dim = c(nsite      = 26,
                            nvisit     = 2,
                            nvariables = 4),
                    dimnames = list(site      = site.info$siteID,
                                    visit     = c("visit.1","visit.2"),
                                    variables = visit.vars))


for(ii in 1:nrow(visit.info)){
  site  <- visit.info$siteID[ii]
  visit <- visit.info$Round[ii]

  JAGS.visit[site,visit,] <- unlist(visit.info[ii,visit.vars])
  
}

JAGS.obs <- JAGS.arr


save(bumble.obs,
     site.info,
     visit.info,
     species.names,
     JAGS.obs,
     JAGS.site,
     JAGS.visit,
     file='saved_objects/JAGS.data.RData')


