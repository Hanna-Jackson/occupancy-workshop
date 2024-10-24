## Written by Hanna Jackson hmj2@sfu.ca

my.model <- function() {
  
  ## ~~~~~~~~~~~~ PRIORS ~~~~~~~~~~~~~~ 
  ## Here in the priors we define the parameters that we're going to
  ##    use, and the range of values that we exoect those parameter
  ##    values might be
  
  psi.0                       ~ dnorm(0,0.001) ## Every species shares this intercept
  psi.sitetype                ~ dnorm(0,0.001) ## The effect of a site being burned on occupancy
  psi.canopyopenness          ~ dnorm(0,0.001) ## The effect of canopyopenness on occupancy
  psi.openflowerabundsite     ~ dnorm(0,0.001) ## The effect of floral abundance on occupancy
  psi.flowerSRsite            ~ dnorm(0,0.001) ## The effect of floral specied richness on occupancy
  psi.sitetypexcanopyopenness ~ dnorm(0,0.001) ## An interaction between site type and canopyopenness on occupancy.
                                               ##    This means that burned and unburned sites are free to have
                                               ##    different effects of canopyopenness on occupancy
  
  p.0                    ~ dnorm(0,0.001) ## Every species shares some common intercept of detection probability
  p.sitetype             ~ dnorm(0,0.001) ## The effect of burn status on detection probability
  p.canopyopenness       ~ dnorm(0,0.001) ## The effect of canopyopenness on detection probability
  p.openflowerabundvisit ~ dnorm(0,0.001) ## The effect of floral abundance on detection probability
  p.flowerSRvisit        ~ dnorm(0,0.001) ## The effect of floral species richness on detection probability
  p.hoursout             ~ dnorm(0,0.001) ## The effect of trap effort on detection probability
  p.trapjulian           ~ dnorm(0,0.001) ## The effect of day of the year on detection probability

  

  ## ~~~ Random effect of species on occupancy ~~~
  
  ## This allows different species to have a different "baseline"
  ## occupancy that is due to factors beyond the things we measured in
  ## this study. You can think of it as "the effect of being x
  ## species"
  
  ## This is the standard deviation calculation: 
  sigma.psi.sp   ~ dunif(0,10) 
  tau.psi.sp    <- 1/(sigma.psi.sp*sigma.psi.sp)
  
  ## And this is where we actually draw the values. For each species,
  ## it's unique intercept is a draw from a common distribution with
  ## some estimated standard deviation and a mean of 0. 
  for(species in 1:nspecies) {
    psi.sp[species] ~ dnorm(0, tau.psi.sp) 
  } 
  
  
  ## ~~~~~~~~~~~~~~~ MODEL ~~~~~~~~~~~~~~~

  ## For each site and each species, we get one value of occupancy
  ##    probability that is calculated with the following formula: 
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      logit(psi[site, species]) <- psi.0                        + ## Occupancy intercept + 
        psi.sp[species]                                         + ## Species-specific intercept +
        psi.sitetype                * sitetype           [site] + ## Effect of burn status +
        psi.canopyopenness          * canopyopenness     [site] + ## Effect of canopy openness +
        psi.openflowerabundsite     * openflowerabundsite[site] + ## Effect of floral abund +
        psi.flowerSRsite            * flowerSRsite       [site] + ## Effect of floral spp richness + 
        psi.sitetypexcanopyopenness * sitetype[site]*canopyopenness[site] ## Interaction between sitetype and canopyopenness

      ## For each site, and each species, and each visit, we get one
      ##    value of detection probability that is calculated with the
      ##    following formula: 
      for(visit in 1:nvisit){ 
        logit(p[site,visit, species]) <- p.0                       + ## Detection intercept
         p.sitetype             * sitetype            [site]       + ## Effect of burn status
         p.hoursout             * hoursout            [site,visit] + ## Effect of trap hours out (effort) 
         p.trapjulian           * trapjulian          [site,visit] + ## Effect of day of the year
         p.openflowerabundvisit * openflowerabundvisit[site,visit] + ## Effect of floral abundance
         p.flowerSRvisit        * flowerSRvisit       [site,visit] + ## Effect of floral species richness
         p.canopyopenness       * canopyopenness      [site]         ## Effect of canopy openness
      } 
    } 
  }  

  
  ## ~~~~~~~~~~~~~~ LIKELIHOOD ~~~~~~~~~~~~~~~

  ## For each site and species...
  for(site in 1:nsite) {
    for(species in 1:nspecies){
      ## Occurrence
      Z[site, species] ~ dbern(psi[site, species])  ## The TRUE occupancy status of the site is a weighted coin flip based on
                                                    ##    the occupancy probability of the site

      ## For each site and species and visit... 
      for(visit in 1:nvisit) {
        p.eff[site, visit, species] <- Z[site, species] * p[site,visit,species] ## Effective detection is the detection
                                                                                ##    probability times whether or not the site
                                                                                ##    is occupied (1 for occupied 0 if not)
        X[site,visit,species] ~ dbern(p.eff[site,visit, species]) ## And whether or not we see the species (our data we
                                                                  ##    collected, X) is a weighted coin flip based on
                                                                  ##    that effective detection probability. 
      }
    }
  }  
} ## End of model 


## Specify the parameters that JAGS will store the values of
my.params <- c('p.0',
               'p.sitetype',
               'p.hoursout', 
               'p.trapjulian',
               'p.canopyopenness',
               "p.openflowerabundvisit", 
               "p.flowerSRvisit",
               'psi.0', 
               'psi.sp',
               'sigma.psi.sp',
               'psi.sitetype',
               'psi.canopyopenness',
               "psi.openflowerabundsite",
               "psi.flowerSRsite",
               'psi.sitetypexcanopyopenness'
)


