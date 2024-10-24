## Written by: Hanna Jackson hmj2@sfu.ca

## ~~~ Figure 4: Effect Sizes (and/or Figure A3-A8) ~~~
make.effect.size.fig <- function(){

  ## Setting up figure labels function
  fig.labs <- c("a",'b')
  add.fig.label <- function(fig.num){
    xmin <- par("usr")[1]
    xmax <- par("usr")[2]
    ymin <- par("usr")[3]
    ymax <- par("usr")[4]
    text(x = xmin + (xmax - xmin) * 0.9,
         y = ymin + (ymax - ymin) * 0.95, 
         fig.labs[fig.num],
         cex = 1.1,
         pos = 4
         )
  }

  ## Make the multipannel layout
  m <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
  layout(m)

  ## Set plotting defaults 
  par(oma=c(3, 2, 1.5, 0.5), 
      mar=c(6.2, 0.5, 1.5, 0.1),
      mgp=c(2,0.2,0),
      tcl=0, 
      cex.axis=1.2,
      pty='s')
  

  
  ## ~~~~~~~~~~~~~~~ Pannel 1: Occupancy ~~~~~~~~~~~~~~~~

  ## ~~~ Solving a frustration of mine ~~~
  
  ## This process was born out of the frustration that came with
  ##   needed to make a new R file that was slightly modified for each
  ##   new model I ran with slightly different combinations of parameters.
  
  ## We dont want to plot all the parameters of the model, as there's
  ##   many that we dont care about, such as the intercept psi.0 or
  ##   the 13 (one for each species) psi.sp parameters that represent
  ##   the random effect of species on occupancy. We'll need a way
  ##   of filtering those out too. 
  
  ## So to solve my frustrations, I wrote this code that would check
  ##    which parameters that we care about visializing were in the
  ##    loaded model run and only plot those. 

  
  ## To start, here is a vector of all the possible occupancy parameters
  ##   that we might want to plot, where each element is named with
  ##   the full name as we would want that variable to be labeled on
  ##   the plot. 
  params.of.interest <- c("Burn status"             = 'psi.sitetype',
                          "Canopy openness"         = 'psi.canopyopenness',
                          'Latitude'                = 'psi.latitude',
                          'Region'                  = 'psi.region' ,
                          "Open flower abundance"   = 'psi.openflowerabundsite',
                          "Flower species richness" = 'psi.flowerSRsite', 
                          "Burn status x canopy openness"       = 'psi.sitetypexcanopyopenness',
                          "Burn status x open flower abundance" = 'psi.sitetypexopenflowerabundsite')

  ## We then grab the names of all parameters in the model run that
  ## contain 'psi.', which indicates that theyre parameters of
  ## occupancy 
  params.psi <- my.params[grepl("psi.", my.params,fixed = TRUE)]

  ## We then ask, which of those are in our params.of.interest vector
  ## we just made? 
  index <- params.of.interest %in% params.psi
  cols  <- params.of.interest[index]

  ## Get the names of those so that we can use it later to label out plot
  labs  <- names(cols)

  ## Get the mean and bayesian credible intervals that we will plot! 
  mu <- summ[cols,'mean']
  lb <- summ[cols,'2.5%']
  ub <- summ[cols,'97.5%']


  
  ## ~~~ Slight, necessary detour ~~~ 
  
  ## In the model run, 0 is burned and 1 is unburned but in the
  ## manuscript we plot it the opposite way around for clairity, so
  ## for plotting we multiply it by negative one 
  if("psi.sitetype" %in% my.params){
    mu['psi.sitetype'] <- mu['psi.sitetype'] * -1
    lb['psi.sitetype'] <- lb['psi.sitetype'] * -1
    ub['psi.sitetype'] <- ub['psi.sitetype'] * -1
  }
  
  if("psi.sitetypexcanopyopenness" %in% my.params){
    mu['psi.sitetypexcanopyopenness'] <- mu['psi.sitetypexcanopyopenness'] * -1
    lb['psi.sitetypexcanopyopenness'] <- lb['psi.sitetypexcanopyopenness'] * -1
    ub['psi.sitetypexcanopyopenness'] <- ub['psi.sitetypexcanopyopenness'] * -1
  }
  
  if("psi.sitetypexopenflowerabundsite" %in% my.params){
    mu['psi.sitetypexopenflowerabundsite'] <- mu['psi.sitetypexopenflowerabundsite'] * -1
    lb['psi.sitetypexopenflowerabundsite'] <- lb['psi.sitetypexopenflowerabundsite'] * -1
    ub['psi.sitetypexopenflowerabundsite'] <- ub['psi.sitetypexopenflowerabundsite'] * -1
  }

  ## Change the plotting options if your model run includes this
  ## parameter with a particularly long name...
  if("psi.sitetypexopenflowerabundsite" %in% my.params){
    par(mar=c(8.7, 5, 2, 2) + 0.1)
  }

  ## If the lower and upper bound of the confidence limits are both
  ## positie or both negative, that means that the confidence interval
  ## doesn't overlap zero and we want to highlight it in our
  ## `fire.col` which is a nice golden yellow.

  ## To figure this out, we just multiply the lower and upper bound
  ## together, and if the result is positive, then either both ub and
  ## lb were positive or both negative! 
  col.vals <- ub * lb
  colours <- ifelse(col.vals>0,fire.col,"darkgrey")

  
  ## ~~~ Okay, now to actually plot it! ~~~

  ## Set up the plotting region with a blank plot 
  plot(NA, 
       xlim=c(1,length(mu)),
       ylim=range(c(lb,ub+10)), 
       ylab=NA, 
       xlab=NA, 
       xaxt='n', 
       xgap.axis=10, 
       las=1,
       main="Occupancy") 

  ## Add all credible intervals from `lb` and `ub` 
  arrows(x0=1:length(mu),
         y0=lb,
         x1=1:length(mu),
         y1=ub,
         lwd=3, code=0, angle=90, length=0.02, 
         col=colours)

  ## Add points from `mu` on top of that 
  points(mu, pch=16, cex=1.5, col='black')

  ## Add a line at zero 
  abline(h=0, col = "darkgrey", lty = 5)
  
  ## Add the variable names using the `labs` object we made
  text(x=seq_along(labs)+0.1,
       y=par('usr')[3],
       labels=labs,
       srt=35,
       adj=c(1.1,1.1),
       xpd=NA,
       cex=1.2)

  ## Y axis label 
  mtext(text="Effect size",
        line=2,
        side=2,
        cex=1.4)

  ## Figure label via the function we wrote earlier 
  add.fig.label(1)



  
  
  ## ~~~~~~~~~~~~~~~ Pannel 2: Detection ~~~~~~~~~~~~~~~~
  
  ## Repeat all that but for detection probability, see above code for
  ## comments explaining how it works  
  
  params.of.interest <- c("Burn status"             = 'p.sitetype',
                          "Canopy openness"         = 'p.canopyopenness',
                          "Open flower abundance"   = 'p.openflowerabundvisit',
                          "Flower species richness" = 'p.flowerSRvisit', 
                          "Trap hours out"          = 'p.hoursout',
                          "Julian day"              = 'p.trapjulian')

  params.p <- my.params[grepl("p.", my.params,fixed = TRUE)]
  
  index <- params.of.interest %in% params.p
  cols  <- params.of.interest[index]
  labs  <- names(cols)

  mu <- summ[cols,'mean']
  lb <- summ[cols,'2.5%']
  ub <- summ[cols,'97.5%']

  if('p.sitetype' %in% cols){
    mu['p.sitetype'] <- mu['p.sitetype'] * -1
    lb['p.sitetype'] <- lb['p.sitetype'] * -1
    ub['p.sitetype'] <- ub['p.sitetype'] * -1
  }


  ## Determine which parameters the BCIs do not overlap zero and
  ## assign them the highlight colour (yellow)
  col.vals <- ub*lb
  colours <- ifelse(col.vals>0,fire.col,"darkgrey")

  
  ## ~~~ Actually plot it ~~~

  if(length(mu)<1){
    plot(NA, 
         xlim=c(1,2),
         ylim=range(0,1), 
         ylab='', 
         xlab=NA, 
         xaxt='n', 
         xgap.axis=10, 
         las=1,
         main="Detection")  
  } else {
    
    plot(NA, 
         xlim=c(1,length(mu)),
         ylim=range(c(lb,ub)), 
         ylab='', 
         xlab=NA, 
         xaxt='n', 
         xgap.axis=10, 
         las=1,
         main="Detection") 
    arrows(x0=1:length(mu),
           y0=lb,
           x1=1:length(mu),
           y1=ub,
           lwd=3, code=0, angle=90, length=0.02, 
           col=colours)
    points(mu, pch=16, cex=1.5, col='black')

    abline(h=0, col = "darkgrey", lty = 5)
    text(x=seq_along(labs)+0.1,
         y=par('usr')[3],
         labels=labs,
         srt=35,
         adj=c(1.1,1.1),
         xpd=NA,
         cex=1.2)
  }
  add.fig.label(2)
  
}



