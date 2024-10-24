## Written by: Hanna Jackson hmj2@sfu.ca

## ~~~ Figure 5: Occupancy vs covariates ~~~~~

make.occ.figure <- function(){

  ## ~~~~~ Setting up ~~~~~ 
  col <- rep(generic.col, 3) # same colour repeated 3 times

  ## Make a function that will add a figure label to each of our pannels consistently!  
  fig.labs <- c("a","b","c","d")
  add.fig.label <- function(fig.num,extra){
    text(x=par('usr')[2]-par('usr')[2]*(0.09)+extra, 
         y=par('usr')[4]*0.95, 
         fig.labs[fig.num],
         cex=1.5, 
         pos=4
         )
  }
  
  ## Make function that will add a point for each site's model
  ## estimated occupancy vs its x-value (for a chosen x variable, `var`)
  add.points <- function(var){
    for (ss in 1:nrow(JAGS.site)){
      chains <- expit(sims.mat[,'psi.0'] +
                      sims.mat[,'psi.sitetype']                * my.data$sitetype[ss]            +
                      sims.mat[,'psi.canopyopenness']          * my.data$canopyopenness[ss]      +
                      sims.mat[,'psi.openflowerabundsite']     * my.data$openflowerabundsite[ss] +
                      sims.mat[,'psi.flowerSRsite']            * my.data$flowerSRsite[ss]        +
                      sims.mat[,'psi.sitetypexcanopyopenness'] * my.data$canopyopenness[ss] * my.data$sitetype[ss]
                      )
      points(y = mean(chains),
             x = JAGS.site[ss,var],
             pch = 16, 
             col = JAGS.site$sitetypecol[ss]
             )
      chains.bci <- quantile(chains, probs=c(0.025,0.975))
      arrows(x0 = JAGS.site[ss,var],
             y0 = chains.bci["2.5%"],
             x1 = JAGS.site[ss,var],
             y1 = chains.bci["97.5%"],
             lwd=1, code=0, angle=90, length=0.02, 
             col=makeTransparent(JAGS.site$sitetypecol[ss],0.5)
             ) 
    }
  }

  ## Layout
  m <- matrix(c(1,2,3,4), ncol=2, nrow=2,byrow=TRUE)
  layout(m) 

  ## Setting plotting options 
  par(oma=c(0.1,4,0.1,0.1),
      mar=c(5, 1.5, 1, 1) + 0.1,
      mgp=c(2,0.2,0),
      tcl=0,
      cex.lab = 1.7,
      cex.axis = 1.4) 
  
  
  ## ~~~~~ Pannel 1: Canopy Openness vs Occupancy ~~~~~

  ## Make a blank plot with the correct dimensions to start 
  plot(NA,
       xlim=range(my.data$canopyopenness),
       ylim=c(0,1),
       xlab='Canopy openness (% open)',
       ylab=NA,
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)

  ## Because we have an interaction term in our model, we are going to
  ## plot 2 lines, one for burned and one for unburned sites. This
  ## doesnt change how we do this, other than that we do all our
  ## calculations on half the data for each line. 
  
  ## ~~~ 1st: FIRE line ~~~

  ## Get the indicies that represent burned sites
  burned <- my.data$sitetype == 0

  ## Make a function that calculates the mean and 95% BCI y values
  ##    (occupancy) for ANY x value (canopyopenness) 
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'psi.0'] +
                    sims.mat[,'psi.sitetype']                * 0 +
                    sims.mat[,'psi.canopyopenness']          * xx +
                    sims.mat[,'psi.openflowerabundsite']     * mean(my.data$openflowerabundsite[burned]) +
                    sims.mat[,'psi.flowerSRsite']            * mean(my.data$flowerSRsite[burned])  +
                    sims.mat[,'psi.sitetypexcanopyopenness'] * xx * 0
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }

  ## Now we need to tell it which x values we want it to calculate
  ##    those y values for 
  x.vals <- seq(from = min(my.data$canopyopenness), ## just tells it to plot x for a range of y's 
                to   = max(my.data$canopyopenness),
                length = 1000)

  ## And now we do the calculation! This will, for each x value,
  ## calculate the mean and the upper and lower confidence limits of
  ## occupancy 
  y.vals <- sapply(x.vals, get.y.val) ## apply get y values to each value in x.vals

  ## Now plot that on the graph 
  ## Mean Line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=fire.col, lwd=4)
  ## Lower CI Line
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=fire.col, lwd=2, lty=2)
  ## Upper CI Line
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=fire.col, lwd=2, lty=2)
  
  ## ~~~ 2nd: FOREST line ~~~
  ## Here we repeat that whole process, but with the forest subset of
  ## sites this time! 
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'psi.0'] +
                    sims.mat[,'psi.sitetype']                * 1 +
                    sims.mat[,'psi.canopyopenness']          * xx +
                    sims.mat[,'psi.openflowerabundsite']     * mean(my.data$openflowerabundsite[!burned]) +
                    sims.mat[,'psi.flowerSRsite']            * mean(my.data$flowerSRsite[!burned])  +
                    sims.mat[,'psi.sitetypexcanopyopenness'] * xx * 1
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from = min(my.data$canopyopenness), #just tells it to plot x for a range of y's 
                to = max(my.data$canopyopenness),
                length = 1000)
  
  y.vals <- sapply(x.vals, get.y.val) 
  
  ## Mean line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=forest.col, lwd=4)
  ## Lower CI line
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=forest.col, lwd=2, lty=2)
  ## Upper CI line 
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=forest.col, lwd=2, lty=2)
  
  ## Final plot arguments 
  par(new=TRUE)
  plot(NA, 
       xlim=range(unscaled$canopyopenness), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)

  ## Call the function that adds the figure labels 
  add.fig.label(1,extra=0)

  ## Call the function we made earlier that will add each site as a
  ## data point for the `canopyopenness` variable
  add.points(var="canopyopenness")

  
  
  ## ~~~~~ Pannel 2: Openflowerabundsite  vs Occupancy ~~~~~

  ## See comments for the first pannel to explain how this works!
  
  plot(NA,
       xlim=range(my.data$openflowerabundsite),
       ylim=c(0,1),
       xlab='Log(Open flower abundance)',
       ylab='',
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'psi.0']                                                                +
                    sims.mat[,'psi.sitetype']                * mean(my.data$sitetype)                 +
                    sims.mat[,'psi.canopyopenness']          * mean(my.data$canopyopenness)           +
                    sims.mat[,'psi.openflowerabundsite']     * xx                                     +
                    sims.mat[,'psi.flowerSRsite']            * mean(my.data$flowerSRsite)             +
                    sims.mat[,'psi.sitetypexcanopyopenness'] * mean(my.data$canopyopenness) * mean(my.data$sitetype)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  x.vals <- seq(from=min(my.data$openflowerabundsite),
                to=  max(my.data$openflowerabundsite),
                length=1000)
  
  y.vals <- sapply(x.vals, get.y.val) 
  
  ## Mean line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=generic.col, lwd=4)
  ## Lower CI line
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=generic.col, lwd=2, lty=2)
  ## Upper CI line 
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=generic.col, lwd=2, lty=2)
  
  ## Final plotting arguments
  par(new=TRUE)
  plot(NA, 
       xlim=range(unscaled$openflowerabundsite), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)
  
  add.fig.label(2,extra=0)
  
  legend(x=1.2,
         y=0.98,
         legend=c("Unburned Site","Burned Site"), 
         col=c(forest.col,fire.col),
         cex = 1.2,
         pt.cex=2.3, 
         pch=15, 
         bty="n",
         ncol = 1)
  
  add.points(var="openflowerabundsite")
  

  
  ## ~~~~~ Pannel 3: Flower SR vs Occupancy ~~~~~

  ## See comments for the first pannel to explain how this works!
  plot(NA,
       xlim=range(my.data$flowerSRsite),
       ylim=c(0,1),
       xlab='Flower species richness',
       ylab=NA, 
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'psi.0']                                                           +
                    sims.mat[,'psi.sitetype']                * mean(my.data$sitetype)            +
                    sims.mat[,'psi.canopyopenness']          * mean(my.data$canopyopenness)      +
                    sims.mat[,'psi.openflowerabundsite']     * mean(my.data$openflowerabundsite) +
                    sims.mat[,'psi.flowerSRsite']            * xx                                +
                    sims.mat[,'psi.sitetypexcanopyopenness'] * mean(my.data$canopyopenness) * mean(my.data$sitetype)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from = min(my.data$flowerSRsite), 
                to   = max(my.data$flowerSRsite),
                length = 1000)
  
  y.vals <- sapply(x.vals, get.y.val) 
  
  ## Mean line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=generic.col, lwd=4)
  ## Lower CI Line
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=generic.col, lwd=2, lty=2)
  ## Upper CI Line
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=generic.col, lwd=2, lty=2)
  
  ## Final plotting arguments
  par(new=TRUE)
  plot(NA, 
       xlim = range(unscaled$flowerSRsite), 
       ylim = c(0,1),
       xlab = '', 
       ylab = '',
       las  = 1)
  
  add.fig.label(3,extra=0)

  add.points(var="flowerSRsite")


  
  ## ~~~ Pannel 4: Sitetype vs Occupancy Boxplot ~~~ 

  ## Okay, here it gets different! In the previous three plots, we had
  ##    continuous variables that we wanted to plot on our x axis. That
  ##    meant that to make a smooth-looking line, we needed to calculate
  ##    occupancy (y value) for many many x values.

  ## SO! In this one, we no longer have to calculate those many many
  ##    diferent values of occuapcny for each possible x value, becasue
  ##    we only have two x values - burned and not burned represented by
  ##    0 and 1.

  ## So instead of making a function that can calculate y values for
  ##    any x value you might possibly want to plug in, we can just
  ##    manually make two objects, one with 0 (burned sites) and one with
  ##    1 (unburned sites).
  burned <- my.data$sitetype==0
  chains.fire <- expit(sims.mat[,'psi.0'] +
                       sims.mat[,'psi.sitetype']                * 0 +
                       sims.mat[,'psi.canopyopenness']          * mean(my.data$canopyopenness       [burned]) +
                       sims.mat[,'psi.openflowerabundsite']     * mean(my.data$openflowerabundsite  [burned]) +
                       sims.mat[,'psi.flowerSRsite']            * mean(my.data$flowerSRsite         [burned]) +
                       sims.mat[,'psi.sitetypexcanopyopenness'] * mean(my.data$canopyopenness       [burned]) * 0  
                       )
  
  chains.forest <- expit(sims.mat[,'psi.0'] +
                         sims.mat[,'psi.sitetype']                * 1 +
                         sims.mat[,'psi.canopyopenness']          * mean(my.data$canopyopenness     [!burned]) +
                         sims.mat[,'psi.openflowerabundsite']     * mean(my.data$openflowerabundsite[!burned]) +
                         sims.mat[,'psi.flowerSRsite']            * mean(my.data$flowerSRsite       [!burned]) +
                         sims.mat[,'psi.sitetypexcanopyopenness'] * mean(my.data$canopyopenness     [!burned]) * 1 
                         )    

  ## Next we'll calculate the interesting quantities from those vectors 
  vals <- c(fire.mean   = mean(chains.fire),
            forest.mean = mean(chains.forest),
            fire        = quantile(chains.fire,   probs=c(0.025,0.975)),
            forest      = quantile(chains.forest, probs=c(0.025,0.975)))

  ## Blank plot 
  plot(NA, 
       xlim=c(0.5,2.5), 
       ylim=c(0,1),
       ylab=NA, 
       xlab="Burn status", 
       xaxt='n', 
       xgap.axis=10, 
       las=1) 

  ## Add the overall mean trends as CIs
  arrows(x0 = 1:2,
         y0 = c(vals[c("fire.2.5%","forest.2.5%")]),
         x1 = 1:2,
         y1 = c(vals[c("fire.97.5%","forest.97.5%")]),
         lwd=7, code=0, angle=90, length=0.02, 
         col=c(fire.col, forest.col)
         )
  
  ## Add the overall mean point 
  points(vals[c("fire.mean","forest.mean")], pch=16, cex=1.8, col='black')

  ## Add labels 
  labs <- c("Burned","Unburned")
  axis(side=1, at=c(1,2), label=labs)

  ## Finally, we want to add in the estimated occupancy of each site
  ##    based on that particular site's characteristics.
  ##    So we're going to iterate through each site, and take the
  ##    variables associated with that site and plug them into our
  ##    equation to estimate occupancy and then plot it on our pannel. 
  for (ss in 1:nrow(JAGS.site)){
    chains <- expit(sims.mat[,'psi.0'] +
                    sims.mat[,'psi.sitetype']                * my.data$sitetype[ss] +
                    sims.mat[,'psi.canopyopenness']          * my.data$canopyopenness[ss] +
                    sims.mat[,'psi.openflowerabundsite']     * my.data$openflowerabundsite[ss] +
                    sims.mat[,'psi.flowerSRsite']            * my.data$flowerSRsite[ss]  +
                    sims.mat[,'psi.sitetypexcanopyopenness'] * my.data$canopyopenness[ss] * my.data$sitetype[ss]
                    )
    makeTransparent(JAGS.site$sitetypecol[ss],0.5)

    ## We'll jitter the points side to side a bit so we can see each
    ##    individual site better. 
    jitt     <- jitter(0, factor=18)
    jittered <- (JAGS.site[ss,"sitetype"]+1) + jitt 

    ## Plot the CIs  
    arrows(x0 = jittered,
           y0 = quantile(chains,probs=c(0.025,0.975))[1],
           x1 = jittered,
           y1 = quantile(chains,probs=c(0.025,0.975))[2],
           lwd=1, code=0, angle=90, length=0.02, 
           col=makeTransparent(JAGS.site$sitetypecol[ss],0.5)
           )
    ## Plot the points 
    points(y = mean(chains),
           x = jittered,
           pch = 16,
           col = JAGS.site$sitetypecol[ss])

  }
  
  add.fig.label(4,extra=0.05)
  
  mtext("Occupancy Probability",
        at = 0.29, 
        line = 1,
        outer = T,
        side = 2,
        cex = 1.4)
  
  mtext("Occupancy Probability",
        at = 0.79, 
        line = 1,
        outer = T,
        side  = 2,
        cex = 1.4)
  
}
