## Written by: Hanna Jackson hmj2@sfu.ca

## ~~~ Figure 6: Detection vs Covariates ~~~

make.det.fig <- function(){

  ## Set up the colours  
  JAGS.site$sitetypecol <- NA
  JAGS.site[which(JAGS.site$sitetype==0),"sitetypecol"] <- fire.col
  JAGS.site[which(JAGS.site$sitetype==1),"sitetypecol"] <- forest.col
  
  ## Make a function that will add a figure label to each of our pannels consistently! 
  fig.labs <-letters[1:6]
  add.fig.label <- function(fig.num,extra){
    text(x=par('usr')[2]-par('usr')[2]*(0.09)+extra, 
         y=par('usr')[4]*0.95, 
         fig.labs[fig.num],
         cex=1.6, 
         pos=4
         )
  }
  
  ## Make function that will add a point for each site-visit's model
  ##    estimated detection probability vs its x-value (for a chosen x
  ##    variable, `var`)
  ## This means there will be 2 times the number of sites = 2 x 26 =
  ##    52 points. This represents the fact that detection can change
  ##    from visit to visit based on the time of day, day of year  or
  ##    changes in the habitat like changing abundance of flowers. 
  
  add.points.site <- function(var){
    for (ss in 1:nrow(JAGS.site)){
      for(vv in 1:2){
        chains <- expit(sims.mat[,'p.0'] +
                        sims.mat[,'p.sitetype']             * my.data$sitetype            [ss]+
                        sims.mat[,'p.hoursout']             * my.data$hoursout            [ss,vv] +
                        sims.mat[,'p.canopyopenness']       * my.data$canopyopenness      [ss] +
                        sims.mat[,'p.openflowerabundvisit'] * my.data$openflowerabundvisit[ss,vv] +
                        sims.mat[,'p.flowerSRvisit']        * my.data$flowerSRvisit       [ss,vv]  +
                        sims.mat[,'p.trapjulian']           * my.data$trapjulian          [ss,vv]
                        )
        points(y = mean(chains),
               x = JAGS.site[ss,var], 
               pch = 16,
               col = JAGS.site$sitetypecol[ss])
        
        chains.bci <- quantile(chains, probs=c(0.025,0.975))
        arrows(x0 = JAGS.site[ss,var],
               y0 = chains.bci["2.5%"],
               x1 = JAGS.site[ss,var],
               y1 = chains.bci["97.5%"],
               lwd = 1, code = 0, angle = 90, length = 0.02, 
               col = makeTransparent(JAGS.site$sitetypecol[ss],0.5)
               )
      }
    }
  }
  
  ## Make function that will add a point for each site-visit's model
  ##    estimated detection probability vs its x-value (for a chosen x
  ##    variable, `var`)
  ## We needed this to be a separate function than the previous one
  ##    becasue we are now taking our x values from JAGS.visit (and
  ##    thus, need to index by an additional dimension), whereas in
  ##    the previous function we took them from JAGS.site.
  
  add.points.visit <- function(var){
    for (ss in 1:nrow(JAGS.site)){
      for(vv in 1:2){
        chains <- expit(sims.mat[,'p.0'] +
                        sims.mat[,'p.sitetype']             * my.data$sitetype            [ss]    +
                        sims.mat[,'p.hoursout']             * my.data$hoursout            [ss,vv] +
                        sims.mat[,'p.canopyopenness']       * my.data$canopyopenness      [ss]    +
                        sims.mat[,'p.openflowerabundvisit'] * my.data$openflowerabundvisit[ss,vv] +
                        sims.mat[,'p.flowerSRvisit']        * my.data$flowerSRvisit       [ss,vv] +
                        sims.mat[,'p.trapjulian']           * my.data$trapjulian          [ss,vv]
                        )
        ## Plot the points 
        points(y=mean(chains),
               x=JAGS.visit[ss,vv,var],
               pch=16, 
               col=JAGS.site$sitetypecol[ss])
        
        ## Make point BCIs
        chains.bci <- quantile(chains, probs=c(0.025,0.975))
        
        ## Then use that to plot the confidence intervals 
        arrows(x0 = JAGS.visit[ss,vv,var],
               y0 = chains.bci["2.5%"],
               x1 = JAGS.visit[ss,vv,var],
               y1 = chains.bci["97.5%"],
               lwd = 1, code = 0, angle = 90, length = 0.02, 
               col = makeTransparent(JAGS.site$sitetypecol[ss],0.5)
               )
      }
    }
  }

  
  ##  ~~~~~~~~~~~~ Actually plotting ~~~~~~~~~~~~

  ## Layout
  m <- matrix(1:6, ncol=3, nrow=2, byrow=TRUE)
  layout(m)

  ## Plotting options
  par(oma=c(0.1, 3, 1.5, 1), 
      mar=c(4.5, 0.5, 0.1, 0.1),
      mgp=c(2,0.4,0),
      tcl=0, 
      cex.axis=1.5,
      cex.lab =1.9,
      pty='s')

  
  ## ~~~~~ Pannel 1: Hoursout vs Detection ~~~~~

  ## Make a blank plot with the correct dimensions to start 
  plot(NA,
       xlim = range(my.data$hoursout),
       ylim = c(0,1),
       xlab = 'Total trap hours (# of traps * hours out)',
       ylab = NA, 
       xaxt = 'n',
       yaxt = 'n',
       bty  = "n",
       las  = 1)

  ## Make a function that calculates the mean and 95% BCI y values
  ##    (detection) for ANY x value (trap # of hours out) 
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'p.0'] + 
                    sims.mat[,'p.sitetype']             * mean(my.data$sitetype)+
                    sims.mat[,'p.hoursout']             * xx +
                    sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness) +
                    sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit) +
                    sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit)  +
                    sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }

  ## Now we need to tell it which x values we want it to calculate
  ##    those y values for 
  x.vals <- seq(from = min(my.data$hoursout), 
                to   = max(my.data$hoursout),
                length = 1000)

  ## And now we do the calculation! This will, for each x value in x.vals,
  ##    calculate the mean, and the upper and lower confidence limits of
  ##    detection probability 
  y.vals <- sapply(x.vals, get.y.val) 

  ## Now plot that on the graph 
  ## Plot mean line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=generic.col, lwd=4)
  ## Plot lower CI 
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=generic.col, lwd=2, lty=2)
  ## Plot upper CI
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=generic.col, lwd=2, lty=2)
  
  ## Final plotting arguments
  par(new=TRUE)
  plot(NA, 
       xlim = range(unscaled$hoursout), 
       ylim = c(0,1),
       xlab = '', 
       ylab = '',
       las  = 1)

  ## Adding the y label 
  mtext(text="Detection Probability",
        line=2.5,
        side=2,
        cex=1.4)

  ## Add a legend to this plot for what colour is burned and which
  ##    colour is unburned
  legend(x=66, 
         y=1.03,
         legend=c("Unburned Site","Burned Site"), 
         col=c(forest.col,fire.col), 
         pt.cex=2.5,
         cex = 1.5,
         pch=15, 
         bty="n",
         ncol = 1)

  ## Call the function we made earlier that will add each site as a
  ##    data point for the `hoursout` variable
  add.points.visit(var="HoursOut")

  ## Call the function that adds the figure labels 
  add.fig.label(fig.num=1, extra=0)
  


  
  ## ~~~~~ Pannel 2: Canopy vs Detection ~~~~~

  ## See pannel 1 code above for comments explaining how this code works
  
  plot(NA,
       xlim=range(my.data$canopyopenness),
       ylim=c(0,1),
       xlab='Canopy openness (% open)',
       ylab='', 
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'p.0'] + 
                    sims.mat[,'p.sitetype']             * mean(my.data$sitetype)+
                    sims.mat[,'p.hoursout']             * mean(my.data$hoursout) +
                    sims.mat[,'p.canopyopenness']       * xx +
                    sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit) +
                    sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit)  +
                    sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from = min(my.data$canopyopenness), 
                to   = max(my.data$canopyopenness),
                length = 1000)
  
  y.vals <- sapply(x.vals, get.y.val) 
  
  ## Plot mean line
  lines(x=x.vals, y=y.vals['mean',], 
        pch=16, col=generic.col, lwd=4)
  ## Plot lower CI 
  lines(x=x.vals, y=y.vals['2.5%',], 
        pch=16, col=generic.col, lwd=2, lty=2)
  ## Plot upper CI
  lines(x=x.vals, y=y.vals['97.5%',],
        pch=16, col=generic.col, lwd=2, lty=2)
  
  ## Final plotting arguments
  par(new=TRUE)
  plot(NA, 
       xlim=range(unscaled$canopyopenness), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)
  
  add.points.site(var='canopyopenness')
  add.fig.label(fig.num=2, extra=0)
  


  
  ## ~~~~~ Pannel 3: Openflowerabundvisit vs Detection ~~~~~

  ## See pannel 1 code above for comments explaining how this code works
  
  plot(NA,
       xlim=range(my.data$openflowerabundvisit),
       ylim=c(0,1),
       xlab='Log(Open flower abundance)',
       ylab='', 
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'p.0'] + 
                    sims.mat[,'p.sitetype']             * mean(my.data$sitetype)+
                    sims.mat[,'p.hoursout']             * mean(my.data$hoursout) +
                    sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness) +
                    sims.mat[,'p.openflowerabundvisit'] * xx +
                    sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit)  +
                    sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from=min(my.data$openflowerabundvisit), 
                to=max(my.data$openflowerabundvisit),
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
       xlim=range(unscaled$openflowerabundvisit), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)
  
  add.points.visit(var="openflowerabundvisit")
  add.fig.label(fig.num=3, extra=0)

  
  
  
  ## ~~~~~ Pannel 4: FlowerSRvisit vs Detection ~~~~~
 
  ## See pannel 1 code above for comments explaining how this code works
  
  plot(NA,
       xlim=range(my.data$flowerSRvisit),
       ylim=c(0,1),
       xlab='Flower species richness',
       ylab=NA, 
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  mtext(text="Detection Probability",
        line=2.5,
        side=2,
        cex=1.4)
  
  get.y.val <- function(xx) {
    chains <- expit(sims.mat[,'p.0'] + 
                    sims.mat[,'p.sitetype']             * mean(my.data$sitetype) +
                    sims.mat[,'p.hoursout']             * mean(my.data$hoursout) +
                    sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness) +
                    sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit) +
                    sims.mat[,'p.flowerSRvisit']        * xx   +
                    sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian)
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from=min(my.data$flowerSRvisit), 
                to=max(my.data$flowerSRvisit),
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
       xlim=range(unscaled$flowerSRvisit), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)
  
  add.points.visit(var="flowerSRvisit")
  add.fig.label(fig.num=4, extra=0)

  

  
  ## ~~~~~ Pannel 5: Julian Day vs Detection ~~~~~

  ## See pannel 1 code above for comments explaining how this code works

  plot(NA,
       xlim=range(my.data$trapjulian),
       ylim=c(0,1),
       xlab='Julian day',
       ylab='', 
       xaxt='n',
       yaxt='n',
       bty="n",
       las=1)
  
  get.y.val <- function(xx){
    chains <- expit(sims.mat[,'p.0'] + 
                    sims.mat[,'p.sitetype']             * mean(my.data$sitetype)+
                    sims.mat[,'p.hoursout']             * mean(my.data$hoursout) +
                    sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness) +
                    sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit) +
                    sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit)  +
                    sims.mat[,'p.trapjulian']           * xx
                    )
    c(mean=mean(chains), quantile(chains, probs=c(0.025,0.975)))
  }
  
  x.vals <- seq(from=min(my.data$trapjulian),  
                to=max(my.data$trapjulian),
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
       xlim=range(unscaled$trapjulian), 
       ylim=c(0,1),
       xlab='', 
       ylab='',
       las=1)
  
  add.points.visit(var="TrapJulian")
  add.fig.label(fig.num=5, extra=14)
  
  
  
  
  
  ## ~~~~~ Pannel 6: Sitetype vs Detection  ~~~~~

  ## Okay, here it gets different! In the previous five plots, we had
  ##    continuous variables that we wanted to plot on our x axis. That
  ##    meant that to make a smooth-looking line, we needed to calculate
  ##    detection (y value) for many many x values.

  ## SO! In this one, we no longer have to calculate those many many
  ##    diferent values of occuapcny for each possible x value, becasue
  ##    we only have two x values - burned and not burned represented by
  ##    0 and 1.

  ## So instead of making a function that can calculate y values for
  ##    any x value you might possibly want to plug in, we can just
  ##    manually make two objects, one with 0 (burned sites) and one with
  ##    1 (unburned sites). We'll call these chains.fire and chains.forest

  ## Index specifying which sites are burned
  burned <- my.data$sitetype == 0
  
  ## Extract chains for the overall trend for both fire and forest 
  chains.fire <- expit(sims.mat[,'p.0'] + 
                       sims.mat[,'p.sitetype']             * 0 +
                       sims.mat[,'p.hoursout']             * mean(my.data$hoursout[burned]) +
                       sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness[burned]) +
                       sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit[burned]) +
                       sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit[burned])  +
                       sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian[burned])
                       )
  
  chains.forest <- expit(sims.mat[,'p.0'] +
                         sims.mat[,'p.sitetype'] * 1 +
                         sims.mat[,'p.hoursout']             *  mean(my.data$hoursout[!burned]) +
                         sims.mat[,'p.canopyopenness']       * mean(my.data$canopyopenness[!burned]) +
                         sims.mat[,'p.openflowerabundvisit'] * mean(my.data$openflowerabundvisit[!burned]) +
                         sims.mat[,'p.flowerSRvisit']        * mean(my.data$flowerSRvisit[!burned]) +
                         sims.mat[,'p.trapjulian']           * mean(my.data$trapjulian[!burned]) )
  
  ## Now from those chains get the mean and sd 
  vals <- c(fire.mean   = mean(chains.fire), 
            fire        = quantile(chains.fire, probs=c(0.025,0.975)),
            forest.mean = mean(chains.forest),
            forest      = quantile(chains.forest, probs=c(0.025,0.975)))
  
  ## Set up the blank plot 
  plot(NA, 
       xlim = c(0.5,2.5),  
       ylim = c(0,1),
       ylab = NA, 
       xlab = "Burn status", 
       xaxt = 'n', 
       xgap.axis = 10, 
       las = 1)
  
  ## Add the overall point and arrow 
  arrows(x0=1:2,
         y0=c(vals[c("fire.2.5%","forest.2.5%")]),
         x1=1:2,
         y1=c(vals[c("fire.97.5%","forest.97.5%")]),
         lwd=7, code=0, angle=90, length=0.02, 
         col=c(fire.col, forest.col)
         )
  points(vals[c("fire.mean","forest.mean")], pch=16, cex=1.8, col='black')
  
  labs <- c("Burned","Unburned")
  axis(side=1, at=c(1,2), label=labs)

  ## Finally, we want to add in the estimated detection of each site-visit
  ##    based on that particular site-visit's characteristics.
  ##    So we're going to iterate through each site and visit and take the
  ##    variables associated with that site-visit and plug them into our
  ##    equation to estimate detection and then plot it on our pannel.
  
  for (ss in 1:nrow(JAGS.site)){
    for (vv in 1:2){
      chains <- expit(sims.mat[,'p.0'] + 
                      sims.mat[,'p.sitetype']             * my.data$sitetype[ss] +
                      sims.mat[,'p.hoursout']             * my.data$hoursout[ss,vv] +
                      sims.mat[,'p.canopyopenness']       * my.data$canopyopenness[ss] +
                      sims.mat[,'p.openflowerabundvisit'] * my.data$openflowerabundvisit[ss,vv] +
                      sims.mat[,'p.flowerSRvisit']        * my.data$flowerSRvisit[ss,vv]  +
                      sims.mat[,'p.trapjulian']           * my.data$trapjulian[ss,vv]
                      )
      
      jitter.val <- 1-jitter(1,factor=19)
      jittered   <- JAGS.site[ss,"sitetype"]+1+jitter.val
      chains.bci <- quantile(chains, probs=c(0.025,0.975))
      arrows(x0=jittered,
             y0=chains.bci["2.5%"],
             x1=jittered,
             y1=chains.bci["97.5%"],
             lwd=1, code=0, angle=90, length=0.02, 
             col=makeTransparent(JAGS.site$sitetypecol[ss],0.5)
             )
      points(y=mean(chains),
             x=jittered,
             pch=16,
             col=JAGS.site$darksitetypecol[ss])
    }
  }
  
  add.fig.label(6,extra=0.05)
}
