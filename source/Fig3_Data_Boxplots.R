## Written by: Hanna Jackson hmj2@sfu.ca

## ~~~ Figure 3: Plot the distribution of the data with boxplots ~~~ 
plot.data.boxplots <- function(){

  ## Setup the multipannel layout 
  m <- matrix(c(1,2,3,4), ncol=2, nrow=2, byrow=TRUE)
  layout(m) 

  ## Set up the plotting options 
  par(oma = c(0.1, 2, 1.5, 1), 
      mar = c(4.2, 0.5, 0.1, 0.1),
      mgp = c(2,0.2,0),
      tcl = 0, 
      cex.axis = 1.3, 
      pty = 's')
  
  fig.labs <- c("a",'b','c','d')

  
  ## ~~~~~ Functions ~~~~~

  ## Function that will add figure labels to each corner consistently 
  add.fig.label <- function(fig.num,extra){
    text(x = par('usr')[1]+(par('usr')[2]-par('usr')[1])*(fig.lab.x)+extra, 
         y = par('usr')[4]*0.95, 
         fig.labs[fig.num],
         cex = 1.5, 
         pos = 4
         )
  }
  
  fig.lab.x <- 0.89

  ## The function that will plot the boxplots that you put into it 
  make.boxplot <- function(forest.data, fire.data, xlabel, fig.num){
    
    boxplot(forest.data, fire.data,
            horizontal=TRUE,
            xlab = xlabel,
            cex.lab = 1.4,
            ylab = "",
            col = c(forest.col, fire.col),
            outline = FALSE
            )
    stripchart(forest.data, add=TRUE, at = 1,
               method = "jitter", jitter = 0.045,
               pch = 22, cex = 1.05,
               bg = dark.forest.col.2
               )
    stripchart(fire.data,   add=TRUE, at = 2,
               method = "jitter", jitter = 0.045,
               pch = 22, cex = 1.05,
               bg =dark.fire.col
               )
    
    ## Only add labels on the y axis if it's a figure on the left
    ##    (pannels 1 and 3) 
    if(fig.num == 1 | fig.num == 3){
      mtext(text="Burned",
            line=1,
            side=2,
            cex = 1.4,
            at=2)
      mtext(text="Unburned",
            line=1,
            side=2,
            cex= 1.4,
            at=1)
    }

    ## Call the function that adds figure labels 
    add.fig.label(fig.num, extra=0)
    
  }

  ## Now we'll use those functions we just made to actually make each
  ##    of the pannels of our figure! 
  
  ## ~~~~~~ PANNEL 1: Canopyopenness ~~~~~~
  
  make.boxplot(forest.data = forest$canopyopenness,
               fire.data   = fire  $canopyopenness,
               xlabel = "Canopy openness (% open)",
               fig.num = 1)
  

  ## ~~~~~~~~ PANNEL 2: Open flower abundance site ~~~~~~~~~~

  make.boxplot(forest.data = forest$openflowerabundsite,
               fire.data   = fire  $openflowerabundsite,
               xlabel = "Log(Open flower abundance)",
               fig.num = 2)
  

  ## ~~~~~~~~~~ PANNEL 3: Flower species richness site ~~~~~~~~~~~
  
  make.boxplot(forest.data = forest$flowerSRsite,
               fire.data   = fire  $flowerSRsite,
               xlabel = "Flower species richness",
               fig.num = 3)
  
  
  ## ~~~~~~~~~~~ PANNEL 4: Hours out ~~~~~~~~~~~~
  make.boxplot(forest.data = hoursout.forest,
               fire.data   = hoursout.fire,
               xlabel = "Total trap hours (# of traps * hrs out)",
               fig.num = 4)
  
}




