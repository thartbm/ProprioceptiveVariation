

plotParadigm <- function(target='inline', xticks='trial', ylab='rotation [°]') {
  
  schedule <- getSchedule()
  schedule <- schedule[which(schedule$type == 'reach'),]
  # add an extra block for easier plotting:
  # schedule <- rbind(schedule, data.frame(trial=961, pair=481, block=31, type='reach', rotation=1,clamp=1))
  
  plot(-1000,-1000, xlim=c(0,481), ylim=c(-30,30),
       main='',xlab='',ylab=ylab,
       bty='n',ax=F)
  
  title(xlab=xticks)
  
  blockstarts <- c(1)
  blockmids   <- c() # do we put block 0 on trial 25?

  for (pair in c(2:480)) {
    
    if (any(schedule[pair,c('rotation','clamp')] != schedule[pair-1,c('rotation','clamp')])) {
      # draw line from previous blockstart to current pair
      # indicating rotation
      # and normal or clamped trials
      lty <- c(1,3)[schedule$clamp[pair-1]+1]
      x <- c(blockstarts[length(blockstarts)],pair)
      y <- rep(schedule$rotation[pair-1],2)
      
      lines(x=x,y=y,lty=lty)
      
      # note block start for x axis ticks
      blockstarts <- c(blockstarts,pair)
      if (schedule$block[pair-1] > 0) {
        blockmids  <- c(blockmids, mean(schedule$pair[which(schedule$block == schedule$block[pair-1])]))
      }

      # draw lines connecting the rotation stuff
      lty <- 1
      x <- rep(pair, 2)
      y <- schedule$rotation[c(pair-1, pair)]
      lines(x=x, y=y, lty=lty)
      
    }
    
  }
  
  lty <- c(1,3)[schedule$clamp[480]+1]
  x <- c(blockstarts[length(blockstarts)],480)
  y <- rep(schedule$rotation[480],2)
  lines(x=x,y=y,lty=lty)
  blockstarts <- c(blockstarts, 480)
  blockmids  <- c(blockmids, mean(schedule$pair[which(schedule$block == schedule$block[pair])]))
  
  if (xticks == 'trial') {
    axis(side=1,at=blockstarts,cex.axis=0.8,las=2)
  }
  if (xticks == 'block') {
    axis(side=1,at=blockmids,labels=sprintf('%d',c(1:length(blockmids))),cex.axis=0.8,las=2)
  }
  axis(side=2,at=seq(-30,30,15),las=2)
  
  
  # if (target %in% c('svg','pdf')) {
  #   dev.off()
  #}
  
}

plotData <- function() {
  
  plotParadigm(xticks='block', ylab='angular deviation [°]')
  
  reaches      <- read.csv('data/training.csv', stringsAsFactors = FALSE)
  localization <- read.csv('data/localization.csv', stringsAsFactors = FALSE)
  
  avg_reaches      <- aggregate(reachdeviation_deg ~ pair, data=reaches,      FUN=mean, na.rm=TRUE)
  avg_localization <- aggregate(taperror_deg ~ pair,       data=localization, FUN=mean, na.rm=TRUE)
  
  FUN = mean
  CI_reaches      <- aggregate(reachdeviation_deg ~ pair, data=reaches,      FUN=function(d) Reach::getConfidenceInterval(d, method='t', FUN=FUN))
  CI_localization <- aggregate(taperror_deg ~ pair,       data=localization, FUN=function(d) Reach::getConfidenceInterval(d, method='t', FUN=FUN))
  
  polygon(
    x = c(CI_reaches$pair, rev(CI_reaches$pair)),
    y = -1 * c(CI_reaches$reachdeviation_deg[,1], rev(CI_reaches$reachdeviation_deg[,2])),
    border=NA,
    col=Reach::colorAlpha('red', alpha=34)
  )
  polygon(
    x = c(CI_localization$pair, rev(CI_localization$pair)),
    y = c(CI_localization$taperror_deg[,1], rev(CI_localization$taperror_deg[,2])),
    border=NA,
    col=Reach::colorAlpha('blue', alpha=34)
  )
  
  lines( x = avg_reaches$pair,
         y = -1 * avg_reaches$reachdeviation_deg,
         col='red')
  lines( x = avg_localization$pair,
         y = avg_localization$taperror_deg,
         col='blue')
  
  legend( x=1,
          y=33,
          legend=c('training reaches','localization'),
          lty=c(1,1),
          col=c('red','blue'),
          bty='n')
  
}