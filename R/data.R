
# getting all the data:
# (need to add processed data to this)

getData <- function() {
  
  Reach::downloadOSFdata(  repository = '5ec7s',
                           filelist   = list('data'=c('raw.zip','demographics.csv', 'processed.zip', 'bootstrap.zip')),
                           folder     = 'data',
                           unzip      = TRUE  )
  
}

# combine data -----

combineLocalizationData <- function() {
  
  demographics <- read.csv('data/demographics.csv', stringsAsFactors = FALSE)
  participants <- demographics$participant
  
  schedule <- getSchedule()
  schedule <- schedule[which(schedule$type == 'localization'),]
  schedule <- schedule[,!names(schedule) %in% c("type")]
  
  locdf <- NA
  
  for (pp in participants) {
    
    # read data
    ppdf <- read.csv(sprintf('data/%s/localization.csv',pp), stringsAsFactors = FALSE)
    
    rm.idx <- c(   which(ppdf$tapx_cm < -10.5),
                   which(ppdf$tapy_cm > 14.5),
                   which(ppdf$tapy_cm < 0),
                   which(sqrt((ppdf$tapx_cm)^2 + (ppdf$tapy_cm)^2) < 5) )
    
    if (length(rm.idx) > 0) {
      rm.idx <- unique(rm.idx)
      ppdf$tapx_cm[rm.idx] <- NA
      ppdf$tapy_cm[rm.idx] <- NA
    }
    
    # merge with schedule info:
    ppdf <- merge(schedule, ppdf)
    ppdf <- ppdf[order(ppdf$trial),]
    ppdf$tapangle_deg <- ((atan2(ppdf$tapy_cm, ppdf$tapx_cm)/pi)*180)
    
    ppdf$taperror_deg <- ppdf$tapangle_deg - ppdf$targetangle_deg
    
    # baseline?
    ppdf$taperror_deg <- ppdf$taperror_deg - mean( ppdf$taperror_deg[c(25:49)] , na.rm=TRUE)
    
    ppdf$tapangle_deg <- round( ppdf$targetangle_deg + ppdf$taperror_deg, digits=4 )
    ppdf$taperror_deg <- round( ppdf$taperror_deg, digits=4 )
        
    if (is.data.frame(locdf)) {
      locdf <- rbind(locdf, ppdf)
    } else {
      locdf <- ppdf
    }
  }
  
  write.csv(locdf, 'data/localization.csv', quote=FALSE, row.names=FALSE)
  # return(locdf)
  
}

getSchedule <- function() {
  
  rotations	<- c(0,	-30,	-15,	0,	-15,	15,	0,	15,	-30,	0,	-30,	0,	30,	0,	30,	-15,	30,	15,	-15,	0,	-15,	-30,	30,	0,	0,	15,	30,	-30,	15,	0,	0)
  clamps	  <- c(0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0)
  pairs	    <- c(49,	12,	12,	12,	12,	12,	24,	24,	24,	12,	12,	12,	12,	24,	12,	24,	24,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	11)
  
  # now expand to full, trial-by-trial schedule
  trial    <- 1:(sum(pairs)*2)
  pair     <- rep(1:(sum(pairs)), each=2)
  block    <- rep(c(1:length(pairs))-1, times=pairs*2)
  type     <- rep(c('reach','localization'),sum(pairs))
  rotation <- rep(rotations, pairs*2)
  clamp    <- rep(clamps, pairs*2)
  
  # remove from localization data schedule?
  # technically yes, but it's better if we have the preceding reach info available
  
  # rotation[which(type == 'localization')] <- NA
  # clamp[which(type == 'localization')] <- NA
  
  return( data.frame( trial,
                      pair,
                      block,
                      type,
                      rotation,
                      clamp      ) )

}

combineTrainingData <- function() {
  
  demographics <- read.csv('data/demographics.csv', stringsAsFactors = FALSE)
  participants <- demographics$participant
  
  schedule <- getSchedule()
  schedule <- schedule[which(schedule$type == 'reach'),]
  schedule <- schedule[,!names(schedule) %in% c("type")]
  
  traindf <- NA
  
  for (pp in participants) {
    
    # read data
    ppdf <- read.csv(sprintf('data/%s/training.csv',pp), stringsAsFactors = FALSE)
    
    # calculate distance from home:
    ppdf$distance <- sqrt(ppdf$handx_cm^2 + ppdf$handy_cm^2)
    # print(ppdf$distance)
    # remove early samples:
    ppdf <- ppdf[which(ppdf$distance > 3),]
    
    # select only the first sample from each trial that's left:
    ppdf$row.idx <- c(1:dim(ppdf)[1])
    row.idx <- aggregate(row.idx ~ trial, data=ppdf, FUN=min)$row.idx
    ppdf <- ppdf[which(ppdf$row.idx %in% row.idx),]
    
    # remove temp columns:
    ppdf <- ppdf[,!names(ppdf) %in% c("row.idx", "distance")]
    
    
    ppdf$reachangle_deg <- NA
    ppdf$reachdeviation_deg <- NA
    
    
    for (targetangle in unique(ppdf$targetangle_deg)) {
      
      # which trials:
      idx <- which(ppdf$targetangle_deg == targetangle)
      
      # reach deviations:
      theta <- -1 * (targetangle / 180) * pi
      R <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),nrow=2,ncol=2)
      XY <- matrix(c(ppdf$handx_cm[idx], ppdf$handy_cm[idx]), byrow=TRUE, nrow=2,ncol=length(idx) )
      XYr <- R %*% XY
      
      # extract reach deviation:
      ppdf$reachdeviation_deg[idx] <- (atan2(XYr[2,],XYr[1,])/pi)*180

    }
    
    # baseline ?
    ppdf$reachdeviation_deg <- ppdf$reachdeviation_deg - mean(ppdf$reachdeviation_deg[c(25:49)], na.rm=TRUE)
    
    ppdf$reachangle_deg     <- round(ppdf$targetangle_deg + ppdf$reachdeviation_deg, digits=4)
    ppdf$reachdeviation_deg <- round(ppdf$reachdeviation_deg,                        digits=4)
    
    
    ppdf <- merge(schedule, ppdf)
    ppdf <- ppdf[order(ppdf$trial),]
    
    
    if (is.data.frame(traindf)) {
      traindf <- rbind(traindf, ppdf)
    } else {
      traindf <- ppdf
    }
    
  }
  
  write.csv(traindf, 'data/training.csv', quote=FALSE, row.names=FALSE)
  # return(traindf)
  
}

participantDescriptors <- function() {
  
  demographics <- read.csv('data/demographics.csv', stringsAsFactors = FALSE)
  
  cat('\n--- age\n\n')
  cat(sprintf('mean: %0.2f\n',mean(demographics$age)))
  cat(sprintf('sd: %0.3f\n',sd(demographics$age)))
  age_range <- range(demographics$age)
  cat(sprintf('range: %d-%d\n',age_range[1],age_range[2]))
  
  for (desc in c('sex','handedness','vision')) {
    cat(sprintf('\n--- %s\n',desc))
    print(table(demographics[,desc]))
  }
  
  
}
