
# getting all the data:
# (need to add processed data to this)

getData <- function() {
  
  Reach::downloadOSFdata(  repository = '5ec7s',
                           filelist   = list('data'=c('raw.zip','demographics.csv')),
                           folder     = 'data',
                           unzip      = TRUE  )
  
}

# combine data -----

combineLocalizationData <- function() {
  
  demographics <- read.csv('data/demographics.csv', stringsAsFactors = FALSE)
  participants <- demographics$participant
  
  for (pp in participants) {
    
    ppdf <- read.csv(sprintf('data/%s/localization.csv',pp), stringsAsFactors = FALSE)
    
    names(ppdf)[which(names(ppdf) == 'block')] <- 'pair'
    
    
    
  }
  
}

getSchedule <- function() {
  
  rotations	<- c(0,	-30,	-15,	0,	-15,	15,	0,	15,	-30,	0,	-30,	0,	30,	0,	30,	-15,	30,	15,	-15,	0,	-15,	-30,	30,	0,	0,	15,	30,	-30,	15,	0,	0)
  clamps	  <- c(0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0)
  pairs	    <- c(49,	12,	12,	12,	12,	12,	24,	24,	24,	12,	12,	12,	12,	24,	12,	24,	24,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12)
  
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

