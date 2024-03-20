
# the scripts in this file bootstrap fits of the exponential functions
# the "model" functions themselves live in the Reach package
# here we apply it to successive blocks
# and bootstrap results across participants


# data in data frame should already be corrected (reach deviations * -1)

# if data has 5 columns:      participant, depvar, block, rotation,  PAIR
# NOT TRIAL!
# the amount of data to pass to cluster nodes is minimized
fitBlockSeries <- function(participants, data, prepend=c()) {
  
  df <- NA
  for (pp in participants) {
    ppdf <- data[which(data$participant == pp),]
    if (is.data.frame(df)) {
      df <- rbind(df, ppdf)
    } else {
      df <- ppdf
    }
  }
  
  avg_df <- aggregate(cbind(depvar, block, rotation) ~ pair, data=df, FUN=mean)
  
  startpoint <- 0
  
  blambdas <- c()
  bN0s     <- c()
  bblocks  <- c()
  bstartpoints <- c()
  bmfs <- c()
  basymptotes <- c()
  btrials <- c()
  
  for (block in c(1:max(avg_df$block))) {
  # for (block in c(1:3)) {
    
    bdf <- avg_df[which(avg_df$block == block),]
    
    asymptote <- bdf$rotation[1]
    signal <- c(prepend, bdf$depvar)
    
    asymptote = asymptote - startpoint # do we even need this?
    signal    = signal    - startpoint
    mf        = 1
    if (asymptote < startpoint) {
      asymptote = -1 * asymptote
      signal    = -1 * signal
      mf        = -1
    }
    
    signal <- signal / asymptote # N0 should be in the range 
    
    # what else?
    bpar <- Reach::exponentialFit(signal=signal)
    

    
    
    bblocks  <- c(bblocks,  block)
    blambdas <- c(blambdas, bpar['lambda'])
    bN0s     <- c(bN0s,     bpar['N0'])
    bstartpoints <- c(bstartpoints, startpoint)
    bmfs <- c(bmfs, mf)
    basymptotes <- c(basymptotes, asymptote)
    btrials <- c(btrials, length(signal))
    
    
    
    # set `startpoint` for the next block:
    predict <- Reach::exponentialModel(par=bpar, timepoints=length(signal)+1)
    startpoint <- startpoint + (mf * predict$output[length(signal)+1] * asymptote)
    
  }
  
  return(data.frame('block'=bblocks,
                    'startpoint'=bstartpoints,
                    'mf'=bmfs,
                    'asymptote'=basymptotes,
                    'trials'=btrials,
                    'lambda'=blambdas,
                    'N0'=bN0s))

}

tryBlockSeriesFit <- function() {
  
  locdf <- read.csv('data/localization.csv', stringsAsFactors = FALSE)

  locdf$depvar <- locdf$taperror_deg
  locdf <- locdf[,c('participant', 'block', 'rotation', 'pair', 'depvar')]
  participants <- unique(locdf$participant)
  
  bs <- fitBlockSeries(  participants = participants,
                         data         = locdf,
                         prepend      = c(0)                 )
  # 
  # write.csv(bs, 'data/bs.csv', quote = F, row.names = F)
  
  # bs <- read.csv('data/bs.csv', stringsAsFactors = FALSE)
  
  # print(bs)
  
  # plot(-1000,-1000,
  #      xlim=c(0,481), ylim=c(-30,30),
  #      main='', xlab='', ylab='',
  #      bty='n', ax=F)
  
  plotParadigm(xticks='block',ylab='angular deviation [°]')
  
  avg_localization <- aggregate(depvar ~ pair, data=locdf, FUN=mean, na.rm=TRUE)
  FUN = mean
  CI_localization <- aggregate(depvar ~ pair, data=locdf, FUN=function(d) Reach::getConfidenceInterval(d, method='t', FUN=FUN))
  
  polygon(
    x = c(CI_localization$pair, rev(CI_localization$pair)),
    y = c(CI_localization$depvar[,1], rev(CI_localization$depvar[,2])),
    border=NA,
    col=Reach::colorAlpha('blue', alpha=34)
  )
  lines( x = avg_localization$pair,
         y = avg_localization$depvar,
         col='blue')
  
  
  trial <- 50
  for (block_idx in c(1:dim(bs)[1])) {
    bpar <- c( 'lambda' = bs$lambda[block_idx],
               'N0'     = bs$N0[block_idx]      )
    
    btrials <- bs$trials[block_idx] - 1 # unprepend a 0?
    
    expfit <- Reach::exponentialModel(par=bpar, timepoints=btrials)$output
    expfit <- (expfit * bs$mf[block_idx] * bs$asymptote[block_idx]) + bs$startpoint[block_idx]
    
    lines(x=c(1:length(expfit))+trial,
          y=expfit,
          col='red')
    
    
    trial <- trial + btrials
  }
  
}


bootstrapAll <- function(bootstraps=5000, bootstrapfile=NULL, tasks=c('localization','training')) {
  
  ncores   <- parallel::detectCores()
  # clust    <- parallel::makeCluster(max(c(1,floor(ncores*0.75))))
  
  backend <- parabar::start_backend(cores = max(c(1,floor(ncores*0.75))), 
                                    cluster_type = "psock", 
                                    backend_type = "async")
  
  locdf <- read.csv('data/localization.csv', stringsAsFactors = FALSE)
  
  locdf$depvar <- locdf$taperror_deg
  locdf <- locdf[,c('participant', 'block', 'rotation', 'pair', 'depvar')]
  
  if (is.null(bootstrapfile)) {
    
    participants <- unique(locdf$participant)
    
    BSparticipants <- sample(x = participants,
                             size = length(participants) * bootstraps,
                             replace = TRUE)
    
    BSparticipants <- matrix( data = BSparticipants,
                              nrow = bootstraps,
                              ncol = length(participants))
    
    write.csv(BSparticipants, 'data/bootstrap_participants.csv', quote=TRUE, row.names=FALSE)
    
  } else {
    
    BSparticipants <- as.matrix(read.csv('data/bootstrap_participants.csv', stringsAsFactors = FALSE))
    
  }
  
  
  # a <- parallel::parApply(cl = clust,
  #                         X = BSparticipants,
  #                         MARGIN = 1,
  #                         FUN = fitBlockSeries,
  #                         data = locdf,
  #                         prepend = c(0))
  
  if ('localization' %in% tasks) {
    cat('\nLOCALIZATION BOOTSTRAPS\n\n')
    
    a <- parabar::par_apply(backend = backend,
                            x = BSparticipants,
                            margin = 1,
                            fun = fitBlockSeries,
                            data = locdf,
                            prepend = c(0) )
    
    all_bootstraps <- NA
    for (bs in c(1:length(a))) {
      bsdf <- a[[bs]]
      bsdf$bootstrap <- bs
      if (is.data.frame(all_bootstraps)) {
        all_bootstraps <- rbind(all_bootstraps, bsdf)
      } else {
        all_bootstraps <- bsdf
      }
    }
    
    write.csv(all_bootstraps, 'data/localization_bootstraps.csv', quote=F, row.names=F)
    
  }
  
  # TRAINING
  
  if ('training' %in% tasks) {
    
    traindf <- read.csv('data/training.csv', stringsAsFactors = FALSE)
    
    traindf$depvar <- traindf$reachdeviation_deg
    traindf <- traindf[,c('participant', 'block', 'rotation', 'pair', 'depvar')]
    
    cat('\nREACH TRAINING BOOTSTRAPS\n\n')
    
    a <- parabar::par_apply(backend = backend,
                            x = BSparticipants,
                            margin = 1,
                            fun = fitBlockSeries,
                            data = traindf,
                            prepend = c())
    
    all_bootstraps <- NA
    for (bs in c(1:length(a))) {
      bsdf <- a[[bs]]
      bsdf$bootstrap <- bs
      if (is.data.frame(all_bootstraps)) {
        all_bootstraps <- rbind(all_bootstraps, bsdf)
      } else {
        all_bootstraps <- bsdf
      }
    }
    
    write.csv(all_bootstraps, 'data/training_bootstraps.csv', quote=F, row.names=F)
    
  }
  
  # parallel::stopCluster(clust)
  
  parabar::stop_backend(backend)
  
}


# bootstrapTraining <- function(bootstraps=5000) {
#   
#   ncores   <- parallel::detectCores()
#   clust    <- parallel::makeCluster(max(c(1,floor(ncores*0.75))))
#   
# 
#   traindf <- read.csv('data/training.csv', stringsAsFactors = FALSE)
#   
#   traindf$depvar <- traindf$reachdeviation_deg
#   traindf <- traindf[,c('participant', 'block', 'rotation', 'pair', 'depvar')]
#   participants <- unique(traindf$participant)
#   
#   BSparticipants <- sample(x = participants,
#                            size = length(participants) * bootstraps,
#                            replace = TRUE)
#   
#   BSparticipants <- matrix( data = BSparticipants,
#                             nrow = bootstraps,
#                             ncol = length(participants))
#   
#   a <- parallel::parApply(cl = clust,
#                           X = BSparticipants,
#                           MARGIN = 1,
#                           FUN = fitBlockSeries,
#                           data = locdf,
#                           prepend = c(0))
#   
#   parallel::stopCluster(clust)
#   
#   all_bootstraps <- NA
#   for (bs in c(1:length(a))) {
#     bsdf <- a[[bs]]
#     bsdf$bootstrap <- bs
#     if (is.data.frame(all_bootstraps)) {
#       all_bootstraps <- rbind(all_bootstraps, bsdf)
#     } else {
#       all_bootstraps <- bsdf
#     }
#   }
#   
#   write.csv(all_bootstraps, 'data/localization_bootstraps.csv', quote=F, row.names=F)
#   return(all_bootstraps)
#   
# }