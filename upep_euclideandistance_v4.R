# Euclidean distances from Vp_Extract.m parameters

# packages & functions ----------------------------------------------------

library(here)
library(tidyr)
library(dplyr)
library(tibble)
library(data.table)
library(ggplot2)
library(scales)
library(ggbeeswarm)

source('~/Dropbox/phd/framebyframe/genotypeUtilities.R')
source('~/Dropbox/phd/framebyframe/pathUtilities.R')

sem <- function(x) sd(x)/sqrt(length(x))

# get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


# some common settings/choices --------------------------------------------

# days/nights
dani <- c('day1', 'night1', 'day2', 'night2')

# Vp_Extract parameters
pmts <- c('activeBoutLength', 'activeBoutMean', 'activeBoutStd',
          'activeBoutTotal', 'activeBoutMin', 'activeBoutMax',
          'numActiveBouts', 'totTimeActive', 'totActivity',
          'inactiveBoutLength', 'numInactiveBouts', 'totTimeInactive')

# Vp_Extract parameters which we will use
# last s for selection
pmtss <- c('activeBoutLength', 'activeBoutMean', 'activeBoutStd',
          'activeBoutTotal', 'activeBoutMin', 'activeBoutMax',
          'numActiveBouts', 'totTimeActive',
          'inactiveBoutLength')
# deleting
# totActivity, because do not understand plot interpretation
# numInactiveBouts, because same as numActiveBouts
# totTimeInactive, because mirror of totTimeActive

# control group
congrp <- 'WT.Cas9'

# v4: trial with pooling WT (uninjected) & WT+Cas9
poolWt <- TRUE

# function to import files for one experiment into a list -----------------

# ! expects output files from vpExMat2Params.m

# function to import files for one experiment into a dataframe
# each row is a one fish during one time window
# each column is one parameter
# each cell is average for that fish and that parameter during that time window

importParams <- function(fromMatfolder,
                         geno_path) {
  
  
  ### find & import grptags file from input folder
  grptags_path <- list.files(fromMatfolder, full.names=TRUE)[which(list.files(fromMatfolder)=='grptags.csv')]
  
  if (length(grptags_path)==0) stop('\t \t \t \t >>> Cannot find grptags.csv in input folder \n')
  
  grptags <- as.vector(unlist(read.csv(grptags_path, header=FALSE)))
  
  
  ### import genotype file
  geno <- importGenotype(geno_path)
  # we only need the genotype names so we can replace fishids by group names in grptags
  
  # below to replace fishids by group names in grptags
  grptags <- sapply(grptags, function(g) {
    colnames(geno)[g]
  })
  
  
  ### import parameter matrices
  # will import in a list, each slot is a time window
  # namely, day1, night1, day2, night2
  # (will skip night0)
  
  # preallocate list
  pamdn <- vector(mode='list', length=length(dani)) # parameter day/night
  names(pamdn) <- dani
  
  # now fill it in
  # w for time window
  for (w in 1:length(dani)) {
    
    # is w a day or a night?
    # look if first character is d or n
    char1 <- strsplit(dani[w], '')[[1]][1]
    if (char1=='d') {
      daynight <- 'day'
    } else if (char1=='n') {
      daynight <- 'night'
    } else {
      stop('\t \t \t \t >>> Error: one of the days/nights does not start by d or n /n')
    }
    
    # import parameter matrix for that time window
    # filename is e.g. pam_day1.csv
    pa <- read.csv(here(fromMatfolder, paste0('pam_', dani[w], '.csv')), header=FALSE)
    
    # add column names
    colnames(pa) <- pmts
    
    # select only parameter we want
    pa <- pa[, colnames(pa) %in% pmtss]
    
    # add fishid, group, time window, day/night
    # add columns from right to left
    pa <- pa %>%
      add_column(fishid=sprintf('f%i', 1:nrow(pa)), .before=1) %>%
      add_column(grp=grptags, .before=1) %>%
      add_column(win=dani[w], .before=1) %>%
      add_column(daynight=daynight, .before=1)
    
    # ! fishid will not correspond to actual well number as I usually do
    # because Vp_Extract.m deletes data from empty well
    # I find it easier to keep all of the data until late
    
    # add pa to list pamdn
    pamdn[[w]] <- pa
  }
  
  # rbindlist pamdn so makes one dataframe
  pml <- rbindlist(pamdn) # pam long format
  
  # now pml is a sensical, easy-to-read dataframe!
  return(pml)
  
}


# import all experiments --------------------------------------------------

# from now on
# exp1 = 220212_01
# exp2 = 220212_02
# exp3 = 220215_01
# exp4 = 220215_02
# exp5 = 220520_01
# exp6 = 220524_01

### exp1
exp1 <- importParams(fromMatfolder=here('exp220212', '220212_01_fromMat'),
                     geno_path=here('exp220212', '220212_01genotype.txt'))
# add experiment tag
exp1 <- exp1 %>%
  add_column(exp='exp1', .before='daynight')

### exp2
exp2 <- importParams(fromMatfolder=here('exp220212', '220212_02_fromMat'),
                     geno_path=here('exp220212', '220212_02genotype.txt'))
# add experiment tag
exp2 <- exp2 %>%
  add_column(exp='exp2', .before='daynight')

### exp3
exp3 <- importParams(fromMatfolder=here('exp220215', '220215_01_fromMat'),
                     geno_path=here('exp220215', '220215_01genotype.txt'))
# add experiment tag
exp3 <- exp3 %>%
  add_column(exp='exp3', .before='daynight')

### exp4
exp4 <- importParams(fromMatfolder=here('exp220215', '220215_02_fromMat'),
                     geno_path=here('exp220215', '220215_02genotype.txt'))
# add experiment tag
exp4 <- exp4 %>%
  add_column(exp='exp4', .before='daynight')

### exp5
exp5 <- importParams(fromMatfolder=here('exp220520', '220520_01_fromMat'),
                     geno_path=here('exp220520', '220520_01genotype.txt'))
# add experiment tag
exp5 <- exp5 %>%
  add_column(exp='exp5', .before='daynight')
# small correction: in other exps, uninjected group is 'WT', but in exp5 it is 'WT-uninj'
# correct this
exp5[which(exp5$grp=='WT.uninj'), 'grp'] <- 'WT'

### exp6
exp6 <- importParams(fromMatfolder=here('exp220524', '220524_01_fromMat'),
                     geno_path=here('exp220524', '220524_01genotype.txt'))
# add experiment tag
exp6 <- exp6 %>%
  add_column(exp='exp6', .before='daynight')

# now append all experiments into one big dataframe
upep <- rbindlist(list(exp1=exp1, exp2=exp2, exp3=exp3, exp4=exp4, exp5=exp5, exp6=exp6))
# upep for micropeptides

# v4: trial with pooled WT (uninjected) & WT+Cas9 for normalisation
# we leave everything as is, but 'pretend' that any WT (uninjected) larva was WT+Cas9 by replacing its label

if (poolWt) {
  upep[which(upep$grp == 'WT') , 'grp'] <- 'WT.Cas9'
}

# paste exp and grp into one identifier
upep <- upep %>%
  add_column(exp_grp=paste(upep$exp, upep$grp, sep='_'), .after='grp')


# convert to long format --------------------------------------------------
# for normalisation, convert upep in a tall dataframe

# pivot_longer to a tall/long format
upel <- upep %>% # for upep, long format
  pivot_longer(-c(exp, daynight, win, grp, exp_grp, fishid),
               names_to='pam_id',
               values_to='param')

# for each fish and parameter, calculate mean across both days or across both nights
upem <- upel %>% # for upep, mean across days / nights
  group_by(exp, daynight, grp, exp_grp, fishid, pam_id) %>% # ***
  summarise_at(vars(param),
               list(
                 param_mean= ~ mean(.)
               ))
# *** note we do not group_by win, so night1/night2 or day1/day2 ends up in the same 'bag' for mean below

# z-score to controls -----------------------------------------------------

# for each parameter, calculate mean & sd of control fish
# for any parameter/fish, Z-score from control = value - control mean / sd

# first prepare small dataframe mean / sd of control fish across parameters
# i.e. rows = parameters, column = mean / sd of all control fish

con <- upem %>%
  subset(grp==congrp) %>% # congrp defined at the start; group_by exp because we will normalise to controls within each experiment (clutch)
  group_by(exp, daynight, pam_id) %>%
  summarise_at(vars(param_mean),
               list(
                 con_mean= ~ mean(.),
                 con_sd= ~ sd(.),
                 con_n= ~ length(.) # as a sanity check that we have the right number of control fish
               ))
# this creates 'control' dataframe
# i.e. mean / sd of all control larvae within each experiment

# v4 pooling WT & WT+Cas9: smallest control group is N = 22

# small helper function to Z-score one datapoint
zScore2Control <- function(experiment,
                           dayornight,
                           parameter,
                           datapoint) {
  
  # from control dataframe, get mean and sd of controls for that experiment, that parameter
  conmean <- con %>%
    subset(exp==experiment & daynight==dayornight & pam_id==parameter, 'con_mean')
  
  consd <- con %>%
    subset(exp==experiment & daynight==dayornight & pam_id==parameter, 'con_sd')
  
  # now calculate Z-score for the datapoint we were given
  zsco <- (datapoint - conmean) / consd
  
  # make sure it is a simple number and nothing extra
  zsco <- as.double(zsco)
  
  # return z-score
  return(zsco)
  
}


# now go through row by row and Z-score datapoints to controls
# Note, we also Z-score datapoints from control fish, i.e. to their own mean/sd
upem$z2con <- apply(upem, 1, function(row) {
  
  zScore2Control(experiment=as.character(row['exp']),
                 dayornight=as.character(row['daynight']),
                 parameter=as.character(row['pam_id']),
                 datapoint=as.double(row['param_mean']))
  
})


# prepare fingerprints ----------------------------------------------------

# create new column daynight_pamid, e.g. day_activeBoutLength
# easier, as we are essentially dealing with parameters during day and during night as different parameters
upem <- upem %>%
  add_column(daynight_pamid=paste(upem$daynight, upem$pam_id, sep='_'), .before='param_mean')
# Note, I checked vs code from Marcus, same values at this stage

# we have Z-scores normalised to controls
# summarise (mean, sd, sem, n) across fish in each group
fgp <- upem %>%
  group_by(exp, daynight, grp, exp_grp, daynight_pamid) %>%
  summarise_at(vars(z2con),
               list(
                 mean= ~ mean(.),
                 sd= ~ sd(.),
                 sem= ~ sem(.),
                 nfis= ~ length(.)
               ))
# fgp for fingerprints


# plot fingerprints -------------------------------------------------------

# add a column daynight_grp will allow us to disconnect day parameters / night parameters
fgp <- fgp %>%
  add_column(exp_grp_daynight=paste(fgp$exp_grp, fgp$daynight, sep='_'), .before='daynight_pamid')

# function to plot fingerprints of specific group

ggFingerprint <- function(fgp,
                          whichexps='all',
                          whichgrps='all',
                          connectOrNo=TRUE,
                          removeControls=FALSE,
                          controlGroup=NA,
                          colours=NA,
                          legendOrNo=TRUE,
                          exportOrNo=FALSE,
                          exportname,
                          width=150,
                          height=100) {
  
  fgptmp <- fgp # fgp, temporary copy for plotting
  
  ### if we are not plotting every experiment, only keep the ones we want
  if (whichexps[1] != 'all') {
    fgptmp <- subset(fgptmp, exp %in% whichexps)
  }
  
  ### if we are not plotting every group, only keep the ones we want
  if (whichgrps[1] != 'all') {
    fgptmp <- subset(fgptmp, grp %in% whichgrps)
  }
  
  ### if we are not plotting controls, remove them
  if (removeControls) {
    fgptmp <- subset(fgptmp, ! grp %in% controlGroup) # i.e. only keep groups which are not control groups
  }
  
  ### if user did provide any colours, we use automatic ggplot colours
  if(is.na(colours[1])) {
    colours <- hue_pal()(length(unique(fgptmp$exp_grp)))
  }
  
  ### plot
  gg_fgp <- ggplot(fgptmp, aes(x=daynight_pamid, y=mean, colour=exp_grp)) +
    geom_hline(yintercept=0, linetype=1, colour='#a7a7a7', size=0.5) +
    geom_pointrange(aes(ymin=mean-sem, ymax=mean+sem), size=0.3) +
    {if(connectOrNo) geom_line(aes(group=exp_grp_daynight), size=0.2)} +
    scale_colour_manual(values=colours) + 
    theme_minimal() +
    theme(
      panel.grid.minor.y=element_blank(),
      axis.title.x=element_text(size=9),
      # axis.title.y=element_text(size=9),
      axis.title.y=element_blank(),
      axis.text.x=element_text(size=7),
      axis.text.y=element_text(size=7)
    ) +
    
    {if(!legendOrNo) theme(legend.position='none')} +
    
    coord_cartesian(ylim=c(-2.4, 2.4)) +
    ylab('Deviation from controls') +
    xlab('Parameter') +
    scale_x_discrete(labels=c(seq(1, 9), seq(1, 9)))
  
  ### save plot, if required
  if (exportOrNo) {
    ggsave(filename=exportname, gg_fgp, width=width, height=height, units='mm')
  }
  
  ### tell user the Ns
  cat('\t \t \t \t >>> Sample sizes \n')
  fgptmp %>%
    group_by(exp_grp) %>%
    distinct(nfis) %>%
    print()
  
  ### pairwise correlation between fingerprints plotted
  
  # need to pivot to wide format, where each column is an exp_grp, each row a parameter, and each cell the mean Z-score
  fgw <- fgptmp %>% # fgptmp, wider format
    group_by(exp_grp) %>%
    pivot_wider(-c(exp, grp, exp_grp_daynight, sd, sem, nfis), names_from=exp_grp, values_from=mean)
  
  # now calculate pairwise correlations between columns (data for each exp_grp)
  fcor <- as.data.frame(cor(fgw[3:ncol(fgw)])) # default is Pearson
  fcor <- get_upper_tri(fcor)
  
  # print the correlation matrix
  cat('\t \t \t \t >>> Pairwise correlations between fingerprints plotted \n')
  print(fcor)

  # gives a summary of the correlation coefficients as mean + sd
  rs <- as.vector(unlist(fcor)) # rs for r coefficients
  # remove any 1.0 (correlation of fingerprint vs itself) and any NA
  rs <- rs[-which(rs==1)]
  rs <- rs[-which(is.na(rs))]
  # tell user mean + sd of correlation coefficients left
  cat('\t \t \t \t >>> Mean & sd of correlation coefficients: ', round(mean(rs),2), round(sd(rs),2), '\n')

  ### return plot
  # so displays in RStudio
  return(gg_fgp)
  
}

dir.create(here('fingerprints_v4'))

# WT.Cas9 only to show normalisation
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='WT.Cas9',
              legendOrNo=TRUE,
              exportOrNo=TRUE,
              exportname=here('fingerprints_v4', 'fgp_WTCas9.pdf'),
              width=150,
              height=100)

# all WT to compare uninjected vs injected
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='WT',
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_WTuninj.pdf'),
              width=150,
              height=100)

###  when available, replication between boxes/clutches

# dreammist
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='dreammist',
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_dreammist.pdf'),
              width=150,
              height=100)

# libra.F0
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='libra.F0',
              connectOrNo=TRUE,
              colours=c('#fcb505', '#fcd98a'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_libraf0.pdf'),
              width=60,
              height=46)


# linc.cd74
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='linc.cd74',
              connectOrNo=TRUE,
              colours=c('#78ac63', '#bbd4ae'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_linccd74.pdf'),
              width=60,
              height=46)

# linc.foxp2
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='linc.foxp2',
              connectOrNo=TRUE,
              colours=c('#cb2a20', '#e4928f'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_lincfoxp2.pdf'),
              width=60,
              height=46)

# linc.loc485.F0
ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps='linc.loc485.F0',
              connectOrNo=TRUE,
              colours=c('#982150', '#cb8fa8'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_lincloc485.pdf'),
              width=60,
              height=46)

# Linc.112305.F0
ggFingerprint(fgp=fgp,
              whichexps='exp4',
              whichgrps='Linc.112305.F0',
              connectOrNo=TRUE,
              colours=c('#ec6354'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_linc112305.pdf'),
              width=60,
              height=46)

# Linc.239.F0
ggFingerprint(fgp=fgp,
              whichexps='exp4',
              whichgrps='Linc.239.F0',
              connectOrNo=TRUE,
              colours=c('#417dcd'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_linc239.pdf'),
              width=60,
              height=46)

# linc.mettl3.F0
ggFingerprint(fgp=fgp,
              whichexps='exp4',
              whichgrps='linc.mettl3.F0',
              connectOrNo=TRUE,
              colours=c('#7e7ba4'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_lincmettl3.pdf'),
              width=60,
              height=46)

# linc.oncecut.F0
ggFingerprint(fgp=fgp,
              whichexps='exp4',
              whichgrps='linc.oncecut.F0',
              connectOrNo=TRUE,
              colours=c('#db5072'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_linconcecut.pdf'),
              width=60,
              height=46)

# exp5
ggFingerprint(fgp=fgp,
              whichexps='exp5',
              whichgrps='linc.mipep1.ex3inj',
              connectOrNo=TRUE,
              colours=c('#23aa4c'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_mipep1.pdf'),
              width=60,
              height=46)

# exp6
ggFingerprint(fgp=fgp,
              whichexps='exp6',
              whichgrps='linc.wrb.F0',
              connectOrNo=TRUE,
              colours=c('#e94537'),
              legendOrNo=FALSE,
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_wrb.pdf'),
              width=60,
              height=46)


### experiment by experiment
ggFingerprint(fgp=fgp,
              whichexps='exp1',
              whichgrps=c('dreammist', 'libra.F0', 'linc.cd74', 'linc.foxp2', 'linc.loc485.F0'), # ***
              removeControls=TRUE,
              controlGroup='WT.Cas9',
              exportOrNo=TRUE,
              exportname=here('fingerprints_v4', 'fgp_exp1.pdf'),
              width=150,
              height=100)
# *** all groups except WT, for fair comparison vs analysis v4 where we normalise to all WT

ggFingerprint(fgp=fgp,
              whichexps='exp2',
              whichgrps=c('dreammist', 'libra.F0', 'linc.cd74', 'linc.foxp2', 'linc.loc485.F0'), # ***
              removeControls=TRUE,
              controlGroup='WT.Cas9',
              exportOrNo=TRUE,
              exportname=here('fingerprints_v4', 'fgp_exp2.pdf'),
              width=150,
              height=100)
# *** all groups except WT, for fair comparison vs analysis v4 where we normalise to all WT

ggFingerprint(fgp=fgp,
              whichexps='exp3',
              whichgrps=c('Linc.112305.F0', 'Linc.239.F0', 'linc.mettl3.F0', 'linc.oncecut.F0'), # ***
              removeControls=TRUE,
              controlGroup='WT.Cas9',
              exportOrNo=TRUE,
              exportname=here('fingerprints_v4', 'fgp_exp3.pdf'),
              width=150,
              height=100)
# *** all groups except WT, for fair comparison vs analysis v4 where we normalise to all WT

ggFingerprint(fgp=fgp,
              whichexps='exp4',
              whichgrps=c('Linc.112305.F0', 'Linc.239.F0', 'linc.mettl3.F0', 'linc.oncecut.F0'), # ***
              removeControls=TRUE,
              controlGroup='WT.Cas9',
              exportOrNo=TRUE,
              exportname=here('fingerprints_v4', 'fgp_exp4.pdf'),
              width=150,
              height=100)
# *** all groups except WT, for fair comparison vs analysis v4 where we normalise to all WT

ggFingerprint(fgp=fgp,
              whichexps='exp5',
              whichgrps='all',
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_exp5.pdf'),
              width=150,
              height=100)

ggFingerprint(fgp=fgp,
              whichexps='exp6',
              whichgrps='all',
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_exp6.pdf'),
              width=150,
              height=100)

ggFingerprint(fgp=fgp,
              whichexps='all',
              whichgrps=c('linc.mipep1.ex3inj', 'linc.wrb.F0'),
              exportOrNo=TRUE,
              exportname=here('fingerprints', 'fgp_wrbmipep.pdf'),
              width=150,
              height=100)

fgp=fgp
whichexps='all'
whichgrps=c('linc.mipep1.ex3inj', 'linc.wrb.F0')
exportOrNo=TRUE
exportname=here('fingerprints', 'fgp_wrbmipep.pdf')
width=150
height=100

fgw2 <- fgw %>%
  subset(daynight=='day')

as.data.frame(cor(fgw2[3:ncol(fgw2)]))

fgptmp2 <- fgptmp %>%
  subset(daynight=='day')

gg_fgp <- ggplot(fgptmp, aes(x=daynight_pamid, y=mean, colour=exp_grp)) +
  geom_hline(yintercept=0, linetype=1, colour='#a7a7a7', size=0.5) +
  geom_pointrange(aes(ymin=mean-sem, ymax=mean+sem), size=0.3) +
  {if(connectOrNo) geom_line(aes(group=exp_grp_daynight), size=0.2)} +
  scale_colour_manual(values=c('#23aa4c', '#e94537')) + 
  theme_minimal() +
  theme(
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_text(size=9),
    # axis.title.y=element_text(size=9),
    axis.title.y=element_blank(),
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=7)
  ) +
  
  {if(!legendOrNo) theme(legend.position='none')} +
  
  coord_cartesian(ylim=c(-2.4, 2.4)) +
  ylab('Deviation from controls') +
  xlab('Parameter') +
  scale_x_discrete(labels=c(seq(1, 9), seq(1, 9)))

ggsave(here('fingerprints', 'fpg_mipep1wrb.pdf'), width=60, height=46, units='mm')



# for each fish, euclidean distance from controls -------------------------

# convert back to wide
upew <- upem %>% # upew for upep wide
  pivot_wider(id_cols=c(fishid, exp, grp, exp_grp),
              names_from=daynight_pamid,
              values_from=z2con)

# illustrate two dimensions of N-dimension space (N = number of parameters)
# e.g. day_activeBoutLengtth and day_activeBoutMean

upew2 <- upew %>%
  subset(exp_grp %in% c('exp1_WT.Cas9', 'exp2_WT.Cas9', 'exp1_linc.loc485.F0', 'exp2_linc.loc485.F0'))

twodim_eg <- ggplot(upew2, aes(x=day_activeBoutLength, y=night_numActiveBouts, colour=exp_grp)) +
  geom_point() +
  coord_cartesian(ylim=c(-2.5, 2.5), xlim=c(-2.5, 2.5)) +
  theme_minimal()
twodim_eg
ggsave(filename=here('twodim_example.pdf'), twodim_eg, width=150, height=100, units='mm')


# each fish has a fingerprint, which is its Z-scores from controls
# i.e. each row of upew

### sanity check:
# compute the mean fingerprint of control fish
# take fingerprints of control fish
upew_con <- upew %>%
  subset(grp==congrp)
# compute mean by column
fgp_con <- as.vector(apply(upew_con[5:ncol(upew_con)], 2, mean))
# this indeed gives ~ 0, i.e. in the multidim space, the center of the control fish is at coordinates 0, 0, 0, ...

# round the reference point's coordinates to 0
fpg_con <- round(fgp_con) # gives 18 times zero

# now, Euclidean distance for each fish is
# distance between its point (whose coordinates are a vector of 18 values, one per parameter) and the reference point at 0, 0, 0, ...

# fip for fish's point, i.e. the coordinates of that point; in practice row of upew
euc <- apply(upew[, 5:ncol(upew)], 1, function(fip) {
  dist(rbind(fgp_con, fip) , 'euclidean')
})
# euc for Euclidean distances from controls
# one value per larva

# put back fishid and grp from pmz
eud <- as.data.frame(cbind(upew$exp_grp, upew$exp, upew$grp, upew$fishid, euc))
colnames(eud) <- c('exp_grp', 'exp', 'grp', 'fishid', 'euclid')
eud$euclid <- as.numeric(eud$euclid)


# plot Euclidean distances ------------------------------------------------

# summarise Euclidean distances (mean, sd, sem, n) across fish in each group
euds <- eud %>%
  group_by(exp_grp, exp, grp) %>%
  summarise_at(vars(euclid),
               list(
                 mean= ~ mean(.),
                 sd= ~ sd(.),
                 sem= ~ sem(.),
                 nfis= ~ length(.)
               ))


# do not plot experiment 3 or dreammist
euds2 <- euds %>%
  subset(exp != 'exp3') %>%
  subset(grp != 'dreammist')

# also remove from individual Euclidean distances
eud2 <- eud %>%
  subset(exp != 'exp3') %>%
  subset(grp != 'dreammist')


### colours

# first define the exp_grp levels order, this is the order in the plot
expgrp_order <- c('exp1_WT.Cas9', 'exp1_libra.F0', 'exp1_linc.cd74', 'exp1_linc.foxp2', 'exp1_linc.loc485.F0',
                  'exp2_WT.Cas9', 'exp2_libra.F0', 'exp2_linc.cd74', 'exp2_linc.foxp2', 'exp2_linc.loc485.F0',
                  'exp4_WT.Cas9', 'exp4_Linc.112305.F0', 'exp4_Linc.239.F0', 'exp4_linc.mettl3.F0', 'exp4_linc.oncecut.F0',
                  'exp5_WT.Cas9', 'exp5_linc.mipep1.ex3inj',
                  'exp6_WT.Cas9', 'exp6_linc.wrb.F0')
eud2$exp_grp <- factor(eud2$exp_grp, levels=expgrp_order)
euds2$exp_grp <- factor(euds2$exp_grp, levels=expgrp_order)

# as this will be the order in the plot, we can also write the labels properly

expgrp_labels <- c('controls exp1',
                   expression(paste(italic('libra'), ' exp1')), 
                   expression(paste(italic('linc.cd74'), ' exp1')),
                   expression(paste(italic('linc.foxp2'), ' exp1')),
                   expression(paste(italic('linc.loc485'), ' exp1')),
                   
                   'controls exp2',
                   expression(paste(italic('libra'), ' exp2')), 
                   expression(paste(italic('linc.cd74'), ' exp2')),
                   expression(paste(italic('linc.foxp2'), ' exp2')),
                   expression(paste(italic('linc.loc485'), ' exp2')),
                   
                   'controls exp3',
                   expression(italic('linc.112305')), 
                   expression(italic('linc.239')),
                   expression(italic('linc.mettl3')),
                   expression(italic('linc.oncecut')),
                   
                   'controls exp4',
                   expression(italic('linc.mipep1')),
                   
                   'controls exp5',
                   expression(italic('linc.wrb')))

# second define the grp order
grp_order <- c('WT.Cas9', 'libra.F0', 'linc.cd74', 'linc.foxp2', 'linc.loc485.F0',
               'Linc.112305.F0', 'Linc.239.F0', 'linc.mettl3.F0', 'linc.oncecut.F0',
               'linc.mipep1.ex3inj',
               'linc.wrb.F0')
eud2$grp <- factor(eud2$grp, levels=grp_order)
euds2$grp <- factor(euds2$grp, levels=grp_order)

# this is now the order in which we should give the colours
grp_cols = c('#697a87', '#fcb505', '#78ac63', '#cb2a20', '#982150', '#ec6354', '#417dcd', '#7e7ba4', '#db5072', '#23aa4c', '#e94537')

gg_eucli <- ggplot(euds2, aes(x=grp, y=mean, colour=grp)) +
  geom_linerange(aes(ymin=mean-sem, ymax=mean+sem), size=1, colour='black') +
  geom_point(aes(colour=grp), fill='black', shape=21, stroke=1.5) +
  geom_quasirandom(data=eud2, aes(x=grp, y=euclid, fill=grp), width=0.25, alpha=0.15, shape=21, size=3, stroke=0) +
  facet_grid(~exp, scales='free_x', space='free_x') +
  scale_color_manual(values=grp_cols) +
  scale_fill_manual(values=grp_cols) +
  theme_minimal() +
  theme(
    panel.spacing=unit(12, 'pt'),
    strip.text.x=element_blank(), # this removes facet titles
    panel.grid.minor.y=element_blank(),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=9),
    #axis.text.x=element_text(size=7, angle=45, hjust=1, vjust=1),
    axis.text.x=element_blank(),
    legend.position='none'
  ) +
  coord_cartesian(ylim=c(0, 11.2)) +
  ylab('Euclidean distance from controls')
gg_eucli

ggsave(here('euclidean.pdf'), gg_eucli, width=110, height=75, units='mm')


# PCA ---------------------------------------------------------------------

upe_pca <- prcomp(upew[,4:ncol(upew)], scale=FALSE)

var_ex <- upe_pca$sdev / sum(upe_pca$sdev)
var_ex <- as.data.frame(cbind(sprintf('pc%i', 1:length(var_ex)), var_ex))
colnames(var_ex) <- c('pc', 'var_explained')

var_ex$pc <- factor(var_ex$pc, levels=sprintf('pc%i', 1:nrow(var_ex)))
var_ex$var_explained <- as.numeric(var_ex$var_explained)

# scree plot
ggplot(var_ex, aes(x=pc, y=var_explained)) +
  geom_point()


library(ggfortify)

autoplot(upe_pca,
         data=upew,
         colour='exp_grp')

autoplot(kmeans(upew[, 4:ncol(upew)], 10),
         data=upew)

kmeans(upew[, 4:ncol(upew)], 3)