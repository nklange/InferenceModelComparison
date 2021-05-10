library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("ggforce")
library("purrr")
library("gridExtra")
library("grid")
library(magick)
library(cowplot)
library(ungeviz)


options(pillar.sigfig = 5)


# Load data --------------------------------------------------------------------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],
                            data_vdb14,data_pratte17)

models <- c("MK_RNminus","MK_RNplus","J_RNminus","J_RNplus")

# Load Fits --------------------------------------------------------------------

FittedFull <- readRDS("Fits/FitFull.rds") %>% filter(model %in% models)

Fitted10Fold <- readRDS("Fits/CV10Fold_TestLL.rds") %>% filter(model %in% models)

Fitted5x2 <- readRDS("Fits/CV5x2_TestLL.rds") %>% filter(model %in% models)

FittedSs <- readRDS("Fits/CVLOSzO_TestLL.rds") %>% filter(model %in% models)

# Rename models -----------

oldmodellabels <- c("J_RNminus","MK_RNminus","J_RNplus","MK_RNplus")
newmodellabels <- c("VP(J)A-","VP(K)A-","VP(J)A+","VP(K)A+")
manuscriptmodelorder <- c("VP(J)A-","VP(K)A-","VP(J)A+","VP(K)A+")

# Exclude Participants -------------

minDelta <- FittedFull %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  mutate(mindiffLL = outer(LL,LL, `-`)[2,1],
         bestLL = LL[[1]]) %>%
  select(model,exp,cvid,mindiffLL,bestLL) %>%
  distinct()

# Tolerance 0.5LL = 1 Deviance point
excludeDeltaLL <- unique(minDelta %>%
                           distinct() %>%
                           filter(mindiffLL > 0.5) %>%
                           .$cvid)


Failed5x2Prediction <- unique(Fitted5x2 %>% filter(CVLL > 1e4) %>% .$id)
length(Failed5x2Prediction)
Failed10FPrediction <- unique(Fitted10Fold %>% filter(CVLL > 1e4) %>% .$id)
length(Failed10FPrediction)
FailedSzPrediction <- unique(FittedSs %>% filter(CVLL > 1e4) %>% .$id)
length(FailedSzPrediction)

# global exclusion:

# Instability + Isues in crossvalidation
  exclall <- unique(c(excludeDeltaLL,
               Failed10FPrediction,Failed5x2Prediction,FailedSzPrediction))
# Issues in crossvalidation only
exclall <- unique(c(Failed10FPrediction,Failed5x2Prediction,FailedSzPrediction))
#do not exclude globally
#exclall <- c()

195 - length(exclall)

allsetsizes <- unique(FittedFull$exp) # include all experiments

# include only experiments with Setsize 1 - 8
#allsetsizes <- c("VSCGM21a","VSCGM21b","VSCGM21c","ol17_e1")


# Identify best fit for all models/approaches  ---------------------------------

ntrials <- experimentfile %>%  filter(exp %in% allsetsizes) %>%
  group_by(id) %>% count() %>%
  slice(rep(1:n(),each=length(models))) %>%
  filter(!id %in% c(exclall)) %>%

  arrange(id)

BestFull <- FittedFull %>%
  filter(model %in% models) %>%
  filter(!cvid %in% c(exclall)) %>%

  filter(exp %in% allsetsizes) %>%

  group_by(model,cvid) %>%
  filter(Dev == min(Dev, na.rm = TRUE)) %>% #finds best fit of 20 runs
  arrange(cvid) %>%
  bind_cols(trials = as.numeric(ntrials$n)) %>% # add number of trials
  mutate(BIC = log(trials) * nparams + 2 * Dev) %>%
  distinct() %>%
  group_by(cvid) %>%
  mutate(delLL = LL - min(LL),
         delDev = Dev - min(Dev),
         delAIC = AIC - min(AIC),
         delBIC = BIC - min(BIC))

BestFull$model <- mapvalues(BestFull$model,
                            from = oldmodellabels,
                            to = newmodellabels)

## globally determines order of bars in graphs

modelsinorder <- manuscriptmodelorder
# modelsinorder <- BestFull %>%
#   group_by(model) %>%
#   summarize(means = mean(AIC)) %>%
#   arrange(-means) %>% .$model

BestFull$model <- factor(BestFull$model,levels = modelsinorder)


BestCV10Fold <- Fitted10Fold %>%
  filter(model %in% models) %>%
  filter(!id %in% exclall) %>%
  filter(!id %in% Failed10FPrediction) %>%
  filter(exp %in% allsetsizes) %>%

  distinct() %>%
  mutate(holdoutCVDev = CVLL * 2) %>%
  group_by(model,exp,id) %>%
  summarize(CVDev = sum(holdoutCVDev)) %>%  # sum over holdout sets
  ungroup() %>%
  mutate(model = mapvalues(model,from = oldmodellabels, to = newmodellabels)) %>%
  mutate(model = factor(model,modelsinorder))

BestCV5x2 <- Fitted5x2 %>%
  filter(model %in% models) %>%
  filter(!id %in% exclall) %>%
  filter(!id %in% Failed5x2Prediction) %>%
  filter(exp %in% allsetsizes) %>%

  mutate(CVDev = CVLL * 2) %>%
  mutate(leftoutrep = substr(leftout, start = 1, stop = 1)) %>%
  group_by(model,exp,id,leftoutrep) %>%
  summarize(CVDevsz = sum(CVDev)) %>%
  group_by(model,exp,id) %>%
  summarize(CVDev = mean(CVDevsz)) %>%
  ungroup() %>%
  mutate(model = mapvalues(model,from = oldmodellabels, to = newmodellabels)) %>%
  mutate(model = factor(model,modelsinorder))

BestCVSs <- FittedSs %>%
  filter(model %in% models) %>%
  filter(!id %in% exclall) %>%
  filter(!id %in% FailedSzPrediction) %>%
  filter(exp %in% allsetsizes) %>%

  #filter(leftout != 1) %>%  # == 1 only setsize 1; != bar setsize 1; commented out = all set sizes
  mutate(holdoutCVDev = CVLL * 2) %>%
  group_by(model,exp,id) %>%
  summarize(CVDev = sum(holdoutCVDev)) %>%
  ungroup() %>%
  mutate(model = mapvalues(model,from = oldmodellabels, to = newmodellabels)) %>%
  mutate(model = factor(model,modelsinorder))

# Set globally winning model as reference AIC/Deviance and calculate deltas ----

# The same happens for all approaches with possible difference in exclusions
# this could probably be in a neat function rather than repeating code...

AICmeans <- BestFull %>%
  group_by(model) %>%
  summarize(means = mean(AIC)) %>%
  arrange(means)

#sum all aics
AICtot <- BestFull %>%
  group_by(cvid,model) %>%
  mutate(deltaAIC = AIC - AICmeans %>% slice(1) %>% .$means) %>%
  group_by(model) %>%
  summarize(mdeltaAIC = mean(deltaAIC)) %>%
  mutate(mdeltaAICadj = ifelse(mdeltaAIC > 20,20,mdeltaAIC)) %>%
  arrange(model)

WinningFull1 <- BestFull %>%
  filter(delAIC == 0) %>% # only 1 model can win
  group_by(exp,model) %>%
  summarize(numBest = length(delAIC)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

WinningFull2 <- BestFull %>%
  filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
  group_by(exp,model) %>%
  summarize(numBest = length(delAIC)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

AICtot$BestModelstrict <- WinningFull1$wu
AICtot$BestModellenient<- WinningFull2$wu



BICmeans <- BestFull %>%
  group_by(model) %>%
  summarize(means = mean(BIC)) %>%
  arrange(means)

#sum all aics
BICtot <- BestFull %>%
  group_by(cvid,model) %>%
  mutate(deltaBIC = BIC - BICmeans %>% slice(1) %>% .$means) %>%
  group_by(model) %>%
  summarize(mdeltaBIC = mean(deltaBIC)) %>%
  mutate(mdeltaBICadj = ifelse(mdeltaBIC > 20,20,mdeltaBIC)) %>%
  arrange(model)

WinningFull1 <- BestFull %>%
  #filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
  filter(delBIC == 0) %>% # only 1 model can win
  group_by(exp,model) %>%
  summarize(numBest = length(delAIC)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

WinningFull2 <- BestFull %>%
  filter(delBIC < 1) %>% # assumes equivalence of models within 1 Deviance point
  #filter(delAIC == 0) %>% # only 1 model can win
  group_by(exp,model) %>%
  summarize(numBest = length(delAIC)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

BICtot$BestModelstrict <- WinningFull1$wu
BICtot$BestModellenient<- WinningFull2$wu

## CV 10F

Dev10means <- BestCV10Fold %>%
  group_by(model) %>%
  summarize(means = mean(CVDev)) %>%
  arrange(means)

CV10Foldtot <- BestCV10Fold %>%
  group_by(id,model) %>%
  mutate(deltaDev = CVDev - Dev10means %>% slice(1) %>% .$means) %>%
  group_by(model) %>%
  summarize(mdeltaDev = mean(deltaDev)) %>%
  arrange(model) %>%
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 20,20,mdeltaDev))

Winning10F1 <-  BestCV10Fold %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  #filter(delDevCV < 1) %>%
  filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

Winning10F2 <-  BestCV10Fold %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  filter(delDevCV < 1) %>%
  #filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

CV10Foldtot$BestModelstrict <- Winning10F1$wu
CV10Foldtot$BestModellenient <- Winning10F2$wu


## CV 5x2F

Dev5x2means <- BestCV5x2 %>%
  group_by(model) %>%
  summarize(means = mean(CVDev)) %>%
  arrange(means)

CV5x2Foldtot <- BestCV5x2 %>%
  group_by(id,model) %>%
  mutate(deltaDev = CVDev - Dev5x2means %>% slice(1) %>% .$means) %>%
  arrange(id) %>%
  group_by(model) %>%
  summarize(mdeltaDev = mean(deltaDev)) %>%
  arrange(model) %>%
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 20,20,mdeltaDev))


Winning5x2F1 <- BestCV5x2 %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  #filter(delDevCV < 1) %>%
  filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

Winning5x2F2 <- BestCV5x2 %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  filter(delDevCV < 1) %>%
  #filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/(195 - length(exclall))) %>%
  arrange(model)

CV5x2Foldtot$BestModelstrict <- Winning5x2F1$wu
CV5x2Foldtot$BestModellenient <- Winning5x2F2$wu

## CV LOSsO

DevSsmeans <- BestCVSs %>%
  group_by(model) %>%
  summarize(means = mean(CVDev)) %>%
  arrange(means)

CVSstot <- BestCVSs %>%
  group_by(id,model) %>%
  mutate(deltaDev = CVDev - DevSsmeans %>% slice(1) %>% .$means) %>%
  group_by(model) %>%
  summarize(mdeltaDev = mean(deltaDev)) %>%
  arrange(model) %>%
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 20,20,mdeltaDev))

inclid <- length(unique(BestCVSs$id))
WinningSs1 <- BestCVSs %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  #filter(delDevCV < 1) %>%
  filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/inclid) %>%
  arrange(model)

WinningSs2 <- BestCVSs %>%
  group_by(id) %>%
  mutate(delDevCV = CVDev - min(CVDev)) %>%
  filter(delDevCV < 1) %>%
  #filter(delDevCV == 0) %>%
  group_by(model) %>%
  summarize(numBest = length(delDevCV)) %>%
  complete(model, fill = list(0)) %>%
  group_by(model) %>%
  summarize(wu = sum(numBest,na.rm=T)) %>%
  mutate(wu = wu/inclid) %>%
  arrange(model)

CVSstot$BestModelstrict <- WinningSs1$wu
CVSstot$BestModellenient <- WinningSs2$wu

# Response Noise, Parameterization of memory precision -------------------------

Total <- bind_rows(AICtot %>% mutate(method = "AIC"),
                   BICtot %>% mutate(method = "BIC"),
                   CV10Foldtot %>% mutate(method = "10Fold-CV"),
                   CV5x2Foldtot %>% mutate(method = "5x2-CV"),
                   CVSstot %>% mutate(method = "LOSsO-CV")) %>%
  mutate(mdeltaDev = ifelse(method == "AIC",mdeltaAIC,
                            ifelse(method == "BIC",mdeltaBIC,mdeltaDev)),
         mdeltaDevadj = ifelse(method == "AIC",mdeltaAICadj,
                               ifelse(method == "BIC",mdeltaBICadj,
                                      mdeltaDevadj))) %>%
  mutate(labeldelta = ifelse(mdeltaDev > 20,mdeltaDev,NA))

Total$method <- factor(Total$method,
                       levels = c("AIC","BIC","10Fold-CV","5x2-CV","LOSsO-CV"),
                       labels = c("AIC","BIC","10-F-CV","5x2-F-CV","LOSsO-CV"))


ylim.prim <- c(0,20)
ylim.sec <- c(0,1)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])


ModelCompTotal <- ggplot(Total %>% filter(!method == "BIC"),
       aes(y = mdeltaDevadj,x = model,fill = mdeltaDevadj )) +
  geom_rect(aes(xmin=0, xmax=4.5, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
  geom_point(aes(y = a + BestModelstrict*b),
             color = "black",size=4,shape=21,fill="transparent") +
  geom_point(aes(y = a + BestModellenient*b),
             color = "black",size=4) +
  coord_flip(ylim = c(0,20)) +
  scale_fill_gradientn(name = expression(paste(Delta)),
                       colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
                       limits=c(0,20))+
  geom_text(aes(label=round(labeldelta), y = 18.5), hjust = -0.5,size=3)+
  scale_y_continuous(name = expression(paste(Delta," AIC/Deviance")),
                     breaks = seq(0,20,4),
                     labels = c(seq(0,16,4),"..."),
                     sec.axis = sec_axis(~ (. - a)/b,
                                         name = "Proportion datasets best fit",
                                         labels = c("0",".25",".50",".75","1"))) +

  scale_x_discrete(name = "Model", labels = c(
    expression(paste("VP(",italic(J),")A-")),
    expression(paste("VP(",kappa,")A-")),
    expression(paste("VP(",italic(J),")A+")),
    expression(paste("VP(",kappa,")A+"))
  )) +
  facet_wrap(.~method)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


# LOSsO: Response Noise, Parameterization of memory precision ------------------


TotsSs <- NULL
for (Ss in c(1:8)){

  BestCVSs <- FittedSs %>%
    filter(model %in% models) %>%
    filter(exp %in% allsetsizes) %>%
    #filter(!id %in% exclall) %>%
    filter(!id %in% FailedSzPrediction) %>%
    filter(leftout == Ss) %>%
    mutate(holdoutCVDev = CVLL * 2) %>%
    group_by(model,exp,id) %>%
    summarize(CVDev = sum(holdoutCVDev)) %>%
    ungroup() %>%
    mutate(model = mapvalues(model,
                             from = oldmodellabels,
                             to = newmodellabels)) %>%
    mutate(model = factor(model,modelsinorder))

  DevSsmeans <- BestCVSs %>%
    group_by(model) %>%
    summarize(means = mean(CVDev)) %>%
    arrange(means)

  tots <- length(unique(BestCVSs$id))
  CVSstot <- BestCVSs %>%
    group_by(id,model) %>%
    mutate(deltaDev = CVDev - DevSsmeans %>% slice(1) %>% .$means) %>%
    group_by(model) %>%
    summarize(mdeltaDev = mean(deltaDev)) %>%
    arrange(model) %>%
    mutate(mdeltaDevadj = ifelse(mdeltaDev > 20,20,mdeltaDev))  %>%
    mutate(leftout = Ss,
           totalppl = tots)


  WinningSs1 <- BestCVSs %>%
    group_by(id) %>%
    mutate(delDevCV = CVDev - min(CVDev)) %>%
    #filter(delDevCV < 1) %>%
    filter(delDevCV == 0) %>%
    group_by(model) %>%
    summarize(numBest = length(delDevCV)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/tots) %>%
    arrange(model)

  WinningSs2 <- BestCVSs %>%
    group_by(id) %>%
    mutate(delDevCV = CVDev - min(CVDev)) %>%
    filter(delDevCV < 1) %>%
    #filter(delDevCV == 0) %>%
    group_by(model) %>%
    summarize(numBest = length(delDevCV)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/tots) %>%
    arrange(model)

  CVSstot$BestModelstrict <- WinningSs1$wu
  CVSstot$BestModellenient <- WinningSs2$wu

  TotsSs <- bind_rows(TotsSs,CVSstot)
}

TotsSs$leftout <- factor(TotsSs$leftout, levels = c(1:8),
                         labels = c("Left out: Ss 1",
                                    "Left out: Ss 2",
                                    "Left out: Ss 3",
                                    "Left out: Ss 4",
                                    "Left out: Ss 5",
                                    "Left out: Ss 6",
                                    "Left out: Ss 7",
                                    "Left out: Ss 8"))

Ppl <- TotsSs %>% group_by(leftout) %>% slice(1)
TotsSs <- TotsSs %>% mutate(labeldelta = ifelse(mdeltaDev > 20,mdeltaDev,NA))


ModelCompSs<- ggplot(TotsSs,aes(y = mdeltaDevadj,x = model,fill = mdeltaDevadj )) +
  geom_rect(aes(xmin=0, xmax=4.5, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
  geom_point(aes(y = a + BestModelstrict*b),
             color = "black",size=3,shape=21,fill="transparent") +
  geom_point(aes(y = a + BestModellenient*b),color = "black",size=4) +
  coord_flip(ylim = c(0,20),xlim=c(0,5)) +
  scale_fill_gradientn(name = expression(paste(Delta," Dev (out-of-sample)")),
                       colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
                       limits=c(0,20))+
  geom_text(aes(y = 19,label = round(labeldelta)),size = 3)+
  geom_text(data = Ppl, aes(x = 5, y = 17,label = paste("N =",totalppl)),size=3.5)+
  #geom_text(aes(label=round(mdeltaDev,1)), hjust = -0.5,size=3.5)+
 # scale_x_discrete(name = "Model") +
  scale_y_continuous(name = expression(paste(Delta," Dev (out-of-sample)")),
                     breaks = seq(0,20,4),
                     labels = c(seq(0,16,4),"..."),
                     sec.axis = sec_axis(~ (. - a)/b,
                                         name = "Proportion datasets best fit",
                                         labels = c("0",".25",".50",".75","1"))) +

  scale_x_discrete(name = "Model", labels = c(
    expression(paste("VP(",italic(J),")A-")),
    expression(paste("VP(",kappa,")A-")),
    expression(paste("VP(",italic(J),")A+")),
    expression(paste("VP(",kappa,")A+"))
  )) +
  facet_wrap(.~leftout)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        panel.border = element_rect(colour = "black",fill = "transparent"))


# # Graphs for supplementary material: -------------------------------------------
#
# # Excluding unstable participants
# ModelCompTotalAll <- ModelCompTotal #using exclAll for instability
# ModelCompTotalLim <- ModelCompTotal #using only experiments with limited set sizes
# # LOSSO when only including experiments with all set sizes
#
# ModelCompSsAll <- ModelCompSs #using exclAll for instability
# ModelCompSsLim <- ModelCompSs #using only experiments with limited set sizes
#
# # By experiment
#
# experimentfile$exp <- as.character(experimentfile$exp)
# pptperexp <- experimentfile %>%
#   group_by(exp) %>%
#   summarize(ntot = length(unique(id))) %>%
#   arrange(exp) %>%
#   slice(rep(1:n(),each=length(models)))
#
#
#
# AICmeans <- BestFull %>%
#   group_by(exp,model) %>%
#   summarize(means = mean(AIC),
#             lengthppt = length(AIC)) %>%
#   arrange(exp,means) %>%
#   group_by(exp) %>% slice(1)
#
# pptperexp <- pptperexp %>%
#   filter(exp %in% unique(AICmeans$exp))
#
# BestAIC = rep(AICmeans$means, AICmeans$lengthppt*length(models))
#
# AICbyExp <- BestFull
# AICbyExp <- AICbyExp %>% ungroup() %>%
#   arrange(exp) %>%
#   mutate(BestAIC = BestAIC) %>%
#   mutate(delAIC = AIC - BestAIC) %>%
#   group_by(exp,model) %>%
#   summarize(deltaAIC = mean(delAIC),
#             nppt = length(delAIC))
# AICbyExp$ntot <- pptperexp$ntot
#
#
# Dev5x2means <- BestCV5x2 %>%
#   group_by(exp,model) %>%
#   summarize(means = mean(CVDev),
#             lengthppt = length(CVDev)) %>%
#   arrange(exp,means) %>%
#   group_by(exp) %>%
#   slice(1)
#
# BestDev5x2 = rep(Dev5x2means$means, Dev5x2means$lengthppt*length(models))
# Dev5x2byExp <- BestCV5x2
# Dev5x2byExp <- Dev5x2byExp %>% ungroup() %>%
#   arrange(exp) %>%
#   mutate(BestDev = BestDev5x2) %>%
#   mutate(delDev = CVDev - BestDev) %>%
#   group_by(exp,model) %>%
#   summarize(deltaDev = mean(delDev),
#             nppt = length(delDev))
# Dev5x2byExp$ntot <- pptperexp$ntot
#
#
# Dev10Fmeans <- BestCV10Fold %>% group_by(exp,model) %>% summarize(means = mean(CVDev),
#                                                                   lengthppt = length(CVDev)) %>% arrange(exp,means) %>% group_by(exp) %>% slice(1)
#
# BestDev10F = rep(Dev10Fmeans$means, Dev10Fmeans$lengthppt*length(models))
# Dev10FbyExp <- BestCV10Fold
# Dev10FbyExp <- Dev10FbyExp %>% ungroup() %>%
#   arrange(exp) %>%
#   mutate(BestDev = BestDev10F) %>%
#   mutate(delDev = CVDev - BestDev) %>%
#   group_by(exp,model) %>%
#   summarize(deltaDev = mean(delDev),
#             nppt = length(delDev))
# Dev10FbyExp$ntot <- pptperexp$ntot
#
#
# DevSsmeans <- BestCVSs %>% group_by(exp,model) %>% summarize(means = mean(CVDev),
#                                                              lengthppt = length(CVDev)) %>% arrange(exp,means) %>% group_by(exp) %>% slice(1)
#
# BestDevSs = rep(DevSsmeans$means, DevSsmeans$lengthppt*length(models))
# DevSsbyExp <- BestCVSs
# DevSsbyExp <- DevSsbyExp %>% ungroup() %>%
#   arrange(exp) %>%
#   mutate(BestDev = BestDevSs) %>%
#   mutate(delDev = CVDev - BestDev) %>%
#   group_by(exp,model) %>%
#   summarize(deltaDev = mean(delDev),
#             nppt = length(delDev))
# DevSsbyExp$ntot <- pptperexp$ntot
#
#
# Total <- bind_rows(AICbyExp %>% mutate(method = "AIC"),
#                    Dev10FbyExp  %>% mutate(method = "10-F-CV"),
#                    Dev5x2byExp %>% mutate(method = "5x2-F-CV"),
#                    DevSsbyExp  %>% mutate(method = "LOSsO-CV")) %>%
#   mutate(model = factor(model, levels = modelsinorder))
# Total <- Total %>% mutate(deltaDev = ifelse(method == "AIC",deltaAIC,deltaDev),
#                           deltadjDev = ifelse(deltaDev > 20,20,deltaDev))
# Total$method <- as.character(Total$method)
# Total$method <- factor(Total$method,
#                        levels = c("AIC","10-F-CV","5x2-F-CV","LOSsO-CV"))
#
# Nleft <- AICbyExp %>%
#   group_by(exp) %>%
#   slice(1) %>%
#   slice(rep(1:n(),4))
# Nleft$method <- rep(c("AIC","10-F-CV","5x2-F-CV","LOSsO-CV"),length(unique(AICbyExp$exp)))
# Nleft$method <- factor(Nleft$method,
#                        levels = c("AIC","10-F-CV","5x2-F-CV","LOSsO-CV"))
#
#
# ExpG <- ggplot(data = Total,aes(y = model,x = exp)) +
#   geom_tile(data = Total, aes(fill = deltadjDev)) +
#
#   geom_text(data = Nleft,aes(y = 5.2, x = exp,
#                              label = paste0(nppt,"/",ntot)),vjust=1.5,size=3) +
#   geom_text(aes(label = round(deltaDev)),size = 3.5) +
#   scale_fill_gradientn(name = expression(paste(Delta)),
#                        colors = rev(brewer.pal(n = 9, name = "YlOrRd")))+
#   scale_x_discrete(position = "top",name = "Experiment") +
#
#
#   scale_y_discrete(name = "Model", labels = c(
#     expression(paste("VP(",italic(J),")A-")),
#     expression(paste("VP(",kappa,")A-")),
#     expression(paste("VP(",italic(J),")A+")),
#     expression(paste("VP(",kappa,")A+"))
#
#
#
# )) +
#
#   # below: only exclude for extreme set sizes
#   # scale_y_discrete(name = "Model", labels = c(
#   #   expression(paste("VP(",kappa,")U-")),
#   #   expression(paste("VP(",kappa,")A-")),
#   #   expression(paste("VP(",kappa,")P-")),
#   #   expression(paste("VP(",kappa,")F+")),
#   #   expression(paste("VP(",kappa,")F-")),
#   #   expression(paste("VP(",kappa,")A+")),
#   #   expression(paste("VP(",kappa,")P+")),
#   #   expression(paste("VP(",kappa,")U+"))
#   # )) +
#   facet_wrap(.~method) +
#   theme(axis.text.x = element_text(angle = 90,size=10),
#         axis.title = element_blank(),
#         axis.text.y = element_text(size=12),
#         axis.ticks = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         panel.border = element_rect(fill = "transparent",colour="black"),
#         strip.placement = "outside",
#         legend.position = "none",
#         strip.background = element_rect(fill = "transparent",
#                                         colour="transparent"),
#         strip.text = element_text(size = 12))
#
# ByExp <- plot_grid(ExpG,labels=c("E"))
#
# Others <- plot_grid(ModelCompTotalAll,ModelCompSsAll,ModelCompTotalLim,ModelCompSsLim,
#                    labels=c("A","B","C","D"),nrow=2,
#                    rel_widths = c(1,1),align = 'h', axis = 't',scale = 0.98)
#
# plot_grid(Others,ByExp,nrow=2,rel_heights = c(2,0.9))
#
# ggsave("SuppRQ1.png", units="cm", width=30, height=35, dpi=600)
