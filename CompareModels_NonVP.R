library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("ggforce")
library("purrr")

options(pillar.sigfig = 5)


# Load data --------------------------------------------------------------------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],
                            data_vdb14,data_pratte17)


models <- c("MK_RNminus",
            "MK_FM_RNplus","MK_RNplus","MK_U_RNplus","MK_P_RNplus",
            "MK_P_RNminus","MK_U_RNminus",
            "MK_FM_RNminus",
            
            "SA_F_RNplus","SA_F_RNminus",
            "SA_P_RNplus", "SA_U_RNplus","EP_P_RNplus","EP_U_RNplus",
            "EP_RNplus","SA_RNplus",
            "EP_RNminus", "EP_FM_RNplus","EP_FM_RNminus",
            "EP_P_RNminus","EP_U_RNminus","SA_P_RNminus","SA_U_RNminus")


# Load Fits --------------------------------------------------------------------

FittedFull <- readRDS("Fits/FitFull.rds") %>% filter(model %in% models)

Fitted10Fold <- readRDS("Fits/CV10Fold_TestLL.rds") %>% filter(model %in% models)

Fitted5x2 <- readRDS("Fits/CV5x2_TestLL.rds") %>% filter(model %in% models)

FittedSs <- readRDS("Fits/CVLOSzO_TestLL.rds") %>% filter(model %in% models)

# Rename models -----------


oldmodellabels <- c("MK_RNminus",
                    "MK_FM_RNplus","MK_RNplus","MK_U_RNplus","MK_P_RNplus",
                    "MK_P_RNminus","MK_U_RNminus",
                    "MK_FM_RNminus",
                    
                    "SA_F_RNplus","SA_F_RNminus",
                    "SA_P_RNplus", "SA_U_RNplus","EP_P_RNplus","EP_U_RNplus",
                    "EP_RNplus","SA_RNplus",
                    "EP_RNminus", "EP_FM_RNplus","EP_FM_RNminus",
                    "EP_P_RNminus","EP_U_RNminus","SA_P_RNminus","SA_U_RNminus")
                   
newmodellabels <- c("VP(K)A-",
                     "VP(K)FM+","VP(K)A+","VP(K)U+","VP(K)P+",
                     "VP(K)P-","VP(K)U-",
                     "VP(K)FM-",
                     
                     "SA(K)F+","SA(K)F-",
                     "SA(K)P+","SA(K)U+","EP(K)P+","EP(K)U+",
                     
                     "EP(K)A+","SA(K)A+",
                     "EP(K)A-","EP(K)FM+","EP(K)FM-",
                     "EP(K)P-","EP(K)U-","SA(K)P-","SA(K)U-")


manuscriptmodelorder <- c("SA(K)A+","EP(K)A+","EP(K)A-","SA(K)U+","SA(K)U-",
                           "EP(K)FM+","EP(K)FM-","EP(K)U-","SA(K)F-","SA(K)F+","EP(K)U+", 
                           "SA(K)P+","SA(K)P-","VP(K)U-","EP(K)P-","VP(K)A-","EP(K)P+",
                           "VP(K)P-","VP(K)FM+","VP(K)FM-","VP(K)A+","VP(K)P+","VP(K)U+")
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

# Isues in crossvalidation + AIC fit unstable (tolerance 1 deviance point)
#exclall <- unique(c(excludeDeltaLL,Failed10FPrediction,Failed5x2Prediction,FailedSzPrediction))
# Issues in crossvalidation only
exclall <- unique(c(Failed10FPrediction,Failed5x2Prediction,FailedSzPrediction))
#do not exclude globally
#exclall <- c()

195 - length(exclall)

allsetsizes <- unique(FittedSs$exp)

# uncomment below to only include data sets that contain data for all set sizes 1 - 8
#allsetsizes <- c("VSCGM21a","VSCGM21b","VSCGM21c","ol17_e1")

# Identify best fit for all models/approaches  ---------------------------------

ntrials <- experimentfile %>% group_by(id) %>%  filter(exp %in% allsetsizes) %>% 
 count() %>% 
  slice(rep(1:n(),each=length(models))) %>% filter(!id %in% c(exclall)) %>% 
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
  mutate(mdeltaAICadj = ifelse(mdeltaAIC > 50,50,mdeltaAIC)) %>% 
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
  mutate(mdeltaBICadj = ifelse(mdeltaBIC > 50,50,mdeltaBIC)) %>% 
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
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 50,50,mdeltaDev))

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
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 50,50,mdeltaDev))


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
  mutate(mdeltaDevadj = ifelse(mdeltaDev > 50,50,mdeltaDev))

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
  mutate(labeldelta = ifelse(mdeltaDev > 50,mdeltaDev,NA)) %>% 
  mutate(modeldist = case_when(grepl("VP", model, fixed = TRUE) ~ "VP",
                               grepl("EP",model,fixed=TRUE) ~ "EP",
                               TRUE ~ "SA"))

Total$method <- factor(Total$method,
                       levels = c("AIC","BIC","10Fold-CV","5x2-CV","LOSsO-CV"),
                       labels = c("AIC","BIC","10-F-CV","5x2-F-CV","LOSsO-CV"))

ylim.prim <- c(0,50)
ylim.sec <- c(0,1)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])



ModelCompTotal <- ggplot(Total %>% filter(!method == "BIC"),
       aes(y = mdeltaDevadj,x = model,fill = modeldist )) +
  geom_rect(aes(xmin=0, xmax=24, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
  geom_point(aes(y = a + BestModelstrict*b),
             color = "black",size=2,shape=21,fill="transparent") +
  geom_point(aes(y = a + BestModellenient*b),
             color = "black",size=2) +
  coord_flip(ylim = c(0,55)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  # scale_fill_gradientn(name = expression(paste(Delta)),
  #                      colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
  #                      limits=c(0,50))+
  geom_text(aes(label=round(labeldelta), y = 52), hjust = 0.5,size=3)+
  scale_y_continuous(name = expression(paste(Delta," AIC/Deviance")),
                     breaks = seq(0,50,10), 
                     labels = c(seq(0,40,10),"..."),
                     sec.axis = sec_axis(~ (. - a)/b,
                                         name = "Proportion datasets best fit",
                                         labels = c("0",".25",".50",".75","1"))) +

  scale_x_discrete(name = "Model", labels = c(
    expression(paste("SA(",kappa,")A+")),
    expression(paste("EP(",kappa,")A+")),
    expression(paste("EP(",kappa,")A-")),
    expression(paste("SA(",kappa,")U+")),
    expression(paste("SA(",kappa,")U-")),
    expression(paste("EP(",kappa,")F+")),
   
    expression(paste("EP(",kappa,")F-")),
    expression(paste("EP(",kappa,")U-")),
    expression(paste("SA(",kappa,")F-")),

    expression(paste("SA(",kappa,")F+")),
    expression(paste("EP(",kappa,")U+")),
    expression(paste("SA(",kappa,")P+")),
    expression(paste("SA(",kappa,")P-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("EP(",kappa,")P-")),
    expression(paste("VP(",kappa,")A-")),
    expression(paste("EP(",kappa,")P+")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")A+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+"))
  )) +
  facet_wrap(.~method)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12,hjust=0),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
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
    filter(!id %in% exclall) %>%
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
    mutate(mdeltaDevadj = ifelse(mdeltaDev > 50,50,mdeltaDev))  %>% 
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

Ppl <- TotsSs %>% group_by(leftout) %>% slice(1) %>% mutate(modeldist = NA)
TotsSs <- TotsSs %>% mutate(labeldelta = ifelse(mdeltaDev > 50,mdeltaDev,NA)) %>% 
  mutate(modeldist = case_when(grepl("VP", model, fixed = TRUE) ~ "VP",
                               grepl("EP",model,fixed=TRUE) ~ "EP",
                               TRUE ~ "SA"))


ModelCompSs<-ggplot(TotsSs,aes(y = mdeltaDevadj,x = model,fill = modeldist)) +
  geom_rect(aes(xmin=0, xmax=25, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
  geom_point(aes(y = a + BestModelstrict*b),
             color = "black",size=3,shape=21,fill="transparent") +
  geom_point(aes(y = a + BestModellenient*b),color = "black",size=3) +
  coord_flip(ylim = c(0,50),xlim=c(0,25)) +
  # scale_fill_gradientn(name = expression(paste(Delta," Dev")),
  #                      colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
  #                      limits=c(0,50),guide = "colorbar")+
  geom_text(aes(y = 45,label = round(labeldelta)),size = 3)+
  geom_text(data = Ppl, aes(x = 25, y = 24,label = paste("N = ",totalppl)))+
  #geom_text(aes(label=round(mdeltaDev,1)), hjust = -0.5,size=3.5)+
 # scale_x_discrete(name = "Model") +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  
  scale_y_continuous(name = expression(paste(Delta," Dev (out-of-sample)")),
                     breaks = seq(0,50,10), 
                     labels = c(seq(0,40,10),"..."),
                     sec.axis = sec_axis(~ (. - a)/b,
                                         name = "Proportion datasets best fit",
                                         labels = c("0",".25",".50",".75","1"))) +
  scale_x_discrete(name = "Model", labels = c(
    expression(paste("SA(",kappa,")A+")),
    expression(paste("EP(",kappa,")A+")),
    expression(paste("EP(",kappa,")A-")),
    expression(paste("SA(",kappa,")U+")),
    expression(paste("SA(",kappa,")U-")),
    expression(paste("EP(",kappa,")F+")),
    
    expression(paste("EP(",kappa,")F-")),
    expression(paste("EP(",kappa,")U-")),
    expression(paste("SA(",kappa,")F-")),
    
    expression(paste("SA(",kappa,")F+")),
    expression(paste("EP(",kappa,")U+")),
    expression(paste("SA(",kappa,")P+")),
    expression(paste("SA(",kappa,")P-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("EP(",kappa,")P-")),
    expression(paste("VP(",kappa,")A-")),
    expression(paste("EP(",kappa,")P+")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")A+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+"))
  )) +
  facet_wrap(.~leftout)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=9,hjust=0),
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
# ModelCompTotalAll <- ModelCompTotal
# ModelCompTotalLim <- ModelCompTotal
# # LOSSO when only including experiments with all set sizes
# 
# ModelCompSsAll <- ModelCompSs #using exclAll for instability
# ModelCompSsLim <- ModelCompSs #using only experiments with limited set sizes
# 
# # By experiment 

experimentfile$exp <- as.character(experimentfile$exp)
pptperexp <- experimentfile %>%
  group_by(exp) %>%
  summarize(ntot = length(unique(id))) %>%
  arrange(exp) %>%
  slice(rep(1:n(),each=length(models)))



AICmeans <- BestFull %>%
  group_by(exp,model) %>%
  summarize(means = mean(AIC),
            lengthppt = length(AIC)) %>%
  arrange(exp,means) %>%
  group_by(exp) %>% slice(1)

pptperexp <- pptperexp %>%
  filter(exp %in% unique(AICmeans$exp))

BestAIC = rep(AICmeans$means, AICmeans$lengthppt*length(models))

AICbyExp <- BestFull
AICbyExp <- AICbyExp %>% ungroup() %>%
  arrange(exp) %>%
  mutate(BestAIC = BestAIC) %>%
  mutate(delAIC = AIC - BestAIC) %>%
  group_by(exp,model) %>%
  summarize(deltaAIC = mean(delAIC),
            nppt = length(delAIC))
AICbyExp$ntot <- pptperexp$ntot



Nleft <- AICbyExp %>%
  group_by(exp) %>%
  slice(1) %>%
  slice(rep(1:n(),4)) %>% mutate(modeldist=NA)


Total <- AICbyExp %>%
  mutate(model = factor(model, levels = modelsinorder))
Total <- Total %>% mutate(deltaDev = deltaAIC,
                          deltadjDev = ifelse(deltaDev > 50,50,deltaDev)) %>% 
  mutate(modeldist = case_when(grepl("VP", model, fixed = TRUE) ~ "VP",
                               grepl("EP",model,fixed=TRUE) ~ "EP",
                               TRUE ~ "SA"))



ylim.prim <- c(0,50)
ylim.sec <- c(0,1)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])



ExpG <- ggplot(Total,
                         aes(y = deltadjDev,x = model,fill = modeldist )) +
  geom_rect(aes(xmin=0, xmax=24, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
    geom_text(data = Nleft,aes(x = 22, y = 45,
                               label = paste0(nppt,"/",ntot)),vjust=1.5,size=3) +
  # geom_point(aes(y = a + BestModelstrict*b),
  #            color = "black",size=2,shape=21,fill="transparent") +
  # geom_point(aes(y = a + BestModellenient*b),
  #            color = "black",size=2) +
  coord_flip(ylim = c(0,55)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+
  # scale_fill_gradientn(name = expression(paste(Delta)),
  #                      colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
  #                      limits=c(0,50))+
  geom_text(aes(label=ifelse(round(deltaDev) > 50,round(deltaDev),NA), y = 52), hjust = 0.5,size=3)+
  scale_y_continuous(name = expression(paste(Delta," AIC/Deviance")),
                     breaks = seq(0,50,10), 
                     labels = c(seq(0,40,10),"...")) +
  #scale_x_discrete() +
  
  #include all data sets
   scale_x_discrete(name = "Model", labels = c(
     expression(paste("SA(",kappa,")A+")),
     expression(paste("EP(",kappa,")A+")),
     expression(paste("EP(",kappa,")A-")),
     expression(paste("SA(",kappa,")U+")),
     expression(paste("SA(",kappa,")U-")),
     expression(paste("EP(",kappa,")F+")),
     
     expression(paste("EP(",kappa,")F-")),
     expression(paste("EP(",kappa,")U-")),
     expression(paste("SA(",kappa,")F-")),
     
     expression(paste("SA(",kappa,")F+")),
     expression(paste("EP(",kappa,")U+")),
     expression(paste("SA(",kappa,")P+")),
     expression(paste("SA(",kappa,")P-")),
     expression(paste("VP(",kappa,")U-")),
     expression(paste("EP(",kappa,")P-")),
     expression(paste("VP(",kappa,")A-")),
     expression(paste("EP(",kappa,")P+")),
     expression(paste("VP(",kappa,")P-")),
     expression(paste("VP(",kappa,")F+")),
     expression(paste("VP(",kappa,")F-")),
     expression(paste("VP(",kappa,")A+")),
     expression(paste("VP(",kappa,")P+")),
     expression(paste("VP(",kappa,")U+"))
)) +
facet_wrap(.~exp,nrow=4)+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=10,hjust=0),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=12))


# plot_grid(Others,ByExp,nrow=2,rel_heights = c(2,0.9))
# 
# ggsave("Test2.png", units="cm", width=30, height=40, dpi=600)
