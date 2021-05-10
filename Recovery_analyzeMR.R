#library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("ggforce")
library("purrr")
library("vwmvp")
library("stringr")
library("grid")
library("gridExtra")

# models ------------------------------

modellevels <- c("MK_RNminus","J_RNminus",
                 "MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
                 "MK_RNplus", "J_RNplus","MK_FM_RNplus",
                 "MK_P_RNplus","MK_U_RNplus")

modellabels <- c('"VP("~kappa~")A-"','"VP("~J~")A-"',
                 '"VP("~kappa~")F-"','"VP("~kappa~")P-"',
                 '"VP("~kappa~")U-"',
                 '"VP("~kappa~")A+"','"VP("~J~")A+"',
                 '"VP("~kappa~")F+"','"VP("~kappa~")P+"',
                 '"VP("~kappa~")U+"')

graphmodellabels <- c(
  expression(paste("VP(",kappa,")A-")),
  expression(paste("VP(",italic(J),")A-")),
  expression(paste("VP(",kappa,")F-")),
  expression(paste("VP(",kappa,")P-")),
  expression(paste("VP(",kappa,")U-")),
  expression(paste("VP(",kappa,")A+")),
  expression(paste("VP(",italic(J),")A+")),
  expression(paste("VP(",kappa,")F+")),
  expression(paste("VP(",kappa,")P+")),
  expression(paste("VP(",kappa,")U+"))
)


# Load files ---------------------------

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = TRUE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}

# Median parameter estimates ---------------------------------------------------

modelrecovery <- load_files("ModelRecovery/","rec")
parameterrecovery <- load_files("ParameterRecovery","trialsSs120")


modrec <- modelrecovery %>%
  group_by(genmodel,model,id) %>%
  arrange(objective) %>%
  slice(1) %>%
  group_by(model,genmodel) %>%
  mutate(nparam = length(vwmvp::get_start_vp(model)))
parrec <- parameterrecovery %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  slice(1) %>%
  group_by(model) %>%
  mutate(nparam = length(vwmvp::get_start_vp(model))) %>%
  mutate(genmodel = model)

# combine model and parameter recovery

fitted <- bind_rows(modrec,parrec) %>%
  arrange(genmodel) %>%
  mutate(AIC = 2 * objective + 2 * nparam) %>%
  group_by(genmodel,model)

mr <- fitted %>%
  summarize(avAIC = mean(AIC)) %>%
  mutate(delAIC = avAIC - min(avAIC)) %>%
  mutate(deladjAIC = ifelse(delAIC > 6,6,delAIC))



mr$model <- factor(mr$model,
                         levels = modellevels,
                         labels = modellabels)
mr$genmodel <- factor(mr$genmodel,
                   levels = modellevels,
                   labels = modellabels)

modeldevmedian <- ggplot(mr,aes(x = genmodel,y = model,fill = deladjAIC)) +
  geom_tile() +
 # geom_text(data = Nleft,aes(y = 11, x = exp,label = paste0(nppt,"/",ntot)),vjust=1.5,size=4) +
  geom_text(aes(label = round(delAIC)),size = 3.5) +
  scale_fill_gradientn(name = expression(paste(Delta," AIC")),colors = brewer.pal(n = 5, name = "Spectral"))+
  #scale_x_discrete(position = "top",name = "Experiment") +
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))


# alternative graph for model preference

BestFull <- fitted %>%
  group_by(genmodel,id) %>%
  mutate(delAIC = AIC - min(AIC)) %>% ungroup() %>%
  #mutate(model = mapvalues(model,from = oldmodellabels, to = newmodellabels)) %>%
  mutate(model = factor(model,levels=modellevels))


TotsSs <- NULL
for (Ss in models){

  AICmeans <- fitted %>%
    filter(genmodel %in% Ss) %>%
    group_by(model) %>%
    summarize(means = mean(AIC)) %>%
    arrange(means)


  CVSstot<- BestFull %>%
    filter(genmodel %in% Ss) %>%
    mutate(deltaAIC = AIC - AICmeans %>% slice(1) %>% .$means) %>%
    group_by(model) %>%
    summarize(mdeltaAIC = mean(deltaAIC)) %>%
    mutate(mdeltaAICadj = ifelse(mdeltaAIC > 20,20,mdeltaAIC)) %>%
    arrange(model) %>%
    mutate(genmodel = Ss,
           totalppl = 40)



  WinningSs1 <-BestFull %>%
    filter(genmodel %in% Ss) %>%
    #filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
    filter(delAIC == 0) %>% # only 1 model can win
    group_by(model) %>%
    summarize(numBest = length(delAIC)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/40) %>%
    arrange(model)

  WinningSs2 <- BestFull %>%
    filter(genmodel %in% Ss) %>%
    filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
    #filter(delAIC == 0) %>% # only 1 model can win
    group_by(model) %>%
    summarize(numBest = length(delAIC)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/40) %>%
    arrange(model)

  CVSstot$BestModelstrict <- WinningSs1 %>% filter(model %in% unique(CVSstot$model)) %>% .$wu
  CVSstot$BestModellenient <- WinningSs2 %>% filter(model %in% unique(CVSstot$model)) %>% .$wu

  TotsSs <- bind_rows(TotsSs,CVSstot)
}

TotsSs <- TotsSs %>%
  mutate(genmodel = factor(genmodel,levels = modellevels)) %>%
  arrange(genmodel,model) %>%
  mutate(whichmodel = rep(c(1:length(models)),each = length(models)))


ylim.prim <- c(0,20)
ylim.sec <- c(0,1)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

TotsSs$genmodel <- factor(TotsSs$genmodel,
                      levels = modellevels,
                      labels = modellabels)

# ggplot(TotsSs,aes(y = mdeltaAICadj,x = model,fill = mdeltaAICadj )) +
#
#   geom_rect(aes(xmin=0, xmax=11, ymin=0, ymax=1),fill="grey") +
#   geom_rect(aes(xmin=whichmodel-0.5,
#                 xmax=whichmodel+0.5,ymin = 0,ymax=20),
#             fill = "transparent",colour="black",linetype="dashed")+
#   geom_bar(stat = "identity") +
#   geom_point(aes(y = a + BestModelstrict*b),
#              color = "black",size=4,shape=21,fill="transparent") +
#   geom_point(aes(y = a + BestModellenient*b),
#              color = "black",size=4) +
#   coord_flip(ylim = c(0,20)) +
#   scale_fill_gradientn(name = expression(paste(Delta," AIC")),
#                        colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
#                        limits=c(0,20))+
#   scale_y_continuous(name = expression(paste(Delta," AIC")),
#                      breaks = seq(0,20,4),
#                      labels = c(seq(0,16,4),"..."),
#                      sec.axis = sec_axis(~ (. - a)/b,
#                                          name = "Proportion datasets best fit",
#                                          labels = c("0",".25",".50",".75","1"))) +
#   scale_x_discrete(name = "Recovering Model",
#                    labels = graphmodellabels) +
#   facet_wrap(.~genmodel,ncol = 5,labeller= label_parsed)+
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.ticks.y = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         legend.background = element_rect(fill="transparent"),
#         strip.background = element_rect(fill = "white"),
#         strip.placement = "outside",
#         panel.border = element_rect(colour = "black",fill = "transparent"))


# Proportion of data sets recovered by models

modelprop_strict_median <- ggplot(TotsSs,aes(x = genmodel,y = model,fill = BestModelstrict)) +
  geom_tile() +
  geom_text(aes(label = sub("^[0]+", "", round(BestModelstrict,2))),size = 3.5) +
  scale_fill_gradientn(name = "Proportion\ndatasets\nbest fit",
                       colors = rev(brewer.pal(n = 5, name = "Spectral")),
                       limits = c(0,1))+
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))

modelprop_lenient_median <- ggplot(TotsSs,aes(x = genmodel,y = model,fill = BestModellenient)) +
  geom_tile() +
  geom_text(aes(label = sub("^[0]+", "", round(BestModellenient,2))),size = 3.5) +
  scale_fill_gradientn(name = "Proportion\ndatasets\nbest fit",
                       colors = rev(brewer.pal(n = 5, name = "Spectral")),
                       limits = c(0,1))+
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
      axis.text.y = element_text(size=12),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white"))



# multivariate samples ---------------------------------------------------------

modelrecovery <- load_files("ModelRecoveryrmv/","rec")
parameterrecovery <- load_files("ParameterRecoveryrmv","trialsSs120")


modrec <- modelrecovery %>%
  group_by(genmodel,model,id) %>%
  arrange(objective) %>%
  slice(1) %>%
  group_by(model,genmodel) %>%
  mutate(nparam = length(vwmvp::get_start_vp(model)))
parrec <- parameterrecovery %>%
  group_by(model,id) %>%
  arrange(objective) %>%
  slice(1) %>%
  group_by(model) %>%
  mutate(nparam = length(vwmvp::get_start_vp(model))) %>%
  mutate(genmodel = model)

# combine model and parameter recovery

fitted <- bind_rows(modrec,parrec) %>%
  arrange(genmodel) %>%
  mutate(AIC = 2 * objective + 2 * nparam) %>%
  group_by(genmodel,model)

mr <- fitted %>%
  summarize(avAIC = mean(AIC)) %>%
  mutate(delAIC = avAIC - min(avAIC)) %>%
  mutate(deladjAIC = ifelse(delAIC > 6,6,delAIC))



mr$model <- factor(mr$model,
                   levels = modellevels,
                   labels = modellabels)
mr$genmodel <- factor(mr$genmodel,
                      levels = modellevels,
                      labels = modellabels)

modeldevRMV <- ggplot(mr,aes(x = genmodel,y = model,fill = deladjAIC)) +
  geom_tile() +
  # geom_text(data = Nleft,aes(y = 11, x = exp,label = paste0(nppt,"/",ntot)),vjust=1.5,size=4) +
  geom_text(aes(label = round(delAIC)),size = 3.5) +
  scale_fill_gradientn(name = expression(paste(Delta," AIC")),colors = brewer.pal(n = 5, name = "Spectral"))+
  #scale_x_discrete(position = "top",name = "Experiment") +
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))



# alternative graph for model preference

BestFull <- fitted %>%
  group_by(genmodel,id) %>%
  mutate(delAIC = AIC - min(AIC)) %>% ungroup() %>%
  #mutate(model = mapvalues(model,from = oldmodellabels, to = newmodellabels)) %>%
  mutate(model = factor(model,levels=modellevels))


TotsSs <- NULL
for (Ss in models){

  AICmeans <- fitted %>%
    filter(genmodel %in% Ss) %>%
    group_by(model) %>%
    summarize(means = mean(AIC)) %>%
    arrange(means)


  CVSstot<- BestFull %>%
    filter(genmodel %in% Ss) %>%
    mutate(deltaAIC = AIC - AICmeans %>% slice(1) %>% .$means) %>%
    group_by(model) %>%
    summarize(mdeltaAIC = mean(deltaAIC)) %>%
    mutate(mdeltaAICadj = ifelse(mdeltaAIC > 20,20,mdeltaAIC)) %>%
    arrange(model) %>%
    mutate(genmodel = Ss,
           totalppl = 40)



  WinningSs1 <-BestFull %>%
    filter(genmodel %in% Ss) %>%
    #filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
    filter(delAIC == 0) %>% # only 1 model can win
    group_by(model) %>%
    summarize(numBest = length(delAIC)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/40) %>%
    arrange(model)

  WinningSs2 <- BestFull %>%
    filter(genmodel %in% Ss) %>%
    filter(delAIC < 1) %>% # assumes equivalence of models within 1 Deviance point
    #filter(delAIC == 0) %>% # only 1 model can win
    group_by(model) %>%
    summarize(numBest = length(delAIC)) %>%
    complete(model, fill = list(0)) %>%
    group_by(model) %>%
    summarize(wu = sum(numBest,na.rm=T)) %>%
    mutate(wu = wu/40) %>%
    arrange(model)

  CVSstot$BestModelstrict <- WinningSs1 %>% filter(model %in% unique(CVSstot$model)) %>% .$wu
  CVSstot$BestModellenient <- WinningSs2 %>% filter(model %in% unique(CVSstot$model)) %>% .$wu

  TotsSs <- bind_rows(TotsSs,CVSstot)
}

TotsSs <- TotsSs %>%
  mutate(genmodel = factor(genmodel,levels = modellevels)) %>%
  arrange(genmodel,model) %>%
  mutate(whichmodel = rep(c(1:length(models)),each = length(models)))


ylim.prim <- c(0,20)
ylim.sec <- c(0,1)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

TotsSs$genmodel <- factor(TotsSs$genmodel,
                          levels = modellevels,
                          labels = modellabels)
#
# ggplot(TotsSs,aes(y = mdeltaAICadj,x = model,fill = mdeltaAICadj )) +
#
#   geom_rect(aes(xmin=0, xmax=11, ymin=0, ymax=1),fill="grey") +
#   geom_rect(aes(xmin=whichmodel-0.5,
#                 xmax=whichmodel+0.5,ymin = 0,ymax=20),
#             fill = "transparent",colour="black",linetype="dashed")+
#   geom_bar(stat = "identity") +
#   geom_point(aes(y = a + BestModelstrict*b),
#              color = "black",size=4,shape=21,fill="transparent") +
#   geom_point(aes(y = a + BestModellenient*b),
#              color = "black",size=4) +
#   coord_flip(ylim = c(0,20)) +
#   scale_fill_gradientn(name = expression(paste(Delta," AIC")),
#                        colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
#                        limits=c(0,20))+
#   scale_y_continuous(name = expression(paste(Delta," AIC")),
#                      breaks = seq(0,20,4),
#                      labels = c(seq(0,16,4),"..."),
#                      sec.axis = sec_axis(~ (. - a)/b,
#                                          name = "Proportion datasets best fit",
#                                          labels = c("0",".25",".50",".75","1"))) +
#   scale_x_discrete(name = "Recovering Model",
#                    labels = graphmodellabels) +
#   facet_wrap(.~genmodel,ncol = 5,labeller= label_parsed)+
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12),
#         axis.ticks.y = element_blank(),
#         panel.background = element_rect(fill = "white"),
#         legend.background = element_rect(fill="transparent"),
#         strip.background = element_rect(fill = "white"),
#         strip.placement = "outside",
#         panel.border = element_rect(colour = "black",fill = "transparent"))
#

# Proportion of data sets recovered by models

modelprop_strict_rmv <- ggplot(TotsSs,aes(x = genmodel,y = model,fill = BestModelstrict)) +
  geom_tile() +
  geom_text(aes(label = sub("^[0]+", "", round(BestModelstrict,2))),size = 3.5) +
  scale_fill_gradientn(name = "Proportion\ndatasets\nbest fit",
                       colors = rev(brewer.pal(n = 5, name = "Spectral")),
                       limits = c(0,1))+
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))

modelprop_lenient_rmv <- ggplot(TotsSs,aes(x = genmodel,y = model,fill = BestModellenient)) +
  geom_tile() +
  geom_text(aes(label = sub("^[0]+", "", round(BestModellenient,2))),size = 3.5) +
  scale_fill_gradientn(name = "Proportion\ndatasets\nbest fit",
                       colors = rev(brewer.pal(n = 5, name = "Spectral")),
                       limits = c(0,1))+
  scale_y_discrete(name = "Recovering Model", labels = graphmodellabels) +
  scale_x_discrete(name = "Generating Model", labels = graphmodellabels) +
  theme(axis.text.x = element_text(size=12,angle=90,hjust=0.5),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"))
