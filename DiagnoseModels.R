#library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library(vwmvp)
library("grid")
library("gridExtra")

library(reshape2)
options(pillar.sigfig = 5)


# Speed of Fitting -------------------------------------------------------------

FittedFull <- readRDS("Fits/FitFull.rds")
models <- c("MK_RNminus","MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")
models <- c("J_RNminus","J_RNplus","MK_RNminus","MK_RNplus")
Speeds <- FittedFull %>% filter(model %in% models) %>%
  group_by(model) %>%
  summarize(meanT = mean(time),
            medianT = median(time),
            iqrT = IQR(time),
            sdT = sd(time),
            meanIt = mean(iterations),
            medianIt = median(iterations),
            iqrIt = IQR(iterations),
            sdIt = sd(iterations))


# Parameters prep   -----------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,
                            data_pratte17)


models <- c("J_RNminus","J_RNplus","MK_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_FM_RNminus",
            "MK_RNplus","MK_P_RNplus","MK_U_RNplus","MK_FM_RNplus")
# models <- c("J_RNminus","J_RNplus","MK_RNminus",
#             "MK_RNplus")
numbermodels <- length(models)

FittedFull <- readRDS("Fits/FitFull.rds")
numberid <- length(unique(FittedFull$cvid))

exclall <- c()

ntrials <- experimentfile %>% group_by(id) %>% count() %>%
  slice(rep(1:n(),each=length(models))) %>% filter(!id %in% c(exclall)) %>%
  mutate(id = as.character(id)) %>%
  arrange(id)

BestFull <- FittedFull %>%
  filter(model %in% models) %>%
  filter(!cvid %in% c(exclall)) %>%
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


# Parameter estimates table ------------------------
maxsetsizes <- experimentfile  %>%
  mutate_if(is.factor,as.character)%>%
  group_by(id) %>%
  summarize(maxSs = max(set_size)) %>%
  arrange(id) %>%
  # best two runs for number of models
  slice(rep(1:n(),each=numbermodels))


parameters <- c("alpha","tau","kappa_r","K","J1bar","mkappa1")

RecodeParameters <- BestFull %>% ungroup() %>%
  select(model,cvid,all_of(parameters)) %>%
  bind_cols(maxSs = maxsetsizes$maxSs) %>%
  mutate(K = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                           K > maxSs),maxSs,K)) %>%
  mutate(K = ifelse((model == "MK_U_RNminus" | model == "MK_U_RNplus"),K/2,K)) %>%
  gather(key = "Parameter",value = "value",parameters) %>%
  group_by(model,Parameter) %>%
  mutate(Parameter = recode(Parameter, "mkappa1" = "phi",
                            "J1bar" = "phi"))


#MK_FM_RNminus, MK_FM_RNplus models: proportion where K >= max set size:

RecodeParameters %>%
  filter(model %in% c("MK_FM_RNminus","MK_FM_RNplus") &
           Parameter == "K") %>%
  group_by(model) %>%
  summarize(abovemaxSs = length(value[value==maxSs])/195)


ParEstimates <- RecodeParameters %>%
  # Remove outliers
  # where any parameter value outside 2.5,97.5 for any model
  filter(value < quantile(value,.975,na.rm=T) & value > quantile(value,.025,na.rm=T))


ParEstimates %>% group_by(model,Parameter) %>%
  summarize(meanPar = mean(value,na.rm=T),
            lower = meanPar - (qnorm(0.975) * sd(value,na.rm=T)/sqrt(length(value))),
            upper = meanPar + (qnorm(0.975) * sd(value,na.rm=T)/sqrt(length(value))),
            medianPar = median(value,na.rm=T)) %>%
  mutate(Parameter = factor(Parameter,
                            levels=c("phi","alpha","tau","kappa_r","K"))) %>%
  arrange(model,Parameter) %>%
  filter(model == "MK_U_RNminus")


# Correlation and stability of parameter estimates  --------------
# Note MK_FM_RNminus, MK_FM_RNplus are corrected for K > max Ss

reparameter <- c("phi","alpha","tau","kappa_r","K")

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

corfun <- function(x,y) {

  result <- cor.test(forCor[[x]],forCor[[y]])
  tibble(V1 = x,
         V2 = y,
         cor = result$estimate[[1]],
         sig = result$p.value < .05)
}


Correlations <- NULL
for (modeln in models){


forCor <- ParEstimates %>%
  spread(Parameter,value) %>%
  filter(model == modeln) %>%
  # only correlate complete estimates,
  # i.e. where not 1 of the values is an outlier relative to whole set
  drop_na(any_of(c("phi",names(vwmvp::get_start_vp(modeln)))))


numberpar <- length(names(vwmvp::get_start_vp(modeln)))
selectee <- as.vector(get_lower_tri(matrix(c(1:(numberpar * numberpar)),
                               nrow=numberpar,ncol=numberpar)))
selectee <- selectee[!is.na(selectee)]

parameters <- reparameter[reparameter %in% c("phi",names(vwmvp::get_start_vp(modeln)))]
combinations <- tibble(V1 = parameters,
                       V2 = parameters) %>% expand.grid() %>%
  arrange(V1,V2)


CorrelationResults <- NULL
for (i in selectee){

  combined <- corfun(combinations[[i,1]],combinations[[i,2]])
  CorrelationResults <- bind_rows(CorrelationResults,combined)
}

CorrelationResults$model <- modeln
Correlations <- bind_rows(Correlations,CorrelationResults)

}


# Format correlation reporting

Correlations <- Correlations %>% mutate(labels = sub("^[0]+","",ifelse(cor < .9999, round(cor,2),NA)))  %>%
  mutate(corlabels = sub("^[-0]+","-",ifelse(cor < .9999, labels, NA))) %>%
  mutate(corlabels = ifelse(round(cor,2) == 0,".00",corlabels)) %>%
  mutate(model = factor(model,
                          levels = c(
                            "J_RNminus","MK_RNminus",
                            "MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
                            "J_RNplus","MK_RNplus","MK_FM_RNplus",
                            "MK_P_RNplus",
                            "MK_U_RNplus"),
                          labels = c(
                            '"VP("~J~")A-"','"VP("~kappa~")A-"',
                            '"VP("~kappa~")F-"', '"VP("~kappa~")P-"','"VP("~kappa~")U-"',
                            '"VP("~J~")A+"','"VP("~kappa~")A+"','"VP("~kappa~")F+"',
                           '"VP("~kappa~")P+"',
                            '"VP("~kappa~")U+"')))


Corgraph <- ggplot(data = Correlations, aes(V1, V2, fill = cor))+
  geom_tile(color = "transparent")+
  geom_text(aes(V1, V2,label =  corlabels),size=2.5)+
  geom_text(aes(V1, V2,label = ifelse(cor < .9999 & sig, "*",NA)),size=3,vjust = -0.5,hjust=-1.5)+
  facet_wrap(.~model,ncol=5,labeller= label_parsed) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  scale_x_discrete(labels = c(expression(paste(phi)),
                              expression(paste(alpha)),
                              expression(paste(tau)),
                              expression(paste(kappa[r])),
                              expression(paste(K)))) +
  scale_y_discrete(labels = c(expression(paste(phi)),
                              expression(paste(alpha)),
                              expression(paste(tau)),
                              expression(paste(kappa[r])),
                              expression(paste(K)))) +
  theme_minimal()+
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill = "#e8e8e8"),
        axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  coord_fixed()






# Stability of Parameter estimates (AIC)

maxsetsizes <- experimentfile  %>%
  mutate_if(is.factor,as.character)%>%
  group_by(id) %>%
  summarize(maxSs = max(set_size)) %>%
  arrange(id) %>%
  # best two runs for number of models
  slice(rep(1:n(),each=numbermodels * 2))

Gatherparameters <- FittedFull %>% filter(model %in% models) %>%
  group_by(model,exp,cvid) %>%
  arrange(model,exp,cvid,objective) %>%
  slice(1:2) %>% ungroup() %>%
  select(alpha,J1bar,tau,mkappa1,kappa_r,K,model,cvid,objective) %>%
  mutate(phi = ifelse(model %in% c("J_RNminus","J_RNplus"),J1bar,mkappa1)) %>%
  select(-J1bar,-mkappa1) %>%
  arrange(cvid,model) %>%
  bind_cols(maxss = maxsetsizes$maxSs) %>%
  mutate(K = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                       K > maxss),maxss,K)) %>%
  gather(key = "Parameters", value = "value", -model,-cvid,-maxss) %>%
  filter(!(value == Inf | is.na(value))) %>%
  group_by(model,cvid,Parameters) %>%
  mutate(diffpercent = (value[2]-value[1])/value[1],
         diff = value[2] - value[1]) %>%
  group_by(model,cvid,Parameters) %>%
  slice(1)  %>%
  mutate(Parameters =
           factor(Parameters,
                  levels=c("phi","alpha","tau","kappa_r","K","objective"),
                  labels=c("phi","alpha","tau","kappa[r]","K","objective"))) %>%
  mutate(model =
           factor(model,
                  levels = c("MK_RNminus","J_RNminus","MK_FM_RNminus",
                             "MK_P_RNminus","MK_U_RNminus",
                             "MK_RNplus","J_RNplus","MK_FM_RNplus",
                             "MK_P_RNplus","MK_U_RNplus"))) %>%
  mutate(modelnum = as.numeric(model))

AboveCrit <- Gatherparameters %>%
  mutate(above = ifelse(abs(diffpercent) > 2, 1, 0)) %>%
  group_by(model,modelnum,Parameters) %>%
  filter(!(above == Inf | is.na(above))) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(model,modelnum,Parameters) %>% filter(!notplotted %in% c(0,100) )


Parestg <- ggplot(Gatherparameters %>% filter(Parameters != "objective"),
                  aes(x =modelnum,y = diffpercent)) +

  geom_abline(slope=0,intercept = 0) +
  geom_rect(aes(xmin = 0.5,xmax = 10.5,ymin = -.2,ymax=.2),fill = "transparent",
            colour="grey",
            linetype="dotted")+
  geom_jitter(size = 1,shape = 16,alpha = 0.2, fill = "blue",colour="blue",
              width=0.01) +
  scale_x_continuous(name = "Model", breaks = c(1:10),labels = c(
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
  )) +
  geom_text(data = AboveCrit,
            mapping = aes(x = modelnum, y = -2,
                          label = paste(round(notplotted),"%")),size = 2.5) +
  scale_y_continuous(name = "% difference",limits = c(-2,2),
                     breaks = c(-1,-.2,0,.2,1),
                     labels = c("-100%","","0%","","100%")) +

  # geom_point(aes(x = factor(trials), y = 0),size = 1, colour = "black") +
  facet_grid(.~Parameters,labeller = label_parsed) +
  theme(axis.text.x = element_text(size=12,angle=90,vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = c(0.9,0.75))
Parestg

# Stability of parameter estimates versus Stability of fit

DevPar <- Gatherparameters %>%
  group_by(model,cvid) %>%
  #only plot model fits with reasonly stable (Delta Dev < 2 points) fits
  filter(any(Parameters == "objective" & diff < 2))

DevParlength <- as.numeric(DevPar %>%
                             filter(Parameters != "objective") %>%
                             group_by(model,cvid) %>%
                             count() %>% .$n)
DevPargraph <- DevPar %>%
  filter(Parameters != "objective") %>%
  bind_cols(objectivediff =
              rep(DevPar %>% filter(Parameters == "objective") %>% .$diff,
                  DevParlength))


DevParstabg <- ggplot(DevPargraph,
                      aes(x = objectivediff,
                          y = diffpercent,fill = model,shape=model)) +
  geom_point(alpha = 0.7) +
  facet_grid(.~Parameters,labeller = label_parsed) +
  scale_fill_discrete(name = "Model", labels = c(
    expression(paste("VP(",italic(J),")A-")),
    expression(paste("VP(",kappa,")A-")),

    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("VP(",italic(J),")A+")),
    expression(paste("VP(",kappa,")A+")),

    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+"))
  )) +
  scale_shape_manual(name="Model",values = c(rep(21,5),rep(22,5)),
                     labels = c(
                       expression(paste("VP(",italic(J),")A-")),
                       expression(paste("VP(",kappa,")A-")),

                       expression(paste("VP(",kappa,")F-")),
                       expression(paste("VP(",kappa,")P-")),
                       expression(paste("VP(",kappa,")U-")),
                       expression(paste("VP(",italic(J),")A+")),
                       expression(paste("VP(",kappa,")A+")),

                       expression(paste("VP(",kappa,")F+")),
                       expression(paste("VP(",kappa,")P+")),
                       expression(paste("VP(",kappa,")U+"))))+

  scale_y_continuous(name = "abs (% difference)",limits = c(0,2),
                     breaks = c(0,1,2),labels = c("0%","100%","200%")) +
  scale_x_continuous(name = expression(paste(Delta,"Deviance (in points)")),
                     breaks = c(0,1,2)) +

  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size = 12),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(fill = "transparent",color="transparent"),
        legend.title = element_blank())
DevParstabg


plot_grid(Parestg,DevParstabg,Corgraph,nrow=3,rel_heights = c(1,1,1.3),labels = c("A","B","C"))
#ggsave("ParameterStability.png", units="cm", width=30, height=35, dpi=600)

# Stability of Fit (AIC) numerical integration figures -------------------------

models <- c("MK_RNminus","J_RNminus","MK_RNplus","J_RNplus")

numbermodels <- length(models)


AllDiff <- FittedFull %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  mutate(diffBest = Dev - Dev[[1]]) %>%
  mutate(model = factor(model,
                        levels = c("MK_RNminus","J_RNminus","MK_FM_RNminus",
                                   "MK_P_RNminus","MK_U_RNminus",
                                   "MK_RNplus","J_RNplus","MK_FM_RNplus",
                                   "MK_P_RNplus","MK_U_RNplus"))) %>%
  mutate(modelnum = as.numeric(model)) %>%
  filter(model %in% models) %>%
  mutate(modelposition = model) %>%
  mutate(modelposition = factor(modelposition,
                                levels = models)) %>%
  mutate(position = as.numeric(modelposition)) %>%
  select(model,exp,cvid,rep,objective,diffBest)

ggplot(AllDiff, aes(x =cvid,y=diffBest,color=exp)) +
  facet_wrap(.~model) +
  scale_y_continuous(limits = c(0,50))+
  geom_jitter()

# Best 2 runs

models <- c("MK_RNminus", "MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")
models <- c("MK_RNminus","J_RNminus","MK_RNplus","J_RNplus")
numbermodels <- length(models)

ntrials <- experimentfile %>% group_by(id) %>% count() %>%
  slice(rep(1:n(),each=length(models))) %>% filter(!id %in% c(exclall)) %>%
  mutate(id = as.character(id)) %>%
  arrange(id)

minDelta <- FittedFull %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  mutate(mindiffLL = outer(LL,LL, `-`)[2,1],
         bestLL = LL[[1]],
         diffDev = mindiffLL * 2,
         bestDev = bestLL * 2) %>%
  select(model,exp,cvid,diffDev,bestDev) %>%
  distinct() %>%
  arrange(cvid) %>%
  ungroup() %>%
  mutate(diffObj = diffDev/bestDev * 100,
         #trials = ntrials$n,
         id = cvid) %>%
  mutate(model = factor(model,
                        levels = c("MK_RNminus","J_RNminus","MK_FM_RNminus",
                                   "MK_P_RNminus","MK_U_RNminus",
                                   "MK_RNplus","J_RNplus","MK_FM_RNplus",
                                   "MK_P_RNplus","MK_U_RNplus"))) %>%
  mutate(modelnum = as.numeric(model)) %>%
  filter(model %in% models) %>%
  mutate(modelposition = model) %>%
  mutate(modelposition = factor(modelposition,
                                levels = models)) %>%
  mutate(position = as.numeric(modelposition))



# Deviance in Points figure:

## Text: data points per model not plotted in the graph

PtsAbove10 <- minDelta %>% mutate(above = ifelse(diffDev > 10, 1, 0)) %>%
  group_by(model,modelnum,position) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(model) %>% filter(!notplotted %in% c(0,100) )


## Text: proportion of data sets underneath a threshold per model

thresholdmodel <- minDelta %>%
  group_by(model,modelnum,position) %>%
  # in deviance terms...
  summarize(belowDev05 = length(diffDev[diffDev < 0.5]),
            belowDev1 = length(diffDev[diffDev < 1]),
            belowDev2 = length(diffDev[diffDev < 2]),
            belowDev10 = length(diffDev[diffDev < 10])) %>%
  gather(key = "crit",value="value",-model,-modelnum,-position) %>%
  group_by(model,modelnum,position,crit) %>%
  summarize(propbelow = value/numberid) %>%
  # y axis graph placement
  mutate(yplace = case_when(crit == "belowDev05" ~ 0.25,
                            crit == "belowDev1" ~ 0.75,
                            crit == "belowDev2" ~ 1.5,
                            crit == "belowDev10" ~ 6))


## Text: proportion of data sets underneath a treshold for any of the models

thresholdtotal <- minDelta %>% ungroup() %>%
  mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
         belowDev1 = ifelse(diffDev < 1,1,0),
         belowDev2 = ifelse(diffDev < 2,1,0),
         belowDev10 = ifelse(diffDev < 10,1,0)) %>%
  group_by(id) %>%
  summarize(modelsbelowDev05 = sum(belowDev05),
            modelsbelowDev1 = sum(belowDev1),
            modelsbelowDev2 = sum(belowDev2),
            modelsbelowDev10 = sum(belowDev10)) %>%
  gather(key = "crit",value="value",-id) %>%
  group_by(crit) %>%
  summarize(propbelow = length(value[value == numbermodels])/length(value)) %>%
  mutate(yplace = case_when(crit == "modelsbelowDev05" ~ 0.25,
                            crit == "modelsbelowDev1" ~ 0.75,
                            crit == "modelsbelowDev2" ~ 1.5,
                            crit == "modelsbelowDev10" ~ 6))

TotalPtsAbove10 <- (1 - thresholdtotal %>%
                      filter(crit == "modelsbelowDev10") %>%
                      .$propbelow) * 100

ggplotColours <- function(n = 10, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

colorin <- ggplotColours(n = 10)

plotcolors <- colorin[unique(thresholdmodel$modelnum)]

datastabpoints <- ggplot(minDelta %>% filter(diffDev <= 10),
                         aes(x = modelnum,y = diffDev )) +

  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0,ymax=0.5),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 0.5,ymax=1),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 1,ymax=2),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 2,ymax=10),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = position,xmax = position + 0.4,ymin = 10,ymax=10.5),
            fill = "transparent",colour="grey")+

  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 0,ymax=0.5),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 0.5,ymax=1),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 1,ymax=2),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 2,ymax=10),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(position) + 1,xmax = max(position)+ 1.5,
                ymin = 10,ymax=10.5),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = 0.75,xmax = max(position) + 1.5,
                ymin = 10.4,ymax=10.5),fill = "white",colour="white")+


  geom_abline(aes(intercept=0.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=2,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=10,slope = 0),linetype="dotted",colour="grey") +

  geom_text(data = thresholdmodel,
            mapping = aes(x = position + 0.2, y = yplace,
                          color = factor(modelnum),
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            size = 3.5) +
  geom_text(data = thresholdtotal,
            mapping = aes(x = numbermodels + 1.25, y = yplace,
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            color = "#A9A9A9",size = 3.5) +

  geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x= position - 0.2,
                                                           color = factor(modelnum)),
              size = 2,alpha=0.6,
              width=0.15) +
  geom_jitter(data = minDelta %>% filter(diffDev <10), aes(x = max(position) + 0.8,
                                                           y = diffDev),
              color = "#A9A9A9",size = 2,alpha=0.6,width=0.2) +
  scale_color_manual(values = plotcolors)+
  geom_text(data = PtsAbove10,
            mapping = aes(x = position - 0.3, y = 10.5,
                          label = paste(round(notplotted),"%")),size = 2.5) +
  annotate("text",x = numbermodels + .75, y = 10.5,
           label = paste(round(TotalPtsAbove10),"%"),size = 2.5) +

  scale_y_continuous(name = expression(paste(Delta,"Deviance (in points)")),
                     limits = c(0,10.5),
                     breaks = c(c(0,0.5,1),seq(2,10,2))) +
  # scale_x_continuous(name = "Model", limits = c(.75,numbermodels + 1.5),
  #                    breaks = c(c(1:numbermodels),numbermodels + 1),
  #                    labels = c(expression(paste("VP(",kappa,")A-")),
  #                               expression(paste("VP(",italic(J),")A-")),
  #                               expression(paste("VP(",kappa,")A+")),
  #                               expression(paste("VP(",italic(J),")A+")),
  #                               "total")
  #                    ) +
  scale_x_continuous(name = "Model", limits = c(.75,numbermodels + 1.5),
                     breaks = c(c(1:numbermodels),numbermodels + 1),
                     labels = c(

                       expression(paste("VP(",kappa,")A-")),
                       expression(paste("VP(",kappa,")F-")),
                       expression(paste("VP(",kappa,")P-")),
                       expression(paste("VP(",kappa,")U-")),
                       expression(paste("VP(",kappa,")A+")),
                       expression(paste("VP(",kappa,")F+")),
                       expression(paste("VP(",kappa,")P+")),
                       expression(paste("VP(",kappa,")U+")),
                       "total"
                     )
  ) +
  theme(axis.text.x = element_text(size=12),#,angle=90,vjust=0.5
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent",colour="transparent"))
datastabpoints

# Deviance in Percentage terms figure:

## Text: data points per model not plotted in the graph

PtsAbove <- minDelta %>% mutate(above = ifelse(diffObj > 0.5, 1, 0)) %>%
  group_by(model,modelnum) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(model) %>% filter(!notplotted %in% c(0,100) )


## Text: proportion of data sets underneath a threshold per model

thresholdmodelPer <- minDelta %>%
  group_by(model,modelnum) %>%
  # in deviance terms...
  summarize(belowPer025 = length(diffObj[diffObj < 0.025]),
            belowPer05 = length(diffObj[diffObj < .05]),
            belowPer1 = length(diffObj[diffObj < .1]),
            belowPer5 = length(diffObj[diffObj < .5])) %>%
  gather(key = "crit",value="value",-model,-modelnum) %>%
  group_by(model,modelnum,crit) %>%
  summarize(propbelow = value/numberid) %>%
  # y axis graph placement
  mutate(yplace = case_when(crit == "belowPer025" ~ 0.0125,
                            crit == "belowPer05" ~ 0.0375,
                            crit == "belowPer1" ~ 0.075,
                            crit == "belowPer5" ~ 0.3))


## Text: proportion of data sets underneath a treshold for any of the models

thresholdtotalPer <- minDelta %>% ungroup() %>%
  mutate(belowPer025 = ifelse(diffObj < 0.025,1,0),
         belowPer05 = ifelse(diffObj < .05,1,0),
         belowPer1 = ifelse(diffObj < .1,1,0),
         belowPer5 = ifelse(diffObj < .5,1,0)) %>%
  group_by(id) %>%
  summarize(modelsbelowPer025 = sum(belowPer025),
            modelsbelowPer05 = sum(belowPer05),
            modelsbelowPer1 = sum(belowPer1),
            modelsbelowPer5 = sum(belowPer5)) %>%
  gather(key = "crit",value="value",-id) %>%
  group_by(crit) %>%
  summarize(propbelow = length(value[value == numbermodels])/length(value)) %>%
  mutate(yplace = case_when(crit == "modelsbelowPer025" ~0.0125,
                            crit == "modelsbelowPer05" ~ 0.0375,
                            crit == "modelsbelowPer1" ~ 0.075,
                            crit == "modelsbelowPer5" ~ 0.3))

TotalPtsAbove <- (1 - thresholdtotalPer %>%
                    filter(crit == "modelsbelowPer5") %>%
                    .$propbelow) * 100

datastabpercent <- ggplot(minDelta %>% filter(diffObj <= 0.5),
                          aes(x = modelnum,y = diffObj )) +

  geom_rect(aes(xmin = modelnum,xmax = modelnum + 0.4,ymin = 0,ymax=0.025),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = modelnum,xmax = modelnum + 0.4,ymin = 0.025,ymax=.05),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = modelnum,xmax = modelnum + 0.4,ymin = .05,ymax=.1),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = modelnum,xmax = modelnum + 0.4,ymin = .1,ymax=.5),
            fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = modelnum,xmax = modelnum + 0.4,ymin = .5,ymax=.55),
            fill = "transparent",colour="grey")+

  geom_rect(aes(xmin = max(modelnum) + 1.5,xmax = max(modelnum)+ 1.5 + 1,
                ymin = 0,ymax=0.025),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(modelnum) + 1.5,xmax = max(modelnum)+ 1.5 + 1,
                ymin = 0.025,ymax=.05),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(modelnum) + 1.5,xmax = max(modelnum)+ 1.5 + 1,
                ymin = .05,ymax=.1),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(modelnum) + 1.5,xmax = max(modelnum)+ 1.5 + 1,
                ymin = .1,ymax=.5),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = max(modelnum) + 1.5,xmax = max(modelnum)+ 1.5 + 1,
                ymin = .5,ymax=.55),fill = "transparent",colour="grey")+
  geom_rect(aes(xmin = 0.75,xmax = max(modelnum) + 1.5 + 1,ymin = .53,ymax=.55),
            fill = "white",colour="white")+
  #
  #
  geom_abline(aes(intercept=0.025,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0.05,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +

  geom_text(data = thresholdmodelPer,
            mapping = aes(x = modelnum + 0.2, y = yplace,
                          color = model,
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            size = 3.5) +
  geom_text(data = thresholdtotalPer,
            mapping = aes(x = numbermodels + 2, y = yplace,
                          label = sub("^0+", "", paste(round(propbelow,2)))),
            color = "#A9A9A9",size = 3.5) +

  geom_jitter(aes(x= modelnum - 0.2, color = factor(model)),size = 2,alpha=0.6,
              width=0.15) +
  geom_jitter(data = minDelta %>% filter(diffObj < 0.5),
              aes(x = max(modelnum) + 1,y = diffObj),
              color = "#A9A9A9",size = 2,alpha=0.6) +
  geom_text(data = PtsAbove,
            mapping = aes(x = modelnum - 0.3, y = 0.52,
                          label = paste(round(notplotted),"%")),size = 2.5) +
  annotate("text",x = numbermodels + 1, y = 0.52,
           label = paste(round(TotalPtsAbove),"%"),size = 2.5) +

  scale_y_continuous(name = expression(paste(Delta,"Deviance (in %)")),
                     limits = c(0,0.55),
                     breaks = c(c(0,0.025,.05),seq(.1,.5,.1)),
                     labels = c("0%","0.025%","0.05%","0.10%","0.20%",
                                "0.30%","0.40%","0.50%")) +
  scale_x_continuous(name = "Model", limits = c(.75,numbermodels + 2.5),
                     breaks = c(c(1:numbermodels),numbermodels + 1.5),
                     labels = c(expression(paste("VP(",kappa,")A-")),
                                expression(paste("VP(",italic(J),")A-")),
                                expression(paste("VP(",kappa,")F-")),
                                expression(paste("VP(",kappa,")P-")),
                                expression(paste("VP(",kappa,")U-")),
                                expression(paste("VP(",kappa,")A+")),
                                expression(paste("VP(",italic(J),")A+")),
                                expression(paste("VP(",kappa,")F+")),
                                expression(paste("VP(",kappa,")P+")),
                                expression(paste("VP(",kappa,")U+")),
                                "total")) +
  theme(axis.text.x = element_text(size=12,angle=90,vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent",colour="transparent"))
datastabpercent


# Stability of Fit (AIC) GA, GA + numeric integration table --------------------

GAonly <- readRDS("Fits/GA_FitFull.rds")
GAnlminb <- readRDS("Fits/GAnlminb_FitFull.rds")


GADelta <- GAonly %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  arrange(model,exp,cvid,-objective) %>%
  slice(1:2) %>%
  mutate(bestDev = -2 * objective,
         diffDev = diff(-objective),
         id = cvid) %>%
  arrange(model,exp,cvid) %>%
  slice(1) %>%
  mutate(model = factor(model,
                        levels = c("MK_RNminus","J_RNminus",
                                   "J_RNplus","MK_RNplus"))) %>%
  mutate(modelnum = as.numeric(model))

## Table Simulation: proportion of data sets that exceed a threshold

thresholdmodelGA <- GADelta %>%
  group_by(model,modelnum) %>%
  # in deviance terms...
  summarize(belowDev05 = length(diffDev[diffDev < 0.5]),
            belowDev1 = length(diffDev[diffDev < 1]),
            belowDev2 = length(diffDev[diffDev < 2]),
            belowDev10 = length(diffDev[diffDev < 10])) %>%
  gather(key = "crit",value="value",-model,-modelnum) %>%
  group_by(model,modelnum,crit) %>%
  summarize(propbelow = value/numberid)

thresholdtotalGA <- GADelta %>% ungroup() %>%
  mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
         belowDev1 = ifelse(diffDev < 1,1,0),
         belowDev2 = ifelse(diffDev < 2,1,0),
         belowDev10 = ifelse(diffDev < 10,1,0)) %>%
  group_by(id) %>%
  summarize(modelsbelowDev05 = sum(belowDev05),
            modelsbelowDev1 = sum(belowDev1),
            modelsbelowDev2 = sum(belowDev2),
            modelsbelowDev10 = sum(belowDev10)) %>%
  gather(key = "crit",value="value",-id) %>%
  group_by(crit) %>%
  summarize(propbelow = length(value[value == 4])/length(value))

GAnlminbDelta <- GAnlminb %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  arrange(model,exp,cvid,objective) %>%
  slice(1:2) %>%
  mutate(bestDev = -2 * -objective,
         diffDev = diff(objective),
         id = cvid) %>%
  arrange(model,exp,cvid) %>%
  slice(1) %>%
  mutate(model = factor(model,
                        levels = c("MK_RNminus","J_RNminus",
                                   "J_RNplus","MK_RNplus"))) %>%
  mutate(modelnum = as.numeric(model))

## Table Simulation + numerical integration:
## proportion of data sets underneath a treshold for any of the models

thresholdmodelGAnlminb <- GAnlminbDelta %>%
  group_by(model,modelnum) %>%
  # in deviance terms...
  summarize(belowDev05 = length(diffDev[diffDev < 0.5]),
            belowDev1 = length(diffDev[diffDev < 1]),
            belowDev2 = length(diffDev[diffDev < 2]),
            belowDev10 = length(diffDev[diffDev < 10])) %>%
  gather(key = "crit",value="value",-model,-modelnum) %>%
  group_by(model,modelnum,crit) %>%
  summarize(propbelow = value/numberid)


thresholdtotalGAnlminb <- GAnlminbDelta %>% ungroup() %>%
  mutate(belowDev05 = ifelse(diffDev < 0.5,1,0),
         belowDev1 = ifelse(diffDev < 1,1,0),
         belowDev2 = ifelse(diffDev < 2,1,0),
         belowDev10 = ifelse(diffDev < 10,1,0)) %>%
  group_by(id) %>%
  summarize(modelsbelowDev05 = sum(belowDev05),
            modelsbelowDev1 = sum(belowDev1),
            modelsbelowDev2 = sum(belowDev2),
            modelsbelowDev10 = sum(belowDev10)) %>%
  gather(key = "crit",value="value",-id) %>%
  group_by(crit) %>%
  summarize(propbelow = length(value[value == 4])/length(value))

# Stability of Fit (AIC) numerical integration all models

FittedFull %>%
  arrange(model,exp,cvid,objective) %>%
  group_by(model,exp,cvid) %>%
  mutate(mindiffLL = outer(LL,LL, `-`)[2,1],
         bestLL = LL[[1]],
         diffDev = mindiffLL * 2,
         bestDev = bestLL * 2) %>%
  select(model,exp,cvid,diffDev,bestDev) %>%
  distinct() %>%
  arrange(cvid) %>%
  ungroup() %>%
  mutate(diffObj = diffDev/bestDev * 100,
         trials = ntrials$n,
         id = cvid) %>%
  mutate(model = factor(model,
                        levels = c("MK_RNminus","J_RNminus","MK_FM_RNminus",
                                   "MK_P_RNminus","MK_U_RNminus",
                                   "J_RNplus","MK_RNplus","MK_FM_RNplus",
                                   "MK_P_RNplus","MK_U_RNplus"))) %>%
  mutate(modelnum = as.numeric(model))

unique(FittedFull$model)

