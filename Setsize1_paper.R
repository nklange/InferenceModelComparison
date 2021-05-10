library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("grid")
library(magick)
library(cowplot)
library(ungeviz)


# COWPLOT: Make Set size 1 graphs p1 --------------------------------------------


models <- c("MK_RNminus","MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")
models2 <- c("MK_U_RNminus","MK_RNminus","MK_P_RNminus",
             "MK_FM_RNplus","MK_FM_RNminus","MK_RNplus","MK_P_RNplus","MK_U_RNplus")

# behavioral signature only for aggregate fit of 1 experiment
exps <- "ol17_e1"
experiment <- "ol17_e1"

# Load empirical data

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size") %>%
  filter(exp %in% exps)

expsetsize <- sort(unique(experimentfile %>% .$setsize))
setsizelabels <- paste("Set size:",expsetsize)
leftoutlabels <- c("Full data set fitted", paste("Left out: Ss",expsetsize))


experimentfile <- experimentfile %>%
  filter(exp == experiment)%>%
  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels),
         id = as.character(factor(id)),
         Predicted = NA)

# Make AIC

leftoutlabels <- c("Full data set fitted", paste("Left out: Ss",expsetsize))
leftoutlabels2 <- c(NA, paste("Ss",expsetsize))



# LOSSO out-of-sample prediction

LOTest <- readRDS("Fits/CVLOSzO_TestLL.rds") %>%
  filter(model %in% models2)

FailedSzPrediction <- unique(LOTest %>% filter(CVLL > 1e4) %>% .$id)

excludeppt <- length(FailedSzPrediction)

# LOSSO in-sample prediction

LOTraining <- readRDS("Fits/CVLOSzO_Trainingfits.rds")  %>%
  filter(model %in% models2) %>%
  filter(!id %in% FailedSzPrediction)

LOTest <- LOTest %>% filter(!id %in% FailedSzPrediction)

# AIC in-sample prediction

FullFit <- readRDS("Fits/FitFull.rds")%>%
  filter(model %in% models2) %>%
  filter(!id %in% FailedSzPrediction)



experimentfile <- experimentfile %>%
  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels)) %>%
  filter(!id %in% FailedSzPrediction)

participants <- unique(experimentfile$id)

DeltaAICFull <- FullFit %>%
  group_by(id,model) %>%
  arrange(AIC) %>% slice(1) %>%
  group_by(model) %>%
  summarize(mAIC = mean(AIC),
            lengths = length(AIC)) %>%
  mutate(deltaAIC = mAIC - min(mAIC)) %>%
  mutate(leftout = "Full data set fitted") %>%
  mutate(model = factor(model, levels = models2))

DeltaAICLOSsO <- LOTraining %>%
  group_by(model,id,leftout) %>%
  mutate(AIC = objective * 2 + 2 * nparams) %>%
  arrange(AIC) %>% slice(1) %>%
  group_by(leftout,model) %>%
  summarize(mAIC = mean(AIC),
            lengths = length(AIC)) %>%
  mutate(deltaAIC = mAIC - min(mAIC)) %>%
  mutate(leftout = factor(leftout,  levels = c(0,expsetsize),
                          labels =  leftoutlabels))  %>%
  mutate(model = factor(model, levels = models2))


AICDEV <- bind_rows(DeltaAICFull,DeltaAICLOSsO) %>%
  mutate(deltaAICadj = ifelse(deltaAIC > 20, 20,deltaAIC)) %>%
  mutate(type = "In-sample\nprediction (AIC)")

DeltaCVDev <- LOTest %>%
  filter(model %in% models) %>%
  group_by(model,id,leftout) %>%
  mutate(CVDev = CVLL * 2 ) %>%
  arrange(CVDev) %>% slice(1) %>%
  group_by(leftout,model) %>%
  summarize(mAIC = mean(CVDev),
            lengths = length(CVDev)) %>%
  mutate(deltaAIC = mAIC - min(mAIC)) %>%
  mutate(leftout = factor(leftout,  levels = c(0,expsetsize),
                          labels =  leftoutlabels2)) %>%
  mutate(deltaAICadj = ifelse(deltaAIC > 20, 20,deltaAIC)) %>%
  mutate(type = "Out-of-sample\nprediction (Dev)")  %>%
  mutate(model = factor(model, levels = models2))


Devs <- bind_rows(AICDEV,DeltaCVDev) %>%
  mutate(model = factor(model, levels = models2)) %>%
  mutate(type  = factor(type,
                        levels = c(
                                   "Out-of-sample\nprediction (Dev)",
                                   "In-sample\nprediction (AIC)"
                                   )))


Annotation <- Devs %>% select(type,leftout,lengths) %>% distinct() %>%
  #filter(leftout %in% leftoutlabels[selectleftouts + 1] ) %>%

  mutate(lengths = paste("N =",lengths)) %>%
  add_row(type = "Out-of-sample\nprediction (Dev)",
          leftout = "-",
          lengths = NA) %>%
  mutate(type  = factor(type,
                        levels = c(
                          "Out-of-sample\nprediction (Dev)",
                          "In-sample\nprediction (AIC)"
                        )))

Devs2 <- Devs %>% filter(type=="In-sample\nprediction (AIC)") %>% mutate(type = factor(type))
Annotation2 <- Annotation %>% filter(type=="In-sample\nprediction (AIC)") %>% mutate(type = factor(type))

AICDev2 <- ggplot(data=Devs2,aes(y = deltaAICadj,x = model)) +
  # geom_text(data = Annotation2, aes(y = 15, x = 4,label = paste("Delta")),parse=TRUE) +
  geom_text(data = Annotation2, aes(y = 10, x = 8.1,label =lengths)) +
  #
  geom_bar(data=Devs2,stat = "identity",aes(fill = deltaAICadj)) +

  coord_flip(ylim = c(0,20)) +
  scale_fill_gradientn(name = expression(paste(Delta)),
                       colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
                       limits=c(0,20))+
  geom_text(data = Devs2,
            aes(label=round(deltaAIC), y = 20), hjust = 1,size=3)+

  scale_y_continuous(name = expression(paste(Delta)), breaks = seq(0,20,5),
                     labels = c(seq(0,15,5),"...") ) +
  scale_x_discrete(name = "Model"
                   ,labels = c(
                     expression(paste("VP(",kappa,")U-")),
                     expression(paste("VP(",kappa,")A-")),
                     expression(paste("VP(",kappa,")P-")),
                     expression(paste("VP(",kappa,")F+")),
                     expression(paste("VP(",kappa,")F-")),
                     expression(paste("VP(",kappa,")A+")),
                     expression(paste("VP(",kappa,")P+")),
                     expression(paste("VP(",kappa,")U+"))
                   )
  ) +

  facet_grid(type~leftout,drop=FALSE) +
  theme(axis.text.x = element_text(size=11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),

        axis.text.y = element_text(size=11,hjust=0),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill ="transparent",colour ="transparent"),

        legend.position = "none",

        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=11))

AICDev2

Devs1 <- Devs %>% filter(type=="Out-of-sample\nprediction (Dev)") %>% mutate(type = factor(type))
Annotation1 <- Annotation %>% filter(type=="Out-of-sample\nprediction (Dev)") %>% mutate(type = factor(type))


AICDev1 <- ggplot(data=Devs1,aes(y = deltaAICadj,x = model)) +
  # geom_text(data = Annotation2, aes(y = 15, x = 4,label = paste("Delta")),parse=TRUE) +
  geom_text(data = Annotation1, aes(y = 10, x = 8.1,label =lengths)) +
  #
  geom_bar(data=Devs1,stat = "identity",aes(fill = deltaAICadj)) +

  coord_flip(ylim = c(0,20)) +
  scale_fill_gradientn(name = expression(paste(Delta)),
                       colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
                       limits=c(0,20))+
  geom_text(data = Devs1,
            aes(label=round(deltaAIC), y = 20), hjust = 1,size=3)+

  scale_y_continuous(name = expression(paste(Delta)), breaks = seq(0,20,5),
                     labels = c(seq(0,15,5),"...") ) +
  scale_x_discrete(name = "Model"
                   ,labels = c(
                     expression(paste("VP(",kappa,")U-")),
                     expression(paste("VP(",kappa,")A-")),
                     expression(paste("VP(",kappa,")P-")),
                     expression(paste("VP(",kappa,")F+")),
                     expression(paste("VP(",kappa,")F-")),
                     expression(paste("VP(",kappa,")A+")),
                     expression(paste("VP(",kappa,")P+")),
                     expression(paste("VP(",kappa,")U+"))
                   )
  ) +

  facet_grid(type~leftout,drop=FALSE) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=11,hjust=0),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill ="transparent",colour ="transparent"),

        legend.position = "none",

        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.placement = "outside",
        strip.text = element_text(size=11))

AICDev1

AICDev <- plot_grid(AICDev1,AICDev2,ncol=1,rel_heights = c(0.85,1))

# Errordistribution


# select elements for manuscript graph
exps <- "ol17_e1"
experiment <- exps
selectleftouts <- c(0,
                    1,
                    2,
                    8)
selectsetsizes <- c("Set size 1",
                    "Set size 2",
                    "Set size 8")

Prediction <- readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds")) %>%
  filter(model %in% models) %>%
  mutate(model = factor(model,levels = models)) %>%

  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels),
         leftout = factor(leftout,
                          levels = c(0,expsetsize),
                          labels = leftoutlabels))



Prediction<- Prediction  %>%
  mutate(Predicted = ifelse((leftout == "Left out: Ss 1" & setsize == "Set size: 1") |
                              (leftout == "Left out: Ss 2" & setsize == "Set size: 2") |
                              (leftout == "Left out: Ss 3" & setsize == "Set size: 3") |
                              (leftout == "Left out: Ss 4" & setsize == "Set size: 4")|
                              (leftout == "Left out: Ss 5" & setsize == "Set size: 5")|
                              (leftout == "Left out: Ss 6" & setsize == "Set size: 6")|
                              (leftout == "Left out: Ss 7" & setsize == "Set size: 7") |
                              (leftout == "Left out: Ss 8" & setsize == "Set size: 8"),
                            "Out-of-sample prediction","In-sample prediction"))

Annotation <- Prediction %>% select(setsize,leftout,Predicted) %>%
  distinct() %>% mutate(error_0 = Inf ) %>%
  filter(setsize %in% setsizelabels[selectleftouts[c(2:length(selectleftouts))]]) %>%
  filter(leftout %in% leftoutlabels[selectleftouts + 1])




load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size") %>%
  filter(exp %in% exps)

expsetsize <- sort(unique(experimentfile %>% .$setsize))
setsizelabels <- paste("Set size:",expsetsize)
leftoutlabels <- c("Full data set fitted", paste("Left out: Ss",expsetsize))


experimentfile <- experimentfile %>%
  filter(exp == experiment)%>%
  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels),
         id = as.character(factor(id)),
         Predicted = NA)

# Make error distribution

errordistribution <- ggplot(data = experimentfile %>%
                              filter(setsize %in% setsizelabels[selectleftouts[c(2:length(selectleftouts))]]),
                            aes(x = error_0)) +
  geom_rect(data = Annotation,xmin = -Inf,ymin = -Inf, xmax = Inf,ymax = Inf,aes(fill=Predicted),color="black") +
  scale_fill_manual(name=exps,values = c("white","#dedede"))+
  geom_histogram(aes(y = ..density..),color="darkgrey",fill="#edeaea",bins=60) +
  coord_cartesian(xlim = c(-pi/2,pi/2))+

  scale_y_continuous(limits=c(0,max(Prediction %>%
                                      filter(setsize %in% setsizelabels[selectleftouts[c(2:length(selectleftouts))]]) %>%
                                      filter(leftout %in% leftoutlabels[selectleftouts + 1] )  %>%

                                      .$prediction,na.rm=T) + 0.2)) +
  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(-pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_line(data = Prediction %>%
              filter(setsize %in% setsizelabels[selectleftouts[c(2:length(selectleftouts))]]) %>%
              filter(leftout %in% leftoutlabels[selectleftouts + 1]),

            aes(x = data,y=prediction,linetype=model,color=model,group=model), size = 0.8) +
  facet_grid(setsize~leftout) +

  scale_linetype_manual(name=experiment,
                        values = c("dashed","dashed","dashed","dashed",
                                   "solid","solid","solid","solid"),
                        labels = c(

                          expression(paste("VP(",kappa,")A-")),
                          expression(paste("VP(",kappa,")F-")),
                          expression(paste("VP(",kappa,")P-")),
                          expression(paste("VP(",kappa,")U-")),
                          expression(paste("VP(",kappa,")A+")),
                          expression(paste("VP(",kappa,")F+")),
                          expression(paste("VP(",kappa,")P+")),
                          expression(paste("VP(",kappa,")U+"))
                        )
  ) +

  scale_color_manual(name = experiment,
                     values =  c("#332288","#E69F00","#44AA99","#CC79A7",
                                 "#332288","#E69F00","#44AA99","#CC79A7"),
                     labels = c(
                       expression(paste("VP(",kappa,")A-")),
                       expression(paste("VP(",kappa,")F-")),
                       expression(paste("VP(",kappa,")P-")),
                       expression(paste("VP(",kappa,")U-")),
                       expression(paste("VP(",kappa,")A+")),
                       expression(paste("VP(",kappa,")F+")),
                       expression(paste("VP(",kappa,")P+")),
                       expression(paste("VP(",kappa,")U+"))
                     )
  )+
guides(linetype = guide_legend(override.aes = list(size = 1.2),ncol=4 ) ) +

theme(axis.text.x = element_text(size=11),
      plot.title = element_text(size = 14,hjust=0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size=12),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      strip.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent"),
      strip.text = element_text(size = 12),
      legend.position = c(0.3,0.2),
      legend.direction = "horizontal",
      legend.text = element_text(size = 13),
      legend.background = element_rect(fill="white",color="lightgrey"),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      legend.title = element_blank(),
      legend.key.width = unit(3,"line"))

errordistribution

# Summary statistics across set size

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size") %>%
  filter(exp %in% exps)

expsetsize <- sort(unique(experimentfile %>% .$setsize))
setsizelabels <- paste("Set size:",expsetsize)
leftoutlabels <- c("Full data set fitted", paste("Left out: Ss",expsetsize))


experimentfile <- experimentfile %>%
  filter(exp == experiment)%>%
  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels),
         id = as.character(factor(id)),
         Predicted = NA)


Prediction <- Prediction %>% filter(model %in% models)

CVtot <- NULL

for(leftouts in unique(Prediction %>% .$leftout)){

  CV <- NULL
  for(SS in unique(Prediction %>% .$setsize)){

    CircVars <- NULL

    for(modeln in unique(Prediction$model)){
      testNA <- any(is.na(Prediction %>%
                            filter(setsize == SS) %>%
                            filter(leftout == leftouts) %>%
                            filter(model == modeln) %>%
                            .$normprediction) == TRUE)
      if (testNA){

        CircVar <- NA
        CircSD <- NA
        MAE <- NA

      } else {

        samples <- sample(size=6e4,Prediction %>%
                            filter(setsize == SS) %>%
                            filter(leftout == leftouts) %>%
                            filter(model == modeln) %>%
                            .$data,replace=TRUE,
                          prob =  Prediction %>%
                            filter(setsize == SS) %>%
                            filter(leftout == leftouts) %>%
                            filter(model == modeln) %>%
                            .$normprediction)
        CircVar <- circular::var.circular(samples)
        CircSD <- circular::sd.circular(samples)
        MAE <- mean(abs(samples))
        kurtosis <- vwmvp::circkurtosis(samples)

      }

      out <- as_tibble(data.frame(setsize = SS,
                                  leftout = leftouts,
                                  model = modeln,
                                  Var = CircVar,
                                  SD = CircSD,
                                  meanErr = MAE,
                                  kurt = kurtosis))

      CircVars <- CircVars %>%  bind_rows(out)
    }


    expdata <- experimentfile %>%
      filter(exp == "ol17_e1") %>%
      filter(setsize == SS) %>%
      .$error_0



    cv <- circular::var.circular(expdata)
    csd <- circular::sd.circular(expdata)
    mae <- mean(abs(expdata))
    kurtosis_data <- vwmvp::circkurtosis(expdata)

    out2 <- as_tibble(data.frame(setsize = SS,
                                 leftout = leftouts,
                                 model = "Data",
                                 Var = cv,
                                 SD = csd,
                                 meanErr = mae,
                                 kurt = kurtosis_data))
    CircVars <- CircVars %>%  bind_rows(out2)
    CV <- CV %>% bind_rows(CircVars)
  }
  CVtot <- CVtot %>% bind_rows(CV)
}

CV <- CVtot



CV <- CV %>% mutate(model = factor(model,
                                   levels = c("Data",models)))


CVplot <- CV %>% select(-SD) %>%
  filter(model %in% c("Data","MK_RNminus","MK_RNplus")) %>%
  pivot_longer(!c(model,setsize,leftout),names_to = "measure",values_to="value") %>%
  mutate(measure = factor(measure, levels = c("Var","meanErr","kurt"),
                          labels = c(expression(paste("circular ",sigma^2)),
                                     expression(paste("mean abs. deviation (radian)")),
                                     expression(paste("kurtosis"))))) %>%
  filter(leftout %in% c("Full data set fitted","Left out: Ss 1","Left out: Ss 2")) %>%
  mutate(leftout = factor(leftout, levels = c("Full data set fitted","Left out: Ss 1","Left out: Ss 2"),
                          labels = c(expression(paste("Full data set fitted")),
                                     expression(paste("Left out: Ss 1")),
                                     expression(paste("Left out: Ss 2"))))) %>%
  mutate(model2 = paste0(model,"_",leftout)) %>%
  mutate(model2 = factor(model2, levels = c("MK_RNplus_paste(\"Full data set fitted\")",
                                            "MK_RNminus_paste(\"Full data set fitted\")",
                                            "MK_RNplus_paste(\"Left out: Ss 1\")",
                                            "MK_RNminus_paste(\"Left out: Ss 1\")",
                                            "MK_RNplus_paste(\"Left out: Ss 2\")",
                                            "MK_RNminus_paste(\"Left out: Ss 2\")",
                                            "Data_paste(\"No Ss left out\")",
                                            "Data_paste(\"Left out: Ss 1\")",
                                            "Data_paste(\"Left out: Ss 2\")")))


Vars1 <- ggplot(CVplot %>%
                  filter(model == "Data"),
                aes(x = setsize, y = value)) +
  geom_hpline(aes(linetype = model),size=1.5,width=1.2) +
  geom_point(data = CVplot %>% filter(model != "Data"),
             aes(x = setsize, y = value,shape = model2, fill = model2),
             position = position_dodge(0.6),size=4,alpha=0.7,color="black") +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_grid(~measure,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25))+
  scale_linetype_manual(values = "solid") +
  scale_shape_manual(values = rep(c(22,21),3),
                     labels = c(
                       expression(paste("VP(",kappa,")A+, full data set fitted")),
                       expression(paste("VP(",kappa,")A-, full data set fitted")),
                       expression(paste("VP(",kappa,")A+, left out: Ss 1")),
                       expression(paste("VP(",kappa,")A-, left out: Ss 1")),
                       expression(paste("VP(",kappa,")A+, left out: Ss 2")),
                       expression(paste("VP(",kappa,")A-, left out: Ss 2"))
                       )
  ) +
  scale_fill_manual(values = c("#D55E00","#ffc599","#0072B2","#99daff","#CC79A7", "#f1dae7"),
                    labels = c(
                      expression(paste("VP(",kappa,")A+, full data set fitted")),
                      expression(paste("VP(",kappa,")A-, full data set fitted")),
                      expression(paste("VP(",kappa,")A+, left out: Ss 1")),
                      expression(paste("VP(",kappa,")A-, left out: Ss 1")),
                      expression(paste("VP(",kappa,")A+, left out: Ss 2")),
                      expression(paste("VP(",kappa,")A-, left out: Ss 2"))
                    )
  ) +
  theme(axis.text = element_text(size=14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=14),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 14),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.1,0.6),
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
Vars1

# NRMSD

CVdat <- CV %>% filter(model == "Data") %>%  slice(rep(1:n(),each=8)) %>% arrange(leftout,setsize,model)


# ALL
CVdat_meandat <- CVdat %>% group_by(model) %>% summarize(meanE = mean(meanErr),
                                                         meanV = mean(Var),
                                                         meank = mean(kurt)) %>%
  pivot_longer(!c(model),names_to = "measure",values_to="value") %>%
  slice(rep(1:n(),8*9))



CVdat2 <- CV %>% filter(model != "Data") %>%
  arrange(leftout,setsize,model) %>%
  mutate(MeanErrodiff = (meanErr - CVdat$meanErr)^2,
         Vardiff = (Var - CVdat$Var)^2,
         kurtdiff = (kurt - CVdat$kurt)^2) %>%
  group_by(model,leftout) %>%
  summarize(ssqe = sqrt(sum(MeanErrodiff)/length(MeanErrodiff)),
            ssqvar = sqrt(sum(Vardiff)/length(Vardiff)),
            ssqkurt = sqrt(sum(kurtdiff)/length(kurtdiff))) %>%
  mutate(measure2 = "RMSE") %>%
  pivot_longer(!c(model,measure2,leftout),names_to = "measure",values_to="value") %>%
  bind_cols(meanval = CVdat_meandat %>% select(value)) %>%
  mutate(vals = value...5/value...6)
models3 <- c("MK_U_RNplus","MK_U_RNminus",
             "MK_P_RNplus","MK_P_RNminus",
             "MK_FM_RNplus","MK_FM_RNminus",
             "MK_RNplus","MK_RNminus")
SumStat2_all <-CVdat2 %>%
  mutate(measure = factor(measure, levels = c("ssqvar","ssqe","ssqkurt"),
                          labels = c(expression(paste("circular ",sigma^2)),
                                     expression(paste("mean abs. deviation (radian)")),
                                     expression(paste("kurtosis"))))) %>%
  filter(model %in% models) %>%
  mutate(model = factor(model, levels = c("MK_RNminus","MK_RNplus","MK_FM_RNminus","MK_FM_RNplus",
                                          "MK_P_RNminus","MK_P_RNplus","MK_U_RNminus","MK_U_RNplus"))) %>%
  # mutate(model = factor(model, levels = models)) %>%
  filter(leftout %in% c("Full data set fitted","Left out: Ss 1","Left out: Ss 2")) %>%
  mutate(leftout = factor(leftout, levels = c("Left out: Ss 2","Left out: Ss 1","Full data set fitted"),
                          labels = c(
                            expression(paste("Left out: Ss 2")),
                                     expression(paste("Left out: Ss 1")),
                                     expression(paste("No Ss left out"))))) %>%
  mutate(combs = paste0(model,"_",leftout)) %>%
  mutate(modellimit = case_when(model %in% c("MK_RNminus","MK_RNplus") ~ "A",
                                model %in% c("MK_FM_RNminus","MK_FM_RNplus") ~ "F",
                                model %in% c("MK_U_RNminus","MK_U_RNplus") ~ "u",
                                model %in% c("MK_P_RNminus","MK_P_RNplus") ~ "P",
                                TRUE ~ "NA")) %>%
  mutate(modelRN = case_when(model %in% models[c(1:4)] ~ "Minus",
                             model %in% models[c(5:8)] ~ "Plus",
                             TRUE ~ "NA")) %>%
  mutate(type = "all")


SumStat <- SumStat2_all  %>%  mutate(model = factor(model, levels = models3))



Vars2_SumStat <- ggplot(SumStat %>% filter(measure2 != "Fit"), aes(x = leftout, y = vals)) +
  geom_bar(data = SumStat %>% filter(measure2 != "Fit"),stat="identity",
           position = "dodge",
           aes(x = leftout, y = vals, fill=model,linetype=model),
           color="black",width=0.9) +
  coord_flip()+
  scale_x_discrete(labels = c("Left out:\nSs 2", "Left out:\nSs 1", "Full data\nset fitted"))+
  facet_wrap(measure~.,ncol=3,labeller=label_parsed) +
  scale_y_continuous(limits  = c(0,0.15),name="normalized RMSD")+



  scale_fill_manual(values = c("#CC79A7","#f1dae7",
                               "#44AA99","#9ee6da",
                               "#E69F00","#f2d188",
                               "#332288","#b2a8e0"),
                     labels = c(
                      expression(paste("VP(",kappa,")U+")),
                      expression(paste("VP(",kappa,")U-")),
                      expression(paste("VP(",kappa,")P+")),

                      expression(paste("VP(",kappa,")P-")),
                      expression(paste("VP(",kappa,")F+")),

                      expression(paste("VP(",kappa,")F-")),
                      expression(paste("VP(",kappa,")A+")),
                      expression(paste("VP(",kappa,")A-"))




                    ), guide = guide_legend(reverse = TRUE)) +
  scale_linetype_manual(values = rep(c("solid","dashed"),4),
                        labels = c(
                          expression(paste("VP(",kappa,")U+")),
                          expression(paste("VP(",kappa,")U-")),
                          expression(paste("VP(",kappa,")P+")),

                          expression(paste("VP(",kappa,")P-")),
                          expression(paste("VP(",kappa,")F+")),

                          expression(paste("VP(",kappa,")F-")),
                          expression(paste("VP(",kappa,")A+")),
                          expression(paste("VP(",kappa,")A-"))




                        ), guide = guide_legend(reverse = TRUE)
                        )+
  #ggtitle("(N)RMSD")+
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=14),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text = element_text(size = 13),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.95,0.5),
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
Vars2_SumStat

# CP: Finalize

plot_grid(AICDev,errordistribution,Vars1,Vars2_SumStat,ncol=1,rel_heights = c(0.7,1.5,0.7,0.7), labels = "AUTO",
          align = 'v', axis = 'l')



ggsave("SetSize1_graphs_p1.png", units="cm", width=35, height=50, dpi=600)


# COWPLOT: Make Set size 1 graphs p2--------------------------------------------

# get all the data: Separate fits

SepSS <- readRDS("Fits/FitFull_sepSS.rds")


DeltaSepSS <- SepSS %>%
  filter(model %in% c("VPnosetsize","VPplusnosetsize")) %>%
  mutate(nparam = ifelse(model == "VPnosetsize",2,3)) %>%
  mutate(AIC = 2 * nparam + 2 * objective) %>%
  group_by(id,setsize,model) %>%
  arrange(objective) %>%
  slice(1) %>%
  mutate(deviance = 2 * objective) %>%
  select(id,setsize,model,deviance,leftout,nparam,AIC) %>%
  group_by(setsize,model) %>%
  summarize(mDev = mean(deviance),
            lengths = length(deviance)) %>%
  mutate(deltaDev = mDev - min(mDev)) %>%
  mutate(model = ifelse(model == "VPnosetsize","MK_RNminus","MK_RNplus")) %>%
  mutate(leftout = "1SepSS")

# get all the data: AIC for set sizes from full data set fitted

FullFitSS <- readRDS("Fits/FitFull_simSS.rds")


DeltaFullSS <- FullFitSS %>%
  filter(model %in% c("MK_RNminus","MK_RNplus")) %>%
  mutate(nparam = ifelse(model == "MK_RNminus",3,4)) %>%
  mutate(AIC = 2 * nparam + 2 * objective) %>%
  group_by(id,setsize,model) %>%
  arrange(objective) %>%
  slice(1) %>%
  mutate(deviance = 2 * objective) %>%

  group_by(setsize,model) %>%
  summarize(mDev = mean(deviance),
            lengths = length(deviance)) %>%
  mutate(deltaDev = mDev - min(mDev)) %>%
  mutate(leftout = "2FullSS")


dub <- bind_rows(DeltaSepSS,DeltaFullSS) %>%
  mutate(setsize = paste0("Ss ", setsize)) %>%
  mutate(leftout = factor(leftout,levels = c("2FullSS","1SepSS"),
                          labels = c("Full data\nset fitted","Set sizes\nfitted separately")))

Ppl <- dub %>% group_by(leftout,setsize) %>% slice(1)

onedub <- ggplot(dub %>% filter(leftout == "Full data\nset fitted"),
                 aes(x = model,y=deltaDev,group=model)) +
  geom_bar(aes(fill = model,linetype=model), color="black",stat="identity",  position = "dodge")+
  geom_text(data = Ppl %>% filter(leftout == "Full data\nset fitted"),
            aes(x = 2.5, y = 2,label = paste("N = ",lengths)))+

  facet_grid(leftout~setsize) +
  scale_y_continuous(limits=c(0,4),name=expression(paste(Delta, " Deviance")))+
  scale_x_discrete(labels = c(expression(paste("VP(",kappa,")A-")),
                              expression(paste("VP(",kappa,")A+"))))+
  coord_flip() +



  scale_fill_manual(values = c("#f1dae7","#CC79A7"),
                    guide = guide_legend(FALSE)) +
  scale_linetype_manual(values = c("dashed","solid"),
                        guide = guide_legend(FALSE)) +
  #ggtitle("(N)RMSD")+
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text = element_text(size = 13),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position ="none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")

twodub <- ggplot(dub %>% filter(leftout == "Set sizes\nfitted separately"),aes(x = model,y=deltaDev,group=model)) +
  geom_bar(aes(fill = model,linetype=model), color="black",stat="identity",  position = "dodge")+
  facet_grid(leftout~setsize) +
  geom_text(data = Ppl %>% filter(leftout == "Set sizes\nfitted separately"),
            aes(x = 2.5, y = 2,label = paste("N = ",lengths)))+

  scale_y_continuous(limits=c(0,4),name=expression(paste(Delta, " Deviance")))+
  scale_x_discrete(labels = c(expression(paste("VP(",kappa,")A-")),
                              expression(paste("VP(",kappa,")A+"))))+
  coord_flip() +



  scale_fill_manual(values = c("#afd8f0","#56B4E9"),
                    guide = guide_legend(FALSE)) +
  scale_linetype_manual(values = c("dashed","solid"),
                        guide = guide_legend(FALSE)) +
  #ggtitle("(N)RMSD")+
  theme(axis.text.y = element_text(size=14),
        axis.text.x = element_text(size=13),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=14),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size = 13),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position ="none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
DeltaDev <- plot_grid(onedub,twodub,ncol=1,rel_heights = c(0.96,1))



# # Agg data
#
# load("Data/ol17_prepared.rda")
# load("Data/vdb14_prepared.rda")
# load("Data/pratte17_prepared.rda")
#
# experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
#   filter(exp %in% "ol17_e1") %>%
#   mutate(id = exp) %>%
#   nest_by(id)
#
# # aggregate fit:
#
# FitsLOSSO <- readRDS("Fits/FitFull_AggLOSSOTraining_firstbatch.rds") %>% filter(exp == "ol17_e1") %>%
#   filter(model %in% c("MK_RNminus","MK_RNplus")) %>%
#   filter(leftout == 1) %>%
#   group_by(model) %>%
#   arrange(objective) %>%
#   slice(1)
#
# Fits <-  readRDS("Fits/FitFull_agg.rds") %>% filter(exp == "ol17_e1") %>%
#   filter(model %in% c("MK_RNminus","MK_RNplus")) %>%
#   group_by(model) %>%
#   arrange(objective) %>%
#   slice(1)
# allmodel <- NULL
# for (modeln in c("MK_RNminus","MK_RNplus")){
#
#   pars <- as.numeric(Fits %>%
#                        filter(model == modeln) %>%
#                        select(names(vwmvp::get_start_vp(modeln))) %>%
#                        ungroup() %>%
#                        select(-model))
#
#
#   parm <- vwmvp::prep_parameters(pars,modeln,c(1:8))
#
#   error_list <- prep_data(experimentfile %>% .$data %>% .[[1]] %>% mutate(id = exp)) %>%
#     .$datalist
#
#
#   out <- vector("numeric", length(c(1:8)))
#
#
#
#   for (i in c(1:8)){
#
#
#     out[[i]] <- vp_routine(precision = parm[[1]][i], parscont = c(parm[[2]],parm[[3]]),
#                            errors = error_list[[i]],
#                            model = modeln, predictOrLL = "LL")
#
#   }
#
#   fin <- data.frame(setsize = c(1:8),
#                     objective = out,
#                     model = modeln)
#
#   allmodel <- allmodel %>% bind_rows(fin)
#
# }
#
#
# VPA_ss <- bind_rows(readRDS("FitVPnosetsize/Agg/ol17_e1VPnosetsize_Full20Random2.rds"),
#                     readRDS("FitVPnosetsize/Agg/ol17_e1VPplusnosetsize_Full20Random2.rds"),
#                     readRDS("FitVPnosetsize/Agg/ol17_e1_VPnosetsize_Full20Random.rds"),
#                     readRDS("FitVPnosetsize/Agg/ol17_e1_VPplusnosetsize_Full20Random.rds")
# )
#
# DeltaAIC <- VPA_ss %>%
#   mutate(nparam = ifelse(model == "VPnosetsize",2,3)) %>%
#   mutate(Deviance = 2 * objective) %>%
#   group_by(setsize,model) %>%
#   arrange(objective) %>%
#   slice(1) %>%
#
#   select(setsize,model,objective,leftout,nparam,Deviance) %>%
#   mutate(leftout = 2) %>%
#   mutate(model = ifelse(model == "VPnosetsize","MK_RNminus","MK_RNplus"))
#
#
#
# variouswaysofpredicting <- bind_rows(as_tibble(allmodel) %>% mutate(leftout = 0)) %>%
#   mutate(nparam = ifelse(model == "MK_RNminus",3,4)) %>%
#   mutate(Deviance = 2 * objective) %>%
#   bind_rows(DeltaAIC)
#
#
# dubs <- variouswaysofpredicting %>%
#   group_by(leftout,setsize) %>%
#   mutate(diffDev = Deviance - min(Deviance)) %>%
#   mutate(leftout = factor(leftout,levels = c(0,2),
#                           labels = c("Simultaneous",
#                                      "Separate"))) %>%
#   mutate(setsize = paste0("Set size ", setsize))
#
#
#
#
# onedub <- ggplot(dubs %>% filter(leftout == "Simultaneous"),aes(x = model,y=diffDev,group=model)) +
#   geom_bar(aes(fill = model,linetype=model), color="black",stat="identity",  position = "dodge")+
#   geom_text(data = Ppl, aes(x = 8, y = 17,label = paste("N = ",totalppl)))+
#   facet_grid(leftout~setsize) +
#   scale_y_continuous(limits=c(0,30),name=expression(paste(Delta, " Deviance")))+
#   scale_x_discrete(labels = c(expression(paste("VP(",kappa,")A-")),
#                               expression(paste("VP(",kappa,")A+"))))+
#   coord_flip() +
#
#
#
#   scale_fill_manual(values = c("#f1dae7","#CC79A7"),
#                     guide = guide_legend(FALSE)) +
#   scale_linetype_manual(values = c("dashed","solid"),
#                         guide = guide_legend(FALSE)) +
#   #ggtitle("(N)RMSD")+
#   theme(axis.text.y = element_text(size=14),
#         axis.text.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         strip.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent"),
#         plot.title = element_text(hjust=0.5,face="plain"),
#         strip.text = element_text(size = 13),
#         legend.key = element_rect(colour = "transparent", fill = "transparent"),
#         legend.background = element_rect(fill = "transparent"),
#         legend.text = element_text(size = 12,face="plain"),
#         legend.position ="none",
#         legend.text.align = 0,
#         legend.key.height = unit(0.5,"line"),
#         legend.title = element_blank(),
#         legend.box = "vertical")
#
# twodub <- ggplot(dubs %>% filter(leftout == "Separate"),aes(x = model,y=diffDev,group=model)) +
#   geom_bar(aes(fill = model,linetype=model), color="black",stat="identity",  position = "dodge")+
#   facet_grid(leftout~setsize) +
#   scale_y_continuous(limits=c(0,30),name=expression(paste(Delta, " Deviance")))+
#   scale_x_discrete(labels = c(expression(paste("VP(",kappa,")A-")),
#                               expression(paste("VP(",kappa,")A+"))))+
#   coord_flip() +
#
#
#
#   scale_fill_manual(values = c("#afd8f0","#56B4E9"),
#                     guide = guide_legend(FALSE)) +
#   scale_linetype_manual(values = c("dashed","solid"),
#                         guide = guide_legend(FALSE)) +
#   #ggtitle("(N)RMSD")+
#   theme(axis.text.y = element_text(size=14),
#         axis.text.x = element_text(size=13),
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size=14),
#         panel.background = element_rect(fill = "transparent"),
#         panel.border = element_rect(colour = "black", fill=NA, size=0.5),
#         strip.background = element_rect(fill = "transparent"),
#         plot.background = element_rect(fill = "transparent"),
#         plot.title = element_text(hjust=0.5,face="plain"),
#         strip.text.x = element_blank(),
#         strip.text.y = element_text(size = 13),
#         legend.key = element_rect(colour = "transparent", fill = "transparent"),
#         legend.background = element_rect(fill = "transparent"),
#         legend.text = element_text(size = 12,face="plain"),
#         legend.position ="none",
#         legend.text.align = 0,
#         legend.key.height = unit(0.5,"line"),
#         legend.title = element_blank(),
#         legend.box = "vertical")
#  plot_grid(onedub,twodub,ncol=1,rel_heights = c(0.96,1))




exps <- "ol17_e1"
experiment <- exps
experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size") %>%
  filter(exp %in% exps)

expsetsize <- sort(unique(experimentfile %>% .$setsize))
setsizelabels <- paste("Set size ",expsetsize)
leftoutlabels <- c("Full data set fitted", paste("Left out: Ss",expsetsize))


experimentfile <- experimentfile %>%
  filter(exp == experiment)%>%
  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels),
         id = as.character(factor(id)),
         Predicted = NA)


expanded <- experimentfile %>% expand(exp,id,setsize)


PredictionSS<- bind_rows(readRDS(paste0("Prediction/prediction_sepSS_full_",experiment,"_agg.rds")),
                         readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds"))) %>%
                           filter(leftout == 0) %>%
                           mutate(leftout = "Full data set fitted")%>%
  filter(model %in% c("MK_RNminus","MK_RNplus","VPnosetsize","VPplusnosetsize")) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_RNplus",
                                         "VPnosetsize","VPplusnosetsize"
  ))) %>%

  mutate(setsize = factor(setsize,
                          levels = expsetsize,
                          labels =  setsizelabels))


errordistribution <- ggplot(data = experimentfile,
                            aes(x = error_0)) +

  geom_histogram(aes(y = ..density..),color="darkgrey",fill="#edeaea",bins=60) +
  coord_cartesian(xlim = c(-pi/2,pi/2))+

  scale_y_continuous(limits=c(0,max(PredictionSS %>%

                                      .$prediction,na.rm=T) + 0.2)) +
  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(-pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_line(data = PredictionSS,

            aes(x = data,y=prediction,linetype=model,color=model,group=model), size = 0.75) +
  facet_wrap(setsize~.,ncol=4) +

  scale_linetype_manual(name = NULL,values = c("dashed","solid","dashed","solid"),
                        labels = c(


                          expression(paste("VP(",kappa,")A-, full data set fitted")),
                          expression(paste("VP(",kappa,")A+, full data set fitted")),
                          expression(paste("VP(",kappa,")A-, set sizes fitted separately")),
                          expression(paste("VP(",kappa,")A+, set sizes fitted separately"))
                        ), guide = guide_legend(ncol=2)) +
  scale_color_manual(name=NULL,
                     values =  c("#CC79A7","#CC79A7",
                                 "#56B4E9","#56B4E9"),
                     labels = c(


                       expression(paste("VP(",kappa,")A-, full data set fitted")),
                       expression(paste("VP(",kappa,")A+, full data set fitted")),
                       expression(paste("VP(",kappa,")A-, set sizes fitted separately")),
                       expression(paste("VP(",kappa,")A+, set sizes fitted separately"))

                     )
  )+

theme(axis.text.x = element_text(size=11),
      plot.title = element_text(size = 14,hjust=0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size=12),
      axis.title.y = element_blank(),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      strip.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent"),
      strip.text = element_text(size = 12),
      legend.position = c(0.65,0.3),
      legend.text = element_text(size = 13),
      legend.direction = "vertical",
      legend.background = element_rect(fill="white",color="lightgrey"),
      legend.key = element_rect(colour = "transparent", fill = "white"),
      legend.title = element_blank(),
      legend.key.width = unit(3,"line"))


# Make summary stat

PredictionSS<- bind_rows(readRDS(paste0("Prediction/prediction_sepSS_full_",experiment,"_agg.rds")),
                         readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds"))) %>%
  filter(leftout == 0) %>%
  mutate(leftout = "Full data set fitted")%>%
  filter(model %in% c("MK_RNminus","MK_RNplus","VPnosetsize","VPplusnosetsize")) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_RNplus",
                                         "VPnosetsize","VPplusnosetsize"
  )))

Prediction <- PredictionSS
models <- unique(PredictionSS$model)
CV <- NULL
for(SS in unique(Prediction %>% .$setsize)){

  CircVars <- NULL

  for(modeln in models){
    testNA <- any(is.na(Prediction %>%
                          filter(setsize == SS) %>%
                          filter(model == modeln) %>%
                          .$normprediction) == TRUE)
    if (testNA){

      CircVar <- NA
      CircSD <- NA
      MAE <- NA

    } else {

      samples <- sample(size=6e4,Prediction %>%
                          filter(setsize == SS) %>%
                          filter(model == modeln) %>%
                          .$data,replace=TRUE,
                        prob =  Prediction %>%
                          filter(setsize == SS) %>%
                          filter(model == modeln) %>%
                          .$normprediction)
      CircVar <- circular::var.circular(samples)
      CircSD <- circular::sd.circular(samples)
      MAE <- mean(abs(samples))
      kurtosis <- vwmvp::circkurtosis(samples)

    }

    out <- as_tibble(data.frame(setsize = SS,
                                model = modeln,
                                Var = CircVar,
                                SD = CircSD,
                                meanErr = MAE,
                                kurt = kurtosis))

    CircVars <- CircVars %>%  bind_rows(out)
  }


  exps <- "ol17_e1"
  experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
    rename(setsize = "set_size") %>%
    filter(exp %in% exps)

  expdata <- experimentfile %>%
    filter(exp == experiment) %>%
    filter(setsize == SS) %>%
    .$error_0



  cv <- circular::var.circular(expdata)
  csd <- circular::sd.circular(expdata)
  mae <- mean(abs(expdata))
  kurtosis_data <- vwmvp::circkurtosis(expdata)

  out2 <- as_tibble(data.frame(setsize = SS,
                               model = "Data",
                               Var = cv,
                               SD = csd,
                               meanErr = mae,
                               kurt = kurtosis_data))
  CircVars <- CircVars %>%  bind_rows(out2)
  CV <- CV %>% bind_rows(CircVars)
}


# CV <- CV %>% mutate(model = factor(model,
#                                    levels = c("Data",models)))


CVdat <- CV %>% filter(model == "Data") %>%  slice(rep(1:n(),each=4))
CVdat_meandat <- CVdat %>% group_by(model) %>% summarize(meanE = mean(meanErr),
                                                         meanV = mean(Var),
                                                         meank = mean(kurt)) %>%
  pivot_longer(!c(model),names_to = "measure",values_to="value") %>%
  slice(rep(1:n(),4))


CVdat2 <- CV %>% filter(model != "Data") %>%
  arrange(setsize,model) %>%
  mutate(MeanErrodiff = (meanErr - CVdat$meanErr)^2,
         Vardiff = (Var - CVdat$Var)^2,
         kurtdiff = (kurt - CVdat$kurt)^2) %>%
  group_by(model) %>%
  summarize(ssqe = sqrt(sum(MeanErrodiff)/length(MeanErrodiff)),
            ssqvar = sqrt(sum(Vardiff)/length(Vardiff)),
            ssqkurt = sqrt(sum(kurtdiff)/length(kurtdiff))) %>%
  mutate(measure2 = "RMSE") %>%
  pivot_longer(!c(model,measure2),names_to = "measure",values_to="value") %>%
  mutate(value = value/CVdat_meandat$value)




CVplot <- CV %>% select(-SD) %>%
  pivot_longer(!c(model,setsize),names_to = "measure",values_to="value") %>%
  mutate(measure = factor(measure, levels = c("Var","meanErr","kurt"),
                          labels = c(expression(paste("circular ",sigma^2)),
                                     expression(paste("mean abs. deviation (radian)")),
                                     expression(paste("kurtosis"))))) %>%
  mutate(model = factor(model, levels = c("Data",c("MK_RNminus","MK_RNplus","VPnosetsize","VPplusnosetsize"))))

Vars1 <- ggplot(CVplot %>% filter(model == "Data") , aes(x = setsize, y = value)) +
  geom_hpline(aes(linetype = model),size=1.5,width=1) +
  geom_point(data = CVplot %>% filter(model != "Data"),
             aes(x = setsize, y = value,shape = model, fill = model),
             position = position_dodge(0.6),size=3.5,alpha=0.7,color="black") +
  scale_x_continuous(breaks=c(1:8),labels=c(c(1:8)),name="Set size") +
  facet_grid(.~measure,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25))+
  scale_shape_manual(values = rep(c(21,22),2),
                     labels = c(
                       expression(paste("VP(",kappa,")A-, full data set fitted")),
                       expression(paste("VP(",kappa,")A+, full data set fitted")),
                       expression(paste("VP(",kappa,")A-, set sizes fitted separately")),
                       expression(paste("VP(",kappa,")A+, set sizes fitted separately"))

                     ), guide = guide_legend(ncol=2,order=1)
  ) +

  scale_fill_manual(values =  c("#f1dae7","#CC79A7", "#afd8f0","#56B4E9"),
                    labels = c(
                      expression(paste("VP(",kappa,")A-, full data set fitted")),
                      expression(paste("VP(",kappa,")A+, full data set fitted")),
                      expression(paste("VP(",kappa,")A-, set sizes fitted separately")),
                      expression(paste("VP(",kappa,")A+, set sizes fitted separately"))),
                    guide = guide_legend(ncol=2,order=1)
  )+
  theme(axis.text = element_text(size=11),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 14),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="lightgrey"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.25,0.72),
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
Vars1


plot_grid(DeltaDev,errordistribution,Vars1,ncol=1,labels="AUTO",rel_heights = c(1,1.25,0.8))

#ggsave("SetSize1_graphs_p2.png", units="cm", width=30, height=30, dpi=600)