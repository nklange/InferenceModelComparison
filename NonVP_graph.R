
source("CompareModels_NonVP.R") # for quantitative fit plots

library(vwmvp)
library("grid")
library("gridExtra")
library("cowplot")
library(ungeviz)


# Stability by points ----------------------------------------------------------


load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,
                            data_pratte17)



FittedFull <- readRDS("Fits/FitFull.rds")
numberid <- length(unique(FittedFull$cvid))
exclall <- c()


models <- c("MK_RNminus",
            "MK_FM_RNplus","MK_RNplus","MK_U_RNplus","MK_P_RNplus",
            "MK_P_RNminus","MK_U_RNminus",
            "MK_FM_RNminus",
            "SA_F_RNplus","SA_F_RNminus",
            "SA_P_RNplus", "SA_U_RNplus","EP_P_RNplus","EP_U_RNplus",
            "EP_RNplus","SA_RNplus",
            "EP_RNminus", "EP_FM_RNplus","EP_FM_RNminus",
            "EP_P_RNminus","EP_U_RNminus","SA_P_RNminus","SA_U_RNminus")

modelsStaborder <- c("MK_RNminus","MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_FM_RNplus","MK_RNplus","MK_P_RNplus","MK_U_RNplus",

            "EP_RNminus","EP_FM_RNminus","EP_P_RNminus","EP_U_RNminus",
            "EP_RNplus","EP_FM_RNplus","EP_P_RNplus","EP_U_RNplus",


            "SA_F_RNminus","SA_P_RNminus","SA_U_RNminus",
            "SA_RNplus","SA_F_RNplus","SA_P_RNplus", "SA_U_RNplus")

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
                        levels = modelsStaborder)) %>%
  mutate(modelnum = as.numeric(model)) %>%
  filter(model %in% models) %>%
  mutate(modelposition = model) %>%
  mutate(modelposition = factor(modelposition,
                                levels =  modelsStaborder)) %>%
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

colorin <- ggplotColours(n = 23)

plotcolors <- c(rep("#009E73",8),rep("#E69F00",8), rep("#56B4E9",7))#colorin[unique(thresholdmodel$modelnum)]

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

                        expression(paste("EP(",kappa,")A-")),
                        expression(paste("EP(",kappa,")F-")),
                        expression(paste("EP(",kappa,")P-")),
                        expression(paste("EP(",kappa,")U-")),

                        expression(paste("EP(",kappa,")A+")),
                        expression(paste("EP(",kappa,")F+")),
                        expression(paste("EP(",kappa,")P+")),
                        expression(paste("EP(",kappa,")U+")),

                        expression(paste("SA(",kappa,")F-")),
                        expression(paste("SA(",kappa,")P-")),
                        expression(paste("SA(",kappa,")U-")),

                        expression(paste("SA(",kappa,")A+")),
                        expression(paste("SA(",kappa,")F+")),
                        expression(paste("SA(",kappa,")P+")),
                        expression(paste("SA(",kappa,")U+")),

                       "total")

  ) +
  theme(axis.text.x = element_text(size=11,angle=90,vjust=0.5),
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

# AIC, LOSSO-CV ----------------------------------------------------------------

Quant <- plot_grid(ModelCompTotal,ModelCompSs,labels=c("B","C"),nrow=1,rel_widths = c(1,1),align = 'h', axis = 't',scale = 0.98)

# Summary statistics -----------------------------------------------------------


manuscriptmodelorder2 <- c("SA_RNplus", "EP_RNplus","EP_RNminus", "SA_U_RNplus","SA_U_RNminus",
                    "EP_FM_RNplus","EP_FM_RNminus","EP_U_RNminus","SA_F_RNminus","SA_F_RNplus",
                    "EP_U_RNplus", "SA_P_RNplus","SA_P_RNminus","MK_U_RNminus","EP_P_RNminus",
                    "MK_RNminus","EP_P_RNplus",
                    "MK_P_RNminus","MK_FM_RNplus","MK_FM_RNminus","MK_RNplus","MK_P_RNplus",
                    "MK_U_RNplus")

TotCVind1 <- NULL
TotCVind2 <- NULL
for (experiment in c("ol17_e1","WM4", "RTT12", "ZL8", "BCH9", "AVA11",
                       "AA12a", "AA12b", "VSCGM21a", "VSCGM21b", "VSCGM21c","pratte17")){


CVind1 <- readRDS(paste0("SummaryStat/SumStat_altmod_",experiment,"_ind.rds"))
CVind2 <- readRDS(paste0("SummaryStat/SumStat_",experiment,"_ind.rds"))

TotCVind1 <- TotCVind1 %>% bind_rows(CVind1)
TotCVind2 <- TotCVind2 %>% bind_rows(CVind2)

}



Sslevels <- sort(as.character(levels(unique(TotCVind2$setsize))))

CVind <- bind_rows(TotCVind1 %>% mutate(setsize = factor(setsize,levels=c(1:8),labels=c(1:8))) %>% filter(model != "Data"),
TotCVind2 %>% mutate(setsize = factor(setsize,levels=Sslevels,labels = c(1:8)))) %>%
  filter(model %in% c("Data",models))



CVdat_meandat <- CVind %>% filter(id != "Av") %>% filter(model == "Data") %>%
  group_by(exp,id,model) %>% summarize(ssqe = mean(meanErr),
                                       ssqvar = mean(Var),
                                       ssqkurt = mean(kurt)) %>%
  pivot_longer(!c(exp,id,model),names_to = "measure",values_to="value") %>%
  slice(rep(1:n(),each=length(models))) %>%
  arrange(exp,id,measure) %>% ungroup()


CVdat <- CVind %>% filter(model == "Data") %>%  slice(rep(1:n(),each=length(models)))


predCV <- CVind %>% filter(model != "Data") %>%  filter(id != "Av") %>%
  arrange(id,setsize,model)
obsCV <- CVdat %>% arrange(id,model,setsize)



NRMSD_ind <- predCV %>%
  mutate(MeanErrodiff = (meanErr - obsCV$meanErr)^2,
         Vardiff = (Var - obsCV$Var)^2,
         kurtdiff = (kurt - obsCV$kurt)^2) %>%
  group_by(exp,id,model) %>%
  summarize(ssqe = sqrt(sum(MeanErrodiff)/length(MeanErrodiff)),
            ssqvar = sqrt(sum(Vardiff)/length(Vardiff)),
            ssqkurt = sqrt(sum(kurtdiff)/length(kurtdiff))) %>%
  mutate(measure2 = "RMSE") %>%
  pivot_longer(!c(id,exp,model,measure2),names_to = "measure",values_to="value")  %>%
  arrange(exp,id,measure) %>%
  ungroup() %>%
  mutate(value = value/CVdat_meandat$value)




  BP <- NRMSD_ind %>%  filter(model %in% models) %>%
    mutate(measure = factor(measure, levels = c("ssqvar","ssqe","ssqkurt"),
                            labels = c(expression(paste("circular ",sigma^2)),
                                       expression(paste("mean abs. deviation (radian)")),
                                       expression(paste("kurtosis"))))) %>%
    mutate(model = factor(model, levels = manuscriptmodelorder2)) %>%
    mutate(modeldist = case_when(grepl("MK", model, fixed = TRUE) ~ "VP",
                                 grepl("EP",model,fixed=TRUE) ~ "EP",
                                 TRUE ~ "SA"))

NRMSD <-  ggplot(BP , aes(x = model, y = value)) +
    geom_boxplot(data = BP,
                 aes(x = model, y = value,fill=modeldist),outlier.shape = 1) +

    coord_flip()+
    scale_fill_manual(values=c("#E69F00", "#56B4E9", "#009E73"))+

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


    facet_wrap(.~measure,ncol=3,labeller=label_parsed) +
    scale_y_continuous(limits  = c(0,1),name="normalized RMSD",expand = c(0.002, 0.002))+
    # scale_fill_manual(values =  c("#CC79A7","#f1dae7",
    #                               "#44AA99","#9ee6da",
    #                               "#E69F00","#f2d188",
    #                               "#332288","#b2a8e0"
    #
    #
    #
    # )) +
    #scale_linetype_manual(values = rep(c("solid","dashed"),4))+
    #ggtitle("(N)RMSD")+
    theme(axis.text.y = element_text(size=10,hjust=0),
          axis.text.x = element_text(size=12),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=12),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          plot.title = element_text(hjust=0.5,face="plain"),
          strip.text = element_text(size = 13),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12,face="plain"),
          legend.position = "none",
          legend.text.align = 0,
          legend.key.height = unit(0.5,"line"),
          legend.title = element_blank(),
          legend.box = "vertical")


# CP: Finalize -----------------------------------------------------------------

plot_grid(datastabpoints,Quant,labels = c("A",NA),ncol=1,
          rel_heights = c(0.7,1))

#ggsave("AltMod1.png", units="cm", width=35, height=40, dpi=600)

plot_grid(NRMSD,ExpG,labels = c("A","B"),ncol=1,
          rel_heights = c(0.4,1))

#ggsave("AltMod2.png", units="cm", width=35, height=45, dpi=600)

# Shen & Ma (2019) Factorial Importance factors --------------------------------

capLimit <- "_FM_"

FIMVPEP <- function(capLimit){
if(capLimit == "_FM_"){
  LABELCAP <- c("A","B")
  limLabel<- "Fixed memory capacity"
} else if (capLimit == "_P_") {
  limLabel <- "Variable capacity limit (Poisson-dist.)"
  LABELCAP <- c("C","D")

} else {
  LABELCAP <- c("E","F")
  limLabel <- "Variable capacity limit (Unif.-dist.)"

}

Capmodels <- c("MK_RNminus","MK_RNplus",
                "MK_U_RNminus","MK_U_RNplus",
                "MK_P_RNminus","MK_P_RNplus",
                "MK_FM_RNminus", "MK_FM_RNplus",
                "EP_P_RNminus","EP_P_RNplus",
                "EP_RNplus", "EP_RNminus",
                "EP_FM_RNplus","EP_FM_RNminus",
                "EP_U_RNplus","EP_U_RNminus"
               )
models <- c("MK_RNminus","MK_RNplus", "EP_RNplus", "EP_RNminus",Capmodels[grepl(capLimit,Capmodels)])

FittedFull <- readRDS("Fits/FitFull.rds") %>%
  filter(model %in% models) %>%
  group_by(exp,id,model) %>%
  arrange(objective) %>%
  slice(1) %>%
  select(exp,id,model,AIC) %>%
  mutate(FactorRN = ifelse(grepl("RNminus",model),"Absent","Present")) %>%
  mutate(FactorVP = ifelse(grepl("MK",model),"Present","Absent")) %>%
  mutate(FactorCap = case_when(grepl("_P_",model) ~ "Present",
                               grepl("_U_",model) ~ "Present",
                               grepl("_FM_",model) ~ "Present",
                               TRUE ~ "Absent")) %>%
  mutate(type = case_when(model == "EP_RNminus" ~ 1,
                          model %in% c("MK_U_RNplus","MK_FM_RNplus","MK_P_RNplus") ~ 3,
                          TRUE ~ 2)) %>%
  group_by(exp,id) %>% arrange(exp,id,type,FactorRN,FactorVP,FactorCap)

KIDKOD <- FittedFull %>% group_by(exp,id) %>% mutate(KID = AIC[1] - AIC,
                                  KOD = AIC - AIC[length(AIC)],
                                  best = AIC - min(AIC))


LFLR <- KIDKOD %>% mutate(expon = exp(-best)) %>%
  summarize(RNFactor = log(mean(expon[5:8])/mean(expon[1:4])),
            VPFactor = log(mean(expon[c(3,4,7,8)])/mean(expon[c(1,2,5,6)])),
            CapFactor = log(mean(expon[c(2,4,6,8)])/mean(expon[c(1,3,5,7)])),
            RNVP =  log(mean(expon[c(7,8)])/mean(expon[c(1:2)])),
            VPCap = log(mean(expon[c(4,8)])/mean(expon[c(1,5)])),
            RNCap = log(mean(expon[c(6,8)])/mean(expon[c(1,3)]))) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value")%>%
  mutate(FIM = "LFLR")

KID <- KIDKOD %>% summarize(RNFactor = KID[5],
                     VPFactor = KID[3],
                     CapFactor = KID[2],
                     RNVP = KID[7],
                     VPCap = KID[4],
                     RNCap = KID[6]) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
  mutate(FIM = "KID")


KOD <- KIDKOD %>% summarize(RNFactor = KOD[4],
                            VPFactor = KOD[6],
                            CapFactor = KOD[7],
                            RNVP = KOD[2],
                            VPCap = KOD[5],
                            RNCap = KOD[3]) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
  mutate(FIM = "KOD")

FIM <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,measure,FIM) %>%
  summarize(meanFIM = mean(value),
            seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "VPFactor",
                                             "CapFactor",
                                             "RNVP",
                                             "VPCap",
                                             "RNCap"),
                          labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)")))

FIMind <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,id,measure,FIM) %>%
  summarize(meanFIM = mean(value),
            seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "VPFactor",
                                             "CapFactor",
                                             "RNVP",
                                             "VPCap",
                                             "RNCap"),
                          labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)")))


expgroupFULL<-ggplot(data=FIMind, aes(x = exp,y=meanFIM,fill=measure)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
  geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
  geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
  geom_boxplot(data = FIMind,aes(x= exp,y=meanFIM,fill=measure),color="darkgrey",
               outlier.size=1,position=position_dodge(.8),alpha=0.25) +


  geom_errorbar(data = FIM,aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                width=0,color="black",size=0.5)+
  geom_point(data=FIM, aes(x = exp,y=meanFIM,group=measure,fill=measure),
             stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +

  #scale_fill_manual(values = c("#0072B2", "#CC79A7"))+
  #scale_color_manual(values = c("#0072B2", "#CC79A7"))+
  scale_y_continuous(limits = c(-20,55),breaks=seq(-10,50,10))+
  #  geom_text(aes(label = starred))+
  facet_grid(FIM~Lim) +

  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_text(size = 10),
        strip.text.y = StripLabel,
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = LegendInfo,
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")

full <- ggplot(data=FIMind %>%  filter(measure%in%c("RN","Cap","RNCap")),
                       aes(x = exp,y=meanFIM,fill=measure)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
  geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
  geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
  geom_boxplot(data = FIMind %>%  filter(measure%in%c("RN","Cap","RNCap")) ,
               aes(x= exp,y=meanFIM,fill=measure),color="darkgrey",
               outlier.size=1,position=position_dodge(.8),alpha=0.25) +


  geom_errorbar(data = FIM %>%  filter(measure%in%c("RN","Cap","RNCap")),
                aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                width=0,color="black",size=0.5)+
  geom_point(data=FIM %>%   filter(measure%in%c("RN","Cap","RNCap")),
             aes(x = exp,y=meanFIM,group=measure,fill=measure),
             stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +

  scale_fill_manual(values = c("#0072B2", "#CC79A7","#009E73"),
                    labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
  #scale_color_manual(values = c("#0072B2", "#CC79A7"))+
  scale_y_continuous(limits = c(-20,200),breaks=seq(0,200,50))+
  #  geom_text(aes(label = starred))+
  facet_grid(FIM~Lim) +

  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_text(size = 10),
        strip.text.y = StripLabel,
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = LegendInfo,
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")


filtered <- ggplot(data=FIMind %>%  filter(measure%in%c("RN","Cap","RNCap")),
       aes(x = exp,y=meanFIM,fill=measure)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
  geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
  geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
  geom_boxplot(data = FIMind %>%  filter(measure%in%c("RN","Cap","RNCap")) ,
               aes(x= exp,y=meanFIM,fill=measure),color="darkgrey",
               outlier.size=1,position=position_dodge(.8),alpha=0.25) +


  geom_errorbar(data = FIM %>%  filter(measure%in%c("RN","Cap","RNCap")),
                aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                width=0,color="black",size=0.5)+
  geom_point(data=FIM %>%   filter(measure%in%c("RN","Cap","RNCap")),
             aes(x = exp,y=meanFIM,group=measure,fill=measure),
             stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +

  scale_fill_manual(values = c("#0072B2", "#CC79A7","#009E73"),
                    labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
  #scale_color_manual(values = c("#0072B2", "#CC79A7"))+
  scale_y_continuous(limits = c(-20,55),breaks=seq(-10,50,10))+
  #  geom_text(aes(label = starred))+
  facet_grid(FIM~Lim) +

  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_text(size = 10),
        strip.text.y = StripLabel,
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = LegendInfo,
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")

FIMtotal <- bind_rows(KID,KOD,LFLR) %>% group_by(measure,FIM) %>%
  summarize(meanFIM = mean(value),
            seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "VPFactor",
                                             "CapFactor",
                                             "RNVP",
                                             "VPCap",
                                             "RNCap"),
                          labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)")))

totfiltered <- ggplot(FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")), aes(x = 1,y=meanFIM,group=measure,fill=measure)) +
  scale_x_continuous(breaks=1,label="         overall") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
  geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
  geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
  geom_boxplot(data = FIMind%>%  filter(measure%in%c("RN","Cap","RNCap")),aes(x= 1,y=meanFIM,fill=measure),color="darkgrey",
               outlier.size=1,position=position_dodge(.8),alpha=0.25) +


  geom_errorbar(data = FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")),aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                width=0,color="black",size=0.5)+
  geom_point(data=FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")), aes(x = 1,y=meanFIM,group=measure,fill=measure),
             stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00","#009E73"),
                    labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
  scale_y_continuous(limits = c(-20,55),breaks=seq(-10,50,10))+
  #  geom_text(aes(label = starred))+
  facet_grid(FIM~Lim) +

  theme(axis.text.y = element_text(size=10),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9,color="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_text(size = 10,color="white"),
        strip.text.y = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")

totfull <- ggplot(FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")), aes(x = 1,y=meanFIM,group=measure,fill=measure)) +
  scale_x_continuous(breaks=1,label="         overall") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
  geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
  geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
  geom_boxplot(data = FIMind%>%  filter(measure%in%c("RN","Cap","RNCap")),aes(x= 1,y=meanFIM,fill=measure),color="darkgrey",
               outlier.size=1,position=position_dodge(.8),alpha=0.25) +


  geom_errorbar(data = FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")),aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                width=0,color="black",size=0.5)+
  geom_point(data=FIMtotal%>%  filter(measure%in%c("RN","Cap","RNCap")), aes(x = 1,y=meanFIM,group=measure,fill=measure),
             stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +
  scale_fill_manual(values = c("#0072B2", "#D55E00","#009E73"),
                    labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
  scale_y_continuous(limits = c(-20,200),breaks=seq(0,200,50))+
  #  geom_text(aes(label = starred))+
  facet_grid(FIM~Lim) +

  theme(axis.text.y = element_text(size=10),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9,color="black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text.x = element_text(size = 10,color="white"),
        strip.text.y = element_blank(),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")


full <- plot_grid(totfull,full,rel_widths=c(0.3,1),nrow=1)
filtered <-  plot_grid(totfiltered,filtered,rel_widths=c(0.3,1),nrow=1)
plot_grid(full,filtered,nrow=1,rel_widths = c(1,1),labels=LABELCAP,scale=.95)

}
plot_grid(FIMVPEP("_FM_"),FIMVPEP("_P_"),FIMVPEP("_U_"),nrow=3,rel_heights=c(1,1,1))

#ggsave("FIM_InteractionCapacityRN.png", units="cm", width=35, height=45, dpi=600)


FIMVPEP2 <- function(capLimit){
  if(capLimit == "_FM_"){
    LABELCAP <- c("A","B")
    limLabel<- "Fixed memory capacity"
  } else if (capLimit == "_P_") {
    limLabel <- "Variable capacity limit (Poisson-dist.)"
    LABELCAP <- c("C","D")

  } else {
    LABELCAP <- c("E","F")
    limLabel <- "Variable capacity limit (Unif.-dist.)"

  }

  Capmodels <- c("MK_RNminus","MK_RNplus",
                 "MK_U_RNminus","MK_U_RNplus",
                 "MK_P_RNminus","MK_P_RNplus",
                 "MK_FM_RNminus", "MK_FM_RNplus",
                 "EP_P_RNminus","EP_P_RNplus",
                 "EP_RNplus", "EP_RNminus",
                 "EP_FM_RNplus","EP_FM_RNminus",
                 "EP_U_RNplus","EP_U_RNminus"
  )
  models <- c("MK_RNminus","MK_RNplus", "EP_RNplus", "EP_RNminus",Capmodels[grepl(capLimit,Capmodels)])

  FittedFull <- readRDS("Fits/FitFull.rds") %>%
    filter(model %in% models) %>%
    group_by(exp,id,model) %>%
    arrange(objective) %>%
    slice(1) %>%
    select(exp,id,model,AIC) %>%
    mutate(FactorRN = ifelse(grepl("RNminus",model),"Absent","Present")) %>%
    mutate(FactorVP = ifelse(grepl("MK",model),"Present","Absent")) %>%
    mutate(FactorCap = case_when(grepl("_P_",model) ~ "Present",
                                 grepl("_U_",model) ~ "Present",
                                 grepl("_FM_",model) ~ "Present",
                                 TRUE ~ "Absent")) %>%
    mutate(type = case_when(model == "EP_RNminus" ~ 1,
                            model %in% c("MK_U_RNplus","MK_FM_RNplus","MK_P_RNplus") ~ 3,
                            TRUE ~ 2)) %>%
    group_by(exp,id) %>% arrange(exp,id,type,FactorRN,FactorVP,FactorCap)

  KIDKOD <- FittedFull %>% group_by(exp,id) %>% mutate(KID = AIC[1] - AIC,
                                                       KOD = AIC - AIC[length(AIC)],
                                                       best = AIC - min(AIC))


  LFLR <- KIDKOD %>% mutate(expon = exp(-best)) %>%
    summarize(RNFactor = log(mean(expon[5:8])/mean(expon[1:4])),
              VPFactor = log(mean(expon[c(3,4,7,8)])/mean(expon[c(1,2,5,6)])),
              CapFactor = log(mean(expon[c(2,4,6,8)])/mean(expon[c(1,3,5,7)])),
              RNVP =  log(mean(expon[c(7,8)])/mean(expon[c(1:2)])),
              VPCap = log(mean(expon[c(4,8)])/mean(expon[c(1,5)])),
              RNCap = log(mean(expon[c(6,8)])/mean(expon[c(1,3)]))) %>%
    pivot_longer(!c(id,exp),names_to = "measure",values_to="value")%>%
    mutate(FIM = "LFLR")

  KID <- KIDKOD %>% summarize(RNFactor = KID[5],
                              VPFactor = KID[3],
                              CapFactor = KID[2],
                              RNVP = KID[7],
                              VPCap = KID[4],
                              RNCap = KID[6]) %>%
    pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
    mutate(FIM = "KID")


  KOD <- KIDKOD %>% summarize(RNFactor = KOD[4],
                              VPFactor = KOD[6],
                              CapFactor = KOD[7],
                              RNVP = KOD[2],
                              VPCap = KOD[5],
                              RNCap = KOD[3]) %>%
    pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
    mutate(FIM = "KOD")

  FIM <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,measure,FIM) %>%
    summarize(meanFIM = mean(value),
              seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
    mutate(measure = factor(measure,levels = c("RNFactor",
                                               "VPFactor",
                                               "CapFactor",
                                               "RNVP",
                                               "VPCap",
                                               "RNCap"),
                            labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
    mutate(Lim = limLabel) %>%
    mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                        labels = c("KID\n(base model + factor)",
                                   "KOD\n(full model - factor)",
                                   "2*LFLR\n(contrasting models with\nand without factor)")))

  FIMind <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,id,measure,FIM) %>%
    summarize(meanFIM = mean(value),
              seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
    mutate(measure = factor(measure,levels = c("RNFactor",
                                               "VPFactor",
                                               "CapFactor",
                                               "RNVP",
                                               "VPCap",
                                               "RNCap"),
                            labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
    mutate(Lim = limLabel) %>%
    mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                        labels = c("KID\n(base model + factor)",
                                   "KOD\n(full model - factor)",
                                   "2*LFLR\n(contrasting models with\nand without factor)")))



  full <- ggplot(data=FIMind ,
                 aes(x = exp,y=meanFIM,fill=measure)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
    geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
    geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
    geom_boxplot(data = FIMind  ,
                 aes(x= exp,y=meanFIM,fill=measure),color="darkgrey",
                 outlier.size=1,position=position_dodge(.8),alpha=0.25) +


    geom_errorbar(data = FIM ,
                  aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                  width=0,color="black",size=0.5)+
    geom_point(data=FIM ,
               aes(x = exp,y=meanFIM,group=measure,fill=measure),
               stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +

    scale_fill_manual(values = c("#0072B2","#D55E00", "#CC79A7",
                                 "#F0E442","#999999",
                                 "#009E73"),
                      labels = c("RN","VP","Cap Limit",
                                 "RN + VP", "VP + Cap Limit",
                                 "RN + Cap Limit"))+
    # #scale_color_manual(values = c("#0072B2", "#CC79A7"))+
    scale_y_continuous(limits = c(-20,200),breaks=seq(0,200,50))+
    #  geom_text(aes(label = starred))+
    facet_grid(FIM~Lim) +

    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          plot.title = element_text(hjust=0.5,face="plain"),
          strip.text.x = element_text(size = 10),
          strip.text.y = StripLabel,
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12,face="plain"),
          legend.position = "right",
          legend.text.align = 0,
          legend.key.height = unit(0.5,"line"),
          legend.title = element_blank(),
          legend.box = "vertical")


  filtered <- ggplot(data=FIMind,
                     aes(x = exp,y=meanFIM,fill=measure)) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
    geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
    geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
    geom_boxplot(data = FIMind %>%  filter(measure%in%c("RN","Cap","RNCap")) ,
                 aes(x= exp,y=meanFIM,fill=measure),color="darkgrey",
                 outlier.size=1,position=position_dodge(.8),alpha=0.25) +


    geom_errorbar(data = FIM ,
                  aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                  width=0,color="black",size=0.5)+
    geom_point(data=FIM,
               aes(x = exp,y=meanFIM,group=measure,fill=measure),
               stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +

    # scale_fill_manual(values = c("#0072B2", "#CC79A7","#009E73"),
    #                   labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
    #scale_color_manual(values = c("#0072B2", "#CC79A7"))+
    scale_y_continuous(limits = c(-20,55),breaks=seq(-10,50,10))+
    #  geom_text(aes(label = starred))+
    facet_grid(FIM~Lim) +

    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          plot.title = element_text(hjust=0.5,face="plain"),
          strip.text.x = element_text(size = 10),
          strip.text.y = StripLabel,
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12,face="plain"),
          legend.position = LegendInfo,
          legend.text.align = 0,
          legend.key.height = unit(0.5,"line"),
          legend.title = element_blank(),
          legend.box = "vertical")

  FIMtotal <- bind_rows(KID,KOD,LFLR) %>% group_by(measure,FIM) %>%
    summarize(meanFIM = mean(value),
              seFIM =qnorm(0.975) * (sd(value)/sqrt(length(value))))  %>%
    mutate(measure = factor(measure,levels = c("RNFactor",
                                               "VPFactor",
                                               "CapFactor",
                                               "RNVP",
                                               "VPCap",
                                               "RNCap"),
                            labels = c("RN","VP","Cap","RNVP","VPCap","RNCap"))) %>%
    mutate(Lim = limLabel) %>%
    mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                        labels = c("KID\n(base model + factor)",
                                   "KOD\n(full model - factor)",
                                   "2*LFLR\n(contrasting models with\nand without factor)")))

  totfiltered <- ggplot(FIMtotal, aes(x = 1,y=meanFIM,group=measure,fill=measure)) +
    scale_x_continuous(breaks=1,label="         overall") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
    geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
    geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
    geom_boxplot(data = FIMind,aes(x= 1,y=meanFIM,fill=measure),color="darkgrey",
                 outlier.size=1,position=position_dodge(.8),alpha=0.25) +


    geom_errorbar(data = FIMtotal,aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                  width=0,color="black",size=0.5)+
    geom_point(data=FIMtotal, aes(x = 1,y=meanFIM,group=measure,fill=measure),
               stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +
    # scale_fill_manual(values = c("#0072B2", "#D55E00","#009E73"),
    #                   labels = c("Response Noise","Capacity Limit", "RN + Cap Limit"))+
    scale_y_continuous(limits = c(-20,55),breaks=seq(-10,50,10))+
    #  geom_text(aes(label = starred))+
    facet_grid(FIM~Lim) +

    theme(axis.text.y = element_text(size=10),
          axis.ticks.y = element_line(),
          axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9,color="black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          plot.title = element_text(hjust=0.5,face="plain"),
          strip.text.x = element_text(size = 10,color="white"),
          strip.text.y = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12,face="plain"),
          legend.position = "none",
          legend.text.align = 0,
          legend.key.height = unit(0.5,"line"),
          legend.title = element_blank(),
          legend.box = "vertical")

  totfull <- ggplot(FIMtotal, aes(x = 1,y=meanFIM,group=measure,fill=measure)) +
    scale_x_continuous(breaks=1,label="         overall") +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 9.2,linetype="dashed",color="grey") +
    geom_hline(yintercept = 6.8,linetype="dashed",color="grey") +
    geom_hline(yintercept = 4.6,linetype="dashed",color="grey") +
    geom_boxplot(data = FIMind,aes(x= 1,y=meanFIM,fill=measure),color="darkgrey",
                 outlier.size=1,position=position_dodge(.8),alpha=0.25) +


    geom_errorbar(data = FIMtotal,aes(ymin = meanFIM-seFIM,ymax=meanFIM+seFIM),position=position_dodge(0.8),
                  width=0,color="black",size=0.5)+
    geom_point(data=FIMtotal, aes(x = 1,y=meanFIM,group=measure,fill=measure),
               stat="identity",position=position_dodge(0.8),shape=21,color="black",size=2) +
    scale_fill_manual(values = c("#0072B2","#D55E00", "#CC79A7",
                                 "#F0E442","#999999",
                                 "#009E73"),
                      labels = c("Response Noise","Variable Precision","Capacity Limit",
                                 "RN + VP", "VP + Cap Limit",
                                 "RN + Cap Limit"))+
    scale_y_continuous(limits = c(-20,200),breaks=seq(0,200,50))+
    #  geom_text(aes(label = starred))+
    facet_grid(FIM~Lim) +

    theme(axis.text.y = element_text(size=10),
          axis.ticks.y = element_line(),
          axis.text.x = element_text(size=10,angle=90,vjust=0.5,hjust=0.9,color="black"),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent"),
          plot.title = element_text(hjust=0.5,face="plain"),
          strip.text.x = element_text(size = 10,color="white"),
          strip.text.y = element_blank(),
          legend.key = element_rect(colour = "transparent", fill = "white"),
          legend.background = element_rect(fill = "transparent"),
          legend.text = element_text(size = 12,face="plain"),
          legend.position = "none",
          legend.text.align = 0,
          legend.key.height = unit(0.5,"line"),
          legend.title = element_blank(),
          legend.box = "vertical")


  plot_grid(totfull,full,rel_widths=c(0.3,1),nrow=1)

}
plot_grid(FIMVPEP2("_FM_"),FIMVPEP2("_P_"),FIMVPEP2("_U_"),nrow=3,rel_heights=c(1,1,1))

#ggsave("FIM_InteractionCapacityRN2.png", units="cm", width=35, height=45, dpi=600)
