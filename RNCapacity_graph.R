
source("CompareModels_LimitRN.R") # for quantitative fit plots

library(vwmvp)
library("grid")
library("gridExtra")
library("cowplot")
library(ungeviz)


# Stability by points ----------------------------------------------------------
# taken 1:1 from DiagnoseModels.R

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,
                            data_pratte17)



FittedFull <- readRDS("Fits/FitFull.rds")
numberid <- length(unique(FittedFull$cvid))
exclall <- c()

models <- c("MK_RNminus", "MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")
#models <- c("MK_RNminus","J_RNminus","MK_RNplus","J_RNplus")
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
  theme(axis.text.x = element_text(size=11),#,angle=90,vjust=0.5
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

# Error Prediction -------------------------------------------------------------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,
                            data_pratte17)


experimentfile <- experimentfile %>%
  rename("setsize" = set_size)%>%
  mutate(setsize = paste0("Set size ", setsize))

experiment <- "ol17_e1"
AggPrediction <- readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds")) %>%
  filter(leftout == 0)

models <- c("MK_RNminus", "MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")


AggPrediction <- AggPrediction %>%
  filter(model %in% models) %>%
  #filter(setsize %in% c("Set size: 1","Set size: 2","Set size: 3","Set size: 8")) %>%
  mutate(model = factor(model, levels=models)) %>%
  mutate(setsize = paste0("Set size ", setsize))



aggerrordistribution <- ggplot(data = experimentfile %>% filter(exp==experiment),
                               aes(x = error_0)) +
  geom_histogram(aes(y = ..density..),color="darkgrey",fill="#edeaea",bins=60,alpha=0.8) +
  facet_wrap(.~setsize,ncol=4) +
  # geom_label(data = Annotation, aes(x = -1.2,
  #                                  y = max(AggPrediction %>% .$prediction) + 0.1,
  #                                  label = exp),size=4,fill="#edeaea") +
  #scale_y_continuous(limits=c(0,max(AggPrediction  %>% filter(exp == experiment) %>% .$normprediction) * 45 + 0.2)) +
  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(-pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  coord_cartesian(xlim = c(-pi/2,pi/2))+
  geom_line(data = AggPrediction,
            aes(x = data,y=prediction,linetype=factor(model),color=factor(model),
                group=factor(model)), size = 0.75) +
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
  guides(color=guide_legend(ncol=4,override.aes = list(size = 1.1) ),
         linetype = guide_legend(ncol=4)) +
  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 14,hjust=0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background.x = element_rect(fill = "transparent"),
        strip.background.y = element_rect(fill = "#edeaea",color="black"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="lightgrey"),
        legend.position = c(0.75,0.3),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))
aggerrordistribution
# Summary statistics -----------------------------------------------------------



load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size")

experiment <- "ol17_e1"
models <- c("MK_U_RNminus","MK_P_RNminus","MK_RNminus",
            "MK_FM_RNplus","MK_FM_RNminus",
            "MK_RNplus","MK_P_RNplus","MK_U_RNplus")
models2 <- c("MK_RNminus","MK_RNplus",
             "MK_FM_RNminus","MK_FM_RNplus"
             ,"MK_P_RNminus","MK_P_RNplus",
             "MK_U_RNminus","MK_U_RNplus")

models3 <- c("MK_U_RNplus","MK_U_RNminus",
             "MK_P_RNplus","MK_P_RNminus",
             "MK_FM_RNplus","MK_FM_RNminus",
             "MK_RNplus","MK_RNminus")
Prediction <- readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds")) %>%
  filter(leftout == 0) %>% filter(model %in% models)


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

CV <- CV %>% mutate(model = factor(model,
                                   levels = c("Data",models)))

CVdat <- CV %>% filter(model == "Data") %>%  slice(rep(1:n(),each=length(models)))
CVdat_meandat <- CVdat %>% group_by(model) %>% summarize(meanE = mean(meanErr),
                                                         meanV = mean(Var),
                                                         meank = mean(kurt)) %>%
  pivot_longer(!c(model),names_to = "measure",values_to="value") %>%
  slice(rep(1:n(),length(models)))


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
  mutate(model = factor(model, levels = c("Data",models2)))


Vars1 <- ggplot(CVplot %>% filter(model == "Data") %>% filter(setsize < 9), aes(x = setsize, y = value)) +
  geom_hpline(aes(linetype = model),size=1.5,width=1.2) +
  geom_point(data = CVplot %>% filter(model != "Data") %>% filter(setsize < 9),
             aes(x = setsize, y = value,shape = model, fill = model),
             position = position_dodge(1),size=3.5,alpha=0.7,color="black") +
  scale_x_continuous(breaks=c(1:8),labels=c(c(1:8)),name="Set size") +
  facet_grid(.~measure,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25))+
  scale_shape_manual(values = rep(c(21,22),4),
                     labels = c(
                       expression(paste("VP(",kappa,")A-")),
                       expression(paste("VP(",kappa,")A+")),
                       expression(paste("VP(",kappa,")F-")),
                       expression(paste("VP(",kappa,")F+")),
                       expression(paste("VP(",kappa,")P-")),
                       expression(paste("VP(",kappa,")P+")),
                       expression(paste("VP(",kappa,")U-")),
                       expression(paste("VP(",kappa,")U+"))

                     )
                     ) +

  scale_fill_manual(values =  c("#b2a8e0","#332288",
                                "#f2d188","#E69F00",
                                "#9ee6da","#44AA99",
                                "#f1dae7","#CC79A7" ),
  labels = c(
    expression(paste("VP(",kappa,")A-")),
    expression(paste("VP(",kappa,")A+")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("VP(",kappa,")U+"))

  )
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
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.1,0.7),
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
Vars1



SumStat2 <-CVdat2 %>%
  mutate(measure = factor(measure, levels = c("ssqvar","ssqe","ssqkurt"),
                          labels = c(expression(paste("circular ",sigma^2)),
                                     expression(paste("mean abs. deviation (radian)")),
                                     expression(paste("kurtosis"))))) %>%
  mutate(model = factor(model, levels = models3))


Vars2_SumStat <- ggplot(SumStat2 %>% filter(measure2 != "Fit"), aes(x = model, y = value)) +
  geom_bar(data = SumStat2 %>% filter(measure2 != "Fit"),stat="identity",
           aes(x = model, y = value, fill = model,linetype=model),
           color="black",width=1) +
  coord_flip()+
  scale_x_discrete(
    labels = c(
      expression(paste("VP(",kappa,")U+")),
      expression(paste("VP(",kappa,")U-")),
      expression(paste("VP(",kappa,")P+")),

      expression(paste("VP(",kappa,")P-")),
      expression(paste("VP(",kappa,")F+")),

      expression(paste("VP(",kappa,")F-")),
      expression(paste("VP(",kappa,")A+")),
      expression(paste("VP(",kappa,")A-"))




    )
  )+
  facet_wrap(.~measure,nrow=3,labeller=label_parsed) +
  scale_y_continuous(limits  = c(0,0.155),name="normalized RMSD",expand = c(0.002, 0.002),
                     breaks = c(0,0.05,0.10,0.15),
                     labels =  c("0","0.05","0.10","0.15"))+
  scale_fill_manual(values =  c("#CC79A7","#f1dae7",
                                "#44AA99","#9ee6da",
                                "#E69F00","#f2d188",
                                "#332288","#b2a8e0"



                                )) +
  scale_linetype_manual(values = rep(c("solid","dashed"),4))+
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
Vars2_SumStat


SumStat <- plot_grid(Vars1,Vars2_SumStat,labels=c("E","F"),nrow=1,rel_widths = c(3,1),align = 'h', axis = 't')

# CP: Finalize -----------------------------------------------------------------

plot_grid(datastabpoints,Quant,aggerrordistribution,SumStat,labels = c("A",NA,"D",NA),ncol=1,
          rel_heights = c(0.6,0.9,0.8,0.8))

ggsave("Test2.png", units="cm", width=35, height=45, dpi=600)


# Shen & Ma (2019) Factorial Importance graphs ---------------------------------


models <- c("MK_RNminus", "MK_RNplus",
            "MK_FM_RNminus", "MK_FM_RNplus",
            "MK_P_RNminus","MK_P_RNplus",
            "MK_U_RNminus","MK_U_RNplus"

            )

FIMCapacity <- function(BasisFit,capLimit){


  Capmodels <- c("MK_FM_RNminus", "MK_FM_RNplus",
              "MK_P_RNminus","MK_P_RNplus",
              "MK_U_RNminus","MK_U_RNplus"

  )



  models <- c("MK_RNminus","MK_RNplus",Capmodels[grepl(capLimit,Capmodels)])

  if(capLimit == "_FM_"){
    limLabel <- "Fixed capacity limit"
    LegendInfo <- c(0.25,0.6)
    StripLabel <- element_blank()
    Ticks <- element_line()
    YAxis <- element_text(size=10)
  } else if(capLimit == "_P_"){
    limLabel <- "Variable capacity limit (Poisson-dist.)"
    LegendInfo <- c(0.25,0.6)
    StripLabel <- element_blank()
    Ticks <- element_blank()
    YAxis <- element_blank()
  } else {
    limLabel <- "Variable capacity limit (Unif.-dist.)"
    LegendInfo <- c(0.25,0.6)
    StripLabel <- element_text(size=10)
    Ticks <- element_blank()
    YAxis <- element_blank()

  }

FittedFull <- BasisFit %>%
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
  mutate(type = case_when(model == "MK_RNminus" ~ 1,
                          model %in% c("MK_U_RNplus","MK_FM_RNplus","MK_P_RNplus") ~ 3,
                          TRUE ~ 2)) %>%
  #filter(model %in% CapUmodels) %>%
  group_by(exp,id) %>% arrange(exp,id,type,FactorRN,FactorVP,FactorCap)

KIDKOD <- FittedFull %>% group_by(exp,id) %>% mutate(KID = AIC[1] - AIC,
                                                     KOD = AIC - AIC[length(AIC)],
                                                     best = AIC - min(AIC))


LFLR <- KIDKOD %>% mutate(expon = exp(-best)) %>%
  summarize(RNFactor = log(mean(expon[3:4])/mean(expon[1:2])),
            CapFactor = log(mean(expon[c(2,4)])/mean(expon[c(1,3)]))) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value")%>%
  mutate(FIM = "LFLR")

KID <- KIDKOD %>% summarize(RNFactor = KID[3],

                            CapFactor = KID[2]) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
  mutate(FIM = "KID")


KOD <- KIDKOD %>% summarize(RNFactor = KOD[2],

                            CapFactor = KOD[3]) %>%
  pivot_longer(!c(id,exp),names_to = "measure",values_to="value") %>%
  mutate(FIM = "KOD")

FIM <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,measure,FIM) %>%
  summarize(meanFIM = mean(value),
            seFIM = qnorm(0.975) * (sd(value)/sqrt(length(value)))) %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "CapFactor"),
                          labels = c("Response Noise","Capacity Limit"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)")))

FIMind <- bind_rows(KID,KOD,LFLR) %>% group_by(exp,id,measure,FIM) %>%
  summarize(meanFIM = mean(value),
            seFIM = qnorm(0.975) * (sd(value)/sqrt(length(value)))) %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "CapFactor"),
                          labels = c("Response Noise","Capacity Limit"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)"))) %>% ungroup()


expgroup <- ggplot(data = FIMind, aes(x = exp,y=meanFIM,fill=measure)) +
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

  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
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
        strip.text.y = element_text(size=9),
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
            seFIM = qnorm(0.975) * (sd(value)/sqrt(length(value)))) %>%
  mutate(measure = factor(measure,levels = c("RNFactor",
                                             "CapFactor"),
                          labels = c("Response Noise","Capacity Limit"))) %>%
  mutate(Lim = limLabel) %>%
  mutate(FIM = factor(FIM,levels=c("KID","KOD","LFLR"),
                      labels = c("KID\n(base model + factor)",
                                 "KOD\n(full model - factor)",
                                 "2*LFLR\n(contrasting models with\nand without factor)")))

tot <- ggplot(FIMtotal, aes(x = 1,y=meanFIM,group=measure,fill=measure)) +
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
  scale_fill_manual(values = c("#0072B2", "#D55E00"))+
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

plot_grid(tot,expgroup,rel_widths=c(0.3,1),nrow=1)

}


FittedFull <- readRDS("Fits/FitFull.rds")

LOTraining <- bind_rows(readRDS("Fits/CVLOSzO_Trainingfits.rds") %>% filter(leftout==1) %>%
                          mutate(AIC = Dev + 2 * nparams),
                        FittedFull %>% filter(exp == "RTT12"))

FM <- plot_grid(FIMCapacity(FittedFull,"_FM_"),FIMCapacity(LOTraining,"_FM_"),nrow=1,labels=c("A","B"),
          rel_widths=c(1,1))
P <- plot_grid(FIMCapacity(FittedFull,"_P_"),FIMCapacity(LOTraining,"_P_"),nrow=1,labels=c("C","D"),
                rel_widths=c(1,1))
U <- plot_grid(FIMCapacity(FittedFull,"_U_"),FIMCapacity(LOTraining,"_U_"),nrow=1,labels=c("E","F"),
                rel_widths=c(1,1))

plot_grid(FM,P,U,nrow=3,rel_heights=c(1,1,1))


#ggsave("FIM_CapacityRN.png", units="cm", width=35, height=45, dpi=600)

