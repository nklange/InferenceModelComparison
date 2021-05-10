# Make Figure 2


library("vwmvp")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("ggplot2")
library("stringr")
library("cowplot")
library(ungeviz)
# Load files -------------------------------------------------------------------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)

# Add picture exptrial / model comp--------------------------------------------------------

ExpTrial <- ggdraw() + draw_image("IllustrationVMVP/ExpTrial.png")


ModelComp <- ggdraw() + draw_image("IllustrationVMVP/ModelComparison.png")

# Descriptive data -------------------------------------------------------------

MakeSumStat <- function(data){
  data %>% group_by(set_size) %>%
    mutate(error_0 = circular::circular(error_0),units = "radians") %>%
    summarize(kurtosis = vwmvp::circkurtosis(error_0),
              circvar = circular::var.circular(error_0),
              meanabsdev = as.numeric(mean(abs(error_0)))) %>%
    pivot_longer(!c(set_size),names_to = "measure",values_to = "value") %>%
    mutate(measure = factor(measure,levels = c("circvar","meanabsdev","kurtosis"),
                            labels = c(expression(paste("circular ", sigma^{2})),
                                       expression(paste("mean abs. deviation (in radian)")),
                                       expression(paste("kurtosis")))))
}


Dat_subj <- NULL
for(subj in unique(experimentfile$id)) {


  sumstat <- MakeSumStat(experimentfile %>% filter(id == subj)) %>%
    mutate(id = subj,
           exp = unique(experimentfile %>% filter(id == subj) %>% .$exp))

  Dat_subj <- Dat_subj %>% bind_rows(sumstat)

}

Dat_exp <-Dat_subj %>% group_by(exp,measure,set_size) %>%
  summarize(meanmeasure = mean(value))




# circular variance by Exp
DataSumStat <- ggplot(Dat_subj, aes(x = factor(set_size), y = value, colour=exp,fill=exp)) +

  geom_point(data = Dat_subj,aes(x = set_size, y = value,colour = exp),
             shape=16,position = position_dodge(0.8), alpha = 0.4,size=2) +
  geom_line(data = Dat_exp,aes(x = set_size,y = meanmeasure,group = exp, colour = exp),
            size=1,
            position = position_dodge(0.8))+
  geom_point(data = Dat_exp,aes(x = factor(set_size), y = meanmeasure),
             shape=22,colour="black",position = position_dodge(0.8),size=3)+
  facet_wrap(~measure,ncol=2,labeller=label_parsed,scales="free_y")+
  scale_y_continuous(limits=c(0,NA))+
  scale_fill_discrete(guide= guide_legend(ncol=1))+
  scale_x_discrete(name = "Set size") +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12,hjust = 0.225),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.justification = "top",
        #legend.position = c(0.7,0.2),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.title = element_blank(),
        legend.text.align = 0)

DataSumStat

# Illustrate VP vs VM (fit separately to set sizes) -----------------------------
# VM : SA_RNplus: 1-parameter model where mean-precision kappa is estimated separately for all set sizes
# VP: VPnosetsize: 2-parameter model where mean-precision kappa and tau are estimated separately for all set sizes

PARTICIPANT <- "ol17_e1_6"

models <- c("SA_RNplus","VPnosetsize")
IllustrationPred <- bind_rows(readRDS(paste0("IllustrationVMVP/ol17_e1_6_VPnosetsize_prediction.rds")),
                              readRDS(paste0("IllustrationVMVP/ol17_e1_6_SA_RNplus_prediction.rds"))
) %>%
  mutate(model = factor(model,levels = c("SA_RNplus","VPnosetsize"))) %>%
  mutate(setsize = paste("Set size", setsize))

IllustrationAIC <-  bind_rows(readRDS(paste0("IllustrationVMVP/ol17_e1_6_VPnosetsize_Full20Random.rds")),
                              readRDS(paste0("IllustrationVMVP/ol17_e1_6_SA_RNplus_Full20Random.rds"))) %>%
  mutate(model = factor(model,levels = c("SA_RNplus","VPnosetsize")))


# Errordistribution
experimentfile <- experimentfile %>% rename(setsize = "set_size") %>%
  mutate(setsize = paste("Set size", setsize))


errordistribution <- ggplot(data = experimentfile %>%

                              filter(id == PARTICIPANT),
                            aes(x = error_0)) +
  geom_histogram(aes(y = ..density..),color="darkgrey",fill="#edeaea",bins=60) +
  facet_wrap(.~setsize,ncol=4) +
  coord_cartesian(xlim = c(-pi/2,pi/2)) +
  #geom_text(data = Annotation, aes(x = 0.75,y = max(Prediction$prediction) + 0.1,label = Predicted),size=3) +

  # scale_y_continuous(limits=c(0,max(IndPrediction$normprediction) + 0.2)) +
  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_line(data = IllustrationPred,
            aes(x = data,y=prediction,color=model,group=model,linetype=model), size = 1) +
  #
  scale_color_manual(values = c("blue","red"),
                     labels = c(expression(paste("Von Mises")),
                                expression(paste("Variable Precision")))) +
  scale_linetype_manual(values = c("dashed","solid"),
                        labels = c(expression(paste("Von Mises")),
                                   expression(paste("Variable Precision")))) +

  theme(axis.text.x = element_text(size=11),
        plot.title = element_text(size = 14,hjust=0.5),
        axis.text.y = element_text(size=18,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="transparent"),
        legend.position = c(0.13,0.35),
        #legend.box.background = element_rect(fill="white"),
        legend.text = element_text(size=12),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))

errordistribution

# Make circular variacne graph -------------------------------------------------
Prediction <- IllustrationPred

CV <- NULL
for(SS in unique(Prediction %>% .$setsize)){

  CircVars <- NULL

  for(modeln in unique(Prediction$model)){
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
    filter(id == PARTICIPANT) %>%
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
                                   levels = c("Data","SA_RNplus","VPnosetsize")))
CVdat_meandat <- CVdat %>% group_by(model) %>% summarize(meanE = mean(meanErr),
                                                         meanV = mean(Var),
                                                         meank = mean(kurt)) %>%
  pivot_longer(!c(model),names_to = "measure",values_to="value") %>%
  slice(rep(1:n(),2))


CVdat <- CV %>% filter(model == "Data") %>%  slice(rep(1:n(),each=2))
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
                                     expression(paste("kurtosis")))))


Vars1 <- ggplot(CVplot %>% filter(model == "Data"), aes(x = setsize, y = value)) +
  geom_hpline(aes(linetype = model),size=1.5,width=1) +
  geom_point(data = CVplot %>% filter(model != "Data"),
             aes(x = setsize, y = value,shape = model, fill = model),
             position = position_dodge(0.5),size=4,alpha=0.6,color="black") +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_grid(.~measure,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25))+
  scale_shape_manual(values = c(21,22),
                     labels = c(
                       expression(paste("Von Mises")),
                       expression(paste("Variable Precision")))) +
  scale_fill_manual(values = c("blue","red"),
                    labels = c(
                      expression(paste("Von Mises")),
                      expression(paste("Variable Precision")))) +

  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.1,0.75),
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
  mutate(model = factor(model, levels = rev(models)))


Vars2_SumStat <- ggplot(SumStat2 %>% filter(measure2 != "Fit"), aes(x = model, y = value)) +
  geom_bar(data = SumStat2 %>% filter(measure2 != "Fit"),stat="identity",
           aes(x = model, y = value, fill = model,linetype=model),
           color="black",width=1) +
  coord_flip()+
  scale_x_discrete(labels = c("VP","VM"

  ))+
  facet_wrap(.~measure,nrow=3,labeller=label_parsed) +
  scale_y_continuous(limits  = c(0,0.4),name="normalized RMSD")+
  scale_fill_manual(values = c("red","blue")) +
  scale_linetype_manual(values = c("solid","dashed"))+
  #ggtitle("(N)RMSD")+
  theme(axis.text.y = element_text(size=11),
        axis.text.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        plot.title = element_text(hjust=0.5,face="plain"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
Vars2_SumStat



#part1<-plot_grid(ExpTrial,pIllu,ncol=2,rel_widths = c(1,4),labels =c('A',NA))



# Add picture modelcomparison --------------------------------------------------------


p1 <- plot_grid(ExpTrial,errordistribution,ncol=2, rel_widths=c(1,4),labels=c("A","B"),scale=.975)
p2 <- plot_grid(Vars1,Vars2_SumStat,ncol=2,rel_widths=c(3,1),labels=c("C","D"))


p3 <- ggdraw(DataSumStat) +
draw_plot(ModelComp, .5, -.025, .5, .5) +
  draw_plot_label(
    c("E", "F"),
    c(0, 0.49),
    c(1, 0.455),
    size = 14
  )


plot_grid(p1,p2,p3,ncol=1,rel_heights=c(0.9,0.6,1.2),align = 'v', axis = 'l')

#ggsave("CaseStudyIntro.png", units="cm", width=35, height=40, dpi=600)














#mean absolute error by Exp
ggplot(dats, aes(x = factor(set_size), y = meanErr, colour=exp,fill=exp)) +

  geom_point(data = dats,aes(x = set_size, y = meanErr, colour = exp),
             shape=16,position = position_dodge(0.8), alpha = 0.4,size=2) +
  geom_line(data = datsgroup,aes(x = set_size,y = mErr,group = exp, colour = exp),
            size=1,
            position = position_dodge(0.8))+
  geom_point(data = datsgroup,aes(x = factor(set_size), y = mErr),
             shape=22,colour="black",position = position_dodge(0.8),size=3)+

  scale_y_continuous(limits = c(0,1.7),name="Mean absolute error (radian)")+

  scale_x_discrete(name = "Set size") +
  theme(axis.text = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        #legend.position = c(0.85,0.2),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.title = element_blank(),
        legend.text.align = 0)

# Split by other-preferences

dats <- experimentfile %>%
  group_by(exp,id,set_size) %>%
  summarize(meanErr = mean(abs(error_0)),
            circvar = circular::var.circular(error_0)) %>%
  mutate(category = ifelse(exp %in% c("pratte17"),"pratte17",
                           ifelse(exp %in% c("WM4"),"WM4",
                                  ifelse(exp %in% c("VSCGM21a"),"VSCGM21a",
                                         ifelse(exp %in% c("VSCGM21b"),"VSCGM21b","others"))))) %>%
  mutate(exp = factor(exp))

datsout <- dats %>% filter(category != "others")

datsgroup <- experimentfile %>%
  group_by(exp,id,set_size) %>%
  summarize(meanErr = mean(abs(error_0)),
            circvar = circular::var.circular(error_0)) %>%
  group_by(exp,set_size) %>%
  summarize(mErr = mean(meanErr),
            mcircvar = mean(circvar),
            ci = qnorm(.975) * sd(circvar)/length(circvar)) %>%
  mutate(category = ifelse(exp %in% c("pratte17"),"pratte17",
                           ifelse(exp %in% c("WM4"),"WM4",
                                  ifelse(exp %in% c("VSCGM21a"),"VSCGM21a",
                                         ifelse(exp %in% c("VSCGM21b"),"VSCGM21b","others"))))) %>%
  mutate(exp = factor(exp))

##mean absolute error by Exp
ggplot(dats, aes(x = set_size, y = meanErr)) +
  geom_point(data = dats %>% filter(category == "others"),
             aes(x = set_size, y = meanErr,group=id),
             colour = "#cccccc",shape=16,
             position = position_dodge(0.6),alpha = 0.8,size=1)+
  geom_line(data = datsgroup %>% filter(category == "others"),
            aes(y = mErr,group = exp),size=1.5,colour="#cccccc",alpha = 0.4,
            position = position_dodge(0.6))+
  geom_point(data = datsout,
             aes(x = set_size, y = meanErr,colour = factor(category),group=id),
             shape=16,position = position_dodge(0.6), alpha = 0.5,size=2) +
  geom_line(data = datsgroup %>% filter(category != "others"),
            aes(y = mErr,group = exp, colour = factor(category)),size=1.5,
            position = position_dodge(0.6))+
  geom_point(data = datsgroup %>% filter(category != "others"),
             aes(x = set_size, y = mErr,colour = factor(category),group=exp),
             size = 4,shape=16,position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(0,1.7),name="Mean absolute error (radian)")+
  scale_colour_manual(values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))  +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8), name = "Set size") +
  theme(
    axis.text = element_text(size=12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(0.85,0.2),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.title = element_blank(),
    legend.text.align = 0)

##circular variance by Exp
ggplot(dats, aes(x = set_size, y = circvar)) +

  geom_point(data = dats %>% filter(category == "others"),
             aes(x = set_size, y =circvar,group=id),
             colour = "#cccccc",shape=16,
             position = position_dodge(0.6),alpha = 0.8,size=1)+
  geom_line(data = datsgroup %>% filter(category == "others"),
            aes(y = mcircvar,group = exp),size=1.5,colour="#cccccc",alpha = 0.4,
            position = position_dodge(0.6))+
  geom_point(data = datsout,
             aes(x = set_size, y = circvar,colour = factor(category),group=id),
             shape=16,position = position_dodge(0.6), alpha = 0.5,size=2) +
  geom_line(data = datsgroup %>% filter(category != "others"),
            aes(y = mcircvar,group = exp, colour = factor(category)),size=1.5,
            position = position_dodge(0.6))+
  geom_point(data = datsgroup %>% filter(category != "others"),
             aes(x = set_size, y = mcircvar,colour = factor(category),group=exp),
             size = 4,shape=16,position = position_dodge(0.6)) +
  scale_y_continuous(limits = c(0,1),name="Circular variance (radian)")+
  scale_colour_manual(values=c("#009E73", "#0072B2", "#D55E00", "#CC79A7"))  +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8), name = "Set size") +
  theme(
    axis.text = element_text(size=12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = c(0.85,0.2),
    legend.background = element_rect(fill = "transparent"),
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.title = element_blank(),
    legend.text.align = 0)
