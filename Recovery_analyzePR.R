library("plyr")
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
library("tmvtnorm")
library("grid")
library("gridExtra")
library(cowplot)

options(pillar.sigfig = 5)
source("Recovery_analyzeMR.R")

# models -----------------------------------------------------------------------

models <- c("MK_RNminus", "MK_FM_RNplus","MK_RNplus","MK_U_RNplus","MK_P_RNplus",
            "MK_P_RNminus","MK_U_RNminus","MK_FM_RNminus","J_RNminus","J_RNplus")

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

# Load data --------------------------------------------------------------------

load_files <- function(path,pattern) {

  files <- dir(path, pattern = pattern,full.names = TRUE,recursive = TRUE)
  # search for test data (dataTEST) file in all subdirectories (recursive) of path
  tables <- lapply(files, readRDS)
  do.call(bind_rows, tables)
  # collate all of them into 1 dataframe
}



# Multivariate generating parameters -------------------------------------------

# for analysis of generating parameter estimates
# identical for 420 and 820 trials per set size
generating <- readRDS("RecoveryData/ParameterRecovery_data_120_rmv.rds")

# Analyse generating parameters in MR/PR

# group by data sets where K is below/above max set size
# estimate of K irrelevant for FM model if K > 8
generating %>%
  .$Parameters %>%
  bind_rows() %>%
  filter(model == "MK_U_RNminus") %>%

  # for MK_U_RNplus,MK_U_RNminus to report K not 2*K
  # comment out for other models
  mutate(K = K/2) %>%

  mutate(Kbelow = ifelse(K < 8,1,2)) %>%
  gather(key = "Parameter",
         value = "value","kappa", "tau","alpha", "kappa_r","K")  %>%


  group_by(Kbelow,Parameter) %>%
  summarize(means = mean(value),
            sds = sd(value),
            mins = min(value),
            max = max(value),
            n = length(Kbelow))


# Parameter Recovery

PRrecovery <- load_files("ParameterRecoveryrmv/","trialsSs")

recovery <- PRrecovery %>%
  group_by(trials,model,id) %>%
  arrange(objective) %>%
  slice(1)

recovery <- recovery %>%
  select(alpha,J1bar,tau,mkappa1,kappa_r,K,model,id,trials) %>%
  mutate(phi = ifelse(model %in% c("J_RNminus","J_RNplus"),J1bar,mkappa1)) %>%
  select(-J1bar,-mkappa1) %>%
  gather(key = "Parameters", value = "value", -model,-id,-trials) %>%
  filter(!(value == Inf | is.na(value))) %>%
  arrange(trials,model,Parameters)


generating <- generating %>%
  .$Parameters %>%
  bind_rows() %>%
  rename(phi = kappa) %>%
  gather(key = "Parameters",value = "value",-model,-parid) %>%
  filter(!(value == Inf | is.na(value))) %>%
  arrange(model,parid,Parameters) #%>% slice(rep(1:n(), each = 40))

# combine generating parameters and estimates by data set
recoverycom <- NULL
recovery$genvalue <- NA
for (subj in unique(generating %>% .$parid)){

  for (par in c("phi","alpha","tau","kappa_r","K")){


    singlepar <- generating %>%
      filter(parid == subj) %>%
      filter(Parameters == par) %>%
      .$value
    part <- recovery %>%
      filter(id == subj) %>%
      filter(Parameters == par) %>%
      mutate(genvalue = singlepar)
    recoverycom <- bind_rows(recoverycom,part)
  }
}

# Ensure that irrelevant K estimates do not penalize MKFM models
recovery <- recoverycom %>%
  mutate(value = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                           Parameters == "K" & value > 8),8,value),
         genvalue = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                              Parameters == "K" & genvalue > 8),8,genvalue))

# Calculate percent difference of generating and recovered value
recovery$percentvalue <- (recovery$value - recovery$genvalue)/recovery$genvalue

# For order in graphs
recovery$Parameters <- factor(recovery$Parameters,
                              levels=c("phi","alpha","tau","kappa_r","K"),
                              labels=c("phi","alpha","tau","kappa[r]","K"))
recovery$model <- factor(recovery$model,
                         levels = modellevels)

# percentag of data points not plotted because they exceed (abs)200% difference
AboveCrit <- recovery %>% mutate(above = ifelse(percentvalue > 2, 1, 0)) %>%
  group_by(trials,model,Parameters) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(model,Parameters) %>% filter(!notplotted %in% c(0,100) )

PRRecRMV <- ggplot(recovery, aes(x = model,y = percentvalue)) +

  geom_abline(slope=0,intercept = 0) +
  geom_rect(aes(xmin = 0.5,xmax = 10.5,ymin = -.2,ymax=.2),
            fill = "transparent",colour="grey",
            linetype="dotted")+
  geom_jitter(size = 1,shape = 16,alpha = 0.2, fill = "red",colour="red",width=0.01) +

  #geom_abline(intercept = 0.2,slope=0)+
  #geom_abline(intercept = -0.2,slope = 0)+
  geom_text(data = AboveCrit,
            mapping = aes(x = model, y = -1.5, label = paste0(round(notplotted),"%")),size = 2.5) +
  scale_y_continuous(name = "% difference",limits = c(-2,2),
                     breaks = c(-1,-.2,0,.2,1),labels = c("-100%","","0%","","100%")) +
  scale_x_discrete(name = "Model",labels = c(
    expression(paste("VP(",kappa,")A-")),
    expression(paste("VP(",italic(J),")A-")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("VP(",kappa,")A+")),
    expression(paste("VP(",italic(J),")A+")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+")),
    "total"
  )) +
  facet_grid(trials~Parameters,labeller = label_parsed) +
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

# Stability of fits (as with observed data)

diffObj <- PRrecovery %>%
  group_by(trials,model,id) %>%
  arrange(objective) %>%
  slice(1:2) %>%
  summarize(mindiffLL = (diff(objective) * 2),
            bestLL = objective[1]) %>%
  mutate(diffObj = mindiffLL/bestLL * 100)

Total <- bind_rows(diffObj) %>%
  mutate(model = factor(model,levels = modellevels,labels = modellabels))


AboveCrit <- Total %>%
  mutate(above = ifelse(diffObj > 0.5, 1, 0)) %>%
  group_by(trials,model) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(trials,model) %>%
  filter(!notplotted %in% c(0,100) )



pr_stabpercent <- ggplot(Total %>% filter(diffObj < 0.5),
                         aes(x = model,y = diffObj)) +

  geom_jitter(aes(colour = model),size = 2,alpha=0.4,width=0.05) +
  geom_abline(aes(intercept=0.025,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0.05,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_text(data = AboveCrit,
            mapping = aes(x = model, y = 0.52,
                          label = paste(round(notplotted),"%")),size = 2.5) +
  scale_x_discrete(name = "Model", labels = graphmodellabels) +
  scale_y_continuous(name = expression(paste(Delta,"Deviance (in %)")),
                     limits = c(0,0.53),
                     breaks = c(c(0,0.025,0.05),seq(0.1,0.5,0.1)),
                     labels = c("0%","0.025%","0.05%","0.1%",
                                "0.2%","0.3%","0.4%","0.5%")) +
  scale_colour_discrete(name = "Data") +
  facet_wrap(.~factor(trials),ncol=5) +
  theme(axis.text.x = element_text(size=12,angle=90,vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title= element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent",colour="transparent"))
pr_stabpercent

AboveCritpt <- Total %>% mutate(above = ifelse(mindiffLL > 10, 1, 0)) %>%
  group_by(trials,model) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(trials,model) %>%
  filter(!notplotted %in% c(0,100) )

pr_stabpoints_RMV <- ggplot(Total %>% filter(mindiffLL < 10),
                        aes(x = model,y = mindiffLL)) +

    geom_jitter(aes(color = model),size = 2,alpha=0.4,width=0.05) +
  geom_abline(aes(intercept=0.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=2,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=10,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_text(data = AboveCritpt,
            mapping = aes(x = model, y = 10.3,
                          label = paste0(round(notplotted),"%")),size = 2.5) +
  scale_x_discrete(name = "Model", labels = graphmodellabels) +
  scale_y_continuous(name = expression(paste(Delta,"Deviance (in points)")),
                     limits = c(0,10.5),breaks = c(c(0,0.5,1),seq(2,10,2))) +
  facet_wrap(.~factor(trials),labeller=label_parsed,ncol = 5) +
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
pr_stabpoints_RMV



# Median generating parameters -------------------------------------------

generating <- readRDS("RecoveryData/ParameterRecovery_data_120_median.rds")

# Analyse generating parameters in MR/PR

# group by data sets where K is below/above max set size
# estimate of K irrelevant for FM model if K > 8
generating %>%
  .$Parameters %>%
  bind_rows() %>%
  filter(model == "MK_FM_RNminus") %>%
  mutate(Kbelow = ifelse(K < 8,1,2)) %>%
  gather(key = "Parameter",
         value = "value","kappa", "tau","alpha", "kappa_r","K")  %>%
  group_by(Kbelow,Parameter) %>%
  summarize(means = mean(value),
            sds = sd(value),
            mins = min(value),
            max = max(value),
            n = length(Kbelow))


# Parameter Recovery

PRrecovery <- load_files("ParameterRecovery/","trialsSs")

recovery <- PRrecovery %>%
  group_by(trials,model,id) %>%
  arrange(objective) %>%
  slice(1)

recovery <- recovery %>%
  select(alpha,J1bar,tau,mkappa1,kappa_r,K,model,id,trials) %>%
  mutate(phi = ifelse(model %in% c("J_RNminus","J_RNplus"),J1bar,mkappa1)) %>%
  select(-J1bar,-mkappa1) %>%
  gather(key = "Parameters", value = "value", -model,-id,-trials) %>%
  filter(!(value == Inf | is.na(value))) %>%
  arrange(trials,model,Parameters)


generating <- generating %>%
  .$Parameters %>%
  bind_rows() %>%
  rename(phi = kappa) %>%
  gather(key = "Parameters",value = "value",-model,-parid) %>%
  filter(!(value == Inf | is.na(value))) %>%
  arrange(model,parid,Parameters) %>%
  slice(rep(1:n(), length(unique(recovery$trials))))


recovery$genvalue <- generating$value

# Ensure that irrelevant K estimates do not penalize MKFM models
recovery <- recovery %>%
  mutate(value = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                           Parameters == "K" & value > 8),8,value),
         genvalue = ifelse(((model == "MK_FM_RNminus" | model == "MK_FM_RNplus") &
                              Parameters == "K" & genvalue > 8),8,genvalue))

# Calculate percent difference of generating and recovered value
recovery$percentvalue <- (recovery$value - recovery$genvalue)/recovery$genvalue

# For order in graphs
recovery$Parameters <- factor(recovery$Parameters,
                              levels=c("phi","alpha","tau","kappa_r","K"),
                              labels=c("phi","alpha","tau","kappa[r]","K"))
recovery$model <- factor(recovery$model,
                         levels = modellevels)

# percentag of data points not plotted because they exceed (abs)200% difference
AboveCrit <- recovery %>% mutate(above = ifelse(percentvalue > 2, 1, 0)) %>%
  group_by(trials,model,Parameters) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(model,Parameters) %>% filter(!notplotted %in% c(0,100) )

PRmedian <- ggplot(recovery, aes(x = model,y = percentvalue)) +

  geom_abline(slope=0,intercept = 0) +
  geom_rect(aes(xmin = 0.5,xmax = 10.5,ymin = -.2,ymax=.2),
            fill = "transparent",colour="grey",
            linetype="dotted")+
  geom_jitter(size = 1,shape = 16,alpha = 0.2, fill = "red",colour="red",width=0.01) +

  #geom_abline(intercept = 0.2,slope=0)+
  #geom_abline(intercept = -0.2,slope = 0)+
  geom_text(data = AboveCrit,
            mapping = aes(x = model, y = -1.5, label = paste0(round(notplotted),"%")),size = 2.5) +
  scale_y_continuous(name = "% difference",limits = c(-2,2),
                     breaks = c(-1,-.2,0,.2,1),labels = c("-100%","","0%","","100%")) +
  scale_x_discrete(name = "Model",labels = c(
    expression(paste("VP(",kappa,")A-")),
    expression(paste("VP(",italic(J),")A-")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("VP(",kappa,")A+")),
    expression(paste("VP(",italic(J),")A+")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+")),
    "total"
  )) +
  facet_grid(trials~Parameters,labeller = label_parsed) +
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

# Stability of fits (as with observed data)

diffObj <- PRrecovery %>%
  group_by(trials,model,id) %>%
  arrange(objective) %>%
  slice(1:2) %>%
  summarize(mindiffLL = (diff(objective) * 2),
            bestLL = objective[1]) %>%
  mutate(diffObj = mindiffLL/bestLL * 100)

Total <- bind_rows(diffObj) %>%
  mutate(model = factor(model,levels = modellevels,labels = modellabels))


AboveCrit <- Total %>%
  mutate(above = ifelse(diffObj > 0.5, 1, 0)) %>%
  group_by(trials,model) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(trials,model) %>%
  filter(!notplotted %in% c(0,100) )



pr_stabpercent <- ggplot(Total %>% filter(diffObj < 0.5),
                         aes(x = model,y = diffObj)) +

  geom_jitter(aes(colour = model),size = 2,alpha=0.4,width=0.05) +
  geom_abline(aes(intercept=0.025,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0.05,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_text(data = AboveCrit,
            mapping = aes(x = model, y = 0.52,
                          label = paste0(round(notplotted),"%")),size = 2.5) +
  scale_x_discrete(name = "Model", labels = graphmodellabels) +
  scale_y_continuous(name = expression(paste(Delta,"Deviance (in %)")),
                     limits = c(0,0.53),
                     breaks = c(c(0,0.025,0.05),seq(0.1,0.5,0.1)),
                     labels = c("0%","0.025%","0.05%","0.1%",
                                "0.2%","0.3%","0.4%","0.5%")) +
  scale_colour_discrete(name = "Data") +
  facet_wrap(.~factor(trials),ncol=5) +
  theme(axis.text.x = element_text(size=12,angle=90,vjust=0.5),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        legend.title= element_blank(),
        axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.position = "none",
        legend.key = element_rect(fill = "transparent",colour="transparent"))
pr_stabpercent

AboveCritpt <- Total %>% mutate(above = ifelse(mindiffLL > 10, 1, 0)) %>%
  group_by(trials,model) %>%
  summarize(notplotted = sum(above)/length(above) * 100) %>%
  arrange(trials,model) %>%
  filter(!notplotted %in% c(0,100) )

pr_stabpoints_median <- ggplot(Total %>% filter(mindiffLL < 10),
                        aes(x = model,y = mindiffLL)) +

  geom_jitter(aes(color = model),size = 2,alpha=0.4,width=0.05) +
  geom_abline(aes(intercept=0.5,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=1,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=2,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=10,slope = 0),linetype="dotted",colour="grey") +
  geom_abline(aes(intercept=0,slope = 0),linetype="dotted",colour="grey") +
  geom_text(data = AboveCritpt,
            mapping = aes(x = model, y = 10.3,
                          label = paste0(round(notplotted),"%")),size = 2.5) +
  scale_x_discrete(name = "Model", labels = graphmodellabels) +
  scale_y_continuous(name = expression(paste(Delta,"Deviance (in points)")),
                     limits = c(0,10.5),breaks = c(c(0,0.5,1),seq(2,10,2))) +
  facet_wrap(.~factor(trials),labeller=label_parsed,ncol = 5) +
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



plot_grid(modeldevmedian,modeldevRMV,
          modelprop_strict_median,modelprop_strict_rmv,
          modelprop_lenient_median,modelprop_lenient_rmv,
                      nrow=3,labels=c("A","D","B","E","C","F"))

#ggsave("ModelRecovery.png", units="cm", width=30, height=40, dpi=600)


first <- plot_grid(pr_stabpoints_median,pr_stabpoints_RMV,nrow=1,rel_widths=c(5,3),labels=c("A","B"))
plot_grid(first,PRmedian,PRRecRMV,nrow=3,rel_heights=c(1,1.5,1.3),labels=c(NA,"C","D"))
#ggsave("ParameterRecovery.png", units="cm", width=35, height=40, dpi=600)
