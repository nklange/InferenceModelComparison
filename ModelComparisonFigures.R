library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("grid")
library("gridExtra")
library(cowplot)
linesize = 0.75


global1 <- tibble(dimension = rep(c(1:5),each = 6),
                 assumptions = rep(c(1:6),5),
                 chosen = c(
                   c(0,0,0,0,1,0),
                   c(0,0,0,2,0,0),
                   c(0,0,0,0,0,1),
                   c(0,0,0,0,0,1),
                   c(0,1,0,0,0,0))) %>% mutate(names = "Model 1")
global2 <- tibble(dimension = rep(c(1:5),each = 6),
                  assumptions = rep(c(1:6),5),
                  chosen = c(c(0,0,0,1,0,0),
                             c(0,0,0,2,0,0),
                             c(0,0,0,1,0,0),
                             c(0,1,0,0,0,0),
                             c(1,0,0,0,0,0))) %>% mutate(names = "Model 2")

globalgraph <- bind_rows(global1,global2)
globalgraph <- global %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,2,chosen))

globalg <- ggplot(data = globalgraph,aes(y = factor(dimension),x = factor(assumptions))) +
  geom_tile(data = globalgraph, aes(fill = factor(chosen)),color = "black",size=linesize) +
  scale_x_discrete(name = "Assumptions") +
  scale_y_discrete(name = "Dimension", labels = c(
    expression(paste(italic("E"))),
    expression(paste(italic("D"))),
    expression(paste(italic("C"))),
    expression(paste(italic("B"))),
    expression(paste(italic("A")))
  )) +
  facet_wrap(.~names)+
  scale_fill_manual(values = c("white","red","black")) +
 # facet_wrap(.~method) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour="transparent"),
        strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill = "transparent",
                                        colour="transparent"),
        strip.text = element_text(size = 12))


local1 <- global1 %>%  mutate(chosen = ifelse(dimension == 2 & assumptions == 4,1,chosen)) 

local2 <- local1 %>%
  mutate(chosen = ifelse(dimension == 5 & assumptions == 1,2,chosen))  %>% 
  mutate(chosen = ifelse(dimension == 5 & assumptions == 2,0,chosen)) %>%  
  mutate(names = "Model 3")
local3 <- local1 %>% mutate(chosen = ifelse(dimension == 5 & assumptions == 3,2,chosen)) %>% 
  mutate(chosen = ifelse(dimension == 5 & assumptions == 2,0,chosen)) %>% 
  mutate(names = "Model 4")
local1 <- local1 %>% mutate(chosen = ifelse(dimension == 5 & assumptions == 2,2,chosen))
  

localgraph <- bind_rows(local1,local2,local3)
localg <- ggplot(data = localgraph,aes(y = factor(dimension),x = factor(assumptions))) +
  geom_tile(data = localgraph, aes(fill = factor(chosen)),color = "black",size=linesize) +
  scale_x_discrete(name = "Assumptions") +
  scale_y_discrete(name = "Dimension", labels = c(
    expression(paste(italic("E"))),
    expression(paste(italic("D"))),
    expression(paste(italic("C"))),
    expression(paste(italic("B"))),
    expression(paste(italic("A")))
  )) +
  facet_wrap(.~names)+
  scale_fill_manual(values = c("white","black","red")) +
  # facet_wrap(.~method) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "transparent",colour="transparent"),
        strip.placement = "outside",
        legend.position = "none",
        strip.background = element_rect(fill = "transparent",
                                        colour="transparent"),
        strip.text = element_text(size = 12))

factorial1 <- local1 %>%  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,3,chosen)) %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 
factorial2 <- local2 %>%  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,3,chosen)) %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 
factorial3 <- local3 %>%  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,3,chosen)) %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 


factorial4 <- factorial1  %>% mutate(chosen = ifelse(dimension == 4 & assumptions == 5,3,chosen))  %>% 
  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,0,chosen)) %>%  mutate(names = "Model 5") %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 
factorial5 <- factorial2  %>% mutate(chosen = ifelse(dimension == 4 & assumptions == 5,3,chosen))  %>% 
  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,0,chosen)) %>%  mutate(names = "Model 6") %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 
factorial6 <- factorial3  %>% mutate(chosen = ifelse(dimension == 4 & assumptions == 5,3,chosen))  %>% 
  mutate(chosen = ifelse(dimension == 4 & assumptions == 6,0,chosen)) %>%  mutate(names = "Model 7") %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,4,chosen)) 




factorial7 <- factorial1 %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 8") %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen))
factorial8 <- factorial2  %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 9")  %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen))
factorial9 <- factorial3  %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 10")  %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen))
factorial10 <- factorial4  %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 11")  %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen)) 
factorial11 <- factorial5  %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 12")  %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen))
factorial12 <- factorial6  %>% mutate(chosen = ifelse(dimension == 2 & assumptions == 4,0,chosen))  %>%  mutate(names = "Model 13")  %>% 
  mutate(chosen = ifelse(dimension == 2 & assumptions == 3,4,chosen))



factorialgraph <- bind_rows(factorial1,factorial2,factorial3,
                            factorial4,factorial5,factorial6,
                            factorial7,factorial8,factorial9,
                            factorial10,factorial11,factorial12) %>% 
  mutate(names = factor(names,levels = c("Model 1","Model 3","Model 4", "Model 8","Model 9","Model 10",
                                         "Model 5","Model 6","Model 7","Model 11",
                                        "Model 12","Model 13")))

  factorialg <- ggplot(data = factorialgraph,aes(y = factor(dimension),x = factor(assumptions))) +
    geom_tile(data = factorialgraph, aes(fill = factor(chosen)),color = "black",size=linesize) +
    scale_x_discrete(name = "Assumptions") +
    scale_y_discrete(name = "Dimension", labels = c(
      expression(paste(italic("E"))),
      expression(paste(italic("D"))),
      expression(paste(italic("C"))),
      expression(paste(italic("B"))),
      expression(paste(italic("A")))
    )) +
    facet_wrap(.~names,ncol=6)+
    scale_fill_manual(values = c("white","black","red","blue","green")) +
    # facet_wrap(.~method) +
    theme(axis.text.x = element_text(size =12),
          axis.text.y = element_text(size=12),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = "transparent",colour="transparent"),
          strip.placement = "outside",
          legend.position = "none",
          strip.background = element_rect(fill = "transparent",
                                          colour="transparent"),
          strip.text = element_text(size = 12))

p1 <- plot_grid(globalg,localg,ncol=2,labels='AUTO',rel_widths=c(2.2,3),scale=.95)
plot_grid(p1,factorialg,ncol=1,labels=c(NA,"C"),rel_heights = c(1.3,2),align='l')

ggsave("Test1.png", units="cm", width=22, height=15, dpi=600)
