#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2024-07-11
# Facetted figure:
# - annual temporal and spatial variability
# - monthly spatial variability 
# - number of months of litter fall
# with composition as variable included
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)

#============================ Dataset ===============================

# df.all.wide.info <- read.csv

# NO! average across traps ####
df.all.wide.trap <- df.all.wide.info 

df.all.wide.trap$div<-as.factor(df.all.wide.trap$tree_species_richness)
df.all.wide.trap$blk<-as.factor(df.all.wide.trap$block)
df.all.wide.trap$myc<-as.factor(df.all.wide.trap$mycorrhizal_type)
df.all.wide.trap$sr<-df.all.wide.trap$tree_species_richness
df.all.wide.trap$sr_myc<-paste(df.all.wide.trap$sr,df.all.wide.trap$myc,sep="_")#interaction terms
df.all.wide.trap$myc <- recode_factor(df.all.wide.trap$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
df.all.wide.trap$block <- as.factor(df.all.wide.trap$block)


# long format ####
df.all.long <- df.all.wide.trap %>% 
  pivot_longer(cols=c(20:29),
               names_to="species",
               values_to="litterfall")


#### 1) Annual temporal variability  ####
df.year.temp.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1, composition) |>  # initially "composition" was not included here! but we included it when calculation annual, cumulative and monthly litterfall
  summarise(litterfall_mean = mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID, composition) |> 
  summarise(litterfall_CV = sd(litterfall_mean, na.rm=T)/mean(litterfall_mean, na.rm=T))

temp.var.yearly <- ggplot()+
  geom_jitter(data = df.year.temp.variability , 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(data = df.year.temp.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Annual temporal variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
temp.var.yearly


#### 2) Annual spatial variability (Coefficient of variation CV) ####
df.year.spat.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1, composition) |>   # initially "composition" was not included here! but we included it when calculation annual, cumulative and monthly litterfall
  summarise(litterfall_CV = sd(litterfall, na.rm=T)/mean(litterfall, na.rm=T))|>  # mean between traps and sd between traps
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID, composition) |> 
  summarise(litterfall_CV = mean(litterfall_CV, na.rm=T))

spa.var.yearly <- ggplot()+
  geom_jitter(data = df.year.spat.variability , 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(data = df.year.spat.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Annual spatial variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  #scale_y_log10() +
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")
spa.var.yearly


#### 3) Monthly spatial variability (CV coefficient of variation) ####
df.month.spat.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, month1, composition) |>  # initially "composition" was not included here! but we included it when calculation annual, cumulative and monthly litterfall
  summarise(litterfall_CV = sd(litterfall, na.rm=T)/mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]),labels=month.abb))

df.month.spat.variability$month2 <- dplyr::recode_factor(df.month.spat.variability$month1, 
                                             "March"="Mar", 
                                             "April"="Apr",
                                             "May"="May",
                                             "June"="Jun",
                                             "July"="Jul",
                                             "August"="Aug",
                                             "September"="Sep",
                                             "October"="Oct",
                                             "November"="Nov",
                                             "December"="Dec",
                                             "January"="Jan",
                                             "February"="Feb")

df.month.spat.variability$month1 = factor(df.month.spat.variability$month2, 
                              levels = c(month.abb[3:12], month.abb[1:2]))

M = map_df( .x = unique(df.month.spat.variability$month1),
            .f = ~ {
              mod = 
                lme(log(litterfall_CV) ~ log2(sr) * myc, 
                    random= ~1|composition,
                    data= df.month.spat.variability |>
                      filter(month1 == .x), 
                    na.action = na.omit) |> 
                anova(type="sequential") |> 
                data.frame() |> 
                mutate(month = .x)
              mod$exp = rownames(mod)
              rownames(mod) <- NULL
              mod |> 
                select(month, explanatory = exp, 
                       numDF, denDF, F.value, p.value)
            }) |> 
  mutate(sign = if_else(p.value < 0.001, '***', 
                        if_else(p.value < 0.01, "**", 
                                if_else(p.value < 0.05, '*',
                                        if_else(p.value<0.1, '.', 'n.s.')))))
M

M$sign[M$explanatory == 'log2(sr)' & M$month == 'March'] = 
  paste0("log2(sr) = ",M$sign[M$explanatory == 'log2(sr)' & M$month == 'March'])
M$sign[M$explanatory == 'myc' & M$month == 'March'] = 
  paste0("myc = ",M$sign[M$explanatory == 'myc' & M$month == 'March'])
M$sign[M$explanatory == 'log2(sr):myc' & M$month == 'March'] = 
  paste0("log2(sr):myc = ",M$sign[M$explanatory == 'log2(sr):myc' & M$month == 'March'])

spa.var.monthly <- ggplot()+
  geom_jitter(data = df.month.spat.variability, 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(data = df.month.spat.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Monthly spatial variability"), 
       x = "Tree species richness")+
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1,10,100),labels=c(0.0001,0.001,0.01,0.1,1,10,100))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr)'),
            aes(x=1, y=20,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=15,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr):myc'),
            aes(x=1, y=10,
                label = sign,
                hjust = 0),
            color = 'black') +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(d)")

spa.var.monthly


#### 4) Number of months of litterfall ####

# average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac_mean = mean(Ac, na.rm = TRUE),	
                   Ae_mean = mean(Ae, na.rm = TRUE),
                   Be_mean = mean(Be, na.rm = TRUE),	
                   Ca_mean = mean(Ca, na.rm = TRUE),
                   Fa_mean = mean(Fa, na.rm = TRUE),	
                   Fr_mean = mean(Fr, na.rm = TRUE),	
                   Pr_mean = mean(Pr, na.rm = TRUE),	
                   Qu_mean = mean(Qu, na.rm = TRUE), 
                   So_mean = mean(So, na.rm = TRUE),	
                   Ti_mean = mean(Ti, na.rm = TRUE), 
                   cont_mean = mean(cont, na.rm = TRUE))

df.all.wide.mean$div<-as.factor(df.all.wide.mean$tree_species_richness)
df.all.wide.mean$blk<-as.factor(df.all.wide.mean$block)
df.all.wide.mean$myc<-as.factor(df.all.wide.mean$mycorrhizal_type)
df.all.wide.mean$sr<-df.all.wide.mean$tree_species_richness
df.all.wide.mean$sr_myc<-paste(df.all.wide.mean$sr,df.all.wide.mean$myc,sep="_")#interaction terms
df.all.wide.mean$myc <- recode_factor(df.all.wide.mean$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")

# sum across species per plot
df.litter.sum = df.all.long |> 
  group_by(block, plotID, div, sr, myc, month1, composition) |>  # initially "composition" was not included here! but we included it when calculation annual, cumulative and monthly litterfall
  summarise(litterfall = sum(dryweight, na.rm = T),
            litterfall_mean = mean(dryweight, na.rm =T),
            litterfall_sd = sd(dryweight, na.rm = T))

# monthly coverage - sum per plot, if litter dryweight larger then 0 ####
df.litter.cover = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall >0) |>
  summarise(number_month_litterfall = n())

temp.cov <- ggplot(df.litter.cover, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(b)")

temp.cov


library(gridExtra)
plot2 <-grid.arrange(layout_matrix = rbind(c(1,2,3), 
                                           c(4,4,4)),
                     grobs= list(temp.var.yearly, temp.cov, spa.var.yearly, spa.var.monthly))

plot2





############# test different threshold for temporal coverage 
# monthly coverage - sum per plot, if litter dryweight larger then 0 ####
df.litter.cover.0 = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall >0) |>
  summarise(number_month_litterfall = n())

temp.cov.0 <- ggplot(df.litter.cover.0, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall (>0g leaf litter)"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")

temp.cov.0


# monthly coverage - sum per plot, if litter dryweight larger then 0.5 ####
df.litter.cover.0.5 = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall >0.5) |>
  summarise(number_month_litterfall = n())

temp.cov.0.5 <- ggplot(df.litter.cover.0.5, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall (>0.5g leaf litter)"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(b)")

temp.cov.0.5

# monthly coverage - sum per plot, if litter dryweight larger then 1 ####
df.litter.cover.1 = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall >1) |>
  summarise(number_month_litterfall = n())

temp.cov.1 <- ggplot(df.litter.cover.0.5, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall (>1g leaf litter)"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")

temp.cov.1


# monthly coverage - sum per plot, if litter dryweight larger then 1 ####
df.litter.cover.1.5 = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall >1.5) |>
  summarise(number_month_litterfall = n())

temp.cov.1.5 <- ggplot(df.litter.cover.1.5, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall (>1.5g leaf litter)"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 12),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(d)")

temp.cov.1.5


library(cowplot)

plot.cov <- plot_grid(temp.cov.0, temp.cov.0.5, temp.cov.1, nrow = 1, align = "h")
plot.cov


# 5) Model ####
library(nlme)

# model 0
mod.coverage.0 =
  lme(number_month_litterfall ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.litter.cover.0)

# Anova (Type III SS) ####
anova(object = mod.coverage.0, type = 'sequential')

# model 0.5
mod.coverage.0.5 =
  lme(number_month_litterfall ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.litter.cover.0.5)

# Anova (Type III SS) ####
anova(object = mod.coverage.0.5, type = 'sequential')


# model 1
mod.coverage.1 =
  lme(number_month_litterfall ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.litter.cover.1)

# Anova (Type III SS) ####
anova(object = mod.coverage.1, type = 'sequential')


### end ###