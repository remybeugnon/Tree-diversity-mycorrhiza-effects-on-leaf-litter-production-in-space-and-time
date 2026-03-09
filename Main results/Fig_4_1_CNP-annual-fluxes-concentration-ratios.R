#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-06-26
# annual effects (tree diversity and mycorrhizal type) on leaf litter C, N, P concentration, flux and CN, CP, NP
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)

#============================ Packages ===============================

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(ggpubr)

#============================ Datasets ===============================

# df.CNP.info <- read_csv
# df.mass.info <- read_csv

#============================ Data sorting ============================

# g/m2 in kg/ha

df.mass.info.kgha <- df.mass.info %>%
  dplyr::mutate(Ac_kgha = Ac_m2 * 10,	
                Ae_kgha = Ae_m2 * 10,
                Be_kgha = Be_m2 * 10,
                Ca_kgha = Ca_m2 * 10,
                Fa_kgha = Fa_m2 * 10,
                Fr_kgha = Fr_m2 * 10,
                Pr_kgha = Pr_m2 * 10,	
                Qu_kgha = Qu_m2 * 10, 
                So_kgha = So_m2 * 10,
                Ti_kgha = Ti_m2 * 10,
                cont_kgha = cont_m2 * 10)

df.mass.info.kgha <- df.mass.info.kgha %>%
  dplyr::select(cols=-c(8:19))


# 1) average across traps ####
df.mass.info.mean <- df.mass.info.kgha %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac = mean(Ac_kgha, na.rm = TRUE),	
                   Ae = mean(Ae_kgha, na.rm = TRUE),
                   Be = mean(Be_kgha, na.rm = TRUE),	
                   Ca = mean(Ca_kgha, na.rm = TRUE),
                   Fa = mean(Fa_kgha, na.rm = TRUE),	
                   Fr = mean(Fr_kgha, na.rm = TRUE),	
                   Pr = mean(Pr_kgha, na.rm = TRUE),	
                   Qu = mean(Qu_kgha, na.rm = TRUE), 
                   So = mean(So_kgha, na.rm = TRUE),	
                   Ti = mean(Ti_kgha, na.rm = TRUE), 
                   cont = mean(cont_kgha, na.rm = TRUE))%>%
  dplyr::select(col=-c(cont))


# long format ####
df.mass.long <- df.mass.info.mean %>% 
  pivot_longer(cols=c("Ac":"Ti"),
               names_to="species",
               values_to="mass")

df.mass.long <- na.omit(df.mass.long)

df.mass.long <- df.mass.long %>% 
  dplyr::mutate(month2 = recode(month1,
                                January = "Jan",
                                February = "Feb",
                                March = "Mar",
                                April = "Apr",
                                May = "May",
                                June = "Jun",
                                July = "Jul",
                                August = "Aug",
                                September = "Sep",
                                October = "Oct",
                                November = "Nov",
                                December = "Dec")) %>%
  ungroup() %>%
  dplyr::select(-c("month","month1","myc","sr","div","blk"))

df.CNP.info <- df.CNP.info |>
  dplyr::select(-c("block","plotName","tree_species_richness","mycorrhizal_type")) |>
  dplyr::rename(month2=month)

# get kg/kg or percentage out of the mg/g concentrations

df.CNP.info.percent <- df.CNP.info %>%
  dplyr::mutate(Cp = C * 0.001,	
                Np = N * 0.001,
                Pp = P * 0.001)

# join both datasets
df.mass.CNP <- full_join(df.mass.long, df.CNP.info.percent, by=c("plotID","month2","species","composition"))


df.mass.CNP$div<-as.factor(df.mass.CNP$tree_species_richness)
df.mass.CNP$blk<-as.factor(df.mass.CNP$block)
df.mass.CNP$myc<-as.factor(df.mass.CNP$mycorrhizal_type)
df.mass.CNP$sr<-df.mass.CNP$tree_species_richness
df.mass.CNP$myc <- recode_factor(df.mass.CNP$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")



# sum(concentration_species * biomass_species)/total biomass

# we need to consider the proportion of litter per species on each plot first
# otherwise we don't know the true backflow of nutrients to the ground

# Cflux and Nflux and Pflux per species
df.mass.CNP$Cflux <- df.mass.CNP$Cp * df.mass.CNP$mass
df.mass.CNP$Nflux <- df.mass.CNP$Np * df.mass.CNP$mass
df.mass.CNP$Pflux <- df.mass.CNP$Pp * df.mass.CNP$mass



### (1) annual effects - sum per plot (no month) ####

# exclude rows with NAs/missing data
df.mass.CNP.nona <- df.mass.CNP %>%
  filter(!(is.na(C) & is.na(N) & is.na(P)))

df.flux.annual <-
  #df.mass.CNP.nona %>%
  df.mass.CNP %>%
  group_by(plotID, blk, sr, myc, composition) |>
  summarise(Cflux.tot = sum(Cflux, na.rm = T),
            Nflux.tot = sum(Nflux, na.rm = T),
            Pflux.tot = sum(Pflux, na.rm = T),
            mass.tot = sum(mass, na.rm = T)) %>%
  mutate(Cc=(Cflux.tot/mass.tot)*1000)%>%   # *1000 so the unit is mg/g
  mutate(Nc=(Nflux.tot/mass.tot)*1000)%>%
  mutate(Pc=(Pflux.tot/mass.tot)*1000)%>%
  mutate(CN=Cc/Nc)%>%
  mutate(CP=Cc/Pc)%>%
  mutate(NP=Nc/Pc)%>%
  mutate(across(c(Cc, Nc, Pc), ~na_if(., 0))) #turn zeros in C, N, and P concentration to NA


# F0 annual biomass
F0 <- ggplot(df.flux.annual, 
             aes(x=sr, y=mass.tot, 
                 color = myc, fill = myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Annual leaf litter biomass"~~(g~m^-2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))
F0

# F1 annual Cflux
F1 <- ggplot(df.flux.annual, 
             aes(x=sr, y=Cflux.tot, 
                 color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=expression(atop("C flux", (kg ~ ha^-1 ~ yr^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(d)")
F1

# F2 annual C concentration
F2 <- ggplot(df.flux.annual, 
             aes(x=sr, y=Cc, 
                 color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= expression(atop("C concentration", (mg ~ g^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
F2

# F3 annual Nflux
F3 <-ggplot(df.flux.annual, 
            aes(x=sr, y=Nflux.tot, 
                color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=expression(atop("N flux", (kg ~ ha^-1 ~ yr^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y=element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(e)")
F3

# F4 annual N concentration
F4 <-ggplot(df.flux.annual, 
            aes(x=sr, y=Nc, 
                color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= expression(atop("N concentration", (mg ~ g^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(b)")
F4 

# F5 annual Pflux
F5 <-ggplot(df.flux.annual, 
            aes(x=sr, y=Pflux.tot, 
                color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=expression(atop("P flux", (kg ~ ha^-1 ~ yr^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(f)")
F5

# F6 annual P concentration
F6 <-ggplot(df.flux.annual, 
            aes(x=sr, y=Pc, 
                color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= expression(atop("P concentration", (mg ~ g^-1))),
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type"
                    ,
                    guide="none"
                    )+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type"
                     ,
                     guide="none"
                     )+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "bottom",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")
F6

# F7 annual CN 
F7 <-ggplot(df.flux.annual, 
            aes(x=sr, y=CN, 
                color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= "C : N", 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(g)")
F7

# F8 annual CP
F8 <-ggplot(df.flux.annual, 
            aes(x=sr, y=CP, 
                color = myc, fill = myc))+
  #geom_point(shape =21, size = 1, alpha=0.5)+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= "C : P", 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(h)")

F8

# F9 annual NP 
F9 <-ggplot(df.flux.annual, 
            aes(x=sr, y=NP, 
                color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
 # geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y= "N : P", 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none"
                    )+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none"
                     )+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(i)")
F9


legend <- get_legend(F2)

annual.plot <- ggarrange(F2, F4, F6, F1, F3, F5, F7, F8, F9,
          ncol = 3, nrow = 3, align ="hv", common.legend = T, legend="bottom")

annual.plot




#### TESTING ANNUAL EFFECTS ####
# concentrations
mod.Cc=
  lme(Cc ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Cc)
summary(mod.Cc)
anova(mod.Cc, type = "sequential")

mod.Nc =
  lme(Nc ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Nc)
summary(mod.Nc)
anova(mod.Nc, type = "sequential")

mod.Pc =
  lme(Pc ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Pc)
summary(mod.Pc)
anova(mod.Pc, type = "sequential")

# fluxs
mod.Cflux =
  lme(Cflux.tot ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Cflux)
summary(mod.Cflux)
anova(mod.Cflux, type = "sequential")

mod.Nflux =
  lme(Nflux.tot ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Nflux)
summary(mod.Nflux)
anova(mod.Nflux, type = "sequential")

mod.Pflux =
  lme(Pflux.tot ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.Pflux)
summary(mod.Pflux)
anova(mod.Pflux, type = "sequential")


# ratios
mod.CN=
  lme(CN ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.CN)
summary(mod.CN)
anova(mod.CN, type = "sequential")

mod.CP=
  lme(CP ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.CP)
summary(mod.CP)
anova(mod.CP, type = "sequential")

mod.NP=
  lme(NP ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.flux.annual)

performance::check_model(mod.NP)
summary(mod.NP)
anova(mod.NP, type = "sequential")








### (2) sum per plot and month (across species) & cumulative effects ####

df.flux.plotmonth <-
  df.mass.CNP %>% 
  group_by(plotID, blk, month2, sr, myc, composition) |>
  summarise(Cflux.tot = sum(Cflux, na.rm = T),
            Nflux.tot = sum(Nflux, na.rm = T),
            Pflux.tot = sum(Pflux, na.rm = T),
            mass.tot = sum(mass, na.rm = T)) %>%
  mutate(Cc=Cflux.tot/mass.tot)%>%
  mutate(Nc=Nflux.tot/mass.tot)%>%
  mutate(Pc=Pflux.tot/mass.tot)%>%
  mutate(CN=Cc/Nc)%>%
  mutate(CP=Cc/Pc)%>%
  mutate(NP=Nc/Pc)


df.flux.plotmonth$month2 = factor(df.flux.plotmonth$month2, levels = 
                                    c( month.abb[3:12],  month.abb[1:2]))

df.flux.plotmonth$myc <- factor(df.flux.plotmonth$myc, levels = c("AM", "EM", "AM + EM"))

df.flux.plotmonth =
  df.flux.plotmonth |>
  group_by(plotID, blk, sr, myc, composition) |>
  arrange(month2) |>
  mutate(Cflux.cs = cumsum(Cflux.tot),
         Nflux.cs = cumsum(Nflux.tot),
         Pflux.cs = cumsum(Pflux.tot)) |>
  ungroup()


#### plot cumulative effects  ###
cs.Cflux<- 
  ggplot()+
  geom_point(data = df.flux.plotmonth, 
             aes(x=sr, y=Cflux.cs, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.flux.plotmonth, 
              aes(x=sr, y=Cflux.cs, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month2)+
  labs(y=bquote("Cumulative sum - Cflux"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "none",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))
cs.Cflux

cs.Nflux<- 
  ggplot()+
  geom_point(data = df.flux.plotmonth, 
             aes(x=sr, y= Nflux.cs, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.flux.plotmonth, 
              aes(x=sr, y=Nflux.cs, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month2)+
  labs(y=bquote("Cumulative sum - Nflux"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "none",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))
cs.Nflux


cs.Pflux<- 
  ggplot()+
  geom_point(data = df.flux.plotmonth, 
             aes(x=sr, y= Pc, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.flux.plotmonth, 
              aes(x=sr, y=Pc, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month2)+
  labs(y=bquote("Cumulative sum - Pflux"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "none",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))
cs.Pflux

