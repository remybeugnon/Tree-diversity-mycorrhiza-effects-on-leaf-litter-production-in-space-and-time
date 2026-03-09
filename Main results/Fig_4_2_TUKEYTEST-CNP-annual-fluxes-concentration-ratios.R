#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-11-06
# annual effects (tree diversity and mycorrhizal type) on leaf litter C, N, P concentration, fluxes and CN, CP, NP
# Tukey Test
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

library(emmeans)
emmeans(mod.NP, list(pairwise ~ "myc"), adjust = "tukey")


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



# Get average conc/flux/ratio per mycorrhizal type
mean_conc_myc <- df.flux.annual %>%
  group_by(myc) %>%
  summarise(meanC = mean(NP, na.rm = TRUE), .groups = "drop")

mean_conc_myc_wide <- mean_conc_myc %>%
  pivot_wider(names_from = myc, values_from = meanC)

mean_conc_myc_wide <- mean_conc_myc_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per tree species richness
mean_conc_sr <- df.flux.annual %>%
  group_by(sr) %>%
  summarise(meanC = mean(NP, na.rm = TRUE), .groups = "drop")

mean_conc_sr_wide <- mean_conc_sr %>%
  pivot_wider(names_from = sr, values_from = meanC)

mean_conc_sr_wide <- mean_conc_sr_wide %>%
  mutate(
    `1_vs_2` = (`1` - `2`) / `2` * 100,
    `1_vs_4` = (`1` - `4`) / `4` * 100,
    `2_vs_4` = (`2` - `4`) / `4` * 100)





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


#### plot cumulative effects - line ###
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
  labs(y=bquote("Cumulative sum - C flux"), 
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
  labs(y=bquote("Cumulative sum - N flux"), 
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
  labs(y=bquote("Cumulative sum - P flux"), 
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




##############################################################################
########## PLOT INVERSE OF P-VALUES AS DOTS BASED ON MODEL OUTPUT ############

# Define the abbreviated months in your custom order (March to February)
months <- c("Mar", "Apr", "May", "Jun", "Jul", "Aug", 
            "Sep", "Oct", "Nov", "Dec", "Jan", "Feb")

response_vars <- c("Cc", "Nc", "Pc", "Cflux.tot", "Nflux.tot", "Pflux.tot", "Cflux.cs", "Nflux.cs", "Pflux.cs", "CN", "CP", "NP")

# Initialize an empty data frame to store results
p_values <- data.frame(
  Month = character(),
  Response = character(),
  Variable = character(),
  P_Value = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

table(df.flux.plotmonth$month2) 

df.test <- df.flux.plotmonth %>%
  group_by(month2, plotID)%>%
  filter(month2%in%c("Jun","Jan"))%>%    # in June plot 39, January 17 missing
  summarise()
 
# Loopmodel# Loop over each response variable
for (response in response_vars) {
  for (month in months) {
    dd_month <- subset(df.flux.plotmonth, month2 == month) %>%
      filter(!is.na(get(response)))
    
    if (nrow(dd_month) > 0) {
      # Fit model with tryCatch to handle errors
      model <- tryCatch({
        lme(as.formula(paste(response, "~ I(log2(sr)) * myc")), 
            random = ~1 | composition, 
            data = dd_month)
      }, error = function(e) {
        warning(paste("Error fitting model for month:", month, "and response:", response, ":", e$message))
        return(NULL)
      })
      
      if (is.null(model)) next  # Skip iteration if model fitting failed
      
      # Perform Type I (sequential) ANOVA (default for lme models)
      anova_res <- anova(model)
      
      terms <- c("I(log2(sr))", "myc", "I(log2(sr)):myc")
      
      for (term in terms) {
        if (term %in% rownames(anova_res)) {
          p_val <- anova_res[term, "p-value"]  # Correct column name for p-values in anova.lme()
          
          if (!is.na(p_val) && length(p_val) > 0) {  # Ensure p_val is valid
            sig_level <- ifelse(p_val < 0.05, "Significant", 
                                ifelse(p_val < 0.1, "Marginal", "Not Significant"))
            
            # Append results to p_values
            new_row <- data.frame(
              Month = month,
              Response = response,
              Variable = term,
              P_Value = p_val,
              Significance = sig_level,
              stringsAsFactors = FALSE
            )
            
            p_values <- rbind(p_values, new_row)
          }
        }
      }
    }
  }
}

# Check results
head(p_values)

#-----

# Print the output data frame to see all p-values and significance
print(p_values)

#saveRDS(p_values, file = "1-data/2-1-data-handling/p_values_monthlyCN.rds")

# Convert significance to a factor for plotting
p_values$Significance <- factor(p_values$Significance, 
                                levels = c("Significant", "Marginal", "Not Significant"))

# Rename and reorder the response variables
resp.var_names <- c("Cflux.cs" = "C flux CS",
                    "Cflux.tot" = "C flux",
                    "Cc" = "C concentration",
                    "Nflux.cs" = "N flux CS",
                    "Nflux.tot" = "N flux",
                    "Nc" = "N concentration",
                    "Pflux.cs" = "P flux CS",
                    "Pflux.tot" = "P flux",
                    "Pc" = "P concentration",
                    "CN" = "C : N",
                    "CP" = "C : P",
                    "NP" = "N : P")

# Update the Response column with new names
p_values$Response <- recode(p_values$Response, !!!resp.var_names)

unique(p_values$Response)  

# Rename and reorder the explanatory variables
variable_names <- c("myc" = "Myc", 
                    "I(log2(sr))" = "Sr", 
                    "I(log2(sr)):myc" = "Sr : Myc")

# Update the Variable column with new names
p_values$Variable <- recode(p_values$Variable, !!!variable_names)

# Check unique values to ensure the names match
unique_variables <- unique(p_values$Variable)
print(unique_variables)  # Print to check what levels exist in your data


# Set up levels for month and variable
p_values$Month <- factor(p_values$Month, levels = months)  # Ensure month is in the correct order
p_values$Variable <- factor(p_values$Variable, levels = c("Sr", "Myc", "Sr : Myc"))  # Ensure variable order
p_values$Response <- factor(p_values$Response, levels = c("C flux CS", "C flux", "C concentration", 
                                                          "N flux CS", "N flux", "N concentration", 
                                                          "P flux CS", "P flux", "P concentration",
                                                          "C : N",
                                                          "C : P",
                                                          "N : P"))

# Map significance levels to colors based on response variable
color_mapping <- list(
  'C flux CS' = c("Significant" = "darkred", "Marginal" = "pink", "Not Significant" = "lightgrey"),
  'C flux' = c("Significant" = "darkred", "Marginal" = "pink", "Not Significant" = "lightgrey"),
  'C concentration' = c("Significant" = "darkred", "Marginal" = "pink", "Not Significant" = "lightgrey"),
  'N flux CS' = c("Significant" = "#55828B", "Marginal" = "lightblue", "Not Significant" = "lightgrey"),
  'N flux' = c("Significant" = "#55828B", "Marginal" = "lightblue", "Not Significant" = "lightgrey"),
  'N concentration' = c("Significant" = "#55828B", "Marginal" = "lightblue", "Not Significant" = "lightgrey"),
  'P flux CS' = c("Significant" = "orange", "Marginal" = "peachpuff", "Not Significant" = "lightgrey"),
  'P flux' = c("Significant" = "orange", "Marginal" = "peachpuff", "Not Significant" = "lightgrey"),
  'P concentration' = c("Significant" = "orange", "Marginal" = "peachpuff", "Not Significant" = "lightgrey"),
  'C : N' = c("Significant" = "#874E71", "Marginal" = "#F1CCE3", "Not Significant" = "lightgrey"),
  'C : P' = c("Significant" = "#874E71", "Marginal" = "#F1CCE3", "Not Significant" = "lightgrey"),
  'N : P' = c("Significant" = "#874E71", "Marginal" = "#F1CCE3", "Not Significant" = "lightgrey")
)

# Create a new column for colors based on significance and response variable
p_values <- p_values %>%
  mutate(Color = case_when(
    Significance == "Significant" ~ map_chr(Response, ~ color_mapping[[.x]][["Significant"]] %||% "black"),  # Default color 'black'
    Significance == "Marginal" ~ map_chr(Response, ~ color_mapping[[.x]][["Marginal"]] %||% "grey"),  # Default color 'grey'
    Significance == "Not Significant" ~ map_chr(Response, ~ color_mapping[[.x]][["Not Significant"]] %||% "lightgrey")  # Default color 'lightgrey'
  ))


# Filter out the last three variables (C:N, C:P, N:P ratios)
filtered_data <- p_values %>% 
  filter(!Response %in% c("C : N", "C : P", "N : P"))

# Initialize an empty list for plots
plot_list <- list()

# Iterate over the remaining response variables
for (response in unique(filtered_data$Response)) {
  # Subset data for the current response variable
  plot_data <- subset(filtered_data, Response == response)
  
  if (nrow(plot_data) > 0) {
    # Get color mapping for the current response
    response_colors <- color_mapping[[response]]
    
    # Reclassify p-values into significance levels
    plot_data$Significance <- cut(
      plot_data$P_Value,
      breaks = c(-Inf, 0.05, 0.1, Inf),
      labels = c("Significant", "Marginal", "Not Significant")
    )
    
    # Map size categories to the same significance levels
    plot_data$SizeCategory <- plot_data$Significance
    
    # Create the dot plot
    plot <- ggplot(plot_data, aes(x = Month, y = Variable)) +
      geom_point(
        aes(size = Significance, fill = Significance), 
        shape = 21, color = "black", stroke = 0.2, alpha = 0.8
      ) +
      scale_size_manual(
        values = c("Significant" = 4, "Marginal" = 3, "Not Significant" = 2),
        name = "Significance Level"
      ) +
      scale_fill_manual(
        values = response_colors,
        name = "Significance Level"
      ) +
      facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
      labs(x = "Month", y = "") +
      theme_minimal() +
      theme(
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12),
        panel.spacing.y = unit(0.2, "lines"),
        axis.title.x = element_text(color = "black", size = 12),
        strip.placement = "outside",
        strip.text.y.right = element_text(angle = 0, hjust = 1, vjust = 0.5, color = "black", size = 12),
        strip.background = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        legend.position = ifelse(response == unique(filtered_data$Response)[1], "right", "none")
      )
    
    # Save the plot to the list
    plot_list[[response]] <- plot
  }
}


# Reorder plot_list manually by the desired response order
plot_order <- c("C flux CS", "C flux", "C concentration", 
                "N flux CS", "N flux", "N concentration", 
                "P flux CS", "P flux", "P concentration",
                "C : N",
                "C : P",
                "N : P")

plot_list <- plot_list[plot_order]


library(gridExtra)
F.raster <- grid.arrange(
  grobs = plot_list,  # Use grobs to supply the list of plots
  nrow = 12, ncol = 1  # Arrange with the desired rows and columns
)


plot_list$`C concentration` <- plot_list$`C concentration`  + theme(axis.text.x = element_blank(),
                                                                    axis.title.x = element_blank(),
                                                                    axis.ticks = element_blank(),
                                                                    strip.text.y.right = element_text(size = 8))
plot_list$`N concentration` <- plot_list$`N concentration`  + theme(axis.text.x = element_blank(),
                                                                    axis.title.x = element_blank(),
                                                                    axis.ticks = element_blank(),
                                                                    strip.text.y.right = element_text(size = 8))
plot_list$`P concentration` <- plot_list$`P concentration` + theme(axis.text.x = element_blank(),
                                                                   axis.title.x = element_blank(),
                                                                   axis.ticks = element_blank(),
                                                                   strip.text.y.right = element_text(size = 8))

plot_list$`C flux` <- plot_list$`C flux` + theme(axis.text.x = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 strip.text.y.right = element_text(size = 8))
plot_list$`N flux` <- plot_list$`N flux` + theme(axis.text.x = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 strip.text.y.right = element_text(size = 8))
plot_list$`P flux` <- plot_list$`P flux` + theme(axis.text.x = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 strip.text.y.right = element_text(size = 8))


plot_list$`C flux CS` <- plot_list$`C flux CS` + theme(axis.text.x = element_blank(),
                                                 axis.title.x = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 strip.text.y.right = element_text(size = 8))
plot_list$`N flux CS` <- plot_list$`N flux CS` + theme(axis.text.x = element_blank(),
                                                       axis.title.x = element_blank(),
                                                       axis.ticks = element_blank(),
                                                       strip.text.y.right = element_text(size = 8))
plot_list$`P flux CS` <- plot_list$`P flux CS` + theme(axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
                                                       axis.title.x = element_text(size = 12),
                                                       axis.ticks = element_blank(),
                                                       strip.text.y.right = element_text(size = 8))


library(ggplot2)
library(patchwork)


all_plots <- (
    plot_list$`C concentration` / 
    plot_list$`N concentration` / 
    plot_list$`P concentration` / 
    
    plot_list$`C flux` / 
    plot_list$`N flux` / 
    plot_list$`P flux` /  
    
    plot_list$`C flux CS` /
    plot_list$`N flux CS` /
    plot_list$`P flux CS` 
)

all_plots  # excluding ratios, since they have missing months


### plot to check the direction of treatment effects 

ggplot()+
  geom_point(data = df.flux.plotmonth, 
             aes(x=sr, y=Pflux.tot, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.flux.plotmonth, 
              aes(x=sr, y=Pflux.tot, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month2)+
  labs(y= "Leaf litter P flux", 
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
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

