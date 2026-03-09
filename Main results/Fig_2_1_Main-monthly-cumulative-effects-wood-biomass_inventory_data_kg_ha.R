#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-05-21
# Fig. 2: Annual, monthly and cumulative litter fall + wood biomass increment (productivity)
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

df.all.wide.info <- read_csv

# g/m2 in kg/ha

df.all.wide.info.kgha <- df.all.wide.info %>%
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


# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info.kgha %>%
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

df.all.wide.mean$div<-as.factor(df.all.wide.mean$tree_species_richness)
df.all.wide.mean$blk<-as.factor(df.all.wide.mean$block)
df.all.wide.mean$myc<-as.factor(df.all.wide.mean$mycorrhizal_type)
df.all.wide.mean$sr<-df.all.wide.mean$tree_species_richness
df.all.wide.mean$sr_myc<-paste(df.all.wide.mean$sr,df.all.wide.mean$myc,sep="_")#interaction terms
df.all.wide.mean$myc <- recode_factor(df.all.wide.mean$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")


#### !!!! check if the selcetd columns are correct
# 2) long format ####
# long format ####
df.mass.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c("Ac":"Ti"),
               names_to="species",
               values_to="dryweight")

# exclude rows with NAs/missing data
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
  ungroup() 

###### main effect ######
df.annual.litter = 
  df.mass.long |> 
  group_by(plotID, myc, tree_species_richness, block, composition) |>
  summarise(litter.prod = sum(dryweight, na.rm = T), 
            sd.litter.prod = sd(dryweight, na.rm = T)) 

main.eff <- ggplot(df.annual.litter, 
                   aes(x=tree_species_richness, y=litter.prod, 
                       color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1, height = 0)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(
    y = expression(atop("Annual litterfall", (kg ~ ha^-1))),  # Use atop() for line break
    x = "Tree species richness"
  ) +
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
        strip.text = element_text(color="black", size = 13),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 13),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text(size = 13),
        legend.text = element_text(size=13),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")

main.eff

mod.litter.bio=
  lme(litter.prod ~ log2(tree_species_richness) * myc,
      random = ~1|composition,
      data = df.annual.litter)

performance::check_model(mod.litter.bio)
summary(mod.litter.bio)
anova(mod.litter.bio, type = "sequential")





###### monthly effect ######

df.monthly.litter = df.mass.long |> 
  group_by(block, plotID, sr, div, myc, composition, month1) |> 
  summarise(litter.prod = sum(dryweight, na.rm = T)) 


df.monthly.litter$month2 <- dplyr::recode_factor(df.monthly.litter$month1, 
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

df.monthly.litter$month1 = factor(df.monthly.litter$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))

M = map_df( .x = unique(df.monthly.litter$month1),
            .f = ~ {
              mod = 
                lme(litter.prod ~ log2(sr) * myc, 
                    random= ~1|composition,
                    data= df.monthly.litter |>
                      filter(month1 == .x)) |> 
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


fig.month <- ggplot()+
  geom_jitter(data = df.monthly.litter, 
             aes(x=sr, y=litter.prod, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5, width=0.1, height = 0)+
  geom_smooth(data = df.monthly.litter, 
              aes(x=sr, y=litter.prod, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(
    y = expression(atop("Monthly litterfall", (kg~ha^-1))),  # Use atop() for line break
    x = "Tree species richness"
  ) +
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 2700),
                     breaks= c(500,1000,1500,2000,2500,3000))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide = "none")+
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr)'),
            aes(x=1, y=2700,
                label = sign,
                hjust = 0),
            size=5,
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=2600,
                label = sign,
                hjust = 0),
            size=5,
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr):myc'),
            aes(x=1, y=2500,
                label = sign,
                hjust = 0),
            size=5,
            color = 'black') +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 13),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 13),
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
fig.month


##### cumulative effect ######

df.cum.litter = df.mass.long |> 
  group_by(block, plotID, div, sr, myc, composition, month1) |> 
  summarise(litter.prod = sum(dryweight, na.rm = T))

df.cum.litter$month2 <- dplyr::recode_factor(df.cum.litter$month1, 
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

df.cum.litter$month1 = factor(df.cum.litter$month2, 
                                  levels = c(month.abb[3:12], month.abb[1:2]))
df.cum.litter.1 = 
  df.cum.litter |> 
  group_by(block, plotID, div, sr, myc, composition) |> 
  arrange(month1) |> 
  mutate(cs = cumsum(litter.prod)) 

M = map_df( .x = unique(df.cum.litter.1$month1),
            .f = ~ {
              mod = 
                lme(cs ~ log2(sr) * myc, 
                    random= ~1|composition,
                    data= df.cum.litter.1 |>
                      filter(month1 == .x)) |> 
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

# 5) plot cumulative effects - facets
cum_facet<- ggplot()+
  geom_jitter(data = df.cum.litter.1, 
             aes(x=sr, y=cs, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5, width=0.1, height = 0)+
  geom_smooth(data = df.cum.litter.1, 
              aes(x=sr, y=cs, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(
    y = expression(atop("Litterfall (cumulative)", (kg ~ ha^-1))),  # Use atop() for line break
    x = "Tree species richness"
  ) +
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_y_continuous(limits = c(0, 5800),
                     breaks=c(0, 1000, 2000, 3000, 4000, 5000))+
  #coord_cartesian(ylim = c(0, 140), xlim = c(1,4))+ 
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr)'),
            aes(x=1, y=5800,
                label = sign,
                #fontface = "bold",
                hjust = 0),
            size=5,
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=5650,
                label = sign,

                #fontface = "bold",
                hjust = 0),
            size=5,
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr):myc'),
            aes(x=1, y=5500,
                label = sign,
                #fontface = "bold",
                hjust = 0),
            size=5,
            color = 'black') +
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(color="black", size = 13),
        axis.line = element_line(color="black"),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks.x = element_line(color="black"),
        strip.text.x = element_text(color="black", size = 13),
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
  labs(tag = "(e)")
cum_facet





### wood ###########

#============================ Packages ===============================

# from MyDiv database
library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()

#============================ Datasets ===============================

wood <- read_csv

wood2 <- wood %>%
  dplyr::rename(wood_bio_22 = "2022",
                wood_bio_23 = "2023")

# t/ha in kg/ha

wood2 <- wood2 %>%
  dplyr::mutate(wood_bio_22 = wood_bio_22*1000,
                wood_bio_23 = wood_bio_23*1000,
                increment = increment*1000)

plotinfo <- rBExIS::bexis.GetDatasetById(id =43)

# join the datasets
wood3 <- full_join(wood2, plotinfo, by=c("plotID","plotName","block","tree_species_richness","mycorrhizal_type","composition"))

wood3$myc <- recode_factor(wood3$mycorrhizal_type, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

# join the datasets   (!!!!! SOMETHING DOESNT WORK HERE)
wood_litter <- full_join(wood3, df.annual.litter, by=c("plotID","block","tree_species_richness","myc","composition"))

wood_litter$div<-as.factor(wood_litter$tree_species_richness)
wood_litter$blk<-as.factor(wood_litter$block)
wood_litter$myc<-as.factor(wood_litter$mycorrhizal_type)
wood_litter$sr<-wood_litter$tree_species_richness
wood_litter$myc <- recode_factor(wood_litter$mycorrhizal_type, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")


#============================ Plotting ===============================

# 2022
F1_2022 <- ggplot(wood_litter, 
                  aes(x=sr, y=wood_bio_22, 
                      color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(
    y = expression(atop("Wood biomass 2022", (kg ~ ha^-2))),  # Use atop() for line break
    x = "Tree species richness"
  ) +
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                    guide="none"
                    )+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=13),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(13),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
F1_2022


mod.tree.bio.22 =
  lme(wood_bio_22 ~ log2(sr) * myc,
      random = ~1|composition,
      data = wood_litter)

#performance::check_model(mod.tree.bio.22)
summary(mod.tree.bio.22)
anova(mod.tree.bio.22, type = "sequential")



# 2023
F1_2023 <- ggplot(wood_litter, 
                  aes(x=sr, y=wood_bio_23, 
                      color = myc, fill = myc))+
  geom_jitter(shape =21, size = 1, alpha=0.5, width= 0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(
    y = expression(atop("Wood biomass 2023", (kg ~ ha^-2))),  # Use atop() for line break
    x = "Tree species richness"
  ) +
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none"
  )+
  theme_minimal()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=13),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(13),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =13),
        plot.subtitle = element_text(size=13),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(b)")
F1_2023


mod.tree.bio.23 =
  lme(wood_bio_23 ~ log2(sr) * myc,
      random = ~1|composition,
      data = wood_litter)

#performance::check_model(mod.tree.bio.22)
summary(mod.tree.bio.23)
anova(mod.tree.bio.23, type = "sequential")


library(patchwork)

comb <- (F1_2022 | F1_2023) +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom") +
  plot_layout(heights = c(2, 5))  # Allocate more height to the bottom rows

comb


### plotting the increment ### 

incr <- ggplot(wood_litter, 
                    aes(x=sr, y=increment, 
                        color = myc, fill = myc)) +
  geom_jitter(shape = 21, size = 1, alpha=0.5, width=0.1, height=0) +
  geom_smooth(method="lm", alpha=0.3) +
  labs(y=expression(atop("Annual wood increment", (kg~ha^-1))), #(kg~ha^-1~yr^-1)
       x = "Tree species richness")+  #(kg~ha^-1)
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4)) +
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none") + 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none") +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=13),
    axis.line = element_line(color='black'),
    axis.text.y = element_text(color="black", size = 12),
    axis.text.x = element_text(color="black", size = 12),
    axis.title.y = element_text(color="black", size = 13, angle = 90),
    axis.title.x = element_text(color="black", size = 13),
    axis.ticks = element_line(color="black"),
    strip.text.x = element_text(size=13),
    panel.border = element_rect(colour="black", fill=NA),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 13),
    plot.subtitle = element_text(size = 13),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key = element_rect(color="transparent"),   
    legend.title = element_text(size = 13),
    legend.text = element_text(size=13),
    legend.background = element_rect(colour=NA),
    legend.box.background = element_rect(color="transparent")
  ) +
  labs(tag = "(b)")

incr

mod.incr =
  lme(increment ~ log2(sr) * myc,
      random = ~1|composition,
      data = wood_litter)

#performance::check_model(mod.incr)
summary(mod.incr)
anova(mod.incr, type = "sequential")





###### litter and wood biomass ####

# annual litter biomass + annual wood increment
wood_litter2 <- wood_litter |>
  # mutate(bio_increment_plot_kgha = bio_increment_plot*10000,
  #        bio_increment_tree_kgha = bio_increment_tree*10000) |>
  mutate(total_bio_tree = increment + litter.prod)

Fig_lit_wood <-
ggplot(wood_litter2, 
       aes(x=sr, y=total_bio_tree, color=myc, fill=myc))+  # per plot or per tree???
  geom_jitter(shape =21, size = 1, alpha=0.5, width=0.1)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(
    y = expression(atop("Annual biomass production", "(wood + litter) (kg " * ha^-1 * ")")),
    x = "Tree species richness"
  )+
  # labs(
  #   y = expression(atop("Total biomass", atop("(wood + litter)", (kg ~ m^-3 ~ yr^-1)))),  # Use atop() for line break
  #   x = "Tree species richness"
  # ) +
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
        strip.text = element_text(size=13),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(color="black", size = 13),
        axis.title.x = element_text(color="black", size = 13),
        axis.ticks = element_line(color="black"),
        strip.text.x = element_text(13),
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
  labs(tag = "(a)")

Fig_lit_wood

mod.total.bio =
  lme(total_bio_tree ~ log2(sr) * myc,
      random = ~1|composition/block,
      data = wood_litter2)

#performance::check_model(mod.total.bio)
summary(mod.total.bio)
anova(mod.total.bio, type = "sequential")





### all plots together ####
library(patchwork)

top_plots <- (Fig_lit_wood | incr | main.eff)

bottom_plots <- (fig.month / cum_facet)


combined_plot <- (top_plots / bottom_plots) +
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom") +
  plot_layout(heights = c(2, 5))  # Allocate more height to the bottom rows

combined_plot

### end ###


