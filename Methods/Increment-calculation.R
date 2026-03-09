# MyDiv Litterfall experiment 
# wood biomass increment 
# Elisabeth Bönisch (boenisch.elisabeth@idiv.de)

# plot size
# 64 m2 

# treatments, plot names, 
library(tidyverse)
library(readxl)
data_inventory <- read_excel #use inventory data from 2022, 2023 (wood biomass)

# from MyDiv database
library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()

plotinfo <- rBExIS::bexis.GetDatasetById(id =43)


data.density = data.frame(
  species = c("Prunus avium", "Tilia platyphyllos", "Fagus sylvatica", "Fraxinus excelsior", 
              "Sorbus aucuparia", "Quercus petraea", "Betula pendula", "Acer pseudoplatanus",
              "Aesculus hippocastanum", "Carpinus betulus"),
  density = c(0.620, 0.475, 0.710, 0.680, 0.610, 0.710, 0.640, 0.615, 0.500, 0.735)
)

df = data_inventory |> 
  left_join(data.density) |> 
  mutate(BA = 3.1415926535 * ((as.numeric(d0)/2)^2)) |>
  mutate(biomass = 0.42 * (as.numeric(d0)^2) * as.numeric(height) * 100 * density / 1000000)

df.1 = 
  df |> 
  group_by(year, block, plotID) |> 
  summarise(biomass.ha = sum(biomass, na.rm = T) / 64 * 10000) |> 
  pivot_wider(names_from = year, values_from = biomass.ha, id_cols = plotID) |> 
  mutate(increment = `2023` - `2022`)


df.final = 
  plotinfo |> 
  left_join(df.1 |> select(plotID, `2022`, `2023`, increment))


ggpubr::ggarrange(
  ggplot(data = df.final,
         aes(x = tree_species_richness, color = mycorrhizal_type, y = `2022`)) + 
    geom_jitter(width = .05) + 
    geom_smooth(method = 'lm') + 
    scale_x_continuous(trans = 'log'),
  ggplot(data = df.final,
         aes(x = tree_species_richness, color = mycorrhizal_type, y = `2023`)) + 
    geom_jitter(width = .05) + 
    geom_smooth(method = 'lm') + 
    scale_x_continuous(trans = 'log'),
  ggplot(data = df.final,
         aes(x = tree_species_richness, color = mycorrhizal_type, y = increment)) + 
    geom_jitter(width = .05) + 
    geom_smooth(method = 'lm') + 
    scale_x_continuous(trans = 'log'),
  nrow = 3
)

write.csv(df.final, 'C:/Users/eb64wupi/Documents/MyDiv-Litterfall-project/1-data/Increment June 2025/increment-data-updated.csv', row.names = F)
