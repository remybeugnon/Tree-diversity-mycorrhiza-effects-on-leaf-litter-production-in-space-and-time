#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-09-05
# dead trees, average tree height, tree dbh
# tree inventory data from MyDiv in 2023
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)
# anova Type  I

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)

#============================ Dataset ===============================


tree_inv_23 <- read_csv2

# Shows structure with variable names and types
str(tree_inv_23)
levels(tree_inv_23$vitality)
unique(tree_inv_23$vitality)

library(devtools)
library(httr)
library(jsonlite)
library(XML)
library(rBExIS)
load_all("rBExIS")
bexis.options("base_url" = "https://mydivdata.idiv.de")
bexis.get.datasets()

plotinfo <- rBExIS::bexis.GetDatasetById(id =43)



dead_trees <- tree_inv_23 %>%
  count(vitality) %>%
  mutate(percent = n / sum(n) * 100)




# 1. Total number of dead trees
spc_tbl <- tree_inv_23  %>%
  filter(vitality == "dead") %>%
  summarise(total_dead = n())

spc_tbl <- tree_inv_23  %>%
  summarise(
    total_trees = n(),
    dead_trees = sum(vitality == "dead"),
    percent_dead = dead_trees / total_trees * 100
  )


# 2. Dead trees split by treatments (after merging treatment info)
# Assuming you have treatment_tbl with plotID, sr, myc

inven_join <- full_join(plotinfo, tree_inv_23, by=c("plotID","block"))

inven_join$div<-as.factor(inven_join$tree_species_richness)
inven_join$blk<-as.factor(inven_join$block)
inven_join$myc<-as.factor(inven_join$mycorrhizal_type)
inven_join$sr<-inven_join$tree_species_richness

spc_tbl_dead <- inven_join  %>%
  filter(vitality == "dead") %>%
  summarise(total_dead = n())


spc_tbl_treat <- inven_join  %>%
  group_by(sr, myc) %>%
  summarise(
    total_trees = n(),
    dead_trees = sum(vitality == "dead"),
    percent_dead = dead_trees / total_trees * 100,
    .groups = "drop"
  )

# 3. Dead trees split by species
spc_tbl_treat_spec <- inven_join  %>%
  group_by(species) %>%
  summarise(
    total_trees = n(),
    dead_trees = sum(vitality == "dead"),
    percent_dead = dead_trees / total_trees * 100,
    .groups = "drop"
  )


install.packages("writexl")   # run once if not installed
library(writexl)

# save plot-level table
write_xlsx(spc_tbl_treat_spec, "C:/Users/eb64wupi/Desktop/Litterfall Journal of Ecology/First Revision Oct 2025/Tree_species_dead_2025-09-05.xlsx")
write_xlsx(spc_tbl_treat, "C:/Users/eb64wupi/Desktop/Litterfall Journal of Ecology/First Revision Oct 2025/Tree_species_per_treatment_dead_2025-09-05.xlsx")



####################
# Reviewer comment: It would be great to know some information about the sizes of the trees during your experiment. Maximum or average height or diameter, for instance, would give the readers a sense of the forest.

spc_tbl_sizes <- inven_join %>%
  mutate(
    height_num = as.numeric(height),
    stem01_dbh_num = as.numeric(stem01_dbh),
    stem02_dbh_num = as.numeric(stem02_dbh),
    stem03_dbh_num = as.numeric(stem03_dbh)
  ) %>%
  rowwise() %>%
  mutate(
    # combine stems into a list per row
    stem_values = list(c(stem01_dbh_num, stem02_dbh_num, stem03_dbh_num)),
    # remove NA and 0 for min and max calculation
    stem_nonzero = list(stem_values[[1]][!is.na(stem_values[[1]]) & stem_values[[1]] > 0]),
    max_dbh = ifelse(length(stem_nonzero[[1]]) > 0, max(stem_nonzero[[1]]), NA),
    min_dbh = ifelse(length(stem_nonzero[[1]]) > 0, min(stem_nonzero[[1]]), NA)
  ) %>%
  ungroup() %>%
  select(-stem_values, -stem_nonzero)


# 1. Overall summary
spc_tbl_sizes %>%
  summarise(
    avg_height = mean(height_num, na.rm = TRUE),
    min_height = min(height_num, na.rm = TRUE),
    max_height = max(height_num, na.rm = TRUE),
    avg_dbh = mean(max_dbh, na.rm = TRUE),
    min_dbh = max(min_dbh, na.rm = TRUE),
    max_dbh = max(max_dbh, na.rm = TRUE)
  )

# 2. Summary per species
species_size <- spc_tbl_sizes %>%
  group_by(species) %>%
  summarise(
    avg_height = mean(height_num, na.rm = TRUE),
    min_height = min(height_num, na.rm = TRUE),
    max_height = max(height_num, na.rm = TRUE),
    avg_dbh = mean(max_dbh, na.rm = TRUE),
    min_dbh = max(min_dbh, na.rm = TRUE),
    max_dbh = max(max_dbh, na.rm = TRUE),
    .groups = "drop"
  )

# save plot-level table
write_xlsx(species_size, "C:/Users/eb64wupi/Desktop/Litterfall Journal of Ecology/First Revision Oct 2025/Tree_species_height-dbh_2025-09-05.xlsx")

# 3. Summary per treatment (if you have sr & myc)
# Join sizes with treatment info
spc_tbl_sizes_treat <- spc_tbl_sizes %>%
  left_join(inven_join, by = "plotID", relationship = "many-to-many")

# Now summarise per treatment
spc_tbl_sizes_treat2 <- spc_tbl_sizes_treat %>%
  group_by(sr, myc) %>%
  summarise(
    avg_height = mean(height_num, na.rm = TRUE),
    max_height = max(height_num, na.rm = TRUE),
    avg_dbh = mean(max_dbh, na.rm = TRUE),
    max_dbh = max(max_dbh, na.rm = TRUE),
    .groups = "drop")



