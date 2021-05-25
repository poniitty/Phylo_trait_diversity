
library(tidyverse)
library(ape)
library(V.PhyloMaker)
library(readxl)

spn <- read_csv("comm_data/Species_list_Rissanen.csv")
all_names <- read_excel("nomenclature/NamesFromDatabasesFINAL.xlsx")

# Remove genus level taxa
spn %>% filter(grepl("_",species)) -> spn

all_names %>% mutate(species = gsub(" ","_",final)) %>% 
  select(species, OldName, final, genus, family) -> all_names

spn$species[!which(spn$species %in% all_names$species)]

spn <- left_join(spn, all_names %>% filter(!duplicated(all_names %>% select(species))))

data(GBOTB.extended)

spn$species[which(!spn$species %in% GBOTB.extended$tip.label)]

# GBOTB.extended$tip.label[grepl("glacialis",GBOTB.extended$tip.label)]
# td$species[grepl("glacialis",td$species)]
# GBOTB.extended$tip.label[grepl("Woodsia_", GBOTB.extended$tip.label)]
spn %>% mutate(species = ifelse(species == "Anthoxanthum_monticola", "Hierochloe_alpina", species),
              species = ifelse(species == "Anthoxanthum_nipponicum", "Anthoxanthum_alpinum", species),
              species = ifelse(species == "Cherleria_biflora", "Minuartia_biflora", species),
              species = ifelse(species == "Omalotheca_norvegica", "Gnaphalium_norvegicum", species),
              species = ifelse(species == "Omalotheca_supina", "Gnaphalium_supinum", species),
              species = ifelse(species == "Athyrium_distentifolium", "Athyrium_alpestre", species),
              species = ifelse(species == "Avenella_flexuosa", "Deschampsia_flexuosa", species),
              species = ifelse(species == "Harrimanella_hypnoides", "Cassiope_hypnoides", species),
              species = ifelse(species == "Kalmia_procumbens", "Loiseleuria_procumbens", species),
              species = ifelse(species == "Spinulum_annotinum", "Lycopodium_annotinum", species),
              species = ifelse(species == "Vahlodea_atropurpurea", "Deschampsia_atropurpurea", species),
              species = ifelse(species == "Bistorta_vivipara", "Persicaria_vivipara", species),
              species = ifelse(species == "Ranunculus_glacialis", "Beckwithia_glacialis", species)) -> spn


taxonom <- data.frame(species = spn$species, 
                      genus = unlist(lapply(spn$species, function(x) strsplit(x, "_")[[1]][1])), 
                      family = spn$family)

tree1 <- phylo.maker(taxonom, tree = GBOTB.extended, nodes = nodes.info.1, scenarios=c("S3"), r = 1)

plot(tree1$scenario.3)

spn$species[!spn$species %in% tree1$scenario.3$tip.label]

save(tree1, file = "phylo_trees/phylo_tree_rissanen.Rdata")



