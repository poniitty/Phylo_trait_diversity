
library(tidyverse)
library(ape)
library(V.PhyloMaker)
library(readxl)

# Species names
spn <- read_csv("comm_data/Species_list_Rissanen.csv")
# Large nomenclature database for Fennoscandian flora
all_names <- read_excel("nomenclature/NamesFromDatabasesFINAL.xlsx")

# Remove genus level taxa
spn %>% filter(grepl("_",species)) -> spn

all_names %>% mutate(species = gsub(" ","_",final)) %>% 
  select(species, OldName, final, genus, family) %>% 
  mutate(OldName = gsub(" ","_",OldName)) -> all_names

# Check if all the study species in the nomenclature data
spn$species[!which(spn$species %in% all_names$species)]

# Join datasets
spn <- left_join(spn, all_names %>% filter(!duplicated(all_names %>% select(species))))

# Function to replace second occurence of specified character string
# e.g. "Silene_flos_cuculi" to "Silene_flos-cuculi"
replace_second_ <- function(x, strch, repl_chr){
  unlist(lapply(x, function(x){
    strs <- strsplit(x, strch)[[1]]
    
    if(length(strs) > 2){
      if(length(strs) > 3){
        return(paste0(strs[1],strch,paste(strs[2:3], collapse = repl_chr),strch,paste(strs[4:length(strs)], collapse = strch)))
      } else {
        return(paste0(strs[1],strch,paste(strs[2:3], collapse = repl_chr)))
      }
    } else {
      return(x)
    }
  }))
}
spn %>% mutate(species = replace_second_(species, strch = "_", repl_chr = "-")) -> spn

#  This is the vacular plant mega three from V.PhyloMaker
data(GBOTB.extended)

# Which study species names not found in the phylogenetic tree
spn$species[which(!spn$species %in% GBOTB.extended$tip.label)]


# Try to find macthes from known synonyms
spn %>% filter(!species %in% GBOTB.extended$tip.label) %>% 
  pull(species) -> missing_sp

all_names %>% filter(species %in% missing_sp) -> missing

missing_sp2 <- missing$OldName[which(missing$OldName %in% GBOTB.extended$tip.label)]

for(i in missing_sp2){
  print(i)
  
  all_names %>% filter(OldName == i) %>% 
    pull(species) -> temp_sp
  
  if(length(temp_sp) == 1){
    
    spn %>% mutate(species = ifelse(species == temp_sp,
                                    i, species)) -> spn
    
  } else {
    
    print("Juuh elikkäs?!?!?!")
    stop()
    
  }
  
}

# Which one are still missing?
# Consider going these through manually
spn$species[which(!spn$species %in% GBOTB.extended$tip.label)]

# GBOTB.extended$tip.label[grepl("glacialis",GBOTB.extended$tip.label)]
# td$species[grepl("glacialis",td$species)]
GBOTB.extended$tip.label[grepl("Dryopteris_", GBOTB.extended$tip.label)] %>% sort()

# This is how you can change a single species name in the data
spn %>% mutate(species = ifelse(species == "Silene_flos_cuculi", "Silene_flos-cuculi", species)) -> spn

# Create dataset with study species' taxonomies to create the subtree from mega phylogeny
taxonom <- data.frame(species = spn$species, 
                      genus = unlist(lapply(spn$species, function(x) strsplit(x, "_")[[1]][1])), 
                      family = spn$family)

# Create the subtree
tree1 <- phylo.maker(taxonom, tree = GBOTB.extended, nodes = nodes.info.1, scenarios=c("S3"), r = 1)

plot(tree1$scenario.3)

# Check if all study species in the subtree
spn$species[!spn$species %in% tree1$scenario.3$tip.label]

# Save as Rdata
save(tree1, file = "phylo_trees/phylo_tree_rissanen.Rdata")



