library(tidyverse)
# devtools::install_github("ibartomeus/fundiv")
library(fundiv)
library(readxl)
library(doParallel)
library(scales)


td <- read_csv("trait_data/median_traits_raw.csv")
id <- read_csv("trait_data/imputed_median_traits_phylopars.csv")
all_names <- read_excel("nomenclature/NamesFromDatabasesFINAL.xlsx")
all_names %>% mutate(species = gsub(" ","_",final)) %>% 
  select(species, OldName, final, genus, family) %>% 
  mutate(OldName = gsub(" ","_",OldName)) -> all_names


# Prepare trait data
td %>% select(species:`Leaf phosphorus`) %>% 
  pivot_longer(Height:`Leaf phosphorus`,names_to = "trait") %>% 
  filter(!is.na(value)) %>% 
  mutate(type = "median") -> td

id %>% select(species:`Leaf phosphorus`) %>% 
  pivot_longer(Height:`Leaf phosphorus`,names_to = "trait") %>% 
  filter(!is.na(value)) %>% 
  mutate(type = "imputed")  -> id

d <- bind_rows(td, id) %>% 
  arrange(species)

# Extract species names not found in final names of nomenclature table
confl_trait_names <- unique(d$species)[!unique(d$species) %in% all_names$species]

# Function to check nomenclature and find the most recent name if synonyms
pull_newname <- function(xx){(unlist(lapply(xx, function(x){
  all_names %>% filter(OldName == x) %>% pull(species) -> temp_sp
  if(length(temp_sp)){
    return(temp_sp)
  } else {
    return("JUUH ELIKKÄSSSSSSSSSSSSSSS")
  }
})))}

# Change conflicting names
for(i in confl_trait_names){
  print(i)
  
  d %>% mutate(species = ifelse(species == i, pull_newname(i), species)) -> d
  
}

# Should be now length of zero
unique(d$species)[!unique(d$species) %in% all_names$species]


###############################################################
# Read community data

cd <- read_csv("comm_data/Community_data_aa.csv")

# There are both Deschampsia_atropurpurea and Vahlodea_atropurpurea in the original comm data
#  Only one of those should be selected
sum(cd$Vahlodea_atropurpurea)
sum(cd$Deschampsia_atropurpurea)
# Vahlodea_atropurpurea selected
# cd %>% select(-Deschampsia_atropurpurea) -> cd

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
names(cd)[5:ncol(cd)] <- replace_second_(names(cd)[5:ncol(cd)], strch = "_", repl_chr = "-")

# Exctract names of the study species
sp_names <- names(cd)[2:ncol(cd)]

# Check which one are missing from the trait data
missing_traits_sp <- sp_names[!sp_names %in% d$species]
missing_traits_sp[missing_traits_sp %in% all_names$OldName]

# Function to check nomenclature and find the most recent name if synonyms
pull_newname <- function(xx){(unlist(lapply(xx, function(x){
  all_names %>% filter(OldName == x) %>% pull(species) -> temp_sp
  if(length(temp_sp)){
    return(temp_sp)
  } else {
    return("JUUH ELIKKÄSSSSSSSSSSSSSSS")
  }
})))}

for(i in missing_traits_sp){
  
  names(cd)[which(names(cd) == i)] <- pull_newname(i)
  
}

sp_names <- names(cd)[2:ncol(cd)]

missing_traits_sp <- sp_names[!sp_names %in% d$species] 
missing_traits_sp
# 19 species has no trait data, these will be dropped from further trait based analyses

cd %>% select(-missing_traits_sp) -> cd


#####################################################################################
# FUNCTIONAL DIVERSITY INDICES

d %>% filter(species %in% names(cd)) -> d

# Check the percentage of species level trait values that have been imputed
d %>% group_by(trait, type) %>% 
  summarise(n = n()) %>% 
  mutate(n = n/length(unique(d$species))*100)

# Log transfrom some of the traits
d %>% select(-type) %>% 
  mutate(trait = gsub(" ","_", trait)) %>% 
  pivot_wider(id_cols = species, names_from = trait, values_from = value) %>% 
  mutate(Height = log(Height),
         Seed_mass = log(Seed_mass),
         Leaf_area = log(Leaf_area)) -> traits

# Rescale all the traits
traits %>% mutate(across(Height:Leaf_phosphorus, ~as.numeric(scale(.x)))) -> traits

# Check how the traits are interrelated
summary(prcomp(traits[,2:ncol(traits)]))
cor(traits[,2:NCOL(traits)])

# These steps are to make sure that the community and trait data are at similar order
sp_names <- names(cd)[2:ncol(cd)]
sp_names <- sort(sp_names)

traits %>% arrange(species) %>% 
  as.data.frame() -> traits

# FD functions need trait names as row names
rownames(traits) <- traits$species

# TEST IF IDENTICAL and in the same order
identical(sp_names, traits$species) # Should be TRUE
min(sp_names == traits$species) # Should be 1

# TRAIT DENDROGRAM
trait.clust <- hclust(dist(traits[,2:NCOL(traits)]), method="average")
# plot(trait.clust)

# Function to calculate multiple FD indices with multiple CPUs
fd_parallel <- function(com, trait){
  
  com <- com[,rownames(trait)]
  
  num_splits<-cores*1
  split_testing<-sort(rank(1:nrow(com))%%num_splits)
  
  jo <- foreach(ii=unique(split_testing), .combine=rbind,
                .packages = c("fundiv")) %dopar% {
                  FD_dendro(trait, com[split_testing == ii,],
                            Cluster.method = "average", ord = "podani",
                            stand.x = F)
                }
  jo2 <- foreach(ii=unique(split_testing), .combine=rbind,
                 .packages = c("fundiv")) %dopar% {
                   tempdf <- com[split_testing == ii,]
                   tempdf <- tempdf[,colSums(tempdf > 0) > 0]
                   temptrait <- trait[rownames(trait) %in% names(tempdf),]
                   as.data.frame(dbFD(temptrait, tempdf, w.abun = TRUE, stand.x = F,
                                      ord = c("podani"), asym.bin = NULL,
                                      corr = c("cailliez"), calc.FRic = TRUE, 
                                      m = 3L, stand.FRic = FALSE,
                                      scale.RaoQ = FALSE, clust.type = "ward",
                                      calc.CWM = T, calc.FDiv = TRUE, dist.bin = 2, 
                                      print.pco = F, messages = F))
                   
                 }
  
  all <- bind_cols(jo, jo2)
  
  return(as.data.frame(all))
}

# Set the number of CPU cores for calculations
cores <- 1

# Calculate indices in parallel
cl <- parallel::makeCluster(cores)
doParallel::registerDoParallel(cl)
true <- fd_parallel(com = cd[rowSums(cd[,sp_names] > 0) > 0,sp_names],
                    trait = traits[,2:NCOL(traits)])
parallel::stopCluster(cl)


# Variables of special interest
keep <- c("FDpg","FDw","FRic","FEve","FDis","RaoQ")

# Combine with plot ids etc...
true <- bind_cols(cd[rowSums(cd[,sp_names] > 0) > 0,1], true)

# test which are strongly correlated
round(cor(true[,c("n_sp", keep)], use = "pairwise.complete.obs"),3)


write_csv(true, "outputs/FD_indices_AA.csv")


