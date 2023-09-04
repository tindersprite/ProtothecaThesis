#### Setup #####
#### Packages ####
library(tidyverse)
library(ComplexUpset)

#### Functions ####
remove_spaces <- function(x){
  return(gsub(" ", "", x))
  }

#### Misc ####
data_files <- Sys.glob("Data/busco_lists/*")

strains <- read.csv("Data/busco_scores.csv")
strains <- strains[1:2]
strains$Strain <- remove_spaces(strains$Strain)

genera <- c("Auxenochlorella",
            "Chlorella",
            "Helicosporidium",
            "Micractinium",
            "Prototheca")
names(genera) <- substring(genera,1,1)

lineage <- c("AHP",
             "CHAMP",
             "AHP",
             "CHAMP",
             "AHP")
names(lineage) <- genera
#####



#####
#### Load all BUSCO data ####
for(i in 1:length(data_files)){
   buscos_for_strain <- read.delim(data_files[i], 
           header = T,
           skip = 2)
   buscos_for_strain <- buscos_for_strain[1:2]
   
   names(buscos_for_strain)[1] <- "BuscoID"
   
   buscos_for_strain <- unique(buscos_for_strain)
   
   strain <- strsplit(data_files[i], split = "/")[[1]][3]
   strain <- strsplit(strain, split = "_")[[1]][1]
   
   buscos_for_strain$Strain <- strain
   
   if(i == 1){
     all_buscos <- buscos_for_strain
   } else {
     all_buscos <- rbind(all_buscos, buscos_for_strain)
   }
}


#### Formatting for species upset ####  
### Identify species for each strain
all_buscos <- left_join(all_buscos, strains)

### Choose which genes count as being present
logical_buscos <- all_buscos
logical_buscos$Status <- logical_buscos$Status == "Complete" | logical_buscos$Status == "Duplicated" | logical_buscos$Status == "Fragmented" 

#### Investigations of Species ####
### Identify the species for which there is only one genome, and thus cannot be 
### used to find a union/intersect
unique_species <- unique(logical_buscos$Species)
n.species <- length(unique_species)
species_to_combine <- c()
consolidated_logical_buscos <- data.frame()
for(i in 1:n.species){
  unique_strains <- unique(subset(logical_buscos, Species == unique_species[i])$Strain)
  if(length(unique_strains) < 2){
    if(nrow(consolidated_logical_buscos) < 1){
      consolidated_logical_buscos <- subset(logical_buscos, Species == unique_species[i])
    } else {
      consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                           subset(logical_buscos, Species == unique_species[i]))
    }
  } else {
    species_to_combine <- c(species_to_combine, unique_species[i])
  }
}

consolidated_logical_buscos <- select(consolidated_logical_buscos,
                                      -Strain)

species_to_combine # Need to find union of BUSCO genes of these species

### P. wickerhamii ###
wickerhamii_buscos <- subset(logical_buscos, Species == "P. wickerhamii")
# wickerhamii_buscos <- subset(wickerhamii_buscos, !grepl("HP",Strain))
consolidated_logical_wickerhamii <- wickerhamii_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

wickerhamii_buscos <- select(wickerhamii_buscos, !c("Species"))
wickerhamii_buscos_wide <- pivot_wider(wickerhamii_buscos, names_from = "Strain", values_from = "Status") 
wickerhamii_buscos_wide <- select(wickerhamii_buscos_wide, -BuscoID)

(wickerhamii_upset <- upset(data = wickerhamii_buscos_wide, intersect = names(wickerhamii_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 10,
      width_ratio = 0.125,
      set_sizes = FALSE))

ggsave("Outputs/BUSCO_intersections/wickerhamii_upset.png",
       wickerhamii_upset,
       width = 20,
       height = 12,
       units = "cm")

sum(consolidated_logical_wickerhamii$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_wickerhamii)
# Best possible score is 1377 

# w/o fragmented genes, using all six genomes, 145 genes are missing
# w/o fragmented genes, using just my genomes, 151 genes are missing.
# w/o fragmented genes, using just published genomes, 151 genes are missing.

# w/ fragmented genes, using all six genomes, 142 genes are missing
# w/ fragmented genes, using just my genomes, 147 genes are missing.
# w/ fragmented genes, using just published genomes, 148 genes are missing.

### P. cutis ###
cutis_buscos <- subset(logical_buscos, Species == "P. cutis")
consolidated_logical_cutis <- cutis_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

cutis_buscos <- select(cutis_buscos, !c("Species"))
cutis_buscos_wide <- pivot_wider(cutis_buscos, names_from = "Strain", values_from = "Status") 
cutis_buscos_wide <- select(cutis_buscos_wide, -BuscoID)

upset(data = cutis_buscos_wide, intersect = names(cutis_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 0,
      width_ratio = 0.125)

sum(consolidated_logical_cutis$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_cutis)

### P. bovis ###
bovis_buscos <- subset(logical_buscos, Species == "P. bovis")
consolidated_logical_bovis <- bovis_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

bovis_buscos <- select(bovis_buscos, !c("Species"))
bovis_buscos_wide <- pivot_wider(bovis_buscos, names_from = "Strain", values_from = "Status") 
bovis_buscos_wide <- select(bovis_buscos_wide, -BuscoID)

(bovis_upset <- upset(data = bovis_buscos_wide, intersect = names(bovis_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 10,
      width_ratio = 0.125,
      set_sizes = FALSE))

ggsave("Outputs/BUSCO_intersections/bovis_upset.png",
       bovis_upset,
       width = 20,
       height = 12,
       units = "cm")

sum(consolidated_logical_bovis$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_bovis)

### P. ciferrii ###
ciferrii_buscos <- subset(logical_buscos, Species == "P. ciferrii")
consolidated_logical_ciferrii <- ciferrii_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

ciferrii_buscos <- select(ciferrii_buscos, !c("Species"))
ciferrii_buscos_wide <- pivot_wider(ciferrii_buscos, names_from = "Strain", values_from = "Status") 
ciferrii_buscos_wide <- select(ciferrii_buscos_wide, -BuscoID)

upset(data = ciferrii_buscos_wide, intersect = names(ciferrii_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 0,
      width_ratio = 0.125)

sum(consolidated_logical_ciferrii$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_ciferrii)

### P. xanthoriae 
xanthoriae_buscos <- subset(logical_buscos, Species == "P. xanthoriae")
consolidated_logical_xanthoriae <- xanthoriae_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

xanthoriae_buscos <- select(xanthoriae_buscos, !c("Species"))
xanthoriae_buscos_wide <- pivot_wider(xanthoriae_buscos, names_from = "Strain", values_from = "Status") 
xanthoriae_buscos_wide <- select(xanthoriae_buscos_wide, -BuscoID)

upset(data = xanthoriae_buscos_wide, intersect = names(xanthoriae_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 0,
      width_ratio = 0.125)

sum(consolidated_logical_xanthoriae$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_xanthoriae)

### A. Protothecoides
protothecoides_buscos <- subset(logical_buscos, Species == "A. protothecoides")
consolidated_logical_protothecoides <- protothecoides_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

protothecoides_buscos <- select(protothecoides_buscos, !c("Species"))
protothecoides_buscos_wide <- pivot_wider(protothecoides_buscos, names_from = "Strain", values_from = "Status") 
protothecoides_buscos_wide <- select(protothecoides_buscos_wide, -BuscoID)

upset(data = protothecoides_buscos_wide, intersect = names(protothecoides_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 0,
      width_ratio = 0.125)

sum(consolidated_logical_protothecoides$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_protothecoides)

### C. sorokiniana ###
sorokiniana_buscos <- subset(logical_buscos, Species == "C. sorokiniana")
consolidated_logical_sorokiniana <- sorokiniana_buscos %>%
  group_by(BuscoID, Species) %>%
  summarise(Status = any(Status))

sorokiniana_buscos <- select(sorokiniana_buscos, !c("Species"))
sorokiniana_buscos_wide <- pivot_wider(sorokiniana_buscos, names_from = "Strain", values_from = "Status") 
sorokiniana_buscos_wide <- select(sorokiniana_buscos_wide, -BuscoID)

upset(data = sorokiniana_buscos_wide, intersect = names(sorokiniana_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 0,
      width_ratio = 0.125)

sum(consolidated_logical_sorokiniana$Status)
consolidated_logical_buscos <- rbind(consolidated_logical_buscos,
                                     consolidated_logical_sorokiniana)

#### Summarise individual species ####

n.buscos.in.species <- consolidated_logical_buscos %>%
  group_by(Species) %>%
  summarise(TotalPresent = sum(Status)) %>%
  mutate(TotalMissing = 1519 - TotalPresent)

n.buscos.in.species

#### Investigation of sublineages ####


consolidated_buscos_wide <- pivot_wider(consolidated_logical_buscos, names_from = "Species", values_from = "Status") 
consolidated_buscos_wide <- select(consolidated_buscos_wide, -BuscoID)

species_upset <- (upset(data = consolidated_buscos_wide, intersect = names(consolidated_buscos_wide), 
      name="Combination of Genomes", 
      min_size = 10,
      width_ratio = 0.125,
      set_sizes = FALSE))



ggsave("Outputs/BUSCO_intersections/species_upset.png",
       species_upset,
       width = 20,
       height = 22,
       units = "cm")

#### Provide list of buscos not present in wickerhamii or bovis ####
absent_wickerhamii <- consolidated_logical_wickerhamii %>%
  subset(!Status) %>%
  select(-Status)

absent_bovis <- consolidated_logical_bovis %>%
  subset(!Status) %>%
  select(-Status)

absent_combined <- full_join(absent_wickerhamii, absent_bovis, by = "BuscoID", suffix = c("w", "b"))
absent_combined[is.na(absent_combined)] <- ""
absent_combined <- mutate(absent_combined, 
       AbsentFrom = paste(Speciesw, Speciesb, sep = ";"))
absent_combined$AbsentFrom <- gsub("^;", "", absent_combined$AbsentFrom)
absent_combined$AbsentFrom <- gsub(";$", "", absent_combined$AbsentFrom)
absent_combined$AbsentFrom <- gsub(";", " and ", absent_combined$AbsentFrom)

absent_combined <- select(absent_combined, BuscoID, AbsentFrom)

buscos_function <- read.delim(data_files[1], 
                                header = T,
                                skip = 2)
buscos_function <- buscos_function[c(1,10)]
names(buscos_function)[1] <- "BuscoID"

absent_combined_function <- left_join(absent_combined, buscos_function)


subset(absent_combined_function, AbsentFrom == "P. wickerhamii")
####
