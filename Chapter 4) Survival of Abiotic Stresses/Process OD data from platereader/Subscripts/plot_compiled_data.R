##### Setup #####
#### Packages ####
require(ggplot2)
require(tidyverse)

#### Misc ####
plate.reader.link <- '~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/'
plate.reader.data <- 'Plate_Reader_Data/Normalised_Data/'
plate.reader.plots<- 'Plate_Reader_Plots/Normalised_Plots/'
setwd(plate.reader.link)


# created with $echo {0..100}_h_{0,30}_min | tr " " "." | tr _ " " > time_point_levels.txt
time.point.levels <- read.delim('~/Documents/Main_Project/Scripts/update_platereader_data/time_point_levels.txt', header = F, sep = ".")
time.point.levels <- unlist(time.point.levels)
names(time.point.levels) <- NULL
time.point.levels <- time.point.levels[-length(time.point.levels)]

strains <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")
strain.levels <- strains$H_number
original.species.level <- NA

species.levels.pathology <- c("Prototheca bovis","Prototheca wickerhamii","Prototheca blaschkeae",
                              "Prototheca cutis","Prototheca miyajii","Prototheca ciferrii",
                              "Prototheca pringsheimii","Prototheca zopfii",
                              "Prototheca cerasi","Prototheca cookei",
                              "Prototheca tumulicola","Prototheca stagnora",
                              "Prototheca moriformis","Prototheca paracutis",
                              "Prototheca xanthoriae","Auxenochlorella protothecoides",
                              "Helicosporidium spp.")
species.levels.descent   <- c("Prototheca bovis","Prototheca pringsheimii","Prototheca zopfii",
                              "Prototheca cerasi","Prototheca ciferrii","Prototheca cookei",
                              "Prototheca blaschkeae","Prototheca tumulicola","Prototheca stagnora",
                              "Prototheca moriformis","Prototheca cutis","Prototheca paracutis",
                              "Prototheca miyajii","Prototheca xanthoriae","Auxenochlorella protothecoides",
                              "Prototheca wickerhamii","Helicosporidium spp.")
species.levels.alphabet <- c('Auxenochlorella protothecoides','Prototheca blaschkeae','Prototheca bovis',
                             'Prototheca ciferrii','Prototheca cookei','Prototheca cutis',
                             'Prototheca miyajii','Prototheca paracutis','Prototheca tumulicola',
                             'Prototheca wickerhamii','Prototheca xanthoriae')

species.levels <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")

#####

##### Check All Strains 

###### Temperature Growth Curves ######
##### Data Processing #####
#### Import data ####
GC25.1 <- read.csv(paste(plate.reader.data,"normalised_GC25_OD750_Plate1_2021-05-28.csv", sep = ""))
GC30.1 <- read.csv(paste(plate.reader.data,"normalised_GC30_OD750_Plate1_2021-06-18.csv", sep = ""))
GC37.1 <- read.csv(paste(plate.reader.data,"normalised_GC37_OD750_Plate1_2021-07-16.csv", sep = ""))
GC42.1 <- read.csv(paste(plate.reader.data,"normalised_GC42_OD750_Plate1_2021-08-06.csv", sep = ""))
GC25.2 <- read.csv(paste(plate.reader.data,"normalised_GC25_OD750_Plate1_2021-08-13.csv", sep = ""))
GC30.2 <- read.csv(paste(plate.reader.data,"normalised_GC30_OD750_Plate1_2021-09-03.csv", sep = ""))
GC37.2 <- read.csv(paste(plate.reader.data,"normalised_GC37_OD750_Plate1_2021-09-07.csv", sep = ""))
GC42.2 <- read.csv(paste(plate.reader.data,"normalised_GC42_OD750_Plate1_2021-09-12.csv", sep = ""))
GC30.3 <- read.csv(paste(plate.reader.data,"normalised_GC30_OD750_Plate1_2022-02-06.csv", sep = ""))
GC37.3 <- read.csv(paste(plate.reader.data,"normalised_GC37_OD750_Plate1_2022-02-11.csv", sep = ""))
GC30.4 <- read.csv(paste(plate.reader.data,"normalised_GC30_OD750_Plate1_2022-02-25.csv", sep = ""))
GC37.4 <- read.csv(paste(plate.reader.data,"normalised_GC37_OD750_Plate1_2022-03-04.csv", sep = ""))

#### Compile data ####

GC_batch1 <- rbind(GC25.1, GC30.1, GC37.1, GC42.1)
GC_batch2 <- rbind(GC25.2, GC30.2, GC37.2, GC42.2)
GC_batch3 <- rbind(GC30.3, GC37.3)
GC_batch4 <- rbind(GC30.4, GC37.4)

GC_batch1 <- mutate(GC_batch1, Batch = 1)
GC_batch2 <- mutate(GC_batch2, Batch = 2)
GC_batch3 <- mutate(GC_batch3, Batch = 3)
GC_batch4 <- mutate(GC_batch4, Batch = 4)

GC <- rbind(GC_batch1, GC_batch2, GC_batch3, GC_batch4)
GC$Batch <- factor(GC$Batch,
                   levels = sort(unique(as.numeric(GC$Batch))))

#### Process Data ####
GC$Infection <- NA
GC$Host <- NA
for(i in 1:nrow(strains)){
  strain.ID <- strains$H_number[i]
  infect <- strains$From_infection[i]
  host <- strains$Host[i]
  
  rows <- which(GC$StrainID == strain.ID)
  GC$Infection[rows] <- infect
  GC$Host[rows]      <- host
}

GC$TimePoint <- factor(GC$TimePoint, levels = time.point.levels)
GC$StrainID <- factor(GC$StrainID, levels = strain.levels)

write.csv(GC,
          file = "Plate_Reader_Data/Compiled_Data/compiled_temperature_OD750.csv",
          row.names = FALSE)

##### Plotting #####
#### Check replicates #### Not yet done ####
#### Plot Total Data ####
GC$SpeciesName <- factor(GC$SpeciesName, 
                         levels = species.levels$Descent) #descent, alphabet, or pathology

(GC_byspecies <- ggplot(data = GC, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = StrainID)) +
  facet_grid(SpeciesName ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))) +
  ggtitle("Data from the first replicate of Prototheca grown under temperature stress")

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/temperature_OD750_bystrain_",
                        Sys.Date(),".png", sep = ""),
      plot = GC_byspecies,
      device = "png",
      width = 20,
      height = 20,
      units = "cm")

# Colour by infection
(GC_byinfection <- ggplot(data = GC, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
  facet_grid(SpeciesName ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))+
  ggtitle("Data from the first replicate of Prototheca grown under temperature stress",
          subtitle = "Coloured by whether the isolate was derived from an infection"))

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/temperature_OD750_byinfection_",
                        Sys.Date(),".png", sep = ""),
       plot = GC_byinfection,
       device = "png",
       width = 20,
       height = 20,
       units = "cm")

# Colour by host
(GC_byhost <- ggplot(data = GC, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Host)) +
  facet_grid(SpeciesName ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0)) + 
  ggtitle("Data from the first replicate of Prototheca grown under temperature stress",
          subtitle = "Coloured by which host the isolate is associated with"))

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/temperature_OD750_byhost_",
                        Sys.Date(),".png", sep = ""),
       plot = GC_byhost,
       device = "png",
       width = 20,
       height = 20,
       units = "cm")

#### Facet on strain instead of temperature and species ####
# Including all the strains is too much to keep track of
GC <- mutate(GC, StrTempDate = paste(StrainID, CultureTemperature, Date, sep = "_"))
GC$CultureTemperature <- factor(GC$CultureTemperature,
                                levels = unique(sort(GC$CultureTemperature)))
single.species <- subset(GC, SpeciesName == "Prototheca blaschkeae")

ggplot(data = single.species, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrTempDate, col = CultureTemperature)) +
  facet_wrap(StrainID ~ ., scales = "fixed", nrow = ceiling(sqrt(length(unique(single.species$StrainID)))))+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(data = single.species, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrTempDate, col = Infection)) +
  facet_wrap(CultureTemperature ~ ., scales = "fixed", nrow = 2)+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(data = single.species, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = StrainID)) +
  facet_grid(StrainID ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0)) 

#### Plot Individual Strains ####

single.strain <- subset(GC, StrainID %in% c("HP51"))
single.strain <- mutate(single.strain, StrDate = paste(StrainID, Date, collapse = "_"))

### Plot each experiment for a single strain
ggplot(data = single.strain, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = xy, col = Date)) +
  facet_grid(CultureTemperature  ~ Batch, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0)) + 
  ggtitle(label = paste("Comparison of the growth of individual wells of ",
                        single.strain$StrainID[1], 
                        " (", single.strain$SpeciesName[1], ")",
  #                      "(", single.strain$Host[1], ")",
  #                      "(", single.strain$Infection[1], ")",
                        "\n",
                        "across different experimental replicates", 
                        sep = ""))

### Plot strains from different batches on the same graphs. 
### Capable of plotting multiple strains
single.column <- ggplot(data = single.strain, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = xy, col = Date)) +
  facet_grid(CultureTemperature ~ StrainID , scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

if(length(unique(single.strain$StrainID))==1){
  single.column + 
    ggtitle(label = paste("Comparison of the growth of individual wells of ",
                          single.strain$StrainID[1], 
                          " (", single.strain$SpeciesName[1], ")",
                          #                      "(", single.strain$Host[1], ")",
                          #                      "(", single.strain$Infection[1], ")",
                          "\n",
                          "across different experimental replicates", 
                          sep = ""))
} else {
  single.column + 
    ggtitle(label = paste("Comparison of the growth of individual wells of strains",
                          "\n",
                          "across different experimental replicates", 
                          sep = ""))
}
#### Exploratory Plots ####

## I observed that in GC37.2, some bovis strains appeared to be reaching
## ODs of about 3. I remembered this was higher than expected, and 
## wondered if infectious strains might grow better at higher temperature
# subset <- subset(GC, SpeciesName == "Prototheca bovis")
ggplot(data = subset, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = StrainID)) +
  facet_grid(Infection ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

ggplot(data = subset, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
  facet_grid(Batch ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

## It seems there was no clear difference between infectious strains 
## growing at higher temperatures, though some bovis does seem to grow
## better at 37°C. Is this a batch effect?

ggplot(data = subset, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Batch)) +
  facet_grid(Infection ~ CultureTemperature, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0))

## This difference in growth might be a batch effect. Batch 1 bovis 
## seems to grow slightly better at 37°

#####

###### pH plots ######
#### Import data ####
GCacid1.1 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-03-11.csv", sep = ""))
GCacid1.2 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-03-18.csv", sep = ""))
GCacid1.3 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-04-01.csv", sep = ""))
GCacid1.4 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-04-08.csv", sep = ""))
GCacid1.5 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-04-15.csv", sep = ""))

GCacid2.1 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-06-11.csv", sep = ""))
GCacid2.2 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-06-17.csv", sep = ""))
GCacid2.3 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-07-01.csv", sep = ""))
GCacid2.4 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-07-08.csv", sep = ""))
GCacid2.5 <- read.csv(paste(plate.reader.data,"normalised_ACID_OD750_Plate1_2022-07-15.csv", sep = ""))

#### Compile data ####

GC_batch1 <- rbind(GCacid1.1, GCacid1.2, GCacid1.3, GCacid1.4, GCacid1.5)
GC_batch1 <- mutate(GC_batch1, Batch = 1)

GC_batch2 <- rbind(GCacid2.1, GCacid2.2, GCacid2.3, GCacid2.4, GCacid2.5)
GC_batch2 <- mutate(GC_batch2, Batch = 2)

GC <- rbind(GC_batch1, GC_batch2)
GC$Batch <- factor(GC$Batch,
                   levels = sort(unique(as.numeric(GC$Batch))))

#### Process Data ####
GC$Infection <- NA
GC$Host <- NA
for(i in 1:nrow(strains)){
  strain.ID <- strains$H_number[i]
  infect <- strains$From_infection[i]
  host <- strains$Host[i]
  
  rows <- which(GC$StrainID == strain.ID)
  GC$Infection[rows] <- infect
  GC$Host[rows]      <- host
}

GC$TimePoint <- factor(GC$TimePoint, levels = time.point.levels)
GC$StrainID <- factor(GC$StrainID, levels = strain.levels)

write.csv(GC,
          file = "Plate_Reader_Data/Compiled_Data/compiled_acid-25deg_OD750.csv",
          row.names = FALSE)

##### Plotting #####
#### Plot Total Data ####
GC$SpeciesName <- factor(GC$SpeciesName, 
                         levels = species.levels$Descent) 
#descent, alphabet, or pathology

## Show variation by strain, across the two experiments
GC <- mutate(GC, Ident = paste(pH, Batch, Replicate, sep = "_"))
for(i in unique(GC$Batch)){
  (GC_bystrainbatch <- ggplot(data = subset(GC, Batch == i), 
                              aes(x = (as.numeric(TimePoint)-1)/2)) +
     geom_line(aes(y = NormalisedOD, group = Ident, col = as.character(pH))) +
     facet_wrap(StrainID ~ ., 
                nrow = ceiling(sqrt(length(unique(GC$StrainID)))),
                scales = "free")+
     scale_x_continuous(name = "Time (hours)",
                        breaks = seq(from = 0, 
                                     to = max(as.numeric(GC$TimePoint)), 
                                     by = 12)) +
     scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
     coord_cartesian(xlim = c(0, 66),
                     ylim = NULL) +
     theme(strip.text.y = element_text(angle = 0)) +
     ggtitle(label = paste("Data from a replicate of Prototheca grown ",
                           "under acid stress",
                           "\ndata from batch ", i, sep = "")))
  
  ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/acid_OD750_varCheck_batch",
                          i,"_",
                          Sys.Date(),".png", sep = ""),
         plot = GC_bystrainbatch,
         device = "png",
         width = 25,
         height = 20,
         units = "cm")
}

(GC_bystrainbatch <- ggplot(data = GC, 
                            aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_line(aes(y = NormalisedOD, group = Ident, col = as.character(pH))) +
    facet_wrap(StrainID ~ ., 
               nrow = ceiling(sqrt(length(unique(GC$StrainID)))),
               scales = "free")+
    scale_x_continuous(name = "Time (hours)",
                       breaks = seq(from = 0, 
                                    to = max(as.numeric(GC$TimePoint)), 
                                    by = 12)) +
    scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
    coord_cartesian(xlim = c(0, 66),
                    ylim = NULL) +
    theme(strip.text.y = element_text(angle = 0)) +
    ggtitle(label = paste("Data from all replicate of Prototheca grown ",
                          "under acid stress", sep = "")))

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/acid_OD750_varCheck_",
                        Sys.Date(),".png", sep = ""),
       plot = GC_bystrainbatch,
       device = "png",
       width = 25,
       height = 20,
       units = "cm")

## Show variation by species
(GC_byspecies <- ggplot(data = GC, aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_smooth(aes(y = NormalisedOD, group = StrainID, col = StrainID)) +
    facet_grid(SpeciesName ~ pH, scales = "free")+
    scale_x_continuous(name = "Time (hours)",
                       breaks = seq(from = 0, 
                                    to = max(as.numeric(GC$TimePoint)), 
                                    by = 12)) +
    scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
    coord_cartesian(xlim = c(0, 66),
                    ylim = NULL) +
    theme(strip.text.y = element_text(angle = 0))) +
  ggtitle("Data from the first replicate of Prototheca grown under acid stress")

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/acid_OD750_bystrain_",
                        Sys.Date(),".png", sep = ""),
       plot = GC_byspecies,
       device = "png",
       width = 20,
       height = 20,
       units = "cm")

### Colour by infection

(GC_byinfection <- ggplot(data = GC, aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
    facet_grid(SpeciesName ~ pH, scales = "free")+
    scale_x_continuous(name = "Time (hours)",
                       breaks = seq(from = 0, 
                                    to = max(as.numeric(GC$TimePoint)), 
                                    by = 12)) +
    scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
    coord_cartesian(xlim = c(0, 66),
                    ylim = NULL) +
    theme(strip.text.y = element_text(angle = 0))) +
  ggtitle("Data from the first replicate of Prototheca grown under acid stress")

ggsave(filename = paste("Plate_Reader_Plots/Compiled_Plots/acid_OD750_byinfection_",
                        Sys.Date(),".png", sep = ""),
       plot = GC_byinfection,
       device = "png",
       width = 20,
       height = 20,
       units = "cm")

#### Exploratory plots ####

single.species <- subset(GC, SpeciesName == "Prototheca wickerhamii")

ggplot(data = single.species, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = xy, col = Date)) +
  facet_grid(pH ~ StrainID, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 66),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0)) + 
  ggtitle(label = paste("Comparison of the growth of individual wells of ",
                        single.species$SpeciesName[1],
                        #                      "(", single.strain$Host[1], ")",
                        #                      "(", single.strain$Infection[1], ")",
                        "\n",
                        "across different experimental replicates", 
                        sep = ""))

#####






