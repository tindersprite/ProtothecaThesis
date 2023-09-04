##### Setup #####
#### Packages ####
require(openxlsx)
require(tidyverse)
require(ggplot2)
require(metafolio)

#### Misc ####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data")

strains <- read_csv(file = "strain_numbers.csv")
strain.names <- c("NEG", strains$H_number)
strain.names <- strain.names[-length(strain.names)]

setwd("plate_reader_link/Plate_Reader_Data/")

# plate.reader.dir   <- "plate_reader_link/"
plate.reader.data  <- "plate_reader_link/Plate_Reader_Data/"
plate.reader.plots <- "plate_reader_link/Plate_Reader_Plots/Mapped_QC_Plots/"

#####

#### Import data ####
data <- read.xlsx(xlsxFile = "Mapping_Data/GC30_1.0.1_Plate1_2021-09-03.xlsx")

#### Assign species names ####

data$SpeciesName <-NA

for(i in 1:nrow(strains)){
  hx <- strains$H_number[i]
  data$SpeciesName[data$StrainID == hx] <- strains$Species[i]
}

#### Set factors ####
data$x <- factor(data$x, levels = 1:12)
data$y <- factor(data$y, levels = LETTERS[26:1])
data$StrainID <- factor(data$StrainID,
                        levels = strain.names)
data$SpeciesName <- factor(data$SpeciesName, 
                           levels = unique(strains$Species))

#### Create Plots ####
map <- ggplot(data = data, aes(x = x, y = y)) + 
  ggtitle(paste(data$Experiment[1], 
                data$Date[1], 
                sep = "_")) +
  theme(plot.title = element_text(hjust = 0.5))


rep.map <- map + # Optional plot, to show the position of replicates
  geom_point(aes(shape = Replicate), size = 8)

(test.map <- map +
  geom_point(aes(colour = SpeciesName), size = 16)+
  scale_shape_manual(values = c(15:19)) +
  scale_colour_manual(values = strainpal) + 
  ggtitle("Map showing the position of different species in \na 48 well plate for the 30°C growth curve run\non the 3rd of September 2021"))

if(length(data$Strain) > 8){
  factor.of.interest <- "SpeciesName"
  strainpal <- gg_color_hue(length(unique(data[[factor.of.interest]]))-1)
  strainpal <- c(strainpal, '#000000')
  
  (strain.map <- map +
    geom_text(aes_string(colour = factor.of.interest, label = "StrainID"))+
      scale_colour_manual(values = strainpal))
}

