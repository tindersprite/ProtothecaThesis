##### Setup #####
#### Packages ####
require(ggplot2)
require(tidyverse)

#### Misc ####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/")
# "/Users/u1984259/Documents/Main_Project/Data/Plate_Reader/Plate_Reader_Data/Mapped_Data/"

plate.reader.data  <- "Plate_Reader_Data/"
plate.reader.plots <- "Plate_Reader_Plots/Normalised_QC_Plots/"

# created with $echo {0..100}_h_{0,30}_min | tr " " "." | tr _ " " > time_point_levels.txt
time.point.levels <- read.delim('~/Documents/Main_Project/Scripts/update_platereader_data/time_point_levels.txt', header = F, sep = ".")
time.point.levels <- unlist(time.point.levels)
names(time.point.levels) <- NULL
time.point.levels <- time.point.levels[-length(time.point.levels)]

strain.numbers <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")


#####

mapped.data <- list.files("Plate_Reader_Data/Mapped_Data/")

for(a in 1:length(mapped.data)){
  ### Skip analysis for files that already exist
  base.filename <- mapped.data[a]
  base.filename <- substring(base.filename, 8,nchar(base.filename)-4)
  
  qc.plot.name     <- paste(plate.reader.plots, 
                            base.filename, 
                            ".png", sep = "")
  normalised.data.name <- paste(plate.reader.data, 
                                "Normalised_Data/normalised_",
                                base.filename,
                                ".csv", sep = "")
  
  if(file.exists(qc.plot.name) && file.exists(normalised.data.name)) next
  
  ### Force processing of files that were updated by previous scripts
  
  ### Load mapped data
  total.data <- read.csv(paste("Plate_Reader_Data/Mapped_Data/",
                               mapped.data[a],
                               sep = ""))
  total.data$Cond <- NA
  cond.columns <- names(total.data)
  columns.to.remove <- c("StrainID", 'x', 'y', 'OD', 'Replicate', 'xy',  # Columns that vary per well
                         'InoculatedDensity', 'InocDensWavelength', 'PrecultureMedia', 'PrecultureTemperature', # Columns that don't exist for NEG
                         'PrecultureShaking', 'StockMedia', 'StockTemperature', # Columns that don't exist for NEG
                         'Notes', 'MapVersion', 'Lamp', 'Ncol', 'Nrow', 'PlateSize', 'PlateNumber') # Columns that aren't worth including
  cond.columns <- cond.columns[!(cond.columns %in% columns.to.remove)]
  for(i in 1:nrow(total.data)){
    total.data$Cond[i] <- paste(total.data[cond.columns][i,], collapse = "")
  }
  
  ### Identify and average blank wells
  blank.wells <- subset(total.data, StrainID == "NEG")
  blanks <- c()
  
  for (condition in unique(blank.wells$Cond)){
    temp <- blank.wells[which(blank.wells$Cond == condition),]
    blank.OD <- mean(temp$OD)
    new.blank <- data.frame(Cond = condition, OD = blank.OD)
    blanks <- rbind(blanks, new.blank)
  }
  
  ### Normalise experimental wells
  normalised.data <- subset(total.data, toupper(StrainID) != "NEG")
  normalised.data <- mutate(normalised.data, NormalisedOD = NA)
  
  for(well in 1:nrow(normalised.data)){
    row <- which(blanks$Cond == normalised.data$Cond[well])
    OD.norm <- normalised.data$OD[well] - blanks$OD[row]
    normalised.data$NormalisedOD[well] <- OD.norm
  }
  
  if(any(is.na(normalised.data$NormalisedOD))){
    print(paste("NAs are present after normalising OD. Normalising OD has failed."))
  }
  
  #### Add species ####
  normalised.data <- mutate(normalised.data, SpeciesName = NA)
  
  for(i in 1:nrow(strain.numbers)){
    strain <- which(strain.numbers$H_number[i] == normalised.data$StrainID)
    normalised.data$SpeciesName[strain] <- strain.numbers$Species[i]
  }
 
  #### Plot normalised data ####
  normalised.data$TimePoint <- factor(normalised.data$TimePoint, levels = time.point.levels)
  normalised.data$StrainID  <- factor(normalised.data$StrainID, 
                                      levels = strain.numbers$H_number[-length(strain.numbers$H_number)])
  
  number.of.species <- length(unique(normalised.data$SpeciesName))
  normalised.data <- mutate(normalised.data, xy = paste(x,y,sep = ""))
  
  quality.plot <- ggplot(data = normalised.data, 
                         aes(y = NormalisedOD, 
                             x = (as.numeric(TimePoint)/2+1))) +
    geom_line(aes(group = xy, col = StrainID)) +
    scale_x_continuous(name = "Time (hours)",
                       breaks = floor(seq(from = 0, 
                                    to = max(as.numeric(normalised.data$TimePoint)/2+1), 
                                    length.out = 8))) +
    scale_y_continuous(name = paste("Normalised OD",
                                    unique(normalised.data$Wavelength),
                                    sep = "")) +
    facet_wrap(~ SpeciesName,
               nrow = ceiling(sqrt(number.of.species)),
               ncol = ceiling(sqrt(number.of.species)),
               scales = "free") 
  
  ggsave(filename = qc.plot.name,
         plot = quality.plot,
         device = "png",
         width = 18,
         height = 10,
         units = "cm")
  
  #### Export normalised data ####
  normalised.data <- normalised.data[!names(normalised.data)=="Cond"]
  
  write.csv(normalised.data,
            file = normalised.data.name,
            row.names = FALSE)
}



