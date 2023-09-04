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

#####

# Each section should be unique
# Collapsable titles for the experiments and gc
# use a function to process the data
# could I set up standard plots to fill in? e.g. plot.temp() for growth curve data




##### Inputs #####
data.file <- file.choose()
#####

##### Import the data #####
normalised.data <- read.csv(data.file)
##### Set up function to process data
# 
##### Tidy the data #####
### Set factors

normalised.data$StrainID <- factor(normalised.data$StrainID, 
                                   levels = strain.levels)
normalised.data$TimePoint <- factor(normalised.data$TimePoint,
                                    levels = time.point.levels)

#####

##### plot pH data #####
ggplot(data = test.wells, aes(x = Time)) +
  geom_line(aes(y = OD.norm, group = xy, col = Strain)) +
  facet_grid(Species ~ pH, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(data = sub.wells, aes(x = Time)) +
  geom_line(aes(y = OD.norm, group = xy, col = xy)) +
  facet_grid(Strain ~ pH, scales = "fixed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sub.wells <- subset(test.wells, pH == 2)

### Neater times

ggplot(data = test.wells, aes(x = (as.numeric(Time)-1)/2)) +
  geom_line(aes(y = OD.norm, group = xy, col = Strain)) +
  facet_grid(Species ~ pH, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(test.wells$Time)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", test.wells$Wavelength[1], sep = "")) +
  #scale_y_continuous(name = "Normalised OD600") +
  ggtitle(paste("pH growth curve", test.wells$Date[1]))


ggplot(data = test.wells, aes(x = (as.numeric(Time)-1)/2)) +
  geom_smooth(aes(y = OD.norm, group = Strain, col = Strain)) +
  facet_grid(Species ~ H2O2Conc, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(test.wells$Time)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", test.wells$Wavelength[1], sep = "")) +
  #scale_y_continuous(name = "Normalised OD600") +
  ggtitle(paste("pH growth curve", test.wells$Date[1]))

##### plot ROS data #####
ggplot(data = test.wells, aes(x = Time)) +
  geom_line(aes(y = OD.norm, group = xy, col = Strain)) +
  facet_grid(Species ~ h2o2.conc, scales = "free")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggplot(data = sub.wells, aes(x = Time)) +
  geom_line(aes(y = OD.norm, group = xy, col = xy)) +
  facet_grid(Strain ~ h2o2.conc, scales = "fixed")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sub.wells <- subset(test.wells, H2O2Conc == 0)

### Neater times


ggplot(data = test.wells, aes(x = (as.numeric(Time)-1)/2)) +
  geom_line(aes(y = OD.norm, group = xy, col = Strain)) +
  facet_grid(Species ~ pH, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(test.wells$Time)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", test.wells$Wavelength[1], sep = "")) +
  #scale_y_continuous(name = "Normalised OD600") +
  ggtitle(paste("ROS growth curve", test.wells$Date[1]))


ggplot(data = test.wells, aes(x = (as.numeric(Time)-1)/2)) +
  geom_smooth(aes(y = OD.norm, group = Strain, col = Strain)) +
  facet_grid(Species ~ H2O2Conc, scales = "free")+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(test.wells$Time)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", test.wells$Wavelength[1], sep = "")) +
  #scale_y_continuous(name = "Normalised OD600") +
  ggtitle(paste("ROS growth curve", test.wells$Date[1]))

##### plot temperature data #####

ggplot(data = normalised.data, aes(x = TimePoint)) +
  geom_line(aes(y = NormalisedOD, group = xy, col = StrainID)) +
  facet_wrap(SpeciesName ~ ., scales = "fixed") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Neater times

ggplot(data = normalised.data, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_line(aes(y = NormalisedOD, group = xy, col = StrainID)) +
  facet_wrap(SpeciesName ~ ., scales = "fixed") +
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(normalised.data$TimePoint)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", normalised.data$Wavelength[1], sep = "")) +
  ggtitle(paste(normalised.data$CultureTemperature[1], 
                "degree growth curve", normalised.data$Date[1]))

ggplot(data = normalised.data, aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_line(aes(y = NormalisedOD, group = xy, col = StrainID)) +
  facet_wrap(SpeciesName ~ ., scales = "free") +
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(normalised.data$TimePoint)), 
                                  by = 4)) +
  scale_y_continuous(name = paste("Normalised OD", normalised.data$Wavelength[1], sep = "")) +
  ggtitle(paste(normalised.data$CultureTemperature[1], 
                "degree growth curve", normalised.data$Date[1]))


#####



