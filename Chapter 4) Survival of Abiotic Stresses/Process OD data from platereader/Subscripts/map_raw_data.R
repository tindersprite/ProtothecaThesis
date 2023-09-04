##### Setup #####
#### Packages ####
require(openxlsx)
require(tidyverse)
require(ggplot2)

#### Misc ####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data")

strain.names <- read_csv(file = "strain_numbers.csv")
strain.names <- c("NEG", strain.names$H_number)
strain.names <- strain.names[-length(strain.names)]



# plate.reader.dir   <- "plate_reader_link/"
plate.reader.data  <- "plate_reader_link/Plate_Reader_Data/"
plate.reader.plots <- "plate_reader_link/Plate_Reader_Plots/Mapped_QC_Plots/"

#####

### Import mapping pairs

mapping.pairs <- read.csv("tmp/mapping_pairs.csv")

raw.data.files     <- mapping.pairs$raw_data_file
mapping.data.files <- mapping.pairs$mapping_data_file


# "plate_reader_link/Plate_Reader_Data/Mapping_Data/"
for(a in 1:length(raw.data.files)){
  
  
  #### Skip data that has already been processed
  base.filename <- raw.data.files[a]
  base.filename <- strsplit(base.filename, split = "[.]")[[1]][1]
  
  qc.plot.name     <- paste(plate.reader.plots, 
                            base.filename, 
                            ".png", sep = "")
  mapped.data.name <- paste(plate.reader.data, 
                            "Mapped_Data/mapped_",
                            base.filename,
                            ".csv", sep = "")
  if(file.exists(mapped.data.name) & file.exists(qc.plot.name)){
    next
    }
  
  read.data    <- paste(plate.reader.data,"Raw_Data/",
                        raw.data.files[a], sep = "")
  mapping.data <- paste(plate.reader.data,"Mapping_Data/",
                        mapping.data.files[a], sep = "")
  
  ### Find useful intervals
  sheets <- getSheetNames(read.data)
  intervals <- sheets[substring(sheets, 1,1) == "I"]
  
  if("Interval 1 (0 h  )" %in% intervals){
    spare <- which(intervals == "Interval 1 (0 h  )")
    intervals <- intervals[-spare]
  }
  
  ### There is occasionally variation in the intervals.
  ### Some files have x h 60 minutes, instead of x+1 h
  # e.g. acidCheck_OD600_Plate1_2020-11-13
  ### Some will vary sightly around the half hour time points
  ### Remove this variation
  split.intervals <- strsplit(intervals, split = " ")
  length.of.split <- length(split.intervals)
  
  minutes <- hours <- vector(length = length.of.split)
  for(i in 1:length.of.split){
  hours[i]   <- split.intervals[[i]][3]
  minutes[i] <- split.intervals[[i]][5]
  }
  
  hours <- substring(hours, 2, nchar(hours))
  minutes[is.na(minutes)] <- "0"
  
  hours <- as.numeric(hours)
  minutes <- as.numeric(minutes)
  
  minutes[minutes <= 10] <- 0
  minutes[minutes >= 20 & minutes <= 40] <- 30
  minutes[minutes >= 50] <- 60
  
  rounding.errors <- which(minutes == 60)
  minutes[rounding.errors] <- 0
  hours[rounding.errors] <- hours[rounding.errors] + 1
  # Is it worth recording which values are modified?
  
  time.points <- paste(hours, " h ", minutes, " min", sep = "")
  
  
  ### Check what wavelength was was used to measure OD ###
  OD.wavelength <- raw.data.files[a]
  OD.wavelength <- strsplit(OD.wavelength, split = "_")[[1]][2]
  OD.wavelength <- substring(OD.wavelength, 3,5)
  
  info.sheet <- read.xlsx(read.data, 
                          sheet = "Protocol Information")
  excitation <- info.sheet[which(info.sheet[,1] == "Excitation:"), 2]
  
  if(excitation != OD.wavelength){
    print(paste("An error has occurred. The OD wavelength of a file name does not match the excitation ",
                "wavelength in the corresponding raw data file. The file in question is: ",
                raw.data.files[a],
                sep = "", collapse = ""))
  }
  
  
  ### Import the raw data and position data ###
  ## Read data from the intervals
  treatment.data <- read.xlsx(mapping.data, 
                              sheet = "Condition Mapping")
  
  ## The data can be in different places, depending on the options chosen in exporting
  ## Find the cell that contains "Raw Data (ODxxx)"
  startrow <- 0
  label.raw.data <- "Placeholder"
  
  while(substring(label.raw.data,1,8) != "Raw Data"){
    startrow <- startrow +1
    label.raw.data <- read.xlsx(xlsxFile = read.data,
                                sheet = intervals[1],
                                rows = startrow,
                                cols = 2, 
                                colNames = F)
    # There will be warnings when trying to load cells with nothing in them
    # The warnings will be "no data found on worksheet"
    if(is.null(label.raw.data)){
      label.raw.data <- "Placeholder"
    } else {
      label.raw.data <- label.raw.data[1,1]
    }
  }
  
  ### Match measurements to the relevant mapping data, and save time points
  total.data <- c()
  treatment.data <- arrange(treatment.data, x, y) # Allows row position to be used to assign values to OD measurements
  for(i in 1:length(intervals)){
    sheet <- intervals[i]
    ### This should work for larger or smaller plates
    interval.data <- read.xlsx(xlsxFile = read.data,
                               sheet = sheet,
                               startRow = startrow+1,
                               colNames = T,
                               rowNames = T)
    
    interval.data <- unlist(interval.data)
    names(interval.data) <- NULL
    
    time <- time.points[i]
    
    treatment.data$TimePoint  <- time 
    treatment.data$OD         <- interval.data
    treatment.data$Wavelength <- excitation
    total.data <- rbind(total.data, treatment.data)
  }
  
  ### Plot quality checking graphs
  total.data <- mutate(total.data, xy = paste(x, y, sep = '')) # to allow proper grouping
  
  total.data$StrainID  <- factor(total.data$StrainID, levels = strain.names)
  ### If files are combined manually, typos can introduce duplicate time points, which will cause errors here
  ### Update: a script has now been written to automatically combine fragmented data,
  ### which can result in duplicated time points. Possible need for FIXME
  total.data$TimePoint <- factor(total.data$TimePoint, levels = unique(time.points))
  number.of.strains    <- length(unique(treatment.data$StrainID))
  
  ### Automatically combining different files seems to classify
  ### the OD values as character data.
  total.data$OD <- as.numeric(total.data$OD)
  
  quality.plot <- ggplot(data = total.data, 
                         aes(y = OD, 
                             x = (as.numeric(TimePoint)/2+1))) +
    geom_line(aes(group = xy, col = StrainID)) +
    scale_x_continuous(name = "Time (hours)",
                       breaks = floor(seq(from = 0, 
                                          to = max(as.numeric(total.data$TimePoint)/2+1), 
                                          length.out = 4))) +
    facet_wrap(~ StrainID,
               nrow = ceiling(sqrt(number.of.strains)),
               ncol = ceiling(sqrt(number.of.strains)),
               scales = "free") 
  
  ggsave(filename = qc.plot.name,
         plot = quality.plot,
         device = "png",
         width = 18,
         height = 10,
         units = "cm")
  
  ### Export mapped data
  
  write.csv(x = total.data, 
            file = mapped.data.name,
            row.names = FALSE)
  
  #### Remove carryover variables
  intervals <- time.points <- OD.wavelength <- excitation <- base.filename <- qc.plot.name <- mapped.data.name <- total.data <- NA
}
#####

#### Custom QC Plots ####
