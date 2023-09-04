##### Setup #####
#### Packages ####
require(openxlsx)
require(lubridate)

#### Misc ####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/Plate_Reader_Data")

#####

##### Code #####
#### Import fragmented data files ####
f <- Sys.glob(file.path(getwd(),"Fragmented_Raw_Data/*"))

f <- strsplit(f, split = "/")
files <- c()

for(i in 1:length(f)){
  file.name.length <- length(f[[i]])
  files[i] <- f[[i]][file.name.length]
}



#### Identify fragmented experiments ####
experiments <- c()

for(i in 1:length(files)){
  s <- strsplit(files[i], split = "_")[[1]]
  s <- s[1:(length(s)-1)]
  experiments[i] <- paste(s, collapse = "_")
}

experiments <- unique(experiments)

#### Loop through fragmented experiments ####

for(i in 1:length(experiments)){
  #### Identify which files need to be used ####
  experiment <- experiments[i]
  matching.files <- files[grep(experiment, files)]
  
  
  #### Skip experiments that have already been combined ####
  output.file.name <- paste("Raw_Data/",
                            experiment,
                            ".xlsx", 
                            sep = "")
  
  if(file.exists(output.file.name)) next
  

  #### Identify the relevant sheet names in each file ####
  sheet.names <- list()
  
  for(j in 1:length(matching.files)){
    ### Load all sheet names
    sheet.names.from.file <- getSheetNames(file.path("Fragmented_Raw_Data", 
                                                     matching.files[j]))
    
    ### Remove redundant sheets 
    redundant.sheets <- c("Interval 1 (0 h  )", 
                          'Interval 1 (0 min  )',
                          "Protocol Information")
    sheet.names.from.file <- sheet.names.from.file[!(sheet.names.from.file %in% 
                                                       redundant.sheets)]
    
    ### Export sheet names in a single object
    # It is necessary to convert the sheet.names.from.file object into a list
    # to pass it to the sheet.names object 
    sheet.names[paste("Session", j, sep = "_")] <- list(sheet.names.from.file)
  }
  
  #### Set new sheet names ####
  ### Work out when the fragments started, relative to each other
  start.times <- c()
  start.dates <- c()
  for(j in 1:length(matching.files)){
    data <- read.xlsx(file.path("Fragmented_Raw_Data", 
                                matching.files[j]))
    time <- unlist(data[grep("Time", data)])
    time <- time[grep("Time", time)]
    time <- strsplit(time, split = " ")[[1]][2]
    start.times[j] <- time
    
    date <- unlist(data[grep("Date", data)])
    date <- date[grep("Date", date)]
    date <- strsplit(date, split = " ")[[1]][2]
    start.dates[j] <- date
  }
  
  start.dates <- dmy(start.dates)
  starts <- ymd_hms(paste0(start.dates, start.times))
  
  starts.offset <- difftime(starts[1:5], starts[1], units = "hours")
  # The first few data sets were always in intervals of half an hour, so now
  # the script expects them to be in those intervals. 
  starts.offset <- as.numeric(round(starts.offset*2)/2)
  
  ### Work out how many intervals are in each fragment
  n.intervals <- vector(mode = "numeric", length = length(sheet.names))
  for(j in 1:length(sheet.names)){
    n.intervals[j] <- length(sheet.names[[j]])
  }
  n.intervals
  
  ### Generate new sheet names
  new.sheet.names <- c()
  a <- 1
  for(j in 1:length(n.intervals)){
    start <- starts.offset[j]
    intervals <- seq(from = start, by = 0.5, length.out = n.intervals[j])
    
    ## Calculate hours
    h <- floor(intervals)
    
    ## Calculate minutes
    min <- rep(0, n.intervals[j])
    min[intervals != round(intervals)] <- 30
    
    
    timing <- paste("(",h," h", sep = "")
    for(k in 1:length(timing)){
      timing[k] <- paste(a, timing[k], sep = " ")
      a <- a+1
      if(min[k] == 30) timing[k] <- paste(timing[k], " 30 min", sep = "")
    }
    timing <- paste(timing, ")", sep = "")
    
    
    new.sheet.names <- c(new.sheet.names, timing)
    ### FIXME: Sheetnames will always be unique, but may refer to the same time point
    ### i.e. if file 1 has 2 sheets, and file 2 starts 35 minutes after file 1
    ### this may not cause a problem for ggplot
  }
  new.sheet.names <- paste("Interval", new.sheet.names, sep = " ")
  
  #### Compile fragmented data files ####
  ### Create a new workbook to contain the data ###
  wb <- createWorkbook()
  
  a <- 1
  for(j in 1:length(matching.files)){
    ### Choose the relevant file and sheet names
    file <- file.path("Fragmented_Raw_Data", 
                     matching.files[j])
    sheets <- sheet.names[[j]]
    
    ### Copy each sheet to the new workbook, with the new sheet name
    for(k in 1:length(sheets)){
      addWorksheet(wb = wb, 
                   sheet=new.sheet.names[a])
      writeData(wb = wb, 
                sheet = new.sheet.names[a],
                colNames = FALSE,
                rowNames = FALSE,
                startRow = 3, # not sure if is required, but other files start on row 3
                read.xlsx(xlsxFile = file,
                          sheet = sheets[k],
                          colNames = FALSE,
                          rowNames = FALSE))
      a <- a+1
    }
  }
  
  ### Add Protocol information 
  # Protocol information is used by the map_raw_data script to ensure the 
  # file name matches the protocol measured by the script. 
  file <- file.path("Fragmented_Raw_Data", 
                    matching.files[1])
  
  addWorksheet(wb = wb,
               sheet = "Protocol Information")
  
  writeData(wb = wb,
            sheet = "Protocol Information",
            colNames = FALSE,
            rowNames = FALSE,
            read.xlsx(xlsxFile = file,
                      sheet = "Protocol Information",
                      colNames = FALSE,
                      rowNames = FALSE))
  
  #### Export workbook for later analysis ####
  saveWorkbook(wb,
               output.file.name, overwrite = T)
}

#####