##### Setup #####
#### Misc ####
plate.reader.link <- '~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/'
plate.reader.data <- 'Plate_Reader_Data/Normalised_Data/'
plate.reader.plots<- 'Plate_Reader_Plots/Normalised_Plots/'
setwd(paste(plate.reader.link,
            "Plate_Reader_Data/Normalised_Data/",
            sep = ""))


#####


# Find all data files
# ensure they are all represented in thee lists.
# combine all files 

### FIXME: some experiments don't have OD750 files: normalised_acid_OD600_Plate1_2020-11-20.csv
### FIXME: One raw data file is missing half the values. These NAs are messing up the analysis

##### Code #####

### Define which files belong to which data set

experiments <- list(SAB_vardeg_pH5.5_750 = c(
  # Batch 1
  "normalised_GC25_OD750_Plate1_2021-05-28.csv",
  "normalised_GC30_OD750_Plate1_2021-06-18.csv",
  "normalised_GC37_OD750_Plate1_2021-07-16.csv",
  "normalised_GC42_OD750_Plate1_2021-08-06.csv",
  # Batch 2
  "normalised_GC25_OD750_Plate1_2021-08-13.csv",
  "normalised_GC30_OD750_Plate1_2021-09-03.csv",
  "normalised_GC37_OD750_Plate1_2021-09-07.csv",
  "normalised_GC42_OD750_Plate1_2021-09-12.csv",
  # Batch 3
  "normalised_GC30_OD750_Plate1_2022-02-06.csv",
  "normalised_GC37_OD750_Plate1_2022-02-11.csv",
  # Batch 4
  "normalised_GC30_OD750_Plate1_2022-02-25.csv",
  "normalised_GC37_OD750_Plate1_2022-03-04.csv",
  # Batch 5
  "normalised_GC37_OD750_Plate1_2022-10-07.csv",
  "normalised_GC30_OD750_Plate1_2022-10-14.csv"),
  SAB_vardeg_pH5.5_600 = c(
    # Pre-standardisation
    "normalised_GC37_OD600_Plate1_2020-10-09.csv", # has shaking, multiple media
    "normalised_GC42_OD600_Plate1_2020-10-30.csv", # has shaking, multiple media
    "normalised_GC30_OD600_Plate1_2020-11-06.csv", # has shaking, multiple media
    # Batch 1
    "normalised_GC25_OD600_Plate1_2021-05-28.csv",
    "normalised_GC30_OD600_Plate1_2021-06-18.csv",
    "normalised_GC37_OD600_Plate1_2021-07-16.csv",
    "normalised_GC42_OD600_Plate1_2021-08-06.csv",
    # Batch 2
    "normalised_GC25_OD600_Plate1_2021-08-13.csv",
    "normalised_GC30_OD600_Plate1_2021-09-03.csv",
    "normalised_GC37_OD600_Plate1_2021-09-07.csv",
    "normalised_GC42_OD600_Plate1_2021-09-12.csv",
    # Batch 3
    "normalised_GC30_OD600_Plate1_2022-02-06.csv",
    "normalised_GC37_OD600_Plate1_2022-02-11.csv",
    # Batch 4
    "normalised_GC30_OD600_Plate1_2022-02-25.csv",
    "normalised_GC37_OD600_Plate1_2022-03-04.csv",
    # Batch 5
    "normalised_GC37_OD600_Plate1_2022-10-07.csv",
    "normalised_GC30_OD600_Plate1_2022-10-14.csv"
    
  ),
  YPD_vardeg_pH5.5_750 = c(
    "normalised_GC37_OD750_Plate1_2022-09-30.csv",
    "normalised_GC42_OD750_Plate1_2022-09-16.csv",
    "normalised_GC30_OD750_Plate1_2022-09-09.csv"
  ),
  YPD_vardeg_pH5.5_600 = c(
    "normalised_GC37_OD600_Plate1_2022-09-30.csv",
    "normalised_GC42_OD600_Plate1_2022-09-16.csv",
    "normalised_GC30_OD600_Plate1_2022-09-09.csv"
  ),
  SAB_skintemp_pH5.5_750 = c(
    "normalised_GC33_OD750_Plate1_2023-01-27.csv",
    "normalised_GC33_OD750_Plate1_2023-01-20.csv",
    "normalised_GC33_OD750_Plate1_2023-01-13.csv"
  ),
  SAB_skintemp_pH5.5_600 = c(
    "normalised_GC33_OD600_Plate1_2023-01-27.csv",
    "normalised_GC33_OD600_Plate1_2023-01-20.csv",
    "normalised_GC33_OD600_Plate1_2023-01-13.csv"
  ),
  SAB_25deg_pHvar_750 = c(
    # Batch 1
    "normalised_ACID_OD750_Plate1_2022-03-11.csv",
    "normalised_ACID_OD750_Plate1_2022-03-18.csv",
    "normalised_ACID_OD750_Plate1_2022-04-01.csv",
    "normalised_ACID_OD750_Plate1_2022-04-08.csv",
    "normalised_ACID_OD750_Plate1_2022-04-15.csv",
    # Batch 2 
    "normalised_ACID_OD750_Plate1_2022-06-11.csv",
    "normalised_ACID_OD750_Plate1_2022-06-17.csv",
    "normalised_ACID_OD750_Plate1_2022-07-01.csv",
    "normalised_ACID_OD750_Plate1_2022-07-08.csv",
    "normalised_ACID_OD750_Plate1_2022-07-15.csv",
    # Batch 3
    "normalised_ACID_OD750_Plate1_2022-11-11.csv",
    "normalised_ACID_OD750_Plate1_2022-12-02.csv"
  ),
  SAB_25deg_pHvar_600 = c(
    # Batch 1
    "normalised_ACID_OD600_Plate1_2022-03-11.csv",
    "normalised_ACID_OD600_Plate1_2022-03-18.csv",
    "normalised_ACID_OD600_Plate1_2022-04-01.csv",
    "normalised_ACID_OD600_Plate1_2022-04-08.csv",
    "normalised_ACID_OD600_Plate1_2022-04-15.csv",
    # Batch 2 
    "normalised_ACID_OD600_Plate1_2022-06-11.csv",
    "normalised_ACID_OD600_Plate1_2022-06-17.csv",
    "normalised_ACID_OD600_Plate1_2022-07-01.csv",
    "normalised_ACID_OD600_Plate1_2022-07-08.csv",
    "normalised_ACID_OD600_Plate1_2022-07-15.csv",
    # Batch 3
    "normalised_ACID_OD600_Plate1_2022-11-11.csv"
    # "normalised_ACID_OD600_Plate1_2022-12-02.csv" # Only half the plate had data collected, 
                                                    # which caused problems for the pipeline.
                                                    # Removal was the most expedient solution
  ),
  SAB_37deg_pHvar_600 = c( # Need to check experimental shaking
    "normalised_acid_OD600_Plate1_2020-11-20.csv", # experimental shaking
    "normalised_acid_OD600_Plate1_2021-03-26.csv", # no shaking
    "normalised_acid_OD600_Plate1_2021-03-19.csv", # no shaking
    "normalised_acid_OD600_plate1_2021-03-12.csv", # no shaking
    "normalised_acidCheck_OD600_Plate1_2020-11-13.csv", # experimental shaking
    "normalised_acidSAB_OD600_Plate1_2020-10-16.csv", # experimental shaking
    "normalised_acidSAB_OD600_Plate1_2020-10-23.csv" # experimental shaking
    ),
  YPD_25deg_pHvar_750 = c( # Only one replicate done
    "normalised_ACID_OD750_Plate1_2022-09-02.csv",
    "normalised_ACID_OD750_Plate1_2022-08-19.csv"
    ),
  YPD_25deg_pHvar_600 = c( # Only one replicate done
    "normalised_ACID_OD600_Plate1_2022-09-02.csv",
    "normalised_ACID_OD600_Plate1_2022-08-19.csv"
  )
)






### Compile data

for(i in 1:length(experiments)){
  
  ### Choose the files from a given experiemnt
  files <- experiments[[i]]
  
  if(length(files) < 1) next
  
  ### Load the first data file, 
  compiled.data <- read.csv(files[1])
  
  ### Add all subsequent data files
  if(length(files) > 1){
    for(j in 2:length(files)){
      compiled.data <- rbind(compiled.data, read.csv(files[j]))
    }
  }
  
  ### Create name, with relevant information
  filename <- paste("../Compiled_Data/",
                    names(experiments)[i], 
                    ".csv",
                    sep = "")
  
  ### Export compiled data
  write.csv(compiled.data,
            filename,
            row.names = FALSE)
}



### Finding mismatching files
# setdiff(names(compiled.data), names(read.csv(files[j])))
# setdiff(names(read.csv(files[j])), names(compiled.data))
# files[1]
