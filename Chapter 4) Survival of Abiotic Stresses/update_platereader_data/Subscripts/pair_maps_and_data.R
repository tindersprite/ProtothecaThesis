##### Setup #####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data/")

substrRight <- function(x, n){
  # produces a substring of n characters from the end of x
  substr(x, nchar(x)-n+1, nchar(x))
}



#####

### Import the data frame of data files
mapping_pairs <- read.csv("tmp/mapping_pairs.csv")
mapping_pairs$raw_data_file <- paste(mapping_pairs$raw_data_file, ".xlsx", sep = "")

### Import the names of the map files
maps <- read.csv("tmp/cond_maps.txt", header = F)
maps <- as.character(maps$V1)

trimmed.maps <- substring(maps, 1, nchar(maps)-5)

### Identify the date and plate number of each map
#  strsplit(trimmed.maps, split = "_")[[1]]
position.of.date <- length(strsplit(trimmed.maps, split = "_")[[1]])
position.of.plate.number  <- length(strsplit(trimmed.maps, split = "_")[[1]]) - 1

map.date  <- c()
map.plate <- c()
for(i in 1:length(maps)){
  map.date[i]  <- strsplit(trimmed.maps, split = "_")[[i]][position.of.date]
  map.plate[i] <- strsplit(trimmed.maps, split = "_")[[i]][position.of.plate.number]
}

###
map.column  <- which(colnames(mapping_pairs) == "mapping_data_file")

for(i in 1:length(map.date)){
  for(j in 1:length(map.plate)){
    date  <- map.date[i]
    plate <- map.plate[j]
    data_for_map <- date == mapping_pairs$date & plate == mapping_pairs$plate_number
    
    mapping_pairs[date == mapping_pairs$date, map.column] <- maps[i]
  }
}

write.csv(mapping_pairs,
          file = "tmp/mapping_pairs.csv",
          row.names = FALSE)

if(any(is.na(mapping_pairs[,map.column]))){
  write(x = "An Error Has Ocurred. 
NAs remain in the mapping_pairs.csv file.
A condition map appears to be missing.", 
        file = "error_message")
}
# Define columns that should exist

# find all the mapping info

# iterate through columns, ensure that each map has the correct column

# If a map is missing a column, add it and provide that information
## Keep note of every map that is updated

# Once all the maps have the correct columns, regenerate raw data

# Once raw data is combined with maps, create a normalised OD column

# Save normalised data

# Compile all the normalised data together 