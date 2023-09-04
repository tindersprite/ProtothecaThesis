##### Setup #####
#### Packages ####
require(ggplot2)
require(dplyr)
require(metafolio) # For gg_color_hue
require(openxlsx)

#### Functions ####
df_combine <- function(df, new_col_name, new_vals){
  in_df  <- df
  out_df <- data.frame()
  for(i in 1:length(new_vals)){
    in_df[new_col_name] <- new_vals[i]
    out_df <- rbind(out_df, in_df)
  }
  return(out_df)
}

#### Misc ####
setwd("~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/")

SAB <- "SAB"
PIM <- "PIM"

photosynthetic.genera <- c('Auxenochlorella', 'Chlorella')

possible.plate.size <- c(6,12,24,48,96)
possible.cols <- c(3,4,6,8,12) # the numbers
possible.rows <- c(2,3,4,6,8)  # the letters

plate.reader.filepath <- '~/Documents/Main_Project/Data/Plate_Reader/'
map.plot.filepath <- paste(plate.reader.filepath, 
                            "Plate_Reader_Plots/Plate_Maps/", sep = "")
map.data.filepath<- paste(plate.reader.filepath, 
                          "Plate_Reader_Data/Mapping_Data/", sep = "")

#####

##### Inputs #####
#### Highly variable inputs ####
map.version    <- "1.0.5"
prototheca.strains      <- c(1,2,3,4,18,27,28,29,30,31,32,37,40,45,51,52,53,54) # Just provide numbers. e.g. 1:3 for HP1, HP2, HP3
auxenochlorella.strains <- c(1,2,5,6,7) # Just provide numbers. e.g. 1:3 for HA1, HA2, HA3
helicosporidium.strains <- c() # Just provide numbers. e.g. 1:3 for HH1, HH2, HH3

date           <- Sys.Date() # Format: "YYYY-MM-DD". Set as NA for today
experiment     <- "GC33"
test.condition <- "CultureTemperature" # Options: 
  # StrainID, CultureMedia, PrecultureTemperature, 
  # CultureTemperature, pH, H2O2Concentration

culture.temp   <- 33
ph             <- c("N") # N for not controlled
ph.media.batch <- c(NA)  # one for each pH; 
                    # Format MmmYYYY, e.g. Mar2020 ## paired with ph
h2o2.conc      <- c(0)


save.rep.map <- FALSE

#### Less variable inputs ####

culture.media       <- c("SAB") # SAB or PIM
preculture.media    <- culture.media ## Paired with culture.media
culture.shaking     <- 0 # speed (RPM) the cultures are shaken at
measurement.shaking <- culture.shaking # Speed cultures are shaken before measurement
stock.temp          <- c("4°C")
plate.number        <- 1
rows                <- 6 # the letters
cols                <- 8 # the numbers

#### Fixed inputs ####

culture.vol          <- c(500) # in microlitres
inoculated.density   <- c(0.05) # OD wells are inoculated to
n.replicates         <- 2
preculture.temp      <- c(25)
preculture.shaking   <- c(20) # 150 in shaking incubator, 20 in rotator
stock.media          <- c("PIM")

plate.size           <- 48 # number of wells in the plate
lamp                 <- 2
notes                <- NA #

photosynthetic.OD.wavelength <- 600 # 750 would be better, but we can't do.
# Is 405 better than 600?
nonphotosynthetic.OD.wavelength <- 600

#####

##### Code #####
#### Tidy inputs ####
### Collect strains into one vector ###
strain.IDs <- c("NEG")

if(length(prototheca.strains) > 0){
  prototheca.strains <- paste("HP", sort(prototheca.strains), sep = "")
  strain.IDs <- c(strain.IDs, prototheca.strains)
}

if(length(auxenochlorella.strains) > 0){
  auxenochlorella.strains <- paste("HA", sort(auxenochlorella.strains), sep = "")
  strain.IDs <- c(strain.IDs, auxenochlorella.strains)
}
if(length(helicosporidium.strains) > 0){
  helicosporidium.strains <- paste("HH", sort(helicosporidium.strains), sep = "")
  strain.IDs <- c(strain.IDs, helicosporidium.strains)
}

### Capitalise everything ###
experiment       <- toupper(experiment)
ph               <- toupper(ph)
ph.media.batch   <- toupper(ph.media.batch)
culture.media    <- toupper(culture.media)
preculture.media <- toupper(preculture.media)
stock.media      <- toupper(stock.media)

replicate        <- LETTERS[1:n.replicates]

### Fix missing inputs, if required ###
if(length(ph) == 1){
  if(ph == "N"){
    pH.media.batch  <- NA
  }
}

if(is.na(date)){
  date <- Sys.Date()
}

#### Generate Data Frame ####
df <- data.frame(StrainID = strain.IDs)

### Add vector variables
df <- df_combine(df, "CultureMedia", culture.media)
df <- df_combine(df, 'pH', ph)
df <- df_combine(df, 'H2O2Concentration', h2o2.conc)
df <- df_combine(df, 'CultureMedia', culture.media)
df <- df_combine(df, 'CultureVolume', culture.vol)
df <- df_combine(df, 'InoculatedDensity', inoculated.density)
df <- df_combine(df, 'Replicate', replicate)
df <- df_combine(df, 'PrecultureTemperature', preculture.temp)
df <- df_combine(df, 'PrecultureShaking', preculture.shaking)
df <- df_combine(df, 'StockMedia', stock.media)
df <- df_combine(df, 'StockTemperature', stock.temp)

### Add scalar variables
df$Date <- as.character(date)
df$Experiment <- experiment
df$MapVersion <- map.version    
df$CultureTemperature <- culture.temp   
df$CultureShaking <- culture.shaking   
df$MeasurementShaking <- measurement.shaking 
df$PlateNumber <- plate.number      
df$Nrow <- rows             
df$Ncol <- cols             
df$PlateSize<- plate.size     
df$Lamp <- lamp           
df$Notes <- notes           

### Add paired variables ###

df$pHMediaBatch <- NA
for(i in 1:length(ph)){
  df$pHMediaBatch[which(df$pH == ph[i])] <- ph.media.batch[i]
}

df$PrecultureMedia <- NA
for(i in 1:length(culture.media)){
  df$PrecultureMedia[which(df$CultureMedia == culture.media[i])] <- preculture.media[i]
}

df$InocDensWavelength <- NA
photosynthetic.strain.ids <- substring(photosynthetic.genera,1,1)
  # identify the strain identifiers that correspond with photosynthetic organisms
checking.strain.ids <- substring(df$StrainID,2,2)
  # identify relevant elements to identify photosynthetic organisms
photos <- checking.strain.ids %in% photosynthetic.strain.ids
df$InocDensWavelength[photos] <- photosynthetic.OD.wavelength
df$InocDensWavelength[!photos] <- nonphotosynthetic.OD.wavelength

### Remove irrelevant variables from negative controls ###
df$InoculatedDensity[df$StrainID == 'NEG']     <- NA
df$InocDensWavelength[df$StrainID == 'NEG']    <- NA
df$PrecultureMedia[df$StrainID == 'NEG']       <- NA
df$PrecultureTemperature[df$StrainID == 'NEG'] <- NA
df$PrecultureShaking[df$StrainID == 'NEG']     <- NA
df$StockMedia[df$StrainID == 'NEG']            <- NA
df$StockTemperature[df$StrainID == 'NEG']      <- NA

#### Generate Positions ####
### Define coordinates
# This should allow for not all of the wells to be filled, but for the wells that are
# filled to be placed randomly within the plate
max.cols <- possible.cols[possible.plate.size == plate.size]
max.rows <- possible.rows[possible.plate.size == plate.size]
total.wells  <- plate.size
filled.wells <- cols * rows

if(!filled.wells == nrow(df)){
  print("Length of data frame does not match number of wells. 
        Something is wrong with the inputs or the dataframe.")
}

possible.x <- rep(1:max.cols, length = total.wells) 
possible.y <- rep(LETTERS[max.rows:1], each = max.cols)

### Set the seed
seed <- paste(strsplit(as.character(date), split = "-")[[1]][1:3], collapse = "")
if(plate.number > 1){
  seed <- paste(plate.number, seed, sep = "")
}
seed <- as.numeric(seed)

### Randomise coordinates
set.seed(seed)
coords <- sample(1:total.wells, total.wells)

df$x <- possible.x[coords[1:filled.wells]]
df$y <- possible.y[coords[1:filled.wells]]

#### Generate Mapping Plot ####
### Set factors ###
df$y <- factor(df$y, levels = LETTERS[max.rows:1])
df$x <- factor(df$x, levels = 1:max.cols)
df$StrainID <- factor(df$StrainID, levels = strain.IDs)

# To be able to set the shape of the test condition, it must not be a continuous
# variable. Two square brackets are needed, because one will return a df which
# cannot be converted into factors
df[[test.condition]] <- factor(df[[test.condition]], 
                             levels = unique(df[[test.condition]]))

### Set Colours ###
# To evenly space colours for all the strains, and then add NEG as black
strainpal <- gg_color_hue(length(strain.IDs)-1) 
strainpal <- c('#000000', strainpal)

cond.factor <- max(c(length(culture.media), 
                      length(preculture.temp), 
                      length(culture.temp), 
                      length(ph)))
condpal <- gg_color_hue(cond.factor) 

### Create graph ###

map <- ggplot(data = df, aes(x = x, y = y)) + 
  ggtitle(paste(experiment, date, sep = "_"),
          subtitle = paste("Save data with the file name: ",
                           experiment, "_",
                           "ODxxx", '_',
                           "Plate", plate.number, '_',
                           date, sep = "")) +
  theme(plot.title = element_text(hjust = 0.5))


rep.map <- map + # Optional plot, to show the position of replicates
  geom_point(aes(shape = Replicate), size = 8)

test.map <- map +
  geom_point(aes_string(colour = "StrainID", shape = test.condition), size = 8)+
  scale_shape_manual(values = c(15:19)) +
  scale_colour_manual(values = strainpal)

if(length(unique(df$Strain)) > 9){
strain.map <- map +
  geom_text(aes(colour = StrainID, label = StrainID))+
  scale_colour_manual(values = strainpal)
  strain.map
} else {
  test.map
}

#### Save Map Plot ####
base.filename       <- paste(experiment, plate.number, date, sep = "_")
rep.map.filename    <- paste("rep-map", base.filename, sep = "_")
test.map.filename   <- paste("test-map", base.filename, sep = "_")
strain.map.filename <- paste("strain-map", base.filename, sep = "_")

setwd(map.plot.filepath)

if(save.rep.map){
  ggsave(filename = paste(rep.map.filename,"png", sep = "."),
         plot = rep.map,
         device = "png",
         width = 2*max.cols,
         height = 2*max.rows,
         units = "cm")
}

ggsave(filename = paste(test.map.filename,"png", sep = "."),
       plot = test.map,
       device = "png",
       width = 2*max.cols,
       height = 2*max.rows,
       units = "cm")

if(exists("strain.map")){
  ggsave(filename = paste(strain.map.filename,"png", sep = "."),
         plot = strain.map,
         device = "png",
         width = 2*max.cols,
         height = 2*max.rows,
         units = "cm")
}

#### Save Map Data ####
setwd(map.data.filepath)

df[is.na(df)] <- "NA"

col.order <- c('StrainID', 'Date', 'Experiment', 'x', 'y',
               'CultureMedia', 'CultureTemperature', 'CultureVolume', 'CultureShaking',
               'InoculatedDensity', 'InocDensWavelength', 'MeasurementShaking',
               'PrecultureMedia', 'PrecultureTemperature', 'PrecultureShaking', 
               'StockMedia', 'StockTemperature', 
               'Replicate', 'PlateNumber',
               'pH', 'pHMediaBatch', 'H2O2Concentration',
               'PlateSize', 'Nrow', 'Ncol', 'Lamp',
               'MapVersion', 'Notes')

df <- arrange(df, StrainID, Replicate,
              CultureMedia, CultureTemperature, CultureVolume, CultureShaking, 
              InoculatedDensity, 
              pH, pHMediaBatch, H2O2Concentration, 
              Lamp, MapVersion, Notes)

df <- df[,col.order]

wb <- createWorkbook()

addWorksheet(wb,
             sheetName = "Condition Mapping")

writeData(wb,
          "Condition Mapping",
          df,
          startCol = 1,
          startRow = 1,
          xy = NULL,
          colNames = TRUE,
          rowNames = FALSE)


saveWorkbook(wb, file = paste(experiment,"_",
                              map.version, "_",
                              "Plate", plate.number, "_", 
                              date,".xlsx", sep = ""), 
             overwrite = T)

#saveWorkbook(wb, file = paste("../../Plate_Reader_Backups/Plate_Reader_Data/Mapping_Data/",
#                              experiment,"_",
#                              map.version, "_",
#                              "Plate", plate.number, "_", 
#                              date,".xlsx", sep = ""), 
#             overwrite = T)
#####


revcomp <- function(seq){
  base <- c("A", "T", "G", "C")
  pair <- c("T", "A", "C", "G")
  
  seq.split <- strsplit(seq, split = "*")[[1]]
  seq.pair <- character(length = length(seq.split))
  for(i in 1:length(seq.split)){
    seq.pair[i] <- pair[seq.split[i]==base]
  }
  seq.pair.rev <- rev(seq.pair)
  return(paste(seq.pair.rev, collapse = ""))
  }

xanthor <- read.delim("~/Documents/Main_Project/Data/SangerSequencing/CytB/Pxanthor SAG263-11 cytb.fasta")
names(xanthor) <- NULL
xanthor <-  paste(unlist(xanthor), collapse = "")
revcomp(xanthor)
