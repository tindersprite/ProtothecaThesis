#####

setwd("~/Documents/Main_Project/Thesis/Figures/Ch1 Identification/")


rev.comp <- function(input.seq){
  bases <- c("A", "T", "G", "C", "N", "M", 'K', "R", 'Y', 'W', 'S', 'V', 'B', 'H', 'D')
  pairs <- c("T", "A", "C", "G", "N", 'K', "M", 'Y', "R", 'W', 'S', 'B', 'V', 'D', 'H')
  
  input.seq <- toupper(input.seq)
  
  input.seq <- strsplit(input.seq, split = "")[[1]]
  output.seq <- character(length = length(input.seq))
  
  for(i in 1:length(bases)){
    output.seq[input.seq == bases[i]] <- pairs[i]
  }
  
  output.seq <- rev(output.seq)
  
  return(paste(output.seq, collapse = ""))
}

###### UK maps ######

# library(PostcodesioR) ## requires full postcode

#### Setup ####
library(ggplot2)
library(dplyr)
library(rgdal) # might be retired, but unsure what to replace with
library(maptools) # might replaced by the sp package
library(gpclib)
library(openxlsx)
library(viridis)
library(terra)
library(sf)
library(stars)


gpclibPermit()

keep.letters <- function(string){
  split.string <- strsplit(string, split = "")[[1]]
  split.string <- split.string[split.string %in% LETTERS]
  output <- paste(split.string, collapse= "")
  return(output)
}

setwd("~/Documents/Main_Project/Thesis/Figures/Ch1 Identification/")
#####

##### NML sample mapping #####
#### Load NML data ####

nml <- read.xlsx("~/Documents/Main_Project/Collaborations/NML Prototheca Samples.xlsx")


nml <- select(nml, Postcode, Prototheca.Spec, Prototheca.by.Shave)
nml <- na.omit(nml)

nml[nml=="+"] <- TRUE
nml[nml=="-"] <- FALSE
nml$Postcode


for(i in 1:nrow(nml)) nml$Areacode[i] <- keep.letters(nml$Postcode[i])

nml.freq <- as.data.frame(table(nml$Areacode))
names(nml.freq) <- c("name", "nml.freq")
bham.freq <- as.data.frame(table(subset(nml, Prototheca.by.Shave == T)$Areacode))
names(bham.freq) <- c("name", "bham.freq")
milk.freq <- left_join(nml.freq, bham.freq)

milk.freq[is.na(milk.freq)] <- 0

#### Set up UK map ####
# Download UK postcode polygon Shapefile
# download.file(
#   "http://www.opendoorlogistics.com/wp-content/uploads/Data/UK-postcode-boundaries-Jan-2015.zip",
#   "Geography/postal_shapefile"
# )
# unzip("Geography/postal_shapefile")

# Read the downloaded Shapefile from disk

postal <- maptools::readShapeSpatial("Geography/Distribution/Areas")

# Assign each "region" an unique id
postal.count <- nrow(postal@data)
postal@data$id <- 1:postal.count

# Transform SpatialPolygonsDataFrame to regular data.frame in ggplot format
postal.fort <- ggplot2::fortify(postal, region='id')

# Add "region" id to frequency data
milk.freq <- merge(milk.freq, postal@data)

# Merge frequency data onto geogrphical postal polygons
postal.fort <- merge(postal.fort, milk.freq, by = "id", all.x = T, all.y = F)
postal.fort <- postal.fort[order(postal.fort$order),]

postal.fort$nml.freq[is.na(postal.fort$nml.freq)] <- 0
postal.fort$bham.freq[is.na(postal.fort$bham.freq)] <- 0
#### Make maps ####

## Source of samples that the NML found to be positive
(nml.p <- ggplot(postal.fort) + 
  geom_polygon(aes(x = long, y = lat, group = group, fill=nml.freq), colour='white') + 
  coord_fixed(ratio = 1.5) + 
  scale_fill_viridis() + 
  labs(fill = "Number of\nPrototheca\ncases") +
  theme_void())

ggsave("Geography/NML+ve.png",
       nml.p,
       height = 6.2, width = 5)

## Source of samples I could confirm were positive
bham.p <- ggplot(postal.fort) + 
    geom_polygon(aes(x = long, y = lat, group = group, fill=bham.freq), colour='white') + 
    coord_fixed(ratio = 1.5) + 
    scale_fill_viridis() + 
    labs(fill = "Number of\nConfirmed\nPrototheca\ncases") +
    theme_void()

ggsave("Geography/bham+ve.png",
       bham.p,
       height = 6.2, width = 5)
#####

##### UK dairy mapping #####

### Data carpentries
describe("Geography/20230628_BirminghamUni_LDDGDairyCattle/APHA_LDDG_DairyPopulationDensity_2021.tif")
dairy <- rast("Geography/20230628_BirminghamUni_LDDGDairyCattle/APHA_LDDG_DairyPopulationDensity_2021.tif")

summary(dairy)
names(dairy) <- "Cattle"

dairy.df <- as.data.frame(dairy, xy = TRUE)
str(dairy.df)
crs(dairy)

(cow.dens <- ggplot() +
  geom_raster(data = dairy.df , aes(x = x, y = y, fill = Cattle)) +
  scale_fill_viridis_c() +
  coord_quickmap() +
  #  coord_cartesian(xlim = c(0,700000), ylim = c(-100000,1250000))+
  labs(fill = "Cattle per\nsquare Km")+ 
    theme_void())

ggsave("Geography/Dairy cattle density.png",
       cow.dens,
       height = 6.2, width = 5)

