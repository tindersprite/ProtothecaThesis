##### Setup #####
#### Packages ####
require(tidyverse)

#### Functions ####
sem <- function(x, length = NA){
  if(is.na(length)) length <- length(x)
  sd(x)/sqrt(length)
}

wrapper <- function(x, ...) {
  paste(strwrap(x, ...), collapse = "\n")
}

shorten.species.name <- function(long.spp.names){
  long.spp.names <- as.character(long.spp.names)
  full.names   <- strsplit(long.spp.names, " ")
  full.names   <- do.call(rbind, full.names)
  
  genus.abbrev <- substring(full.names[,1],1,1)
  species <- full.names[,2]
  
  short.spp.names <- paste(genus.abbrev, ". ", species, sep = "")
  return(short.spp.names)
}

#### Levels  ####
strains <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")
species.order <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")
timepoint.levels <- read.delim('~/Documents/Main_Project/Scripts/update_platereader_data/time_point_levels.txt', 
                               header = F, sep = ".")

#### Labels ####
dilution.labs <- c('1:1',"1:10", "1:100", "1:1000", "1:10,000")
dilution.labs.alt <-c('1e0',"1e-1", "1e-2", '1e-3', '1e-4')
names(dilution.labs) <- c(1,0.1, 0.01, 0.001, 0.0001)
names(dilution.labs.alt) <- c(1,0.1, 0.01, 0.001, 0.0001)

temperature.labs <- c("37°C", "30°C", "25°C")
names(temperature.labs) <- c(37,30,25)

#### Misc ####
strain.to.species.df <- strains[c(1,2)]
names(strain.to.species.df)[1] <- "Strain"
#####



#### Settings ####

rm_NAs <- T
rm_zero_vals <- T
rm_low_conf_vals <- T
rm_late_vals <- T # "late" values are measurements taken after 96 hours
rm_strain_1replicate <- F

##### Analysis #####
#### Load data ####
cfu_data <- read.csv("Data/RawCFUs.csv")

### Identify strains to species level
cfu_data <- merge(cfu_data, strain.to.species.df)

#### Apply Settings ####
### Remove missing values ###
# missing values cause errors when calculating statistics: mean(c(1,2,NA)) == NA
NAs_in_data <- sum(is.na(cfu_data$CFUs))
if(rm_NAs){
  cfu_data <- cfu_data[!(is.na(cfu_data$CFUs)),]
}

### Remove zero values ###
# zero values artificially reduce means:
## mean(c(10,20)) == 15
## mean(c(10,20, 0)) == 10
## Need to remove zero values from measurements where there are non-zero 
## values, as these will skew means, but keep zero values where there are 
## no non-zero values, as some timepoints genuinely have no CFUs. 

if(rm_zero_vals){
  for(st in unique(cfu_data$Strain)){
    for(ti in unique(cfu_data$Timepoint)){
      for(te in unique(cfu_data$ExperimentalTemperature)){
        for(da in unique(cfu_data$Date)){
          rows <- cfu_data$Strain == st &
            cfu_data$Timepoint == ti &
            cfu_data$ExperimentalTemperature == te &
            cfu_data$Date == da
          
          if(sum(rows) == 0) next # Skip any combinations with no data
          
          if(any(cfu_data[rows,]$CFUs == 0)){
            if(all(cfu_data[rows,]$CFUs == 0)){next 
            } else {
              cfu_data$CFUs[cfu_data$CFUs == 0 & rows] <- "Zero"
            }
          }
        }
      }
    }
  }
}

zeros_in_data <- sum(cfu_data$CFUs == "Zero")
if(rm_zero_vals){
  cfu_data <- cfu_data[cfu_data$CFUs != "Zero",]
  cfu_data$CFUs <- as.numeric(cfu_data$CFUs)
}

### Remove low confidence values ###
# While counting, there will be some values I was uncertain of while counting
low_conf_in_data <- sum(!cfu_data$ConfidentInCount)
if(rm_low_conf_vals){
  cfu_data <- cfu_data[cfu_data$ConfidentInCount,]
}

### Remove values after 4 days incubation
if(rm_late_vals){
  cfu_data <- cfu_data[cfu_data$Timepoint<100,]
}

### Remove values from species with only one replicate
if(rm_strain_1replicate){
  for(st in unique(cfu_data$Strain)){
    if(length(unique(subset(cfu_data, Strain == st)$Date)) < 2){
      cfu_data <- cfu_data[(cfu_data$Strain != st),]
    }
  }
#  
}


#### Process Data ####
### Remove unwanted rows:
## don't want rows with CFUs that are too high to reliably count
## This should be covered by removing values which are low confidence and
## removing rows which have zeros in them

filtered_cfu_data <- data.frame()
for(strain in unique(cfu_data$Strain)){
  for(timepoint in unique(cfu_data$Timepoint)){
    for(temperature in unique(cfu_data$ExperimentalTemperature)){
      for(date in unique(cfu_data$Date)){
        sub.data <- subset(cfu_data,
                           Strain == strain &
                             Timepoint == timepoint &
                             ExperimentalTemperature == temperature &
                             Date == date)
        if(nrow(sub.data) == 0) next
        
        # Remove counts less than ten, if there are larger values, as these will
        # contain more error
        if(any(sub.data$CFUs > 10)){
          sub.data <- sub.data[sub.data$CFUs>=10,]
        } else {
          # If there are no counts less than ten, but there are counts that
          # aren't zero, remove the zeros, as they will skew the mean
          if(any(sub.data$CFUs > 0)){
            sub.data <- sub.data[sub.data$CFUs > 0,]
          }
        }
        
        # Remove counts greater than 250, if there are smaller values, as these will
        # contain more error
        if(any(sub.data$CFUs < 250)) sub.data <- sub.data[sub.data$CFUs<=250,]
        
        filtered_cfu_data <- rbind(sub.data, filtered_cfu_data)
      }
    }
  }
}

cfu_data <- filtered_cfu_data

### Calculate cells per ml for each measurement
cfu_data$CellsPerMl <- ((cfu_data$CFUs/cfu_data$VolSpotted)*1000)/cfu_data$Dilution

### Combine measurements from the same time point and date
cfus_grouped <- group_by(cfu_data, Strain, Date, Timepoint, ExperimentalTemperature) %>%
  summarise(CellsPerMl = mean(CellsPerMl))

### Calculate average cells per ml for each time point

cfus_grouped %>%
  group_by(Strain, Timepoint, ExperimentalTemperature) %>%
  summarise(MeanCFUs = mean(CellsPerMl),
            SDCFUs = sd(CellsPerMl),
            SEMCFUs = sem(CellsPerMl, length = length(unique(Date))),
            Replicates = length(unique(Date))) ->
  cfu_summary

### NAs will be introduced where there is only one replicate, because there is no
### SD or SEM for a single number
cfu_summary[is.na(cfu_summary)] <- 0

# An alternative way to do the same thing:
#cfu_data %>%
#  group_by(Strain, Timepoint) %>%
#  summarise_at(vars(CellsPerMl), list(MeanCFUs = mean,
#                                      SDCFUs = sd))

cfu_data %>%
  inner_join(cfu_summary) ->
  cfu_data

#### Set Levels ####

cfu_data$Dilution <- factor(cfu_data$Dilution,
                            levels = c(1,0.1,0.01,0.001,0.0001))
cfu_data$Species <- factor(cfu_data$Species,
                           levels = species.order$Descent)
cfu_data$Strain <- factor(cfu_data$Strain,
                          levels = strains$H_number)

#### Provide summary statistics ####
data_summary <- data.frame(
  PossibleError = c("Too high CFU", "0 CFU", "Low Confidence", "Late Timepoint"),
  NumberPresent = c(NAs_in_data, zeros_in_data, low_conf_in_data, NA),
  IsRemoved = c(rm_NAs, rm_zero_vals, rm_low_conf_vals, rm_late_vals)
)

#### Plot ####
title.width <- 40
nstrains <-  length(unique(cfu_data$Strain))
ntemps <- length(unique(cfu_data$ExperimentalTemperature))
ndilutions <- length(unique(cfu_data$Dilution))

(deathcurve.individual <- ggplot(data = cfu_data, aes(x = Timepoint, y = CellsPerMl)) +
  geom_point(aes(col = Strain, pch = as.character(ExperimentalTemperature))) +
  geom_smooth(aes(col = Strain), method = "lm") +
  facet_grid(Species ~ Dilution, scales = "free",
             labeller = labeller(Dilution = dilution.labs)) +
  scale_x_continuous(name = "Hours", breaks = seq(0, 
                                                  max(cfu_data$Timepoint), 
                                                  by = 24)) + 
  scale_y_continuous(name = "Cells per ml", labels = scales::comma) +
  ggtitle(wrapper("37°C 'Growth' Curve for temperature - sensitive strains",
                  width = title.width),
          subtitle = "Measurements plotted individually") )

ggsave("Outputs/37°C-deathcurve_individual.png",
       plot = deathcurve.individual,
       device = "png",
       width = ndilutions *2+1,
       height = nstrains *2 + 1)



(deathcurve.grid.titleless <- ggplot(data = cfu_data, aes(x = Timepoint, y = MeanCFUs)) +
  geom_pointrange(aes(ymin = MeanCFUs-(2*SEMCFUs), 
                      ymax = MeanCFUs+(2*SEMCFUs),
                      col = Strain)) + 
  geom_line(aes(col = Strain)) + 
  facet_grid(Species ~ ExperimentalTemperature, scales = "free",
             labeller = labeller(ExperimentalTemperature = temperature.labs)) + 
  scale_x_continuous(name = "Hours", breaks = seq(0, 
                                                  max(cfu_data$Timepoint), 
                                                  by = 24)) + 
  scale_y_continuous(name = "Cells per ml", labels = scales::comma) +
    expand_limits(y = 0)
  )



(deathcurve.wrapped.titleless <- ggplot(data = cfu_data, aes(x = Timepoint, y = MeanCFUs)) +
    geom_pointrange(aes(ymin = MeanCFUs-(2*SEMCFUs), 
                        ymax = MeanCFUs+(2*SEMCFUs),
                        col = Strain)) + 
    geom_line(aes(col = Strain)) + 
    facet_wrap( ~ paste(Species, " ", 
                        ExperimentalTemperature, "°C",
                        sep = ""), 
               scales = "free",
               nrow = length(unique(cfu_data$Species))/2,
               labeller = labeller(ExperimentalTemperature = temperature.labs)) + 
    scale_x_continuous(name = "Hours", breaks = seq(0, 
                                                    max(cfu_data$Timepoint), 
                                                    by = 24)) + 
    scale_y_continuous(name = "Cells per ml", labels = scales::comma) +
    expand_limits(y = 0)
  )


(deathcurve.stats.grid <- deathcurve.grid.titleless + 
    geom_text(aes(label = Replicates)) +
  ggtitle(wrapper("37°C 'Growth' Curve for temperature - sensitive strains",
                  width = title.width),
          subtitle = "Summary of measurements plotted. Bars ±2SE.\nNumbers indicate replicates."))


(deathcurve.stats.wrapped <- deathcurve.wrapped.titleless + 
    geom_text(aes(label = Replicates)) +
    ggtitle(wrapper("37°C 'Growth' Curve for temperature - sensitive strains",
                    width = title.width),
            subtitle = "Summary of measurements plotted. Bars ±2SE.\nNumbers indicate replicates."))

ggsave("Outputs/37°C-deathcurve_grid_notitle.png",
       plot = deathcurve.grid.titleless,
       device = "png",
       width = 8,
       height = nstrains *1.6)

ggsave("Outputs/37°C-deathcurve_wrapped_notitle.png",
       plot = deathcurve.wrapped.titleless,
       device = "png",
       width = nstrains *2.1 + 1,
       height = 5.5)

ggsave("Outputs/37°C-deathcurve_summary_grid.png",
       plot = deathcurve.stats.grid,
       device = "png",
       width = 8,
       height = nstrains *1.6 + 1)

ggsave("Outputs/37°C-deathcurve_summary_wrapped.png",
       plot = deathcurve.stats.wrapped,
       device = "png",
       width = nstrains *2.1 + 1,
       height = 5.5)

(deathcurve.wrapped2.titleless <- ggplot(data = cfu_data, aes(x = Timepoint, y = MeanCFUs)) +
    geom_pointrange(aes(ymin = MeanCFUs-(2*SEMCFUs), 
                        ymax = MeanCFUs+(2*SEMCFUs),
                        col = Strain, shape = as.character(ExperimentalTemperature))) + 
    geom_line(aes(col = Strain, group = paste(ExperimentalTemperature, Strain))) + 
    facet_wrap( ~ Species, 
                scales = "free",
                nrow = ceiling(sqrt(length(unique(cfu_data$Species)))),
                labeller = labeller(ExperimentalTemperature = temperature.labs)) + 
    scale_x_continuous(name = "Hours", breaks = seq(0, 
                                                    max(cfu_data$Timepoint), 
                                                    by = 24)) + 
    scale_y_continuous(name = "Cells per ml", labels = scales::comma) +
    guides(shape = guide_legend(title = "Temperature",
                                order = 1),
           colour = guide_legend(title = "Strain",
                               order = 2))+
    expand_limits(y = 0) +
    theme(strip.text = element_text(face = "italic"))
)

nspecies <- length(unique(cfu_data$Species))
ggsave("Outputs/37°C-deathcurve_wrapped2_notitle.png",
       plot = deathcurve.wrapped2.titleless,
       device = "png",
       width = ceiling(sqrt(nspecies)) *2.5 + 1,
       height = ceiling(sqrt(nspecies)) *2)


(deathcurve.wrapped3.titleless <- ggplot(data = cfu_data, aes(x = Timepoint, y = MeanCFUs)) +
    geom_pointrange(aes(ymin = MeanCFUs-(2*SEMCFUs), 
                        ymax = MeanCFUs+(2*SEMCFUs),
                        col = Strain, shape = paste0(as.character(ExperimentalTemperature), "°C"))) + 
    geom_line(aes(col = Strain, group = paste(ExperimentalTemperature, Strain))) + 
    facet_wrap( ~ shorten.species.name(Species), 
                scales = "free",
                nrow = 2) + 
    scale_x_continuous(name = "Hours", breaks = seq(0, 
                                                    max(cfu_data$Timepoint), 
                                                    by = 24)) + 
    scale_y_continuous(name = "Cells per ml", labels = scales::comma) +
    guides(shape = guide_legend(title = "Temperature",
                                order = 1),
           colour = guide_legend(title = "Strain",
                                 order = 2))+
    expand_limits(y = 0) +
    theme(strip.text = element_text(face = "italic"))
)

ggsave("~/Documents/Main_Project/Admin/Presentations/PostdocDylan/TemperatureDeathCurves.png",
       deathcurve.wrapped3.titleless,
        width = 30,
       height = 12, 
       units = "cm")

#####
data_summary


#### Statistical comparison ####
### Identify data for comparison
head(cfus_grouped)
cfus_grouped 
cfus_grouped <- left_join(cfus_grouped, strain.to.species.df)
arrange(cfus_grouped, Strain, Timepoint, ExperimentalTemperature)

### Create empty data frame for results
n.rows <- length(unique(cfus_grouped$Species)) * length(unique(cfus_grouped$Timepoint))
cfus_compared <- data.frame(Species = character(length = n.rows), 
                            Timepoint = numeric(length = n.rows),
                            Pvalue = numeric(length = n.rows))

### Compare between temperatures at each time point, for each species
row <- 0
for(sp in unique(cfus_grouped$Species)){
  for(tp in unique(cfus_grouped$Timepoint)){
    row <- row + 1
    
    cfus_compared$Species[row] <- sp
    cfus_compared$Timepoint[row] <- tp
  
    control <- subset(cfus_grouped,  
                      ExperimentalTemperature == 25 & 
                        Species == sp & 
                        Timepoint == tp)$CellsPerMl
    experiment <- subset(cfus_grouped, ExperimentalTemperature == 37& 
                           Species == sp & 
                           Timepoint == tp)$CellsPerMl
    
    if(length(control) < 2 | length(experiment) < 2){
      cfus_compared$Pvalue[row] <- NA
    } else {
      cfus_compared$Pvalue[row] <- t.test(control, experiment)$p.value
      # cfus_compared$Pvalue[row] <- wilcox.test(control, experiment)$p.value
    }
    
  }
}
### Adjust for multiple comparisons 
cfus_compared$PAdj <- p.adjust(cfus_compared$Pvalue, method = "BH")

### Indicate significance with text
cfus_compared$Significance <- "NS"
cfus_compared$Significance[cfus_compared$PAdj <= 0.05] <- "*"
cfus_compared$Significance[cfus_compared$PAdj <= 0.01] <- "**"
cfus_compared$Significance[cfus_compared$PAdj <= 0.001] <- "***"


head(cfu_data)

### Choose heights for the text
unique(select(cfu_data, Species, MeanCFUs,SEMCFUs, Timepoint, ExperimentalTemperature)) %>%
  mutate(UpperPoint = (MeanCFUs + 2*SEMCFUs) * 1.2) %>%
  group_by(Species, Timepoint) %>%
  summarise(InitialTextHeight = max(UpperPoint)) %>%
  group_by(Species) %>%
  summarise(InitialTextHeight = InitialTextHeight,
            Timepoint = Timepoint,
            HighestPoint = max(InitialTextHeight)) %>%
  mutate(TextHeight = InitialTextHeight + 0.05*HighestPoint) %>%
  select(Species, Timepoint, TextHeight) %>%
  right_join(cfus_compared) -> cfus_compared


### Set levels to fix panel order
cfus_compared$Species <- factor(cfus_compared$Species, levels = species.order$Descent)

### Produce plot
deathcurve.wrapped2.annotated.titleless <- deathcurve.wrapped2.titleless +
  geom_text(data = cfus_compared, aes(x = Timepoint, label = Significance, y = TextHeight)) 
  
ggsave("Outputs/37°C-deathcurve_wrapped2-annotated_notitle.png",
       plot = deathcurve.wrapped2.annotated.titleless,
       device = "png",
       width = ceiling(sqrt(nspecies)) *2.5 + 1,
       height = ceiling(sqrt(nspecies)) *2)
