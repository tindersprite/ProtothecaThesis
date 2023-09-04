##### Setup #####
#### Packages ####
require(tidyverse)
require(msir) # For loess prediction intervals
require(ggpubr)
require(growthcurver)
# require(growthcurver) # requires a format my data doesn't have, 
#### Functions ####
shorten.species.name <- function(long.spp.names){
  long.spp.names <- as.character(long.spp.names)
  full.names   <- strsplit(long.spp.names, " ")
  full.names   <- do.call(rbind, full.names)
  
  genus.abbrev <- substring(full.names[,1],1,1)
  species <- full.names[,2]
  
  short.spp.names <- paste(genus.abbrev, ". ", species, sep = "")
  return(short.spp.names)
}

#### Filepaths ####
plate.reader.link <- '~/Documents/Main_Project/Scripts/update_platereader_data/plate_reader_link/'
plate.reader.data <- paste(plate.reader.link, 'Plate_Reader_Data/Compiled_Data/', sep = "")


#### Levels ####
# created with $echo {0..100}_h_{0,30}_min | tr " " "." | tr _ " " > time_point_levels.txt
time.point.levels <- read.delim('~/Documents/Main_Project/Scripts/update_platereader_data/time_point_levels.txt', header = F, sep = ".")
time.point.levels <- unlist(time.point.levels)
names(time.point.levels) <- NULL
time.point.levels <- time.point.levels[-length(time.point.levels)]

strains <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")
strain.levels <- strains$H_number

species.levels <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")

#species.levels.pathology <- c("Prototheca bovis","Prototheca wickerhamii","Prototheca blaschkeae",
#                              "Prototheca cutis","Prototheca miyajii","Prototheca ciferrii",
#                              "Prototheca pringsheimii","Prototheca zopfii",
#                              "Prototheca cerasi","Prototheca cookei",
#                              "Prototheca tumulicola","Prototheca stagnora",
#                              "Prototheca moriformis","Prototheca paracutis",
#                              "Prototheca xanthoriae","Auxenochlorella protothecoides",
#                              "Auxenochlorella symbiontica","Helicosporidium spp.")
#species.levels.descent   <- c("Prototheca bovis","Prototheca pringsheimii","Prototheca zopfii",
#                              "Prototheca cerasi","Prototheca ciferrii","Prototheca cookei",
#                              "Prototheca blaschkeae","Prototheca tumulicola","Prototheca stagnora",
#                              "Prototheca moriformis","Prototheca cutis","Prototheca paracutis",
#                              "Prototheca miyajii","Prototheca xanthoriae","Auxenochlorella protothecoides",
#                              "Auxenochlorella symbiontica","Prototheca wickerhamii","Helicosporidium spp.")
#species.levels.alphabet <- c('Auxenochlorella protothecoides',"Auxenochlorella symbiontica",
#                             'Prototheca blaschkeae','Prototheca bovis',
#                             'Prototheca ciferrii','Prototheca cookei','Prototheca cutis',
#                             'Prototheca miyajii','Prototheca paracutis','Prototheca tumulicola',
#                             'Prototheca wickerhamii','Prototheca xanthoriae')
#

original.species.level <- NA

infection.levels <- c("TRUE", "FALSE", "UNKNOWN", "SYMBIOSIS")

#### Labels ####
species.labels <- shorten.species.name(species.levels$Alphabet) 
names(species.labels) <- species.levels$Alphabet
species.labels[species.labels == "A. sp.1"] <- "Auxenochlorella sp. 1"
species.labels[species.labels == "A. sp.2"] <- "Auxenochlorella sp. 2"

#### Working directory ####
setwd("~/Documents/Main_Project/Thesis/Figures/Ch2 Abiotic Stress/Pt1 Temperature stress/")

#### Misc ####
strain.to.species.df <- strains[c(1,2)]
names(strain.to.species.df)[1] <- "Strain"
#####

##### Plate reader OD plots #####

### Identify files that contain relevant data:
data.files <- c("SAB_vardeg_pH5.5_750.csv",
                "YPD_vardeg_pH5.5_750.csv",
                "SAB_vardeg_pH5.5_600.csv",
                "YPD_vardeg_pH5.5_600.csv")

#####
for(dataset in data.files){
  #### Import data ####
  ### Load dataset
  # dataset <- data.files[1]
  GC <- read.csv(paste(plate.reader.data, dataset, sep = ""))
  
  ### Remove any entries which were shaken
  GC <- subset(GC, MeasurementShaking == 0)
  
  #### Process Data ####
  ### Take constants from GC data 
  wavelength <- unique(GC$Wavelength)
  media <- unique(GC$CultureMedia)
  
  
  if(media == "SAB"){
    ### Identify species with only one replicate, to be removed
    ### Only to be done for SAB cultures - I only did one YPD replicate
    species_replicates <- GC %>%
      group_by(StrainID, CultureTemperature) %>%
      summarise(N = length(unique(Date))) %>%
      left_join(strain.to.species.df, by = c("StrainID" = "Strain")) %>%
      group_by(Species, CultureTemperature) %>%
      summarise(Replicates = sum(N)) %>%
      pivot_wider(names_from = CultureTemperature, values_from = Replicates)
    
    single.replicate.species <- species_replicates$Species[species_replicates$`30` < 2 | species_replicates$`37` < 2]
    
    GC <- subset(GC, !(SpeciesName %in% single.replicate.species))
  }
  
  ### Count number of strains in experiment
  
  length(unique((GC$StrainID)))
  
  ### Count the number of replicates per species
  species_replicates <- GC %>%
    group_by(StrainID, CultureTemperature, CultureMedia, Wavelength) %>%
    summarise(N = length(unique(Date))) %>%
    left_join(strain.to.species.df, by = c("StrainID" = "Strain")) %>%
    group_by(Species, CultureTemperature, CultureMedia, Wavelength) %>%
    summarise(Replicates = sum(N)) %>%
    pivot_wider(names_from = CultureTemperature, values_from = Replicates)
  
  
  names(species_replicates)[4:ncol(species_replicates)] <- paste(names(species_replicates)[4:ncol(species_replicates)], 
                                                                 "degC", sep = "")
  write.csv(species_replicates,
            paste("replicates_per_species_",
                  "OD", wavelength, "_",
                  media, 
                  ".csv",
                  sep = ""),
            row.names = FALSE)
  
  ### Count the number of replicates per strain
  strain_replicates <- GC %>%
    group_by(StrainID, CultureTemperature, CultureMedia, Wavelength) %>%
    summarise(Replicates = length(unique(Date))) %>%
    pivot_wider(names_from = CultureTemperature, values_from = Replicates)
  
  strain_replicates[is.na(strain_replicates)] <- 0
  
  names(strain_replicates)[4:ncol(strain_replicates)] <- paste(names(strain_replicates)[4:ncol(strain_replicates)], 
                                                               "degC", sep = "")
  
  strain_replicates$StrainID <- factor(strain_replicates$StrainID,
                                       levels = strains$H_number)
  strain_replicates <- strain_replicates[order(strain_replicates$StrainID),]
  
  write.csv(strain_replicates,
            paste("replicates_per_strain_",
                  "OD", wavelength, "_",
                  media, 
                  ".csv",
                  sep = ""),
            row.names = FALSE)
  
  
  
  
  ### 
  GC$Infection <- NA
  GC$Host <- NA
  for(i in 1:nrow(strains)){
    strain.ID <- strains$H_number[i]
    infect <- strains$From_infection[i]
    host <- strains$Host[i]
    species <- strains$Species[i]
    
    rows <- which(GC$StrainID == strain.ID)
    GC$Infection[rows] <- infect
    GC$Host[rows]      <- host
    GC$SpeciesName[rows]<- species
  }
  
  GC$TimePoint <- factor(GC$TimePoint, levels = time.point.levels)
  GC$StrainID  <- factor(GC$StrainID, levels = strain.levels)
  GC$Infection <- factor(GC$Infection, levels = infection.levels)
  GC$Source <- case_when(GC$Infection == "TRUE" ~ "Infection",
                         GC$Infection == "FALSE" ~ "Environment",
                         GC$Infection == "UNKNOWN" ~ "Unknown",
                         GC$Infection == "SYMBIOSIS" ~ "Symbiosis")
  
  #### Set species order ####
  GC$SpeciesName <- factor(GC$SpeciesName, 
                           levels = species.levels$Descent) #Descent, Alphabet, or Pathology
  
  #### Set grid labels ####
  factor.labels <- paste(unique(GC$CultureTemperature), "°C", sep = "")
  names(factor.labels) <- unique(GC$CultureTemperature)
  factor.labels[which(names(factor.labels) == 25)] <- "Room Temperature"
  
  #### (Optional) Subset data ####
  #GC <- subset(GC, CultureTemperature %in% c(30,37))
  #GC <- subset(GC, !(Infection %in% c("UNKNOWN")))
  #GC <- subset(GC, !(Infection %in% c("UNKNOWN", "SYMBIOSIS", "COMMENSAL")))
  
  #### Make Temperature plots by strain ####
  
  GC_bystrain<- ggplot(data = GC, 
                       aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_smooth(aes(y = NormalisedOD, group = StrainID, col = StrainID),
                method = "loess") +
    facet_grid(SpeciesName ~ CultureTemperature, 
               scales = "free", 
               labeller = labeller(SpeciesName = species.labels,
                                   CultureTemperature = factor.labels))+
    scale_x_continuous(name = "Time (hours)",
                       breaks = seq(from = 0, 
                                    to = max(as.numeric(GC$TimePoint)), 
                                    by = 12)) +
    scale_y_continuous(name = paste("Normalised OD", 
                                    GC$Wavelength[1], sep = "")) + 
    coord_cartesian(xlim = c(0, 60),
                    ylim = NULL) +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(col = "Strain ID") 
  
  GC_bystrain + ggtitle("Growth of Prototheca (or related) strains under temperature stress")
  
  # GC_byinfection <- ggplot(data = GC, 
  #                          aes(x = (as.numeric(TimePoint)-1)/2)) +
  #     geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection),
  #                 method = "loess") +
  #     facet_grid(SpeciesName ~ CultureTemperature, 
  #                scales = "free", 
  #                labeller = labeller(SpeciesName = species.labels,
  #                                    CultureTemperature = factor.labels))+
  #     scale_x_continuous(name = "Time (hours)",
  #                        breaks = seq(from = 0, 
  #                                     to = max(as.numeric(GC$TimePoint)), 
  #                                     by = 12)) +
  #     scale_y_continuous(name = paste("Normalised OD", 
  #                                     GC$Wavelength[1], sep = "")) + 
  #     coord_cartesian(xlim = c(0, 60),
  #                     ylim = NULL) +
  #     theme(strip.text.y = element_text(angle = 0)) +
  #     scale_color_manual(values = c(`TRUE`="red", 
  #                                   `FALSE` = "black", 
  #                                   UNKNOWN = "blue", 
  #                                   SYMBIOSIS = "green")) +
  #     labs(col = "Strain isolated\nfrom infection") 
  # 
  # GC_byinfection +
  #     ggtitle("Growth of Prototheca (or related) strains under temperature stress")
  
  GC_byinfection <- ggplot(data = GC, 
                           aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Source),
                method = "loess") +
    facet_grid(SpeciesName ~ CultureTemperature, 
               scales = "free", 
               labeller = labeller(SpeciesName = species.labels,
                                   CultureTemperature = factor.labels))+
    scale_x_continuous(name = "Time (hours)",
                       breaks = seq(from = 0, 
                                    to = max(as.numeric(GC$TimePoint)), 
                                    by = 12)) +
    scale_y_continuous(name = paste("Normalised OD", 
                                    GC$Wavelength[1], sep = "")) + 
    coord_cartesian(xlim = c(0, 60),
                    ylim = NULL) +
    theme(strip.text.y = element_text(angle = 0)) +
    scale_color_manual(values = c(Infection="red", 
                                  Environment = "black", 
                                  Symbiosis = "green",
                                  Unknown = "blue")) +
    labs(col = "Origin of strain", shape = "Origin of strain") 
  
  GC_byinfection +
    ggtitle("Growth of Prototheca (or related) strains under temperature stress")
  
  #### Calculate error bars by strain ####
  
  # A proper calculation of predictive intervals can use the msir package
  # and is explained here: 
  #https://stats.stackexchange.com/questions/141552/how-to-calculate-prediction-intervals-for-loess
  # However, I can't just plug my data in. I would probably need to subset
  # for each strain, calculate the intervals for each strain, and then 
  # combine them to somehow include in ggplot. 
  
  eight.hour.timepoints <- (0:10*8)*2+1
  
  prediction.bars <-  data.frame(TimePoint = numeric(0), 
                                 NormalisedOD = numeric(0), 
                                 upper = numeric(0), 
                                 lower = numeric(0),
                                 StrainID = character(0), 
                                 CultureTemperature = numeric(0),
                                 Infection = logical(0),
                                 SpeciesName = character(0),
                                 Source = character(0))
  
  # strain <- "HP51"
  # temperature <- 37
  
  for(strain in unique(GC$StrainID)){
    for(temperature in unique(GC$CultureTemperature)){
      strain.temp.data <- subset(GC, StrainID == strain & CultureTemperature == temperature)
      
      if(nrow(strain.temp.data) == 0) next 
      
      l <- loess.sd(x = strain.temp.data$TimePoint, 
                    y = strain.temp.data$NormalisedOD, nsigma = 1.96)
      
      #plot(l$x, l$y, ylim = c(0, 0.4))
      #lines(l$x, l$upper, lty=2)
      #lines(l$x, l$lower, lty=2)
      
      loess.vals <- data.frame(TimePoint = l$x, NormalisedOD = l$y, 
                               upper = l$upper, lower = l$lower,
                               StrainID = strain, CultureTemperature = temperature,
                               Infection = strain.temp.data$Infection[1],
                               SpeciesName = strain.temp.data$SpeciesName[1],
                               Source = strain.temp.data$Source[1])
      
      #ggplot(data = strain.temp.data, aes(x = TimePoint)) +
      #  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
      #  coord_cartesian(ylim = c(0,0.4))
      #
      #ggplot(data = strain.temp.data, aes(x = (as.numeric(TimePoint)-1)/2)) +
      #  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
      #  geom_line(data = loess.vals, aes(y = upper))  +
      #  geom_line(data = loess.vals, aes(y = lower))  +
      #  coord_cartesian(ylim = c(0,0.35))
      
      
      # to add bars, try:
      
      loess.vals <- subset(loess.vals, TimePoint %in% eight.hour.timepoints)
      
      prediction.bars <- rbind(prediction.bars, loess.vals)
    }
  }
  
  #### Plot with prediction bars by strain ####
  ggplot(data = strain.temp.data, aes(x = (as.numeric(TimePoint)-1)/2)) +
    geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Infection)) +
    geom_pointrange(data = prediction.bars, 
                    aes(y = NormalisedOD, ymin = lower, ymax = upper, col = Infection))  +
    coord_cartesian(ylim = c(0,0.35))
  
  
  #https://stackoverflow.com/questions/43956574/using-position-dodge-with-geom-pointrange
  
  GC_bystrain_predbars <- GC_bystrain +
    geom_pointrange(data = prediction.bars, 
                    aes(y = NormalisedOD, ymin = lower, ymax = upper, 
                        col = StrainID, shape = Source),
                    size = 0.1,
                    position = position_dodge(width = 1))
  
  (GC_byinfection_predbars <- GC_byinfection +
      geom_pointrange(data = prediction.bars, 
                      aes(y = NormalisedOD, ymin = lower, ymax = upper, group=StrainID, 
                          col = Source, shape = Source),
                      size = 0.4,
                      position = position_dodge(width = 1))
  )
  
  
  ggsave(filename = paste0("TemperatureGC_platereader_bars_",
                           wavelength, media, ".png"),
         plot = GC_byinfection_predbars,
         device = "png",
         width = 20,
         height = 30,
         units = "cm")



##### Growthcurver #####
  
  ####
  if(media == "YPD") next # skip growthcurver for YPD
  if(wavelength == 600) next # skip growthcurver for OD600

  #### Generate Growthcurver models for all wells ####
  # Growthcurver appears unable to handle curves in which no growth occurred,
  # resulting in absurdly high carrying capacities (k) or growth rates (r)
  # for individual wells
  
  ### Convert data into the format required by growthcurver
  growthdata <- GC[c("TimePoint", "StrainID", "CultureTemperature","Date","NormalisedOD", "Replicate")]
  
  d <- pivot_wider(growthdata, 
                   names_from=c(StrainID, CultureTemperature, Date, Replicate), 
                   values_from = NormalisedOD)
  ### Set the time points to be 0.5 hours
  d$TimePoint <- (as.numeric(d$TimePoint)-1)/2
  
  # names(d)
  ### Create output data frame
  num_analyses <- sum(grepl("^H",names(d)))
  d_gc_a <- data.frame(sample = character(num_analyses),
                       k = numeric(num_analyses),
                       n0  = numeric(num_analyses),
                       r = numeric(num_analyses),
                       t_mid = numeric(num_analyses),
                       t_gen = numeric(num_analyses),
                       auc_l = numeric(num_analyses),
                       auc_e = numeric(num_analyses),
                       sigma = numeric(num_analyses),
                       stringsAsFactors = FALSE)
  
  ### Set the upper time limit
  trim_at_time <- 60
  
  ### Run model fitting for each well
  # par(mfrow = c(4,3))
  j <- 1
  n.unfitted.alldata <- 0
  unfitted.samples <- c()
  for(i in 1:ncol(d)){
    col_name <- names(d)[i]
    if(!grepl("^H", col_name)) next
    
    gc_fit <- SummarizeGrowth(data_t = d[, "TimePoint"], 
                              data_n = d[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    d_gc_a$sample[j] <- col_name
    d_gc_a[j, 2:9] <- c(gc_fit$vals$k,
                        gc_fit$vals$n0,
                        gc_fit$vals$r,
                        gc_fit$vals$t_mid,
                        gc_fit$vals$t_gen,
                        gc_fit$vals$auc_l,
                        gc_fit$vals$auc_e,
                        gc_fit$vals$sigma)
    j <- j+1
    
    if(gc_fit$model[1] == ""){
      unfitted.samples <- c(unfitted.samples, col_name)
      n.unfitted.alldata <- n.unfitted.alldata +1
    }
  }
  
  ## Many replicates fail to have a reasonable model
  
  ## The number that fail to fit at all are:
  n.unfitted.alldata
  
  ### Generate figure to show how model fitting can fail
  ## Some example failures include
  failures<- c(NoFit = "HP54_37_2022-02-11_B", # no model
               PoorFit = "HP31_37_2021-09-07_B", # Sigma = 0.3
               HighCarryingCapacity = "HP4_25_2021-05-28_B", # k = 6.83e7
               HighGrowthRate = "HP33_42_2021-09-12_B", # rate = 40.7
               HighGenerationTime = "HA1_42_2021-09-12_A", # gen time = 285
               NegativeInflectionTime = "HP39_42_2021-09-12_A") # t_mid = -1419.17
  #"HP18_42_2021-08-06_A"
  ## Produce the figure
  png("FittingGrowthcurverModels.png",
      width = 9, height = 6, units = "in", res = 300)
  par(mfrow = c(2,3))
  for(i in 1:length(failures)){
    col_name <- failures[i]
    gc_fit <- SummarizeGrowth(data_t = d[, "TimePoint"], 
                              data_n = d[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    title <- LETTERS[i]
      # paste(paste(strsplit(col_name, split = "_")[[1]][1:2], collapse = " at "), 
      #              "°C \n",
      #              names(failures)[i], 
      #              sep = "")
    
    plot(gc_fit$data$t, gc_fit$data$N, 
      #   main = title,
         xlab = "Time Point", 
         ylab = paste("Normalised OD", GC$Wavelength[1], sep = ""))
    title(main = title, 
          adj = 0,
          cex.main = 2)
    if(gc_fit$vals$note != "cannot fit data"){
      lines(gc_fit$data$t, 
            predict(gc_fit$model), 
            col = "red")
    }
  }
  
  dev.off()
  
  d_gc_a[d_gc_a$sample %in% failures,]
  
  #### Identify which models failed to be fitted ####
  ### Identify carrying capacity threshold for failure
  par(mfrow = c(1,2))
  hist(d_gc_a$k, 
       main = "Histogram of\ncarrying capacities",
       xlab = "Carrying Capacity (k)")
  hist(d_gc_a$k[d_gc_a$k < 200], 
       main = "Histogram of\ncarrying capacities\nunder 200",
       xlab = "Carrying Capacity (k)")
  # export at 900 wide, 600 tall
  
  png("DeterminingThresholds1.png", width = 600, height = 300)
  par(mfrow = c(1,2))
  hist(d_gc_a$k, 
       main = "",
       xlab = "Carrying Capacity (k)")
  title(main = "A", adj = 0, cex.main = 2)
  hist(d_gc_a$k[d_gc_a$k < 200], 
       main = '', # "Histogram of\ncarrying capacities\nunder 200",
       xlab = "Carrying Capacity (k)")
  title(main = "B", adj = 0, cex.main = 2)
   dev.off()
  
  kthreshold <- 30
  
  # number of carrying capacity outliers, with threshold judged by eye
  n.k.failures <- sum(d_gc_a$k>kthreshold) 
  # 1-nrow(d_gc_a[d_gc_a$k>kthreshold,])/nrow(d_gc_a)
  
  sort(d_gc_a$k, decreasing = T)[1:(n.k.failures+3)]
  
  par(mfrow = c(1,1))
  
  ### Identify growth rate threshold for failure
  par(mfrow = c(2,2))
  hist(d_gc_a$r,
       main = "Histogram of\ngrowth rates",
       xlab = "Growth Rate (r)")
  hist(d_gc_a$r[d_gc_a$r<10],
       main = "Histogram of\ngrowth rates\nunder 10",
       xlab = "Growth Rate (r)")
  hist(d_gc_a$r[d_gc_a$r<4],
       main = "Histogram of\ngrowth rates\nunder 4",
       xlab = "Growth Rate (r)")
  # export at 900 wide, 600 tall
  
  png("DeterminingThresholds2.png", width = 600, height = 300)
  par(mfrow = c(1,3))
  hist(d_gc_a$r,
       main = '', #"Histogram of\ngrowth rates",
       xlab = "Growth Rate (r)")
  title(main = "C", adj = 0, cex.main = 2)
  hist(d_gc_a$r[d_gc_a$r<10],
       main = '', #"Histogram of\ngrowth rates\nunder 10",
       xlab = "Growth Rate (r)")
  title(main = "D", adj = 0, cex.main = 2)
  hist(d_gc_a$r[d_gc_a$r<4],
       main = '', #"Histogram of\ngrowth rates\nunder 4",
       xlab = "Growth Rate (r)")
  title(main = "E", adj = 0, cex.main = 2)
  dev.off()
  
  rthreshold <- 1.5
  
  # number of growth rate outliers, as judged by eye
  n.r.failures <- sum(d_gc_a$r>rthreshold) 
  #1-nrow(d_gc_a[d_gc_a$r>rthreshold,])/nrow(d_gc_a)
  
  sort(d_gc_a$r, decreasing = T)[1:(n.r.failures+3)]
  par(mfrow = c(1,1))
  
 
  
  ### Identify inflection point threshold for failure
  # inflection point is half way to the carrying capacity, so should be positive
  par(mfrow = c(1,2))
  hist(d_gc_a$t_mid,
       main = "Histogram of\ninflection points",
       xlab = "Inflection point (t_mid)")
  
  tmidthreshold <- 0
  
  n.tmid.failures <- sum(d_gc_a$t_mid<tmidthreshold) 
  sort(d_gc_a$t_mid)[1:(n.tmid.failures+3)]
  #head(d_gc_a[d_gc_a$t_mid<tmidthreshold,])
  
  par(mfrow = c(3,3))
  tmid.samples <- d_gc_a$sample[d_gc_a$t_mid < tmidthreshold]
  
  tmid.samples <- tmid.samples[order(d_gc_a$t_mid[d_gc_a$t_mid < tmidthreshold])]
  for(i in tmid.samples){
    col_name <- i
    gc_fit <- SummarizeGrowth(data_t = d[, "TimePoint"], 
                              data_n = d[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    title <- paste(paste(strsplit(col_name, split = "_")[[1]][c(1,2,4)], collapse = " at "), 
                   "°C\n", round(d_gc_a$t_mid[d_gc_a$sample == i], digits = 3),
                   sep = "")
    
    plot(gc_fit$data$t, gc_fit$data$N, 
         main = title,
         xlab = "Time Point", 
         ylab = paste("Normalised OD", GC$Wavelength[1], sep = ""))
    if(gc_fit$vals$note != "cannot fit data"){
      lines(gc_fit$data$t, 
            predict(gc_fit$model), 
            col = "red")
    }
  }
  
  par(mfrow = c(1,1))
  
  ### Identify sigma for failure
  
  png("SigmaThreshold.png",
      width = 600, height = 600)
  
  par(mfrow = c(3,3))
  hist(d_gc_a$sigma)
  
  for(i in d_gc_a$sample[d_gc_a$sigma > 0.15]){
    col_name <- i
    gc_fit <- SummarizeGrowth(data_t = d[, "TimePoint"], 
                              data_n = d[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    title <- paste(paste(strsplit(col_name, split = "_")[[1]][c(1,2,4)], collapse = " at "), 
                   "°C\nSigma = ", round(d_gc_a$sigma[d_gc_a$sample == i], digits = 3),
                   sep = "")
    
    plot(gc_fit$data$t, gc_fit$data$N, 
         main = title,
         xlab = "Time Point", 
         ylab = paste("Normalised OD", GC$Wavelength[1], sep = ""))
    if(gc_fit$vals$note != "cannot fit data"){
      lines(gc_fit$data$t, 
            predict(gc_fit$model), 
            col = "red")
    }
  }
  dev.off()
  
  png("SigmaThreshold1.png",
      width = 300, height = 200)
  
  hist(d_gc_a$sigma,
       main = "",
       xlab = "Sigma values")
  title(main = "A", adj = 0, cex.main = 2)
  dev.off()
  
  png("SigmaThreshold2.png",
      width = 800, height = 400)
  par(mfrow = c(2,4))

  j <- 2
  for(i in d_gc_a$sample[d_gc_a$sigma > 0.15]){
    col_name <- i
    gc_fit <- SummarizeGrowth(data_t = d[, "TimePoint"], 
                              data_n = d[, col_name],
                              t_trim = trim_at_time,
                              bg_correct = "none")
    
    
    
    plot(gc_fit$data$t, gc_fit$data$N, 
         main = '',
         xlab = "Time Point", 
         ylab = paste("Normalised OD", GC$Wavelength[1], sep = ""))
    title(main = LETTERS[j], adj = 0, cex.main = 2)
    j <- j+1
    if(gc_fit$vals$note != "cannot fit data"){
      lines(gc_fit$data$t, 
            predict(gc_fit$model), 
            col = "red")
    }
  }
  dev.off()
  paste(LETTERS[2:8], round(d_gc_a$sigma[d_gc_a$sigma > 0.15], digits = 3), sep = " = ", collapse = ", ")
  
  sigmathreshold <- 0.2
  n.sigma.failures <- sum(d_gc_a$sigma>sigmathreshold) 
  
  ### calculate number of total failures 
  
  failed.models.rows <- d_gc_a$r>rthreshold|
    d_gc_a$k>kthreshold|
    #  d_gc_a$t_mid<tmidthreshold| # I can't justify excluding models based on tmid
    d_gc_a$sigma>sigmathreshold|
    d_gc_a$n0==0
  
  failed.models <- d_gc_a[failed.models.rows,]
  
  n.tot.failures <- sum(failed.models.rows) 
  
  ### 
  summary.of.failures <- data.frame(Statistic = c("Carrying Capacity",
                                                  "Growth Rate", 
                                                  #    "Time of Inflection Point",
                                                  "Sigma",
                                                  "Failure to Fit a Model",
                                                  "All"),
                                    Threshold = c(kthreshold, rthreshold, 
                                                  #tmidthreshold, 
                                                  sigmathreshold,
                                                  NA, "All"),
                                    N.Failures = 
                                      c(n.k.failures,
                                        n.r.failures,
                                        #    n.tmid.failures,
                                        n.sigma.failures,
                                        n.unfitted.alldata,
                                        n.tot.failures))
  
  summary.of.failures$PercentFailed <- round(summary.of.failures$N.Failures/nrow(d_gc_a)*100,
                                             digits = 2)
  
  summary.of.failures
  
  ### Identify which models fail, and how many fail
  #### Split sample names so each variable has its own column
  models <- do.call(rbind, strsplit(d_gc_a$sample, split = "_"))
  models <- as.data.frame(models)
  names(models) <- c("Strain", "Temperature", "Date", "Replicate")
  
  models <- cbind(models, d_gc_a)
  
  #### Provide species names
  models <- merge(models, strain.to.species.df)
  models$Species <- factor(models$Species,
                           levels = species.levels$Descent)
  
  #### Group models according to factors of interest: stress and species
  models <- group_by(models, Species, Temperature)
  
  #### Filter for models that meet thresholds set above
  valid.models <- models[models$r <= rthreshold &
                           models$k <= kthreshold &
                           #  models$t_mid >= tmidthreshold &
                           models$sigma <= sigmathreshold &
                           models$n0 != 0,] 
  
  #### Determine the number of models from each species that fail each threshold
  summary.models <- summarise(models, 
                              n.wells.total = length(sample),
                              k.failed = sum(k>kthreshold),
                              r.failed = sum(r>rthreshold),
                              # tmid.failed = sum(t_mid<tmidthreshold),
                              sigma.failed = sum(sigma>sigmathreshold),
                              model.failed = sum(n0==0))
  
  #### Determine the number of models that remain after removing outlier models
  summary.valid.models <- summarise(valid.models, 
                                    n.wells.valid = length(sample))
  
  #### Combine the summary statistics
  calculate.valid.models <- left_join(summary.models, summary.valid.models, 
                                      by = c("Species", "Temperature"))
  
  #### Remove any NAs, which can emerge if there are no samples of a given species/stress level
  #### (might happen if all models are outliers, thus that species/stress level doesn't exist
  #### in the summary.valid.models data frame)
  calculate.valid.models$n.wells.valid[is.na(calculate.valid.models$n.wells.valid)] <- 0
  
  #### Calculate proportions of models that fail
  calculate.valid.models$n.wells.invalid <- 
    calculate.valid.models$n.wells.total - 
    calculate.valid.models$n.wells.valid
  
  calculate.valid.models$prop.wells.invalid <-
    round(calculate.valid.models$n.wells.invalid/
            calculate.valid.models$n.wells.total,
          digits = 4)
  sum(calculate.valid.models$n.wells.invalid)
  
  #### Filter the summary, to only show species where models have been removed
  problematic.models <- calculate.valid.models[calculate.valid.models$n.wells.invalid > 0,]
  
  
  problematic.models$Species <- shorten.species.name(problematic.models$Species)
  #### Summarise models that do not fit ####
  summary.of.failures
  100-summary.of.failures[nrow(summary.of.failures), "PercentFailed"]
  problematic.models <- as.data.frame(problematic.models)
  problematic.models <- arrange(problematic.models, Temperature, Species)
  
  subset(problematic.models, Temperature == 25)
  subset(problematic.models, Temperature == 30)
  subset(problematic.models, Temperature == 37)
  subset(problematic.models, Temperature == 42)
  subset(problematic.models, n.wells.valid == 0)
  
  ### Export summaries
  write.csv(summary.of.failures,
            "SummaryOfFailuresToFitModels.csv",
            row.names = FALSE)
  
  write.csv(problematic.models,
            "BreakdownOfFailuresToFitModels.csv",
            row.names = FALSE)
  
  
  #### Subset models using thresholds ####
  d_gc <- subset(d_gc_a,
                 k <= kthreshold &
                   r <= rthreshold &
                   t_mid >= tmidthreshold &
                   sigma <= sigmathreshold &
                   n0 != 0)
  #### Identify lost species ####
  ## In removing data, we might lose too much data from some strains to do
  ## statistical comparisons. Check to see if there are any temperatures
  ## with only one well's worth of data
  sample.names <- do.call(rbind,strsplit(d_gc$sample, split = "_"))
  sample.names <- as.data.frame(sample.names)
  names(sample.names) <- c("Strain",
                           "Temperature",
                           "Date",
                           "Replicate")
  
  ### Identify number of acceptable technical replicates on a given date
  n.wells.remaining <- table(sample.names[,c("Strain", "Date", "Temperature")])
  n.wells.remaining <- as.data.frame(n.wells.remaining)
  n.wells.remaining <- n.wells.remaining[n.wells.remaining$Freq != 0,]
  
  ### Identify number of dates with at least one technical replicate (i.e.
  ### number of biological replicates)
  l <- length(unique(n.wells.remaining$Strain))
  n.dates.remaining <- data.frame(Strain = character(length = l),
                                  nDates30 = numeric(length = l),
                                  nDates37 = numeric(length = l))
  
  for(i in 1:length(unique(n.wells.remaining$Strain))){
    s <- unique(n.wells.remaining$Strain)[i]
    sub <- subset(n.wells.remaining, Strain == s)
    n.dates.remaining$Strain[i] <- as.character(s)
    n.dates.remaining$nDates30[i] <- sum(sub$Temperature==30)
    n.dates.remaining$nDates37[i] <- sum(sub$Temperature==37)
  }
  
  ### Identify strains with fewer than 2 biological replicates under relevant conditions
  strains.lost <- n.dates.remaining$nDates30 < 2 |
    n.dates.remaining$nDates37 < 2
  
  lost.strains <- n.dates.remaining[strains.lost,]
  lost.strains <- merge(lost.strains,strain.to.species.df)
  lost.strains <- lost.strains[c("Strain", 'Species', "nDates30", "nDates37")]
  lost.strains
  
  ### Identify species with fewer than 2 biological replicates under relevant conditions
  lost.species <- n.dates.remaining %>%
    left_join(strain.to.species.df) %>%
    group_by(Species) %>%
    summarise(nMeasurements30 = sum(nDates30), 
              nMeasurements37 = sum(nDates37), 
              Strains = list(Strain))
  #### Show lost species #### FINDME
  lost.species[lost.species$nMeasurements30<2 |
                 lost.species$nMeasurements37<2,]
  
  paste(shorten.species.name(lost.species[lost.species$nMeasurements30<2 |
                                            lost.species$nMeasurements37<2,]$Species), collapse = ", ")
  #gridExtra::grid.table(lost.strains,
  #                      rows = NULL)
  
  #### Remove lost species ####
  ### remove strains where a given species doesn't have at least two observations
  strains.to.remove <- unlist(lost.species[lost.species$nMeasurements30<2 |
                                             lost.species$nMeasurements37<2,]$Strains)
  strains.to.remove <- paste(strains.to.remove,"_", sep = "")
  
  models.to.remove <- vector("logical", length = nrow(d_gc))
  for(i in 1:length(strains.to.remove)){
    models.to.remove <- grepl(strains.to.remove[i], d_gc$sample) | models.to.remove
  }
  
  d_gc <- d_gc[!models.to.remove,]
  
  
  #### Combine technical replicates ####
  d_gc$sample <- substring(d_gc$sample,1,nchar(d_gc$sample)-2)
  d_gc <- summarise_all(group_by(d_gc, sample), mean)
  
  ### How many data points do I have left?
  
  dp <- strsplit(d_gc$sample, split = "_")
  
  data.points <- do.call(rbind, dp)
  data.points <- as.data.frame(data.points)
  names(data.points) <- c("Strain", "Temperature", "Date")
  data.points <- merge(data.points, strain.to.species.df)
  data.points.summary <- group_by(data.points, Species, Temperature) %>%
    summarise(n = length(Date))
  
  data.points.summary <- pivot_wider(data.points.summary, names_from = Temperature, values_from = n)
  names(data.points.summary)[-1] <- paste0(names(data.points.summary)[-1], "°C")
  
  data.points.summary[is.na(data.points.summary) ]<-0
  
  ### Order the resulting table
  data.points.summary$Species <- factor(data.points.summary$Species, levels = species.levels$Descent)
  data.points.summary <- data.points.summary[order(data.points.summary$Species),]
  data.points.summary$Species <- shorten.species.name(data.points.summary$Species)
  
  ### Export table
  write.csv(data.points.summary,
            file = "NumberOfDataPointsFromConditions.csv",
            row.names = F)
  
  #### Plot growth rate against carrying capacity ####
  ### Set d_gc to be a data frame, so new columns can be added
  d_gc <- as.data.frame(d_gc)
  
  ### Split sample names to provide information on strain, stress, and date
  for(i in 1:nrow(d_gc)){
    s <- strsplit(d_gc$sample[i], split = "_")
    d_gc$Strain[i] <- s[[1]][1]
    d_gc$Temperature[i] <- s[[1]][2]
    d_gc$Date[i] <-s[[1]][3]
  }
  
  ### Assign species names to strain names, again
  d_gc <- merge(d_gc,strain.to.species.df)
  
  
  scatter.KvsR.str <- ggplot(data = d_gc, aes(x = r, y = k)) +
    geom_point(aes(colour = factor(Strain, levels = strain.levels),
                   shape = Temperature)) +
    labs(x = "Growth Rate",
         y = "Carrying Capacity") +
    guides(shape = guide_legend(order = 1),
           colour=guide_legend(title = "Strain",
                               order = 2))
  
  scatter.KvsR.str + ggtitle("Scatter plot of carrying capacity against growth rate") 
  
  ggsave(paste0("ScatterplotOfKvsRStrain_",
                wavelength, media, ".png"),
         scatter.KvsR.str,
         height = 6,
         width = 9)
  
  scatter.KvsR.spp <- ggplot(data = d_gc, aes(x = r, y = k)) +
    geom_point(aes(colour = factor(Species, levels = species.levels$Descent),
                   shape = Temperature)) +
    labs(x = "Growth Rate",
         y = "Carrying Capacity") +
    guides(shape = guide_legend(order = 1),
           colour=guide_legend(title = "Species",
                               order = 2))
  
  scatter.KvsR.spp + ggtitle("Scatter plot of carrying capacity against growth rate") 
  
  ggsave(paste0("ScatterplotOfKvsRSpecies_",
                wavelength, media, ".png"),
         scatter.KvsR.spp,
         height = 6,
         width = 9)
  
  ggplot(data = d_gc_a, aes(x = r, y = k)) +
    geom_point() +
    #scale_x_log10()+
    scale_y_log10()+
    ggtitle("Scatter plot of log carrying capacity against log growth rate") +
    labs(x = "Growth Rate",
         y = "Carrying Capacity")# export h = 500, w = 900
  #### Statistical comparisons of growthcurver models within species ####
  species.present <- unique(d_gc$Species)

  ## Checking carrying capacity
  for(i in 1:length(species.present)){
    d_gc_strain <- subset(d_gc, Species == species.present[i])
    print(species.present[i])
    if(length(unique(d_gc_strain$Strain))== 1){
      print("only 1 strain")
      next
    }
    k <- aov(k ~ Temperature*Strain, data = d_gc_strain)
    ks <- summary(k)

    print(ks)
  }
  
  d_gc_strain <- subset(d_gc, Species == "Prototheca bovis")
  k <- aov(k ~ Temperature*Strain, data = d_gc_strain)
  thsd_k_strain <- TukeyHSD(k)
  sig.dif.strains <- (thsd_k_strain$Strain)[,4]
  sig.dif.strains <- sig.dif.strains[sig.dif.strains<0.05]
  significant.comparisons <- table(unlist(strsplit(names(sig.dif.strains), split = "-")))
  n.strains <- length(unique(d_gc_strain$Strain))
  different.strains.k <- names(significant.comparisons)[significant.comparisons >= n.strains/2]
  
 
  
  ## Checking growth rate
  dif.species <- c()
  different.strains.r <- c()
  for(i in 1:length(species.present)){
    d_gc_strain <- subset(d_gc, Species == species.present[i])
    print(species.present[i])
    if(length(unique(d_gc_strain$Strain))== 1){
      print("only 1 strain")
      next
    }
    r <- aov(r ~ Temperature*Strain, data = d_gc_strain)
    rs <- summary(r)
    
    print(rs)
    
    if(rs[[1]][[5]][2] < 0.05){
      n.strains <- length(unique(d_gc_strain$Strain))
      dif.species <- c(dif.species, species.present[i])
    #  break
      if(n.strains > 2){
      thsd_r_strain <- TukeyHSD(r)
      sig.dif.strains <- (thsd_r_strain$Strain)[,4]
      sig.dif.strains <- sig.dif.strains[sig.dif.strains<0.05]
      significant.comparisons <- table(unlist(strsplit(names(sig.dif.strains), split = "-")))
      different.strains.r <- c(different.strains.r,
                               names(significant.comparisons)[significant.comparisons 
                                                              >= n.strains/2])
      }
    }
  }
  
  ## Check whether clinical ciferrii strains grow faster/slower than environmental
  d <- left_join(subset(d_gc, Species == "Prototheca ciferrii"), 
                 strains, by = c("Strain" = "H_number"))
  d %>% group_by(Strain, From_infection) %>% summarise(Mean = mean(r)) %>% arrange(Mean)
  ## They do not
  
  
  ## Remove non-representative strains
  d_gc <- subset(d_gc, !(Strain %in% different.strains.k))
  
  #d_gc <- subset(d_gc, !(Strain %in% different.strains.r)) 
  
  ## Removing strains based on different growth rates would remove almost all
  ## ciferrii data.
  
  #### Statistical comparisons of growthcurver models between species ####
  ### anova of species data, carrying capacity ###
  k <- aov(k ~ Temperature*Species, data = d_gc)
  summary(k)
  #gridExtra::grid.table(round(as.data.frame(summary(k)[[1]]),digits = 4)) 
  # Export h = 150, w = 500
  
  comparisons_k <- list()
  for(i in 1:length(unique(d_gc$Species))){
    species_ <- unique(d_gc$Species)[i]
    data <- subset(d_gc, Species == species_)
    x <- data$k[data$Temperature == 30]
    y <- data$k[data$Temperature == 37]
    if(length(x)>1 & length(y)>1) {
      comparisons_k[species_] <- list(t.test(x,y))
    } else{
      comparisons_k[species_] <- "insufficient data"
    }
  }
  
  results_k <- data.frame(species = character(length(comparisons_k)),
                          pval = numeric(length(comparisons_k)))
  
  for(i in 1:length(comparisons_k)){
    results_k$species[i] <- names(comparisons_k)[i]
    results_k$pval[i] <- comparisons_k[[i]]$p.value # cannot pull values that don't exist
    results_k$mean30[i] <- comparisons_k[[i]]$estimate["mean of x"]
    results_k$mean37[i] <- comparisons_k[[i]]$estimate["mean of y"]
  }
  
  #results_k$padj_Bon <- p.adjust(results_k$pval, method = "bonferroni", n = length(results_k$pval)) 
  results_k$padj_BH  <- p.adjust(results_k$pval, method = "BH", n = length(results_k$pval)) 
  results_k$difference <- results_k$mean30 - results_k$mean37
  results_k$prop.difference <- results_k$difference/results_k$mean30
  
  results_k <- relocate(results_k, species, mean30, mean37, difference, prop.difference)
  
  ### Reorder results
  results_k$species <- factor(results_k$species, levels = species.levels$Descent)
  results_k <- results_k[order(results_k$species),]
  
  results_k$species <- shorten.species.name(results_k$species)
  #sort(results_k[results_k$padj_Bon < 0.05,]$species)
  names(results_k) <- c("Species",
                        "Mean Carrying Capacity at 30°C",
                        "Mean Carrying Capacity at 37°C",
                        "Difference between means",
                        "Proportional difference",
                        "P value",
                        "Adjusted P value")
  
  results_k$Significance <- "NS"
  results_k$Significance[results_k$`Adjusted P value` <= 0.05] <- "*"
  results_k$Significance[results_k$`Adjusted P value` <= 0.01] <- "**"
  results_k$Significance[results_k$`Adjusted P value` <= 0.001] <- "***"
  
  write.csv(cbind(results_k[1],
                  signif(results_k[-c(1, ncol(results_k))], digits = 3),
                  results_k[ncol(results_k)]),
            "CarryingCapacityStats.csv",
            row.names = FALSE)
  
  ### anova of species data, growth rate ###
  r <- aov(r ~ Temperature*Species, data = d_gc)
  summary(r)
  # gridExtra::grid.table(round(as.data.frame(summary(r)[[1]]),digits = 4)) 
  # Export w = 500, h = 150
  
  comparisons_r <- list()
  for(i in 1:length(unique(d_gc$Species))){
    species_ <- unique(d_gc$Species)[i]
    data <- subset(d_gc, Species == species_)
    x <- data$r[data$Temperature == 30]
    y <- data$r[data$Temperature == 37]
    if(length(x)>1 & length(y)>1) {
      comparisons_r[species_] <- list(t.test(x,y))
    } else{
      comparisons_r[species_] <- "insufficient data"
      
    }
  }
  
  results_r <- data.frame(species = character(length(comparisons_r)),
                          pval = numeric(length(comparisons_r)))
  
  for(i in 1:length(comparisons_r)){
    results_r$species[i] <- names(comparisons_r)[i]
    results_r$pval[i] <- comparisons_r[[i]]$p.value # cannot pull values that don't exist
    results_r$mean30[i] <- comparisons_r[[i]]$estimate["mean of x"]
    results_r$mean37[i] <- comparisons_r[[i]]$estimate["mean of y"]
  }
  
  results_r$difference <- results_r$mean30 - results_r$mean37
  # results_r$padj_Bon <- p.adjust(results_r$pval, method = "bonferroni", n = length(results_r$pval)) 
  results_r$padj_BH  <- p.adjust(results_r$pval, method = "BH", n = length(results_r$pval)) 
  #sort(results_r[results_r$padj_Bon < 0.05,]$species)
  
  sort(results_r[results_r$padj < 0.05,]$species)
  
  results_r <- relocate(results_r, species, mean30, mean37, difference)
  
  ### Reorder results
  results_r$species <- factor(results_r$species, levels = species.levels$Descent)
  results_r <- results_r[order(results_r$species),]
  
  results_r$species <- shorten.species.name(results_r$species)
  
  #sort(results_k[results_k$padj_Bon < 0.05,]$species)
  names(results_r) <- c("Species",
                        "Mean Growth Rate at 30°C",
                        "Mean Growth Rate at 37°C",
                        "Difference between means",
                        "P value",
                        "Adjusted P value")
  
  results_r$Significance <- "NS"
  results_r$Significance[results_r$`Adjusted P value` <= 0.05] <- "*"
  results_r$Significance[results_r$`Adjusted P value` <= 0.01] <- "**"
  results_r$Significance[results_r$`Adjusted P value` <= 0.001] <- "***"
  
  write.csv(cbind(results_r[1],
                  signif(results_r[-c(1, ncol(results_r))], digits = 3),
                  results_r[ncol(results_r)]),
            "GrowthRateStats.csv",
            row.names = FALSE)
  #comparisons_r["Prototheca cookei"]
  
  #### Results of growthcurver statistical comparisons ####
  # Which strains were different from other strains in the same species?
  different.strains.k
  different.strains.r
  
  data.points.summary
  
  results_k
  results_r
  # comparisons_k["Prototheca cookei"]
  
  #### Compare start and end ODs ####
  ## Compare individual wells
  ### Subset to only have 37°C data from desired time points
  se.data <- subset(GC, TimePoint %in% c("0 h 0 min", "60 h 0 min"))
  se.data <- subset(se.data, CultureTemperature == 37)
  se.data <- subset(se.data, StrainID != "HP18") # removed due to different carrying capacity
  
  ### Create single identifier for wells
  se.data <- mutate(se.data, 
                    WellID = paste(Date, PlateNumber, xy, sep = "_"))
  
  ### Discard unnecessary columns
  se.data <- select(se.data, StrainID, SpeciesName, 
                    WellID, 
                    NormalisedOD, TimePoint)
  
  ### Split the data into two columns, for the T tests
  se.data <- pivot_wider(se.data, names_from = TimePoint, values_from = NormalisedOD)
  
  ## Compare strains between plates
  #se.data <- subset(GC, TimePoint %in% c("0 h 0 min", "60 h 0 min"))
  #se.data <- subset(se.data, CultureTemperature == 37)
  #se.data <- select(se.data, StrainID, SpeciesName, Date, NormalisedOD, TimePoint)
  #se.data <- group_by(se.data, StrainID, Date, TimePoint)
  #se.data <- unique(summarise(se.data, MeanOD = mean(NormalisedOD), SpeciesName = SpeciesName))
  #se.data <- pivot_wider(se.data, names_from = TimePoint, values_from = MeanOD)
  
  
  comparisons_se <- list()
  for(i in 1:length(unique(se.data$SpeciesName))){
    species_ <- as.character(unique(se.data$SpeciesName)[i])
    data_ <- subset(se.data, SpeciesName == species_)
    comparisons_se[species_] <- list(t.test(data_$`60 h 0 min`, data_$`0 h 0 min`, paired = T))
  }
  
  results_se <- data.frame(species = character(length(comparisons_se)),
                           growth = numeric(length(comparisons_se)),
                           conf.int.upper = numeric(length(comparisons_se)),
                           conf.int.lower = numeric(length(comparisons_se)),
                           pval = numeric(length(comparisons_se)))
  
  for(i in 1:length(comparisons_se)){
    results_se$species[i] <- names(comparisons_se)[i]
    results_se$growth[i] <- comparisons_se[[i]]$estimate
    results_se$conf.int.upper[i] <- comparisons_se[[i]]$conf.int[2]
    results_se$conf.int.lower[i] <- comparisons_se[[i]]$conf.int[1]
    results_se$pval[i] <- comparisons_se[[i]]$p.value # cannot pull values that don't exist
  }
  
  
  results_se$padj_BH  <- p.adjust(results_se$pval, method = "BH", n = length(results_se$pval)) 
  
  (results_se <- cbind(results_se[1],signif(results_se[2:6], digits = 3)))
  
  names(results_se) <- c("Species",
                         "Growth at 37°C",
                         "Upper 95% Confidence Limit",
                         "Lower 95% Confidence Limit",
                         "P value",
                         "Adjusted P value")
  
  results_se$Species <- factor(results_se$Species, levels = species.levels$Descent)
  results_se <- results_se[order(results_se$Species),]
  results_se$Species <- shorten.species.name(results_se$Species)
  
  results_se$Significance <- "NS"
  results_se$Significance[results_se$`Adjusted P value` <= 0.05] <- "*"
  results_se$Significance[results_se$`Adjusted P value` <= 0.01] <- "**"
  results_se$Significance[results_se$`Adjusted P value` <= 0.001] <- "***"
  
  write.csv(results_se,
            "PairedComparisonsBetweenTimepoints.csv",
            row.names = FALSE)
  # gridExtra::grid.table(results_se, rows = NULL) 
  # export w = 800, h = 360
  #####
}

##### Skin temperature plots #####

#### Load data ####
### Choose dataset
data.files <- c("SAB_vardeg_pH5.5_750.csv", "SAB_skintemp_pH5.5_750.csv")
# data.files <- gsub("750", "600", data.files)

### Load data
GC <- rbind(
  read.csv(paste(plate.reader.data, data.files[1], sep = "")),
  read.csv(paste(plate.reader.data, data.files[2], sep = "")))

### Remove any entries which were shaken
GC <- subset(GC, MeasurementShaking == 0)

### Remove data from HA6, as it doesn't grow in the plate reader
GC <- subset(GC, StrainID != "HA6")

#### Process Data ####
### Take constants from GC data 
wavelength <- unique(GC$Wavelength)
media <- unique(GC$CultureMedia)




### Count the number of replicates per species
species_replicates <- GC %>%
  group_by(StrainID, CultureTemperature, CultureMedia, Wavelength) %>%
  summarise(N = length(unique(Date))) %>%
  left_join(strain.to.species.df, by = c("StrainID" = "Strain")) %>%
  group_by(Species, CultureTemperature, CultureMedia, Wavelength) %>%
  summarise(Replicates = sum(N)) %>%
  pivot_wider(names_from = CultureTemperature, values_from = Replicates)


names(species_replicates)[4:ncol(species_replicates)] <- paste(names(species_replicates)[4:ncol(species_replicates)], 
                                                               "degC", sep = "")

species_replicates[is.na(species_replicates)] <- 0

write.csv(species_replicates,
          paste("replicates_per_species_",
                "OD", wavelength, "_",
                media, 
                "_skintemp.csv",
                sep = ""),
          row.names = FALSE)




### 
GC$Infection <- NA
GC$Host <- NA
for(i in 1:nrow(strains)){
  strain.ID <- strains$H_number[i]
  infect <- strains$From_infection[i]
  host <- strains$Host[i]
  species <- strains$Species[i]
  
  rows <- which(GC$StrainID == strain.ID)
  GC$Infection[rows] <- infect
  GC$Host[rows]      <- host
  GC$SpeciesName[rows]<- species
}

GC$TimePoint <- factor(GC$TimePoint, levels = time.point.levels)


#### Set grid labels ####
factor.labels <- paste(unique(GC$CultureTemperature), "°C", sep = "")
names(factor.labels) <- unique(GC$CultureTemperature)
factor.labels[which(names(factor.labels) == 25)] <- "Room\nTemperature"

#### Calculate prediction intervals ####
eight.hour.timepoints <- (0:10*8)*2+1

prediction.bars <-  data.frame(TimePoint = numeric(0), 
                               NormalisedOD = numeric(0), 
                               upper = numeric(0), 
                               lower = numeric(0),
                               SpeciesName = character(0), 
                               CultureTemperature = numeric(0),
                               Infection = logical(0),
                               Source = character(0))

# strain <- "HP51"
# temperature <- 37

for(species in unique(GC$SpeciesName)){
  for(temperature in unique(GC$CultureTemperature)){
    species.temp.data <- subset(GC, SpeciesName == species & CultureTemperature == temperature)
    
    if(nrow(species.temp.data) == 0) next 
    
    l <- loess.sd(x = species.temp.data$TimePoint, 
                  y = species.temp.data$NormalisedOD, nsigma = 1.96)
    
    loess.vals <- data.frame(TimePoint = l$x, NormalisedOD = l$y, 
                             upper = l$upper, lower = l$lower,
                             SpeciesName = species, CultureTemperature = temperature,
                             Infection = strain.temp.data$Infection[1],
                             Source = strain.temp.data$Source[1])
    
    loess.vals <- subset(loess.vals, TimePoint %in% eight.hour.timepoints)
    
    prediction.bars <- rbind(prediction.bars, loess.vals)
  }
}
#### Set species order ####
GC$SpeciesName <- factor(GC$SpeciesName, 
                         levels = species.levels$Descent) #Descent, Alphabet, or Pathology

prediction.bars$SpeciesName <- factor(prediction.bars$SpeciesName,
                                      levels = species.levels$Descent)

#### Make Temperature plot by species ####

GC_byspecies<- ggplot(data = GC, 
                     aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = SpeciesName, col = SpeciesName),
              method = "loess") +
  facet_grid(SpeciesName ~ CultureTemperature, 
             scales = "free", 
             labeller = labeller(SpeciesName = species.labels,
                                 CultureTemperature = factor.labels))+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", 
                                  GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 60),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(col = "Species") 

# GC_byspecies + ggtitle("Growth of Prototheca (or related) species under temperature stress")

GC_byspecies_predbars <- GC_byspecies +
  geom_pointrange(data = prediction.bars, 
                  aes(y = NormalisedOD, ymin = lower, ymax = upper, 
                      col = SpeciesName),
                  size = 0.1)

# GC_byspecies_predbars + ggtitle("Growth of Prototheca (or related) species under temperature stress")

ggsave(filename = paste0("SkinTemperatureGC_platereader_bars_",
                         wavelength, media, ".png"),
       plot = GC_byspecies_predbars,
       device = "png",
       width = 22,
       height = 30,
       units = "cm")

#####

##### Produce graphs for powerpoint ####
require(ggpubr)

f.half <- sort(unique(GC$SpeciesName))[1:(length(unique(GC$SpeciesName))/2)]
s.half <- base::setdiff( sort(unique(GC$SpeciesName)),f.half)

GC$Source <- factor(GC$Source,
                    levels = c('Infection', 
                               'Environment',
                               'Symbiosis' ,
                               'Unknown'))

left <- ggplot(data = subset(GC, SpeciesName %in% f.half), aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Source), method = "loess")+
  geom_pointrange(data = subset(prediction.bars, SpeciesName %in% f.half), 
                  aes(y = NormalisedOD, ymin = lower, ymax = upper, group=StrainID, col = Source),
                  size = 0.1,
                  position = position_dodge(width = 1)) +
  facet_grid(SpeciesName ~ CultureTemperature, 
             scales = "free", 
             labeller = labeller(SpeciesName = species.labels,
                                 CultureTemperature = factor.labels))+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", 
                                  GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 60),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0, face = "italic")) +
  scale_color_manual(values = c(Infection="red", 
                                Environment = "black", 
                                Symbiosis = "green",
                                Unknown = "blue"),
                     drop = FALSE) +
  labs(col = "Origin of strain")

right <- ggplot(data = subset(GC, SpeciesName %in% s.half), aes(x = (as.numeric(TimePoint)-1)/2)) +
  geom_smooth(aes(y = NormalisedOD, group = StrainID, col = Source), method = "loess") +
  geom_pointrange(data = subset(prediction.bars, SpeciesName %in% s.half), 
                  aes(y = NormalisedOD, ymin = lower, ymax = upper, group=StrainID, col = Source),
                  size = 0.1,
                  position = position_dodge(width = 1))+
  facet_grid(SpeciesName ~ CultureTemperature, 
             scales = "free", 
             labeller = labeller(SpeciesName = species.labels,
                                 CultureTemperature = factor.labels))+
  scale_x_continuous(name = "Time (hours)",
                     breaks = seq(from = 0, 
                                  to = max(as.numeric(GC$TimePoint)), 
                                  by = 12)) +
  scale_y_continuous(name = paste("Normalised OD", 
                                  GC$Wavelength[1], sep = "")) + 
  coord_cartesian(xlim = c(0, 60),
                  ylim = NULL) +
  theme(strip.text.y = element_text(angle = 0,
                                    face = "italic")) +
  scale_color_manual(values = c(Infection="red", 
                                Environment = "black", 
                                Symbiosis = "green",
                                Unknown = "blue"),
                     drop = FALSE) +
  labs(col = "Origin of strain")

combo <- ggarrange(left, right,
                   common.legend = TRUE, legend = "right")

ggsave("~/Documents/Main_Project/Admin/Presentations/PostdocDylan/TempGrowthCurves.png",
       plot = combo,
       width = 36,
       height = 20,
       units = "cm")

#####