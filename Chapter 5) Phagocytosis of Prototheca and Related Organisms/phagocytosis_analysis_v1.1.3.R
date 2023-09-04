##### Setup #####
#### Packages ####
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(lubridate)
#### Functions ####
shorten.species.name <- function(s){
  split.name <- strsplit(s, split = " ")[[1]]
  genus.shortened <- paste0(substring(split.name[1],1,1), ".")
  short.name <- paste(genus.shortened, split.name[2])
  return(short.name)
}

std.err <- function(x, na.rm = TRUE) {
  if(na.rm){
   x <- na.omit(x) 
  }
  sd(x) / sqrt(length(x))
  }

significance.stars <- function(v){
  o <- rep("NS", length(v))
  o[v <= 0.05] <- "*"
  o[v <= 0.01] <- "**"
  o[v <= 0.001] <- "***"
  return(o)
}

#### Misc ####

options(scipen=-3)
strain_numbers <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")

pathogenic_species <- strain_numbers[2:3]
pathogenic_species <- pathogenic_species[!is.na(pathogenic_species$From_infection),]
pathogenic_species$From_infection[pathogenic_species$From_infection == "SYMBIOSIS"] <- TRUE
pathogenic_species$From_infection[pathogenic_species$From_infection == "UNKNOWN"] <- FALSE
pathogenic_species <- pathogenic_species %>%
  group_by(Species) %>%
  summarise(Infectious = any(as.logical(From_infection)))

pathogenic_rows <- pathogenic_species$Infectious
non_pathogenic_rows <- !pathogenic_species$Infectious
pathogenic_species$Classification <- NA
pathogenic_species$Classification[pathogenic_rows] <- "Pathogenic"
pathogenic_species$Classification[non_pathogenic_rows] <- "Non-Pathogenic"
pathogenic_species$Classification[pathogenic_species$Species == "Auxenochlorella symbiontica"] <-
  "Symbiotic"

lineage_group <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")
lineage_group <- lineage_group["Descent"]
names(lineage_group) <- "Species"
lineage_group$Group <- c("Helico",
                         rep("Cattle", 7),
                         rep("Other", 3),
                         rep("Human", 4),
                         rep("Auxeno", 6),
                         rep("Chlorella", 4))

ccap_strains <- data.frame(Species = c("Keratococcus bicaudatus",
                                       "Crucigeniella apiculata",
                                       "Coccomyxa galuniae",
                                       "Trebouxia decolorans",
                                       "Chlamydomonas reinhardtii"),
                           Group = c("Genus Keratococcus",
                                     "Family Oocystaceae",
                                     "Order incertae sedis",
                                     "Order Trebouxiales",
                                     "Class Chlorophyceae"))

lineage_group <- rbind(lineage_group, ccap_strains)

species.order <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")
species.order <- species.order["Descent"]
species.order <- data.frame(LongFormat = c(species.order$Descent, ccap_strains$Species))


species.order$ShortFormat <- unlist(lapply(species.order$LongFormat, shorten.species.name))
species.order$ShortFormat[species.order$ShortFormat=="A. sp.1"] <- "Auxenochlorella sp. 1"
species.order$ShortFormat[species.order$ShortFormat=="A. sp.2"] <- "Auxenochlorella sp. 2"

strain_numbers <- strain_numbers[1:2]
names(strain_numbers)[1] <- "Strain"

cone_numbers <- data.frame(donor = c(1,2),
                           start.date = c("2023-04-26", "2023-05-03"),
                           end.date = c("2023-04-28", "2023-05-04"))


age.experiments <- c("S5P007","S5P009")
community.experiments <- "S5P019"

cell_types <- c("Human Primary Cells", "Murine Cell Line")
names(cell_types) <- c("GM-CSF Primary", "J774.1")
#####

# Make sure to check assumptions of ANOVAs
# README: https://stats.stackexchange.com/questions/11296/using-anova-on-percentages 
# May want contingency tables for some of the comparisons (e.g. comparing cell types or strains of a species)


##### Full analysis for thesis #####
#### Identify files to use ####

f <- Sys.glob("Data/Simplified_Phagocytosis_Events_Raw/*")

#### Load all the data ####
data <- read.csv(f[1])

data$Page <- strsplit(f[1], split = "_")[[1]][6]

for(i in 2:length(f)){
  new_data <- read.csv(f[i])
  new_data$Page <- strsplit(f[i], split = "_")[[1]][6]
  data <- rbind(data, new_data)
}

#### Assign species ####

data <- left_join(data, strain_numbers)
data$Species <- factor(data$Species,
                       levels = species.order$LongFormat)

#### Identify number of replicates ####
## This is primarily because counts were done non-chronologically,
## so I had to make sure I didn't forget to count any of my replicates

data <- mutate(data, ID = paste0(Page, "_", Date))
total.files <- data %>%
  select(Strain, ID, Page, Cell_Type) %>%
  group_by(Strain, Cell_Type) %>%
  summarise(NFiles = length(unique(ID)),
            Pages = paste(unique(Page), collapse = ", "))

counted.files <- data %>%
  subset(!is.na(N_Phagocytosis_events)) %>%
  group_by(Strain, Cell_Type) %>%
  summarise(NCounted = length(unique(ID)),
            CountedPages = paste(unique(Page), collapse = ", "))


n.replicates <- full_join(total.files, counted.files)
n.replicates$NCounted[is.na(n.replicates$NCounted)] <- 0
n.replicates$Missing <- n.replicates$NFiles != n.replicates$NCounted
n.replicates <- relocate(n.replicates, Strain,Cell_Type, Missing, NFiles, NCounted)

### reorder rows for easy human reading
n.replicates$Strain <- factor(n.replicates$Strain,
                               levels = strain_numbers$Strain)
  
n.replicates <- arrange(n.replicates, Strain)
#### Produce summaries for each well ####
## These summaries are needed because different numbers of macrophages were
## counted for different days

seconds_per_frame <- read.xlsx("Data/secondsperframe.xlsx")
## Fix the format of dates in the seconds_per_frame data frame, because I was inconsistent and excel does funny things to it
seconds_per_frame$Date <- as.character(as.Date(seconds_per_frame$Date, origin = "1899-12-30"))
seconds_per_frame$Date <-unlist(lapply(lapply(strsplit(seconds_per_frame$Date, "-"), rev), paste, collapse = "/"))

data <- data %>%
  left_join(seconds_per_frame) %>%
  mutate(Time_of_first_event = (Frame_of_first_event-1)*Seconds_Per_Frame)

exp_data <- data %>%
  group_by(Strain, Page, Date, Species, Cell_Type) %>%
  summarise(PhagocytosisCount = sum(N_Phagocytosis_events != 0),
            NMacrophages = length(Macrophage),
            IngestedProtothecaMean = mean(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            IngestedProtothecaSD = sd(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            FirstEventMean = mean(Time_of_first_event, na.rm = T),
            FirstEventSD = sd(Time_of_first_event, na.rm = T)) %>%
  mutate(PhagocyticIndex = PhagocytosisCount/NMacrophages)

exp_data <- relocate(exp_data, Strain, Species,
                     Cell_Type, Page,
                     PhagocytosisCount, NMacrophages, PhagocyticIndex)


### Split data to experiment 

## Some experiments were done with additional factors included, and thus should
## not be part of the main analysis

exp_data <- subset(exp_data,!( Page %in% c(age.experiments, community.experiments)))

age_metadata <- read.csv("Data/inducedphagocytosis_age.csv")
age_exp_data <- data %>%
  subset(Strain == "HA6" & Cell_Type == "J774.1") %>%
  left_join(age_metadata) %>%
  group_by(Strain, Page, Date, Species, Cell_Type, Type) %>%
  summarise(PhagocytosisCount = sum(N_Phagocytosis_events != 0),
            NMacrophages = length(Macrophage),
            IngestedProtothecaMean = mean(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            IngestedProtothecaSD = sd(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            FirstEventMean = mean(Time_of_first_event, na.rm = T),
            FirstEventSD = sd(Time_of_first_event, na.rm = T)) %>%
  mutate(PhagocyticIndex = PhagocytosisCount/NMacrophages,
         Age = as.numeric(dmy(Date) - dmy(Type)))



community_metadata <- read.csv("Data/inducedphagocytosis_community.csv")
community_exp_data <- data %>%
  subset(Page %in% community.experiments) %>%
  left_join(community_metadata) %>%
  group_by(Strain, Page, Date, Species, Cell_Type, Type) %>%
  summarise(PhagocytosisCount = sum(N_Phagocytosis_events != 0),
            NMacrophages = length(Macrophage),
            IngestedProtothecaMean = mean(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            IngestedProtothecaSD = sd(N_Phagocytosis_events[N_Phagocytosis_events!=0]),
            FirstEventMean = mean(Time_of_first_event, na.rm = T),
            FirstEventSD = sd(Time_of_first_event, na.rm = T)) %>%
  mutate(PhagocyticIndex = PhagocytosisCount/NMacrophages)
#### Classify ####
## Want to have record of:
### Which lineage species are part of  
exp_data <- left_join(exp_data, lineage_group)

### Which species are pathogenic/symbiotic
exp_data <- left_join(exp_data, pathogenic_species)

### Which wells were from which donor?
exp_data$Date <- as.Date(exp_data$Date, format = "%d/%m/%Y")

exp_data$Donor <- NA
for(experiment in unique(exp_data$Page)){
  ## Skip experiment if it has the wrong ell type
  cell_type <- unique(exp_data$Cell_Type[exp_data$Page == experiment])
  if(!("GM-CSF Primary" %in% cell_type)) next 
  
  ## Identify the date
  experiment_date <- unique(exp_data$Date[exp_data$Page == experiment])
  
  ## Identify which donor's cells were used on that date
  for(i in 1:nrow(cone_numbers)){
    if(experiment_date >= cone_numbers$start.date[i] & experiment_date <= cone_numbers$end.date[i]){
      donor <- cone_numbers$donor[i]
    }
  }
  
  ## Record which donor this page used
  exp_data$Donor[exp_data$Page == experiment] <- donor
}

exp_data$Donor <- as.character(exp_data$Donor)


##### Comparisons #####
#### Overall comparisons ####
## Desired statistics include:
## - N macrophages phagocytosed
## - Average N Prototheca per macrophage
#### Likely to be problematic as some experiments had more algal cells
## - Time of first phagocytosis (if applicable)
#### Likely to be problematic, as some focusing took longer. 
#### Phagocytosis events that occurred before the footage started was classified
#### as occurring in frame 1
subset(exp_data, !is.na(PhagocytosisCount))

overall.results <- exp_data %>%
  group_by(Species, Cell_Type) %>%
  summarise(Replicates = length(Page),
            PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE)) %>% 
  left_join(pathogenic_species)

overall.results <- select(overall.results, Species, Cell_Type, Infectious,
                          Replicates,
                  starts_with("PhagocyticIndex"),
                  starts_with("InternalisedAlgae"),
                  starts_with("FirstPhagocytosis"))
as.data.frame(overall.results)
## FIXME: should probably arrange in species relatedness order
## 



#### Pathogenic vs non-pathogenic species? ####
### Differences in phagocytosis by murine macrophages in pathogenic species
### some are virtually phagocytosed, while others are almost all phagocytosed#run 2023-06-01
# as.data.frame(subset(exp_data, Cell_Type == "J774.1" & Infectious == TRUE & grepl("Prototheca", Species)))

mouse_path_data <- subset(exp_data, Cell_Type == "J774.1" & Infectious == TRUE & grepl("Prototheca", Species))
anova_mouse_path <- aov(PhagocyticIndex ~ Species, mouse_path_data)
# 
summary(anova_mouse_path)
TukeyHSD(anova_mouse_path)
TukeyHSD(anova_mouse_path)$Species[TukeyHSD(anova_mouse_path)$Species[,4] <= 0.05,]
a <- TukeyHSD(anova_mouse_path)$Species[TukeyHSD(anova_mouse_path)$Species[,4] <= 0.05,]
a <- TukeyHSD(anova_mouse_path)$Species
a <- as.data.frame(a)
a$Comparison <- rownames(a)
rownames(a) <- NULL
colnames(a) <- c("Difference",
                 "95% CI lower bound",
                 "95% CI upper bound",
                 "Adjusted p value",
                 "Comparison")
a[1:4] <- signif(a[1:4], digits = 3)

a$Significance <- significance.stars(a$`Adjusted p value`)
write.csv(relocate(a, Comparison),
          file = "Outputs/Pathogenic_Comparisons/THSD_pathogenic_MouseCells.csv",
          row.names = F)
### Differences in phagocytosis by human macrophages in pathogenic species
### some are not phagocytosed at all.#run 2023-06-01
#subset(exp_data, Cell_Type == "GM-CSF Primary" & Infectious == TRUE & grepl("Prototheca", Species))

human_path_data <- subset(exp_data, Cell_Type == "GM-CSF Primary" & Infectious == TRUE & grepl("Prototheca", Species))
anova_human_path <- aov(PhagocyticIndex ~ Species, human_path_data)
summary(anova_human_path)
TukeyHSD(anova_human_path)
TukeyHSD(anova_human_path)$Species[TukeyHSD(anova_human_path)$Species[,4] <= 0.05,4]

a <- TukeyHSD(anova_human_path)$Species[TukeyHSD(anova_human_path)$Species[,4] <= 0.05,4]
names(a) <- NULL
signif(a, digits = 3)

###
path_data <- subset(exp_data, Infectious == TRUE & grepl("Prototheca", Species))

results <- path_data %>%
  group_by(Species, Cell_Type, Group) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))
results$Species <- factor(results$Species,
                          levels = species.order$LongFormat)
results <- arrange(results, Cell_Type, Species)

write.csv(cbind(results[1:2],
                signif(results[4:10], digits = 3)),
          file = "Outputs/Pathogenic_Comparisons/summary_pathogenic_AllCells.csv",
          row.names = F)

(path_data_bar <- ggplot(data = results, aes(x = Species)) + 
  geom_col(aes(y = PhagocyticIndexMean,
               fill = Group))+
    scale_fill_discrete(labels = c("Cattle-\nassociated", "Human-\nassociated"))+
  geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                    ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                width = 0.3)+
  facet_wrap(.~Cell_Type, ncol = 1, labeller = as_labeller(cell_types)) +
   labs(y= "Phagocytic Index") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1))) 
  

ggsave("Outputs/Pathogenic_Comparisons/bar_pathogenic_AllCells.png",
       path_data_bar,
       height = 6,
       width = 7.5,
       units = "in")

as.data.frame(results)

anova_pathogens <- aov(PhagocyticIndex ~ Species*Cell_Type, data = path_data)
summary(anova_pathogens) # significant differences between the cow pathogenic species and all others #run 2023-06-01
path_comparisons <- TukeyHSD(anova_pathogens)$Species
path_comparisons <- as.data.frame(path_comparisons)
path_comparisons$round_Pajd <- round(path_comparisons$`p adj`, digits = 4)
path_comparisons
path_comparisons[path_comparisons$round_Pajd<=0.05,]

anova_groups <- aov(PhagocyticIndex ~ Group*Cell_Type, data = exp_data)
summary(anova_groups) # significant differences only between the cow pathogenic group and all others #run 2023-06-01
groups_comparisons <- TukeyHSD(anova_groups)$Group
groups_comparisons <- as.data.frame(groups_comparisons)
groups_comparisons$round_Pajd <- round(groups_comparisons$`p adj`, digits = 4)
groups_comparisons

#### within/between Lineage/group, including non-pathogens ####
### cattle species
cow_data <- subset(exp_data, Group == "Cattle")
# cow_data <- subset(exp_data, Group == "Cattle" & Cell_Type == "J774.1")

anova_cattle <- aov(PhagocyticIndex ~ Species*Cell_Type, data = cow_data)
# anova_cattle <- aov(PhagocyticIndex ~ Species, data = cow_data)
summary(anova_cattle) 
cattle_comparisons <- TukeyHSD(anova_cattle)$Species
cattle_comparisons <- as.data.frame(cattle_comparisons)
cattle_comparisons$round_Pajd <- round(cattle_comparisons$`p adj`, digits = 4)
cattle_comparisons

cow_results <- cow_data %>%
  group_by(Species, Cell_Type, Infectious, Classification) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))

cow_results$Species <- factor(cow_results$Species,
                          levels = species.order$LongFormat)

(cow_data_bar <- ggplot(data = cow_results, aes(x = Species)) + 
  geom_col(aes(y = PhagocyticIndexMean,
               fill = Classification))+
    scale_fill_manual(values = c("Pathogenic" = "red", 
                                 "Non-Pathogenic" = "darkgrey",
                                 "Symbiotic" = "green"))+
  geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                    ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                width = 0.3)+
  facet_wrap(.~Cell_Type, ncol = 1, labeller = as_labeller(cell_types)) +
  labs(y= "Phagocytic Index") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1)))

ggsave("Outputs/Cattle_Comparisons/bar_cow_AllCells.png",
       cow_data_bar,
       height = 6,
       width = 6,
       units = "in")

### human species
human_data <- subset(exp_data, Group == "Human")
# human_data <- subset(exp_data, Group == "Human" & Cell_Type == "J774.1")

anova_human <- aov(PhagocyticIndex ~ Species*Cell_Type, data = human_data)
# anova_human <- aov(PhagocyticIndex ~ Species, data = human_data)
summary(anova_human) # significant difference by species within human group, but no significant pairwise difference #run 2023-06-01
human_comparisons <- TukeyHSD(anova_human)$Species
human_comparisons <- as.data.frame(human_comparisons)
human_comparisons$round_Pajd <- round(human_comparisons$`p adj`, digits = 4)
human_comparisons

if(all(human_data$Cell_Type == "J774.1")){
  a <- TukeyHSD(anova_human)$Species
  a <- as.data.frame(a)
  a$Comparison <- rownames(a)
  rownames(a) <- NULL
  colnames(a) <- c("Difference",
                   "95% CI lower bound",
                   "95% CI upper bound",
                   "Adjusted p value",
                   "Comparison")
  a[1:4] <- signif(a[1:4], digits = 3)
  
  a$Significance <- significance.stars(a$`Adjusted p value`)
  write.csv(relocate(a, Comparison),
            file = "Outputs/Human_Comparisons/THSD_human_MouseCells.csv",
            row.names = F)
}


human_results <- human_data %>%
  group_by(Species, Cell_Type, Infectious, Classification) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))
human_results$Species <- factor(human_results$Species,
                          levels = species.order$LongFormat)
options(scipen = 0)
human_data$PhagocyticIndex <- as.numeric(format(human_data$PhagocyticIndex, scientific = F))

(human_data_bar <- ggplot(data = human_results, aes(x = Species)) + 
    geom_col(aes(y = PhagocyticIndexMean,
                 fill = Classification))+
    scale_fill_manual(values = c("Pathogenic" = "red", 
                                 "Non-Pathogenic" = "darkgrey",
                                 "Symbiotic" = "green"))+
    geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                      ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                  width = 0.3)+
    facet_wrap(.~Cell_Type, ncol = 1, labeller = as_labeller(cell_types)) +
    labs(y= "Phagocytic Index") +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1)))

ggsave("Outputs/Human_Comparisons/bar_human_AllCells.png",
       human_data_bar,
       height = 6,
       width = 6,
       units = "in")


human_data_1 <- subset(exp_data, Group == "Human" & Cell_Type == "J774.1")
anova_human_1 <- aov(PhagocyticIndex ~ Species, data = human_data)
summary(anova_human_1) 
human_comparisons_1 <- TukeyHSD(anova_human_1)$Species
human_comparisons_1 <- as.data.frame(human_comparisons_1)
human_comparisons_1$Comparison <- rownames(human_comparisons_1)
rownames(human_comparisons_1) <- NULL
colnames(human_comparisons_1) <- c("Difference",
                 "95% CI lower bound",
                 "95% CI upper bound",
                 "Adjusted p value",
                 "Comparison")
human_comparisons_1[1:4] <- signif(human_comparisons_1[1:4], digits = 3)


write.csv(relocate(human_comparisons_1, Comparison),
          file = "Outputs/Human_Comparisons//THSD_human_MouseCells.csv",
          row.names = F)


### STU lineage/other
other_data <- subset(exp_data, Group == "Other")
# other_data <- subset(exp_data, Group == "Other" & Cell_Type != "J774.1")

anova_other <- aov(PhagocyticIndex ~ Species*Cell_Type, data = other_data)
# anova_other <- aov(PhagocyticIndex ~ Species, data = other_data)
summary(anova_other) 
other_comparisons <- TukeyHSD(anova_other)$Species
other_comparisons <- as.data.frame(other_comparisons)
other_comparisons$round_Pajd <- round(other_comparisons$`p adj`, digits = 4)
other_comparisons

other_results <- other_data %>%
  group_by(Species, Cell_Type, Infectious, Classification) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))

other_results$Species <- factor(other_results$Species,
                          levels = species.order$LongFormat)

(other_data_bar <- ggplot(data = other_results, aes(x = Species)) + 
    geom_col(aes(y = PhagocyticIndexMean,
                 fill = Classification))+
    scale_fill_manual(values = c("Pathogenic" = "red", 
                                 "Non-Pathogenic" = "darkgrey",
                                 "Symbiotic" = "green"))+
    geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                      ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                  width = 0.3)+
    facet_wrap(.~Cell_Type, ncol = 2, labeller = as_labeller(cell_types)) +
    labs(y= "Phagocytic Index") +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1)))

ggsave("Outputs/Environmental_Comparisons/bar_other_AllCells.png",
       other_data_bar,
       height = 3,
       width = 6,
       units = "in")



#### Auxenochlorella data ####
full_auxeno_data <- subset(exp_data, Strain %in% c("HP53", "HP54", "HA1", "HA2", "HA5", "HA6", "HA7"))
auxeno_data <- subset(full_auxeno_data, as.Date(Date, format = "%d/%m/%Y") < "2023-05-05")
# auxeno_data <- subset(auxeno_data, Cell_Type != "J774.1")

## Possible outliers

# auxeno_data <- subset(auxeno_data, !(Strain == "HA6" & Page == "S4P141")) # this culture was contaminated

# auxeno_data <- subset(exp_data, Group == "Auxeno" & Cell_Type != "J774.1")

auxeno_species_results <- auxeno_data %>%
  group_by(Species, Cell_Type, Infectious, Classification) %>% 
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))
as.data.frame(auxeno_species_results) 

auxeno_species_results$Species <- factor(auxeno_species_results$Species,
                          levels = species.order$LongFormat)

(auxeno_species_data_bar <- ggplot(data = auxeno_species_results, aes(x = Species)) + 
    geom_col(aes(y = PhagocyticIndexMean,
                 fill = Classification))+
    scale_fill_manual(values = c("Pathogenic" = "red", 
                                 "Non-Pathogenic" = "darkgrey",
                                 "Symbiotic" = "green"))+
    geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                      ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                  width = 0.3)+
    facet_wrap(.~Cell_Type, nrow = 2, labeller = as_labeller(cell_types)) +
    labs(y= "Phagocytic Index") +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1))) 

ggsave("Outputs/Auxenochlorella_Comparisons/bar_auxenoSp_AllCells.png",
       auxeno_species_data_bar,
       height = 6,
       width = 7.5,
       units = "in")

anova_auxeno <- aov(PhagocyticIndex ~ Species*Cell_Type, data = auxeno_data)
# anova_auxeno <- aov(PhagocyticIndex ~ Species, data = auxeno_data)
summary(anova_auxeno) # no significant differences within Auxeno group #run 2023-06-01
auxeno_comparisons <- TukeyHSD(anova_auxeno)$Species
auxeno_comparisons <- as.data.frame(auxeno_comparisons)
auxeno_comparisons$round_Pajd <- round(auxeno_comparisons$`p adj`, digits = 4)
auxeno_comparisons


## Looks like Auxeno is usually low except for HA1. But there are a few data 
## points that are quite separate. It's hard to tell if it's really that 
## different from other groups/strains
auxeno_Jdata <- subset(auxeno_data, Cell_Type == "J774.1")
auxeno_Jdata$Strain <- factor(auxeno_Jdata$Strain, levels = strain_numbers$Strain)
auxeno_Jdata$Species <- factor(auxeno_Jdata$Species, levels = species.order$LongFormat)
replicate_diversity <- ggplot(auxeno_Jdata, aes(y = PhagocyticIndex, x = Strain))+
  geom_point(aes(colour = Species))  + coord_cartesian(ylim = 0:1)
ggsave("Outputs/Auxenochlorella_Comparisons/point_ReplicateDiversity_MouseCells.png",
       replicate_diversity,
       width = 7.5,
       height = 3,
       units = "in")



### Induced phagocytosis: age

# S5P007 went terribly out of focus. Usable video starts 1.5 hours after infection
# S5P013 was also unusably unfocused. No other video available

# backup_age_exp_data <- age_exp_data
# age_exp_data <- backup_age_exp_data
age_exp_data$Type <- as.character(dmy(age_exp_data$Type))

age_exp_data<- left_join(age_exp_data, data.frame(Type = unique(age_exp_data$Type),
                                   Culture = LETTERS[1:length(unique(age_exp_data$Type))]) )

# the two cultures that are not phagocytosed are the youngest, but one culture
# is phagocytosed at 16 days and another at 16 days is not
arrange(age_exp_data, Age)

options(scipen = 0)
inducible_age <- ggplot(age_exp_data, aes(x = Culture, y = PhagocyticIndex)) + 
  geom_point(aes(col = Age)) + 
  coord_cartesian(ylim = 0:1)

ggplot(age_exp_data, aes(x = Age, y = PhagocyticIndex)) + 
  geom_point(aes(shape = Culture)) + 
  coord_cartesian(ylim = 0:1, xlim = c(0, max(age_exp_data$Age))) 

  
ggsave("Outputs/Auxenochlorella_Comparisons/point_Inducible-HA6-Age_MouseCells.png",
       inducible_age, 
       height = 3,
       width = 3,
       units = "in")

ggplot(age_exp_data, aes(x = Age, y = PhagocyticIndex)) +
  geom_point(aes(colour = Culture)) + 
  geom_smooth(method= "lm") +
  coord_cartesian(ylim = 0:1)

cor(age_exp_data$Age, age_exp_data$PhagocyticIndex, method = "pearson")
## not great correlation
summary(lm(age_exp_data$Age ~ age_exp_data$PhagocyticIndex))
## not linear relationship

### Induced phagocytosis: community
inducible_community <- ggplot(community_exp_data, aes(x = Strain, y = PhagocyticIndex)) + 
  geom_point(aes(col = Type)) + coord_cartesian(ylim = 0:1)

# phagocytosis not obviously caused by community, though only one replicate. 
ggsave("Outputs/Auxenochlorella_Comparisons/point_Inducible-HA1,2,5-Community_MouseCells.png",
       inducible_community, 
       height = 3,
       width = 3,
       units = "in")

### Include all auxeno data
## format the different data frames so they are compatible
auxeno_data$Type <- NA
age_exp_data$Date <- as.Date(age_exp_data$Date, format = "%d/%m/%Y")
community_exp_data$Date <- as.Date(community_exp_data$Date, format = "%d/%m/%Y")
names(community_exp_data)
auxeno_data <- select(auxeno_data, !Infectious)

combined_auxeno_data <- rbind(auxeno_data,age_exp_data,community_exp_data)
combined_auxeno_data$Strain <- factor(combined_auxeno_data$Strain,
                                      levels = strain_numbers$Strain)

combined_auxeno_Jdata <- subset(combined_auxeno_data, Cell_Type == "J774.1")
combined_auxeno_Jdata$Species <- factor(combined_auxeno_Jdata$Species,
                                        levels = species.order$LongFormat)

replicate_diversity2 <- ggplot(combined_auxeno_Jdata, aes(x = Strain, y = PhagocyticIndex)) +
  geom_point(aes(colour = Species)) + coord_cartesian(ylim = 0:1)

ggsave("Outputs/Auxenochlorella_Comparisons/point_ReplicateDiversityFullData_MouseCells.png",
       replicate_diversity2,
       width = 7.5,
       height = 3,
       units = "in")

combined_auxeno_species_results <- combined_auxeno_data %>%
  subset(select = -Classification) %>%
  left_join(pathogenic_species) %>%
  group_by(Species, Cell_Type, Infectious, Classification) %>% 
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))
combined_auxeno_species_results$Species <- factor(combined_auxeno_species_results$Species,
                          levels = species.order$LongFormat)



(auxeno_species_fulldata_bar <- ggplot(data = combined_auxeno_species_results, aes(x = Species)) + 
    geom_col(aes(y = PhagocyticIndexMean,
                 fill = Classification))+
    scale_fill_manual(values = c("Pathogenic" = "red", 
                                 "Non-Pathogenic" = "darkgrey",
                                 "Symbiotic" = "green"))+
    geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                      ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                  width = 0.3)+
    facet_wrap(.~Cell_Type, nrow = 2, labeller = as_labeller(cell_types)) +
    labs(y= "Phagocytic Index") +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5)) + 
    coord_cartesian(ylim = 0:1))

ggsave("Outputs/Auxenochlorella_Comparisons/bar_auxenoSpFulldata_AllCells.png",
       auxeno_species_fulldata_bar,
       height = 6,
       width = 7.5,
       units = "in")




anova_auxenoCombined <- aov(PhagocyticIndex ~ Species*Cell_Type, data = combined_auxeno_data)
# combined_auxeno_data_cell <- subset(combined_auxeno_data, Cell_Type == "J774.1")
# anova_auxenoCombined <- aov(PhagocyticIndex ~ Species, data = combined_auxeno_data_cell)
summary(anova_auxenoCombined) 
auxeno_comparisons <- TukeyHSD(anova_auxenoCombined)$Species
auxeno_comparisons <- as.data.frame(auxeno_comparisons)
auxeno_comparisons$round_Pajd <- round(auxeno_comparisons$`p adj`, digits = 4)
auxeno_comparisons


if(all(combined_auxeno_data$Cell_Type == "J774.1")){
  a <- TukeyHSD(anova_auxenoCombined)$Species
  a <- as.data.frame(a)
  a$Comparison <- rownames(a)
  rownames(a) <- NULL
  colnames(a) <- c("Difference",
                   "95% CI lower bound",
                   "95% CI upper bound",
                   "Adjusted p value",
                   "Comparison")
  a[1:4] <- signif(a[1:4], digits = 3)
  
  a$Significance <- significance.stars(a$`Adjusted p value`)
  write.csv(relocate(a, Comparison),
            file = "Outputs/Auxenochlorella_Comparisons/THSD_Auxeno_MouseCells.csv",
            row.names = F)
}


#### Timing and counts #### 
groups <- c("Cattle", 
            "Other", # true environmentals
            "Human",
            "Auxeno")

global_data <- exp_data
global_data$Type <- NA

global_data <- rbind(global_data,age_exp_data,community_exp_data)
global_data$Strain <- factor(global_data$Strain,
                                      levels = strain_numbers$Strain)
global_data <- select(global_data, !Group)
global_data <- left_join(global_data, lineage_group)
global_data <- subset(global_data, grepl("H[AP]", Strain))

global_results <- global_data %>%
  group_by(Species, Cell_Type, Group) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSE = std.err(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSE = std.err(FirstEventMean, na.rm = TRUE))
global_results$Species <- factor(global_results$Species,
                          levels = species.order$LongFormat)


global_results$Grouping <- NA
for(i in 1:length(global_results$Group)){
  global_results$Grouping[i] <- switch(global_results$Group[i],
        Auxeno = "Auxenochlorella",
        Cattle = "Cattle-Associated",
        Human = "Human-Associated",
        Other = "True Environmentals")
}


ggplot(global_results, aes(x = Species, y = PhagocyticIndexMean)) +
  geom_col(aes(fill = Grouping)) + 
  geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                    ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                width = 0.3) +
  facet_wrap(Cell_Type~., nrow = 2, labeller = as_labeller(cell_types)) 



### Identify differences between the number of algae internalised
(bar_internalised <- ggplot(global_results, aes(x = Species, y = InternalisedAlgaeMean)) +
  geom_col(aes(fill = Grouping)) + 
  geom_errorbar(aes(ymin = InternalisedAlgaeMean - 2*InternalisedAlgaeSE,
                    ymax = InternalisedAlgaeMean + 2*InternalisedAlgaeSE),
                width = 0.3) +
  facet_wrap(.~Cell_Type, ncol = 1, labeller = as_labeller(cell_types)) +
    labs(y = "Number of Internalised Algal Cells")+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5)))
  

ggsave("Outputs/Other_Comparisons/bar_InternalisedAlgae_allcells.png",
       bar_internalised,
       height = 6,
       width = 12,
       units = "in")

## Check within each group for differences
for(i in 1:length(groups)){
  print(groups[i])
internalised <- aov(IngestedProtothecaMean~Species*Cell_Type,
                    subset(global_data, Group == groups[i]))
print(summary(internalised))
}

## Differeneces withing groups were only found for human strains
internalised <- aov(IngestedProtothecaMean~Species,
                    subset(global_data, 
                           Group == "Human" & Cell_Type == "J774.1"))
summary(internalised)
TukeyHSD(internalised)
## But no pairwise differences were identified

## Check between groups for differences
internalised <- aov(IngestedProtothecaMean~Group*Cell_Type,global_data)
summary(internalised)

internalised_mouse <- aov(IngestedProtothecaMean~Group,
                    subset(global_data, 
                           Cell_Type == "J774.1"
                           ))
summary(internalised_mouse)
internalised_mouse_df <- as.data.frame(TukeyHSD(internalised_mouse)$Group)
internalised_mouse_df$celltype <- "J774.1"
internalised_mouse_df$comparison <- rownames(internalised_mouse_df)

internalised_human <- aov(IngestedProtothecaMean~Group,
                    subset(global_data, 
                           Cell_Type != "J774.1"
                    ))
summary(internalised_human)
internalised_human_df <- as.data.frame(TukeyHSD(internalised_human)$Group)
internalised_human_df$celltype <- "GM-CSF Primary"
internalised_human_df$comparison <- rownames(internalised_human_df)

internalised_df <- rbind(internalised_human_df,
                         internalised_mouse_df)


rownames(internalised_df) <- NULL
internalised_df<- relocate(internalised_df, comparison, celltype)
internalised_df[,3:6]<- signif(internalised_df[,3:6], digits = 4)
internalised_df$comparison <- gsub("Other", "Environmental", internalised_df$comparison)
colnames(internalised_df) <- c("Group Comparison",
                               "Cell Type",
                               "Difference", 
                               "95% CI lower bound",
                               "95% CI upper bound",
                               "Adjusted p value")

write.csv(internalised_df,
          "Outputs/Other_Comparisons/THSD_InternalisedGroups_AllCells.csv",
          row.names = F)

## Cattle significantly different from human and Auxeno in murine cells
## Cattle significantly different from all groups in human cells

### Identify difference between time taken to complete first phagocytosis 
options(scipen = 10)
(bar_timing <- ggplot(global_results, aes(x = Species, y = PhagocytosisRateMean)) +
    geom_col(aes(fill = Grouping)) + 
    geom_errorbar(aes(ymin = PhagocytosisRateMean - 2*PhagocytosisRateSE,
                      ymax = PhagocytosisRateMean + 2*PhagocytosisRateSE),
                  width = 0.3) +
    facet_wrap(.~Cell_Type, ncol = 1, labeller = as_labeller(cell_types)) +    
    labs(y = "Time to Complete First Phagocytosis Event (Seconds)")+
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5)))

ggsave("Outputs/Other_Comparisons/bar_PhagocyticRate_allcells.png",
       bar_timing,
       height = 6,
       width = 12,
       units = "in")


for(i in 1:length(groups)){
  print(groups[i])
  timing <- aov(FirstEventMean~Species*Cell_Type,
                      subset(global_data, Group == groups[i]))
  print(summary(timing))
}

## Differences within groups were seen for both human and auxeno groups
timing <- aov(FirstEventMean~Species,
                    subset(global_data, 
                           Group == "Auxeno" & Cell_Type == "J774.1"))
summary(timing)
TukeyHSD(timing) 



timing_auxeno_df <- as.data.frame(TukeyHSD(timing)$Species)
timing_auxeno_df$comparison <- rownames(timing_auxeno_df)

rownames(timing_auxeno_df) <- NULL
timing_auxeno_df<- relocate(timing_auxeno_df, comparison)
timing_auxeno_df[,2:5]<- signif(timing_auxeno_df[,2:5], digits = 4)

colnames(timing_auxeno_df) <- c("Species Comparison",
                               "Difference", 
                               "95% CI lower bound",
                               "95% CI upper bound",
                               "Adjusted p value")

timing_auxeno_df[2:5] <- signif(timing_auxeno_df[2:5], digits = 3)

timing_auxeno_df$Significance <- significance.stars(timing_auxeno_df$`Adjusted p value`)

write.csv(timing_auxeno_df,
          "Outputs/Other_Comparisons/THSD_TimingAuxeno_MouseCells.csv",
          row.names = F)


## no pairwise differences in human strains
## Significant difference between Auxenochlorella sp1 and all other Auxenochlorella strains,
## as well as between xanthoriae and symbiontica

## Check between groups for differences
timing <- aov(FirstEventMean~Group*Cell_Type,global_data)
summary(timing)

timing <- aov(FirstEventMean~Group,
                    subset(global_data, 
                           Cell_Type != "J774.1"
                    ))
summary(timing)
TukeyHSD(timing)
## Cattle significantly different from all other groups in murine
## and significant difference between Auxeno and true environmentals in murine
## Significant differences between cattle and both Auxeno and Human in human cells

## Check between groups for differences
timing <- aov(FirstEventMean~Group*Cell_Type,global_data)
summary(timing)

timing_mouse <- aov(FirstEventMean~Group,
                          subset(global_data, 
                                 Cell_Type == "J774.1" 
                          ))
summary(timing_mouse)
timing_mouse_df <- as.data.frame(TukeyHSD(timing_mouse)$Group)
timing_mouse_df$celltype <- "J774.1"
timing_mouse_df$comparison <- rownames(timing_mouse_df)

timing_human <- aov(FirstEventMean~Group,
                          subset(global_data, 
                                 Cell_Type != "J774.1"
                          ))
summary(timing_human)
timing_human_df <- as.data.frame(TukeyHSD(timing_human)$Group)
timing_human_df$celltype <- "GM-CSF Primary"
timing_human_df$comparison <- rownames(timing_human_df)

timing_df <- rbind(timing_human_df,
                         timing_mouse_df)


rownames(timing_df) <- NULL
timing_df<- relocate(timing_df, comparison, celltype)
timing_df[,3:6]<- signif(timing_df[,3:6], digits = 4)
timing_df$comparison <- gsub("Other", "Environmental", timing_df$comparison)
colnames(timing_df) <- c("Group Comparison",
                               "Cell Type",
                               "Difference", 
                               "95% CI lower bound",
                               "95% CI upper bound",
                               "Adjusted p value")

write.csv(timing_df,
          "Outputs/Other_Comparisons/THSD_TimingGroups_AllCells.csv",
          row.names = F)

#####
#### Non AHP algae ####
ccap_data <- subset(exp_data, grepl("CCAP", Strain))

ccap_data$Species <- factor(ccap_data$Species, 
                            levels = species.order$LongFormat)

ccap_results <- ccap_data %>%
  group_by(Species) %>%
  summarise(PhagocyticIndexMean = mean(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSD = sd(PhagocyticIndex, na.rm = TRUE),
            PhagocyticIndexSE = std.err(PhagocyticIndex),
            InternalisedAlgaeMean = mean(IngestedProtothecaMean, na.rm = TRUE),
            InternalisedAlgaeSD = sd(IngestedProtothecaMean, na.rm = TRUE),
            PhagocytosisRateMean = mean(FirstEventMean, na.rm = TRUE),
            PhagocytosisRateSD = sd(FirstEventMean, na.rm = TRUE))

ccap_bar <- ggplot(data = ccap_results, aes(x = Species, y = PhagocyticIndexMean))+
  geom_col(aes(fill = Species))  + 
  geom_errorbar(aes(ymin = PhagocyticIndexMean - 2*PhagocyticIndexSE,
                    ymax = PhagocyticIndexMean + 2*PhagocyticIndexSE),
                width = 0.3)+
  labs(y= "Phagocytic Index") +
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
  scale_y_continuous(breaks=seq(from = -0.25, to = 1, by = 0.25)) +
  coord_cartesian(ylim = c(-0.25,1))


ggsave("Outputs/CCAP/bar_ccap_MouseCells.png",
       ccap_bar,
              height = 4,
              width = 7.5,
              units = "in")
ccap_data$Date <- as.character(ccap_data$Date)
(ccap_points <- ggplot(ccap_data, aes(x = Species, y = PhagocyticIndex))+
    geom_point(aes(colour = Date))+ 
    labs(y= "Phagocytic Index") +
    theme(axis.text.x = element_text(angle = 30, vjust = 0.6, hjust=0.5))+
    coord_cartesian(ylim = c(0,1)))

ggsave("Outputs/CCAP/point_ccap_MouseCells.png",
       ccap_points,
       height = 4,
       width = 7.5,
       units = "in")

# will want additional replicates to show means
# might want some prototheca for context

exp_data <- exp_data[!is.na(exp_data$Species),]
ggplot(data = exp_data, aes(x = Group, y = PhagocyticIndex)) +
  geom_boxplot() + coord_cartesian(ylim = 0:1)

ccap_species <- unique(exp_data$Species[grepl("CCAP", exp_data$Strain)])
for(i in 1:5){
  s <- ccap_species[i]
  exp_data$Group[exp_data$Species == s] <- s
}



#### Difference between primary donors ####

primary.data <- subset(exp_data, Cell_Type == "GM-CSF Primary")
primary.data$Species <- unlist(lapply(primary.data$Species, shorten.species.name))
primary.data$Species[primary.data$Species=="A. sp.1"] <- "Auxenochlorella sp. 1"
primary.data$Species[primary.data$Species=="A. sp.2"] <- "Auxenochlorella sp. 2"

primary.data$Species <- factor(primary.data$Species, levels = species.order$ShortFormat)
primary.data <- arrange(primary.data, Species)

ggplot(data = primary.data, aes(y = PhagocyticIndex, x = Group))+
  geom_boxplot(aes(colour = Donor)) + coord_cartesian(ylim = 0:1) # Should this be a bar chart?
  

n.species <-length(unique(primary.data$Species))
n.graph.rows <- ceiling(n.species/10)
n.species.per.row <- ceiling(n.species/n.graph.rows)
graph.species.rows <- data.frame(Species = unique(primary.data$Species),
                                 Row = rep(1:n.graph.rows, each = n.species.per.row)[1:n.species])
primary.data <- left_join(primary.data, graph.species.rows)

donor_box <- ggplot(data = subset(primary.data, !is.na(PhagocyticIndex)), aes(y = PhagocyticIndex, x = Species))+
  geom_boxplot(aes(colour = Donor)) +
  facet_wrap(Row~., ncol = 1, scales = "free") +
  theme(strip.text.x = element_blank(),
  axis.text.x=element_text(angle = 15)) +
  coord_cartesian(ylim = 0:1)

donor_points <- ggplot(data = subset(primary.data, !is.na(PhagocyticIndex)), aes(y = PhagocyticIndex, x = Species))+
  geom_point(aes(colour = Donor)) +
  facet_wrap(Row~., ncol = 1, scales = "free") +
  theme(strip.text.x = element_blank(),
        axis.text.x=element_text(angle = 15)) +
  coord_cartesian(ylim = 0:1)

ggsave("Outputs/Other_Comparisons/boxplot_donorComp.png",
       donor_box,
       height = 7, 
       width = 9,
       units = "in")
ggsave("Outputs/Other_Comparisons/points_donorComp.png",
       donor_points,
       height = 7, 
       width = 9,
       units = "in")

anova_donor_comparisons <- aov(PhagocyticIndex ~ Species*Donor, data = primary.data)
summary(anova_donor_comparisons)
TukeyHSD(anova_donor_comparisons)


#### Differences between strains of the same species ####
# On reflection, I can't make a great argument for the different not being that
# one strain had all three replicates with donor 1, while the other had two
# replicates with donor 2. I did initially say that wickerhamii not showing a
# difference between the donors suggested that the difference between bovis
# strains was genuine, but then I remembered that wickerhamii is one of the
# species that isn't different between the donors.

# I do still think there might be strain-specific differences, I just don't
# think this data is sufficient to show that

### difference between bovis strains
subset(exp_data, Species == "Prototheca bovis"  & Cell_Type == "J774.1")
HP3_data  <- subset(exp_data, Strain == "HP3"  & Cell_Type != "J774.1")
HP40_data <-  subset(exp_data, Strain == "HP40"& Cell_Type != "J774.1")

p1 <- ggplot(data = subset(data, Species == "Prototheca bovis" & Cell_Type == "J774.1"),
             aes(x = N_Phagocytosis_events), binwidth = 3) +
  geom_histogram(aes(fill = Strain)) +
  facet_wrap(~Strain) + 
  labs(x = "Number of Algal Cells Internalised")

p2 <- ggplot(data = subset(data, Species == "Prototheca bovis" & Cell_Type == "J774.1"), 
             aes(x = Time_of_first_event/60)) + 
  geom_histogram(aes(fill = Strain), binwidth = 3) +
  facet_wrap(~Strain) +
  labs(x = "Time of First Phagocytosis Event (Minutes)")

p3 <- ggarrange(p1,p2,
          nrow = 2,
          common.legend = T,
          legend = "none")

ggsave("Outputs/Other_Comparisons/bar_bovisStrainDiff_humanCells.png",
       p3,
       height = 9,
       width = 9,
       units = "in")

p.vals <- numeric(length = 3L)
p.vals[1] <- t.test(HP3_data$PhagocytosisCount, HP40_data$PhagocytosisCount)$p.value
p.vals[2] <- t.test(HP3_data$IngestedProtothecaMean, HP40_data$IngestedProtothecaMean)$p.value
p.vals[3] <- t.test(HP3_data$FirstEventMean, HP40_data$FirstEventMean)$p.value
p.vals
p.vals.adj <- p.adjust(p.vals, method = "BH")
p.vals.adj
# probably worth including

### difference between wickerhamii
HP50_data <- subset(exp_data, Strain == "HP50" & Cell_Type != "J774.1")
HP52_data <-  subset(exp_data, Strain == "HP52"& Cell_Type != "J774.1") 

p1 <- ggplot(data = subset(data, Species == "Prototheca wickerhamii"& Cell_Type != "J774.1"), 
             aes(x = N_Phagocytosis_events)) +
  geom_histogram(aes(fill = Strain), binwidth = 3) +
  facet_wrap(~Strain) + 
  labs(x = "Number of Algal Cells Internalised")

p2 <- ggplot(data = subset(data, Species == "Prototheca wickerhamii"& Cell_Type != "J774.1"),
             aes(x = Time_of_first_event/60)) + # FIXME: use time, not frames
  geom_histogram(aes(fill = Strain), binwidth = 3) +
  facet_wrap(~Strain) +
  labs(x = "Time of First Phagocytosis Event (Minutes)")

ggarrange(p1,p2,
          nrow = 2,
          common.legend = T,
          legend = "none")

p.vals <- numeric(length = 3L)
p.vals[1] <- t.test(HP50_data$PhagocytosisCount, HP52_data$PhagocytosisCount)$p.value
p.vals[2] <- t.test(HP50_data$IngestedProtothecaMean, HP52_data$IngestedProtothecaMean)$p.value
p.vals[3] <- t.test(HP50_data$FirstEventMean, HP52_data$FirstEventMean)$p.value
p.vals
p.vals.adj <- p.adjust(p.vals, method = "BH")
p.vals.adj


##### Outputs #####

as.data.frame(n.replicates)
# HP4 S4P180 is uncountable

as.data.frame(subset(exp_data, Strain == "HA6"))
as.data.frame(subset(exp_data, Cell_Type != "J774.1" &
                       Strain == "HP27"))

n.replicates.species <- n.replicates %>%
  left_join(strain_numbers) %>%
  group_by(Species, Cell_Type) %>%
  summarise(Replicates = sum(NCounted),
            Strains = length(unique(Strain)))


## 
results


#### Export ####

ggsave("Outputs/PhagocytosisHistByStrain.png",
       p1,
       width = 6,
       height = 6)

write.csv(results,
          "Outputs/Results.csv",
          row.names = F)

write.csv(exp_data,
          "Outputs/ExpData.csv",
          row.names = F)

write.csv(n.replicates.species,
          "Outputs/Replicates.csv",
          row.names = F)

#####

