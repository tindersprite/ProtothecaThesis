##### Setup #####
#### Packages ####
library(ggplot2)
library(tidyverse)


#### Functions ####
shorten.species.name <- function(s){
  split.name <- strsplit(s, split = " ")[[1]]
  genus.shortened <- paste0(substring(split.name[1],1,1), ".")
  short.name <- paste(genus.shortened, split.name[2])
  return(short.name)
}

#### Levels ####
strain.levels <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")
strain.levels <- strain.levels$H_number

species.levels <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")
species.order <- unlist(lapply(species.levels$Descent, "shorten.species.name"))
species.order[species.order == "H. sp."] <- "Helicosporidium sp. Sj"
species.order[species.order == "A. sp.1"] <- "Auxenochlorella sp. 1"
species.order[species.order == "A. sp.2"] <- "Auxenochlorella sp. 2"

#####

#### Load data ####
busco.results <- read.csv("Data/busco_scores.csv") 



#### Format data #### 
### Arrange strain IDs across two lines (which looks better)
busco.results <- mutate(busco.results, ID = paste(Species, Strain, sep = "\n"))

### Standardise busco format for text labelling
busco.results <- mutate(busco.results, Written.Score = paste("C:", C,
                                                             " [S:", S,
                                                             ", D:", D,
                                                             "], F:", `F`,
                                                             ", M:", M,
                                                             ", n:", N,
                                                             sep = ""))

### Set scores to same length, which makes it easier to align them on the plot
score.length <- max(nchar(busco.results$Written.Score))
for(i in 1:nrow(busco.results)){
  while(nchar(busco.results$Written.Score[i]) < score.length){
    busco.results$Written.Score[i] <- paste(busco.results$Written.Score[i], " ", sep = "")
  }
}

### Tidy format
busco.results <- pivot_longer(busco.results, 
                              cols = c('S','D','F','M'),
                              names_to = 'Category', values_to = 'Score')

### Set levels
busco.levels <- busco.results[c("Species", "Strain")]
busco.levels$Species <- factor(busco.levels$Species,
                               levels = species.order)

strain.levels <- unique(c(strain.levels, busco.results$Strain))
busco.levels$Strain <- factor(busco.levels$Strain,
                              levels = strain.levels)
busco.levels <- arrange(busco.levels, Species, Strain)
id.levels <- paste(busco.levels$Species, busco.levels$Strain, sep = "\n")
id.levels <- unique(id.levels)

#busco.results$ID <- factor(busco.results$ID, 
#                           levels = rev(unique(busco.results$ID)))
busco.results$ID <- factor(busco.results$ID, 
                           levels = rev(id.levels))
busco.results$Category <- factor(busco.results$Category,
                                 levels = rev(c('S','D','F','M')))

strain.levels <- union(strain.levels, busco.results$Strain)
busco.results$Strain <- factor(busco.results$Strain,
                               levels = rev(strain.levels))
##### Plot data #####
#### All genomes ####
combined <- ggplot(data = busco.results, aes(x = ID, y = Score)) + 
  geom_bar(aes(fill = Category), position = 'fill', stat = "identity")+
  scale_fill_manual(values = c("red", "yellow", "deepskyblue3", "deepskyblue"),
                    labels = c("Missing", "Fragmented", "Duplicated", "Single")) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Species and Strain Identifier for Assembly")

combined.text <- combined + geom_text(data = subset(busco.results, Category == "S"),
                                      aes(label = Written.Score),
                                      y = 0.22, size = 2.5)

ggsave("Outputs/combined_busco_plot.png",
       plot = combined,
       width = 8,
       height = length(unique(busco.results$ID))*0.5+1)

ggsave("Outputs/combined_busco_text_plot.png",
       plot = combined.text,
       width = 8,
       height = length(unique(busco.results$ID))*0.5+1)

#### Internal genomes ####
my.busco.results <- subset(busco.results, Source == "Internal")

my_buscos <- ggplot(data = my.busco.results, aes(x = ID, y = Score)) + 
  geom_bar(aes(fill = Category), position = 'fill', stat = "identity")+
  scale_fill_manual(values = c("red", "yellow", "deepskyblue3", "deepskyblue"),
                    labels = c("Missing", "Fragmented", "Duplicated", "Single")) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Species and Internal Strain Identifier for Assembly")

ggsave("Outputs/internal_busco_plot.png",
       plot = my_buscos,
       width = 8,
       height = length(unique(my.busco.results$ID))*0.5+1)

my_buscos_text <- my_buscos + geom_text(data = subset(my.busco.results, Category == "S"),
                                        aes(label = Written.Score),
                                        y = 0.22, size = 2.5)

ggsave("Outputs/internal_busco_text_plot.png",
       plot = my_buscos_text,
       width = 8,
       height = length(unique(my.busco.results$ID))*0.5+1)

#### Published genomes ####
external.busco.results <- subset(busco.results, Source != "Internal")

external_buscos <- ggplot(data = external.busco.results, aes(x = ID, y = Score)) + 
  geom_bar(aes(fill = Category), position = 'fill', stat = "identity")+
  scale_fill_manual(values = c("red", "yellow", "deepskyblue3", "deepskyblue"),
                    labels = c("Missing", "Fragmented", "Duplicated", "Single")) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Species and Strain Identifier for Assembly")

ggsave("Outputs/external_busco_plot.png",
       plot = external_buscos,
       width = 8,
       height = length(unique(external.busco.results$ID))*0.5+1)

external_buscos_text <- external_buscos + geom_text(data = subset(external.busco.results, Category == "S"),
                                                    aes(label = Written.Score),
                                                    y = 0.22, size = 2.5)

ggsave("Outputs/external_busco_text_plot.png",
       plot = external_buscos_text,
       width = 8,
       height = length(unique(external.busco.results$ID))*0.5+1)

#####