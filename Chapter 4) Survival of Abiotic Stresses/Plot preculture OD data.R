##### Setup #####
#### Packages ####
require(ggplot2)

#### Functions ####
find.counts <- function(d, label.pos = max(OD$`OD600.(calculated)`)){
  return(
    data.frame(
      y = label.pos *1.2 , # Must be named y
      label = paste("n = ", length(d), sep = "") # Why does lenght(df) work?
    )
  )
}

shorten.name <- function(vector.of.names){
  genus.names <- as.data.frame(strsplit(vector.of.names, split = " "))[1,]
  genus.letters <- toupper(substring(genus.names, 1,1))
  species.names <- as.data.frame(strsplit(vector.of.names, split = " "))[2,]
  out <- paste(genus.letters, ". ", species.names, sep = "")
  return(out)
}

#### Factor Levels ####
str.lvls <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/strain_numbers.csv")
sp.lvls <- read.csv("~/Documents/Main_Project/Scripts/update_platereader_data/species_order.csv")

#### Misc ####
setwd("~/Documents/Main_Project/Thesis/Figures/Ch2 Abiotic Stress/Pt3 General Trends/")
OD.data <- "~/Documents/Main_Project/Data/ODplots/ODs.xlsx"
#####

###### 2 day growth ######
#### Load OD data ####
OD <- openxlsx::read.xlsx(OD.data)

OD <- subset(OD, `Incubation.Temp.(°C)` < 30 & 
               Culture.Media == "SAB" &
               `Incubation.Time.(days)` == 2)
#### Assign species names ####

str.lvls <- str.lvls[names(str.lvls) != "Notes"]

names(str.lvls)[names(str.lvls) == "H_number"] <- "Strain"

OD <- merge(OD, str.lvls)

#### Set levels in OD data ####
OD$Strain <- factor(OD$Strain,
                    levels = str.lvls$Strain)
OD$Species <- factor(OD$Species,
                     levels = sp.lvls$Descent)
#####

#### Prepare data for plotting ####
max.strains.per.row <- 18

number.of.rows  <- length(unique(OD$Strain))%/%max.strains.per.row + 1
strains.per.row <- length(unique(OD$Strain))/number.of.rows
strains.per.row <- ceiling(strains.per.row)

ordered.strains <- sort(unique(OD$Strain))
strains.row     <- rep(1:number.of.rows, each = strains.per.row)

for(i in 1:length(ordered.strains)){
  df.rows <- ordered.strains[i] == OD$Strain
  OD$PlotRow[df.rows] <- strains.row[i]
}
#### Plot 2 day growth data ####
upper.y <- max(OD$`OD600.(calculated)`) * 1.3

odplot.b <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
   geom_boxplot(aes(fill = Species)) + # Species, Batch, Infection, Host 
#  geom_point(size = 0.5, position = position_jitter(w = 0.1, h = 0))+
   facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
   # facet_wrap(Species ~ . , scales = "free") +
   stat_summary(
     fun.data = find.counts, 
     geom = "text", 
     hjust = 0.5,
     angle = 45
   ) +
   coord_cartesian(ylim = c(0,upper.y)) +
   theme(
     strip.background = element_blank(),
     strip.text.x = element_blank()
   ) + 
  labs(y = "OD600")

odplot.b +
  ggtitle("ODs of precultures used to inoculate growth curves",
          subtitle = "Precultures grown in: SAB, 25°C, 2 days. Some exceptions apply")

ggsave("day2_PrecultureBox_calcN.png",
       plot = odplot.b,
       height = number.of.rows*5,
       width = 22,
       units = "cm")

odplot.v <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
  geom_violin(aes(fill = Species)) + # Species, Batch, Infection, Host
  facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
  # facet_wrap(Species ~ . , scales = "free") +
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,upper.y)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(y = "OD600")

odplot.v +
  ggtitle("ODs of precultures used to inoculate temperature growth curves",
          subtitle = "Precultures grown in: SAB, 25°C, 2 days. Some exceptions apply")

ggsave("day2_PrecultureViol_calcN.png",
       plot = odplot.v,
       height = number.of.rows*5,
       width = 22,
       units = "cm")

OD$SpeciesSplit <- gsub(" ", "\n", OD$Species)
OD$SpeciesSplit <- factor(OD$SpeciesSplit,
                          levels = gsub(" ", "\n",unique(sort( OD$Species))))

OD$SpeciesShort <- shorten.name(as.character(OD$Species))
OD$SpeciesShort[OD$SpeciesShort=="A. sp.1"] <- "Auxenochlorella sp. 1"
OD$SpeciesShort[OD$SpeciesShort=="A. sp.2"] <- "Auxenochlorella sp. 2"

species.levels <- shorten.name(as.character(sort(unique(OD$Species))))
species.levels[species.levels=="A. sp.1"] <- "Auxenochlorella sp. 1"
species.levels[species.levels=="A. sp.2"] <- "Auxenochlorella sp. 2"
OD$SpeciesShort <- factor(OD$SpeciesShort, 
                          levels = species.levels)

(odplot.species <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = SpeciesShort)) + 
  geom_boxplot(aes(fill = Species)) + # Species, Batch, Infection, Host
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,upper.y)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) + 
  labs(y = "OD600", x = "Species")) 

ggsave("day2_Species_PrecultureBox_calcN.png",
       plot = odplot.species,
       height = 13,
       width = 22,
       units = "cm")
#####

##### Global anova at strain level #####
a.str <- aov(`OD600.(calculated)` ~ Strain, OD)
summary(a.str)

#### Pairwise comparisons, with correction for multiple testing ####
a.str.multi <- TukeyHSD(a.str, which = "Strain")

#### Put comparisons in table ####
### Create matrix to show results ###
strains <- sort(unique(OD$Strain))
rc <- length(strains) -1
out.str <- matrix(nrow = rc,
                  ncol = rc,
                  dimnames = list(strains[1:rc], strains[1:rc+1]),
                  data = NA)

# a.str.multi$Strain[,"p adj"]

# names(a.str.multi$Strain[,"p adj"])[1:5]
### Fill the matrix with data ###
j <- 1
for(i in 1:nrow(out.str)){
  window <- nrow(out.str)-i
  r <- c(rep(NA, i-1), a.str.multi$Strain[,"p adj"][j:(j+window)])
  j <- j+window+1
  out.str[i,] <- r
}

### Subset the matrix, to check data have been positioned correctly ###
n <- 7
choice <- c()
j <- 1
for(i in 1:n){
  window <- 48-i
  choice <- c(choice, j:(j+n-i))
  j <- j+window+1
}
# choice

# out.str[1:n, 1:n]
# a.str.multi$Strain[,"p adj"][choice]

#### Summarise P values ####
### Convert P values to booleans indicating significance ###
out.str.sig <- out.str < 0.05
# out.sig[1:n, 1:n]

#### Present strain P values ####
out.str
out.str.sig
#### Conclusions ####
# Strain has a significant effect on OD (P<2e-16)
# There are significant differences (Padj < 0.05) between some strains. 
#####

##### Anova at strain level for each species #####
species <- unique(OD$Species)

summary.sp.str <- list()
for(sp in species){
  OD.subset <- subset(OD, Species == sp)
  
  if(length(unique(OD.subset$Strain)) == 1) {
    summary.sp.str[[sp]]$Overview <- paste(sp, " has only one strain. ANOVA not performed", sep = "")
    next
    }
  a.sp.str <- aov(`OD600.(calculated)` ~ Strain, OD.subset)
  
  s <- summary(a.sp.str)
  if(s[[1]]$`Pr(>F)`[1] < 0.05){
    summary.sp.str[[sp]]$Overview <- paste("Strains within ", 
                sp, 
                " are statistically different from each other",
                sep = "")
    
    #### Pairwise comparisons, with correction for multiple testing ####
    a.sp.str.multi <- TukeyHSD(a.sp.str, which = "Strain")
    
    
    summary.sp.str[[sp]]$Comparisons <- c("The following comparisons were made:",
                                 attributes(a.sp.str.multi$Strain)$dimnames[[1]])
    
    summary.sp.str[[sp]]$Significant$Comparisons<- c("The following strain comparisons are significant: ",
                                 names(which(a.sp.str.multi$Strain[,4] < 0.05)))
    summary.sp.str[[sp]]$Significant$Padj <- c(a.sp.str.multi$Strain[,4])[a.sp.str.multi$Strain[,4] < 0.05]
    
  } else {
    summary.sp.str[[sp]]$Overview <- paste("Strains within ", 
                sp, 
                " are not statistically different from each other",
                sep = "")
  }
  
}

#### Plot data grouped by species ####
ggplot(data = OD, aes(x = Strain, y = `OD600.(calculated)`))+
  geom_violin(aes(colour = Species)) +
  facet_wrap(Species ~ ., scales = "free")

ggplot(data = OD, aes(x = Strain, y = `OD600.(calculated)`))+
  geom_boxplot(aes(colour = Species)) +
  facet_wrap(Species ~ ., scales = "free_x")
# FIXME: include a counter for how many measurements contribute to each box

#### Look for differences caused by individual strains
pattern <- "HP18"
pattern.species <- str.lvls$Species[str.lvls$Strain == pattern]
tot <- grepl(pattern, summary.sp.str[[pattern.species]]$Comparisons[-1])
sig <- grepl(pattern, summary.sp.str[[pattern.species]]$Significant$Comparisons[-1])
sum(sig)/sum(tot)

#### Conclusions ####
summary.sp.str
## If I only consider a strain to be different from its species if it is different
## in more than half of the comparisons
# HP18 is statistically different from  the other P bovis strains
# the two xanthoriae strains are significantly different from each other
# HP39 is almost statistically significant from the other P ciferrii strains

#### How different is HP18 from the other bovis? ####

subset(OD, Species == "Prototheca bovis") %>%
  group_by(Strain) %>%
  summarise(MeanOD = mean(`OD600.(calculated)`)) %>%
  arrange(MeanOD) -> OD.summary

OD.summary$MeanOD[1]/OD.summary$MeanOD[2] * 100
OD.summary$MeanOD[1]/mean(OD.summary$MeanOD[-1]) * 100


OD %>% group_by(Species) %>%
  summarise(MeanOD = mean(`OD600.(calculated)`))
# FIXME: Need to set a.sp.str.multi to be made with bovis


differences <- a.sp.str.multi$Strain[
  grep("18",dimnames(a.sp.str.multi$Strain)[[1]])&
    a.sp.str.multi$Strain[,4]<0.05,1]
differences[grep("HP18", substring(names(differences),1,5))] <- 
  differences[grep("HP18", substring(names(differences),1,5))] * -1

mean(differences)
#####

##### Anova at species level #####
a.sp <- aov(`OD600.(calculated)` ~ Species, OD)
summary(a.sp)
#### Correct for multiple testing ####
# is the tukey test a pairwise test or a way to correct for multiple testing, or both?
# are there better tests?
a.sp.multi <- TukeyHSD(a.sp, which = "Species")

#pairwise.t.test(OD$`OD600.(calculated)`, OD$Species, p.adjust.method = "bonferro")
## Tukey's is probably better than many t-tests
#### Put results in table ####
### Create matrix to show results ###
species <- sort(unique(OD$Species))
rc <- length(species) -1
out.sp <- matrix(nrow = rc,
                 ncol = rc,
                 dimnames = list(species[1:rc], species[1:rc+1]),
                 data = NA)



# names(a.str.multi$Strain[,"p adj"])[1:5]
### Fill the matrix with data ###

comparisons <- dimnames(a.sp.multi$Species)[[1]]
for(i in 1:length(comparisons)){
  col <- strsplit(comparisons[i], split = "-")[[1]][1]
  row <- strsplit(comparisons[i], split = "-")[[1]][2]
  out.sp[row,col] <- a.sp.multi$Species[,"p adj"][i]
}

#### Summarise P values ####
### Convert P values to booleans indicating significance ###
out.sp.sig <- out.sp < 0.05

### Convert P values to stars indicating significance ###

### Convert P values to stars indicating significance ###
# FIXME : Matrix isn't invertable 
out.sp.star <- matrix(nrow = nrow(out.sp),
                      ncol = ncol(out.sp),
                      data = " ")
out.sp.star[out.sp < 0.05] <- "*"
out.sp.star[out.sp < 0.01] <- "**"
out.sp.star[out.sp < 0.001] <- "***"
out.sp.star[is.na(out.sp)] <- NA
dimnames(out.sp.star) <- dimnames(out.sp)
#### Present species P values ####
out.sp
out.sp.sig
out.sp.star

write.csv(as.data.frame(out.sp.star),
          "OD_SpeciesComparisons_Significance.csv",
          row.names = T)
#####


#### Anova at lineage level ? ####


###### 5 day growth ######
#### Load OD data ####
OD <- openxlsx::read.xlsx(OD.data)

OD <- subset(OD, `Incubation.Temp.(°C)` < 30 & 
               Culture.Media == "SAB" &
               `Incubation.Time.(days)` == 5)
#### Assign species names ####

str.lvls <- str.lvls[names(str.lvls) != "Notes"]

names(str.lvls)[names(str.lvls) == "H_number"] <- "Strain"

OD <- merge(OD, str.lvls)

#### Set levels in OD data ####
OD$Strain <- factor(OD$Strain,
                    levels = str.lvls$Strain)
OD$Species <- factor(OD$Species,
                     levels = sp.lvls$Descent)
#### Prepare data for plotting ####
max.strains.per.row <- 18

number.of.rows  <- length(unique(OD$Strain))%/%max.strains.per.row + 1
strains.per.row <- length(unique(OD$Strain))/number.of.rows
strains.per.row <- ceiling(strains.per.row)

ordered.strains <- sort(unique(OD$Strain))
strains.row     <- rep(1:number.of.rows, each = strains.per.row)

for(i in 1:length(ordered.strains)){
  df.rows <- ordered.strains[i] == OD$Strain
  OD$PlotRow[df.rows] <- strains.row[i]
}
#### Plot 5 day growth data ####
upper.y.5 <- max(OD$`OD600.(calculated)`)*1.3

odplot.b.5 <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
  geom_boxplot(aes(fill = Species)) + # Species, Batch, Infection, Host
  facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
  # facet_wrap(Species ~ . , scales = "free") +
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,upper.y.5)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(y = "OD600")

odplot.b.5 +
  ggtitle("ODs of precultures used to inoculate growth curves",
          subtitle = "Precultures grown in: SAB, 25°C, 5 days")

ggsave("day5_PrecultureBox_calcN.png",
       plot = odplot.b.5,
       height = number.of.rows*5,
       width = 22,
       units = "cm")

odplot.v.5 <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
  geom_violin(aes(fill = Species)) + # Species, Batch, Infection, Host
  facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
  # facet_wrap(Species ~ . , scales = "free") +
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,26)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(y = "OD600")

odplot.v.5 +
  ggtitle("ODs of precultures used to inoculate temperature growth curves",
          subtitle = "Precultures grown in: SAB, 25°C, 5 days")

#ggsave("day5_PrecultureViol_calcN.png",
#       plot = odplot.v.5,
#       height = number.of.rows*5,
#       width = 22,
#       units = "cm")
#####


###### YPD growth ######
#### Load OD data ####
OD <- openxlsx::read.xlsx(OD.data)

OD <- subset(OD, `Incubation.Temp.(°C)` < 30 & 
               Culture.Media == "YPD" &
               `Incubation.Time.(days)` == 2)
#### Assign species names ####

str.lvls <- str.lvls[names(str.lvls) != "Notes"]

names(str.lvls)[names(str.lvls) == "H_number"] <- "Strain"

OD <- merge(OD, str.lvls)

#### Set levels in OD data ####
OD$Strain <- factor(OD$Strain,
                    levels = str.lvls$Strain)
OD$Species <- factor(OD$Species,
                     levels = sp.lvls$Descent)

#### Prepare data for plotting ####
max.strains.per.row <- 18

number.of.rows  <- length(unique(OD$Strain))%/%max.strains.per.row + 1
strains.per.row <- length(unique(OD$Strain))/number.of.rows
strains.per.row <- ceiling(strains.per.row)

ordered.strains <- sort(unique(OD$Strain))
strains.row     <- rep(1:number.of.rows, each = strains.per.row)

for(i in 1:length(ordered.strains)){
  df.rows <- ordered.strains[i] == OD$Strain
  OD$PlotRow[df.rows] <- strains.row[i]
}
#### Plot YPD 2 day growth data ####
upper.y.y <- max(OD$`OD600.(calculated)`)*1.3

odplot.b.y <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
  geom_boxplot(aes(fill = Species)) + # Species, Batch, Infection, Host
  facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
  # facet_wrap(Species ~ . , scales = "free") +
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,upper.y.y)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(y = "OD600")

odplot.b.y +
  ggtitle("ODs of precultures used to inoculate growth curves",
          subtitle = "Precultures grown in: YPD, 25°C, 2 days")

ggsave("day2YPD_PrecultureBox_calcN.png",
       plot = odplot.b.y,
       height = max(c(length(unique(OD$Species))*0.75,number.of.rows*5)),
       width = 22,
       units = "cm")

odplot.v.y <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = Strain)) + 
  geom_violin(aes(fill = Species)) + # Species, Batch, Infection, Host
  facet_wrap(. ~ PlotRow , scales = "free", ncol = 1) +
  # facet_wrap(Species ~ . , scales = "free") +
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,26)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  labs(y = "OD600")

odplot.v.y +
  ggtitle("ODs of precultures used to inoculate temperature growth curves",
          subtitle = "Precultures grown in: YPD, 25°C, 2 days")

#ggsave("day2YPD_PrecultureViol_calcN.png",
#       plot = odplot.v.y,
#       height = number.of.rows*5,
#       width = 22,
#       units = "cm")

OD$SpeciesShort <- shorten.name(as.character(OD$Species))
OD$SpeciesShort <- factor(OD$SpeciesShort, 
                          levels = shorten.name(as.character(sort(unique(OD$Species)))))

odplot.b.y.sp <- ggplot(data = OD, aes(y=`OD600.(calculated)`, x = SpeciesShort)) + 
  geom_boxplot(aes(fill = Species)) + # Species, Batch, Infection, Host
  stat_summary(
    fun.data = find.counts, 
    geom = "text", 
    hjust = 0.5,
    angle = 45
  ) +
  coord_cartesian(ylim = c(0,upper.y.y)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  ) + 
  labs(y = "OD600")

ggsave("day2YPD_Species_PrecultureBox_calcN.png",
       plot = odplot.b.y.sp,
       height = 11,
       width = 22,
       units = "cm")
#####
