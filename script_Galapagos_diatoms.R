## Code "Ecology and diversity of diatoms in Galápagos Islands
## Bachelor final degree (University of Barcelona)
## Authors: Cristina Alcaide and Xavier Benito (xavier.benito.granell@gmail.com)
## September 2024

#Clear workspace
rm(list=ls(all=TRUE))

#Set WD
setwd("")

##loading libraries for functions used
library(vegan) #to perform nonmetric multidimensional analysis, simper and adonis2 
library(ade4) #to perfom PCA analysis
library(ggplot2) #to make nice ordination plots
library(cowplot) #to plot differents ggplots objects in the same grid
library(ggrepel) # to repel overlapping text in ggplots
library(ggpubr) #to draw ellipses
library(goeveg) #allow to select species for vegan ordination objects
library(tidyverse) #allow to manipulate tabulate data
library(scales) #allow to scale variables from 0 to 1 
library(psych) #correlation tests
library(svglite) #to export .svg 
library(egg)

######################################
#### Read and manipulate datasets ####
######################################

# Read in meta data--site information
meta <- read.csv("meta.csv", row.names = 1) %>%
  mutate(island=factor(island))
str(meta)

# Read in Environmental data
env_data <- read.csv("Galapagos_env.csv", row.names = 1)

# Look at the data structure
str(env_data)
colSums(is.na(env_data)) #check how many observations are NAs

# Here we extract the spatial variables to use for later
spatial_var <- env_data[,c(27,28,29)]
names(spatial_var) <- c("area_island","area_habitat","elevation")
spatial_var$area_island <- as.numeric(spatial_var$area_island) #transform area to numeric
spatial_var$elevation <- as.numeric(spatial_var$elevation) #transform area to numeric

# Remove variables we don't want (i.e., variables with all NA, area island and habitat and elevation)
env_data <- env_data[,-c(4,5,12,14,18,19,20,22,23,24,25,27,28,29)] 

# Remove samples with all NA 
env_subset <- env_data %>%
  filter_all(any_vars(!is.na(.))) 
str(env_subset)

# See how many variables have NA observations still
colSums(is.na(env_subset))

#panel correlation plots to assess data distribution
panel.hist <- function(x, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(usr[1:2], 0, 1.5) )     
  h <- hist(x, plot = FALSE)     
  breaks <- h$breaks; nB <- length(breaks)     
  y <- h$counts; y <- y/max(y)     
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) 
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {     
  usr <- par("usr"); on.exit(par(usr))     
  par(usr = c(0, 1, 0, 1))     
  r <- abs(cor(x, y, use = "complete"))   
  txt <- format(c(r, 0.123456789), digits = digits)[1]     
  txt <- paste0(prefix, txt)     
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)     
  text(0.5, 0.5, txt, cex = cex.cor * r) }


# Make Correlation plot
pairs(env_subset, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

# Remove outliers 
env_subset$Water.T[env_subset$Water.T > 39] <- NA

#transform variables to meet assumptions of homogenity of variances
env_transformed <- transform(env_subset, Water.T=log10(Water.T+0.25), 
                             pH=log10(pH+0.25), Cond= log10(Cond+0.25), 
                             Secchi=log10(Secchi+0.25), Alkalinity=log10(Alkalinity+0.25), 
                             Ca=log10(Ca+0.25), Mg=log10(Mg+0.25), K=log10(K+0.25), Na=log10(Na+0.25), 
                             Cl=log10(Cl+0.25), NO3=log(NO3+0.25), SO4=log10(SO4+0.25), PO4=log10(PO4+0.25), 
                             DO=log10(DO+0.25), MaxDepth=log10(MaxDepth+0.25))

# re-run the corrlation plot
pairs(env_transformed, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

# Merge variables and metadata
all_data <- merge(env_transformed, meta, by="row.names")
row.names(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

# boxplots environmental variation across islands
svglite("boxplot_env_island.svg", width = 10, height = 8)
par(oma=c(2,0,0,2))
par(mfrow=c(5,4))
nms <- colnames(all_data[,1:15])
for (i in 1:length(nms)) {
  bp <- boxplot(all_data[, i] ~ island, data = all_data, ylab=NULL, xlab=NULL, main=nms[i])
  tick <- seq_along(bp$names)
  axis(1, at = tick, labels = FALSE)
}
dev.off()

# boxplots environmental variation across habitat types
par(oma=c(2,0,0,2))
par(mfrow=c(5,4))
for (i in 1:length(nms)) {
  bp <- boxplot(all_data[, i] ~ habitat, data = all_data, ylab=NULL, xlab=NULL, main=nms[i])
  tick <- seq_along(bp$names)
  axis(1, at = tick, labels = FALSE)
}

######################################
#### Principal Component Analysis ####
######################################

# re-assign dataframe for PCA
PCA_data <- env_transformed

#Perform PCA with NIPALS algorithm
head(PCA_data)
colSums(is.na(PCA_data)) #check how many observations are NAs

# Remove variables with more than 50% of NAs
PCA_data <- PCA_data %>%
  dplyr::select(-c(8,9,14))

PCA.data <- data.frame(scale(PCA_data)) #scale variables
PCA.nipals <- nipals(PCA.data, nf = 4, rec = FALSE, niter = 1000, tol = 1e-09)

#Save summary matrix results  
PCA.summary <- NULL
PCA.summary <- matrix(data = NA, nrow = 8, ncol = 3, byrow = FALSE, dimnames = NULL)
colnames(PCA.summary) <- c("Parameter", "Value", "Explained variance")

PCA.summary[1, 1] <- c("Number of extracted factors")
PCA.summary[1, 2] <- PCA.nipals$nf

PCA.summary[2, 1] <- c("Number of variables")
PCA.summary[2, 2] <- nrow(PCA.nipals$co)

PCA.summary[3, 1] <- c("Eigenvalue 1")
PCA.summary[3, 2] <- PCA.nipals$eig[1]
PCA.summary[3, 3] <- round((PCA.nipals$eig[1] / nrow(PCA.nipals$co)) * 100,2)

PCA.summary[4, 1] <- c("Eigenvalue 2")
PCA.summary[4, 2] <- PCA.nipals$eig[2]
PCA.summary[4, 3] <- round((PCA.nipals$eig[2] / nrow(PCA.nipals$co)) * 100,2)

PCA.summary[5, 1] <- c("Eigenvalue 3")
PCA.summary[5, 2] <- PCA.nipals$eig[3]
PCA.summary[5, 3] <- round((PCA.nipals$eig[3] / nrow(PCA.nipals$co)) * 100,2)

PCA.summary[6, 1] <- c("Eigenvalue 4")
PCA.summary[6, 2] <- PCA.nipals$eig[4]
PCA.summary[6, 3] <- round((PCA.nipals$eig[4] / nrow(PCA.nipals$co)) * 100,2)

PCA.summary[7, 1] <- c("Number of iterations axis 1")
PCA.summary[7, 2] <- PCA.nipals$nb[1]

PCA.summary[8, 1] <- c("Number of iterations axis 2")
PCA.summary[8, 2] <- PCA.nipals$nb[2]

#Save the column coordinates (Component matrix) #ncol= n variables 
Component.matrix <- NULL
Component.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 5, byrow = FALSE, dimnames = NULL)
Component.matrix[,1] <- rownames(PCA.nipals$co)
Component.matrix[,2] <- (((PCA.nipals$co[,1] ^ 2) * PCA.nipals$eig[1]) / PCA.nipals$eig[1]) + (((PCA.nipals$co[,2] ^ 2) * PCA.nipals$eig[2]) / PCA.nipals$eig[2])

Component.matrix <- cbind(Component.matrix[,1], PCA.nipals$co, Component.matrix[,2])
colnames(Component.matrix) <- c("Variable", "Component 1", "Component 2", "Component 3", "Component 4",  "Power Importance")

#Save the column normed scores (Component scores coefficient matrix)
Component.coefficient.matrix <- NULL
Component.coefficient.matrix <- matrix(data = NA, nrow = ncol(PCA.data), ncol = 1, byrow = FALSE, dimnames = NULL)
Component.coefficient.matrix[,1] <- rownames(PCA.nipals$c1)
Component.coefficient.matrix <- cbind(Component.coefficient.matrix, PCA.nipals$c1)
colnames(Component.coefficient.matrix) <- c("Variable", "Component 1", "Component 2", "Component 3", "Component 4")

#Save the row coordinates (Factor Scores)
colnames(PCA.nipals$li) <- c("Component 1", "Component 2", "Component 3", "Component 4")
Factor.scores <- data.frame(cbind(PCA.data, PCA.nipals$li))

#Create data frame with site scores and regions
PCA.scores <- merge(Factor.scores, meta, by="row.names")
row.names(PCA.scores) <- PCA.scores$Row.names
PCA.scores$Row.names <- NULL

PCA.result <- data.frame(PCA1=PCA.scores$Component.1, PCA2=PCA.scores$Component.2,
                         island=PCA.scores$island, habitat=PCA.scores$habitat)

#extract factor scores (environmental variables)
comp1 <- as.numeric(Component.coefficient.matrix[,2])
comp2 <- as.numeric(Component.coefficient.matrix[,3])

## Plot PCA
# Export plots. First we need to call the function. SVG is a vector file that can be modified with a third party software (e.g., Adobe, Inkscape)
#svglite("PCA.svg", width = 10, height = 8)

par(mfrow=c(1,2))
par(mar=c(5,4,3,3)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text. 

plot(PCA.result$PCA1, PCA.result$PCA2, type = "n", xlab=paste("PCA1","(",PCA.summary[3,3],"%",")"), ylab=paste("PCA2","(",PCA.summary[4,3],"%",")"))

title("Sites")
abline(h=0, col="grey")
abline(v=0, col="grey")

# Transform island and habitat to factor for plotting
PCA.result$island <- as.factor(PCA.result$island)
PCA.result$habitat <- as.factor(PCA.result$habitat)

# Manual plotting
points(PCA.result[(PCA.result$island=="Santiago" & PCA.result$habitat=="Lagoon"), 1:2], col="blue", pch=15, cex=1.5)
points(PCA.result[(PCA.result$island=="Santiago" & PCA.result$habitat=="Lake"), 1:2], col="blue", pch=16,cex=1.5)
points(PCA.result[(PCA.result$island=="Santiago" & PCA.result$habitat=="Pond"), 1:2], col="blue", pch=17,cex=1.5)

points(PCA.result[(PCA.result$island=="Genovesa" & PCA.result$habitat=="Lagoon"), 1:2], col="orange", pch=15,cex=1.5)
points(PCA.result[(PCA.result$island=="Genovesa" & PCA.result$habitat=="Lake"), 1:2], col="orange", pch=16,cex=1.5)
points(PCA.result[(PCA.result$island=="Genovesa" & PCA.result$habitat=="Pond"), 1:2], col="orange", pch=17,cex=1.5)

points(PCA.result[(PCA.result$island=="Isabela" & PCA.result$habitat=="Lagoon"), 1:2], col="red", pch=15,cex=1.5)
points(PCA.result[(PCA.result$island=="Isabela" & PCA.result$habitat=="Lake"), 1:2], col="red", pch=16,cex=1.5)
points(PCA.result[(PCA.result$island=="Isabela" & PCA.result$habitat=="Pond"), 1:2], col="red", pch=17,cex=1.5)

points(PCA.result[(PCA.result$island=="San Cristóbal" & PCA.result$habitat=="Lagoon"), 1:2], col="forestgreen", pch=15,cex=1.5)
points(PCA.result[(PCA.result$island=="San Cristóbal" & PCA.result$habitat=="Lake"), 1:2], col="forestgreen", pch=16,cex=1.5)
points(PCA.result[(PCA.result$island=="San Cristóbal" & PCA.result$habitat=="Pond"), 1:2], col="forestgreen", pch=17,cex=1.5)

points(PCA.result[(PCA.result$island=="Santa Cruz" & PCA.result$habitat=="Lagoon"), 1:2], col="darkviolet", pch=15,cex=1.5)
points(PCA.result[(PCA.result$island=="Santa Cruz" & PCA.result$habitat=="Lake"), 1:2], col="darkviolet", pch=16,cex=1.5)
points(PCA.result[(PCA.result$island=="Santa Cruz" & PCA.result$habitat=="Pond"), 1:2], col="darkviolet", pch=17,cex=1.5)

#legend for monhtly dataset
legend("bottomleft", c("Santiago", "Genovesa", "Isabela", "San Cristóbal", "Santa Cruz",
                       "Lagoon", "Lake", "Pond")
       , cex=0.9, pch=c(NA,NA, NA, NA, NA, 15,16,17),
       text.col=c("blue", "orange", "red","forestgreen","darkviolet","black","black","black"), 
       #col=c("blue", "orange", "red","forestgreen","darkgrey", "bl"),
       ncol = 1, xpd = TRUE)

# Add lagoon names 
#text(PCA.result[,1:2], labels = PCA.result$island, pos = 1, cex = 0.5, offset = 0.2)

#Variables
comp1 <- as.numeric(Component.coefficient.matrix[,2])
comp2 <- as.numeric(Component.coefficient.matrix[,3])

#Labels
labels <- as.character(Component.coefficient.matrix[,1])

#par(mar=c(2,2,4,4)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text.
plot(comp1, comp2, pch=16, col="black", xlab=paste("PCA1","(",PCA.summary[3,3],"%",")"), ylab=paste("PCA2","(",PCA.summary[4,3],"%",")"))
title("Variables")
abline(h=0, col="grey")
abline(v=0, col="grey")

text(comp1, comp2, labels = labels, pos = 3, cex = 0.8, offset = 0.3)

# Run this to save the plot
dev.off()

#Correlations between PCA components and factor variables
# Check if new recoded variables "island" and "habitat" have significant correlations with PCA axes. 
PCA.scores <- PCA.scores %>% mutate(island=dplyr::recode(island,
                                               "Santiago"=1,
                                               "Isabela"=2,
                                               "San Cristóbal"=3,
                                               "Santa Cruz"=4,
                                               "Genovesa"=5)) %>%
  mutate(habitat=dplyr::recode(habitat,"Lagoon"=1, "Lake"=2, "Pond"=3)) 


cor(PCA.scores[,c("Component.1","Component.2")], PCA.scores[,c(1:12,17:18)], use = "pairwise.complete.obs")
cor.r <- round(corr.test(PCA.scores[,c("Component.1","Component.2")],  PCA.scores[,c(1:12,17:18)], method = "spearman")$r,2)
cor.p <- round(corr.test(PCA.scores[,c("Component.1","Component.2")],  PCA.scores[,c(1:12,17:18)], method = "spearman")$p,2)

cor_PCA_results <- rbind(cor.r, cor.p)
row.names(cor_PCA_results) <- c("PCA1r", "PCA2r", "PCA1p", "PCA2p")

# Export the correlation table
#write.table(cor_PCA_results, "cor_PCA_results.txt")

###################
# Diatom analyses #
###################

# Read in diatom data (presence/absence)
diat <- read.csv("Galapagos_diat_PA.csv", row.names = 1)
str(diat)

## NMDS
##Remove rare species
n.occur <- apply(diat>0, 2, sum)
diat <- diat[, n.occur>1] #We keep species present in more than 2 samples

#Run NMDS with presence/absence data and Jaccard dissimilarity index
diat.nmds <- metaMDS(diat, noshare = TRUE, trymax = 500, distance = "jaccard")

#Extract scores from diatom NMDS ordination
diat.nmds.scores <- as.data.frame(scores(diat.nmds, display = "sites"))

## Here we merge nmds scores with env data and metadata
# First we assign row.names (ID) as in environmental data
row.names(diat.nmds.scores) <- row.names(env_data)
diat_env_nmds <- merge(diat.nmds.scores, all_data, by="row.names")
row.names(diat_env_nmds) <- diat_env_nmds$Row.names
diat_env_nmds$Row.names <- NULL

#Plot NMDS
##create a result table to facilitate plotting afterwards
nmds_tbl <- data.frame(NMDS1=diat_env_nmds$NMDS1, NMDS2=diat_env_nmds$NMDS2, 
                       island=as.factor(diat_env_nmds$island), habitat=as.factor(diat_env_nmds$habitat), site=diat_env_nmds$site)

# Make the plot
nmds_plt <- ggplot(nmds_tbl, aes(NMDS1,NMDS2, colour=island, shape=habitat, label=site)) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(size=4) +
  #geom_point(aes(shape = region), size=3) +
  #stat_conf_ellipse(aes(x=NMDS1, y=NMDS2, color=habitat, type="norm")) +
  #scale_colour_manual(values=c("#E69F00", "#999999")) +
  geom_text_repel(colour="black", size=3) +
  geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  theme_article() +
  theme(legend.text = element_text(size=10),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.box = "vertical") 
nmds_plt

# This is to select, for instance, the 10% most abundant species with 90% best environmental fit in NDMS for axes 1 & 2
# the object nmds is your model, and diat is the dataframe used to perform the nmds
selected_nmds <- ordiselect(diat, diat.nmds,  ablim = 0.2, fitlim = 0.7, method = "axes")  

#Extract NMDS species scores. Ideally, the species code should be placed in the row.names to keep these id all the time
scrs.spp <- data.frame(scores(diat.nmds, display = "species", choices = 1:2))

#Select NMDS species scores from the ordiselect
selected_spp_scrs <- scrs.spp[row.names(scrs.spp) %in% selected_nmds, ]

# Here we keep species label that are explaining a significant portion of variability on the NMDS axes, and then drop those from the larger list
keep <- row.names(scrs.spp) %in% row.names(selected_spp_scrs)
scrs.spp.drop <- scrs.spp[!c(keep),]

# Plot diatom NMDS scores
spp_plt <- ggplot(data=selected_spp_scrs, aes(x=NMDS1, y=NMDS2), colour="darkgrey", alpha=0.6) +
  geom_point(size=3) +
  xlab("NMDS1") + ylab("NMDS2") +
  geom_text_repel(data=selected_spp_scrs, aes(x=NMDS1, y=NMDS2, label=row.names(selected_spp_scrs)), colour="black", size=3) +
  geom_vline(aes(xintercept = 0), linetype = "solid", colour="grey") +
  geom_hline(aes(yintercept = 0), linetype = "solid", colour="grey") +
  #coord_fixed() + ## need aspect ratio of 1!
  theme_bw()
spp_plt

#cowplot's plot.grid function to join the two NMDS plots: sites and species
library(cowplot)
nmds_plt <- plot_grid(nmds_plt, spp_plt, 
                      ncol=2, axis = "rlbt",
                      labels=c("a","b")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) 
nmds_plt

# Save the plot
ggsave("NMDS_plot_2.svg", nmds_plt)

# Correlations between NMDS1 and 2 and environmental variables
diat_env_nmds <- diat_env_nmds %>% mutate(island=dplyr::recode(island,
                                                        "Santiago"=1,
                                                        "Isabela"=2,
                                                        "San Cristóbal"=3,
                                                        "Santa Cruz"=4,
                                                        "Genovesa"=5)) %>%
  mutate(habitat=dplyr::recode(habitat,"Lagoon"=1, "Lake"=2, "Pond"=3)) 


cor(diat_env_nmds[,c("NMDS1","NMDS2")], diat_env_nmds[,c(3:19)], use = "pairwise.complete.obs")
cor.r <- round(corr.test(diat_env_nmds[,c("NMDS1","NMDS2")], diat_env_nmds[,c(3:19)], method = "spearman")$r,2)
cor.p <- round(corr.test(diat_env_nmds[,c("NMDS1","NMDS2")], diat_env_nmds[,c(3:19)], method = "spearman")$p,2)

cor_nmds_results <- rbind(cor.r, cor.p)
row.names(cor_nmds_results) <- c("NMDS1r", "NMDS2r", "NMDS1p", "NMDS2p")
cor_nmds_results

## Assign row.names (ID) as in environmental data
row.names(diat) <- row.names(env_data)

## Species richness-area relationships
rich <- apply(diat>0, 1, sum)

# test species richness vs spatial variables vs some environmental variables
rich_space <- cbind(rich,spatial_var,meta[,c("habitat", "island", "site")],
                    env_data[,c("pH","Cond","Secchi","Ca")])
names(rich_space) <- c("rich","area_island","area_habitat", "elevation", "habitat", "island", "site","pH","Cond","Secchi", "Ca")
row.names(rich_space) <- row.names(meta)

habitats_keep <- c("Lake","Pond") #create a vector with the habitats wanted to keep
row.names(meta) <- meta$site
habitats_keep_df <- meta[meta$habitat %in% habitats_keep,] #extract lake and pond habitat meta information

rownames_keep <- row.names(habitats_keep_df)
rich_space <- rich_space[(row.names(rich_space) %in% rownames_keep),]
names(rich_space)

# explorative plots
par(mfrow=c(1,2))
par(mar=c(5,4,3,3)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text. 
plot(log10(rich_space$area_island),log10(rich_space$rich))
text(log10(rich_space$area_island), log10(rich_space$rich), labels = rich_space$site, pos = 3, cex = 0.5, offset = 0.3)
plot(log10(rich_space$area_habitat),log10(rich_space$rich))
text(log10(rich_space$area_habitat), log10(rich_space$rich), labels = rich_space$site, pos = 3, cex = 0.5, offset = 0.3)

rich_space <- rich_space[!(row.names(rich_space) %in% "Genovesa"),] #genovesa site is dragging the dispersion. Let´s try to remove it

#linear regressions with base R
mod1 <- lm(log10(rich)~log10(area_habitat), data = rich_space) #Create the linear regression
plot(log10(rich)~log10(area_habitat), pch = 16, col = "blue", data=rich_space) #Plot the results
abline(mod1) #Add a regression line
summary(mod1)

mod2 <- lm(log10(rich)~log10(area_island), data = rich_space) #Create the linear regression
plot(log10(rich)~log10(area_island), pch = 16, col = "blue", data=rich_space) #Plot the results
abline(mod2) 

plot(mod1$residuals)

# Correlations with Pearson
set.seed(1)
cor.r <- round(corr.test(rich_space[,c("rich")], rich_space[,(c("area_island", "area_habitat", "elevation"))], method = "spearman")$r,2)
cor.p <- round(corr.test(rich_space[,c("rich")], rich_space[,(c("area_island", "area_habitat", "elevation"))], method = "spearman")$p,2)
cor_richness <- rbind(cor.r, cor.p)
row.names(cor_richness) <- c("corr.r", "corr.p")
cor_richness

# plot with ggplot
rich_island <- ggplot(rich_space, aes(x=log10(area_island), y=log10(rich), label=site)) +
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x) +
  geom_point(size=4) +
  stat_cor(label.y = 2) + #here we add Pearson's Rho and p value of the correlation
  xlab("log area island (km2)") + ylab("log species richness") +
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 
rich_island

rich_habitat <- ggplot(rich_space, aes(x=log10(area_habitat), y=log10(rich), label=site)) +
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x) +
  geom_point(size=4) +
  stat_cor(label.x = 1,label.y = 2) + #here we add Pearson's Rho and p value of the correlation
  xlab("log area habitat (km2)") + ylab("log species richness") +
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 
rich_habitat

rich_elevation <- ggplot(rich_space, aes(x=sqrt(elevation), y=log10(rich), label=site)) +
  geom_smooth(method = "lm", se=TRUE, formula = y ~ x) +
  geom_point(size=4) +
  stat_cor(label.y = 2) + #here we add Pearson's Rho and p value of the correlation
  xlab("sqrt (elevation)") + ylab("log species richness") +
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 

# save plots individually
ggsave("richness_island_plot.svg", rich_island,
       width = 10, height = 8)
ggsave("richness_habitat_plot.svg", rich_habitat,
       width = 10, height = 8)
ggsave("richness_elevation_plot.svg", rich_elevation,
       width = 10, height = 8)

# Save plots as a grouped template (optionally)
library(cowplot)
richness_plt <- plot_grid(rich_island, rich_habitat, rich_elevation,
                      ncol=3, axis = "rlbt",
                      labels=c("a","b","c")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(legend.title=element_blank())
richness_plt

ggsave("richess_areas_plot.svg", richness_plt,
       width = 12, height = 8)


#################
## MANTEL TEST ##
#################

#Compute disimilarity distances: geographical distances among sites, diatom speciesd (Jaccard) and environmental (using the four PCA axes)
# Because we need dataframes with the exact number of rows, the primary matrix is the PCA.scores and we join coordinate and diatom dataframe to match it
row.names(PCA.scores) <- PCA.scores$site
diat <- read.csv("Galapagos_diat_PA.csv", row.names = 1)

spp.env <- merge(diat, PCA.scores, by = "row.names")
names(spp.env)

spp <- decostand(spp.env[,2:214], method = "hellinger")
diat.dist <- analogue::distance(spp, method="bray")

# create an indicator for all diagonals in the matrix
d <- row(diat.dist) - col(diat.dist)
# use split to group on these values
diagonal<-split(diat.dist, d)
diat.dist.long<-unlist(diagonal["1"]) # select relevant one (diag = 1)

sites.dist <- analogue::distance(cbind(PCA.scores$lon, PCA.scores$lat), method="euclidean")
# create an indicator for all diagonals in the matrix
d <- row(sites.dist) - col(sites.dist)
# use split to group on these values
diagonal<-split(sites.dist, d)
sites.dist.long<-unlist(diagonal["1"]) # select relevant one (diag = 1)

env.dist <- analogue::distance(PCA.scores[,c("Component.1", "Component.2", "Component.3", "Component.4")], method="euclidean")
# create an indicator for all diagonals in the matrix
d <- row(env.dist) - col(env.dist)
# use split to group on these values
diagonal<-split(env.dist, d)
env.dist.long<-unlist(diagonal["1"]) # select relevant one (diag = 1)

# plot(env.dist,diat.dist)
# plot(sites.dist,diat.dist)

# Create a dataframe for plotting with ggplot
dist_all <- data.frame(diat=diat.dist.long*-1, env=env.dist.long, sites=sites.dist.long, 
                       id=all_data$site[-1], habitat=all_data$habitat[-1], 
                       spp.env[,215:226][-1,])

# plots
env_sites_plt <- ggplot(dist_all, aes(x=sites, y=env, label=id)) +
  geom_smooth(method = "lm", se=TRUE) +
  geom_point(size=3) +
  xlab("Spatial distance (Euclidean)") + ylab("Environmental distance (Euclidean)")+
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 
env_sites_plt

diat_env_plt <- ggplot(dist_all, aes(x=env, y=diat, label=id)) +
  geom_smooth(method = "lm", se=TRUE) +
  geom_point(size=3) +
  xlab("Environmental distance (Euclidean)") + ylab("Diatom community similarity (Jaccard)")+
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 
diat_env_plt

diat_sites_plt <- ggplot(dist_all, aes(x=sites, y=diat, label=id)) +
  geom_smooth(method = "lm", se=TRUE) +
  geom_point(size=3) +
  xlab("Spatial distance (Euclidean)") + ylab("")+
  geom_text_repel(colour="black", size=4) +
  theme_article() +
  theme(axis.text=element_text(size=12),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14)) 
diat_sites_plt

# save plots individually
ggsave("env_sites_mantel.svg", env_sites_plt,
       width = 10, height = 8)
ggsave("diat_env_mantel.svg", diat_env_plt,
       width = 10, height = 8)
ggsave("diat_sites_mantel.svg", diat_sites_plt,
       width = 10, height = 8)

# Save plots as a grouped template (optionally)
ddr_plt <- plot_grid(diat_env_plt, diat_sites_plt,
                          ncol=2, axis = "rlbt",
                          labels=c("a","b")) +
  theme(plot.margin = unit(c(0,0,0,0), "cm")) +
  theme(legend.title=element_blank())

ggsave("ddr_plot.svg", ddr_plt,
       width = 10, height = 8)

ggsave("ddr_plot.png", ddr_plt,
       width = 10, height = 8)

#Mantel test using vegan's mantel  function
set.seed(1) #

mantel(env.dist, sites.dist) #correlation dissimilarities between enviroment and space
mantel(diat.dist, sites.dist) #correlation dissimilarities between diatoms and space
mantel(diat.dist, env.dist) #correlation dissimilarities between diatoms and environment

# Partial Mantel
# What's the effect of environmental distance gradients on diatom community while factoring out spatial distances?
mantel.partial(diat.dist, env.dist, sites.dist, method="spearman")

## Most frequent diatom taxa per island and habitat
df_spp <- merge(diat,meta[,c("habitat","island","site")], by="row.names")
df_spp$Row.names <- NULL

spp_most_frequent <- df_spp %>% 
  gather(taxa, freq, -island, -habitat, -site) %>%
  group_by(taxa,island,habitat) %>%
  summarise(freq = sum(freq)) %>%
  filter(!freq == "0" ) %>% #this is to remove empty samples (rows)
  ungroup() %>%
  group_by(island,habitat) %>%
  mutate(occurrence_perc = freq / sum(freq) * 100) %>%
  slice(which.max(occurrence_perc))

write.csv(spp_most_frequent, "spp_most_frequent.csv")
