library(BIOMASS)
library(brranching)
library(ape)
library(phylobase)
library(picante)
library(phytools)
setwd("//NASLSB/Base de données/BD Cirad")

# data <- read.csv("Extraction_base_phy_meca-Essais meca.csv", header=T, sep=";", dec=".")
# 
# guyane <- data[data$LIB_PAYS == "Guyane Francaise",]
# guyane$taxon <- paste(guyane$genus, guyane$sp, sep=" ")
# taxo_cor <- correctTaxo(guyane$taxon, useCache=F)
# 
# guyane$genus <- taxo_cor$genusCorrected 
# guyane$sp <- taxo_cor$speciesCorrected 
# 
# family_result <- getTaxonomy(guyane$genus)
# guyane$family <- family_result$family
# 
# write.csv(guyane, "Essais meca_guyane.csv")


guyane <- read.csv("Essais meca_guyane.csv", header=T, sep=";", dec=".")
guyane$family <- as.factor(guyane$family)

guyane_m <- aggregate(RB  ~ taxon + family,guyane, mean)
aggD12 <- aggregate(D12  ~ taxon + family,guyane, mean)

guyane_m$D12 <- aggD12$D12

plot(RB~D12, guyane_m)
m <- lm(RB~ D12, guyane_m)
abline(m)
guyane_m$shrinkage_res <- residuals(m)


boxplot(shrinkage_res ~ family, guyane_m, xaxt="n", xlab="") 
axis(side = 1, at=1:nlevels(guyane_m$family),labels = FALSE)
text(x = 1:nlevels(guyane_m$family),
     y = par("usr")[3] - 0.45,
     labels = levels(guyane_m$family),
     xpd = NA,
      adj=c(1,0.5),
     ## Rotate the labels by 35 degrees.
     srt = 35,
     cex = 0.6)

zanne <- read.tree("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity\\Vascular_Plants_rooted.dated.txt")

tip <- gsub(" ", "_", guyane_m$taxon)
 
zanne$tip.label[match(tip, zanne$tip.label)]

in_zanne <- na.omit(match(tip, zanne$tip.label))
not_in_zanne <- which(is.na(match(tip, zanne$tip.label)))


tip[not_in_zanne]


match.phylo.data(zanne, tip)



tree <- keep.tip(zanne, in_zanne)


plot(tree, cex=0.5,no.margin=T)

t1<-add.species.to.genus(tree,"Tetragastris_sp")


plot(t1, cex=0.5,no.margin=T)
