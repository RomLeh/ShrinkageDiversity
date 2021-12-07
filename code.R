library(BIOMASS)
library(brranching)
library(ape)
library(phylobase)
library(picante)
library(phytools)
setwd("//NASLSB/Base de données/BD Cirad")

# data <- read.csv("Extraction_base_phy_meca-Essais meca.csv", header=T, sep=";", dec=".")
# guyane <- data[data$LIB_PAYS == "Guyane Francaise",]
# guyane$taxon <- paste(guyane$genus, guyane$sp, sep=" ")
# taxo_cor <- correctTaxo(guyane$taxon, useCache=F)

# guyane$genus <- taxo_cor$genusCorrected 
# guyane$sp <- taxo_cor$speciesCorrected 
# 
# family_result <- getTaxonomy(guyane$genus)
# guyane$family <- family_result$family
# 
# write.csv(guyane, "Essais meca_guyane.csv")


guyane <- read.csv("Essais meca_guyane.csv", header=T, sep=";", dec=".")
guyane$family <- as.factor(guyane$family)
guyane$sp[guyane$genus == "Tabebuia" & guyane$sp == "serratifolia"] <- "serratifolius"
guyane$genus[guyane$genus == "Tabebuia" & guyane$sp == "serratifolius"] <- "Handroanthus"
guyane$genus[guyane$genus == "Balizia" & guyane$sp == "pedicellaris"] <- "Albizia"
guyane$sp[guyane$genus == "Cedrelinga" & guyane$sp == "catenaeformis"] <- "cateniformis"
guyane$genus[guyane$genus == "Sclerolobium" & guyane$sp == "melinonii"] <- "Tachigali"
guyane$genus[guyane$genus == "Bombacopsis" & guyane$sp == "nervosa"] <- "Pachira"





guyane$taxon  <- paste(guyane$genus, guyane$sp, sep=" ")
guyane <- guyane[guyane$sp !='sp',]


guyane_m <- aggregate(RB  ~ taxon + family,guyane, mean)
aggD12 <- aggregate(D12  ~ taxon + family,guyane, mean)

guyane_m$D12 <- aggD12$D12



plot(RB ~ D12, guyane_m)
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


tip <- gsub(" ", "_", guyane_m$taxon)



zanne <- read.tree("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity\\Vascular_Plants_rooted.dated.txt")
in_zanne <- na.omit(match(tip, zanne$tip.label))
not_in_zanne <- which(is.na(match(tip, zanne$tip.label)))
treeZ <- keep.tip(zanne, in_zanne)
treeZinit <- force.ultrametric(treeZ, method=c("extend"))
plot(treeZinit, cex=0.5, no.margin=T,direction="up")
t1z <- treeZinit
for(i in 1:length(tip[not_in_zanne])){t1z<-add.species.to.genus(t1z, tip[not_in_zanne][i])}
treeZ <- t1z
plot(treeZ, cex=0.5, no.margin=T,direction="up")


sp_not_in_tree <- tip[which(is.na(match(tip, treeZ$tip.label)))]


janssens <- read.nexus("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity\\Janssens.tre")
janssens$tip.label[janssens$tip.label=="OUT_Pinus"] <- 'Pinus_caribaea'
janssens$tip.label[janssens$tip.label=="Xylopia_villosa"] <- 'Xylopia_grandiflora'
janssens$tip.label[janssens$tip.label=="Xylopia_aethiopica "] <- 'Xylopia_sericea'
janssens$tip.label[janssens$tip.label=="Ilex_crenata"] <- 'Ilex_casiquiarensis'
janssens$tip.label[janssens$tip.label=="Tabebuia_rosea"] <- 'Tabebuia_fluviatilis'
janssens$tip.label[janssens$tip.label=="Trattinnickia_demerarae"] <- 'Trattinnickia_rhoifolia'
janssens$tip.label[janssens$tip.label=="Tovomita_longifolia"] <- 'Tovomita_carinata'
janssens$tip.label[janssens$tip.label=="Andira_inermis"] <- 'Andira_coriacea'
janssens$tip.label[janssens$tip.label=="Copaifera_officinalis"] <- 'Copaifera_guianensis'
janssens$tip.label[janssens$tip.label=="Dimorphandra_conjugata"] <- 'Dimorphandra_polyandra'
janssens$tip.label[janssens$tip.label=="Peltogyne_floribunda"] <- 'Peltogyne_venosa'
janssens$tip.label[janssens$tip.label=="Vouacapoua_macropetala"] <- "Vouacapoua_americana"
janssens$tip.label[janssens$tip.label=="Aniba_williamsii"] <- "Aniba_melinonii"
janssens$tip.label[janssens$tip.label=="Albizia_coriaria"] <- "Albizia_pedicellaris"
janssens$tip.label[janssens$tip.label=="Eriotheca_roseorum"] <- "Eriotheca_crassa"
janssens$tip.label[janssens$tip.label=="Eriotheca_discolor"] <- "Eriotheca_globosa"
janssens$tip.label[janssens$tip.label=="Euplassa_duquei"] <- "Euplassa_pinnata"
janssens$tip.label[janssens$tip.label=="Zanthoxylum_ailanthoides"] <- "Zanthoxylum_rhoifolium"
janssens$tip.label[janssens$tip.label=="Simaba_trichilioides"] <- "Simaba_multiflora"
janssens$tip.label[janssens$tip.label=="Ampelocera_hottlei"] <- "Ampelocera_edentula"
janssens$tip.label[janssens$tip.label=="Cecropia_pachystachya"] <- "Cecropia_sciadophylla"
janssens$tip.label[janssens$tip.label=="Pachira_brevipes"] <- "Pachira_nervosa"



in_janssens <- na.omit(match(tip, janssens$tip.label))
not_in_janssens <- which(is.na(match(tip, janssens$tip.label)))
treeJ <- keep.tip(janssens, in_janssens)
treeJinit <- force.ultrametric(treeJ, method=c("extend"))
plot(treeJinit, cex=0.5, no.margin=T,direction="up")

t1j <- treeJinit
for(i in 1:length(tip[not_in_janssens])){
   t1j<-add.species.to.genus(t1j, tip[not_in_janssens][i])
}

treeJ <- t1j

plot(treeJ, cex=0.5, no.margin=T,direction="up")

identify(treeJ)


sp_not_in_treeJ <- tip[which(is.na(match(tip, treeJ$tip.label)))]
