library(BIOMASS)
library(brranching)
library(ape)
library(phylobase)
library(picante)
library(phytools)
library(Taxonstand)
library(ggtree)
library(treeio)
library(ggplot2)
library(tidyverse)

#setwd("//NASLSB/Base de données/BD Cirad")

setwd("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity")
data <- read.csv("Extraction_base_phy_meca-Essais meca.csv", header=T, sep=";", dec=".")
guyane <- data[data$LIB_PAYS == "Guyane Francaise",]
guyane$taxon <- paste(guyane$genus, guyane$sp, sep=" ")

### Correction de la bota par the PlantList
# TPL_cor <- TPL(unique(guyane$taxon))
# write.csv(TPL_cor,"TPL_cor.csv")
 
### Modification des noms 
TPL_cor           <- read.csv("TPL_cor.csv", header=T, sep=";", dec=".")
TPL_cor$taxon_cor <- with(TPL_cor, paste(New.Genus, New.Species, sep=' '))

# Taxon
newTaxon <- rep(NA, dim(guyane)[1])
for(i in 1:dim(guyane)[1]){
newTaxon[i]  <-  TPL_cor$taxon_cor[TPL_cor$Taxon == paste(guyane$taxon[i])] }
names(guyane)[33] <- "taxon_old"
guyane$taxon <- newTaxon

# genre espéce
for(i in 1:dim(guyane)[1]){
   guyane$genus[i]  <-  strsplit(guyane$taxon," ")[[i]][1] 
   guyane$sp[i]     <-  strsplit(guyane$taxon," ")[[i]][2]}
# suppression des sp
guyane <- guyane[guyane$sp !='sp',]

# ajout famille
guyane$family <- rep(NA, dim(guyane)[1])
for(i in 1:dim(guyane)[1]){
      guyane$family[i] <- TPL_cor$Family[TPL_cor$taxon_cor == paste(guyane$taxon[i])]}

guyane$family <- as.factor(guyane$family)
### Correction tachigali
guyane$genus[guyane$genus == "Sclerolobium" & guyane$sp == "melinonii"] <- "Tachigali"
guyane$taxon[guyane$taxon == "Sclerolobium melinonii"] <- "Tachigali melinonii"



### Retrait de cocos et recordo

guyane <- guyane[guyane$genus != "Cocos",]
guyane <- guyane[guyane$genus != "Recordoxylon",]


#### Création data moyénnés
guyane_m <- aggregate(RB  ~ taxon + family,guyane, mean)
aggD12 <- aggregate(D12  ~ taxon + family,guyane, mean)
guyane_m$D12 <- aggD12$D12

### Extaction residus RB vs D12
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



# zanne <- read.tree("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity\\Vascular_Plants_rooted.dated.txt")
# in_zanne <- na.omit(match(tip, zanne$tip.label))
# not_in_zanne <- which(is.na(match(tip, zanne$tip.label)))
# treeZ <- keep.tip(zanne, in_zanne)
# treeZinit <- force.ultrametric(treeZ, method=c("extend"))
# plot(treeZinit, cex=0.5, no.margin=T,direction="up")
# t1z <- treeZinit
# for(i in 1:length(tip[not_in_zanne])){t1z<-add.species.to.genus(t1z, tip[not_in_zanne][i])}
# treeZ <- t1z
# plot(treeZ, cex=0.5, no.margin=T,direction="up")
# sp_not_in_tree <- tip[which(is.na(match(tip, treeZ$tip.label)))]

#### #### #### #### #### #### #### #### #### #### #### #### #### 
#### Creation phylogeny basé sur phylo de Janssens et al.   #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### 
janssens <- read.nexus("E:\\Travaux\\Projets de recherche\\non-financés\\ShrinkageDiversity\\Janssens.tre")
janssens$tip.label[janssens$tip.label=="OUT_Pinus"] <- 'Pinus_caribaea'
janssens$tip.label[janssens$tip.label=="Xylopia_villosa"] <- 'Xylopia_aromatica'
janssens$tip.label[janssens$tip.label=="Xylopia_aethiopica "] <- 'Xylopia_sericea'
janssens$tip.label[janssens$tip.label=="Ilex_crenata"] <- 'Ilex_casiquiarensis'
janssens$tip.label[janssens$tip.label=="Tabebuia_rosea"] <- 'Tabebuia_fluviatilis'
janssens$tip.label[janssens$tip.label=="Trattinnickia_demerarae"] <- 'Trattinnickia_rhoifolia'
janssens$tip.label[janssens$tip.label=="Tovomita_longifolia"] <- 'Tovomita_carinata'
janssens$tip.label[janssens$tip.label=="Andira_inermis"] <- 'Andira_coriacea'
janssens$tip.label[janssens$tip.label=="Copaifera_officinalis"] <- 'Copaifera_guyanensis'
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

### Selection des taxon présents absent dans la phylogenie
in_janssens <- na.omit(match(tip, janssens$tip.label))
not_in_janssens <- which(is.na(match(tip, janssens$tip.label)))
## Elagage de l'arbre
treeJ <- keep.tip(janssens, in_janssens)
treeJinit <- force.ultrametric(treeJ, method=c("extend"))
#### Ajout des taxon absent de la phylogenie : créer des groupes paraphylétique OCOTEA
t1j <- treeJinit
for(i in 1:length(tip[not_in_janssens])){
   t1j<-add.species.to.genus(t1j, tip[not_in_janssens][i])}
treeJ <- t1j

plot(treeJ, cex=0.5, no.margin=T,direction="up")
sp_not_in_treeJ <- tip[which(is.na(match(tip, treeJ$tip.label)))]


rownames(guyane_m) <- guyane_m$taxon
treeJ$tip.label <- gsub("_", " ", treeJ$tip.label)

guyane_m <- guyane_m[treeJ$tip.label,]
guyane_m <- guyane_m[guyane_m$taxon != "Recordoxylon_speciosum",]
names(guyane_m)[1] <- 'label'

tree <- full_join(treeJ,guyane_m, by="label")

p <- ggtree(tree, layout="circular", aes(color=shrinkage_res),size=1) +
         scale_color_gradient2(low='gold', mid='gray', high='forestgreen') +
         geom_tiplab(size=2.5) + xlim(NA, 500)
p



p1 <- p %<+% guyane_m + geom_tiplab(size=2)
