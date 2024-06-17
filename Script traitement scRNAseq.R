###########################################################
# Script pour analyser les données de Single cell RNA Seq #
# Par Ethan Sourasinh stagiaire en M1 de bioinformatique  #
# Fait en 2024                                            #
###########################################################
sessionInfo()

setwd("/Users/osourasinh/Documents/Stage M1/Données SGLRNA")

setwd("C:/Users/osourasinh/WorkingDirectory/Stage M1/Données SGLRNA")

#############################
# Chargement des librairies #
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scales)
library(data.table)
library(clustree)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::version()
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# 
# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
# devtools::install_github("aertslab/SCENIC") 
# packageVersion("SCENIC")
# 
# # Github:
# devtools::install_github("aertslab/AUCell")
# devtools::install_github("aertslab/RcisTarget")
# devtools::install_github("aertslab/GENIE3")

# Bioconductor
# install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/PACKAGENAME.tar.gz", repos=NULL)
# 
# install.packages("remotes")
# remotes::install_github("aertslab/SCopeLoomR")
# remotes::install_github("mojaveazure/seurat-disk")
# remotes::install_github("jmvera255/SCENIC")
library(SCENIC)
library(SeuratDisk)
library(SCopeLoomR)
library(SingleCellExperiment)
library(devtools)
library(RcisTarget)

# Use devtools to install hdf5r and loomR from GitHub
# devtools::install_github(repo = "hhoeflin/hdf5r")
# devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
# 
# devtools::install_github(repo = "satijalab/seurat", ref = "loom")
# 
# devtools::install_github("hhoeflin/hdf5r")

library(Seurat)
library(loomR)
library(hdf5r)

# install.packages("doSNOW")
# 
# install.packages("doParallel") 
# 
# install.packages("doMPI")

library(doSNOW)
library(doParallel)
library(doMPI)
library(AUCell)
library(shiny)
library(KernSmooth)
library(RColorBrewer)
library(arrow)

# devtools::install_github("ncborcherding/scRepertoire")

library(scRepertoire)
# BiocManager::install("slingshot")
library(slingshot)
library(cowplot)
# BiocManager::install("BiocParallel")
# BiocParallel::register(BiocParallel::SerialParam())
# BiocManager::install("tradeSeq")
library(tradeSeq)
# .rs.restartR()
####################
# Données communes #

# # Chargement des listes de gènes signatures des différentes phases de réplication 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Fonction pour modifier les gènes en majuscules en minuscules (Homme -> Souris)
capwords <- function(s, strict = TRUE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

# Modification des listes de gènes signatures

s.genes <-capwords(s.genes)

g2m.genes <-capwords(g2m.genes)

Cytotoxsignature <- c("KLRA10","NCR1","KLRA9","KLRA8","KLRB1C","GZMA","KLRI2","KLRC1","KLRC2","KLRA3","PLEK","GZMB","FCER1G","KLRC3","CCL5","CCR5","IL18RAP","KLRB1F","KLRB1A","KLRK1","ITGAX","KLRB1B","KLRE1","ITGAM","KLRG1","CD244","IL12RB2","ANXA1","FGL2","CCR2","ITGA2","AOAH","OSBPL3","CCL3","TYROBP","SULF2","SAMD3","ADAMTS14","KLRA6","CAR5B","TBX21","SPRY2","ANXA2","SERPINB9","KLRA1","CAR2","PION","SLAMF7","ARSB","GEM","GM11435","LYN","STYK1","IL18R1","BHLHE40","SLC25A24","KHDC1A","FASL","DENND4A","SH2D1B1","CLNK","PRDM1","CTLA2A","SERPINB6B","ID2","SYKB","IL10RA","ITGB1","ATP2B4","SERPINB9B","EOMES","S100A6","PIK3AP1","SLC9A7","ZEB2","ARL4D","RNF216","NRARP","PLCG2","LAT2","RBPMS","KLRI1","KLRA5","ENTPD1","IL2RB","PRF1","CASP1","TCRG-C","KIT","KLRD1","I830127L07RIK","S1PR5","AW112010")

Cytotoxsignature <- capwords(Cytotoxsignature)

# Création de la liste des gènes de cycle cellualaire
CyclingGenes <- c(s.genes, g2m.genes)

# Création de la signature cycling g2m
Cyclingg2msignature <- c("Top2a","Tubb5","Nusap1","H2afx","Cdk1","Hist1h2ap","Spc25","Prc1","Mki67","Hist1h2ae","Ube2s","Kpna2","Cdca3","Spc24","Aurkb","Smc2","Hist1h1b","Ccna2","Esco2","Incenp","Fbxo5","Kif22","H2afv","Kif11","Pbk","Mxd3","Tk1","Tacc3","Asf1b","Racgap1","Tpx2","Tubb4b","Cdca2","Ube2c","Ndc80","Kif15","Ska1","Hmmr","Cit","Mis18bp1","Smc4","Tyms","Cenpf","Cks2","Plk1","Knl1","Dek","Birc5","Ncapg","Cdkn2c","Sgo1","Rrm2","Hist1h2ab","Ezh2","Cdca8","Cenpe","Bub1b","Arl6ip1","Calm2","Kif20b","Hist1h2ak","Shcbp1","Hjurp","Cdca5","Aurka","Hist1h4i","Jpt1","Hist1h2ag","Ccdc34","Tmpo","Hist1h2bj","Kif20a","Lockd","Rrm1","Bub3","Cenpm","Kif23","Rad51","Tuba1c","Rad51ap1","Ncapd2","Hist1h4d","Hist1h1e","Hist1h3c","Ccnb1","Mad2l1","Pmf1","Hnrnpa3","Clspn","Lbr","Rad21","Cmc2","Cenph","Rangap1","Nrm","Prdx4","Gm42031","Hnrnpab","Smc1a","Sumo2","Cenpw","Nsd2","Cbx5","Lsm2","Usp1","Rbbp4","Syce2","Gmnn","Hmgb3","Alyref","Nucks1","Rfc5","Fen1","Dnajc9","Rfc4","Smarca5","Prim1","Lig1","Tceal9","Lsm3","Srsf3","Cdc20","Knstrn","Haus4","Tcf19","Ccne1","Crip1","Gm10282","Dtymk","Hnrnpa2b1","Banf1","Rbbp7","Uhrf1","Dnmt1","Snrpd1","Hpf1","Hnrnpd","Snrpe","Hprt","Orc6","Cycs","Tfdp1","Nasp","Fam111a","Snrpb","Ssrp1","Igfbp4","Tpi1","Srsf7","Nrgn")

# Création de la signature cycling s
cyclingssignature <- c("Ptma","Ran","Ppia","Ybx1","Eif5a","H2afz","Hmgn2","Slc25a5","Npm1","Gapdh","Cdca7","Ldha","Ranbp1","Set","Hmgb2","Dut","Hmgb1","Hells","Srsf2","Ncl","Txn1","Mif","Erh","Hspd1","Ung","Phgdh","Prdx1","Rps27l","Nop58","Mcm6","Tipin","Mcm5","Selenoh","Hint1","Dctpp1","Anp32b","Hsp90b1","Pa2g4","Ccnb2","Nap1l1","Lmnb1","Mcm2","Cks1b","Dbi","Anp32e","Stmn1","Mcm7","Lgals1","Hsp90aa1","Fabp5","Pcna","Mcm3","Cenpa","Pgk1","Slbp","Pclaf","Tuba1b")

# # Création des la signature Mait0
# MAIT0signature <-c ("Lef1","Itm2a","Nrgn","Bcl2","Slamf6","Izumo1r","Gimap6","Malat1","Id3","Cd28","Ccr9","Cd2","Slc29a1","Gimap5","Satb1","Ccr7","Il21r","Tox","Marcksl1","Adk","Cd27","Cldn10","Cd8a","Tuba1a","Cd4","Actn1","Cd81","Atp1b1","Myb","Cd5","Kcnn4","Ubac2","Sox4","Rasgrp1","Ptprc","Hivep3","Ssbp2","Tspan13","Ldhb","Cd8b1","Tnfrsf9","Ly6d”,”Cytip","Mrps6","Gsn","Myl10","Rhoh","Themis","Egr2","Nab2","Dntt","Cd69","Cd24a","Gpr83","Asap1","Tmem108","Frat2","Rnf167","Chst2","Tespa1","Spint2","H2-Q2","Fam169b","Egr1","Lztfl1","Tnfsf8","Zfp330","Tubb2b","Cux1","St6gal1","Prkca","Slc25a17","Arpp21","Zcchc12","Mpp4","Ubash3a","Bach2","Traf4","Pip4k2a","Basp1","Patz1","Aqp11","Nr4a1","Lad1","Zfp281","Siae","Mir142hg","Tsc22d1","Vamp1")

# Création de la signature Mait1
MAIT1signature <- c("Klrd1","Cd8a", "Klrc1", "Samd3", "Plac8", "Il12rb2", "Ifng", "Klrb1c", "Serpina3g", "Fasl", "Klrk1", "Tbx21", "Slamf7", "Zfp683", "Art2b", "Ly6c2", "Xcl1", "Ccl5", "Nkg7")

# Création de la signature Mait17
MAIT17signature <- c("Tmem176a","Tmem176b","Rorc","Slc27a6","Il17rb","Itgae","Il23r","Serpinb1a","Mmp25","1700040d17rik","Pxdc1","Il1r1","Lingo4","Ccr8","Ccr6","Cysltr1","Ltb4r1")

# création de la signature precursor
precursorsignature <- c("Dntt","Rag1","Rag2")

 

# # Création de la signature TCR
# TCRsignature<- c("Gng4","Egr3","Egr2","C2cd4a","Ebi3","Nr4a1","Tnfrsf9","Mir146a","Tnfrsf4","Espl1")

# Création de la signature des cellules en cycle cellulaire
# CyclingGenessignature <- CyclingGenes

# Création de la signature TCRmini
# TCRminisignature <- c("Egr3","Egr2","Nr4a1")

# # Création de la signature Egr2reg
# Egr2regsignature <- c("Bcl2","Birc5","Cd27","Cd69","Ckap2","Cytip","Dusp2","Dusp6","Egr1","Egr3","Id3","Naa25","Nab2","Nr4a1","Nrgn","Pabpn1","Rbms2","Rps2","Slamf6","Stk17b","Tox","Ttc3")

# Création de la signature Mait0 
MAIT0signature <- c("Itm2a","Ccr9","Satb1","Cd8b1","Bcl2","Lef1","Cd8a","Cd27","Cd28","Tox","Gimap6","Rhoh","Tuba1a","Cd2","Sox4","Cd247","Trbc1","Nrgn","Gm43352","Slc29a1","Bcl2","Cd28","Tox","Cldn10","Slamf6","Hivep3","Id3","Myb","Ubac2","Izumo1r","Cytip","Marcksl1","Tubb5","Adk","Egr2","Ldhb","Rasgrp1","H3f3a","Mrps6","B630019A10Rik")


# Interm_Asignature <- c("Igfbp4”,”Ms4a4b”,”Izumo1r”,”Ifi27l2a”,”Klf2”,”Ccr7")
# 
# Interm_Bsignature <- c("Plac8")

# Interm_A_Longsignature <- c("Igfbp4","Ms4a4b","Izumo1r","Ifi27l2a","Klf2","Ccr7","Lef1","Cd2","Sell","Id3","Ms4a6b","Gimap4","Cd4","S1pr1","Gm8369","Tsc22d3","Cytip","Actn1","Gimap6","Dgka","Adk","Kcnn4","Drosha","Klhdc2","Tspan32","Gm12840","Smc4","Ass1","Rasgrp2","Ddit4","Inpp4b","H2-K1","Bcl2","Gimap3","Cd28","Atp1b1","Cd5","Cd247","B630019A10Rik","Gimap1","Dusp2","Itm2a","Slamf6","Foxp1","Satb1","Slfn1")
# 
# Interm_B_Longsignature <- c("Plac8","Ifi27l2a","Izumo1r","Tesc","Ass1","Slamf6","Igfbp4","Ms4a6b","Ms4a4b")
# 
# Immature_A_Longsignature <- c("Itm2a","Ccr9","Satb1","Cd8b1","Bcl2","Lef1","Cd8a","Cd27","Cd28","Tox","Gimap6","Rhoh","Tuba1a","Cd2","Sox4","Cd247","Trbc1","Ldhb","Malat1","Ets1","Cd81","Tcf7","Foxp1","Ssbp2","Lztfl1","Trbc2","Gimap4","Tespa1","Tuba4a","Gsn","Cd24a","Ptprc","Fyb","Trac","Gimap9","Tspan13","Actn1","Acp5","Tmem108","Bcl11b","B630019A10Rik","Tubb5","Themis","Prkca","Kcnn4","Cd4","Il21r","Grap","Lck","Ccr7","Tia1","Sh2d1a","Slc38a2","Dntt","Hnrnpf","Eif4g2","Dnaja1","Stmn1","Pnrc1","Pip4k2a","Tmsb10","Gimap3","Fam169b","Traf4","Fus","H3f3a","Pik3ip1","Tra2b","Supt4a","Ccnd3","Gm43352","Zfp281","Hnrnpk","Arhgap45","Tdrp","Ftl1","Dusp10","H3f3b","Nono","Ephb6","Dynll1","Ube2d3")
# 
# Immature_B_Longsignature <- c("Itm2a","Lef1","Nrgn","Gm43352","Slc29a1","Bcl2","Cd28","Tox","Cldn10","Slamf6","Cd27","Hivep3","Id3","Cd2","Myb","Ubac2","Sox4","Izumo1r","Cytip","Marcksl1","Tubb5","Adk","Sep-07","Egr2","Ldhb","Rasgrp1","H3f3a","Mrps6","B630019A10Rik","Cd53","Pdcd1","Egr1","Ikzf1","Sh2d1a","Tespa1","Frat2","Malat1","Nab2","St6gal1","Rgs10","Rhoh","Rap1a","Slc25a17","H2-Q2","Cd81","Arhgap15","Cd5","Ptprc","Plgrkt","Tuba1a","Cd84","Stmn1","Zfp330","Basp1","Atp1b1","Gimap5","Tmsb10","Hmgn1","Gimap6","Ier3ip1","Prkch","Gpr83","Cd69","Il21r","Gsn","Jarid2","Slc9a9","Rnf167","Gnai3","Pfn1","Ddx19a","Asap1","Cux1","Spint2","Cap2","Supt4a","Mpp4","Ccnd3","Tm6sf1","Med14","Als2","Trac","Rbm4b","Nono","Ttc3","Map1lc3b","Rbmx","Bin2","Tia1","Tubb2b","Arid5b","Kcnn4","Zcchc12","Eif4g2","Siae","Hnrnpdl","Lcp1","Gm11998","Trbc2","Cdh24","Chst2","Gadd45a","Tra2b","Arpc2","Tnfsf8","Srsf3","Gm14718","Elf1","Hnrnpf","Ephx1","Jpt1","Fus")



# Pensez à faire un setwd("DOSSIER") où DOSSIER est le dossier contenant le dossier "filtered_feature_bc_matrix"
# Chargement des données issus du séquençage scRNAseq via 10X :
Data10X <- Read10X("./filtered_feature_bc_matrix/")

# Création de l'objet Seurat principal assay par défaut RNA
mouse = CreateSeuratObject(counts = Data10X$`Gene Expression`)
# Création d'assay supplémentaires pour séparer les données
mouse [["HTO"]]=CreateAssayObject(counts = Data10X$`Antibody Capture`)
mouse [["Protein"]] = CreateAssayObject(counts = Data10X$Protein)
# Normalisation des données selon les Hashtags (HTO)
mouse <- NormalizeData(mouse,assay = "HTO", normalization.method = "CLR")
# Assigne les cellules à leur hashtag, threshold = positive.quantile
mouse <- HTODemux(mouse, assay = "HTO", positive.quantile = 0.999)

# mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")
# VlnPlot(mouse, "percent.mt")
# 
# table(mouse$HTO_classification.global)
# 
# Idents(mouse) <- "HTO_maxID"
# RidgePlot(mouse, assay = "HTO", features = rownames(mouse[["HTO"]])[1:2], ncol = 2)
# 
# FeatureScatter(mouse, feature1 = "hto-Hashtag-1", feature2 = "hto-Hashtag-2")
# Enlève les billes qui contiennent deux cellules.
mouse.subset <- subset(mouse, idents = "Doublet", invert = TRUE)
# Enlève les billes qui ne contient pas de cellule.
mouse.subset <- subset(mouse.subset, idents = "Negative", invert = TRUE)
# table(mouse.subset$HTO_classification.global)

######################
# Séparation Hashtag #

# Séparation de l'objet principal selon les hashtags
mouseHashtag = SplitObject(mouse.subset)

# Isolation du Hashtag 1
mouseHashtag.subset1 <- mouseHashtag$`Hashtag-1`
# Normalisation des données
mouseHashtag.subset1 <- NormalizeData(mouseHashtag.subset1)
# Mise à l'échelle des données
mouseHashtag.subset1 <- ScaleData(mouseHashtag.subset1,verbose = FALSE)

# Isolation du Hashtag 2
mouseHashtag.subset2 <- mouseHashtag$`Hashtag-2`
# Normalisation des données
mouseHashtag.subset2 <- NormalizeData(mouseHashtag.subset2)
# Mise à l'échelle des données
mouseHashtag.subset2 <- ScaleData(mouseHashtag.subset2,verbose = FALSE)

# Isolation du Hashtag 3
mouseHashtag.subset3 <- mouseHashtag$`Hashtag-3`
# Normalisation des données
mouseHashtag.subset3 <- NormalizeData(mouseHashtag.subset3)
# Mise à l'échelle des données
mouseHashtag.subset3 <- ScaleData(mouseHashtag.subset3,verbose = FALSE)

# Isolation du Hashtag 4
mouseHashtag.subset4 <- mouseHashtag$`Hashtag-4`
# Normalisation des données
mouseHashtag.subset4 <- NormalizeData(mouseHashtag.subset4)
# Mise à l'échelle des données
mouseHashtag.subset4 <- ScaleData(mouseHashtag.subset4,verbose = FALSE)

# Isolation du Hashtag 5
mouseHashtag.subset5 <- mouseHashtag$`Hashtag-5`
# Normalisation des données
mouseHashtag.subset5 <- NormalizeData(mouseHashtag.subset5)
# Mise à l'échelle des données
mouseHashtag.subset5 <- ScaleData(mouseHashtag.subset5,verbose = FALSE)

# Isolation du Hashtag 6
mouseHashtag.subset6 <- mouseHashtag$`Hashtag-6`
# Normalisation des données
mouseHashtag.subset6 <- NormalizeData(mouseHashtag.subset6)
# Mise à l'échelle des données
mouseHashtag.subset6 <- ScaleData(mouseHashtag.subset6,verbose = FALSE)


#########################
# Hashtag12 downsampled #
# Fusion des données des Hashtags 1 et 2
mouseHashtag12down = merge(mouseHashtag.subset1,mouseHashtag.subset2)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag12down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag12down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag12down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag12down, pattern = "^Ig")

# VlnPlot(mouseHashtag12down, features = "percent.mt")
 # Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag12down <- subset(mouseHashtag12down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag12down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag12down <- subset(mouseHashtag12down, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag12down<-JoinLayers(mouseHashtag12down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag12down <- CellCycleScoring(mouseHashtag12down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag12down <- AddModuleScore(object = mouseHashtag12down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag12down)<- mouseHashtag12down$old.ident

# Normalisation des données
mouseHashtag12down <- NormalizeData(mouseHashtag12down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag12down$old.ident)

# On choisit au hasard des cellules du Hashtag 1 pour rammener au nombre de cellule du Hashtag 2; ici 284.
mouseHashtag12down<- subset(mouseHashtag12down, downsample= min(table(mouseHashtag12down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag12down <- FindVariableFeatures(mouseHashtag12down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag12down <- ScaleData(mouseHashtag12down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag12down <- RunPCA(mouseHashtag12down, features = rownames(mouseHashtag12down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag12down <- RunTSNE(mouseHashtag12down, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag12down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag12down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag12down <- FindNeighbors(mouseHashtag12down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag12down <- FindClusters(mouseHashtag12down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag12down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag12down <- FindClusters(mouseHashtag12down, resolution = 0.9)

# Umap selon les 15 dimensions de la PCA
mouseHashtag12down <- RunUMAP(mouseHashtag12down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag12down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot12downgenes <- FeaturePlot(mouseHashtag12down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot12downgenes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot hashtag 1 et 2 downsampled avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot12downgenes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot12downsignatures <- FeaturePlot(mouseHashtag12down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot12downsignatures + DimPlot(mouseHashtag12down, reduction = "umap")

# VlnPlot(mouseHashtag12down,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag12down,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag12down,features = c("Xist"), split.by = "old.ident")

# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot hashtag 1 et 2 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12signatures + DimPlot(mouseHashtag12down, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids12down <- c("Mait1", "Mait17a","Mait0", "Mait17b","CyclingS","CyclingG2M","Precursor")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids12down) <- levels(mouseHashtag12down)
mouseHashtag12down <- RenameIdents(mouseHashtag12down, new.cluster.ids12down)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap12downlabel <-DimPlot(mouseHashtag12down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap12downlabel

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/UMAP hashtag 1 et 2 downsampled.png", res = 150, width = 600)
# Umap12downlabel
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot12downsignaturessplithashtag <- FeaturePlot(mouseHashtag12down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "precursorscore1","cyclingsscore1","cyclingg2mscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot12downsignaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12downsignaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag12down <- DimPlot(mouseHashtag12down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag12down
Umapgrouphashtag12down <- DimPlot(mouseHashtag12down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag12down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/UMAP hashtag 1 et Umap hashtag 2downsampled.png", res = 150, width = 1000)
# Umapsplithashtag12down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/UMAP hashtag 1 et hashtag 2 downsampled.png", res = 150, width = 1000)
# Umapgrouphashtag12down
# dev.off()

VlnPlot(mouseHashtag12down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot Egr2 par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/FeaturePlot Irf5 par cluster et par hashtag 12 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12down,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag12down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag12down, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag12down,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag12down,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag12down,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag12down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster12down  <-table(mouseHashtag12down$seurat_clusters,mouseHashtag12down$old.ident)

# Affichage de la table
TableCluster12down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster12down)/2)){
  colonne1 <- c(colonne1 , TableCluster12down[i,1]/(TableCluster12down[i,1]+ TableCluster12down[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster12down[i,2]/(TableCluster12down[i,1]+ TableCluster12down[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster12down [i,1] / min(table(mouseHashtag12down$old.ident)))*100)
  colonne4 <- c(colonne4 , (TableCluster12down [i,2] / min(table(mouseHashtag12down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO12downCluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO12downCluster$`Hashtag 1nb_cell`<-c(TableCluster12down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO12downCluster$`Hashtag 2nb_cell`<-c(TableCluster12down[,2])


# Changement des noms de colonnes
colnames(PercentHTO12downCluster)<- c("Pourcentage intercluster Hashtag 1"," Pourcentage intercluster Hashtag 2","Pourcentage Hashtag1","Pourcentage Hashtag 2","Hashtag 1nb_cell","Hashtag 2nb_cell")

# Changement des noms de lignes
rownames(PercentHTO12downCluster)<- new.cluster.ids12down

PercentHTO12downCluster
# Enregistrement du tableau de données
# write.csv(PercentHTO12downCluster, file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/PercentHTO12downCluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag12down.markers <- FindAllMarkers(mouseHashtag12down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag12down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_12down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap12down <-DoHeatmap(mouseHashtag12down, features = top10_12down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap12down
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap hashtag 1 et 2 downsampled par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap12down
# dev.off()

# Séparation selon les hashtags
mouseHashtag12splitdown <- SplitObject(mouseHashtag12down,split.by = "old.ident")
# Hashtag 1
mouseHashtag12split1down <- mouseHashtag12splitdown$`Hashtag-1`
mouseHashtag12split1down <- JoinLayers(mouseHashtag12split1down)

# Hashtag 2
mouseHashtag12split2down <- mouseHashtag12splitdown$`Hashtag-2`
mouseHashtag12split2down <- JoinLayers(mouseHashtag12split2down)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag12split1down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12split1down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag12split2down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12split2down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag1down.markers <- FindAllMarkers(mouseHashtag12split1down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag1down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_1down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap1down <-DoHeatmap(mouseHashtag12split1down, features = top10_1down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap1down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag2down.markers <- FindAllMarkers(mouseHashtag12split2down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap2down <-DoHeatmap(mouseHashtag12split2down, features = top10_2down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2down

# Affichage des deux heatmaps ensembles
heatmap1down + heatmap2down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap hashtag 1 et hashtag 2 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap1down + heatmap2down
# dev.off()

# Création de nouveaux identifiants pour la suite
mouseHashtag12down$ClusterGroup <- paste(mouseHashtag12down@active.ident, mouseHashtag12down$old.ident, sep = "_")
Idents(mouseHashtag12down) <- "ClusterGroup"

# Recherche des markers du cluster Mait17 du Hashtag 1 différent du Hashtag 2
MarkerMait17adown <- FindMarkers(mouseHashtag12down, ident.1 = "Mait17a_Hashtag-1", ident.2 = "Mait17a_Hashtag-2")
head(MarkerMait17adown, n = 10)

MarkerMait17bdown <- FindMarkers(mouseHashtag12down, ident.1 = "Mait17b_Hashtag-1", ident.2 = "Mait17b_Hashtag-2")
head(MarkerMait17bdown, n = 10)

# Recherche des markers du cluster Mait1 du Hashtag 1 différent du Hashtag 2
MarkerMait1down <- FindMarkers(mouseHashtag12down, ident.1 = "Mait1_Hashtag-1", ident.2 = "Mait1_Hashtag-2")
head(MarkerMait1down, n = 10)

# Recherche des markers du cluster Immature du Hashtag 1 différent du Hashtag 2
MarkerMait0down <- FindMarkers(mouseHashtag12down, ident.1 = "Mait0_Hashtag-1", ident.2 = "Mait0_Hashtag-2")
head(MarkerMait0down, n = 10)

# Recherche des markers du cluster Cycling du Hashtag 1 différent du Hashtag 2
MarkerCyclingSdown <- FindMarkers(mouseHashtag12down, ident.1 = "CyclingS_Hashtag-1", ident.2 = "CyclingS_Hashtag-2")

MarkerCyclingG2Mdown <- FindMarkers(mouseHashtag12down, ident.1 = "CyclingG2M_Hashtag-1", ident.2 = "CyclingG2M_Hashtag-2")

MarkerPrecursordown <- FindMarkers(mouseHashtag12down, ident.1 = "Precursor_Hashtag-1", idents.2 = "Precursor_Hashtag-2")

MarkerMait1down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1down

MarkerMait17adown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17adown

MarkerMait17bdown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17bdown

MarkerCyclingSdown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingSdown

MarkerMait0down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait0down

MarkerCyclingG2Mdown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingG2Mdown

MarkerPrecursordown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Precursordown

Idents(mouseHashtag12down)<- mouseHashtag12down$old.ident
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait1down)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait17adown)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait17bdown)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait0down)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_CyclingSdown)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_CyclingG2Mdown)))
DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Precursordown)))

FeaturePlot(mouseHashtag12down,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag12down, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerCyclingG2Mdown,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerCyclingG2MdownsampledHashtag1.csv")
# 
# write.csv(MarkerCyclingSdown,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerCyclingSdownsampledHashtag1.csv")
# 
# write.csv(MarkerMait0down,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerMait0downsampledHashtag1.csv")
# 
# write.csv(MarkerMait1down,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerMAit1downsampledHashtag1.csv")
# 
# write.csv(MarkerMait17adown,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerMait17adownsampledHashtag1.csv")
# 
# write.csv(MarkerMait17bdown,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerMait17bdownsampledHashtag1.csv")
# 
# write.csv(MarkerPrecursordown,file = "./Dossier Graphiques et Tableaux/Hashtag12 downsampled/GenesmarkerPrecursordownsampledHashtag1.csv")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10Mait1 hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait1down)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10Mait17a hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait17adown)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10Mait17b hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait17bdown)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10Mait0 hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Mait0down)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10CyclingSdown hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_CyclingSdown)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12 downsampled/heatmap top10Precursor hashtag 1 et hashtag 2 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12down, features = c(rownames(top10_Precursordown)))
# dev.off()

#############
# Hashtag12 #
# Fusion des données des Hashtags 1 et 2
mouseHashtag12 = merge(mouseHashtag.subset1,mouseHashtag.subset2)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag12[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag12, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag12[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag12, pattern = "^Ig")

# VlnPlot(mouseHashtag12, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag12 <- subset(mouseHashtag12, subset = percent.mt < 3)

# VlnPlot(mouseHashtag12, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag12 <- subset(mouseHashtag12, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag12<-JoinLayers(mouseHashtag12)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag12 <- CellCycleScoring(mouseHashtag12, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag12 <- AddModuleScore(object = mouseHashtag12, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag12)<- mouseHashtag12$old.ident

# Normalisation des données
mouseHashtag12 <- NormalizeData(mouseHashtag12)


# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag12 <- FindVariableFeatures(mouseHashtag12, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag12 <- ScaleData(mouseHashtag12,verbose = FALSE)

# PCA selon les gènes
mouseHashtag12 <- RunPCA(mouseHashtag12, features = rownames(mouseHashtag12), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag12 <- RunTSNE(mouseHashtag12, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag12)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag12)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag12 <- FindNeighbors(mouseHashtag12, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag12 <- FindClusters(mouseHashtag12, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag12) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag12 <- FindClusters(mouseHashtag12, resolution = 0.9)

# Umap selon les 15 dimensions de la PCA
mouseHashtag12 <- RunUMAP(mouseHashtag12, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag12, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot12genes <- FeaturePlot(mouseHashtag12, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot12genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot hashtag 1 et 2  avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot12genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot12signatures <- FeaturePlot(mouseHashtag12, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot12signatures + DimPlot(mouseHashtag12, reduction = "umap")

# VlnPlot(mouseHashtag12,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag12,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag12,features = c("Xist"), split.by = "old.ident")
# VlnPlot(mouseHashtag12,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 12 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot hashtag 1 et 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12signatures + DimPlot(mouseHashtag12, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids12 <- c("Mait1", "Mait17a","Mait0", "Mait17b","CyclingS","CyclingG2M","Precursor")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids12) <- levels(mouseHashtag12)
mouseHashtag12 <- RenameIdents(mouseHashtag12, new.cluster.ids12)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap12label <-DimPlot(mouseHashtag12, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap12label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/UMAP hashtag 1 et 2 .png", res = 150, width = 600)
# Umap12label
# dev.off()

DotPlot(mouseHashtag12,features = c("Dntt","Rorc","Tbx21","Zbtb16","Mki67"),col.min=0, split.by = "old.ident", cols = c("red","green","blue"))

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot12signaturessplithashtag <- FeaturePlot(mouseHashtag12, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot12signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag12 <- DimPlot(mouseHashtag12, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag12
Umapgrouphashtag12 <- DimPlot(mouseHashtag12, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag12

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/UMAP hashtag 1 et Umap hashtag 2.png", res = 150, width = 1000)
# Umapsplithashtag12
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/UMAP hashtag 1 et hashtag 2 .png", res = 150, width = 1000)
# Umapgrouphashtag12
# dev.off()

VlnPlot(mouseHashtag12,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot Zbtb16 par cluster et par hashtag 12.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot Egr2 par cluster et par hashtag 12.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot Bcl2a1b par cluster et par hashtag 12.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot Nr4a1 par cluster et par hashtag 12.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/FeaturePlot Irf5 par cluster et par hashtag 12.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag12, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag12, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag12,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag12,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag12,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag12,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster12  <-table(mouseHashtag12$seurat_clusters,mouseHashtag12$old.ident)

# Affichage de la table
TableCluster12

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster12)/2)){
  colonne1 <- c(colonne1 , TableCluster12[i,1]/(TableCluster12[i,1]+ TableCluster12[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster12[i,2]/(TableCluster12[i,1]+ TableCluster12[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster12 [i,1] / sum(TableCluster12[,1]))*100)
  colonne4 <- c(colonne4 , (TableCluster12 [i,2] / sum(TableCluster12[,2]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO12Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO12Cluster$`Hashtag 1nb_cell`<-c(TableCluster12[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO12Cluster$`Hashtag 2nb_cell`<-c(TableCluster12[,2])


# Changement des noms de colonnes
colnames(PercentHTO12Cluster)<- c("Pourcentage intercluster Hashtag 1"," Pourcentage intercluster Hashtag 2","Pourcentage Hashtag1","Pourcentage Hashtag 2","Hashtag 1nb_cell","Hashtag 2nb_cell")

# Changement des noms de lignes
rownames(PercentHTO12Cluster)<- new.cluster.ids12

PercentHTO12Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO12Cluster, file="./Dossier Graphiques et Tableaux/Hashtag12/PercentHTO12Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag12.markers <- FindAllMarkers(mouseHashtag12, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag12.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_12

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap12 <-DoHeatmap(mouseHashtag12, features = top10_12$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap12
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap hashtag 1 et 2  par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap12
# dev.off()

# Séparation selon les hashtags
mouseHashtag12split <- SplitObject(mouseHashtag12,split.by = "old.ident")
# Hashtag 1
mouseHashtag12split1 <- mouseHashtag12split$`Hashtag-1`
mouseHashtag12split1 <- JoinLayers(mouseHashtag12split1)

# Hashtag 2
mouseHashtag12split2 <- mouseHashtag12split$`Hashtag-2`
mouseHashtag12split2 <- JoinLayers(mouseHashtag12split2)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag12split1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12split1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag12split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag1.markers <- FindAllMarkers(mouseHashtag12split1, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_1

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap1 <-DoHeatmap(mouseHashtag12split1, features = top10_1$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap1


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag2.markers <- FindAllMarkers(mouseHashtag12split2, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap2 <-DoHeatmap(mouseHashtag12split2, features = top10_2$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2

# Affichage des deux heatmaps ensembles
heatmap1 + heatmap2

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap hashtag 1 et hashtag 2 par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap1 + heatmap2
# dev.off()



# Création de nouveaux identifiants pour la suite
mouseHashtag12$ClusterGroup <- paste(mouseHashtag12@active.ident, mouseHashtag12$old.ident, sep = "_")
Idents(mouseHashtag12) <- "ClusterGroup"

# Recherche des markers du cluster Mait17 du Hashtag 1 différent du Hashtag 2
MarkerMait17a <- FindMarkers(mouseHashtag12, ident.1 = "Mait17a_Hashtag-1", ident.2 = "Mait17a_Hashtag-2")
head(MarkerMait17a, n = 10)

MarkerMait17b <- FindMarkers(mouseHashtag12, ident.1 = "Mait17b_Hashtag-1", ident.2 = "Mait17b_Hashtag-2")
head(MarkerMait17b, n = 10)

# Recherche des markers du cluster Mait1 du Hashtag 1 différent du Hashtag 2
MarkerMait1 <- FindMarkers(mouseHashtag12, ident.1 = "Mait1_Hashtag-1", ident.2 = "Mait1_Hashtag-2")
head(MarkerMait1, n = 10)

# Recherche des markers du cluster Immature du Hashtag 1 différent du Hashtag 2
MarkerMait0 <- FindMarkers(mouseHashtag12, ident.1 = "Mait0_Hashtag-1", ident.2 = "Mait0_Hashtag-2")
head(MarkerMait0, n = 10)


# Recherche des markers du cluster Cycling du Hashtag 1 différent du Hashtag 2
MarkerCyclingS <- FindMarkers(mouseHashtag12, ident.1 = "CyclingS_Hashtag-1", ident.2 = "CyclingS_Hashtag-2")

MarkerCyclingG2M <- FindMarkers(mouseHashtag12, ident.1 = "CyclingG2M_Hashtag-1", ident.2 = "CyclingG2M_Hashtag-2")

MarkerPrecursor <- FindMarkers(mouseHashtag12, ident.1 = "Precursor_Hashtag-1", idents.2 = "Precursor_Hashtag-2")

MarkerMait1 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1

MarkerMait17a %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17a

MarkerMait17b %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17b

MarkerCyclingS %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingS

MarkerMait0 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait0

MarkerCyclingG2M %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingG2M

MarkerPrecursor %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Precursor

Idents(mouseHashtag12)<- mouseHashtag12$old.ident
DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait1)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait17a)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait17b)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_CyclingS)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait0)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_CyclingG2M)))
DoHeatmap(mouseHashtag12, features = c(rownames(top10_Precursor)))

DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait0))) + scale_fill_gradient(low = "blue", high = "red")
DotPlot(mouseHashtag12, features = c(rownames(top10_CyclingS)))

FeaturePlot(mouseHashtag12,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag12, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerCyclingG2M,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerCyclingG2MHashtag1.csv")
# 
# write.csv(MarkerCyclingS,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerCyclingSHashtag1.csv")
# 
# write.csv(MarkerMait0,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerMait0Hashtag1.csv")
# 
# write.csv(MarkerMait1,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerMAit1Hashtag1.csv")
# 
# write.csv(MarkerMait17a,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerMait17aHashtag1.csv")
# 
# write.csv(MarkerMait17b,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerMait17bHashtag1.csv")
# 
# write.csv(MarkerPrecursor,file = "./Dossier Graphiques et Tableaux/Hashtag12/GenesmarkerPrecursorHashtag1.csv")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterMAit0 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait0)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterMAit1 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait1)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterMAit17a par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait17a)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterMAit17b par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_Mait17b)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterCyclingS par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_CyclingS)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag12/heatmap gène clusterPrecursor par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12, features = c(rownames(top10_Precursor)))
# dev.off()


##################
# Hashtag12clean #

mouseHashtag12clean <- subset(mouseHashtag12, subset = seurat_clusters != 6)


# Normalisation des données
mouseHashtag12clean <- NormalizeData(mouseHashtag12clean)


# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag12clean <- FindVariableFeatures(mouseHashtag12clean, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag12clean <- ScaleData(mouseHashtag12clean,verbose = FALSE)

# PCA selon les gènes
mouseHashtag12clean <- RunPCA(mouseHashtag12clean, features = rownames(mouseHashtag12clean), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag12clean <- RunTSNE(mouseHashtag12clean, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag12clean)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag12clean)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag12clean <- FindNeighbors(mouseHashtag12clean, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag12clean <- FindClusters(mouseHashtag12clean, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag12clean) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag12clean <- FindClusters(mouseHashtag12clean, resolution = 0.6)

# Umap selon les 15 dimensions de la PCA
mouseHashtag12clean <- RunUMAP(mouseHashtag12clean, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag12clean, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot12genesclean <- FeaturePlot(mouseHashtag12clean, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot12genesclean
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot hashtag 1 et 2 clean avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot12genesclean
# dev.off()


FeaturePlot(mouseHashtag12clean, features = c("Zbtb16","Ccr9","Itm2a","Id2","Id3"))
# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot12signaturesclean <- FeaturePlot(mouseHashtag12clean, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1"), min.cutoff = 0.2)
Featureplot12signaturesclean + DimPlot(mouseHashtag12clean, reduction = "umap")

# VlnPlot(mouseHashtag12clean,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag12clean,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag12clean,features = c("Xist"), split.by = "old.ident")
# VlnPlot(mouseHashtag12clean,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 12 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot hashtag 1 et 2 clean avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12signaturesclean + DimPlot(mouseHashtag12clean, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids12clean <- c("Mait17a","Mait1","Mait0","Mait17b","CyclingS","CyclingG2M")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids12clean) <- levels(mouseHashtag12clean)
mouseHashtag12clean <- RenameIdents(mouseHashtag12clean, new.cluster.ids12clean)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap12labelclean <-DimPlot(mouseHashtag12clean, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap12labelclean

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/UMAP hashtag 1 et 2 clean.png", res = 150, width = 1000)
# Umap12labelclean
# dev.off()

# DotPlot(mouseHashtag12clean,features = c("Dntt","Rorc","Tbx21","Zbtb16","Mki67"),col.min=0, split.by = "old.ident", cols = c("red","green","blue"))

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot12signaturessplithashtagclean <- FeaturePlot(mouseHashtag12clean, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot12signaturessplithashtagclean

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot12signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
UmapsplitHashtag12clean <- DimPlot(mouseHashtag12clean, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
UmapsplitHashtag12clean
UmapgroupHashtag12clean <- DimPlot(mouseHashtag12clean, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
UmapgroupHashtag12clean

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/UMAP hashtag 1 et Umap hashtag 2clean.png", res = 150, width = 1000)
# UmapsplitHashtag12clean
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/UMAP hashtag 1 et hashtag 2 clean.png", res = 150, width = 1000)
# UmapgroupHashtag12clean
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot Zbtb16 par cluster et par hashtag 12clean.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot Egr2 par cluster et par hashtag 12clean.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot Bcl2a1b par cluster et par hashtag 12clean.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot Nr4a1 par cluster et par hashtag 12clean.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/FeaturePlot Irf5 par cluster et par hashtag 12clean.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag12clean,features = c("Irf5"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag12clean,features = c("Rorc"),split.by = "old.ident")


# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag12clean, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag12clean, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag12clean,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag12clean,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag12clean,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag12clean,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster12clean  <-table(mouseHashtag12clean$seurat_clusters,mouseHashtag12clean$old.ident)

# Affichage de la table
TableCluster12clean

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster12clean)/2)){
  colonne1 <- c(colonne1 , TableCluster12clean[i,1]/(TableCluster12clean[i,1]+ TableCluster12clean[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster12clean[i,2]/(TableCluster12clean[i,1]+ TableCluster12clean[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster12clean [i,1] / sum(TableCluster12clean[,1]))*100)
  colonne4 <- c(colonne4 , (TableCluster12clean [i,2] / sum(TableCluster12clean[,2]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO12Clusterclean <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO12Clusterclean$`Hashtag 1nb_cell`<-c(TableCluster12clean[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO12Clusterclean$`Hashtag 2nb_cell`<-c(TableCluster12clean[,2])


# Changement des noms de colonnes
colnames(PercentHTO12Clusterclean)<- c("Pourcentage intercluster Hashtag 1"," Pourcentage intercluster Hashtag 2","Pourcentage Hashtag1","Pourcentage Hashtag 2","Hashtag 1nb_cell","Hashtag 2nb_cell")

# Changement des noms de lignes
rownames(PercentHTO12Clusterclean)<- new.cluster.ids12clean

PercentHTO12Clusterclean
# Enregistrement du tableau de données
# write.csv(PercentHTO12Clusterclean, file="./Dossier Graphiques et Tableaux/Hashtag12clean/PercentHTO12Clusterclean.csv")

colonne3 = as.numeric(colonne3)
colonne4 = as.numeric(colonne4)

Percentintracluster = as.data.frame(c(colonne3,colonne4))
Percentintracluster$cluster = rownames(PercentHTO12Clusterclean)
Percentintracluster$hashtag = c("Hashtag 1","Hashtag 1","Hashtag 1","Hashtag 1","Hashtag 1","Hashtag 1", "Hashtag 2","Hashtag 2","Hashtag 2","Hashtag 2","Hashtag 2","Hashtag 2")
colnames(Percentintracluster) <- c("Pourcentage","cluster","hashtag")

ggplot(Percentintracluster, aes(x = cluster, y = Pourcentage))+
  geom_col(aes(fill = hashtag), position = position_dodge(0.8), width = 0.7)+
  scale_fill_discrete(name = "Condition", labels = c("Axénique","Monoxénique"))

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/barplot pourcentage par cluster et par hashtag.png", res = 150, height = 1800, width = 1500)
# ggplot(Percentintracluster, aes(x = cluster, y = Pourcentage))+
# geom_col(aes(fill = hashtag), position = position_dodge(0.8), width = 0.7)+
# scale_fill_discrete(name = "Condition", labels = c("Axénique","Monoxénique"))
# dev.off()




# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag12clean.markers <- FindAllMarkers(mouseHashtag12clean, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag12clean.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_12clean

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap12clean <-DoHeatmap(mouseHashtag12clean, features = top10_12clean$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap12clean
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap hashtag 1 et 2  clean par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap12clean
# dev.off()

# Séparation selon les hashtags
mouseHashtag12cleansplit <- SplitObject(mouseHashtag12clean,split.by = "old.ident")
# Hashtag 1
mouseHashtag12cleansplit1 <- mouseHashtag12cleansplit$`Hashtag-1`
mouseHashtag12cleansplit1 <- JoinLayers(mouseHashtag12cleansplit1)

# Hashtag 2
mouseHashtag12cleansplit2 <- mouseHashtag12cleansplit$`Hashtag-2`
mouseHashtag12cleansplit2 <- JoinLayers(mouseHashtag12cleansplit2)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag12cleansplit1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12cleansplit1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag12cleansplit2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag12cleansplit2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag1clean.markers <- FindAllMarkers(mouseHashtag12cleansplit1, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag1clean.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_1clean

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap1clean <-DoHeatmap(mouseHashtag12cleansplit1, features = top10_1clean$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap1clean


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag2clean.markers <- FindAllMarkers(mouseHashtag12cleansplit2, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2clean.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2clean

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap2clean <-DoHeatmap(mouseHashtag12cleansplit2, features = top10_2clean$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2clean

# Affichage des deux heatmaps ensembles
heatmap1clean + heatmap2clean

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap hashtag 1 et hashtag 2 clean par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap1clean + heatmap2clean
# dev.off()


# Création de nouveaux identifiants pour la suite
mouseHashtag12clean$ClusterGroup <- paste(mouseHashtag12clean@active.ident, mouseHashtag12clean$old.ident, sep = "_")
Idents(mouseHashtag12clean) <- "ClusterGroup"


# Recherche des markers du cluster Mait17a du Hashtag 1 différent du Hashtag 2
MarkerMait17aclean <- FindMarkers(mouseHashtag12clean, ident.1 = "Mait17a_Hashtag-1", ident.2 = "Mait17a_Hashtag-2")

# Recherche des markers du cluster Mait17b du Hashtag 1 différent du Hashtag 2
MarkerMait17bclean <- FindMarkers(mouseHashtag12clean, ident.1 = "Mait17b_Hashtag-1", ident.2 = "Mait17b_Hashtag-2")


# Recherche des markers du cluster Mait1 du Hashtag 1 différent du Hashtag 2
MarkerMait1clean <- FindMarkers(mouseHashtag12clean, ident.1 = "Mait1_Hashtag-1", ident.2 = "Mait1_Hashtag-2")

# Recherche des markers du cluster Mait0 du Hashtag 1 différent du Hashtag 2
MarkerMait0clean <- FindMarkers(mouseHashtag12clean, ident.1 = "Mait0_Hashtag-1", ident.2 = "Mait0_Hashtag-2")

# Recherche des markers du cluster CyclingS du Hashtag 1 différent du Hashtag 2
MarkerCyclingSclean <- FindMarkers(mouseHashtag12clean, ident.1 = "CyclingS_Hashtag-1", ident.2 = "CyclingS_Hashtag-2")

# Recherche des markers du cluster CyclingG2M du Hashtag 1 différent du Hashtag 2
MarkerCyclingG2Mclean <- FindMarkers(mouseHashtag12clean, ident.1 = "CyclingG2M_Hashtag-1", ident.2 = "CyclingG2M_Hashtag-2")


MarkerMait1clean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1clean

MarkerMait17aclean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17aclean

MarkerMait17bclean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17bclean

MarkerCyclingSclean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingSclean

MarkerMait0clean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait0clean

MarkerCyclingG2Mclean %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_CyclingG2Mclean


Idents(mouseHashtag12clean)<- mouseHashtag12clean$old.ident
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait1clean)))
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait17aclean)))
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait17bclean)))
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_CyclingSclean)))
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait0clean)))
DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_CyclingG2Mclean)))

# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait0))) + scale_fill_gradient(low = "blue", high = "red")
# DotPlot(mouseHashtag12clean, features = c(rownames(top10_CyclingS)))
# 
# FeaturePlot(mouseHashtag12clean,features = "Xist", split.by = "ident",min.cutoff = 0)
# FeaturePlot(mouseHashtag12clean, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerCyclingG2Mclean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerCyclingG2MHashtag1clean.csv")
# 
# write.csv(MarkerCyclingSclean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerCyclingSHashtag1clean.csv")
# 
# write.csv(MarkerMait0clean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerMait0Hashtag1clean.csv")
# 
# write.csv(MarkerMait1clean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerMAit1Hashtag1clean.csv")
# 
# write.csv(MarkerMait17aclean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerMait17aHashtag1clean.csv")
# 
# write.csv(MarkerMait17bclean,file = "./Dossier Graphiques et Tableaux/Hashtag12clean/GenesmarkerMait17bHashtag1clean.csv")



# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterMAit0clean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait0clean)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterMAit1clean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait1clean)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterMAit17aclean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait17aclean)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterMAit17bclean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_Mait17bclean)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterCyclingSclean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_CyclingSclean)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag12clean/heatmap gène clusterCyclingG2Mclean par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag12clean, features = c(rownames(top10_CyclingG2Mclean)))
# dev.off()


############################
# Hashtag12clean slingshot #

dimred12 <- mouseHashtag12clean@reductions$umap@cell.embeddings
clustering12 <- mouseHashtag12clean$RNA_snn_res.0.8
counts12 <- as.matrix(mouseHashtag12clean@assays$RNA$counts, mouseHashtag12clean@assays$RNA@var.features)


set.seed(1)
lineages12 <- getLineages(data = dimred12,
                        clusterLabels = clustering12,
                        #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                        start.clus = "2") #define where to start the trajectories

lineages12<-SlingshotDataSet(lineages12)

lineages12

par(mfrow=c(1,2))
plot(dimred12[,1:2], col = clustering12,  cex=.5,pch = 16)
for(i in levels(clustering12)){ 
  text( mean(dimred12[clustering12==i,1]),
        mean(dimred12[clustering12==i,2]), labels = i,font = 2) }
plot(dimred12, col = clustering12,  pch = 16)
lines(lineages12, lwd = 3, col = 'black')

curves12<- getCurves(lineages12, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves12<-SlingshotDataSet(curves12)
curves12

plot(dimred12, col = clustering12, asp = 1, pch = 16)
lines(curves12, lwd = 3, col = "black")



filt_counts12 <- counts12[rowSums(counts12 > 5) > ncol(counts12)/100, ]
dim(filt_counts12)

sce12 <- fitGAM(counts = as.matrix(filt_counts12), sds = curves12)

plotGeneCount(curves12, filt_counts12, clusters = clustering12, models = sce12)

# png("./Dossier Graphiques et Tableaux/Hashtag12clean/Slingshot courbes Hashtag12clean.png", res = 95, height = 600, width = 600)
# plotGeneCount(curves12, filt_counts12, clusters = clustering12, models = sce12)
# dev.off()

plot_differential_expression12 <- function(feature_id12) {
  feature_id12 <- pseudotime_association12 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id12)
  cowplot::plot_grid(plotGeneCount(curves12, filt_counts12, gene = "Xcl1", clusters = clustering12, models = sce12) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce12, as.matrix(counts12), gene = "Xcl1"))
}


pseudotime_association12 <- associationTest(sce12)
pseudotime_association12$fdr <- p.adjust(pseudotime_association12$pvalue, method = "fdr")
pseudotime_association12 <- pseudotime_association12[order(pseudotime_association12$pvalue), ]
pseudotime_association12$feature_id12 <- rownames(pseudotime_association12)

feature_id12 <- pseudotime_association12 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id12)
plot_differential_expression12(feature_id12)



Pseudotimes12 <- as.data.frame(slingPseudotime(curves12))
Pseudotimes12$CellName <- rownames(Pseudotimes12)
HashtagPseudotime12 <- as.data.frame(table(rownames(mouseHashtag12clean@assays$RNA@cells), mouseHashtag12clean$old.ident))
HashtagPseudotime12$Hashtag <- ifelse (HashtagPseudotime12$Var2 == "Hashtag-1" & HashtagPseudotime12$Freq == "1" |HashtagPseudotime12$Var2 == "Hashtag-2" & HashtagPseudotime12$Freq == "0", "Hashtag-1", "Hashtag-2") 
HashtagPseudotime12$Cluster <- mouseHashtag12clean@active.ident
HashtagPseudotime12 = HashtagPseudotime12 [,c(-2,-3)]
names(HashtagPseudotime12)[1] <- "CellName"


HashtagPseudotime12 = merge(Pseudotimes12,HashtagPseudotime12, by = "CellName")

ggplot(HashtagPseudotime12,aes(Lineage1,Cluster, colour = Cluster)) +
       geom_violin()


HashtagPseudotime12split <- split(HashtagPseudotime12,HashtagPseudotime12$Hashtag)
HashtagPseudotime12split1 <- HashtagPseudotime12split$`Hashtag-1`
HashtagPseudotime12split2 <- HashtagPseudotime12split$`Hashtag-2`





# Hashtag1clean

dimred1 <- mouseHashtag12cleansplit1@reductions$umap@cell.embeddings
clustering1 <- mouseHashtag12cleansplit1$RNA_snn_res.0.9
counts1 <- as.matrix(mouseHashtag12cleansplit1@assays$RNA$counts, mouseHashtag12cleansplit1@assays$RNA@var.features)


set.seed(1)
lineages1 <- getLineages(data = dimred1,
                        clusterLabels = clustering1,
                        start.clus = "2") #define where to start the trajectories

lineages1<-SlingshotDataSet(lineages1)

lineages1

par(mfrow=c(1,2))
plot(dimred1[,1:2], col = clustering1,  cex=.5,pch = 16)
for(i in levels(clustering1)){ 
  text( mean(dimred1[clustering1==i,1]),
        mean(dimred1[clustering1==i,2]), labels = i,font = 2) }
plot(dimred1, col = clustering1,  pch = 16)
lines(lineages1, lwd = 3, col = 'black')

curves1 <- getCurves(lineages1, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves1<-SlingshotDataSet(curves1)
curves1

plot(dimred1, col = clustering1, asp = 1, pch = 16)
lines(curves1, lwd = 3, col = "black")

filt_counts1 <- counts1[rowSums(counts1 > 5) > ncol(counts1)/100, ]
dim(filt_counts1)

sce1 <- fitGAM(counts = as.matrix(filt_counts1), sds = curves1)

plotGeneCount(curves1, filt_counts1, clusters = clustering1, models = sce1)

plot_differential_expression1 <- function(feature_id1) {
  feature_id1 <- pseudotime_association1 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id1)
  cowplot::plot_grid(plotGeneCount(curves1, filt_counts1, gene = "Zbtb16", clusters = clustering1, models = sce1) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce1, as.matrix(counts1), gene = "Zbtb16"))
}


pseudotime_association1 <- associationTest(sce1)
pseudotime_association1$fdr <- p.adjust(pseudotime_association1$pvalue, method = "fdr")
pseudotime_association1 <- pseudotime_association1[order(pseudotime_association1$pvalue), ]
pseudotime_association1$feature_id1 <- rownames(pseudotime_association1)

feature_id1 <- pseudotime_association1 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id1)
plot_differential_expression1(feature_id1)

# Hashtag2clean

dimred2 <- mouseHashtag12cleansplit2@reductions$umap@cell.embeddings
clustering2 <- mouseHashtag12cleansplit2$RNA_snn_res.0.9
counts2 <- as.matrix(mouseHashtag12cleansplit2@assays$RNA$counts, mouseHashtag12cleansplit2@assays$RNA@var.features)

set.seed(1)
lineages2 <- getLineages(data = dimred2, clusterLabels = clustering2,start.clus = "2")

lineages2<-SlingshotDataSet(lineages2)

lineages2

par(mfrow = c(1, 2))
plot(dimred2[, 1:2], col = clustering2, cex = 0.5, pch = 16)
for (i in levels(clustering2)) {
  text(mean(dimred2[clustering2 == i, 1]), mean(dimred2[clustering2 == i, 2]), labels = i, font = 2)
}
plot(dimred2[, 1:2], col = clustering2, cex = 0.5, pch = 16)
lines(lineages2, lwd = 3, col = "black")

set.seed(1)
lineages2 <- getLineages(data = dimred2,
                        clusterLabels = clustering2,
                        start.clus = "2") #define where to start the trajectories

lineages2<-SlingshotDataSet(lineages2)

lineages2

par(mfrow=c(1,2))
plot(dimred2[,1:2], col = clustering2,  cex=.5,pch = 16)
for(i in levels(clustering2)){ 
  text( mean(dimred2[clustering2==i,1]),
        mean(dimred2[clustering2==i,2]), labels = i,font = 2) }
plot(dimred2, col = clustering2,  pch = 16)
lines(lineages2, lwd = 3, col = 'black')

curves2 <- getCurves(lineages2, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves2 <-SlingshotDataSet(curves2)
curves2

plot(dimred2, col = clustering2, asp = 1, pch = 16)
lines(curves2, lwd = 3, col = "black")


filt_counts2 <- counts2[rowSums(counts2 > 5) > ncol(counts2)/100, ]
dim(filt_counts2)

sce2 <- fitGAM(counts = as.matrix(filt_counts2), sds = curves2)

plotGeneCount(curves2, filt_counts2, clusters = clustering2, models = sce2)

plot_differential_expression2 <- function(feature_id2) {
  feature_id2 <- pseudotime_association2 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id2)
  cowplot::plot_grid(plotGeneCount(curves2, filt_counts2, gene = "Zbtb16", clusters = clustering2, models = sce2) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce2, as.matrix(counts2), gene = "Zbtb16"))
}

pseudotime_association2 <- associationTest(sce2)
pseudotime_association2$fdr <- p.adjust(pseudotime_association2$pvalue, method = "fdr")
pseudotime_association2 <- pseudotime_association2[order(pseudotime_association2$pvalue), ]
pseudotime_association2$feature_id2 <- rownames(pseudotime_association2)

feature_id2 <- pseudotime_association2 %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id2)
plot_differential_expression2(feature_id2)





feature_id

pseudotime_start_end_association <- startVsEndTest(sce, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)

feature_id <- pseudotime_start_end_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)

plot_differential_expression(feature_id)

different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)

branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)

feature_id <- branch_point_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)


mouse.traj.mait17 <- Add.pseudotime.cluster(SeuratObject = mouseHashtag12cleansplit1, slingshot.dataset = curves1, reduction = "umap", start.cluster = 2, seurat.cluster = "seurat_clusters", select.trajectory = 1, bin = 5, stretch = 0, extend =  "n")

mouse.traj.mait1 <- Add.pseudotime.cluster(SeuratObject = mouseHashtag12cleansplit1, slingshot.dataset = curves1, reduction = "umap", start.cluster = 2, seurat.cluster = "seurat_clusters", select.trajectory = 2, bin = 5, stretch = 0, extend =  "n")



list.dataset.separate <- list(mouse.traj.mait17,mouse.traj.mait1)
names(list.dataset.separate) <- c("Hashtag1 MAIT17", "Hashtag1 MAIT1")

megaplot.gene.pseudotime(list.dataset = list.dataset.separate , feature = "Zbtb16", 
                                color = c("#332288", "#882255", "#88CCEE", "#CC6677", "#44AA99", "#DDCC77", "#999933", "#117733"),
                                lwd = 3, bin = 10, inset = c(-0.2,10))


#########################
# Hashtag23 downsampled #
# Fusion des données des Hashtags 2 et 3
mouseHashtag23down = merge(mouseHashtag.subset2,mouseHashtag.subset3)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag23down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag23down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag23down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag23down, pattern = "^Ig")

# VlnPlot(mouseHashtag23down, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag23down <- subset(mouseHashtag23down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag23down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag23down <- subset(mouseHashtag23down, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag23down<-JoinLayers(mouseHashtag23down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag23down <- CellCycleScoring(mouseHashtag23down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag23down <- AddModuleScore(object = mouseHashtag23down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag23down)<- mouseHashtag23down$old.ident

# Normalisation des données
mouseHashtag23down <- NormalizeData(mouseHashtag23down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag23down$old.ident)

# On choisit au hasard des cellules du Hashtag 1 pour rammener au nombre de cellule du Hashtag 2; ici 284.
mouseHashtag23down<- subset(mouseHashtag23down, downsample= min(table(mouseHashtag23down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag23down <- FindVariableFeatures(mouseHashtag23down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag23down <- ScaleData(mouseHashtag23down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag23down <- RunPCA(mouseHashtag23down, features = rownames(mouseHashtag23down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag23down <- RunTSNE(mouseHashtag23down, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag23down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag23down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag23down <- FindNeighbors(mouseHashtag23down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag23down <- FindClusters(mouseHashtag23down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag23down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag23down <- FindClusters(mouseHashtag23down, resolution = 0.6)

# Umap selon les 15 dimensions de la PCA
mouseHashtag23down <- RunUMAP(mouseHashtag23down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag23down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot23downgenes <- FeaturePlot(mouseHashtag23down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot23downgenes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot hashtag 2 et 3 downsampled avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot23downgenes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot23downsignatures <- FeaturePlot(mouseHashtag23down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot23downsignatures + DimPlot(mouseHashtag23down, reduction = "umap")

# VlnPlot(mouseHashtag23down,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag23down,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag23down,features = c("Xist"), split.by = "old.ident")


# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Egr2 par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot hashtag 2 et 3 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot23downsignatures + DimPlot(mouseHashtag23down, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids23down <- c("Mait1", "Mait17","Mait0","Cycling")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids23down) <- levels(mouseHashtag23down)
mouseHashtag23down <- RenameIdents(mouseHashtag23down, new.cluster.ids23down)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap23downlabel <-DimPlot(mouseHashtag23down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap23downlabel

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/UMAP hashtag 2 et 3 downsampled.png", res = 150, width = 600)
# Umap23downlabel
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot23downsignaturessplithashtag <- FeaturePlot(mouseHashtag23down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1","precursorscore1","cyclingsscore1","cyclingg2mscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot23downsignaturessplithashtag

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot hashtag 2 et FeaturePlot hashtag 3 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot23downsignaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag23down <- DimPlot(mouseHashtag23down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag23down
Umapgrouphashtag23down <- DimPlot(mouseHashtag23down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag23down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/UMAP hashtag 2 et Umap hashtag 3d ownsampled.png", res = 150, width = 1000)
# Umapsplithashtag23down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/UMAP hashtag 2 et hashtag 3 downsampled.png", res = 150, width = 1000)
# Umapgrouphashtag23down
# dev.off()

VlnPlot(mouseHashtag23down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Egr2 par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/FeaturePlot Irf5 par cluster et par hashtag 23 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23down,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag23down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag23down, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag23down,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag23down,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag23down,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag23down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster23down  <-table(mouseHashtag23down$seurat_clusters,mouseHashtag23down$old.ident)

# Affichage de la table
TableCluster23down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster23down)/2)){
  colonne1 <- c(colonne1 , TableCluster23down[i,1]/(TableCluster23down[i,1]+ TableCluster23down[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster23down[i,2]/(TableCluster23down[i,1]+ TableCluster23down[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster23down [i,1] / min(table(mouseHashtag23down$old.ident)))*100)
  colonne4 <- c(colonne4 , (TableCluster23down [i,2] / min(table(mouseHashtag23down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO23downCluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO23downCluster$`Hashtag 2nb_cell`<-c(TableCluster23down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO23downCluster$`Hashtag 3nb_cell`<-c(TableCluster23down[,2])


# Changement des noms de colonnes
colnames(PercentHTO23downCluster)<- c("Pourcentage intercluster Hashtag 2"," Pourcentage intercluster Hashtag 3","Pourcentage Hashtag2","Pourcentage Hashtag 3","Hashtag 2nb_cell","Hashtag 3nb_cell")

# Changement des noms de lignes
rownames(PercentHTO23downCluster)<- new.cluster.ids23down

PercentHTO23downCluster
# Enregistrement du tableau de données
# write.csv(PercentHTO23downCluster, file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/PercentHTO23downCluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag23down.markers <- FindAllMarkers(mouseHashtag23down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag23down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_23down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap23down <-DoHeatmap(mouseHashtag23down, features = top10_23down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap23down
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/heatmap hashtag 2 et 3 downsampled par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap23down
# dev.off()

# Séparation selon les hashtags
mouseHashtag23splitdown <- SplitObject(mouseHashtag23down,split.by = "old.ident")
# Hashtag 2
mouseHashtag23split2down <- mouseHashtag23splitdown$`Hashtag-2`
mouseHashtag23split2down <- JoinLayers(mouseHashtag23split2down)

# Hashtag 3
mouseHashtag23split3down <- mouseHashtag23splitdown$`Hashtag-3`
mouseHashtag23split3down <- JoinLayers(mouseHashtag23split3down)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag23split2down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag23split2down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag23split3down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag23split3down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag2down.markers <- FindAllMarkers(mouseHashtag23split2down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap2down <-DoHeatmap(mouseHashtag23split2down, features = top10_2down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag3down.markers <- FindAllMarkers(mouseHashtag23split3down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag3down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_3down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap3down <-DoHeatmap(mouseHashtag23split3down, features = top10_3down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap3down

# Affichage des deux heatmaps ensembles
heatmap2down + heatmap3down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/heatmap hashtag 2 et hashtag 3 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap2down + heatmap3down
# dev.off()


# Création de nouveaux identifiants pour la suite
mouseHashtag23down$ClusterGroup <- paste(mouseHashtag23down@active.ident, mouseHashtag23down$old.ident, sep = "_")
Idents(mouseHashtag23down) <- "ClusterGroup"

# Recherche des markers du cluster Mait17 du Hashtag 2 différents du Hashtag 3
MarkerMait17 <- FindMarkers(mouseHashtag23down, ident.1 = "Mait17_Hashtag-2", ident.2 = "Mait17_Hashtag-3")
head(MarkerMait17, n = 10)

# Recherche des markers du cluster Mait1 du Hashtag 2 différents du Hashtag 3
MarkerMait1 <- FindMarkers(mouseHashtag23down, ident.1 = "Mait1_Hashtag-2", ident.2 = "Mait1_Hashtag-3")
head(MarkerMait1, n = 10)

# Recherche des markers du cluster Immature du Hashtag 1 différent du Hashtag 2
MarkerMait0 <- FindMarkers(mouseHashtag23down, ident.1 = "Mait0_Hashtag-2", ident.2 = "Mait0_Hashtag-3")
head(MarkerMait0, n = 10)


# Recherche des markers du cluster Cycling du Hashtag 1 différent du Hashtag 2
MarkerCycling <- FindMarkers(mouseHashtag23down, ident.1 = "Cycling_Hashtag-2", ident.2 = "Cycling_Hashtag-3")

MarkerMait1 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1

MarkerMait17 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17

MarkerMait0 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait0

MarkerCycling %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cycling

Idents(mouseHashtag23down)<- mouseHashtag23down$old.ident
DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Mait1)))
DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Mait17)))
DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Cycling)))
DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Mait0)))


# write.csv(MarkerCycling,file = "./Dossier Graphiques et Tableaux/Hashtag23 downsampled/GenesmarkerCyclingHashtag2.csv")
# 
# write.csv(MarkerMait0,file = "./Dossier Graphiques et Tableaux/Hashtag23 downsampled/GenesmarkerMait0Hashtag2.csv")
# 
# write.csv(MarkerMait1,file = "./Dossier Graphiques et Tableaux/Hashtag23 downsampled/GenesmarkerMait1Hashtag2.csv")
# 
# write.csv(MarkerMait17,file = "./Dossier Graphiques et Tableaux/Hashtag23 downsampled/GenesmarkerMait17Hashtag2.csv")

 
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/heatmap gène clusterMAit1 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Mait1)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/heatmap gène clusterMAit17 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Mait17)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag23 downsampled/heatmap gène clusterCycling par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23down, features = c(rownames(top10_Cycling)))
# dev.off()


#############
# Hashtag23 #
# Fusion des données des Hashtags 2 et 3
mouseHashtag23 = merge(mouseHashtag.subset2,mouseHashtag.subset3)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag23[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag23, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag23[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag23, pattern = "^Ig")

# VlnPlot(mouseHashtag23, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag23 <- subset(mouseHashtag23, subset = percent.mt < 3)

# VlnPlot(mouseHashtag23, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag23 <- subset(mouseHashtag23, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag23<-JoinLayers(mouseHashtag23)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag23 <- CellCycleScoring(mouseHashtag23, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag23 <- AddModuleScore(object = mouseHashtag23, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag23)<- mouseHashtag23$old.ident

# Normalisation des données
mouseHashtag23 <- NormalizeData(mouseHashtag23)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag23$old.ident)

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag23 <- FindVariableFeatures(mouseHashtag23, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag23 <- ScaleData(mouseHashtag23,verbose = FALSE)

# PCA selon les gènes
mouseHashtag23 <- RunPCA(mouseHashtag23, features = rownames(mouseHashtag23), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag23 <- RunTSNE(mouseHashtag23, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag23)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag23)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag23 <- FindNeighbors(mouseHashtag23, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag23 <- FindClusters(mouseHashtag23, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag23) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag23 <- FindClusters(mouseHashtag23, resolution = 0.6)

# Umap selon les 15 dimensions de la PCA
mouseHashtag23 <- RunUMAP(mouseHashtag23, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag23, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot23genes <- FeaturePlot(mouseHashtag23, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot23genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot hashtag 2 et 3  avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot23genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot23signatures <- FeaturePlot(mouseHashtag23, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot23signatures + DimPlot(mouseHashtag23, reduction = "umap")

# VlnPlot(mouseHashtag23,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag23,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag23,features = c("Xist"), split.by = "old.ident")



# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Egr2 par cluster et par hashtag 23 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot hashtag 2 et 3  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot23signatures
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids23 <- c("Mait1", "Mait17","Mait0", "Cycling")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids23) <- levels(mouseHashtag23)
mouseHashtag23 <- RenameIdents(mouseHashtag23, new.cluster.ids23)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap23label <-DimPlot(mouseHashtag23, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap23label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/UMAP hashtag 2 et 3 .png", res = 150, width = 600)
# Umap23label
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot23signaturessplithashtag <- FeaturePlot(mouseHashtag23, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot23signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot23signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag23 <- DimPlot(mouseHashtag23, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag23
Umapgrouphashtag23 <- DimPlot(mouseHashtag23, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag23

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/UMAP hashtag 2 et Umap hashtag 3.png", res = 150, width = 1000)
# Umapsplithashtag23
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/UMAP hashtag 1 et hashtag 2 .png", res = 150, width = 1000)
# Umapgrouphashtag23
# dev.off()

VlnPlot(mouseHashtag23,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Zbtb16 par cluster et par hashtag 23.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Egr2 par cluster et par hashtag 23.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Bcl2a1b par cluster et par hashtag 23.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Nr4a1 par cluster et par hashtag 23.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag23,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag23/FeaturePlot Irf5 par cluster et par hashtag 23.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag23,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag23, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag23, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag23,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag23,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag23,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag23,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster23  <-table(mouseHashtag23$seurat_clusters,mouseHashtag23$old.ident)

# Affichage de la table
TableCluster23

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster23)/2)){
  colonne1 <- c(colonne1 , TableCluster23[i,1]/(TableCluster23[i,1]+ TableCluster23[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster23[i,2]/(TableCluster23[i,1]+ TableCluster23[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster23 [i,1] / sum(TableCluster23[,1]))*100)
  colonne4 <- c(colonne4 , (TableCluster23 [i,2] / sum(TableCluster23[,2]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO23Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO23Cluster$`Hashtag 1nb_cell`<-c(TableCluster23[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO23Cluster$`Hashtag 2nb_cell`<-c(TableCluster23[,2])


# Changement des noms de colonnes
colnames(PercentHTO23Cluster)<- c("Pourcentage intercluster Hashtag 2"," Pourcentage intercluster Hashtag 3","Pourcentage Hashtag2","Pourcentage Hashtag 3","Hashtag 2nb_cell","Hashtag 3nb_cell")

# Changement des noms de lignes
rownames(PercentHTO23Cluster)<- new.cluster.ids23

PercentHTO23Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO23Cluster, file="./Dossier Graphiques et Tableaux/Hashtag23/PercentHTO23Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag23.markers <- FindAllMarkers(mouseHashtag23, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag23.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_23

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap23 <-DoHeatmap(mouseHashtag23, features = top10_23$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap23
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/heatmap hashtag 1 et 2  par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap23
# dev.off()

# Séparation selon les hashtags
mouseHashtag23split <- SplitObject(mouseHashtag23,split.by = "old.ident")
# Hashtag 2
mouseHashtag23split2 <- mouseHashtag23split$`Hashtag-2`
mouseHashtag23split2 <- JoinLayers(mouseHashtag23split2)

# Hashtag 3
mouseHashtag23split3 <- mouseHashtag23split$`Hashtag-3`
mouseHashtag23split3 <- JoinLayers(mouseHashtag23split3)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag23split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag23split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag23split3, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag23split3, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag2.markers <- FindAllMarkers(mouseHashtag23split2, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap2 <-DoHeatmap(mouseHashtag23split2, features = top10_2$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag3.markers <- FindAllMarkers(mouseHashtag23split3, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_3

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap3 <-DoHeatmap(mouseHashtag23split3, features = top10_3$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap3

# Affichage des deux heatmaps ensembles
heatmap2 + heatmap3

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/heatmap hashtag 2 et hashtag 3 par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap2 + heatmap3
# dev.off()


# Création de nouveaux identifiants pour la suite
mouseHashtag23$ClusterGroup <- paste(mouseHashtag23@active.ident, mouseHashtag23$old.ident, sep = "_")
Idents(mouseHashtag23) <- "ClusterGroup"

# Recherche des markers du cluster Mait17a du Hashtag 2 différents du Hashtag 3
MarkerMait17 <- FindMarkers(mouseHashtag23, ident.1 = "Mait17_Hashtag-2", ident.2 = "Mait17_Hashtag-3")
head(MarkerMait17, n = 10)

# Recherche des markers du cluster Mait1 du Hashtag 2 différents du Hashtag 3
MarkerMait1 <- FindMarkers(mouseHashtag23, ident.1 = "Mait1_Hashtag-2", ident.2 = "Mait1_Hashtag-3")
head(MarkerMait1, n = 10)

# Recherche des markers du cluster Immature du Hashtag 1 différent du Hashtag 2
MarkerMait0 <- FindMarkers(mouseHashtag23, ident.1 = "Mait0_Hashtag-2", ident.2 = "Mait0_Hashtag-3")
head(MarkerMait0, n = 10)


# Recherche des markers du cluster Cycling du Hashtag 1 différent du Hashtag 2
MarkerCycling <- FindMarkers(mouseHashtag23, ident.1 = "Cycling_Hashtag-2", ident.2 = "Cycling_Hashtag-3")

MarkerMait1 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1

MarkerMait17 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17

MarkerMait0 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait0

MarkerCycling %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cycling

Idents(mouseHashtag23)<- mouseHashtag23$old.ident
DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait1)))
DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait17)))
DoHeatmap(mouseHashtag23, features = c(rownames(top10_Cycling)))
DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait0)))

# 
# write.csv(MarkerCycling,file = "./Dossier Graphiques et Tableaux/Hashtag23/GenesmarkerCyclingHashtag2.csv")
# 
# write.csv(MarkerMait0,file = "./Dossier Graphiques et Tableaux/Hashtag23/GenesmarkerMait0Hashtag2.csv")
# 
# write.csv(MarkerMait1,file = "./Dossier Graphiques et Tableaux/Hashtag23/GenesmarkerMait1Hashtag2.csv")
# 
# write.csv(MarkerMait17,file = "./Dossier Graphiques et Tableaux/Hashtag23/GenesmarkerMait17Hashtag2.csv")


# png(file="./Dossier Graphiques et Tableaux/Hashtag23/heatmap gène clusterMAit1 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait1)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/heatmap gène clusterMAit17 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait17)))
# dev.off()
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag23/heatmap gène clusterMait0 par cluster.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag23, features = c(rownames(top10_Mait0)))
# dev.off()


#########################
# Hashtag45 downsampled #
# Fusion des données des Hashtags 4 et 5
mouseHashtag45down = merge(mouseHashtag.subset4,mouseHashtag.subset5)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag45down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag45down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag45down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag45down, pattern = "^Ig")

# VlnPlot(mouseHashtag45down, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag45down <- subset(mouseHashtag45down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag45down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag45down <- subset(mouseHashtag45down, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag45down<-JoinLayers(mouseHashtag45down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag45down <- CellCycleScoring(mouseHashtag45down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Création du score Cytotox
mouseHashtag45down <- AddModuleScore(object = mouseHashtag45down, features = list(Cytotoxsignature), name = "Cytotoxscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag45down)<- mouseHashtag45down$old.ident

# Normalisation des données
mouseHashtag45down <- NormalizeData(mouseHashtag45down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag45down$old.ident)

# On choisit au hasard des cellules du Hashtag 1 pour rammener au nombre de cellule du Hashtag 2; ici 284.
mouseHashtag45down<- subset(mouseHashtag45down, downsample= min(table(mouseHashtag45down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag45down <- FindVariableFeatures(mouseHashtag45down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag45down <- ScaleData(mouseHashtag45down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag45down <- RunPCA(mouseHashtag45down, features = rownames(mouseHashtag45down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag45down <- RunTSNE(mouseHashtag45down, dims = 1:6, perplexity = 50)

# Affichage du t-sne
DimPlot(mouseHashtag45down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag45down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag45down <- FindNeighbors(mouseHashtag45down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag45down <- FindClusters(mouseHashtag45down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag45down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag45down <- FindClusters(mouseHashtag45down, resolution = 0.5)

# Umap selon les 15 dimensions de la PCA
mouseHashtag45down <- RunUMAP(mouseHashtag45down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag45down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot45downgenes <- FeaturePlot(mouseHashtag45down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot45downgenes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot hashtag 4 et 5 downsampled avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot45downgenes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot45downsignatures <- FeaturePlot(mouseHashtag45down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1","Cytotoxscore1"), min.cutoff = 0.2, coord.fixed = TRUE )
Featureplot45downsignatures + DimPlot(mouseHashtag45down, reduction = "umap")


FeaturePlot(mouseHashtag45down, features = c("Gzma","Gzmb"))
# VlnPlot(mouseHashtag45down,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag45down,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag45down,features = c("Xist"), split.by = "old.ident")


# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Graph Lilou 30-04-2024/FeaturePlot hashtag 4 et 5 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot45downsignatures + DimPlot(mouseHashtag45down, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids45down <- c("Mait1", "Mait17","Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids45down) <- levels(mouseHashtag45down)
mouseHashtag45down <- RenameIdents(mouseHashtag45down, new.cluster.ids45down)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap45downlabel <-DimPlot(mouseHashtag45down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap45downlabel

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/UMAP hashtag 4 et 5 downsampled.png", res = 150, width = 600)
# Umap45downlabel
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot45downsignaturessplithashtag <- FeaturePlot(mouseHashtag45down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot45downsignaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot45downsignaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag45down <- DimPlot(mouseHashtag45down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag45down
Umapgrouphashtag45down <- DimPlot(mouseHashtag45down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag45down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/UMAP hashtag 1 et Umap hashtag 2downsampled.png", res = 150, width = 1000)
# Umapsplithashtag45down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/UMAP hashtag 1 et hashtag 2 downsampled.png", res = 150, width = 1000)
# Umapgrouphashtag45down
# dev.off()

VlnPlot(mouseHashtag45down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot Egr2 par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/FeaturePlot Irf5 par cluster et par hashtag 45 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45down,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# FeaturePlot qui superpose les features
FeaturePlot(mouseHashtag45down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
FeaturePlot(mouseHashtag45down, blend = TRUE,features = c("Cd24a","Cd44"))

# ViolinPlot
VlnPlot(mouseHashtag45down,features = c("Cd24a","Cd44"))
VlnPlot(mouseHashtag45down,features = c("Cd24a","Cd44"), split.by = "old.ident")
VlnPlot(mouseHashtag45down,features = c("CD24-prot","CD44-prot"))
VlnPlot(mouseHashtag45down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster45down  <-table(mouseHashtag45down$seurat_clusters,mouseHashtag45down$old.ident)

# Affichage de la table
TableCluster45down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster45down)/2)){
  colonne1 <- c(colonne1 , TableCluster45down[i,1]/(TableCluster45down[i,1]+ TableCluster45down[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster45down[i,2]/(TableCluster45down[i,1]+ TableCluster45down[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster45down [i,1] / min(table(mouseHashtag45down$old.ident)))*100)
  colonne4 <- c(colonne4 , (TableCluster45down [i,2] / min(table(mouseHashtag45down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO45downCluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO45downCluster$`Hashtag 4nb_cell`<-c(TableCluster45down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO45downCluster$`Hashtag 5nb_cell`<-c(TableCluster45down[,2])


# Changement des noms de colonnes
colnames(PercentHTO45downCluster)<- c("Pourcentage intercluster Hashtag 4"," Pourcentage intercluster Hashtag 5","Pourcentage Hashtag 4","Pourcentage Hashtag 5","Hashtag 4nb_cell","Hashtag 5nb_cell")

# Changement des noms de lignes
rownames(PercentHTO45downCluster)<- new.cluster.ids45down

PercentHTO45downCluster
# Enregistrement du tableau de données
# write.csv(PercentHTO45downCluster, file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/PercentHTO45downCluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag45down.markers <- FindAllMarkers(mouseHashtag45down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag45down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_45down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap45down <-DoHeatmap(mouseHashtag45down, features = top10_45down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap45down
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/heatmap hashtag 1 et 2 downsampled par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap45down
# dev.off()

# Séparation selon les hashtags
mouseHashtag45splitdown <- SplitObject(mouseHashtag45down,split.by = "old.ident")
# Hashtag 4
mouseHashtag45split4down <- mouseHashtag45splitdown$`Hashtag-4`
mouseHashtag45split4down <- JoinLayers(mouseHashtag45split4down)

# Hashtag 5
mouseHashtag45split5down <- mouseHashtag45splitdown$`Hashtag-5`
mouseHashtag45split5down <- JoinLayers(mouseHashtag45split5down)

# FeaturePlot de différents gènes du hashtag 4
FeaturePlot(mouseHashtag45split4down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 4 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag45split4down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag45split5down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag45split5down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag4down.markers <- FindAllMarkers(mouseHashtag45split4down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag4down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_4down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap4down <-DoHeatmap(mouseHashtag45split4down, features = top10_4down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap4down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 5 de chaque clusters dont le fold change doit être positif
mouseHashtag5down.markers <- FindAllMarkers(mouseHashtag45split5down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap5down <-DoHeatmap(mouseHashtag45split5down, features = top10_5down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5down

# Affichage des deux heatmaps ensembles
heatmap4down + heatmap5down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/heatmap hashtag 1 et hashtag 2 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap4down + heatmap5down
# dev.off()

# Création de nouveaux identifiants pour la suite
mouseHashtag45down$ClusterGroup <- paste(mouseHashtag45down@active.ident, mouseHashtag45down$old.ident, sep = "_")
Idents(mouseHashtag45down) <- "ClusterGroup"
# Recherche des markers du cluster Mait17 du Hashtag 1 différent du Hashtag 2
MarkerMait17down <- FindMarkers(mouseHashtag45down, ident.1 = "Mait17_Hashtag-4", ident.2 = "Mait17_Hashtag-5")
head(MarkerMait17down, n = 10)

MarkerMait1down<- FindMarkers(mouseHashtag45down, ident.1 = "Mait1_Hashtag-4", ident.2 = "Mait1_Hashtag-5")
head(MarkerMait1down, n = 10)

MarkerCytotoxdown<- FindMarkers(mouseHashtag45down, ident.1 = "Cytotox_Hashtag-4", ident.2 = "Cytotox_Hashtag-5")
head(MarkerCytotoxdown, n = 10)


MarkerMait1down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1down

MarkerMait17down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17down

MarkerCytotoxdown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cytotoxdown


Idents(mouseHashtag45down)<- mouseHashtag45down$old.ident
DoHeatmap(mouseHashtag45down, features = c(rownames(top10_Mait1down)))
DoHeatmap(mouseHashtag45down, features = c(rownames(top10_Mait17down)))
DoHeatmap(mouseHashtag45down, features = c(rownames(top10_Cytotoxdown)))

FeaturePlot(mouseHashtag45down,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag45down, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerMait17down,file = "./Dossier Graphiques et Tableaux/Hashtag45 downsampled/GenesmarkerMait17downsampledHashtag1.csv")
# 
# write.csv(MarkerMait17down,file = "./Dossier Graphiques et Tableaux/Hashtag45 downsampled/GenesmarkerMait17downsampledHashtag1.csv")
# 
# write.csv(MarkerCytotoxdown,file = "./Dossier Graphiques et Tableaux/Hashtag45 downsampled/GenesmarkerCytotoxdownsampledHashtag1.csv")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45 downsampled/heatmap top10Mait1 hashtag 4 et hashtag 5 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag45down, features = c(rownames(top10_Mait1down)))
# dev.off()




#############
# Hashtag45 #
# Fusion des données des Hashtags 4 et 5
mouseHashtag45 = merge(mouseHashtag.subset4,mouseHashtag.subset5)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag45[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag45, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag45[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag45, pattern = "^Ig")

# VlnPlot(mouseHashtag45, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag45 <- subset(mouseHashtag45, subset = percent.mt < 3)

# VlnPlot(mouseHashtag45, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag45 <- subset(mouseHashtag45, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag45<-JoinLayers(mouseHashtag45)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag45 <- CellCycleScoring(mouseHashtag45, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Création du score Cytotox
mouseHashtag45 <- AddModuleScore(object = mouseHashtag45, features = list(Cytotoxsignature), name = "Cytotoxscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag45)<- mouseHashtag45$old.ident

# Normalisation des données
mouseHashtag45 <- NormalizeData(mouseHashtag45)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag45$old.ident)

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag45 <- FindVariableFeatures(mouseHashtag45, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag45 <- ScaleData(mouseHashtag45,verbose = FALSE)

# PCA selon les gènes
mouseHashtag45 <- RunPCA(mouseHashtag45, features = rownames(mouseHashtag45), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag45 <- RunTSNE(mouseHashtag45, dims = 1:6, perplexity = 50)

# Affichage du t-sne
DimPlot(mouseHashtag45)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag45)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag45 <- FindNeighbors(mouseHashtag45, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag45 <- FindClusters(mouseHashtag45, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag45) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag45 <- FindClusters(mouseHashtag45, resolution = 0.5)

# Umap selon les 15 dimensions de la PCA
mouseHashtag45 <- RunUMAP(mouseHashtag45, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag45, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot45genes <- FeaturePlot(mouseHashtag45, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot45genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot hashtag 4 et 5  avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot45genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot45signatures <- FeaturePlot(mouseHashtag45, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1","Cytotoxscore1"), min.cutoff = 0.2)
Featureplot45signatures + DimPlot(mouseHashtag45, reduction = "umap")

# VlnPlot(mouseHashtag45,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag45,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag45,features = c("Xist"), split.by = "old.ident")
# VlnPlot(mouseHashtag45,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Graph Lilou 30-04-2024/FeaturePlot Zbtb16 par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 45 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot hashtag 4 et 5  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot45signatures
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids45 <- c("Mait1a","Mait1b", "Mait17","Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids45) <- levels(mouseHashtag45)
mouseHashtag45 <- RenameIdents(mouseHashtag45, new.cluster.ids45)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap45label <-DimPlot(mouseHashtag45, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap45label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/UMAP hashtag 1 et 2 .png", res = 150, width = 600)
# Umap45label
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot45signaturessplithashtag <- FeaturePlot(mouseHashtag45, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot45signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot45signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag45 <- DimPlot(mouseHashtag45, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag45
Umapgrouphashtag45 <- DimPlot(mouseHashtag45, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag45

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/UMAP hashtag 4 et Umap hashtag 5.png", res = 150, width = 1000)
# Umapsplithashtag45
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/UMAP hashtag 4 et hashtag 5.png", res = 150, width = 1000)
# Umapgrouphashtag45
# dev.off()

VlnPlot(mouseHashtag45,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot Zbtb16 par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot Egr2 par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot Bcl2a1b par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot Nr4a1 par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag45,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/FeaturePlot Irf5 par cluster et par hashtag 45.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag45,features = c("Irf5"),split.by = "old.ident")
# dev.off()


# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag45, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag45, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag45,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag45,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag45,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag45,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster45  <-table(mouseHashtag45$seurat_clusters,mouseHashtag45$old.ident)

# Affichage de la table
TableCluster45

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster45)/2)){
  colonne1 <- c(colonne1 , TableCluster45[i,1]/(TableCluster45[i,1]+ TableCluster45[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster45[i,2]/(TableCluster45[i,1]+ TableCluster45[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster45 [i,1] / sum(TableCluster45[,1]))*100)
  colonne4 <- c(colonne4 , (TableCluster45 [i,2] / sum(TableCluster45[,2]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO45Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO45Cluster$`Hashtag 4nb_cell`<-c(TableCluster45[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO45Cluster$`Hashtag 5nb_cell`<-c(TableCluster45[,2])


# Changement des noms de colonnes
colnames(PercentHTO45Cluster)<- c("Pourcentage intercluster Hashtag 4"," Pourcentage intercluster Hashtag 5","Pourcentage Hashtag 4","Pourcentage Hashtag 5","Hashtag 4nb_cell","Hashtag 5nb_cell")

# Changement des noms de lignes
rownames(PercentHTO45Cluster)<- new.cluster.ids45

PercentHTO45Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO45Cluster, file="./Dossier Graphiques et Tableaux/Hashtag45/PercentHTO45Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag45.markers <- FindAllMarkers(mouseHashtag45, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag45.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_45

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap45 <-DoHeatmap(mouseHashtag45, features = top10_45$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap45
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/heatmap hashtag 4 et 5  par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap45
# dev.off()

# Séparation selon les hashtags
mouseHashtag45split <- SplitObject(mouseHashtag45,split.by = "old.ident")
# Hashtag 4
mouseHashtag45split4 <- mouseHashtag45split$`Hashtag-4`
mouseHashtag45split4 <- JoinLayers(mouseHashtag45split4)

# Hashtag 5
mouseHashtag45split5 <- mouseHashtag45split$`Hashtag-5`
mouseHashtag45split5 <- JoinLayers(mouseHashtag45split5)

# FeaturePlot de différents gènes du hashtag 4
FeaturePlot(mouseHashtag45split4, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 4  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag45split4, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag45split5, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag45split5, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag4.markers <- FindAllMarkers(mouseHashtag45split4, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag4.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_4

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap4 <-DoHeatmap(mouseHashtag45split4, features = top10_4$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap4


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 5 de chaque clusters dont le fold change doit être positif
mouseHashtag5.markers <- FindAllMarkers(mouseHashtag45split5, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap5 <-DoHeatmap(mouseHashtag45split5, features = top10_5$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5

# Affichage des deux heatmaps ensembles
heatmap4 + heatmap5

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag45/heatmap hashtag 4 et hashtag 5 par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap4 + heatmap5
# dev.off()

# Création de nouveaux identifiants pour la suite
mouseHashtag45$ClusterGroup <- paste(mouseHashtag45@active.ident, mouseHashtag45$old.ident, sep = "_")
Idents(mouseHashtag45) <- "ClusterGroup"

MarkerMait17 <- FindMarkers(mouseHashtag45, ident.1 = "Mait17_Hashtag-4", ident.2 = "Mait17_Hashtag-5")
head(MarkerMait17, n = 10)

MarkerMait1a<- FindMarkers(mouseHashtag45, ident.1 = "Mait1a_Hashtag-4", ident.2 = "Mait1a_Hashtag-5")
head(MarkerMait1a, n = 10)

MarkerMait1b <- FindMarkers(mouseHashtag45, ident.1 = "Mait1b_Hashtag-4", ident.2 = "Mait1b_Hashtag-5")
head(MarkerMait1b, n = 10)

MarkerCytotox <- FindMarkers(mouseHashtag45, ident.1 = "Cytotox_Hashtag-4", ident.2 = "Cytotox_Hashtag-5")
head(MarkerCytotox, n = 10)

MarkerMait17 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17

MarkerMait1a %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1a

MarkerMait1b %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1b

MarkerCytotox %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cytotox

Idents(mouseHashtag45)<- mouseHashtag45$old.ident
DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait17)))
DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait1a)))
DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait1b)))
DoHeatmap(mouseHashtag45, features = c(rownames(top10_Cytotox)))

FeaturePlot(mouseHashtag45,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag45, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerMait17,file = "./Dossier Graphiques et Tableaux/Hashtag45/GenesmarkerMait17Hashtag1.csv")

# write.csv(MarkerMait1a,file = "./Dossier Graphiques et Tableaux/Hashtag45/GenesmarkerMait1aHashtag1.csv")
# 
# write.csv(MarkerMait1b,file = "./Dossier Graphiques et Tableaux/Hashtag45/GenesmarkerMait1bHashtag1.csv")
# 
# write.csv(MarkerCytotox,file = "./Dossier Graphiques et Tableaux/Hashtag45/GenesmarkerCytotoxHashtag1.csv")

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/heatmap top10Mait17 hashtag 4 et hashtag 5 .png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait17)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/heatmap top10Mait17 hashtag 4 et hashtag 5 .png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait1a)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag45/heatmap top10Mait17 hashtag 4 et hashtag 5 .png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag45, features = c(rownames(top10_Mait1b)))
# dev.off()


#########################
# Hashtag56 downsampled #
# Fusion des données des Hashtags 5 et 6
mouseHashtag56down = merge(mouseHashtag.subset5,mouseHashtag.subset6)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag56down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag56down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag56down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag56down, pattern = "^Ig")

# VlnPlot(mouseHashtag56down, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag56down <- subset(mouseHashtag56down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag56down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag56down <- subset(mouseHashtag56down, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag56down<-JoinLayers(mouseHashtag56down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag56down <- CellCycleScoring(mouseHashtag56down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score Cytotox
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(Cytotoxsignature), name = "Cytotoxscore")

# Création du score cyclingg2m
mouseHashtag56down <- AddModuleScore(object = mouseHashtag56down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag56down)<- mouseHashtag56down$old.ident

# Normalisation des données
mouseHashtag56down <- NormalizeData(mouseHashtag56down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag56down$old.ident)

# On choisit au hasard des cellules du Hashtag 1 pour rammener au nombre de cellule du Hashtag 2; ici 284.
mouseHashtag56down<- subset(mouseHashtag56down, downsample= min(table(mouseHashtag56down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag56down <- FindVariableFeatures(mouseHashtag56down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag56down <- ScaleData(mouseHashtag56down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag56down <- RunPCA(mouseHashtag56down, features = rownames(mouseHashtag56down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag56down <- RunTSNE(mouseHashtag56down, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag56down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag56down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag56down <- FindNeighbors(mouseHashtag56down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag56down <- FindClusters(mouseHashtag56down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag56down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag56down <- FindClusters(mouseHashtag56down, resolution = 0.3)

# Umap selon les 15 dimensions de la PCA
mouseHashtag56down <- RunUMAP(mouseHashtag56down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag56down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot56downgenes <- FeaturePlot(mouseHashtag56down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot56downgenes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot hashtag 5 et 6 downsampled avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot56downgenes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot56downsignatures <- FeaturePlot(mouseHashtag56down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1","Cytotoxscore1"), min.cutoff = 0.2)
Featureplot56downsignatures + DimPlot(mouseHashtag56down, reduction = "umap")

# VlnPlot(mouseHashtag56down,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag56down,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag56down,features = c("Xist"), split.by = "old.ident")


# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot hashtag 5 et 6 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot56downsignatures + DimPlot(mouseHashtag56down, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids56down <- c("Mait1","Mait17","Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids56down) <- levels(mouseHashtag56down)
mouseHashtag56down <- RenameIdents(mouseHashtag56down, new.cluster.ids56down)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap56downlabel <-DimPlot(mouseHashtag56down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap56downlabel

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/UMAP hashtag 5 et 6 downsampled.png", res = 150, width = 600)
# Umap56downlabel
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot56downsignaturessplithashtag <- FeaturePlot(mouseHashtag56down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot56downsignaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5 et FeaturePlot hashtag 6 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot56downsignaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag56down <- DimPlot(mouseHashtag56down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag56down
Umapgrouphashtag56down <- DimPlot(mouseHashtag56down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag56down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/UMAP hashtag 5 et Umap hashtag 6 downsampled.png", res = 150, width = 1000)
# Umapsplithashtag56down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/UMAP hashtag 5 et hashtag 6 downsampled.png", res = 150, width = 1000)
# Umapgrouphashtag56down
# dev.off()


VlnPlot(mouseHashtag56down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot Egr2 par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/FeaturePlot Irf5 par cluster et par hashtag 56 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56down,features = c("Irf5"),split.by = "old.ident")
# dev.off()


# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag56down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag56down, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag56down,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag56down,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag56down,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag56down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster56down  <-table(mouseHashtag56down$seurat_clusters,mouseHashtag56down$old.ident)

# Affichage de la table
TableCluster56down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster56down)/2)){
  colonne1 <- c(colonne1 , TableCluster56down[i,1]/(TableCluster56down[i,1]+ TableCluster56down[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster56down[i,2]/(TableCluster56down[i,1]+ TableCluster56down[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster56down [i,1] / min(table(mouseHashtag56down$old.ident)))*100)
  colonne4 <- c(colonne4 , (TableCluster56down [i,2] / min(table(mouseHashtag56down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO56downCluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO56downCluster$`Hashtag 5nb_cell`<-c(TableCluster56down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO56downCluster$`Hashtag 6nb_cell`<-c(TableCluster56down[,2])


# Changement des noms de colonnes
colnames(PercentHTO56downCluster)<- c("Pourcentage intercluster Hashtag 5"," Pourcentage intercluster Hashtag 6","Pourcentage Hashtag 5","Pourcentage Hashtag 6","Hashtag 5nb_cell","Hashtag 6nb_cell")

# Changement des noms de lignes
rownames(PercentHTO56downCluster)<- new.cluster.ids56down

PercentHTO56downCluster
# Enregistrement du tableau de données
# write.csv(PercentHTO56downCluster, file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/PercentHTO56downCluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag56down.markers <- FindAllMarkers(mouseHashtag56down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag56down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_56down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap56down <-DoHeatmap(mouseHashtag56down, features = top10_56down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap56down
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/heatmap hashtag 5 et 6 downsampled par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap56down
# dev.off()

# Séparation selon les hashtags
mouseHashtag56splitdown <- SplitObject(mouseHashtag56down,split.by = "old.ident")
# Hashtag 5
mouseHashtag56split5down <- mouseHashtag56splitdown$`Hashtag-5`
mouseHashtag56split5down <- JoinLayers(mouseHashtag56split5down)

# Hashtag 6
mouseHashtag56split6down <- mouseHashtag56splitdown$`Hashtag-6`
mouseHashtag56split6down <- JoinLayers(mouseHashtag56split6down)

# FeaturePlot de différents gènes du hashtag 5
FeaturePlot(mouseHashtag56split5down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag56split5down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag56split6down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 6 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag56split6down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag5down.markers <- FindAllMarkers(mouseHashtag56split5down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap5down <-DoHeatmap(mouseHashtag56split5down, features = top10_5down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 5 de chaque clusters dont le fold change doit être positif
mouseHashtag6down.markers <- FindAllMarkers(mouseHashtag56split6down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag6down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_6down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap6down <-DoHeatmap(mouseHashtag56split6down, features = top10_6down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap6down

# Affichage des deux heatmaps ensembles
heatmap5down + heatmap6down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/heatmap hashtag 5 et hashtag 6 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap5down + heatmap6down
# dev.off()


# Création de nouveaux identifiants pour la suite
mouseHashtag56down$ClusterGroup <- paste(mouseHashtag56down@active.ident, mouseHashtag56down$old.ident, sep = "_")
Idents(mouseHashtag56down) <- "ClusterGroup"

MarkerMait17down <- FindMarkers(mouseHashtag56down, ident.1 = "Mait17_Hashtag-5", ident.2 = "Mait17_Hashtag-6")

MarkerMait1down <- FindMarkers(mouseHashtag56down, ident.1 = "Mait1_Hashtag-5", ident.2 = "Mait1_Hashtag-6")

MarkerCytotoxdown <- FindMarkers(mouseHashtag56down, ident.1 = "Cytotox_Hashtag-5", ident.2 = "Cytotox_Hashtag-6")

MarkerMait1down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1down

MarkerMait17down %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17down

MarkerCytotoxdown %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cytotoxdown

Idents(mouseHashtag56down)<- mouseHashtag56down$old.ident
DoHeatmap(mouseHashtag56down, features = c(rownames(top10_Mait1down)))
DoHeatmap(mouseHashtag56down, features = c(rownames(top10_Mait17down)))
DoHeatmap(mouseHashtag56down, features = c(rownames(top10_Cytotoxdown)))

FeaturePlot(mouseHashtag56down,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag56down, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerMait17down,file = "./Dossier Graphiques et Tableaux/Hashtag56 downsampled/GenesmarkerMait17downsampledHashtag1.csv")
# 
# write.csv(MarkerMait1down,file = "./Dossier Graphiques et Tableaux/Hashtag56 downsampled/GenesmarkerMait1downsampledHashtag1.csv")
# 
# write.csv(MarkerCytotoxdown,file = "./Dossier Graphiques et Tableaux/Hashtag56 downsampled/GenesmarkerCytotoxdownsampledHashtag1.csv")
# 
# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/heatmap top10Mait1 hashtag 5 et hashtag 6 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag56down, features = c(rownames(top10_Mait1down)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag56 downsampled/heatmap top10Mait17 hashtag 5 et hashtag 6 downsampled.png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag56down, features = c(rownames(top10_Mait17down)))
# dev.off()
#############
# Hashtag56 #
# Fusion des données des Hashtags 5 et 6
mouseHashtag56 = merge(mouseHashtag.subset5,mouseHashtag.subset6)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag56[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag56, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag56[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag56, pattern = "^Ig")

# VlnPlot(mouseHashtag56, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag56 <- subset(mouseHashtag56, subset = percent.mt < 3)

# VlnPlot(mouseHashtag56, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag56 <- subset(mouseHashtag56, subset = percent.Ig < 20)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag56<-JoinLayers(mouseHashtag56)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag56 <- CellCycleScoring(mouseHashtag56, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Cytotox
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(Cytotoxsignature), name = "Cytotoxscore")

# Création du score Mait0mini
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag56 <- AddModuleScore(object = mouseHashtag56, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag56)<- mouseHashtag56$old.ident

# Normalisation des données
mouseHashtag56 <- NormalizeData(mouseHashtag56)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag56$old.ident)

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag56 <- FindVariableFeatures(mouseHashtag56, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag56 <- ScaleData(mouseHashtag56,verbose = FALSE)

# PCA selon les gènes
mouseHashtag56 <- RunPCA(mouseHashtag56, features = rownames(mouseHashtag56), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag56 <- RunTSNE(mouseHashtag56, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag56)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag56)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag56 <- FindNeighbors(mouseHashtag56, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag56 <- FindClusters(mouseHashtag56, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag56) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag56 <- FindClusters(mouseHashtag56, resolution = 0.2)

# Umap selon les 15 dimensions de la PCA
mouseHashtag56 <- RunUMAP(mouseHashtag56, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag56, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot56genes <- FeaturePlot(mouseHashtag56, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot56genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot hashtag 5 et 6  avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot56genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot56signatures <- FeaturePlot(mouseHashtag56, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1","Cytotoxscore1"), min.cutoff = 0.2)
Featureplot56signatures + DimPlot(mouseHashtag56, reduction = "umap")

# VlnPlot(mouseHashtag56,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")

# FeaturePlot(mouseHashtag56,features = "Egr2", min.cutoff = 0.2)
# VlnPlot(mouseHashtag56,features = c("Xist"), split.by = "old.ident")


# png(file="./Graph Lilou 30-04-2024/FeaturePlot Egr2 par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Egr2"),split.by = "old.ident")
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot hashtag 1 et 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot56signatures + DimPlot(mouseHashtag56, reduction = "umap")
# dev.off()

# Création des nouveaux nom des clusters manuellement dans l'ordre c(0,1,2,3,4,5,...)
# Attention les clusters peuvent changer d'une fois à l'autre il faut changer l'ordre des noms.
new.cluster.ids56 <- c("Mait1","Mait17","Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids56) <- levels(mouseHashtag56)
mouseHashtag56 <- RenameIdents(mouseHashtag56, new.cluster.ids56)

# Affichage de la Umap avec les noms de cluster afficher sur le graphique
Umap56label <-DimPlot(mouseHashtag56, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap56label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/UMAP hashtag 5 et 6 .png", res = 150, width = 600)
# Umap56label
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot56signaturessplithashtag <- FeaturePlot(mouseHashtag56, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot56signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5 et FeaturePlot hashtag 6  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot56signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag56 <- DimPlot(mouseHashtag56, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag56
Umapgrouphashtag56 <- DimPlot(mouseHashtag56, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag56

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/UMAP hashtag 5 et Umap hashtag 6.png", res = 150, width = 1000)
# Umapsplithashtag56
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/UMAP hashtag 5 et hashtag 6 .png", res = 150, width = 1000)
# Umapgrouphashtag56
# dev.off()

VlnPlot(mouseHashtag56,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot Zbtb16 par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot Egr2 par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot Bcl2a1b par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot Nr4a1 par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag56,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/FeaturePlot Irf5 par cluster et par hashtag 56 .png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag56,features = c("Irf5"),split.by = "old.ident")
# dev.off()


# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag56, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag56, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag56,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag56,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag56,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag56,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster56  <-table(mouseHashtag56$seurat_clusters,mouseHashtag56$old.ident)

# Affichage de la table
TableCluster56

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
for (i in 1:(length(TableCluster56)/2)){
  colonne1 <- c(colonne1 , TableCluster56[i,1]/(TableCluster56[i,1]+ TableCluster56[i,2])*100)
  colonne2 <- c(colonne2 , TableCluster56[i,2]/(TableCluster56[i,1]+ TableCluster56[i,2])*100)
  colonne3 <- c(colonne3 , (TableCluster56 [i,1] / sum(TableCluster56[,1]))*100)
  colonne4 <- c(colonne4 , (TableCluster56 [i,2] / sum(TableCluster56[,2]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]

PercentHTO56Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO56Cluster$`Hashtag 5nb_cell`<-c(TableCluster56[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO56Cluster$`Hashtag 6nb_cell`<-c(TableCluster56[,2])


# Changement des noms de colonnes
colnames(PercentHTO56Cluster)<- c("Pourcentage intercluster Hashtag 5"," Pourcentage intercluster Hashtag 6","Pourcentage Hashtag 5","Pourcentage Hashtag 6","Hashtag 5nb_cell","Hashtag 6nb_cell")

# Changement des noms de lignes
rownames(PercentHTO56Cluster)<- new.cluster.ids56

PercentHTO56Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO56Cluster, file="./Dossier Graphiques et Tableaux/Hashtag56/PercentHTO56Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag56.markers <- FindAllMarkers(mouseHashtag56, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag56.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_56

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap56 <-DoHeatmap(mouseHashtag56, features = top10_56$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap56
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/heatmap hashtag 5 et 6 par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap56
# dev.off()

# Séparation selon les hashtags
mouseHashtag56split <- SplitObject(mouseHashtag56,split.by = "old.ident")
# Hashtag 5
mouseHashtag56split5 <- mouseHashtag56split$`Hashtag-5`
mouseHashtag56split5 <- JoinLayers(mouseHashtag56split5)

# Hashtag 6
mouseHashtag56split6 <- mouseHashtag56split$`Hashtag-6`
mouseHashtag56split6 <- JoinLayers(mouseHashtag56split6)

# FeaturePlot de différents gènes du hashtag 4
FeaturePlot(mouseHashtag56split5, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 5  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag56split5, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag56split6, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 6  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag56split6, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag5.markers <- FindAllMarkers(mouseHashtag56split5, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap5 <-DoHeatmap(mouseHashtag56split5, features = top10_5$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 5 de chaque clusters dont le fold change doit être positif
mouseHashtag6.markers <- FindAllMarkers(mouseHashtag56split6, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag6.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_6

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap6 <-DoHeatmap(mouseHashtag56split6, features = top10_6$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap6

# Affichage des deux heatmaps ensembles
heatmap5 + heatmap6

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag56/heatmap hashtag 5 et hashtag 6 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap5 + heatmap6
# dev.off()

# Création de nouveaux identifiants pour la suite
mouseHashtag56$ClusterGroup <- paste(mouseHashtag56@active.ident, mouseHashtag56$old.ident, sep = "_")
Idents(mouseHashtag56) <- "ClusterGroup"

MarkerMait17 <- FindMarkers(mouseHashtag56, ident.1 = "Mait17_Hashtag-5", ident.2 = "Mait17_Hashtag-6")

MarkerMait1 <- FindMarkers(mouseHashtag56, ident.1 = "Mait1_Hashtag-5", ident.2 = "Mait1_Hashtag-6")

# Recherche des markers du cluster Mait17 du Hashtag 1 différent du Hashtag 2
MarkerCytotox <- FindMarkers(mouseHashtag56, ident.1 = "Cytotox_Hashtag-5", ident.2 = "Cytotox_Hashtag-6")
head(MarkerCytotox, n = 10)

MarkerMait1 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait1

MarkerMait17 %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Mait17

MarkerCytotox %>%
  dplyr::filter(p_val_adj<0.05)%>%
  ungroup() -> top10_Cytotox

Idents(mouseHashtag56)<- mouseHashtag56$old.ident
DoHeatmap(mouseHashtag56, features = c(rownames(top10_Mait1)))
DoHeatmap(mouseHashtag56, features = c(rownames(top10_Mait17)))
DoHeatmap(mouseHashtag56, features = c(rownames(top10_Cytotox)))



FeaturePlot(mouseHashtag56,features = "Xist", split.by = "ident",min.cutoff = 0)
FeaturePlot(mouseHashtag56, features = rownames(top10_Mait1), split.by = "ident", min.cutoff = 0)

# write.csv(MarkerCyclingG2M,file = "./Graph Lilou 30-04-2024/GenesmarkerCyclingG2MHashtag1.csv")

# write.csv(MarkerCyclingS,file = "./Graph Lilou 30-04-2024/GenesmarkerCyclingSHashtag1.csv")

# write.csv(MarkerMait0,file = "./Graph Lilou 30-04-2024/GenesmarkerMait0Hashtag1.csv")

# write.csv(MarkerMait1,file = "./Graph Lilou 30-04-2024/GenesmarkerMAit1Hashtag1.csv")

# write.csv(MarkerMait17a,file = "./Graph Lilou 30-04-2024/GenesmarkerMait17aHashtag1.csv")

# write.csv(MarkerMait17b,file = "./Graph Lilou 30-04-2024/GenesmarkerMait17bHashtag1.csv")

# write.csv(MarkerPrecursor,file = "./Graph Lilou 30-04-2024/GenesmarkerPrecursorHashtag1.csv")

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/heatmap top10Mait1 hashtag 5 et hashtag 6 .png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag56, features = c(rownames(top10_Mait1)))
# dev.off()

# png(file="./Dossier Graphiques et Tableaux/Hashtag56/heatmap top10Mait17 hashtag 5 et hashtag 6 .png", res = 95, height = 1100, width = 1000)
# DoHeatmap(mouseHashtag56, features = c(rownames(top10_Mait17)))
# dev.off()

###########################
# Hashtag 123 downsampled #

mouseHashtag123down = merge(mouseHashtag.subset1,mouseHashtag.subset2)
mouseHashtag123down = merge(mouseHashtag123down,mouseHashtag.subset3)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag123down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag123down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag123down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag123down, pattern = "^Ig")

# VlnPlot(mouseHashtag123down, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag123down <- subset(mouseHashtag123down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag123down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag123down <- subset(mouseHashtag123down, subset = percent.Ig <0.25)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag123down<-JoinLayers(mouseHashtag123down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag123down <- CellCycleScoring(mouseHashtag123down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(precursorsignature), name = "precursorscore")

# # Création du score TCR
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(TCRsignature), name = "TCRscore")

# Création du score TCRmini
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(TCRminisignature), name = "TCRminiscore")

# Création du score Egr2reg
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Egr2regsignature), name = "Egr2regscore")

# Création du score cyclings
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Création du score de la signature des cellules en phase S
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(CyclingGenessignature), name = "Cyclingscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Immature_Asignature), name = "ImmatureAscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Immature_Bsignature), name = "ImmatureBscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Interm_Asignature), name = "IntermAscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Interm_Bsignature), name = "IntermBscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Interm_A_Longsignature), name = "IntermALongscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Interm_B_Longsignature), name = "IntermBLongscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Immature_A_Longsignature), name = "ImmatureALongscore")
# 
# mouseHashtag123down <- AddModuleScore(object = mouseHashtag123down, features = list(Immature_B_Longsignature), name = "ImmatureBLongscore")


# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag123down)<- mouseHashtag123down$old.ident

# Normalisation des données
mouseHashtag123down <- NormalizeData(mouseHashtag123down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag123down$old.ident)

mouseHashtag123down<- subset(mouseHashtag123down, downsample= min(table(mouseHashtag123down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag123down <- FindVariableFeatures(mouseHashtag123down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag123down <- ScaleData(mouseHashtag123down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag123down <- RunPCA(mouseHashtag123down, features = rownames(mouseHashtag123down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag123down <- RunTSNE(mouseHashtag123down, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag123down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag123down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag123down <- FindNeighbors(mouseHashtag123down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag123down <- FindClusters(mouseHashtag123down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag123down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag123down <- FindClusters(mouseHashtag123down, resolution = 0.9)

# Umap selon les 15 dimensions de la PCA
mouseHashtag123down <- RunUMAP(mouseHashtag123down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag123down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot123genesdown <- FeaturePlot(mouseHashtag123down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot123genesdown
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot hashtag 1 et 2 et 3 avec différents gènes downsampled.png", res = 150, height = 1000, width = 1500)
# Featureplot123genesdown
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot123signaturesdown <- FeaturePlot(mouseHashtag123down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot123signaturesdown + DimPlot(mouseHashtag123down, reduction = "umap")

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot hashtag 1 et 2 et 3 avec différentes signatures downsampled.png", res = 150, height = 1000, width = 1500)
# Featureplot123signaturesdown + DimPlot(mouseHashtag123down, reduction = "umap")
# dev.off()

# FeaturePlot(mouseHashtag123down, features = c("Cyclingscore1","cyclingsscore1","cyclingg2mscore1"),min.cutoff = 0)

# FeaturePlot(mouseHashtag123down, features = c("ImmatureAscore1","ImmatureBscore1","IntermAscore1","IntermBscore1","IntermALongscore1","IntermBLongscore1","ImmatureALongscore1","ImmatureBLongscore1"), min.cutoff = 0) + DimPlot(mouseHashtag123down, reduction = "umap")

# FeaturePlot(mouseHashtag123down,features = c("MAIT17score1","MAIT1score1","IntermBLongscore1","MAIT0miniscore1","ImmatureALongscore1","cyclingsscore1", "cyclingg2mscore1")) + DimPlot(mouseHashtag123down, reduction = "umap")

# VlnPlot(mouseHashtag123down,features = "ImmatureAscore1")
# VlnPlot(mouseHashtag123down,features = "TCRscore1", split.by = "old.ident")
# VlnPlot(mouseHashtag123down,features = "TCRminiscore1", split.by = "old.ident")
VlnPlot(mouseHashtag123down,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")
# VlnPlot(mouseHashtag123down,features = "Egr2regscore1", split.by = "old.ident")
# FeaturePlot(mouseHashtag123down, features = "Egr2regscore1")

FeaturePlot(mouseHashtag123down,features = "Egr2", min.cutoff = 0.2)
VlnPlot(mouseHashtag123down,features = c("Xist"), split.by = "old.ident")
VlnPlot(mouseHashtag123down,features = c("Zbtb16"),split.by = "old.ident")

new.cluster.ids123down <- c("Mait1", "Mait17a","Mait17b", "Precursor","Mait0","CyclingG2M","CyclingS")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids123down) <- levels(mouseHashtag123down)
mouseHashtag123down <- RenameIdents(mouseHashtag123down, new.cluster.ids123down)

Umap123labeldown <-DimPlot(mouseHashtag123down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap123labeldown

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/UMAP hashtag 1 et 2 et 3 downsampled .png", res = 150, width = 600)
# Umap123labeldown
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot123signaturessplithashtagdown <- FeaturePlot(mouseHashtag123down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot123signaturessplithashtagdown

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2 et FeaturePlot hashtag 3 avec différentes signatures downsampled.png", res = 150, width = 1000, height = 1000)
# Featureplot123signaturessplithashtagdown
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag123down <- DimPlot(mouseHashtag123down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag123down
Umapgrouphashtag123down <- DimPlot(mouseHashtag123down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag123down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/UMAP hashtag 1 et Umap hashtag 2 et Umap hashtag 3 downsampled.png", res = 150, width = 1000)
# Umapsplithashtag123down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/UMAP hashtag 1 et hashtag 2 et hashtag 3 downsampled .png", res = 150, width = 1000)
# Umapgrouphashtag123down
# dev.off()

VlnPlot(mouseHashtag123down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag123 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot Egr2 par cluster et par hashtag123 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag123 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag123 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/FeaturePlot Irf5 par cluster et par hashtag123 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123down,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag123down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag123down, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag123down,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag123down,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag123down,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag123down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster123down  <-table(mouseHashtag123down$seurat_clusters,mouseHashtag123down$old.ident)

# Affichage de la table
TableCluster123down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
colonne5 <- ""
colonne6 <- ""
for (i in 1:(length(TableCluster123down)/3)){
  colonne1 <- c(colonne1 , TableCluster123down[i,1]/(TableCluster123down[i,1]+ TableCluster123down[i,2] + TableCluster123down[i,3])*100)
  colonne2 <- c(colonne2 , TableCluster123down[i,2]/(TableCluster123down[i,1]+ TableCluster123down[i,2]+ TableCluster123down[i,3])*100)
  colonne3 <- c(colonne3 , TableCluster123down[i,3]/(TableCluster123down[i,1]+ TableCluster123down[i,2]+ TableCluster123down[i,3])*100)
  colonne4 <- c(colonne4 , (TableCluster123down [i,1] / min(table(mouseHashtag123down$old.ident)))*100)
  colonne5 <- c(colonne5 , (TableCluster123down [i,2] / min(table(mouseHashtag123down$old.ident)))*100)
  colonne6 <- c(colonne6 , (TableCluster123down [i,3] / min(table(mouseHashtag123down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]
colonne5 = colonne5[-1]
colonne6 = colonne6[-1]

PercentHTO123Clusterdown <- data.frame(colonne1,colonne2,colonne3,colonne4,colonne5,colonne6)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO123Clusterdown$`Hashtag 1nb_cell`<-c(TableCluster123down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO123Clusterdown$`Hashtag 2nb_cell`<-c(TableCluster123down[,2])

# Ajout du nombre de cellule de hashtag 3 pour chaque cluster
PercentHTO123Clusterdown$`Hashtag 3nb_cell`<-c(TableCluster123down[,3])

# Changement des noms de colonnes
colnames(PercentHTO123Clusterdown)<- c("Pourcentage intercluster Hashtag 1"," Pourcentage intercluster Hashtag 2"," Pourcentage intercluster Hashtag 3","Pourcentage Hashtag 1","Pourcentage Hashtag 2","Pourcentage Hashtag 3","Hashtag 1nb_cell","Hashtag 2nb_cell","Hashtag 3nb_cell")

# Changement des noms de lignes
rownames(PercentHTO123Clusterdown)<- new.cluster.ids123down

PercentHTO123Clusterdown
# Enregistrement du tableau de données
# write.csv(PercentHTO123Clusterdown, file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/PercentHTO123Clusterdownsampled.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag123down.markers <- FindAllMarkers(mouseHashtag123down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag123down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_123down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap123down <-DoHeatmap(mouseHashtag123down, features = top10_123down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap123down
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/heatmap hashtag 1 et 2 et 3 par cluster downsampled.png", res = 150, height = 1800, width = 1500)
# heatmap123down
# dev.off()

# Séparation selon les hashtags
mouseHashtag123downsplit <- SplitObject(mouseHashtag123down,split.by = "old.ident")
# Hashtag 1
mouseHashtag123downsplit1 <- mouseHashtag123downsplit$`Hashtag-1`
mouseHashtag123downsplit1 <- JoinLayers(mouseHashtag123downsplit1)

# Hashtag 2
mouseHashtag123downsplit2 <- mouseHashtag123downsplit$`Hashtag-2`
mouseHashtag123downsplit2 <- JoinLayers(mouseHashtag123downsplit2)

# Hashtag 3
mouseHashtag123downsplit3 <- mouseHashtag123downsplit$`Hashtag-3`
mouseHashtag123downsplit3 <- JoinLayers(mouseHashtag123downsplit3)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag123downsplit1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag123downsplit1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag123downsplit2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag123downsplit2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag1down.markers <- FindAllMarkers(mouseHashtag123downsplit1, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag1down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_1down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap1down <-DoHeatmap(mouseHashtag123downsplit1, features = top10_1down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap1down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag2down.markers <- FindAllMarkers(mouseHashtag123downsplit2, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap2down <-DoHeatmap(mouseHashtag123downsplit2, features = top10_2down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2down

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag3down.markers <- FindAllMarkers(mouseHashtag123downsplit3, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag3down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_3down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap3down <-DoHeatmap(mouseHashtag123downsplit3, features = top10_3down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap3down
# Affichage des deux heatmaps ensembles
heatmap1down + heatmap2down + heatmap3down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123 downsampled/heatmap hashtag 1 et hashtag 2 et hashtag 3 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap1down + heatmap2down + heatmap3down
# dev.off()

###############
# Hashtag 123 #

mouseHashtag123 = merge(mouseHashtag.subset1,mouseHashtag.subset2)
mouseHashtag123 = merge(mouseHashtag123,mouseHashtag.subset3)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag123[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag123, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag123[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag123, pattern = "^Ig")

# VlnPlot(mouseHashtag123, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag123 <- subset(mouseHashtag123, subset = percent.mt < 3)

# VlnPlot(mouseHashtag123, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag123 <- subset(mouseHashtag123, subset = percent.Ig <0.25)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag123<-JoinLayers(mouseHashtag123)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag123 <- CellCycleScoring(mouseHashtag123, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(precursorsignature), name = "precursorscore")

# # Création du score TCR
# mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(TCRsignature), name = "TCRscore")

# Création du score TCRmini
# mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(TCRminisignature), name = "TCRminiscore")

# Création du score Egr2reg
# mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(Egr2regsignature), name = "Egr2regscore")

# Création du score cyclings
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag123 <- AddModuleScore(object = mouseHashtag123, features = list(Cyclingg2msignature), name = "cyclingg2mscore")


# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag123)<- mouseHashtag123$old.ident

# Normalisation des données
mouseHashtag123 <- NormalizeData(mouseHashtag123)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag123$old.ident)

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag123 <- FindVariableFeatures(mouseHashtag123, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag123 <- ScaleData(mouseHashtag123,verbose = FALSE)

# PCA selon les gènes
mouseHashtag123 <- RunPCA(mouseHashtag123, features = rownames(mouseHashtag123), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag123 <- RunTSNE(mouseHashtag123, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag123)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag123)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag123 <- FindNeighbors(mouseHashtag123, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag123 <- FindClusters(mouseHashtag123, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag123) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag123 <- FindClusters(mouseHashtag123, resolution = 0.8)

# Umap selon les 15 dimensions de la PCA
mouseHashtag123 <- RunUMAP(mouseHashtag123, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag123, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot123genes <- FeaturePlot(mouseHashtag123, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot123genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot hashtag 1 et 2 et 3 avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot123genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot123signatures <- FeaturePlot(mouseHashtag123, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","precursorscore1"), min.cutoff = 0.2)
Featureplot123signatures + DimPlot(mouseHashtag123, reduction = "umap")
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot hashtag 1 et 2 et 3 avec différentes signatures.png", res = 150, height = 1000, width = 1500)
# Featureplot123signatures + DimPlot(mouseHashtag123, reduction = "umap")
# dev.off()

# FeaturePlot(mouseHashtag123, features = c("Cyclingscore1","cyclingsscore1","cyclingg2mscore1"),min.cutoff = 0)

# FeaturePlot(mouseHashtag123, features = c("ImmatureAscore1","ImmatureBscore1","IntermAscore1","IntermBscore1","IntermALongscore1","IntermBLongscore1","ImmatureALongscore1","ImmatureBLongscore1"), min.cutoff = 0) + DimPlot(mouseHashtag123, reduction = "umap")

# FeaturePlot(mouseHashtag123,features = c("MAIT17score1","MAIT1score1","IntermBLongscore1","MAIT0miniscore1","ImmatureALongscore1","cyclingsscore1", "cyclingg2mscore1")) + DimPlot(mouseHashtag123, reduction = "umap")

# VlnPlot(mouseHashtag123,features = "ImmatureAscore1")
# VlnPlot(mouseHashtag123,features = "TCRscore1", split.by = "old.ident")
# VlnPlot(mouseHashtag123,features = "TCRminiscore1", split.by = "old.ident")
VlnPlot(mouseHashtag123,features = c("Egr3","Egr2","Nr4a1"),split.by = "old.ident")
# VlnPlot(mouseHashtag123,features = "Egr2regscore1", split.by = "old.ident")
# FeaturePlot(mouseHashtag123, features = "Egr2regscore1")

FeaturePlot(mouseHashtag123,features = "Egr2", min.cutoff = 0.2)
VlnPlot(mouseHashtag123,features = c("Xist"), split.by = "old.ident")
VlnPlot(mouseHashtag123,features = c("Zbtb16"),split.by = "old.ident")

new.cluster.ids123 <- c("Mait1", "Mait17a", "Mait0", "Mait17b0","Precursor","CyclingG2M","CyclingS")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids123) <- levels(mouseHashtag123)
mouseHashtag123 <- RenameIdents(mouseHashtag123, new.cluster.ids123)

Umap123label <-DimPlot(mouseHashtag123, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap123label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/UMAP hashtag 1 et 2 et 3.png", res = 150, width = 600)
# Umap123label
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot123signaturessplithashtag <- FeaturePlot(mouseHashtag123, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot123signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot123signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag123 <- DimPlot(mouseHashtag123, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag123
Umapgrouphashtag123 <- DimPlot(mouseHashtag123, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag123

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/UMAP hashtag 1 et Umap hashtag 2 et Umap hashtag 3.png", res = 150, width = 1000)
# Umapsplithashtag123
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/UMAP hashtag 1 et hashtag 2 et hashtag 3 .png", res = 150, width = 1000)
# Umapgrouphashtag123
# dev.off()

VlnPlot(mouseHashtag123,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot Zbtb16 par cluster et par hashtag123.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot Egr2 par cluster et par hashtag123.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot Bcl2a1b par cluster et par hashtag123.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot Nr4a1 par cluster et par hashtag123.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag123,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag123/FeaturePlot Irf5 par cluster et par hashtag123.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag123,features = c("Irf5"),split.by = "old.ident")
# dev.off()

# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag123, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag123, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag123,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag123,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag123,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag123,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster123  <-table(mouseHashtag123$seurat_clusters,mouseHashtag123$old.ident)

# Affichage de la table
TableCluster123

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
colonne5 <- ""
colonne6 <- ""
for (i in 1:(length(TableCluster123)/3)){
  colonne1 <- c(colonne1 , TableCluster123[i,1]/(TableCluster123[i,1]+ TableCluster123[i,2] + TableCluster123[i,3])*100)
  colonne2 <- c(colonne2 , TableCluster123[i,2]/(TableCluster123[i,1]+ TableCluster123[i,2]+ TableCluster123[i,3])*100)
  colonne3 <- c(colonne3 , TableCluster123[i,3]/(TableCluster123[i,1]+ TableCluster123[i,2]+ TableCluster123[i,3])*100)
  colonne4 <- c(colonne4 , (TableCluster123 [i,1] / sum(TableCluster123[,1]))*100)
  colonne5 <- c(colonne5 , (TableCluster123 [i,2] / sum(TableCluster123[,2]))*100)
  colonne6 <- c(colonne6 , (TableCluster123 [i,3] / sum(TableCluster123[,3]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]
colonne5 = colonne5[-1]
colonne6 = colonne6[-1]

PercentHTO123Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4,colonne5,colonne6)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO123Cluster$`Hashtag 1nb_cell`<-c(TableCluster123[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO123Cluster$`Hashtag 2nb_cell`<-c(TableCluster123[,2])

# Ajout du nombre de cellule de hashtag 3 pour chaque cluster
PercentHTO123Cluster$`Hashtag 3nb_cell`<-c(TableCluster123[,3])

# Changement des noms de colonnes
colnames(PercentHTO123Cluster)<- c("Pourcentage intercluster Hashtag 1"," Pourcentage intercluster Hashtag 2"," Pourcentage intercluster Hashtag 3","Pourcentage Hashtag 1","Pourcentage Hashtag 2","Pourcentage Hashtag 3","Hashtag 1nb_cell","Hashtag 2nb_cell","Hashtag 3nb_cell")

# Changement des noms de lignes
rownames(PercentHTO123Cluster)<- new.cluster.ids123

PercentHTO123Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO123Cluster, file="./Dossier Graphiques et Tableaux/Hashtag123/PercentHTO123Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag123.markers <- FindAllMarkers(mouseHashtag123, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag123.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_123

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap123 <-DoHeatmap(mouseHashtag123, features = top10_123$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap123
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/heatmap hashtag 1 et 2 et 3  par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap123
# dev.off()

# Séparation selon les hashtags
mouseHashtag123split <- SplitObject(mouseHashtag123,split.by = "old.ident")
# Hashtag 1
mouseHashtag123split1 <- mouseHashtag123split$`Hashtag-1`
mouseHashtag123split1 <- JoinLayers(mouseHashtag123split1)

# Hashtag 2
mouseHashtag123split2 <- mouseHashtag123split$`Hashtag-2`
mouseHashtag123split2 <- JoinLayers(mouseHashtag123split2)

# Hashtag 3
mouseHashtag123split3 <- mouseHashtag123split$`Hashtag-3`
mouseHashtag123split3 <- JoinLayers(mouseHashtag123split3)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag123split1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag123split1, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag123split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag123split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag1.markers <- FindAllMarkers(mouseHashtag123split1, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag1.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_1

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap1 <-DoHeatmap(mouseHashtag123split1, features = top10_1$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap1


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag2.markers <- FindAllMarkers(mouseHashtag123split2, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag2.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_2

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap2 <-DoHeatmap(mouseHashtag123split2, features = top10_2$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap2

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag3.markers <- FindAllMarkers(mouseHashtag123split3, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag3.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_3

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap3 <-DoHeatmap(mouseHashtag123split3, features = top10_3$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap3
# Affichage des deux heatmaps ensembles
heatmap1 + heatmap2 + heatmap3

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag123/heatmap hashtag 1 et hashtag 2 et hashtag 3  par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap1 + heatmap2 + heatmap3
# dev.off()
###########################
# Hashtag 456 downsampled #

mouseHashtag456down = merge(mouseHashtag.subset4,mouseHashtag.subset5)
mouseHashtag456down = merge(mouseHashtag456down,mouseHashtag.subset6)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag456down[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag456down, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag456down[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag456down, pattern = "^Ig")

# VlnPlot(mouseHashtag456down, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag456down <- subset(mouseHashtag456down, subset = percent.mt < 3)

# VlnPlot(mouseHashtag456down, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag456down <- subset(mouseHashtag456down, subset = percent.Ig <0.25)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag456down<-JoinLayers(mouseHashtag456down)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag456down <- CellCycleScoring(mouseHashtag456down, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(precursorsignature), name = "precursorscore")

# Création du score cyclings
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Création du score Cytotox
mouseHashtag456down <- AddModuleScore(object = mouseHashtag456down, features = list(Cytotoxsignature), name = "Cytotoxscore")


# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag456down)<- mouseHashtag456down$old.ident

# Normalisation des données
mouseHashtag456down <- NormalizeData(mouseHashtag456down)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag456down$old.ident)

# On choisit au hasard des cellules du Hashtag 1 pour rammener au nombre de cellule du Hashtag 2; ici 284.
mouseHashtag456down<- subset(mouseHashtag456down, downsample= min(table(mouseHashtag456down$old.ident)))

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag456down <- FindVariableFeatures(mouseHashtag456down, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag456down <- ScaleData(mouseHashtag456down,verbose = FALSE)

# PCA selon les gènes
mouseHashtag456down <- RunPCA(mouseHashtag456down, features = rownames(mouseHashtag456down), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag456down <- RunTSNE(mouseHashtag456down, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag456down)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag456down)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag456down <- FindNeighbors(mouseHashtag456down, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag456down <- FindClusters(mouseHashtag456down, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag456down) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag456down <- FindClusters(mouseHashtag456down, resolution = 0.5)

# Umap selon les 15 dimensions de la PCA
mouseHashtag456down <- RunUMAP(mouseHashtag456down, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag456down, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot456genes <- FeaturePlot(mouseHashtag456down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot456genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot hashtag 4 et 5 et 6 avec différents gènes ownsampled.png", res = 150, height = 1000, width = 1500)
# Featureplot456genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot456signatures <- FeaturePlot(mouseHashtag456down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","Cytotoxscore1"), min.cutoff = 0.2)
Featureplot456signatures + DimPlot(mouseHashtag456down, reduction = "umap")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot hashtag 4 et 5 et 6 avec différentes signatures downsampled.png", res = 150, height = 1000, width = 1500)
# Featureplot456signatures + DimPlot(mouseHashtag456down, reduction = "umap")
# dev.off()

FeaturePlot(mouseHashtag456down,features = "Egr2", min.cutoff = 0.2)
VlnPlot(mouseHashtag456down,features = c("Xist"), split.by = "old.ident")
VlnPlot(mouseHashtag456down,features = c("Zbtb16"),split.by = "old.ident")

new.cluster.ids456down <- c("Mait1","MAit17","Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids456down) <- levels(mouseHashtag456down)
mouseHashtag456down <- RenameIdents(mouseHashtag456down, new.cluster.ids456down)

Umap456downlabel <-DimPlot(mouseHashtag456down, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap456downlabel

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/UMAP hashtag 1 et 2 et 3 downsampled.png", res = 150, width = 600)
# Umap456downlabel
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot456downsignaturessplithashtag <- FeaturePlot(mouseHashtag456down, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot456downsignaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2 downsampled avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot456downsignaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag456down <- DimPlot(mouseHashtag456down, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag456down
Umapgrouphashtag456down <- DimPlot(mouseHashtag456down, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag456down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/UMAP hashtag 1 et Umap hashtag 2downsampled.png", res = 150, width = 1000)
# Umapsplithashtag456down
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/UMAP hashtag 1 et hashtag 2 downsampled.png", res = 150, width = 1000)
# Umapgrouphashtag456down
# dev.off()

VlnPlot(mouseHashtag456down,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot Zbtb16 par cluster et par hashtag456 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456down,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456down,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot Egr2 par cluster et par hashtag456 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456down,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456down,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot Bcl2a1b par cluster et par hashtag456 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456down,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456down,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot Nr4a1 par cluster et par hashtag456 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456down,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456down,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/FeaturePlot Irf5 par cluster et par hashtag456 downsampled.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456down,features = c("Irf5"),split.by = "old.ident")
# dev.off()



# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag456down, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag456down, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag456down,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag456down,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag456down,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag456down,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster456down  <-table(mouseHashtag456down$seurat_clusters,mouseHashtag456down$old.ident)

# Affichage de la table
TableCluster456down

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
colonne5 <- ""
colonne6 <- ""
for (i in 1:(length(TableCluster456down)/3)){
  colonne1 <- c(colonne1 , TableCluster456down[i,1]/(TableCluster456down[i,1]+ TableCluster456down[i,2] + TableCluster456down[i,3])*100)
  colonne2 <- c(colonne2 , TableCluster456down[i,2]/(TableCluster456down[i,1]+ TableCluster456down[i,2]+ TableCluster456down[i,3])*100)
  colonne3 <- c(colonne3 , TableCluster456down[i,3]/(TableCluster456down[i,1]+ TableCluster456down[i,2]+ TableCluster456down[i,3])*100)
  colonne4 <- c(colonne4 , (TableCluster456down [i,1] / min(table(mouseHashtag456down$old.ident)))*100)
  colonne5 <- c(colonne5 , (TableCluster456down [i,2] / min(table(mouseHashtag456down$old.ident)))*100)
  colonne6 <- c(colonne6 , (TableCluster456down [i,3] / min(table(mouseHashtag456down$old.ident)))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]
colonne5 = colonne5[-1]
colonne6 = colonne6[-1]

PercentHTO456downCluster <- data.frame(colonne1,colonne2,colonne3,colonne4,colonne5,colonne6)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO456downCluster$`Hashtag 4nb_cell`<-c(TableCluster456down[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO456downCluster$`Hashtag 5nb_cell`<-c(TableCluster456down[,2])

# Ajout du nombre de cellule de hashtag 3 pour chaque cluster
PercentHTO456downCluster$`Hashtag 6nb_cell`<-c(TableCluster456down[,3])

# Changement des noms de colonnes
colnames(PercentHTO456downCluster)<- c("Pourcentage intercluster Hashtag 4"," Pourcentage intercluster Hashtag 5"," Pourcentage intercluster Hashtag 6","Pourcentage Hashtag 4","Pourcentage Hashtag 5","Pourcentage Hashtag 6","Hashtag 4nb_cell","Hashtag 5nb_cell","Hashtag 5nb_cell")

# Changement des noms de lignes
rownames(PercentHTO456downCluster)<- new.cluster.ids456down

PercentHTO456downCluster
# Enregistrement du tableau de données
# write.csv(PercentHTO456downCluster, file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/PercentHTO456downCluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag456down.markers <- FindAllMarkers(mouseHashtag456down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag456down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_456down

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap456down <-DoHeatmap(mouseHashtag456down, features = top10_456down$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap456down
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/heatmap hashtag 1 et 2 et 3 downsampled par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap456down
# dev.off()

# Séparation selon les hashtags
mouseHashtag456splitdown <- SplitObject(mouseHashtag456down,split.by = "old.ident")
# Hashtag 4
mouseHashtag456split4down <- mouseHashtag456splitdown$`Hashtag-4`
mouseHashtag456split4down <- JoinLayers(mouseHashtag456split4down)

# Hashtag 2
mouseHashtag456split5down <- mouseHashtag456splitdown$`Hashtag-5`
mouseHashtag456split5down <- JoinLayers(mouseHashtag456split5down)

# Hashtag 3
mouseHashtag456split6down <- mouseHashtag456splitdown$`Hashtag-6`
mouseHashtag456split6down <- JoinLayers(mouseHashtag456split6down)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag456split4down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 4 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag456split4down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag456split5down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2 downsampled avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag456split2down, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag4down.markers <- FindAllMarkers(mouseHashtag456split4down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag4down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_4down

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap4down <-DoHeatmap(mouseHashtag456split4down, features = top10_4down$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap4down


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag5down.markers <- FindAllMarkers(mouseHashtag456split5down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap5down <-DoHeatmap(mouseHashtag456split5down, features = top10_5down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5down

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag6down.markers <- FindAllMarkers(mouseHashtag456split6down, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag6down.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_6down

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap6down <-DoHeatmap(mouseHashtag456split6down, features = top10_6down$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap6down
# Affichage des deux heatmaps ensembles
heatmap4down + heatmap5down + heatmap6down

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456 downsampled/heatmap hashtag 4 et hashtag 5 et hashtag 6 downsampled par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap4down + heatmap5down + heatmap6down
# dev.off()
###############
# Hashtag 456 #

mouseHashtag456 = merge(mouseHashtag.subset4,mouseHashtag.subset5)
mouseHashtag456 = merge(mouseHashtag456,mouseHashtag.subset6)

# Ajout des pourcentages de gènes mitochondriales par cellule
mouseHashtag456[["percent.mt"]] <- PercentageFeatureSet(mouseHashtag456, pattern = "^mt-")
# Ajout des pourcentages de gènes d'immunoglobuline par cellule
mouseHashtag456[["percent.Ig"]] <- PercentageFeatureSet(mouseHashtag456, pattern = "^Ig")

# VlnPlot(mouseHashtag456, features = "percent.mt")
# Retire les cellules ayant un pourcentage de gène mitochondriales supérieur au threshold; ici 5.
mouseHashtag456 <- subset(mouseHashtag456, subset = percent.mt < 3)

# VlnPlot(mouseHashtag456, features = "percent.Ig")
# Retire les cellules ayant un pourcentage de gène d'immunoglobuline supérieur au threshold; ici 20.
mouseHashtag456 <- subset(mouseHashtag456, subset = percent.Ig <0.25)


# Fusion des différents layers pour lancer les commandes suivantes
mouseHashtag456<-JoinLayers(mouseHashtag456)
# Assigne une phrase de cycle cellulaire à chaque cellule selon les gènes signatures
mouseHashtag456 <- CellCycleScoring(mouseHashtag456, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Création du score Mait0
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(MAIT0signature), name = "MAIT0score")

# Création du score Mait1
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(MAIT1signature), name = "MAIT1score")

# Création du score Mait17
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(MAIT17signature), name = "MAIT17score")

# Création du score Mait0mini
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(precursorsignature), name = "precursorscore")

# # Création du score TCR
# mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(TCRsignature), name = "TCRscore")

# Création du score TCRmini
# mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(TCRminisignature), name = "TCRminiscore")

# Création du score Egr2reg
# mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(Egr2regsignature), name = "Egr2regscore")

# Création du score cyclings
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(cyclingssignature), name = "cyclingsscore")

# Création du score cyclingg2m
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(Cyclingg2msignature), name = "cyclingg2mscore")

# Création du score Cytotox
mouseHashtag456 <- AddModuleScore(object = mouseHashtag456, features = list(Cytotoxsignature), name = "Cytotoxscore")


# Changement des Idents de chaque cellule pour la suite du script
Idents(mouseHashtag456)<- mouseHashtag456$old.ident

# Normalisation des données
mouseHashtag456 <- NormalizeData(mouseHashtag456)

# Affiche le nombre de cellules de chacun des hashtags
table(mouseHashtag456$old.ident)

# On recherche les gènes qui changent parmis les cellules; ici 5000
mouseHashtag456 <- FindVariableFeatures(mouseHashtag456, selection.method = "vst",nfeatures = 5000)

# Mise à l'échelle des données
mouseHashtag456 <- ScaleData(mouseHashtag456,verbose = FALSE)

# PCA selon les gènes
mouseHashtag456 <- RunPCA(mouseHashtag456, features = rownames(mouseHashtag456), approx = FALSE)

# t-sne selon les composantes 1 à 6 de la PCA
mouseHashtag456 <- RunTSNE(mouseHashtag456, dims = 1:6, perplexity = 100)

# Affichage du t-sne
DimPlot(mouseHashtag456)

# Affichage de la déviation standard selon les composantes principales de la PCA
ElbowPlot(mouseHashtag456)
# On choisit le dernier point avant que la déviation standard ne change plus lorsque l'on augmente la composante principale. Pour les données de scRNAseq il ne faut pas prendre une dimension trop petite car on peut écarter des données. Ici, 15.

# Recherche des points voisins basé sur les 15 dimensions 
mouseHashtag456 <- FindNeighbors(mouseHashtag456, dims = 1:15)

# Recherche des clusters sur toutes les résolutions entre 0 et 1
mouseHashtag456 <- FindClusters(mouseHashtag456, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9))

# Constructions de l'arbre de clusterisation
clustree(mouseHashtag456) 
# On garde pour la résolution qui ne change pas au rang suivant; ici 0.4
mouseHashtag456 <- FindClusters(mouseHashtag456, resolution = 0.1)

# Umap selon les 15 dimensions de la PCA
mouseHashtag456 <- RunUMAP(mouseHashtag456, dims = 1:15)

# Affichage de la Umap
DimPlot(mouseHashtag456, reduction = "umap")

# FeaturePlot de différents gènes afin de déterminer les clusters
Featureplot456genes <- FeaturePlot(mouseHashtag456, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0)

Featureplot456genes
# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot hashtag 1 et 2 et 3 avec différents gènes.png", res = 150, height = 1000, width = 1500)
# Featureplot456genes
# dev.off()

# FeaturePlot des différentes signatures afin de déterminer les clusters
Featureplot456signatures <- FeaturePlot(mouseHashtag456, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "cyclingsscore1","cyclingg2mscore1","Cytotoxscore1"), min.cutoff = 0.2)
Featureplot456signatures + DimPlot(mouseHashtag456, reduction = "umap")

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot hashtag 1 et 2 et 3 avec différentes signatures.png", res = 150, height = 1000, width = 1500)
# Featureplot456signatures + DimPlot(mouseHashtag456, reduction = "umap")
# dev.off()


new.cluster.ids456 <- c("Mait1","Mait17", "Cytotox")
# Renomme les clusters selon les nouveaux noms
names(new.cluster.ids456) <- levels(mouseHashtag456)
mouseHashtag456 <- RenameIdents(mouseHashtag456, new.cluster.ids456)

Umap456label <-DimPlot(mouseHashtag456, reduction = "umap", label = TRUE, pt.size = 1.5)
Umap456label

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/UMAP hashtag 1 et 2 .png", res = 150, width = 600)
# Umap456label
# dev.off()

# FeaturePlot des différentes signatures géniques séparé selon les hashtags
Featureplot456signaturessplithashtag <- FeaturePlot(mouseHashtag456, features = c("MAIT17score1", "MAIT1score1", "MAIT0score1", "Cyclingscore1"), min.cutoff = 0, split.by = "old.ident")
Featureplot456signaturessplithashtag

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 1 et FeaturePlot hashtag 2  avec différentes signatures.png", res = 150, width = 1000, height = 1000)
# Featureplot456signaturessplithashtag
# dev.off()

# Affichage de la Umap séparé selon les hashtags
Umapsplithashtag456 <- DimPlot(mouseHashtag456, reduction = "umap", pt.size = 1.5, split.by = "old.ident")
Umapsplithashtag456
Umapgrouphashtag456 <- DimPlot(mouseHashtag456, reduction = "umap", pt.size = 1.5, group.by = "old.ident")
Umapgrouphashtag456

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/UMAP hashtag 4 et Umap hashtag 5 et Umap hashtag 6.png", res = 150, width = 1000)
# Umapsplithashtag456
# dev.off()

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/UMAP hashtag 4 et hashtag 5 et haashtag 6.png", res = 150, width = 1000)
# Umapgrouphashtag456
# dev.off()


VlnPlot(mouseHashtag456,features = c("Zbtb16"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot Zbtb16 par cluster et par hashtag456.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456,features = c("Zbtb16"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456,features = c("Egr2"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot Egr2 par cluster et par hashtag456.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456,features = c("Egr2"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456,features = c("Bcl2a1b"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot Bcl2a1b par cluster et par hashtag456.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456,features = c("Bcl2a1b"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456,features = c("Nr4a1"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot Nr4a1 par cluster et par hashtag456.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456,features = c("Nr4a1"),split.by = "old.ident")
# dev.off()

VlnPlot(mouseHashtag456,features = c("Irf5"),split.by = "old.ident")

# png(file="./Dossier Graphiques et Tableaux/Hashtag456/FeaturePlot Irf5 par cluster et par hashtag456.png", res = 150, width = 1000, height = 800)
# VlnPlot(mouseHashtag456,features = c("Irf5"),split.by = "old.ident")
# dev.off()



# # FeaturePlot qui superpose les features
# FeaturePlot(mouseHashtag456, blend = TRUE,features = c("CD24-prot","CD44-prot"))
# FeaturePlot(mouseHashtag456, blend = TRUE,features = c("Cd24a","Cd44"))
# 
# # ViolinPlot
# VlnPlot(mouseHashtag456,features = c("Cd24a","Cd44"))
# VlnPlot(mouseHashtag456,features = c("Cd24a","Cd44"), split.by = "old.ident")
# VlnPlot(mouseHashtag456,features = c("CD24-prot","CD44-prot"))
# VlnPlot(mouseHashtag456,features = c("CD24-prot","CD44-prot"), split.by = "old.ident")

# Création d'une table entre le numéro du cluster et le hashtag
TableCluster456  <-table(mouseHashtag456$seurat_clusters,mouseHashtag456$old.ident)

# Affichage de la table
TableCluster456

# Création du dataframe contenant les informations du nombre de cellule par cluster et par hashtag
colonne1 <- ""
colonne2 <- ""
colonne3 <- ""
colonne4 <- ""
colonne5 <- ""
colonne6 <- ""
for (i in 1:(length(TableCluster456)/3)){
  colonne1 <- c(colonne1 , (TableCluster456[i,1]/(TableCluster456[i,1]+ TableCluster456[i,2]+ TableCluster456[i,3]))*100)
  colonne2 <- c(colonne2 , (TableCluster456[i,2]/(TableCluster456[i,1]+ TableCluster456[i,2]+ TableCluster456[i,3]))*100)
  colonne3 <- c(colonne3 , (TableCluster456[i,3]/(TableCluster456[i,1]+ TableCluster456[i,2]+ TableCluster456[i,3]))*100)
  colonne4 <- c(colonne4 , (TableCluster456 [i,1] / sum(TableCluster456[,1]))*100)
  colonne5 <- c(colonne5 , (TableCluster456 [i,2] / sum(TableCluster456[,2]))*100)
  colonne6 <- c(colonne6 , (TableCluster456 [i,3] / sum(TableCluster456[,3]))*100)
}
colonne1 = colonne1[-1]
colonne2 = colonne2[-1]
colonne3 = colonne3[-1]
colonne4 = colonne4[-1]
colonne5 = colonne5[-1]
colonne6 = colonne6[-1]

PercentHTO456Cluster <- data.frame(colonne1,colonne2,colonne3,colonne4,colonne5,colonne6)

# Ajout du nombre de cellule de hashtag 1 pour chaque cluster
PercentHTO456Cluster$`Hashtag 4nb_cell`<-c(TableCluster456[,1])

# Ajout du nombre de cellule de hashtag 2 pour chaque cluster
PercentHTO456Cluster$`Hashtag 5nb_cell`<-c(TableCluster456[,2])

# Ajout du nombre de cellule de hashtag 3 pour chaque cluster
PercentHTO456Cluster$`Hashtag 6nb_cell`<-c(TableCluster456[,3])

# Changement des noms de colonnes
colnames(PercentHTO456Cluster)<- c("Pourcentage intercluster Hashtag 4"," Pourcentage intercluster Hashtag 5"," Pourcentage intercluster Hashtag 6","Pourcentage Hashtag 4","Pourcentage Hashtag 5","Pourcentage Hashtag 6","Hashtag 4nb_cell","Hashtag 5nb_cell","Hashtag 5nb_cell")

# Changement des noms de lignes
rownames(PercentHTO456Cluster)<- new.cluster.ids456

PercentHTO456Cluster
# Enregistrement du tableau de données
# write.csv(PercentHTO456Cluster, file="./Dossier Graphiques et Tableaux/Hashtag456/PercentHTO456Cluster.csv")

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 et 2 de chaque clusters dont le fold change doit être positif
mouseHashtag456.markers <- FindAllMarkers(mouseHashtag456, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag456.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_456

# Heatmap des hashtags 1 et 2 des 10 meilleurs gènes de chaque cluster
heatmap456 <-DoHeatmap(mouseHashtag456, features = top10_456$gene, combine = TRUE, angle = 90, group.by = "ident") + NoLegend()

# Affiche la heatmap
heatmap456
# # Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/heatmap hashtag 4 et 5 et 6  par cluster.png", res = 150, height = 1800, width = 1500)
# heatmap456
# dev.off()

# Séparation selon les hashtags
mouseHashtag456split <- SplitObject(mouseHashtag456,split.by = "old.ident")
# Hashtag 4
mouseHashtag456split4 <- mouseHashtag456split$`Hashtag-4`
mouseHashtag456split4 <- JoinLayers(mouseHashtag456split4)

# Hashtag 2
mouseHashtag456split5 <- mouseHashtag456split$`Hashtag-5`
mouseHashtag456split5 <- JoinLayers(mouseHashtag456split5)

# Hashtag 3
mouseHashtag456split6 <- mouseHashtag456split$`Hashtag-6`
mouseHashtag456split6 <- JoinLayers(mouseHashtag456split6)

# FeaturePlot de différents gènes du hashtag 1
FeaturePlot(mouseHashtag456split4, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 4  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag456split4, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# FeaturePlot de différents gènes du hashtag 2
FeaturePlot(mouseHashtag456split5, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)

# Enregistrement de l'image
# png(file="FeaturePlot hashtag 2  avec différents gènes.png", res = 150, height = 1500, width = 2000)
# FeaturePlot(mouseHashtag456split2, features = c("Zbtb16","Ccr7", "Ccr9", "Cd3e","Tbx21","Klrk1","Rorc","Ccr6","Dntt","Rag1","Rag2","Mki67", "Sell","Klf2"), min.cutoff = 0, ncol = 4)
# dev.off()

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 1 de chaque clusters dont le fold change doit être positif
mouseHashtag4.markers <- FindAllMarkers(mouseHashtag456split4, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag4.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_4

# Heatmap des hashtags 1 des 10 meilleurs gènes de chaque cluster
heatmap4 <-DoHeatmap(mouseHashtag456split4, features = top10_4$gene, combine = TRUE,angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap4


# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag5.markers <- FindAllMarkers(mouseHashtag456split5, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag5.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_5

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap5 <-DoHeatmap(mouseHashtag456split5, features = top10_5$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap5

# Recherche des différents markers (gènes les plus marquants) sur les données du hashtag 2 de chaque clusters dont le fold change doit être positif
mouseHashtag6.markers <- FindAllMarkers(mouseHashtag456split6, only.pos = TRUE)

# Parmis tout les marqueurs on les trie par cluster, on conserve ceux dont le fold change est supérieur à 1, on garde les 10 premiers de chaque clusters
mouseHashtag6.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_6

# Heatmap des hashtags 2 des 10 meilleurs gènes de chaque cluster
heatmap6 <-DoHeatmap(mouseHashtag456split6, features = top10_6$gene, combine = TRUE, angle = 90) + NoLegend()

# Affichage de la heatmap
heatmap6
# Affichage des deux heatmaps ensembles
heatmap4 + heatmap5 + heatmap6

# Enregistrement de l'image
# png(file="./Dossier Graphiques et Tableaux/Hashtag456/heatmap hashtag 4 et hashtag 5 et hashtag 6 par cluster.png", res = 95, height = 1100, width = 1000)
# heatmap4 + heatmap5 + heatmap6
# dev.off()
##########
# SCENIC #

# Création des données 
# singleCellMatrix <- Seurat::Read10X(data.dir="./filtered_feature_bc_matrix/")
# mouseHashtag12 <- JoinLayers(mouseHashtag12)
# exprMat <- mouseHashtag12[["RNA"]]$counts
# cellInfo <- data.frame(seuratCluster=(mouseHashtag12$seurat_clusters))

setwd("./Scenic12/")

loomseurat <- as.loom(mouseHashtag12,filename = "./mousehashtag12.loom")
loomseurat

infile <- H5File$new("./mousehashtag12.loom")

loomPath <- system.file(package="SCENIC", "./mouseHashtag.loom")
loom <- open_loom("./mousehashtag12.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

dim(exprMat)
cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)

cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$seurat_clusters))

saveRDS(cellInfo, file="cellInfo.Rds")

colVars <- list(seurat_clusters=c("0"="red", 
                           "1"="darkorange", 
                           "2"="lightgreen", 
                           "3"="blue", 
                           "4"="cyan",
                           "5"="purple",
                           "6"="pink"))
colVars$seurat_clusters <- colVars$seurat_clusters[intersect(names(colVars$seurat_clusters), cellInfo$seurat_clusters)]
saveRDS(colVars, file="./seurat_cluster")
plot.new(); legend(0,1, fill=colVars$seurat_clusters, legend=names(colVars$seurat_clusters))

# Sys.setenv(LIBARROW_MINIMAL = "false")
# Sys.setenv(ARROW_WITH_ZSTD = "ON")
# install_arrow()
# .rs.restartR()

# feather12 <- read_feather("./cistargetmouse/mm9-tss-centered-10kb-7species.mc9nr.feather")

db <- importRankings("./cistargetmouse/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
names(db@rankings)[1] <- "features"
db@org <- "mgi"
db@genome <- "mm10"
arrow::write_feather(db@rankings, "./cistargetmouse/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")



db <- importRankings("./cistargetmouse/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
names(db@rankings)[1] <- "features"
db@org <- "mgi"
db@genome <- "mm10"
arrow::write_feather(db@rankings,
                     "./cistargetmouse/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

# And then you can run the next step to keep genes.




mm10_dbs <- list('500bp'= 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather', 
                 '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather')
db_mcVersion <- 'mc9nr'
db_path <- './cistargetmouse'

scenicOptions <- initializeScenic(org='mgi', dbs = mm10_dbs, dbDir = db_path, datasetTitle='scenicmouse12', nCores=10) 

motifAnnotations_mgi <- motifAnnotations


# scenicOptions <- initializeScenic(org="mgi", dbDir="./cistargetmouse/",dbs = "cistargetmouse/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather", nCores=8) 

exprMat <- as.matrix(mouseHashtag12[["RNA"]]$counts)
genesKept <- geneFiltering(exprMat, scenicOptions,  minCountsPerGene = 3 * 0.01 * ncol(exprMat),  minSamples = ncol(exprMat) * 0.01) 
exprMat_filtered <- exprMat[genesKept, ] #Gene list saved in int/1.1_genesKept.Rds


org <- "mgi"
dbDir <- "./cistargetmouse" # RcisTarget databases location
dbs <- defaultDbNames[[org]]
myDatasetTitle <- "SCENIC mousehashtag12" # choose a name for your analysis
scenicOptions <- initializeScenic(org="mgi", dbDir="./cistargetmouse/", dbs=dbs, datasetTitle=myDatasetTitle)
motifAnnotations_mgi <- motifAnnotations
scenicOptions <- initializeScenic(org="mgi", dbDir="./cistargetmouse/", dbs=dbs, datasetTitle=myDatasetTitle)
scenicOptions@inputDatasetInfo$cellInfo <- "./cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "./colVars.Rds"
scenicOptions@settings$nCores <- 8

saveRDS(scenicOptions, file="./scenicOptions.Rds") 
#################
# Script Scenic #

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

interestingGenes <- c("Zbtb16","Klrk1","Dntt","Sell","fghjdfg")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

rm(exprMat)

runCorrelation(exprMat_filtered, scenicOptions)

# Là c'est giga long genre des heures
exprMat_filtered <- log2(exprMat_filtered+1)
# runGenie3(exprMat_filtered, scenicOptions)
coormat <-readRDS("./int/1.2_corrMat.Rds")

write.csv(coormat,file = "./coordmat.csv")

infile <- H5File$new("./mousehashtag12.loom")

loom <- open_loom("./mousehashtag12.loom")
exprMat <- get_dgem(loom)
close_loom(loom)
# Optional: log expression (for TF expression plot, it does not affect any other calculation)
exprMat_log <- log2(exprMat+1)
dim(exprMat)

scenicOptions <- readRDS("./scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 1 # 1 car sinon runscenic_2_createRegulons marche pas
scenicOptions@settings$seed <- 123
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod="top5perTarget")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="./scenicOptions.Rds") # To save status
# 
# aucellApp <- plotTsne_AUCellHtml(scenicOptions, exprMat_log,"./int/plotTsne_AUCellHtml")
# savedSelections <- runApp(aucellApp)


# # Save the modified thresholds:
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds") 


# # scenicOptions@settings$devType="png"
# scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)

nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="seurat_clusters", cex=.5)

par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="seurat_clusters", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

# DGEM (Digital gene expression matrix)(non-normalized counts)
exprMat <- get_dgem(loom)
dgem <- exprMat
head(colnames(dgem))  #should contain the Cell ID/name

# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/mousehashtagSCENIC.loom"
export2loom(scenicOptions, exprMat)

scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)

# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")

# Show TF expression:
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Rorc","Jun","E2f1")],], plots="Expression")

dens2d <- bkde2D(tSNE_scenic$Y, 2)$fhat
image(dens2d, col=brewer.pal(9, "YlOrBr"), axes=FALSE)
contour(dens2d, add=TRUE, nlevels=5, drawlabels=FALSE)

par(mfrow=c(1,2))

regulonNames <- c( "Rorc","E2f1")
cellCol <- SCENIC::plotEmb_rgb(scenicOptions, regulonNames, aucType="AUC", aucMaxContrast=0.6)


regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Rorc", "Jun")]


cellInfo <- readRDS("cellInfo.Rds")

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")


setwd("..")


#################################
# Ajout des TCRs dans l'analyse #


setwd("/Users/osourasinh/Documents/Stage M1/Données TCR")
TCR <- read.csv ("filtered_contig_annotations.csv")

#########
# TCR12 #
TCRlist12 <- createHTOContigList(TCR, sc = mouseHashtag12, group.by = "ident")

head(TCRlist12[[1]])

combined.TCR12 <- combineTCR(TCRlist12, 
                           samples = c("Mait0","CyclingS","Mait17a","Mait1","CyclingG2M","Mait17b","Precursor"
),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = FALSE)



# write.csv(combined.TCR12$Mait0,"./combinedTCR_Hashtag12_Mait0.csv")
# write.csv(combined.TCR12$CyclingS,"./combinedTCR_Hashtag12_CyclingS.csv")
# write.csv(combined.TCR12$Mait17a,"./combinedTCR_Hashtag12_Mait17a.csv")
# write.csv(combined.TCR12$Mait1,"./combinedTCR_Hashtag12_Mait1.csv")
# write.csv(combined.TCR12$CyclingG2M,"./combinedTCR_Hashtag12_CyclingG2M.csv")
# write.csv(combined.TCR12$Mait17b,"./combinedTCR_Hashtag12_Mait17b.csv")
# write.csv(combined.TCR12$Precursor,"./combinedTCR_Hashtag12_Precursor.csv")



exportClones(combined.TCR12, 
             write.file = TRUE,
             dir = "./",
             file.name = "mousehashtag12clones.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR12, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag12/expressiongeneTRAVhashtag12.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12, 
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",  
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag12/expressiongeneTRAJhashtag12.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag12/expressiongeneTRBVhashtag12.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag12/expressiongeneTRBDhashtag12.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag12/expressiongeneTRBJhashtag12.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR12[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR12, 
          chain = "TRA", 
          aa.length = 13)

# png(file="./mouseHashtag12/percentaaTRAhashtag12.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR12,
#           chain = "TRA",
#           aa.length = 13)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR12, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="./mouseHashtag12/positionalEntropyTRAhashtag12.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR12, 
#                   chain = "TRA", 
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR12, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="./mouseHashtag12/percentTRAVhashtag12.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12, 
#              chain = "TRA", 
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR12, 
             chain = "TRA", 
             gene = "Jgene")

# png(file="./mouseHashtag12/percentTRAJhashtag12.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12, 
#              chain = "TRA", 
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR12, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="./mouseHashtag12/percentTRBVhashtag12.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12, 
#              chain = "TRB", 
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR12, 
             chain = "TRB", 
             gene = "Jgene")

# png(file="./mouseHashtag12/percentTRBJhashtag12.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12, 
#              chain = "TRB", 
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR12, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="./mouseHashtag12/percentTRBDhashtag12.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12, 
#              chain = "TRB", 
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR12, 
                         chain = "TRA", 
                         gene = "Vgene", 
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="./mouseHashtag12/PCA-TRAVhashtag12.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR12,
         chain = "TRA")


df.genes <- percentVJ(combined.TCR12, 
                      chain = "TRA", 
                      exportTable = TRUE)

# png(file="./mouseHashtag12/PCA-TRAVhashtag12.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR12, 
#           chain = "TRA")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic()

percentKmer(combined.TCR12, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag12/percentKmeraaTRAhashtag12.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR12, 
#             cloneCall = "aa",
#             chain = "TRA", 
#             motif.length = 3, 
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR12, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag12/percentKmerntTRAhashtag12.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR12,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR12, 
                cloneCall = "gene", 
                n.boots = 20)

clonalRarefaction(combined.TCR12,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR12,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR12,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR12,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)


clonalSizeDistribution(combined.TCR12, 
                       cloneCall = "aa", 
                       method= "ward.D2")


clonalOverlap(combined.TCR12, 
              cloneCall = "strict", 
              method = "morisita")

# png(file="./mouseHashtag12/clonalOverlaphashtag12.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR12, 
#               cloneCall = "strict", 
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR12, 
              cloneCall = "strict", 
              method = "raw")



sce <- combineExpression(combined.TCR12, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

??plotUMAP

clonalNetwork(TCRlist12, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")

#########
# TCR1 #
mouseHashtag1split = SplitObject(mouseHashtag1, split.by = "old.ident")

# Isolation du Hashtag 1
mouseHashtag1 <- mouseHashtag1split$`Hashtag-1`

TCRlist1 <- createHTOContigList(TCR, sc = mouseHashtag1, group.by = "ident")

head(TCRlist1[[1]])

combined.TCR1 <- combineTCR(TCRlist1, 
                             samples = c("Mait0","CyclingS","Mait17a","Mait1","CyclingG2M","Mait17b","Precursor"
                             ),
                             removeNA = FALSE, 
                             removeMulti = FALSE, 
                             filterMulti = FALSE)




# write.csv(combined.TCR1$Mait0,"./mouseHashtag1/combinedTCR_Hashtag1_Mait0.csv")
# write.csv(combined.TCR1$CyclingS,"./mouseHashtag1/combinedTCR_Hashtag1_CyclingS.csv")
# write.csv(combined.TCR1$Mait17a,"./mouseHashtag1/combinedTCR_Hashtag1_Mait17a.csv")
# write.csv(combined.TCR1$Mait1,"./mouseHashtag1/combinedTCR_Hashtag1_Mait1.csv")
# write.csv(combined.TCR1$CyclingG2M,"./mouseHashtag1/combinedTCR_Hashtag1_CyclingG2M.csv")
# write.csv(combined.TCR1$Mait17b,"./mouseHashtag1/combinedTCR_Hashtag1_Mait17b.csv")
# write.csv(combined.TCR1$Precursor,"./mouseHashtag1/combinedTCR_Hashtag1_Precursor.csv")



exportClones(combined.TCR1, 
             write.file = TRUE,
             dir = "./",
             file.name = "mouseHashtag1clones.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR1, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag1/expressiongeneTRAVHashtag1.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1,
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag1/expressiongeneTRAJHashtag1.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag1/expressiongeneTRBVHashtag1.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag1/expressiongeneTRBDHashtag1.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag1/expressiongeneTRBJHashtag1.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR1[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR1, 
          chain = "TRA", 
          aa.length = 15)

# png(file="./mouseHashtag1/percentaaTRAHashtag1.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR1,
#           chain = "TRA",
#           aa.length = 15)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR1, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="./mouseHashtag1/positionalEntropyTRAHashtag1.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR1,
#                   chain = "TRA",
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR1, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="./mouseHashtag1/percentTRAVHashtag1.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1,
#              chain = "TRA",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR1, 
             chain = "TRA", 
             gene = "Jgene")

# png(file="./mouseHashtag1/percentTRAJHashtag1.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1,
#              chain = "TRA",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR1, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="./mouseHashtag1/percentTRBVHashtag1.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1,
#              chain = "TRB",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR1, 
             chain = "TRB", 
             gene = "Jgene")

# png(file="./mouseHashtag1/percentTRBJHashtag1.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1,
#              chain = "TRB",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR1, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="./mouseHashtag1/percentTRBDHashtag1.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1,
#              chain = "TRB",
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR1, 
                         chain = "TRA", 
                         gene = "Vgene", 
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="./mouseHashtag1/PCA-TRAVHashtag1.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR1,
          chain = "TRA")


df.genes <- percentVJ(combined.TCR1, 
                      chain = "TRA", 
                      exportTable = TRUE)

# png(file="./mouseHashtag1/percentVJ-TRAHashtag1.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR1,
#           chain = "TRA")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic()

percentKmer(combined.TCR1, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag1/percentKmeraaTRAHashtag1.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR1,
#             cloneCall = "aa",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR1, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag1/percentKmerntTRAHashtag1.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR1,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR1, 
                cloneCall = "gene", 
                n.boots = 20)

clonalRarefaction(combined.TCR1,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR1,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR1,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR1,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)


clonalSizeDistribution(combined.TCR1, 
                       cloneCall = "aa", 
                       method= "ward.D2")


clonalOverlap(combined.TCR1, 
              cloneCall = "strict", 
              method = "morisita")

# png(file="./mouseHashtag1/clonalOverlapHashtag1.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR1,
#               cloneCall = "strict",
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR1, 
              cloneCall = "strict", 
              method = "raw")



sce <- combineExpression(combined.TCR1, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

??plotUMAP

clonalNetwork(TCRlist1, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")



#########
# TCR2 #
mouseHashtag2split = SplitObject(mouseHashtag2, split.by = "old.ident")

# Isolation du Hashtag 1
mouseHashtag2 <- mouseHashtag2split$`Hashtag-2`
TCRlist2 <- createHTOContigList(TCR, sc = mouseHashtag2, group.by = "ident")

head(TCRlist2[[1]])

combined.TCR2 <- combineTCR(TCRlist2, 
                            samples = c("Mait1","Mait17b","CyclingG2M","CyclingS","Mait17a","Mait0","Precursor"),
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = FALSE)




# write.csv(combined.TCR2$Mait0,"./mouseHashtag2/combinedTCR_Hashtag2_Mait0.csv")
# write.csv(combined.TCR2$CyclingS,"./mouseHashtag2/combinedTCR_Hashtag2_CyclingS.csv")
# write.csv(combined.TCR2$Mait17a,"./mouseHashtag2/combinedTCR_Hashtag2_Mait17a.csv")
# write.csv(combined.TCR2$Mait1,"./mouseHashtag2/combinedTCR_Hashtag2_Mait1.csv")
# write.csv(combined.TCR2$CyclingG2M,"./mouseHashtag2/combinedTCR_Hashtag2_CyclingG2M.csv")
# write.csv(combined.TCR2$Mait17b,"./mouseHashtag2/combinedTCR_Hashtag2_Mait17b.csv")
# write.csv(combined.TCR2$Precursor,"./mouseHashtag2/combinedTCR_Hashtag2_Precursor.csv")



exportClones(combined.TCR2, 
             write.file = TRUE,
             dir = "./",
             file.name = "mouseHashtag2clones.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR2, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag2/expressiongeneTRAVHashtag2.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2,
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag2/expressiongeneTRAJHashtag2.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag2/expressiongeneTRBVHashtag2.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag2/expressiongeneTRBDHashtag2.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="./mouseHashtag2/expressiongeneTRBJHashtag2.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR2[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR2, 
          chain = "TRA", 
          aa.length = 13)

# png(file="./mouseHashtag2/percentaaTRAHashtag2.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR2,
#           chain = "TRA",
#           aa.length = 13)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR2, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="./mouseHashtag2/positionalEntropyTRAHashtag2.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR2,
#                   chain = "TRA",
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR2, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="./mouseHashtag2/percentTRAVHashtag2.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2,
#              chain = "TRA",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR2, 
             chain = "TRA", 
             gene = "Jgene")

# png(file="./mouseHashtag2/percentTRAJHashtag2.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2,
#              chain = "TRA",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR2, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="./mouseHashtag2/percentTRBVHashtag2.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2,
#              chain = "TRB",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR2, 
             chain = "TRB", 
             gene = "Jgene")

# png(file="./mouseHashtag2/percentTRBJHashtag2.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2,
#              chain = "TRB",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR2, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="./mouseHashtag2/percentTRBDHashtag2.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2,
#              chain = "TRB",
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR2, 
                         chain = "TRA", 
                         gene = "Vgene", 
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="./mouseHashtag2/PCA-TRAVHashtag2.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR2,
          chain = "TRA")


df.genes <- percentVJ(combined.TCR2, 
                      chain = "TRA", 
                      exportTable = TRUE)

# png(file="./mouseHashtag2/percentVJ-TRAVHashtag2.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR2,
#           chain = "TRA")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic()

percentKmer(combined.TCR2, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag2/percentKmeraaTRAHashtag2.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR2,
#             cloneCall = "aa",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR2, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="./mouseHashtag2/percentKmerntTRAHashtag2.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR2,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR2, 
                cloneCall = "gene", 
                n.boots = 20)

clonalRarefaction(combined.TCR2,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR2,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR2,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

clonalRarefaction(combined.TCR2,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)


clonalSizeDistribution(combined.TCR2, 
                       cloneCall = "aa", 
                       method= "ward.D2")


clonalOverlap(combined.TCR2, 
              cloneCall = "strict", 
              method = "morisita")

# png(file="./mouseHashtag2/clonalOverlapHashtag2.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR2,
#               cloneCall = "strict",
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR2, 
              cloneCall = "strict", 
              method = "raw")



sce <- combineExpression(combined.TCR2, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

??plotUMAP

clonalNetwork(TCRlist2, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")




##############
# TCR1 clean #
mouseHashtag1cleansplit = SplitObject(mouseHashtag12cleansplit1, split.by = "old.ident")

# Isolation du Hashtag 1
mouseHashtag1clean <- mouseHashtag1cleansplit$`Hashtag-1`

TCRlist1clean <- createHTOContigList(TCR, sc = mouseHashtag1clean, group.by = "ident")

head(TCRlist1clean[[1]])

combined.TCR1clean <- combineTCR(TCRlist1clean, 
                            samples = c("Mait0","CyclingS","Mait17a","Mait1","CyclingG2M","Mait17b"
                            ),
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = FALSE)
setwd("./MouseHashtag1clean/")



# write.csv(combined.TCR1clean$Mait0,"combinedTCR_Hashtag1_Mait0clean.csv")
# write.csv(combined.TCR1clean$CyclingS,"combinedTCR_Hashtag1_CyclingSclean.csv")
# write.csv(combined.TCR1clean$Mait17a,"combinedTCR_Hashtag1_Mait17aclean.csv")
# write.csv(combined.TCR1clean$Mait1,"combinedTCR_Hashtag1_Mait1clean.csv")
# write.csv(combined.TCR1clean$CyclingG2M,"combinedTCR_Hashtag1_CyclingG2Mclean.csv")
# write.csv(combined.TCR1clean$Mait17b,"combinedTCR_Hashtag1_Mait17bclean.csv")



exportClones(combined.TCR1clean, 
             write.file = TRUE,
             dir = "./",
             file.name = "mouseHashtag1clonesclean.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR1clean, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAVHashtag1clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1clean,
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1clean, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAJHashtag1clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1clean,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1clean, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBVHashtag1clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1clean,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1clean, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBDHashtag1clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1clean,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1clean, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBJHashtag1clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR1clean,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR1clean, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR1clean[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR1clean, 
          chain = "TRA", 
          aa.length = 15)

# png(file="percentaaTRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR1clean,
#           chain = "TRA",
#           aa.length = 15)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR1clean, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="positionalEntropyTRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR1clean,
#                   chain = "TRA",
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR1clean, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="percentTRAVHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1clean,
#              chain = "TRA",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR1clean, 
             chain = "TRA", 
             gene = "Jgene")
 
# png(file="percentTRAJHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1clean,
#              chain = "TRA",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR1clean, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="percentTRBVHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1clean,
#              chain = "TRB",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR1clean, 
             chain = "TRB", 
             gene = "Jgene")

# png(file="percentTRBJHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1clean,
#              chain = "TRB",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR1clean, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="percentTRBDHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR1clean,
#              chain = "TRB",
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR1clean,
                         chain = "TRA",
                         gene = "Vgene",
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="PCA-TRAVHashtag1clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR1clean,
          chain = "TRA")


df.genes <- percentVJ(combined.TCR1clean, 
                      chain = "TRA", 
                      exportTable = TRUE)

# png(file="percentVJ-TRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR1clean,
#           chain = "TRA")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic()

# png(file="PCA-TRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()
 
percentKmer(combined.TCR1clean, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmeraaTRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR1clean,
#             cloneCall = "aa",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR1clean, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmerntTRAHashtag1clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR1clean,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR1clean, 
                cloneCall = "gene", 
                n.boots = 20)

# png(file="clonal Diversity Hashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalDiversity(combined.TCR1clean, 
#                 cloneCall = "gene", 
#                 n.boots = 20)
# dev.off()

clonalRarefaction(combined.TCR1clean,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction1 Hashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR1clean,
#                   plot.type = 1,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR1clean,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction2 Hashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR1clean,
#                   plot.type = 2,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR1clean,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction3 Hashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR1clean,
#                   plot.type = 3,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR1clean,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)

# png(file="clonalRarefaction1-1 Hashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR1clean,
#                   plot.type = 1,
#                   hill.numbers = 1,
#                   n.boots = 2)
# dev.off()

clonalSizeDistribution(combined.TCR1clean, 
                       cloneCall = "aa", 
                       method= "ward.D2")


clonalOverlap(combined.TCR1clean, 
              cloneCall = "strict", 
              method = "morisita")

# png(file="clonalOverlapHashtag1clean.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR1clean,
#               cloneCall = "strict",
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR1clean, 
              cloneCall = "strict", 
              method = "raw")



sce <- combineExpression(combined.TCR1clean, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))


clonalNetwork(TCRlist1clean, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")



##############
# TCR2 clean #
mouseHashtag2cleansplit = SplitObject(mouseHashtag12cleansplit2, split.by = "old.ident")

# Isolation du Hashtag 2
mouseHashtag2clean <- mouseHashtag2cleansplit$`Hashtag-2`

TCRlist2clean <- createHTOContigList(TCR, sc = mouseHashtag2clean, group.by = "ident")

head(TCRlist2clean[[1]])

combined.TCR2clean <- combineTCR(TCRlist2clean, 
                                 samples = c("Mait0","CyclingS","Mait17a","Mait1","CyclingG2M","Mait17b"
                                 ),
                                 removeNA = FALSE, 
                                 removeMulti = FALSE, 
                                 filterMulti = FALSE)

setwd("../MouseHashtag2clean/")



# write.csv(combined.TCR2clean$Mait0,"combinedTCR_Hashtag2_Mait0clean.csv")
# write.csv(combined.TCR2clean$CyclingS,"combinedTCR_Hashtag2_CyclingSclean.csv")
# write.csv(combined.TCR2clean$Mait17a,"combinedTCR_Hashtag2_Mait17aclean.csv")
# write.csv(combined.TCR2clean$Mait1,"combinedTCR_Hashtag2_Mait1clean.csv")
# write.csv(combined.TCR2clean$CyclingG2M,"combinedTCR_Hashtag2_CyclingG2Mclean.csv")
# write.csv(combined.TCR2clean$Mait17b,"combinedTCR_Hashtag2_Mait17bclean.csv")



exportClones(combined.TCR2clean, 
             write.file = TRUE,
             dir = "./",
             file.name = "mouseHashtag2clonesclean.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR2clean, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAVHashtag2clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2clean,
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2clean, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAJHashtag2clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2clean,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2clean, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBVHashtag2clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2clean,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2clean, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBDHashtag2clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2clean,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2clean, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBJHashtag2clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR2clean,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR2clean, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR2clean[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR2clean, 
          chain = "TRA", 
          aa.length = 15)

# png(file="percentaaTRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR2clean,
#           chain = "TRA",
#           aa.length = 15)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR2clean, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="positionalEntropyTRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR2clean,
#                   chain = "TRA",
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR2clean, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="percentTRAVHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2clean,
#              chain = "TRA",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR2clean, 
             chain = "TRA", 
             gene = "Jgene")

# png(file="percentTRAJHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2clean,
#              chain = "TRA",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR2clean, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="percentTRBVHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2clean,
#              chain = "TRB",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR2clean, 
             chain = "TRB", 
             gene = "Jgene")

png(file="percentTRBJHashtag2clean.png", res = 110, height = 1100, width = 1000)
percentGenes(combined.TCR2clean,
             chain = "TRB",
             gene = "Jgene")
dev.off()

percentGenes(combined.TCR2clean, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="percentTRBDHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR2clean,
#              chain = "TRB",
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR2clean,
                         chain = "TRA",
                         gene = "Vgene",
                         exportTable = TRUE)

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="PCA-TRAVHashtag2clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR2clean,
          chain = "TRA")


df.genes <- percentVJ(combined.TCR2clean, 
                      chain = "TRA", 
                      exportTable = TRUE)

# png(file="percentVJ-TRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR2clean,
#           chain = "TRA")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
  theme_classic()

# png(file="PCA-TRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentKmer(combined.TCR2clean, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmeraaTRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR2clean,
#             cloneCall = "aa",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR2clean, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmerntTRAHashtag2clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR2clean,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR2clean, 
                cloneCall = "gene", 
                n.boots = 20)

# png(file="clonal Diversity Hashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalDiversity(combined.TCR2clean, 
#                 cloneCall = "gene", 
#                 n.boots = 20)
# dev.off()

clonalRarefaction(combined.TCR2clean,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction1 Hashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 1,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR2clean,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction2 Hashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 2,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR2clean,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction3 Hashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 3,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR2clean,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)
 
# png(file="clonalRarefaction1-1 Hashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 1,
#                   hill.numbers = 1,
#                   n.boots = 2)
# dev.off()

clonalSizeDistribution(combined.TCR2clean, 
                       cloneCall = "aa", 
                       method= "ward.D2")


clonalOverlap(combined.TCR2clean, 
              cloneCall = "strict", 
              method = "morisita")

# png(file="clonalOverlapHashtag2clean.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR2clean,
#               cloneCall = "strict",
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR2clean, 
              cloneCall = "strict", 
              method = "raw")



sce <- combineExpression(combined.TCR2clean, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))


clonalNetwork(TCRlist2clean, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")
###############
# TCR12 clean #
TCRlist12clean <- createHTOContigList(TCR, sc = mouseHashtag12clean, group.by = "ident")

head(TCRlist12clean[[1]])

combined.TCR12clean <- combineTCR(TCRlist12clean, 
                             samples = c("Mait0_Hashtag-1","Mait17a_Hashtag-1"," Mait1_Hashtag-1"," CyclingS_Hashtag-1"," CyclingG2M_Hashtag-1"," Mait17b_Hashtag-1","Mait0_Hashtag-2","Mait17a_Hashtag-2"," Mait1_Hashtag-2"," CyclingS_Hashtag-2"," CyclingG2M_Hashtag-2"," Mait17b_Hashtag-2"),
                             removeNA = FALSE, 
                             removeMulti = FALSE, 
                             filterMulti = FALSE)

setwd("../mouseHashtag12clean/")

# write.csv(combined.TCR12clean$Mait0,"./combinedTCR_Hashtag12clean_Mait0.csv")
# write.csv(combined.TCR12clean$CyclingS,"./combinedTCR_Hashtag12clean_CyclingS.csv")
# write.csv(combined.TCR12clean$Mait17a,"./combinedTCR_Hashtag12clean_Mait17a.csv")
# write.csv(combined.TCR12clean$Mait1,"./combinedTCR_Hashtag12clean_Mait1.csv")
# write.csv(combined.TCR12clean$CyclingG2M,"./combinedTCR_Hashtag12clean_CyclingG2M.csv")
# write.csv(combined.TCR12clean$Mait17b,"./combinedTCR_Hashtag12clean_Mait17b.csv")
# write.csv(combined.TCR12clean$Precursor,"./combinedTCR_Hashtag12clean_Precursor.csv")



exportClones(combined.TCR12clean, 
             write.file = TRUE,
             dir = "./",
             file.name = "mouseHashtag12cleanclones.csv")

# .rs.restartR()
#ça c'est bien
vizGenes(combined.TCR12clean, 
         x.axis = "TRAV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAVhashtag12clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12clean,
#          x.axis = "TRAV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12clean, 
         x.axis = "TRAJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRAJhashtag12clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12clean,
#          x.axis = "TRAJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12clean, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBVhashtag12clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12clean,
#          x.axis = "TRBV",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12clean, 
         x.axis = "TRBD",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBDhashtag12clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12clean,
#          x.axis = "TRBD",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12clean, 
         x.axis = "TRBJ",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE)

# png(file="expressiongeneTRBJhashtag12clean.png", res = 110, height = 1100, width = 1000)
# vizGenes(combined.TCR12clean,
#          x.axis = "TRBJ",
#          y.axis = NULL,
#          plot = "barplot",
#          scale = TRUE)
# dev.off()

vizGenes(combined.TCR12clean, 
         x.axis = "TRAV",
         y.axis = "TRAJ",
         plot = "heatmap",  
         scale = TRUE)

vizGenes(combined.TCR12clean[c(4,5,6)], 
         x.axis = "TRBV",
         y.axis = "TRBJ",
         plot = "heatmap",  
         scale = TRUE)

#ça c'est bien aussi
percentAA(combined.TCR12clean, 
          chain = "TRA", 
          aa.length = 13)

# png(file="percentaaTRAhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentAA(combined.TCR12clean,
#           chain = "TRA",
#           aa.length = 13)
# dev.off()


#ça c'est bien aussi
positionalEntropy(combined.TCR12clean, 
                  chain = "TRA", 
                  aa.length = 18)

# png(file="positionalEntropyTRAhashtag12clean.png", res = 110, height = 1100, width = 1000)
# positionalEntropy(combined.TCR12clean,
#                   chain = "TRA",
#                   aa.length = 18)
# dev.off()

#ça c'est bien aussi
percentGenes(combined.TCR12clean, 
             chain = "TRA", 
             gene = "Vgene")

# png(file="percentTRAVhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12clean,
#              chain = "TRA",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR12clean, 
             chain = "TRA", 
             gene = "Jgene")

# png(file="percentTRAJhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12clean,
#              chain = "TRA",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR12clean, 
             chain = "TRB", 
             gene = "Vgene")

# png(file="percentTRBVhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12clean,
#              chain = "TRB",
#              gene = "Vgene")
# dev.off()

percentGenes(combined.TCR12clean, 
             chain = "TRB", 
             gene = "Jgene")

# png(file="percentTRBJhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12clean,
#              chain = "TRB",
#              gene = "Jgene")
# dev.off()

percentGenes(combined.TCR12clean, 
             chain = "TRB", 
             gene = "Dgene")

# png(file="percentTRBDhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentGenes(combined.TCR12clean,
#              chain = "TRB",
#              gene = "Dgene")
# dev.off()

df.genes <- percentGenes(combined.TCR12clean, 
                         chain = "TRB", 
                         gene = "Vgene", 
                         exportTable = TRUE)


#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes)))
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) +
  scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
  theme_classic()

# png(file="PCA-TRBVhashtag12clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) +
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) +
#   guides(fill=guide_legend(title="Samples")) +
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) +
#   theme_classic()
# dev.off()

percentVJ(combined.TCR12clean,
          chain = "TRB")


df.genes <- percentAA(combined.TCR12clean, 
                      chain = "TRB", 
                      exportTable = TRUE)

# png(file="PCA-TRBhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentVJ(combined.TCR12clean,
#           chain = "TRB")
# dev.off()

#Performing PCA
pc <- prcomp(df.genes)

#Getting data frame to plot from
df <- as.data.frame(cbind(pc$x[,1:2], rownames(df.genes))) 
df$PC1 <- as.numeric(df$PC1)
df$PC2 <- as.numeric(df$PC2)

#Plotting
ggplot(df, aes(x = PC1, y = PC2)) + 
  geom_jitter(aes(fill =df[,3]), shape = 21, size = 5) + 
  guides(fill=guide_legend(title="Samples")) + 
  scale_fill_manual(values = hcl.colors(nrow(df), "viridis")) + 
  theme_classic()

# png(file="PCA-TRAhashtag12clean.png", res = 110, height = 1100, width = 1000)
# ggplot(df, aes(x = PC1, y = PC2)) + 
#   geom_point(aes(fill =df[,3]), shape = 21, size = 5) + 
#   guides(fill=guide_legend(title="Samples")) + 
#   scale_fill_manual(values = hcl.colors(nrow(df), "inferno")) + 
#   theme_classic()
# dev.off()

percentKmer(combined.TCR12clean, 
            cloneCall = "aa",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmeraaTRAhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR12clean,
#             cloneCall = "aa",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

percentKmer(combined.TCR12clean, 
            cloneCall = "nt",
            chain = "TRA", 
            motif.length = 3, 
            top.motifs = 25)

# png(file="percentKmerntTRAhashtag12clean.png", res = 110, height = 1100, width = 1000)
# percentKmer(combined.TCR12clean,
#             cloneCall = "nt",
#             chain = "TRA",
#             motif.length = 3,
#             top.motifs = 25)
# dev.off()

clonalDiversity(combined.TCR12clean, 
                cloneCall = "gene", 
                n.boots = 20)

# png(file="clonal Diversity Hashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalDiversity(combined.TCR12clean,
#                 cloneCall = "gene",
#                 n.boots = 20)
# dev.off()

clonalRarefaction(combined.TCR12clean,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction1 Hashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 1,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()


clonalRarefaction(combined.TCR12clean,
                  plot.type = 2,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction2 Hashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 2,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR12clean,
                  plot.type = 3,
                  hill.numbers = 0,
                  n.boots = 2)

# png(file="clonalRarefaction3 Hashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 3,
#                   hill.numbers = 0,
#                   n.boots = 2)
# dev.off()

clonalRarefaction(combined.TCR12clean,
                  plot.type = 1,
                  hill.numbers = 1,
                  n.boots = 2)


# png(file="clonalRarefaction1-1 Hashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalRarefaction(combined.TCR2clean,
#                   plot.type = 1,
#                   hill.numbers = 1,
#                   n.boots = 2)
# dev.off()


clonalSizeDistribution(combined.TCR12clean, 
                       cloneCall = "strict", 
                       method= "ward.D2")

?clonalSizeDistribution

clonalOverlap(combined.TCR12clean, 
              cloneCall = "strict", 
              method = "morisita")




# png(file="clonalOverlaphashtag12clean.png", res = 110, height = 1100, width = 1000)
# clonalOverlap(combined.TCR12clean,
#               cloneCall = "strict",
#               method = "morisita")
# dev.off()

clonalOverlap(combined.TCR12clean, 
              cloneCall = "strict", 
              method = "raw")


sce <- combineExpression(combined.TCR12clean, 
                         sce, 
                         cloneCall="gene", 
                         group.by = "sample", 
                         proportion = TRUE)

#Define color palette 
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

plotUMAP(sce, colour_by = "cloneSize") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,5,7)]))

??plotUMAP

clonalNetwork(TCRlist12clean, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              filter.clones = NULL,
              filter.identity = NULL,
              cloneCall = "aa")


##################
# May be usefull #



# mouse[["percent_mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")
# 
# VlnPlot(mouseHashtag[["Hashtag-1"]], features = c("nCount_RNA"))
# VlnPlot(mouseHashtag[["Hashtag-2"]], features = c("nCount_RNA"))
# 
# VlnPlot(mouse, features =c("nFeature_RNA", "nCount_RNA","percent_mt"), alpha = 0.2)
# VlnPlot(mouse, features =c("nFeature_HTO", "nCount_HTO","percent_mt"), alpha = 0.2)
# VlnPlot(mouse, features =c("nFeature_Protein", "nCount_Protein","percent_mt"), alpha = 0.2)
# mouse <- subset(mouse, subset = percent_mt < 10 & percent_C1qc < 0.01 )
# VlnPlot(mouse, features =c("nFeature_RNA", "nCount_RNA","percent_mt"), alpha = 0.2)
# VlnPlot(mouse, features =c("nFeature_HTO", "nCount_HTO","percent_mt"), alpha = 0.2)
# VlnPlot(mouse, features =c("nFeature_Protein", "nCount_Protein","percent_mt"), alpha = 0.2)
# 
# 
# DimHeatmap(mouseHashtag45, dims = 1, cells = 500, balanced = TRUE)

# cluster0.markers <- FindMarkers(mouseHashtag12, ident.1 = "0",only.pos = TRUE)
# row.names(cluster0.markers)[1]
# cluster1.markers <- FindMarkers(mouseHashtag12, ident.1 = "1",only.pos = TRUE)
# row.names(cluster1.markers)[1]
# cluster2.markers <- FindMarkers(mouseHashtag12, ident.1 = "2",only.pos = TRUE)
# row.names(cluster2.markers)[1]
# cluster3.markers <- FindMarkers(mouseHashtag12, ident.1 = "3",only.pos = TRUE)
# row.names(cluster3.markers)[1]
# cluster4.markers <- FindMarkers(mouseHashtag12, ident.1 = "4",only.pos = TRUE)
# row.names(cluster4.markers)[1]
# cluster5.markers <- FindMarkers(mouseHashtag12, ident.1 = "5",only.pos = TRUE)
# row.names(cluster5.markers)[1]
# cluster3.markers <- FindMarkers(mouseHashtag12, ident.1 = 3, only.pos = TRUE)
# head(cluster3.markers, n = 5)
# cluster5.markers <- FindMarkers(mouseHashtag12, ident.1 = 5, ident.2 = c(0, 3))
# head(cluster5.markers, n = 5)
# cluster6.markers <- FindMarkers(mouseHashtag12, ident.1 = 6)
# head(cluster6.markers, n = 5)
#ident.1 choix du cluster ident.2 choisir les clusters qui doivent être différents only.pos donne les gènes qui sont exprimés

# Voir pour faire un programme qui le fait tout seul grâce aux gènes marqueurs des clusters.
# Naming_cluster<- function(objectseurat){
#   old_cluster_names = c("0","1","2","3","4","5")
#   cluster0.markers <- FindMarkers(objectseurat, ident.1 = "0",only.pos = TRUE)
#   Markers0 <- row.names(cluster0.markers)[1]
#   if (Markers0 == MAIT1signature[1]){
#     old_cluster_names[1]="Mait1"}
#   else if (Markers0 == MAIT17signature[1]){
#     old_cluster_names[1]="Mait17"}
#   else if (Markers0 == MAIT0signature[1]){
#     old_cluster_names[1]="Interm"}
#   else if (Markers0 == "Dock2"){
#     old_cluster_names[1]="Immature"}
#   else if (Markers0 == "Pclaf"){
#     old_cluster_names[1]="Cycling G2M"}
#   else if (Markers0 == "Hist1h2af"){
#     old_cluster_names[1]="Cycling S"}
#   
#   cluster1.markers <- FindMarkers(objectseurat, ident.1 = "1",only.pos = TRUE)
#   Markers1 <- row.names(cluster1.markers)[1]
#   if (Markers1 == MAIT1signature[1]){
#     old_cluster_names[2]="Mait1"}
#   else if (Markers1 == MAIT17signature[1]){
#     old_cluster_names[2]="Mait17"}
#   else if (Markers1 == MAIT1signature[1]){
#     old_cluster_names[2]="Interm"}
#   else if (Markers1 == "Dock2"){
#     old_cluster_names[2]="Immature"}
#   else if (Markers1 == "Pclaf"){
#     old_cluster_names[2]="Cycling G2M"}
#   else if (Markers1 == "Hist1h2af"){
#     old_cluster_names[2]="Cycling S"}
#   
#   cluster2.markers <- FindMarkers(objectseurat, ident.1 = "2",only.pos = TRUE)
#   Markers2 <- row.names(cluster2.markers)[1]
#   if (Markers2 == MAIT1signature[1]){
#     old_cluster_names[3]="Mait1"}
#   else if (Markers2 == MAIT17signature[1]){
#     old_cluster_names[3]="Mait17"}
#   else if (Markers2 == MAIT0signature[1]){
#     old_cluster_names[3]="Interm"}
#   else if (Markers2 == "Dock2"){
#     old_cluster_names[3]="Immature"}
#   else if (Markers2 == "Pclaf"){
#     old_cluster_names[3]="Cycling G2M"}
#   else if (Markers2 == "Hist1h2af"){
#     old_cluster_names[3]="Cycling S"}
#   
#   cluster3.markers <- FindMarkers(objectseurat, ident.1 = "3",only.pos = TRUE)
#   Markers3 <- row.names(cluster3.markers)[1]
#   if (Markers3 == MAIT1signature[1]){
#     old_cluster_names[4]="Mait1"}
#   else if (Markers3 == MAIT17signature[1]){
#     old_cluster_names[4]="Mait17"}
#   else if (Markers3 == MAIT0signature[1]){
#     old_cluster_names[4]="Interm"}
#   else if (Markers3 == "Dock2"){
#     old_cluster_names[4]="Immature"}
#   else if (Markers3 == "Pclaf"){
#     old_cluster_names[4]="Cycling G2M"}
#   else if (Markers3 == "Hist1h2af"){
#     old_cluster_names[4]="Cycling S"}
# 
#   cluster4.markers <- FindMarkers(objectseurat, ident.1 = "4",only.pos = TRUE)
#   Markers4 <- row.names(cluster4.markers)[1]
#   if (Markers4 == MAIT1signature[1]){
#     old_cluster_names[5]="Mait1"}
#   else if (Markers4 == MAIT17signature[1]){
#     old_cluster_names[5]="Mait17"}
#   else if (Markers4 == MAIT0signature[1]){
#     old_cluster_names[5]="Interm"}
#   else if (Markers4 == "Dock2"){
#     old_cluster_names[5]="Immature"}
#   else if (Markers4 == "Pclaf"){
#     old_cluster_names[5]="Cycling G2M"}
#   else if (Markers4 == "Hist1h2af"){
#     old_cluster_names[5]="Cycling S"}
#   
#   cluster5.markers <- FindMarkers(objectseurat, ident.1 = "5",only.pos = TRUE)
#   Markers5 <- row.names(cluster5.markers)[1]
#   if (Markers5 == MAIT1signature[1]){
#     old_cluster_names[6]="Mait1"}
#   else if (Markers5 == MAIT17signature[1]){
#     old_cluster_names[6]="Mait17"}
#   else if (Markers5 == MAIT0signature[1]){
#     old_cluster_names[6]="Interm"}
#   else if (Markers5 == "Dock2"){
#     old_cluster_names[6]="Immature"}
#   else if (Markers5 == "Pclaf"){
#     old_cluster_names[6]="Cycling G2M"}
#   else if (Markers5 == "Hist1h2af"){
#     old_cluster_names[6]="Cycling S"}
#   return(old_cluster_names)
#   }

# new.cluster.ids<- Naming_cluster(mouseHashtag12)
# new.cluster.ids

# head(mouseHashtag12[[]])

# RidgePlot(mouseHashtag12, features = c("Dntt","Ccr9","Rag1","Zbtb16","Klf2","S1pr1","Rps27","Rpl37","Rps24", "Rps7", "Rpl18a","Ptma"))

