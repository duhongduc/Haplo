############################################################################################    
# R codes for phylogeographic distribution of mtDNA haplogroups in Southeast Asia

#	A workflow with RStudio: phylogenetic analyses and visualizations using mtDNA sequences
#		
#	Duc Du
#	August 2023

# Device configration
# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 18.04.6 LTS
###########################################################################33

rm(list = ls(all.names = TRUE))

# Install and load packages

InstallPackages = FALSE

if (InstallPackages) {
  if (!requireNamespace("BiocManager", quietly=TRUE)) 
    install.packages("BiocManager")
  BiocManager::install("msa")
  BiocManager::install("ggtreeExtra")
  
  install.packages("adegenet")
  install.packages("ape")
  install.packages("ggtree")
  install.packages("ggplot2")
  install.packages("ips")
  install.packages("bios2mds")
  install.packages("haplotypes")
  install.packages("pegas")
  install.packages("phytools")
  install.packages("stats")
  install.packages("treeio")
}

# install.packages("strap", dependencies=TRUE)
# install.packages("phytools", dependencies=TRUE)
# remotes::install_github("fmichonneau/phyloch")
# install.packages("admixr")
# devtools::install_github("uqrmaie1/admixtools")
# devtools::install_github("uqrmaie1/admixtools", dependencies = TRUE, force = TRUE)
# install.packages("Rcpp")
# install.packages("tidyverse")
# install.packages("igraph")
# install.packages("plotly")
# install.packages("coalescentMCMC")

library(adegenet)
library(ape)
library(BiocManager)
# BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
library(stats)
library(ips)
# BiocManager::install("msa")
library(msa)
library(stringi)
library(stringr)
library(strap)
library(coalescentMCMC)

# Load multiple DNA sequences
fname = "countries.fasta"
fname2 = "Iso_RSRS_RCRS_countriesAlign.fasta"
file <- Biostrings::readDNAStringSet(fname)#for reading multiple DNA sequences from msa package
file2 <- Biostrings::readDNAStringSet(fname2)#for reading multiple DNA sequences from msa package
cb<-file2
nbin<-as.DNAbin(cb) #read aligned data from cb above

an<-as.alignment(nbin)  #converting DNAbin to alignment format
nm<-as.matrix(an)       #converting alignment to matrix
nbinmat<-as.matrix(labels(nbin)) #extraction of the sample names
nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin)

library(ggtree)
ggt<-ggtree::ggtree(tree, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)

library(treeio)
zoom(tree, grep("Vietnam", tree$tip.label, value = TRUE))

groupInfo <- split(tree$tip.label, gsub("_\\w+", "", tree$tip.label))
tree <- groupOTU(tree, groupInfo)
options(ignore.negative.edge=TRUE)
p <- ggtree(tree, aes(color=group)) + geom_tiplab() + xlim(NA, 23)
zoom(tree, grep("Vietnam", tree$tip.label), xmax_adjust=2)

tre<-ladderize(tree)
ggtree(tre, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)+
  geom_highlight(node = 20, fill="purple", alpha=0.2)

# njmsaplot<-msaplot(ggt, nbin, offset = 0.009, width=1, height = 0.5, color = c(rep("coral4", 1), rep("rosybrown", 1), rep("sienna1", 1), rep("lightgoldenrod1", 1), rep("darkseagreen2", 1), rep("lightskyblue1", 1)))
# njmsaplot
# 
# dev.new()
# njmsaplot
# 
# pdf("njmsaplot-sea.pdf", width = 11, height = 9)#save as pdf file
# njmsaplot
# dev.off()

# Phylogenetic tree

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)
library(tidytree)
library(ggstar)
library(TDbook)
library(phytools)
library(ape)

# tree2 <- ape::read.tree("Treefull_root.newick")

tree <- phytools::read.newick("Treefull.newick")
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$tip.label <- gsub("_", "", tree$tip.label)
tree$tip.label <- gsub("-", "", tree$tip.label)
tree$tip.label <- gsub(" ", "", tree$tip.label)
tree$tip.label <- gsub(",", "", tree$tip.label)
tree$tip.label <- trimws(gsub("\\s+", " ", tree$tip.label))

options(max.print=1000000)

ggt<-ggtree::ggtree(tree, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt

library(tidyr)
library(dplyr)
library(data.table)
library(readxl)
library(scales)

pal<-rgb(0,1,0)
show_col(pal)
show_col(hue_pal()(9))
show_col(hue_pal()(16), borders = NA)

show_col(viridis_pal()(16))
show_col(viridis_pal()(16), labels = FALSE)
show_col(viridis_pal()(141))
viridis_pal()(141)

show_col(c("#0099ff", "#163566", "#336699", "#339900", "#660000", "#66ccff", "#9900ff", "#990f80", "#996666", "#99ff99", "#cc0000", "#cc66ff", "#cccc33", "#fa0f0c", "#ff6633", "#ff9999", "#ffcc33", "#ffcc99", "#ffff00"))

dat <- read_excel("IsolateExplanation.xlsx")
dat$name <- gsub("\\+", "", dat$name)
dat$name <- gsub("\\'", "", dat$name)
dat$name <- gsub("\\(", "", dat$name)
dat$name <- gsub("\\)", "", dat$name)
dat$name <- gsub("\\_", "", dat$name)
dat$name <- gsub("\\-", "", dat$name)
dat$name <- gsub(" ", "", dat$name)
dat$name <- gsub("\\*", "", dat$name)
dat$name <- gsub("\\@", "", dat$name)
dat$name <- gsub("\\!", "", dat$name)
dat$name <- gsub("\\,", "", dat$name)
dat$name <- gsub("\\,", "", dat$name)
dat2 <- dat %>% filter(!name %in% tree$tip.label)
dat$name <- trimws(gsub("\\s+", " ", dat$name))
table(tree$tip.label %in% dat$name)
table(dat$name %in% tree$tip.label)
dat <- as.data.frame(dat)

dat <- dat %>%
  mutate(haplogroup2=ifelse(haplogroup2=="A+152"|haplogroup2=="A+152+16362"|haplogroup2=="A+152+16362+200", "A+",
                            ifelse(haplogroup2=="R+16189", "R+", haplogroup2)),
         Country_color=NA,
         Country_color=ifelse(Country=="Brunei", "#ff6633",
                              ifelse(Country=="Cambodia", "#ffff00",
                                     ifelse(Country=="Indonesia", "#9900ff",
                                            ifelse(Country=="Laos", "#0099ff",
                                                   ifelse(Country=="Malaysia", "#990f80",
                                                          ifelse(Country=="Myanmar", "#99ff99",
                                                                 ifelse(Country=="Philippines", "#cc66ff",
                                                                        ifelse(Country=="Singapore", "#ff9999",
                                                                               ifelse(Country=="Thailand", "#339900",
                                                                                      ifelse(Country=="Timor-Leste", "#66ccff",
                                                                                             ifelse(Country=="Vietnam", "#fa0f0c",
                                                                                                    ifelse(Country=="RSRS", "black",
                                                                                                           ifelse(Country=="rCRS", "orange", Country_color))))))))))))),
         `Language family`=ifelse(Ethnicity=="Mon", "Austroasiatic",
                                  ifelse(Ethnicity=="Hmong", "Hmong-Mien",
                                         ifelse(Ethnicity=="Shan", "Tai-Kadai",
                                                ifelse(Ethnicity=="Jehai (or Jahai)", "Austroasiatic",
                                                       ifelse(Ethnicity=="Temuan", "Austronesian",
                                                              ifelse(Ethnicity=="Maranao", "Austronesian",
                                                                     ifelse(Ethnicity=="Semelai", "Austroasiatic",
                                                                            ifelse(Ethnicity=="Bru (Brao)", "Austroasiatic",
                                                                                   ifelse(Ethnicity=="Jarai", "Austronesian",
                                                                                          ifelse(Ethnicity=="Kadazan-Dusun", "Austronesian",
                                                                                                 ifelse(Ethnicity=="Alor", "Austronesian",
                                                                                                        ifelse(Ethnicity=="Arakanese (or Rakhine)", "Sino-Tibetan",
                                                                                                               ifelse(Ethnicity=="Timorese", "Austronesian",
                                                                                                                      ifelse(Ethnicity=="Mang", "Austroasiatic",`Language family`)))))))))))))),
         `Language family`=ifelse(`Language family`=="Austronesian, Austroasiatic" | `Language family`=="Austroasiatic, Austronesian", "Austroasiatic, Austronesian, Sino-Tibetan",
                                  ifelse(`Language family`=="Austronesian, Spanish" | `Language family`=="Austronesian, Trans-New Guinea" | `Language family`=="Papuan, Austronesian, English", "Austronesian",
                                         ifelse(`Language family`=="Hmong-Mien" | `Language family`=="Hmong-Mien, Mongolic", "Hmong-Mien +",
                                                ifelse(`Language family`=="Indo-European, Sino-Tibetan" | `Language family`=="Sino-Tibetan, Austroasiatic, Tai-Kadai" | `Language family`=="Sino-Tibetan, Tai-Kada" | `Language family`=="Tai-Kadai, Sino-Tibetan", "Sino-Tibetan +",
                                                       ifelse(`Language family`=="Trans-New Guinea" | `Language family`=="Trans–New Guinea (Alor-Pantar, Papuan)", "Trans-New Guinea +",
                                                              ifelse(`Language family`=="Austronesian, Austroasiatic, Indo-European", "Austroasiatic, Austronesian + (xSino-Tibetan)",
                                                                     ifelse(`Language family`=="Austroasiatic, Tai-Kadai" | `Language family`=="Austroasiatic, Tai-Kadai, Hmong-Mien, Sino-Tibetan" | `Language family`=="Tai-Kadai, Hmong-Mien, Austroasiatic", "Austroasiatic, Tai-Kadai +",
                                                                            ifelse(`Language family`=="Sino-Tibetan", "Sino-Tibetan +", `Language family`))))))))) %>%
  droplevels()

# df <- df_Candidaauris_data
# tr <- tree_Candidaauris
# library(writexl)
# write_xlsx(dat, "SEA_megadata.xlsx")
# 
# dat <- read_excel("SEA_megadata.xlsx")
countries <- c("Brunei", "Cambodia", "Indonesia", "Laos", "Malaysia", "Myanmar", "Philippines", "rCRS", "RSRS", "Singapore", "Thailand", "Timor-Leste", "Vietnam")

# For the tip points
dat1 <- dat %>% dplyr::select(c("name", "Country", "Country_color", "haplogroup2"))
dat1$Country <- factor(dat1$Country, levels=countries)
Countrycolors <- dat1[match(countries,dat$Country),"Country_color"]

# For the haplogroup layer
dat3 <- dat %>% dplyr::select(c("name", "haplogroup2")) %>%
  melt(id="name", variable.name="haplo", value.name="haplogroup")

# For the clade group
dat4 <- dat %>% dplyr::select(c("name", "Language family"))
dat4 <- aggregate(.~`Language family`, dat4, FUN=paste, collapse=",")
clades <- lapply(dat4$name, function(x){unlist(strsplit(x,split=","))})
names(clades) <- dat4$`Language family`

# ggt<-ggtree::ggtree(tree, cex = 0.8, aes(color=branch.length))+
#   geom_tiplab(align=TRUE, size=2)+
#   geom_treescale(y = - 5, color = "coral4", fontsize = 4)
# ggt
# 
# tree <- groupOTU(tree, clades, "Clade")
# Clade <- NULL
# p <- ggtree(tr=tree, layout="fan", open.angle=15, size=0.2, aes(colour=Clade)) +
#   scale_colour_manual(
#     name="Language",
#     values=c("black", "#fa0f0c", "#FF61CC", "#ED68ED", "#F8766D", "#ABA300", "#C77CFF", "#8494FF", "#E68613", "#0099ff", "#660000", "#163566", "#336699", "#339900", "#66ccff", "lightgrey"),
#     labels=c("RSRS", "Austroasiatic", "Austroasiatic, Austronesian", "Austroasiatic, Austronesian + (xSino-Tibetan)", "Austroasiatic, Austronesian, Sino-Tibetan", "Austroasiatic, Tai-Kadai +", "Austronesian", "Austronesian + (xAustroasiatic)", "rCRS", "Hmong-Mien +", "Mayan", "Sino-Tibetan", "Sino-Tibetan +", "Tai-Kadai", "Trans-New Guinea +", "Unknown"),
#     guide=guide_legend(keywidth=1.5,
#                        keyheight=1.25,
#                        order=1,
#                        override.aes=list(linetype=1, size=6,alpha=1))) +
#   theme(legend.position="right",
#         legend.background=element_rect(fill=NA),
#         legend.title=element_text(size=20),
#         legend.text=element_text(size=15),
#         legend.key.size = unit(40, "cm"),
#         legend.spacing.y = unit(2, "cm")) + 
#   new_scale_colour()
# 
# p1 <- p %<+% dat1 +
#   geom_tippoint(aes(colour=Country), alpha=0) +
#   geom_tiplab(aes(colour=Country),
#               align=TRUE,
#               linetype=3,
#               size=1,
#               linesize=0.2,
#               show.legend=FALSE) +
#   scale_colour_manual(
#     name="Country",
#     values=Countrycolors, 
#     guide=guide_legend(keywidth=1.5,
#                        keyheight=1.25,
#                        order=2,
#                        override.aes=list(size=5,alpha=1))) +
#   theme(legend.position="right",
#         legend.background=element_rect(fill=NA),
#         legend.title=element_text(size=20),
#         legend.text=element_text(size=15),
#         legend.key.size = unit(30, "cm"),
#         legend.spacing.y = unit(2, "cm")) + 
#   new_scale_colour()
# 
# p3 <- p1 +
#   geom_fruit(data=dat3,
#              geom=geom_star,
#              mapping=aes(x=haplogroup, y=name, fill=haplogroup, starshape=haplo),
#              size=3,
#              starstroke=0,
#              pwidth=0.1,
#              inherit.aes = FALSE,
#              grid.params=list(linetype=3, size=0.2)) +
#   scale_fill_discrete(
#     name="Haplogroup",
#     guide=guide_legend(keywidth=1.5, 
#                        keyheight=1.25, 
#                        order=3,
#                        override.aes=list(size=5)),
#     na.translate=FALSE) +
#   scale_starshape_discrete(guide="none") +
#   theme(legend.background=element_rect(fill=NA),
#         legend.title=element_text(size=20), 
#         legend.text=element_text(size=15),
#         legend.key.size = unit(30, "cm"),
#         legend.spacing.y = unit(2, "cm"))
# p3
# ggsave(filename = file.path("figures", "Treefull.png"), width = 30, height = 20)

filename <- "TREEFULL.NEXUS"
tree2 <- ape::read.nexus(filename)
tree2$tip.label <- gsub("'", "", tree2$tip.label)
tree2$tip.label <- gsub("_", "", tree2$tip.label)
tree2$tip.label <- gsub("-", "", tree2$tip.label)
tree2$tip.label <- gsub(" ", "", tree2$tip.label)
tree2$tip.label <- gsub(",", "", tree2$tip.label)
tree2$tip.label <- trimws(gsub("\\s+", " ", tree2$tip.label))

info <- dat
cols <- Countrycolors

metadata <- dat %>%
  mutate(Language_color=NA,
         Language_color=ifelse(`Language family`=="RSRS", "black",
                               ifelse(`Language family`=="Austroasiatic", "#fa0f0c",
                                      ifelse(`Language family`=="Austroasiatic, Austronesian, Sino-Tibetan", "#cccc33",
                                             ifelse(`Language family`=="Austroasiatic, Tai-Kadai +", "#ffcc99",
                                                    ifelse(`Language family`=="Austronesian", "#9900ff",
                                                           ifelse(`Language family`=="rCRS", "white",
                                                                  ifelse(`Language family`=="Hmong-Mien +", "#0099ff",
                                                                         ifelse(`Language family`=="Mayan", "#660000",
                                                                                ifelse(`Language family`=="Sino-Tibetan +", "#336699",
                                                                                       ifelse(`Language family`=="Tai-Kadai", "#339900",
                                                                                              ifelse(`Language family`=="Trans-New Guinea +", "#66ccff",
                                                                                                     ifelse(`Language family`=="Unknown", "lightgrey",
                                                                                                            Language_color)))))))))))),
         Ethnicity=ifelse(Ethnicity=="Lao, Akha, Hmong, Khmu, Yao/Mien, Phuan", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan",
                          ifelse(Ethnicity=="Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan",
                                 ifelse(Ethnicity=="Yao", "Dao", Ethnicity))),
         Ethnicity_color=NA,
         Ethnicity_color=ifelse(Ethnicity %in% c("Abaknon", "Alor", "Ambelau, Ambonese", "Ambonese", "Balantak, Bali Aga, Balinese", "Bali Aga, Balinese", "Banjar", "Banjar, Bantenese, Banyumasan", "Banjar, Dayak, Javanese", "Batak", "Batak, Acehnese", "Batak, Minangkabau, Acehnese, Lampung", "Bicolano", "Bidayuh", "Bruneian Malay", "Bugis", "Bugis, Malay", "Bugkalot (or Ilongot)", "Cebuano", "Cebuano - Filipino", "Cham", "Cuyunin (or Cuyonon)", "Dayak", "Filipino", "Filipino (or Tagalog)", "Ibaloi", "Ifugao", "Igorot", "Indonesian", "Ivatan", "Jarai", "Javanese", "Javanese, Malay", "Javanese, Palembang, Batak, Minangkabau, Komering", "Kadazan-Dusun", "Kankanaey", "Makassarese", "Malay", "Malay, Achehnese", "Malay, Banjar Malay", "Maranao", "Melanau", "Minahasa", "Minangkabau", "Moken", "Palembangese", "Papuan", "Seletar (or Orang Seletar)", "Semende", "SiLa", "Sumbanese", "Sundanese", "Surigaonon", "Tagalog", "Temuan", "Tetum", "The Kalanguya (or Ikalahan)", "Timorese", "Toraja", "UrakLawoi", "Zambal"), "#9900ff",
                                ifelse(Ethnicity %in% c("Achang", "Aini", "Arakanese (or Rakhine)", "Bamar (or Burman)", "Chinese", "Dai", "Deang", "HaNhi", "Hui", "Jingpo", "Karen", "Lahu", "LaHu", "Naga"), "#336699",
                                       ifelse(Ethnicity=="RSRS", "black",
                                              ifelse(Ethnicity %in% c("Akar", "Akar Jambat", "Bru (Brao)", "Ede", "Giarai", "Jehai (or Jahai)", "Khmer", "Khuen", "Kinh", "Kintaq", "Kreung", "Mang", "Mon", "PaThen", "Phnong", "PhuLa", "Semelai", "Stieng", "Tompoun"), "#fa0f0c",
                                                     ifelse(Ethnicity %in% c("Batek", "Jahai, Semang", "Kensiu", "Khmer, Cham, Chinese-Cambodian, Vietnamese", "Lisu", "LoLo", "Orang Asli"), "#cccc33",
                                                            ifelse(Ethnicity %in% c("Bunak", "Fataluku", "Kemak", "Makasae", "Makassai", "Mambai"), "#66ccff",
                                                                   ifelse(Ethnicity %in% c("CoLao", "Isan (or Lao)", "LaChi", "Lao", "Lao Islan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Nung", "Phutai", "Shan", "Tay", "Tay Nung", "Thai"), "#339900",
                                                                          ifelse(Ethnicity %in% c("Dao", "Hmong", "IuMien"), "#0099ff",
                                                                                 ifelse(Ethnicity %in% c("rCRS"), "white",
                                                                                        ifelse(Ethnicity %in% c("Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan"), "#ffcc99",
                                                                                               ifelse(Ethnicity %in% c("Mam"), "#660000",
                                                                                                      ifelse(Ethnicity %in% c("Unknown"), "lightgrey", Ethnicity))))))))))))) %>%
  droplevels()

metadata <- metadata %>%
  left_join(dat1 %>% dplyr::select(name, Country_color)) %>% 
  dplyr::select(c("name", "Country", "Country_color",
           "Language family", "Language_color", "haplogroup1", "Ethnicity", "Ethnicity_color")) %>%
  droplevels()

languages <- c("Austroasiatic", "Austroasiatic, Austronesian, Sino-Tibetan", "Austroasiatic, Tai-Kadai +", "Austronesian", "Hmong-Mien +", "Mayan", "rCRS", "RSRS", "Sino-Tibetan +", "Tai-Kadai", "Trans-New Guinea +", "Unknown")
Languagecolors <- metadata[match(languages, metadata$`Language family`),"Language_color"]

ethnics <- c("Abaknon", "Achang", "Aini", "Akar", "Akar Jambat", "Alor", "Ambelau, Ambonese", "Ambonese", "Arakanese (or Rakhine)", "Balantak, Bali Aga, Balinese", "Bali Aga, Balinese", "Bamar (or Burman)", "Banjar", "Banjar, Bantenese, Banyumasan", "Banjar, Dayak, Javanese", "Batak", "Batak, Acehnese", "Batak, Minangkabau, Acehnese, Lampung", "Batek", "Bicolano", "Bidayuh", "Bru (Brao)", "Bruneian Malay", "Bugis", "Bugis, Malay", "Bugkalot (or Ilongot)", "Bunak", "Cebuano", "Cebuano - Filipino", "Cham", "Chinese", "CoLao", "Cuyunin (or Cuyonon)", "Dai", "Dao", "Dayak", "Deang", "Ede", "Fataluku", "Filipino", "Filipino (or Tagalog)", "Giarai", "HaNhi", "Hmong", "Hui", "Ibaloi", "Ifugao", "Igorot", "Indonesian", "Isan (or Lao)", "IuMien", "Ivatan", "Jahai, Semang", "Jarai", "Javanese", "Javanese, Malay", "Javanese, Palembang, Batak, Minangkabau, Komering", "Jehai (or Jahai)", "Jingpo", "Kadazan-Dusun", "Kankanaey", "Karen", "Kemak", "Kensiu", "Khmer", "Khmer, Cham, Chinese-Cambodian, Vietnamese", "Khuen", "Kinh", "Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "Kintaq", "Kreung", "LaChi", "Lahu", "LaHu", "Lao", "Lao Islan", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Lisu", "LoLo", "Makasae", "Makassai", "Makassarese", "Malay", "Malay, Achehnese", "Malay, Banjar Malay", "Mam", "Mambai", "Mang", "Maranao", "Melanau", "Minahasa", "Minangkabau", "Moken", "Mon", "Naga", "Nung", "Orang Asli", "Palembangese", "Papuan", "PaThen", "Phnong", "PhuLa", "Phutai", "rCRS", "RSRS", "Seletar (or Orang Seletar)", "Semelai", "Semende", "Shan", "SiLa", "Stieng", "Sumbanese", "Sundanese", "Surigaonon", "Tagalog", "Tay", "Tay Nung", "Temuan", "Tetum", "Thai", "The Kalanguya (or Ikalahan)", "Timorese", "Tompoun", "Toraja", "Unknown", "UrakLawoi", "Zambal")
Ethniccolors <- metadata[match(ethnics, metadata$Ethnicity), "Ethnicity_color"]

p <- ggtree(tree2, layout='circular')

p <- p %<+% metadata

p1 <- p +
  geom_tippoint(aes(color=Country),
                size=4) + 
  scale_color_manual(values=cols, 
                     guide=guide_legend(keywidth=2,
                                        keyheight=3,
                                        order=2,
                                        override.aes=list(size=10,alpha=1))) + 
  theme(legend.position="right",
                          legend.background=element_rect(fill=NA),
                          legend.title=element_text(size=30, face="bold"),
                          legend.text=element_text(size=20),
                          legend.key.size = unit(30, "cm"),
                          legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

# p2 <-p1 +
#   geom_fruit(
#     geom=geom_tile,
#     mapping=aes(fill=`Language family`),
#     width=0.01,
#     offset=0.1
#   ) +
#   scale_fill_manual(
#     name="Language",
#     values=Languagecolors,
#     guide=guide_legend(keywidth=2, 
#                        keyheight=1, 
#                        ncol=1, 
#                        order=1,
#                        override.aes=list(size=5,alpha=1))
#   ) +
#   theme(legend.background=element_rect(fill=NA),
#         legend.title=element_text(size=30, face="bold"), 
#         legend.text=element_text(size=20),
#         legend.key.size = unit(30, "cm"),
#         legend.spacing.y = unit(2, "cm")) + 
#   new_scale_colour()
# 
# p3 <- p2 +
#   new_scale_fill() +
#   geom_fruit(geom=geom_star,
#              mapping=aes(fill=haplogroup1),
#              size=10,
#              starstroke=0,
#              pwidth=0.1,
#              inherit.aes = FALSE,
#              grid.params=list(linetype=3, size=0.2)) +
#   scale_fill_discrete(
#     name="Haplogroup",
#     guide=guide_legend(keywidth=2, 
#                        keyheight=1, 
#                        ncol=10,
#                        order=3,
#                        override.aes=list(size=10)),
#     na.translate=FALSE) +
#   scale_starshape_discrete(guide="none") +
#   theme(legend.background=element_rect(fill=NA),
#         legend.title=element_text(size=30, face="bold"), 
#         legend.text=element_text(size=20),
#         legend.key.size = unit(30, "cm"),
#         legend.spacing.y = unit(2, "cm")) + 
#   new_scale_colour()
# 
# p3
# ggsave(filename = file.path("figures", "Treefull(2).png"), width = 30, height = 20)

p4 <-p1 +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=`Language family`),
    width=0.006,
    offset=0.03
  ) +
  scale_fill_manual(
    name="Language",
    values=Languagecolors,
    guide=guide_legend(keywidth=3, 
                       keyheight=3, 
                       ncol=1, 
                       order=1,
                       override.aes=list(size=10,alpha=1))
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"), 
        legend.text=element_text(size=20),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

p5 <- p4 +
  new_scale_fill() +
  geom_fruit(geom=geom_star,
             mapping=aes(fill=haplogroup1),
             size=15,
             starstroke=0,
             pwidth=0.2,
             inherit.aes = FALSE,
             grid.params=list(linetype=3, size=0.2)) +
  scale_fill_discrete(
    name="Haplogroup",
    guide=guide_legend(keywidth=2, 
                       keyheight=2, 
                       ncol=10,
                       order=3,
                       override.aes=list(size=10)),
    na.translate=FALSE) +
  scale_starshape_discrete(guide="none") +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"), 
        legend.text=element_text(size=20),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

p5
ggsave(filename = file.path("figures", "Treefull(3).png"), width = 30, height = 20)

# Final tree

p <- ggtree(tree2, layout='circular')

p <- p %<+% metadata

p1 <- p +
  geom_tippoint(aes(color=Country),
                size=5) + 
  scale_color_manual(values=cols, 
                     guide=guide_legend(keywidth=2,
                                        keyheight=3,
                                        order=2,
                                        override.aes=list(size=15,alpha=1))) + 
  theme(legend.position="right",
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"),
        legend.text=element_text(size=24),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical") + 
  new_scale_colour()

p4 <-p1 +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=`Language family`),
    width=0.005,
    offset=0.03
  ) +
  scale_fill_manual(
    name="Language",
    values=Languagecolors,
    guide=guide_legend(keywidth=3, 
                       keyheight=3, 
                       order=1,
                       override.aes=list(size=10,alpha=1))
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"), 
        legend.text=element_text(size=24),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical") + 
  new_scale_colour()

p6 <- p +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Ethnicity),
    width=0.005,
    offset=0.03
  ) +
  scale_fill_manual(
    name="Ethnicity",
    values=Ethniccolors,
    guide=guide_legend(keywidth=2, 
                       keyheight=2, 
                       ncol = 2,
                       override.aes=list(size=1,alpha=1),
                       title.position="top")
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=36, face="bold"), 
        legend.text=element_text(size=22),
        legend.key.size = unit(10, "cm"),
        legend.spacing.y = unit(2, "cm"),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        legend.position = "right", 
        legend.box.spacing = unit(-30, "pt")) + 
  new_scale_colour()

p7 <- p4 +
  new_scale_fill() +
  geom_fruit(geom=geom_star,
             mapping=aes(fill=haplogroup1),
             size=18,
             starstroke=0,
             pwidth=0.2,
             inherit.aes = FALSE,
             grid.params=list(linetype=3, size=0.2)) +
  scale_fill_discrete(
    name="Haplogroup",
    guide=guide_legend(keywidth=1, 
                       keyheight=1, 
                       order=3,
                       nrow = 4, byrow = TRUE,
                       override.aes=list(size=10)),
    na.translate=FALSE) +
  scale_starshape_discrete(guide="none") +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=40, face="bold"), 
        legend.text=element_text(size=30),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(1, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.position = "left", 
        legend.box.spacing = unit(-50, "pt")) + 
  new_scale_colour()
p7

library(cowplot)   # get_legend() & plot_grid() functions
library(patchwork) # blank plot: plot_spacer()

# Get legends

leg6 <- get_legend(p6)

# create a blank plot for legend alignment 
blank_p <- plot_spacer() + theme_void()

# Make a title

title <- ggdraw() + draw_label("Southeast Asia", fontface='bold') +
  theme(plot.title = element_text(size = 100))
  
# Put everything together

final_p <- plot_grid(p7, leg6, nrow = 1, align = "h", axis = "r", rel_widths = c(1.8,1)) +
  theme(plot.background = element_rect(fill = "white", colour = NA)) +
  ggtitle("Southeast Asia") + 
  theme(plot.title = element_text(size = 100, face = "bold", hjust = 0),
        plot.title.position = "plot")

ggsave(filename = file.path("figures", "Treefull(5).png"), width = 49, height = 33)

## Vietnam

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)
library(tidytree)
library(ggstar)
library(TDbook)
library(phytools)
library(ape)

tree <- phytools::read.newick("Vietnam.newick")
tree$tip.label <- gsub("'", "", tree$tip.label)
tree$tip.label <- gsub("_", "", tree$tip.label)
tree$tip.label <- gsub("-", "", tree$tip.label)
tree$tip.label <- gsub(" ", "", tree$tip.label)
tree$tip.label <- gsub(",", "", tree$tip.label)
tree$tip.label <- trimws(gsub("\\s+", " ", tree$tip.label))

options(max.print=1000000)

library(tidyr)
library(dplyr)
library(data.table)
library(readxl)
library(scales)

dat <- read_excel("IsolateExplanation.xlsx")
dat$name <- gsub("\\+", "", dat$name)
dat$name <- gsub("\\'", "", dat$name)
dat$name <- gsub("\\(", "", dat$name)
dat$name <- gsub("\\)", "", dat$name)
dat$name <- gsub("\\_", "", dat$name)
dat$name <- gsub("\\-", "", dat$name)
dat$name <- gsub(" ", "", dat$name)
dat$name <- gsub("\\*", "", dat$name)
dat$name <- gsub("\\@", "", dat$name)
dat$name <- gsub("\\!", "", dat$name)
dat$name <- gsub("\\,", "", dat$name)
dat$name <- gsub("\\,", "", dat$name)
dat2 <- dat %>% filter(!name %in% tree$tip.label)
dat$name <- trimws(gsub("\\s+", " ", dat$name))
table(tree$tip.label %in% dat$name)
table(dat$name %in% tree$tip.label)
dat <- as.data.frame(dat)

dat <- dat %>%
  mutate(haplogroup2=ifelse(haplogroup2=="A+152"|haplogroup2=="A+152+16362"|haplogroup2=="A+152+16362+200", "A+",
                            ifelse(haplogroup2=="R+16189", "R+", haplogroup2)),
         Country_color=NA,
         Country_color=ifelse(Country=="Brunei", "#ff6633",
                              ifelse(Country=="Cambodia", "#ffff00",
                                     ifelse(Country=="Indonesia", "#9900ff",
                                            ifelse(Country=="Laos", "#0099ff",
                                                   ifelse(Country=="Malaysia", "#990f80",
                                                          ifelse(Country=="Myanmar", "#99ff99",
                                                                 ifelse(Country=="Philippines", "#cc66ff",
                                                                        ifelse(Country=="Singapore", "#ff9999",
                                                                               ifelse(Country=="Thailand", "#339900",
                                                                                      ifelse(Country=="Timor-Leste", "#66ccff",
                                                                                             ifelse(Country=="Vietnam", "#fa0f0c",
                                                                                                    ifelse(Country=="RSRS", "black",
                                                                                                           ifelse(Country=="rCRS", "orange", Country_color))))))))))))),
         `Language family`=ifelse(Ethnicity=="Mon", "Austroasiatic",
                                  ifelse(Ethnicity=="Hmong", "Hmong-Mien",
                                         ifelse(Ethnicity=="Shan", "Tai-Kadai",
                                                ifelse(Ethnicity=="Jehai (or Jahai)", "Austroasiatic",
                                                       ifelse(Ethnicity=="Temuan", "Austronesian",
                                                              ifelse(Ethnicity=="Maranao", "Austronesian",
                                                                     ifelse(Ethnicity=="Semelai", "Austroasiatic",
                                                                            ifelse(Ethnicity=="Bru (Brao)", "Austroasiatic",
                                                                                   ifelse(Ethnicity=="Jarai", "Austronesian",
                                                                                          ifelse(Ethnicity=="Kadazan-Dusun", "Austronesian",
                                                                                                 ifelse(Ethnicity=="Alor", "Austronesian",
                                                                                                        ifelse(Ethnicity=="Arakanese (or Rakhine)", "Sino-Tibetan",
                                                                                                               ifelse(Ethnicity=="Timorese", "Austronesian",
                                                                                                                      ifelse(Ethnicity=="Mang", "Austroasiatic",`Language family`)))))))))))))),
         `Language family`=ifelse(`Language family`=="Austronesian, Austroasiatic" | `Language family`=="Austroasiatic, Austronesian", "Austroasiatic, Austronesian, Sino-Tibetan",
                                  ifelse(`Language family`=="Austronesian, Spanish" | `Language family`=="Austronesian, Trans-New Guinea" | `Language family`=="Papuan, Austronesian, English", "Austronesian",
                                         ifelse(`Language family`=="Hmong-Mien" | `Language family`=="Hmong-Mien, Mongolic", "Hmong-Mien +",
                                                ifelse(`Language family`=="Indo-European, Sino-Tibetan" | `Language family`=="Sino-Tibetan, Austroasiatic, Tai-Kadai" | `Language family`=="Sino-Tibetan, Tai-Kada" | `Language family`=="Tai-Kadai, Sino-Tibetan", "Sino-Tibetan +",
                                                       ifelse(`Language family`=="Trans-New Guinea" | `Language family`=="Trans–New Guinea (Alor-Pantar, Papuan)", "Trans-New Guinea +",
                                                              ifelse(`Language family`=="Austronesian, Austroasiatic, Indo-European", "Austroasiatic, Austronesian + (xSino-Tibetan)",
                                                                     ifelse(`Language family`=="Austroasiatic, Tai-Kadai" | `Language family`=="Austroasiatic, Tai-Kadai, Hmong-Mien, Sino-Tibetan" | `Language family`=="Tai-Kadai, Hmong-Mien, Austroasiatic", "Austroasiatic, Tai-Kadai +",
                                                                            ifelse(`Language family`=="Sino-Tibetan", "Sino-Tibetan +", `Language family`))))))))) %>%
  droplevels() %>%
  filter(Country=="Vietnam")

countries <- c("Vietnam")

# For the tip points
dat1 <- dat %>% dplyr::select(c("name", "Country", "Country_color", "haplogroup2"))
dat1$Country <- factor(dat1$Country, levels=countries)
Countrycolors <- dat1[match(countries,dat$Country),"Country_color"]

filename <- "Vietnam.NEXUS"
tree2 <- ape::read.nexus(filename)
tree2$tip.label <- gsub("'", "", tree2$tip.label)
tree2$tip.label <- gsub("_", "", tree2$tip.label)
tree2$tip.label <- gsub("-", "", tree2$tip.label)
tree2$tip.label <- gsub(" ", "", tree2$tip.label)
tree2$tip.label <- gsub(",", "", tree2$tip.label)
tree2$tip.label <- trimws(gsub("\\s+", " ", tree2$tip.label))

info <- dat
cols <- Countrycolors

metadata <- dat %>%
  mutate(Language_color=NA,
         Language_color=ifelse(`Language family`=="RSRS", "black",
                               ifelse(`Language family`=="Austroasiatic", "#fa0f0c",
                                      ifelse(`Language family`=="Austroasiatic, Austronesian, Sino-Tibetan", "#cccc33",
                                             ifelse(`Language family`=="Austroasiatic, Tai-Kadai +", "#ffcc99",
                                                    ifelse(`Language family`=="Austronesian", "#9900ff",
                                                           ifelse(`Language family`=="rCRS", "white",
                                                                  ifelse(`Language family`=="Hmong-Mien +", "#0099ff",
                                                                         ifelse(`Language family`=="Mayan", "#660000",
                                                                                ifelse(`Language family`=="Sino-Tibetan +", "#336699",
                                                                                       ifelse(`Language family`=="Tai-Kadai", "#339900",
                                                                                              ifelse(`Language family`=="Trans-New Guinea +", "#66ccff",
                                                                                                     ifelse(`Language family`=="Unknown", "lightgrey",
                                                                                                            Language_color)))))))))))),
         Ethnicity=ifelse(Ethnicity=="Lao, Akha, Hmong, Khmu, Yao/Mien, Phuan", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan",
                          ifelse(Ethnicity=="Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan",
                                 ifelse(Ethnicity=="Yao", "Dao", Ethnicity))),
         Ethnicity_color=NA,
         Ethnicity_color=ifelse(Ethnicity %in% c("Abaknon", "Alor", "Ambelau, Ambonese", "Ambonese", "Balantak, Bali Aga, Balinese", "Bali Aga, Balinese", "Banjar", "Banjar, Bantenese, Banyumasan", "Banjar, Dayak, Javanese", "Batak", "Batak, Acehnese", "Batak, Minangkabau, Acehnese, Lampung", "Bicolano", "Bidayuh", "Bruneian Malay", "Bugis", "Bugis, Malay", "Bugkalot (or Ilongot)", "Cebuano", "Cebuano - Filipino", "Cham", "Cuyunin (or Cuyonon)", "Dayak", "Filipino", "Filipino (or Tagalog)", "Ibaloi", "Ifugao", "Igorot", "Indonesian", "Ivatan", "Jarai", "Javanese", "Javanese, Malay", "Javanese, Palembang, Batak, Minangkabau, Komering", "Kadazan-Dusun", "Kankanaey", "Makassarese", "Malay", "Malay, Achehnese", "Malay, Banjar Malay", "Maranao", "Melanau", "Minahasa", "Minangkabau", "Moken", "Palembangese", "Papuan", "Seletar (or Orang Seletar)", "Semende", "SiLa", "Sumbanese", "Sundanese", "Surigaonon", "Tagalog", "Temuan", "Tetum", "The Kalanguya (or Ikalahan)", "Timorese", "Toraja", "UrakLawoi", "Zambal"), "#9900ff",
                                ifelse(Ethnicity %in% c("Achang", "Aini", "Arakanese (or Rakhine)", "Bamar (or Burman)", "Chinese", "Dai", "Deang", "HaNhi", "Hui", "Jingpo", "Karen", "Lahu", "LaHu", "Naga"), "#336699",
                                       ifelse(Ethnicity=="RSRS", "black",
                                              ifelse(Ethnicity %in% c("Akar", "Akar Jambat", "Bru (Brao)", "Ede", "Giarai", "Jehai (or Jahai)", "Khmer", "Khuen", "Kinh", "Kintaq", "Kreung", "Mang", "Mon", "PaThen", "Phnong", "PhuLa", "Semelai", "Stieng", "Tompoun"), "#fa0f0c",
                                                     ifelse(Ethnicity %in% c("Batek", "Jahai, Semang", "Kensiu", "Khmer, Cham, Chinese-Cambodian, Vietnamese", "Lisu", "LoLo", "Orang Asli"), "#cccc33",
                                                            ifelse(Ethnicity %in% c("Bunak", "Fataluku", "Kemak", "Makasae", "Makassai", "Mambai"), "#66ccff",
                                                                   ifelse(Ethnicity %in% c("CoLao", "Isan (or Lao)", "LaChi", "Lao", "Lao Islan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Nung", "Phutai", "Shan", "Tay", "Tay Nung", "Thai"), "#339900",
                                                                          ifelse(Ethnicity %in% c("Dao", "Hmong", "IuMien"), "#0099ff",
                                                                                 ifelse(Ethnicity %in% c("rCRS"), "white",
                                                                                        ifelse(Ethnicity %in% c("Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan"), "#ffcc99",
                                                                                               ifelse(Ethnicity %in% c("Mam"), "#660000",
                                                                                                      ifelse(Ethnicity %in% c("Unknown"), "lightgrey", Ethnicity))))))))))))) %>%
  droplevels()

metadata <- metadata %>%
  left_join(dat1 %>% dplyr::select(name, Country_color)) %>% 
  dplyr::select(c("name", "Country", "Country_color",
                  "Language family", "Language_color", "haplogroup1", "Ethnicity", "Ethnicity_color")) %>%
  droplevels()

languages <- c("Austroasiatic", "Austroasiatic, Tai-Kadai +", "Austronesian", "Hmong-Mien +", "Sino-Tibetan +", "Tai-Kadai")
Languagecolors <- metadata[match(languages, metadata$`Language family`),"Language_color"]

ethnics <- c("Cham", "CoLao", "Dao", "Ede", "Giarai", "HaNhi", "Hmong", "Hui", "Kinh", "Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "LaChi", "LaHu", "LoLo", "Mang", "Nung", "PaThen", "PhuLa", "SiLa", "Stieng", "Tay", "Tay Nung", "Thai")
Ethniccolors <- metadata[match(ethnics, metadata$Ethnicity), "Ethnicity_color"]

p <- ggtree(tree2, layout='circular')

p <- p %<+% metadata

p1 <- p +
  geom_tippoint(aes(color=Country),
                size=5) + 
  scale_color_manual(values=cols, 
                     guide=guide_legend(keywidth=2,
                                        keyheight=3,
                                        order=2,
                                        override.aes=list(size=15,alpha=1))) + 
  theme(legend.position="right",
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"),
        legend.text=element_text(size=24),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

p4 <-p1 +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=`Language family`),
    width=0.0003,
    offset=0.03
  ) +
  scale_fill_manual(
    name="Language",
    values=Languagecolors,
    guide=guide_legend(keywidth=3, 
                       keyheight=3, 
                       order=1,
                       override.aes=list(size=10,alpha=1))
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=30, face="bold"), 
        legend.text=element_text(size=24),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

p6 <- p4 +
  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Ethnicity),
    width=0.0003,
    offset=0.03
  ) +
  scale_fill_manual(
    name="Ethnicity",
    values=Ethniccolors,
    guide=guide_legend(keywidth=2, 
                       keyheight=3, 
                       ncol = 1,
                       order = 4,
                       override.aes=list(size=1,alpha=1),
                       title.position="top")
  ) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=36, face="bold"), 
        legend.text=element_text(size=22),
        legend.key.size = unit(10, "cm"),
        legend.spacing.y = unit(2, "cm")) + 
  new_scale_colour()

p7 <- p6 +
  new_scale_fill() +
  geom_fruit(geom=geom_star,
             mapping=aes(fill=haplogroup1),
             size=30,
             starstroke=0,
             pwidth=0.2,
             inherit.aes = FALSE,
             grid.params=list(linetype=3, size=0.2)) +
  scale_fill_discrete(
    name="Haplogroup",
    guide=guide_legend(keywidth=1, 
                       keyheight=1, 
                       order=3,
                       nrow = 2, byrow = TRUE,
                       override.aes=list(size=10)),
    na.translate=FALSE) +
  scale_starshape_discrete(guide="none") +
  ggtitle("Vietnam") + 
  theme(plot.title = element_text(size = 100, face = "bold")) +
  theme(legend.background=element_rect(fill=NA),
        legend.title=element_text(size=40, face="bold"), 
        legend.text=element_text(size=30),
        legend.key.size = unit(30, "cm"),
        legend.spacing.y = unit(1, "cm"),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.position = "right") + 
  new_scale_colour()
p7

ggsave(filename = file.path("figures", "Tree_Vietnam.png"), width = 49, height = 33)

# Create dataset of sequences and subclades

library(tidyr)
library(dplyr)
library(data.table)
library(readxl)

dat <- read_excel("IsolateExplanation.xlsx")
dat <- as.data.frame(dat)
names <- dat$names

library(janitor)
dat2 <- dat %>%
  mutate(Ethnicity=ifelse(Ethnicity=="Lao, Akha, Hmong, Khmu, Yao/Mien, Phuan", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan",
                          ifelse(Ethnicity=="Yao", "Dao", Ethnicity))) %>%
  dplyr::select(name, Country, `Language family`, Ethnicity, haplo, haplogroup1, haplogroup2, haplogroup3) %>% setDT()
hap_rank <- dat2[, .N, by = .(haplogroup1)] %>% arrange(desc(N))
hap2_rank <- dat2[, .N, by = .(haplogroup2)] %>% arrange(desc(N))
hap3_rank <- dat2[, .N, by = .(haplogroup3)] %>% arrange(desc(N))

hap <- dat2[, .N, by = .(haplogroup1)] %>% arrange(desc(N)) %>% mutate(percent=(N*100)/sum(N)) %>% ungroup()
hap2 <- dat2[, .N, by = .(haplogroup2)] %>% arrange(desc(N)) %>% mutate(percent=(N*100)/sum(N)) %>% ungroup()
hap3 <- dat2[, .N, by = .(haplogroup3)] %>% arrange(desc(N)) %>% mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(pegas)
library(ape)

### Calculate nucleotide diversity, haplotype diversity, MPD and write out all subclades

dat_hap <- NULL
for (i in hap_rank$haplogroup1) {
  cat(i, "\n")
  df_i <- dat2 %>% dplyr::filter(haplogroup1==i)
  n_i <- nrow(df_i)
  nbin_i <- nbin[labels(nbin) %in% df_i$name]
  h_i <- pegas::haplotype(nbin_i)
  n_hi <- length(as.list(h_i))
  fas_i <- file2[labels(nbin) %in% df_i$name]
  writeXStringSet(fas_i, paste0("data/haplo/", i, "Aligned.fasta"))
  # dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
  x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
  dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
  mpd_i <- mean(dist.gene_i)
  hap.div_i <- hap.div(nbin_i, variance = TRUE)
  nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
  dati <- data.frame(df_i$haplogroup1, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
  ## combine
  dat_hap <- rbindlist(l = list(dat_hap, dati)) %>% unique() %>% setDT()
}

## Rename
setnames(x = dat_hap,
         old = c("df_i.haplogroup1", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
         new = c("Haplogroup", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

dat_hap2 <- NULL
for (i in hap2_rank$haplogroup2) {
  cat(i, "\n")
  df_i <- dat2 %>% dplyr::filter(haplogroup2==i)
  n_i <- nrow(df_i)
  nbin_i <- nbin[labels(nbin) %in% df_i$name]
  h_i <- pegas::haplotype(nbin_i)
  n_hi <- length(as.list(h_i))
  fas_i <- file2[labels(nbin) %in% df_i$name]
  writeXStringSet(fas_i, paste0("data/haplo/", i, "Aligned.fasta"))
  # dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
  x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
  dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
  mpd_i <- mean(dist.gene_i)
  hap.div_i <- hap.div(nbin_i, variance = TRUE)
  nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
  dati <- data.frame(df_i$haplogroup2, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
  ## combine
  dat_hap2 <- rbindlist(l = list(dat_hap2, dati)) %>% unique() %>% setDT()
}

## Rename
setnames(x = dat_hap2,
         old = c("df_i.haplogroup2", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
         new = c("Haplogroup", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

dat_hap3 <- NULL
for (i in hap3_rank$haplogroup3) {
  cat(i, "\n")
  df_i <- dat2 %>% dplyr::filter(haplogroup3==i)
  n_i <- nrow(df_i)
  nbin_i <- nbin[labels(nbin) %in% df_i$name]
  h_i <- pegas::haplotype(nbin_i)
  n_hi <- length(as.list(h_i))
  fas_i <- file2[labels(nbin) %in% df_i$name]
  writeXStringSet(fas_i, paste0("data/haplo/", i, "Aligned.fasta"))
  # dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
  x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
  dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
  mpd_i <- mean(dist.gene_i)
  hap.div_i <- hap.div(nbin_i, variance = TRUE)
  nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
  dati <- data.frame(df_i$haplogroup3, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
  ## combine
  dat_hap3 <- rbindlist(l = list(dat_hap3, dati)) %>% unique() %>% setDT()
}

## Rename
setnames(x = dat_hap3,
         old = c("df_i.haplogroup3", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
         new = c("Haplogroup", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

## Combine
dat_hap <- rbindlist(l = list(dat_hap, dat_hap2, dat_hap3)) %>% unique() %>% setDT()

## Rename
setnames(x = dat_hap,
         old = c("df_i.haplogroup1", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
         new = c("Haplogroup", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

library(writexl)
write_xlsx(dat_hap, "SEA_haplogroups_updated.xlsx")

################ SUBCLADES #######################

# Write separate subclade sequences

# Haplogroup G
hap_G <- dat %>% filter(haplogroup1 == "G")
nbin_G <- nbin[labels(nbin) %in% hap_G$name]
G <- file2[hap_G$name]
writeXStringSet(G, "data/GAligned.fasta")

# Haplogroup M
hap_M <- dat %>% filter(haplogroup1 == "M")
nbin_M <- nbin[labels(nbin) %in% hap_M$name]
M <- file2[hap_M$name]
writeXStringSet(M, "data/MAligned.fasta")

# Haplogroup N
hap_N <- dat %>% filter(haplogroup1 == "N")
nbin_N <- nbin[labels(nbin) %in% hap_N$name]
N <- file2[hap_N$name]
writeXStringSet(N, "data/NAligned.fasta")

# Haplogroup P
hap_P <- dat %>% filter(haplogroup1 == "P")
nbin_P <- nbin[labels(nbin) %in% hap_P$name]
P <- file2[hap_P$name]
writeXStringSet(P, "data/PAligned.fasta")

# Haplogroup R
hap_R <- dat %>% filter(haplogroup1 == "R")
nbin_R <- nbin[labels(nbin) %in% hap_R$name]
R <- file2[hap_R$name]
writeXStringSet(R, "data/RAligned.fasta")

# Haplogroup F1
hap_F1 <- dat %>% filter(haplogroup2 == "F1")
nbin_F1 <- nbin[labels(nbin) %in% hap_F1$name]

F1 <- file2[hap_F1$name]
writeXStringSet(F1, "data/F1Aligned.fasta")

# Haplogroup F1(xF1a)

hap_F1xF1a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup2 == "F1" & haplogroup4 != "F1a")
nbin_F1xF1a <- nbin[labels(nbin) %in% hap_F1xF1a$name]

F1xF1a <- file2[hap_F1xF1a$name]
writeXStringSet(F1xF1a, "data/F1xF1aAligned.fasta")

# # Haplogroup F1a
# hap_F1a <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "F1a")
# nbin_F1a <- nbin[labels(nbin) %in% hap_F1a$name]
# 
# F1a <- file[hap_F1a$name]
# writeXStringSet(F1a, "data/F1a.fasta")
# 
# # cat(file="data/F1a.fasta", paste(paste0(">",names(nbin_F1a)),
# #                                  sapply(nbin_F1a, paste, collapse=""), sep="\n"), sep="\n")
# 

# Haplogroup F2
hap_F2 <- dat %>% filter(haplogroup2 == "F2")
nbin_F2 <- nbin[labels(nbin) %in% hap_F2$name]

F2 <- file2[hap_F2$name]
writeXStringSet(F2, "data/F2Aligned.fasta")

# Haplogroup F3
hap_F3 <- dat %>% filter(haplogroup2 == "F3")
nbin_F3 <- nbin[labels(nbin) %in% hap_F3$name]

F3 <- file2[hap_F3$name]
writeXStringSet(F3, "data/F3Aligned.fasta")

# Haplogroup B5
hap_B5 <- dat %>% filter(haplogroup2 == "B5")
nbin_B5 <- nbin[labels(nbin) %in% hap_B5$name]

B5 <- file2[hap_B5$name]
writeXStringSet(B5, "data/B5Aligned.fasta")

# # Haplogroup B5a
# 
# hap_B5a <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="B5", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "B5a")
# 
# B5a <- file[hap_B5a$name]
# writeXStringSet(B5a, "data/B5a.fasta")
# 

# Haplogroup M7
hap_M7 <- dat %>% filter(haplogroup2 == "M7")
nbin_M7 <- nbin[labels(nbin) %in% hap_M7$name]

M7 <- file2[hap_M7$name]
writeXStringSet(M7, "data/M7Aligned.fasta")

# # Haplogroup M7b
# hap_M7b <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "M7b")
# 
# M7b <- file[hap_M7b$name]
# writeXStringSet(M7b, "data/M7b.fasta")
# 
# hap_M7b1a1 <- dat %>%
#   mutate(haplogroup4 = ifelse(haplogroup2=="M7", str_extract(haplo, "^([A-Z])\\d+\\w\\d\\w\\d"), haplogroup3)) %>%
#   filter(haplogroup4 == "M7b1a1")
# 
# M7b1a1 <- file[hap_M7b1a1$name]
# writeXStringSet(M7b1a1, "data/M7b1a1.fasta")
# 

# Haplogroup B4
hap_B4 <- dat %>% filter(haplogroup2 == "B4")
nbin_B4 <- nbin[labels(nbin) %in% hap_B4$name]

B4 <- file2[hap_B4$name]
writeXStringSet(B4, "data/B4Aligned.fasta")

# # Haplogroup B4a
# 
# hap_B4a <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "B4a")
# 
# B4a <- file[hap_B4a$name]
# writeXStringSet(B4a, "data/B4a.fasta")
# 
# # Haplogroup M7c
# hap_M7c <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "M7c")
# 
# M7c <- file[hap_M7c$name]
# writeXStringSet(M7c, "data/M7c.fasta")
# 
# # Haplogroup B4c
# 
# hap_B4c <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "B4c")
# 
# B4c <- file[hap_B4c$name]
# writeXStringSet(B4c, "data/B4c.fasta")
# 
# # Haplogroup F1f
# hap_F1f <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "F1f")
# 
# F1f <- file[hap_F1f$name]
# writeXStringSet(F1f, "data/F1f.fasta")
# 

# Haplogroup N9
hap_N9 <- dat %>% filter(haplogroup2 == "N9")
nbin_N9 <- nbin[labels(nbin) %in% hap_N9$name]

N9 <- file2[hap_N9$name]
writeXStringSet(N9, "data/N9Aligned.fasta")

# # Haplogroup N9a
# hap_N9a <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="N9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "N9a")
# 
# N9a <- file[hap_N9a$name]
# writeXStringSet(N9a, "data/N9a.fasta")
# 

# Haplogroup R9
hap_R9 <- dat %>% filter(haplogroup2 == "R9")
nbin_R9 <- nbin[labels(nbin) %in% hap_R9$name]

R9 <- file2[hap_R9$name]
writeXStringSet(R9, "data/R9Aligned.fasta")

# # Haplogroup R9b
# hap_R9b <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="R9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "R9b")
# 
# R9b <- file[hap_R9b$name]
# writeXStringSet(R9b, "data/R9b.fasta")
# 

# # Haplogroup M74
# hap_M74 <- dat %>%filter(haplogroup3 == "M74")
# 
# M74 <- file[hap_M74$name]
# writeXStringSet(M74, "data/M74.fasta")
# 
# Haplogroup C7
hap_C7 <- dat %>% filter(haplogroup2 == "C7")
nbin_C7 <- nbin[labels(nbin) %in% hap_C7$name]

C7 <- file2[hap_C7$name]
writeXStringSet(C7, "data/C7Aligned.fasta")

# # Haplogroup C7a
# hap_C7a <- dat %>% 
#   mutate(haplogroup4 = ifelse(haplogroup3=="C7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
#   filter(haplogroup4 == "C7a")
# 
# C7a <- file[hap_C7a$name]
# writeXStringSet(C7a, "data/C7a.fasta")

# Haplogroup D4
hap_D4 <- dat %>% filter(haplogroup2 == "D4")
nbin_D4 <- nbin[labels(nbin) %in% hap_D4$name]

D4 <- file2[hap_D4$name]
writeXStringSet(D4, "data/D4Aligned.fasta")

###########################################

names <- as.data.frame(nbinmat)
data <- separate(names, V1, into = c("country", "ethnic", "access_no", "haplo"), sep = "\\.")
data <- cbind(data, names) %>% dplyr::rename(name=V1)

dat <- data %>% mutate(
  haplogroup1 = substr(haplo, 1, 1),
  haplogroup2 = substr(haplo, 1, 2),
  haplogroup3 = str_extract(haplo, "^([A-Z])\\d+"),
  haplogroup3 = ifelse(is.na(haplogroup3), haplo, haplogroup3)
) %>% setDT()

library(writexl)
# dat_write <- data %>% mutate(
#   haplogroup1 = substr(haplo, 1, 1),
#   haplogroup2 = str_extract(haplo, "^([A-Z])\\d+"),
#   haplogroup2 = ifelse(is.na(haplogroup2), haplo, haplogroup2),
#   haplogroup3=ifelse(!(haplo %in% c("A+152", "A+152+16362", " A+152+16362+200", "R+16189")), str_extract(haplo, "^([A-Z])\\d\\w"), haplo),
#   haplogroup3=ifelse(is.na(haplogroup3), haplo, haplogroup3),
#   haplogroup3=case_when(haplogroup3 %in% c("134", "157", "161", "168", "171", "174", "241", "257", "259", "260", "262", "268", "279", "281", "295", "30", "304", 
#                                  "313", "315", "32", "351", "375", "381", "425", "433", "439", "483", "489", "496", "50", "501", "51", "53", "530", 
#                                  "563", "564", "567", "589", "607", "62", "73", "78", "90") ~ haplo,
#                         TRUE ~ haplogroup3),
#   ID=order(name)
#   ) %>% setDT()
# dat_write <- dat_write[, c(9,5,1,2,3,4,6,7,8)]
# write_xlsx(dat_write, "SEA_haplogroups.xlsx")

hap_country <- dat[, .N, by = .(haplo, country)]
hap1_country <- dat[, .N, by = .(haplogroup1, country)]

country_hap <- dat[, .N, by = .(country, haplo)]
country_hap <- country_hap %>%
  group_by(country) %>% arrange(haplo, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

country_hap1 <- dat[, .N, by = .(country, haplogroup1)]
country_hap1 <- country_hap1 %>%
  group_by(country) %>% arrange(haplogroup1, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

country_hap3 <- dat[, .N, by = .(country, haplogroup3)]
country_hap3 <- country_hap3 %>%
  group_by(country) %>% arrange(haplogroup3, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

ethnic <- read.csv("SEAFreq_by_label.csv")
ethnic <- ethnic %>%
  mutate(haplo=str_extract(Haplogroup, "^([A-Z])\\d+"),
         haplo=ifelse(is.na(haplo), Haplogroup, haplo))

library(reshape2)
ethnic_hap <- ethnic %>% select(-Haplogroup) %>%
  melt(., id.vars="haplo")

###################################################

# Country

dat <- read_excel("IsolateExplanation.xlsx")
dat <- as.data.frame(dat)
names <- dat$names

country <- dat %>%
  dplyr::select(name, Country, `Language family`, Ethnicity, haplo, haplogroup1, haplogroup2, haplogroup3) %>% setDT()
country_rank <- country[, .N, by = .(Country)] %>% arrange(desc(N))

country_hap <- country[, .N, by = .(Country, haplo)]
country_hap <- country_hap %>%
  group_by(Country) %>% arrange(haplo, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

country_hap3 <- country[, .N, by = .(Country, haplogroup3)]
country_hap3 <- country_hap3 %>%
  group_by(Country) %>% arrange(haplogroup3, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(janitor)

dat_c <- NULL
for (i in country_rank$Country) {
  cat(i, "\n")
  df_i <- country %>% dplyr::filter(Country==i)
  n_i <- nrow(df_i)
  nbin_i <- nbin[labels(nbin) %in% df_i$name]
  h_i <- pegas::haplotype(nbin_i)
  n_hi <- length(as.list(h_i))
  fas_i <- file2[labels(nbin) %in% df_i$name]
  writeXStringSet(fas_i, paste0("data/country/", i, "_Aligned.fasta"))
  dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
  x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
  dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
  mpd_i <- mean(dist.gene_i)
  hap.div_i <- hap.div(nbin_i, variance = TRUE)
  nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
  dati <- data.frame(df_i$Country, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
  ## combine
  dat_c <- rbindlist(l = list(dat_c, dati)) %>% unique() %>% setDT()
}
# Rename
setnames(x = dat_c,
         old = c("df_i.Country", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
         new = c("Country", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

library(writexl)
write_xlsx(dat_c, "SEA_country.xlsx")

# Ethnicity

ethnic <- read_excel("IsolateExplanation.xlsx")

library(janitor)
ethnic2 <- ethnic %>%
  mutate(Ethnicity=ifelse(Ethnicity=="Lao, Akha, Hmong, Khmu, Yao/Mien, Phuan", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan",
                          ifelse(Ethnicity=="Yao", "Dao", Ethnicity))) %>%
  dplyr::select(name, Country, `Language family`, Ethnicity, haplo, haplogroup1, haplogroup2, haplogroup3) %>% setDT()
ethnic_rank <- ethnic2[, .N, by = .(Ethnicity)] %>% arrange(desc(N))

ethnic_hap <- ethnic2[, .N, by = .(Ethnicity, haplo)]
ethnic_hap <- ethnic_hap %>%
  group_by(Ethnicity) %>% arrange(haplo, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

ethnic_hap3 <- ethnic2[, .N, by = .(Ethnicity, haplogroup3)]
ethnic_hap3 <- ethnic_hap3 %>%
  group_by(Ethnicity) %>% arrange(haplogroup3, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(pegas)
library(ape)

# dat_e <- NULL
# for (i in ethnic_rank$Ethnicity) {
#   cat(i, "\n")
#   df_i <- ethnic2 %>% dplyr::filter(Ethnicity==i)
#   n_i <- nrow(df_i)
#   nbin_i <- nbin[labels(nbin) %in% df_i$name]
#   h_i <- pegas::haplotype(nbin_i)
#   n_hi <- length(as.list(h_i))
#   fas_i <- file2[labels(nbin) %in% df_i$name]
#   writeXStringSet(fas_i, paste0("data/ethnic/", i, ".fasta"))
#   dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
#   x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
#   dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
#   mpd_i <- mean(dist.gene_i)
#   hap.div_i <- hap.div(nbin_i, variance = TRUE)
#   nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
#   dati <- data.frame(df_i$Ethnicity, df_i$`Language family`, df_i$Country, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
#   ## combine
#   dat_e <- rbindlist(l = list(dat_e, dati)) %>% unique() %>% setDT()
# }
# # Rename
# setnames(x = dat_e,
#          old = c("df_i.Ethnicity", "df_i..Language.family.", "df_i.Country", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
#          new = c("Ethnic", "Language", "Country", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))
# 
# library(writexl)
# write_xlsx(dat_e, "SEA_ethnicity.xlsx")

dat_e <- read_excel("SEA_ethnicity.xlsx")

# df <- dat_e %>% na.omit() %>% dplyr::select(1,4,6,8)
df <- dat_e %>%
  mutate(`Language`=ifelse(`Language`=="Austronesian, Austroasiatic" | `Language`=="Austroasiatic, Austronesian", "Austroasiatic, Austronesian, Sino-Tibetan",
                                  ifelse(`Language`=="Austronesian, Spanish" | `Language`=="Austronesian, Trans-New Guinea" | `Language`=="Papuan, Austronesian, English", "Austronesian",
                                         ifelse(`Language`=="Hmong-Mien" | `Language`=="Hmong-Mien, Mongolic", "Hmong-Mien +",
                                                ifelse(`Language`=="Indo-European, Sino-Tibetan" | `Language`=="Sino-Tibetan, Austroasiatic, Tai-Kadai" | `Language`=="Sino-Tibetan, Tai-Kada" | `Language`=="Tai-Kadai, Sino-Tibetan", "Sino-Tibetan +",
                                                       ifelse(`Language`=="Trans-New Guinea" | `Language`=="Trans–New Guinea (Alor-Pantar, Papuan)", "Trans-New Guinea +",
                                                              ifelse(`Language`=="Austronesian, Austroasiatic, Indo-European", "Austroasiatic, Austronesian + (xSino-Tibetan)",
                                                                     ifelse(`Language`=="Austroasiatic, Tai-Kadai" | `Language`=="Austroasiatic, Tai-Kadai, Hmong-Mien, Sino-Tibetan" | `Language`=="Tai-Kadai, Hmong-Mien, Austroasiatic", "Austroasiatic, Tai-Kadai +",
                                                                            ifelse(`Language`=="Sino-Tibetan", "Sino-Tibetan +", `Language`)))))))),
         `Language`=ifelse(Ethnic=="Mon", "Austroasiatic",
                           ifelse(Ethnic=="Hmong", "Hmong-Mien +",
                                  ifelse(Ethnic=="Shan", "Tai-Kadai",
                                         ifelse(Ethnic=="Jehai (or Jahai)", "Austroasiatic",
                                                ifelse(Ethnic=="Temuan", "Austronesian",
                                                       ifelse(Ethnic=="Maranao", "Austronesian",
                                                              ifelse(Ethnic=="Semelai", "Austroasiatic",
                                                                     ifelse(Ethnic=="Bru (Brao)", "Austroasiatic",
                                                                            ifelse(Ethnic=="Jarai", "Austronesian",
                                                                                   ifelse(Ethnic=="Kadazan-Dusun", "Austronesian",
                                                                                          ifelse(Ethnic=="Alor", "Austronesian",
                                                                                                 ifelse(Ethnic=="Arakanese (or Rakhine)", "Sino-Tibetan +",
                                                                                                        ifelse(Ethnic=="Timorese", "Austronesian",
                                                                                                               ifelse(Ethnic=="Mang", "Austroasiatic", `Language`))))))))))))))) %>%
  na.omit() %>% dplyr::select(1,2,3,4,6,8,10) %>% setDF()
df <- df %>% filter(`Sample size`>2 & !`Ethnic` %in% c("Unknown", "Khmer, Cham, Chinese-Cambodian, Vietnamese", "Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Banjar, Bantenese, Banyumasan", "Batak, Minangkabau, Acehnese, Lampung", "Banjar, Dayak, Javanese"))
# dist_matrix <- dist(dat_e[,-1])
# mds_result <- cmdscale(dist_matrix)
# plot(mds_result, col = dat_e$ethnic, pch = 19, xlab = "MDS1", ylab = "MDS2")

library(MASS)
library(magrittr)
library(dplyr)
library(ggpubr)

# Compute MDS
mds <- df %>% na.omit() %>% dplyr::select(-c(1,2,3,4)) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          size = 1,
          repel = TRUE)

# K-means clustering (K=10)
clust <- kmeans(mds, 10)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K10.png"), width = 15, height = 10)

# K-means clustering (K=10) - Language
mds <- df %>% na.omit() %>% dplyr::select(-c(1,2,3,4)) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

clust_lang <- kmeans(mds, 10)$cluster %>%
  as.factor()
mds_lang <- mds %>%
  mutate(groups = clust, language = df$Language)
# Plot and color by groups
ggscatter(mds_lang, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "language",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K10_language.png"), width = 15, height = 10)

# K-means clustering (K=6)
clust <- kmeans(mds, 6)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K6.png"), width = 15, height = 10)

library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
set.seed(123)
df2 <- tibble::column_to_rownames(df, var = "Ethnic") %>% dplyr::select(-c("Language", "Country", "Sample size")) %>% as.data.frame()
distance <- get_dist(df2)

fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
k2 <- kmeans(df2, centers = 2, nstart = 25)
fviz_cluster(k2, data = df2)

df2 %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         ethnic = df$Ethnic) %>%
  ggplot(aes(`Haplotype diversity (H)`, `Nucleotide diveristy (pi)`, color = factor(cluster), label = ethnic)) +
  geom_text()

k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25)
k6 <- kmeans(df2, centers = 6, nstart = 25)
k7 <- kmeans(df2, centers = 7, nstart = 25)
k8 <- kmeans(df2, centers = 8, nstart = 25)
k9 <- kmeans(df2, centers = 9, nstart = 25)
k10 <- kmeans(df2, centers = 10, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5")
p5 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 6")
p6 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 7")
p7 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 8")
p8 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 9")
p9 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 10")

library(gridExtra)
png(filename = file.path("figures", "ethnic_K2-5.png"), width = 1200, height = 800)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

png(filename = file.path("figures", "ethnic_K6-10.png"), width = 1200, height = 800)
grid.arrange(p5, p6, p7, p8, p9, nrow = 2)
dev.off()

fviz_cluster(k5, data = df2)
ggsave(filename = file.path("figures", "ethnic_K5.png"), width = 15, height = 10)

df2 %>%
  as_tibble() %>%
  mutate(cluster = k5$cluster,
         ethnic = df$Ethnic) %>%
  ggplot(aes(`Haplotype diversity (H)`, `Nucleotide diveristy (pi)`, color = factor(cluster), label = ethnic)) +
  geom_text()

fviz_nbclust(df2, kmeans, method = "wss")
fviz_nbclust(df2, kmeans, method = "silhouette")

# Compute k-means clustering with k = 5
set.seed(123)
final <- kmeans(df2, 5, nstart = 25)
print(final)

fviz_cluster(final, data = df2)
ggsave(filename = file.path("figures", "ethnic_K5_final.png"), width = 15, height = 10)

# Ethnicity - FULL

dat_e <- read_excel("SEA_ethnicity.xlsx")

# ethnic <- read_excel("IsolateExplanation.xlsx")
# 
# library(janitor)
# ethnic2 <- ethnic %>%
#   mutate(Ethnicity=ifelse(Ethnicity=="Lao, Akha, Hmong, Khmu, Yao/Mien, Phuan", "Lao, Akha, Hmong, Khmu, Yao, Mien, Phuan", Ethnicity)) %>%
#   dplyr::select(name, Country, `Language family`, Ethnicity, haplo, haplogroup1, haplogroup2, haplogroup3) %>% setDT()
# ethnic_rank <- ethnic2[, .N, by = .(Ethnicity)] %>% arrange(desc(N))
# 
# ethnic_hap <- ethnic2[, .N, by = .(Ethnicity, haplo)]
# ethnic_hap <- ethnic_hap %>%
#   group_by(Ethnicity) %>% arrange(haplo, .by_group = TRUE) %>% 
#   mutate(percent=(N*100)/sum(N)) %>% ungroup()
# 
# ethnic_hap3 <- ethnic2[, .N, by = .(Ethnicity, haplogroup3)]
# ethnic_hap3 <- ethnic_hap3 %>%
#   group_by(Ethnicity) %>% arrange(haplogroup3, .by_group = TRUE) %>% 
#   mutate(percent=(N*100)/sum(N)) %>% ungroup()
# 
# library(pegas)
# library(ape)
# 
# dat_e <- NULL
# for (i in ethnic_rank$Ethnicity) {
#   cat(i, "\n")
#   df_i <- ethnic2 %>% dplyr::filter(Ethnicity==i)
#   n_i <- nrow(df_i)
#   nbin_i <- nbin[labels(nbin) %in% df_i$name]
#   h_i <- pegas::haplotype(nbin_i)
#   n_hi <- length(as.list(h_i))
#   # fas_i <- file[labels(nbin) %in% df_i$name]
#   # writeXStringSet(fas_i, paste0("data/ethnic/", i, ".fasta"))
#   # dnbin_i <- dist.dna(nbin_i, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
#   x_i <- as.matrix.DNAbin(nbin_i)  #converting DNAbin to matrix
#   dist.gene_i <- dist.gene(x_i, method = "pairwise", variance = TRUE, pairwise.deletion = FALSE)
#   mpd_i <- mean(dist.gene_i)
#   hap.div_i <- hap.div(nbin_i, variance = TRUE)
#   nuc.div_i <- nuc.div(nbin_i, variance = TRUE, pairwise.deletion = FALSE) # hap.div = 1 all haplotypes are unique
#   dati <- data.frame(df_i$Ethnicity, df_i$`Language family`, df_i$Country, n_i, n_hi, hap.div_i[1], hap.div_i[2], nuc.div_i[1], nuc.div_i[2], mpd_i) %>% slice(1)
#   ## combine
#   dat_e <- rbindlist(l = list(dat_e, dati)) %>% unique() %>% setDT()
# }
# ## Rename
# setnames(x = dat_e,
#          old = c("df_i.Ethnicity", "df_i..Language.family.", "df_i.Country", "n_i", "n_hi", "hap.div_i.1.", "hap.div_i.2.", "nuc.div_i.1.", "nuc.div_i.2.", "mpd_i"),
#          new = c("Ethnic", "Language", "Country", "Sample size", "Number of haplotypes", "Haplotype diversity (H)", "H variance", "Nucleotide diveristy (pi)", "pi SE", "MPD"))

# df <- dat_e %>% na.omit() %>% dplyr::select(1,4,6,8)
df <- dat_e %>%
  mutate(`Language`=ifelse(`Language`=="Austronesian, Austroasiatic" | `Language`=="Austroasiatic, Austronesian", "Austroasiatic, Austronesian, Sino-Tibetan",
                           ifelse(`Language`=="Austronesian, Spanish" | `Language`=="Austronesian, Trans-New Guinea" | `Language`=="Papuan, Austronesian, English", "Austronesian",
                                  ifelse(`Language`=="Hmong-Mien" | `Language`=="Hmong-Mien, Mongolic", "Hmong-Mien +",
                                         ifelse(`Language`=="Indo-European, Sino-Tibetan" | `Language`=="Sino-Tibetan, Austroasiatic, Tai-Kadai" | `Language`=="Sino-Tibetan, Tai-Kada" | `Language`=="Tai-Kadai, Sino-Tibetan", "Sino-Tibetan +",
                                                ifelse(`Language`=="Trans-New Guinea" | `Language`=="Trans–New Guinea (Alor-Pantar, Papuan)", "Trans-New Guinea +",
                                                       ifelse(`Language`=="Austronesian, Austroasiatic, Indo-European", "Austroasiatic, Austronesian + (xSino-Tibetan)",
                                                              ifelse(`Language`=="Austroasiatic, Tai-Kadai" | `Language`=="Austroasiatic, Tai-Kadai, Hmong-Mien, Sino-Tibetan" | `Language`=="Tai-Kadai, Hmong-Mien, Austroasiatic", "Austroasiatic, Tai-Kadai +",
                                                                     ifelse(`Language`=="Sino-Tibetan", "Sino-Tibetan +", `Language`)))))))),
         `Language`=ifelse(Ethnic=="Mon", "Austroasiatic",
                           ifelse(Ethnic=="Hmong", "Hmong-Mien +",
                                  ifelse(Ethnic=="Shan", "Tai-Kadai",
                                         ifelse(Ethnic=="Jehai (or Jahai)", "Austroasiatic",
                                                ifelse(Ethnic=="Temuan", "Austronesian",
                                                       ifelse(Ethnic=="Maranao", "Austronesian",
                                                              ifelse(Ethnic=="Semelai", "Austroasiatic",
                                                                     ifelse(Ethnic=="Bru (Brao)", "Austroasiatic",
                                                                            ifelse(Ethnic=="Jarai", "Austronesian",
                                                                                   ifelse(Ethnic=="Kadazan-Dusun", "Austronesian",
                                                                                          ifelse(Ethnic=="Alor", "Austronesian",
                                                                                                 ifelse(Ethnic=="Arakanese (or Rakhine)", "Sino-Tibetan +",
                                                                                                        ifelse(Ethnic=="Timorese", "Austronesian",
                                                                                                               ifelse(Ethnic=="Mang", "Austroasiatic", `Language`))))))))))))))) %>%
  na.omit() %>% dplyr::select(1,2,3,4,6,8,10) %>% setDF()
df <- df %>% filter(`Ethnic`!="Unknown")
# dist_matrix <- dist(dat_e[,-1])
# mds_result <- cmdscale(dist_matrix)
# plot(mds_result, col = dat_e$ethnic, pch = 19, xlab = "MDS1", ylab = "MDS2")

# Compute MDS
mds <- df %>% na.omit() %>% dplyr::select(-c(1,2,3,4)) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          size = 1,
          repel = TRUE)

# K-means clustering (K=10)
clust <- kmeans(mds, 10)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K10_full.png"), width = 15, height = 10)

# K-means clustering (K=10) - Language
mds <- df %>% na.omit() %>% dplyr::select(-c(1,2,3,4)) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

clust_lang <- kmeans(mds, 10)$cluster %>%
  as.factor()
mds_lang <- mds %>%
  mutate(groups = clust, language = df$Language)
# Plot and color by groups
ggscatter(mds_lang, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "language",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K10_language_full.png"), width = 15, height = 10)

# K-means clustering (K=6)
clust <- kmeans(mds, 6)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$Ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K6_full.png"), width = 15, height = 10)

library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
set.seed(123)
df2 <- tibble::column_to_rownames(df, var = "Ethnic") %>% dplyr::select(-c("Language", "Country", "Sample size")) %>% as.data.frame()
distance <- get_dist(df2)

fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
k2 <- kmeans(df2, centers = 2, nstart = 25)
fviz_cluster(k2, data = df2)

df2 %>%
  as_tibble() %>%
  mutate(cluster = k2$cluster,
         ethnic = df$Ethnic) %>%
  ggplot(aes(`Haplotype diversity (H)`, `Nucleotide diveristy (pi)`, color = factor(cluster), label = ethnic)) +
  geom_text()

k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25)
k6 <- kmeans(df2, centers = 6, nstart = 25)
k7 <- kmeans(df2, centers = 7, nstart = 25)
k8 <- kmeans(df2, centers = 8, nstart = 25)
k9 <- kmeans(df2, centers = 9, nstart = 25)
k10 <- kmeans(df2, centers = 10, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5")
p5 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 6")
p6 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 7")
p7 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 8")
p8 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 9")
p9 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 10")

library(gridExtra)
png(filename = file.path("figures", "ethnic_K2-5_full.png"), width = 1200, height = 800)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

png(filename = file.path("figures", "ethnic_K6-10_full.png"), width = 1200, height = 800)
grid.arrange(p5, p6, p7, p8, p9, nrow = 2)
dev.off()

fviz_cluster(k5, data = df2)
ggsave(filename = file.path("figures", "ethnic_K5_full.png"), width = 15, height = 10)

df2 %>%
  as_tibble() %>%
  mutate(cluster = k5$cluster,
         ethnic = df$Ethnic) %>%
  ggplot(aes(`Haplotype diversity (H)`, `Nucleotide diveristy (pi)`, color = factor(cluster), label = ethnic)) +
  geom_text()

fviz_nbclust(df2, kmeans, method = "wss")
fviz_nbclust(df2, kmeans, method = "silhouette")

# Compute k-means clustering with k = 5
set.seed(123)
final <- kmeans(df2, 5, nstart = 25)
print(final)

fviz_cluster(final, data = df2)
ggsave(filename = file.path("figures", "ethnic_K5_final_full.png"), width = 15, height = 10)

## Cham

Cham <- VN %>% filter(ethnic=="Cham")
nbin_Cham <- nbin[labels(nbin) %in% Cham$name]
class(nbin_Cham)
dnbin_Cham <- dist.dna(nbin_Cham, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
hap.div(nbin_Cham, variance = TRUE)
nuc.div(nbin_Cham, variance = TRUE, pairwise.deletion = FALSE)
x_Cham <- as.matrix.DNAbin(nbin_Cham)  #converting DNAbin to matrix

disgene_Cham <- dist.gene(x_Cham, method = "pairwise", pairwise.deletion = FALSE, variance = TRUE)

h_Cham <- pegas::haplotype(nbin_Cham)
d_Cham <- dist.dna(h_Cham, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
nt_Cham <- rmst(d_Cham, quiet = TRUE)#constructs the haplotype network
sz_Cham <- summary(h_Cham)
nt.labs_Cham <- attr(nt_Cham, "labels")
sz_Cham <- sz_VN[nt.labs_Cham]

name_Cham <- Cham$name
R <- haploFreq(x_Cham, fac = name_Cham, haplo = h_Cham)
R <- R[nt.labs_Cham, ]
samp<-R %>% t()
phylo_Cham <- ape::as.phylo(nt_Cham)
dis <- cophenetic(phylo_Cham)

dis<-dis[,phylo_Cham$tip.label]
samp<-samp[,phylo_Cham$tip.label]

comdist(samp, dis, abundance.weighted=TRUE)
pd(samp, phylo_Cham, include.root = TRUE)
mpd(samp, dis)
mpdn(samp, dis, abundance.weighted = TRUE, time.output = FALSE)

tree_Cham<-nj(dnbin_Cham)
ggt_Cham<-ggtree::ggtree(tree_Cham, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(size=3)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_Cham

### MDS based on haplogroup proportion (frequency)

library(reshape2)
ehap <- ethnic_hap3 %>%
  dplyr::rename(ethnic=Ethnicity) %>%
  group_by(haplogroup3) %>%
  mutate(ID=order(ethnic)) %>%
  ungroup() %>%
  dplyr::rename(key=haplogroup3, value=N) %>%
  dplyr::select(-percent) %>%
  arrange(key) %>% filter(!ethnic %in% c("Unknown", "Khmer, Cham, Chinese-Cambodian, Vietnamese", "Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung", "Lao, Akha, Hmong, Khmu, Dao, Mien, Phuan", "Lao, Tai Dam, Tai Deng, Tai Yuan, Katang, Phuan", "Banjar, Bantenese, Banyumasan", "Batak, Minangkabau, Acehnese, Lampung", "Banjar, Dayak, Javanese"))

dt_e <- spread(ehap, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_e %>% dplyr::select(-c(ethnic, ID))
m<-as.matrix(DT)
ethnic <- dt_e$ethnic
ID <- dt_e$ID
dt <- aggregate(m, data.frame(ethnic),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
ethnic_hap_all <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100))

library(MASS)
library(magrittr)
library(dplyr)
library(ggpubr)

dist_matrix <- dist(ethnic_hap_all[,-c("ethnic")])
mds_result <- cmdscale(dist_matrix)

# Compute MDS
mds <- ethnic_hap_all %>% dplyr::select(-1) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")
# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = ethnic_hap_all$ethnic,
          size = 1,
          repel = TRUE)

# Compute MDS
df <- ethnic_hap_all
mds <- df %>% na.omit() %>% dplyr::select(-1) %>%
  dist() %>%          
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c("Dim.1", "Dim.2")

# Plot MDS
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$ethnic,
          size = 1,
          repel = TRUE)

# K-means clustering (K=10)
clust <- kmeans(mds, 10)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K10_pro.png"), width = 15, height = 10)

# K-means clustering (K=6)
clust <- kmeans(mds, 6)$cluster %>%
  as.factor()
mds <- mds %>%
  mutate(groups = clust)
# Plot and color by groups
ggscatter(mds, x = "Dim.1", y = "Dim.2", 
          label = df$ethnic,
          color = "groups",
          palette = "jco",
          size = 1, 
          ellipse = TRUE,
          ellipse.type = "convex",
          repel = TRUE)
ggsave(filename = file.path("figures", "ethnic_K6_pro.png"), width = 15, height = 10)

library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
set.seed(123)
df2 <- tibble::column_to_rownames(df, var = "ethnic") %>% as.data.frame()
distance <- get_dist(df2)

fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
k2 <- kmeans(df2, centers = 2, nstart = 25)
fviz_cluster(k2, data = df2)

k3 <- kmeans(df2, centers = 3, nstart = 25)
k4 <- kmeans(df2, centers = 4, nstart = 25)
k5 <- kmeans(df2, centers = 5, nstart = 25)
k6 <- kmeans(df2, centers = 6, nstart = 25)
k7 <- kmeans(df2, centers = 7, nstart = 25)
k8 <- kmeans(df2, centers = 8, nstart = 25)
k9 <- kmeans(df2, centers = 9, nstart = 25)
k10 <- kmeans(df2, centers = 10, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k2, geom = "point", data = df2) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = df2) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = df2) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 5")
p5 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 6")
p6 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 7")
p7 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 8")
p8 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 9")
p9 <- fviz_cluster(k5, geom = "point",  data = df2) + ggtitle("k = 10")

library(gridExtra)
png(filename = file.path("figures", "ethnic_K2-5_pro.png"), width = 1200, height = 800)
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

png(filename = file.path("figures", "ethnic_K6-10_pro.png"), width = 1200, height = 800)
grid.arrange(p5, p6, p7, p8, p9, nrow = 2)
dev.off()

fviz_cluster(k5, data = df2)
ggsave(filename = file.path("figures", "ethnic_K5_pro.png"), width = 15, height = 10)

fviz_nbclust(df2, kmeans, method = "wss")
fviz_nbclust(df2, kmeans, method = "silhouette")

library(viridis)
g1 <- ggplot(country_hap) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplo)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=50, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2, "cm")) +
  coord_flip()
g1
ggsave(filename = file.path("figures", "country_haplo.png"), width = 15, height = 10)

g2 <- ggplot(country_hap1) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplogroup1)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g2
ggsave(filename = file.path("figures", "country_haplo1.png"), width = 15, height = 10)

g3 <- ggplot(country_hap3) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplogroup3)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g3
ggsave(filename = file.path("figures", "country_haplo3.png"), width = 15, height = 10)

g4 <- ggplot(ethnic_hap) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(variable), y=value, fill=factor(haplo)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=5, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Ethinicity") +
  scale_y_continuous(name = "Percent") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  coord_polar()
g4
ggsave(filename = file.path("figures", "ethnic_haplo.png"), width = 20, height = 15)

options(ignore.negative.edge=TRUE)

# Vietnam

VN <- dat %>% filter(country=="Vietnam") %>%
  mutate(ethnic = ifelse(str_detect(ethnic, "Cham") == TRUE, "Cham",
                         ifelse(str_detect(ethnic, "CoLao") == TRUE, "CoLao",
                                ifelse(str_detect(ethnic, "Dao") == TRUE, "Dao",
                                       ifelse(str_detect(ethnic, "DHX") == TRUE, "Unspecified",
                                              ifelse(str_detect(ethnic, "DKK") == TRUE, "Unspecified",
                                                     ifelse(str_detect(ethnic, "DKX") == TRUE, "Unspecified",
                                                            ifelse(str_detect(ethnic, "DOX") == TRUE, "Unspecified",
                                                                   ifelse(str_detect(ethnic, "Ede") == TRUE, "Ede",
                                                                          ifelse(str_detect(ethnic, "Giarai") == TRUE, "Giarai",
                                                                                 ifelse(str_detect(ethnic, "HaNhi") == TRUE, "HaNhi",
                                                                                        ifelse(str_detect(ethnic, "HMong") == TRUE, "HMong",
                                                                                               ifelse(str_detect(ethnic, "Hui") == TRUE, "Unspecified",
                                                                                                      ifelse(str_detect(ethnic, "KDH") == TRUE, "Unspecified",
                                                                                                             ifelse(str_detect(ethnic, "Kinh") == TRUE, "Kinh",
                                                                                                                    ifelse(str_detect(ethnic, "LaChi") == TRUE, "LaChi",
                                                                                                                           ifelse(str_detect(ethnic, "LaHu") == TRUE, "LaHu",
                                                                                                                                  ifelse(str_detect(ethnic, "LoLo") == TRUE, "LoLo",
                                                                                                                                         ifelse(str_detect(ethnic, "Mang") == TRUE, "Mang",
                                                                                                                                                ifelse(str_detect(ethnic, "Nung") == TRUE, "Nung",
                                                                                                                                                       ifelse(str_detect(ethnic, "PaThen") == TRUE, "PaThen",
                                                                                                                                                              ifelse(str_detect(ethnic, "PhuLa") == TRUE, "PhuLa",
                                                                                                                                                                     ifelse(str_detect(ethnic, "SiLa") == TRUE, "SiLa",
                                                                                                                                                                            ifelse(str_detect(ethnic, "Tay") == TRUE, "Tay",
                                                                                                                                                                                   ifelse(str_detect(ethnic, "Thai") == TRUE, "Thai",
                                                                                                                                                                                          ifelse(str_detect(ethnic, "V206") == TRUE, "Viet",
                                                                                                                                                                                                 ifelse(str_detect(ethnic, "VIET") == TRUE, "Viet",
                                                                                                                                                                                                        ifelse(str_detect(ethnic, "VN") == TRUE, "Viet",
                                                                                                                                                                                                               ifelse(str_detect(ethnic, "Yao") == TRUE, "Yao",
                                                                                                                                                                                                                      ethnic)))))))))))))))))))))))))))))

nbin_VN <- nbin[labels(nbin) %in% VN$name]
class(nbin_VN)
dnbin_VN <- dist.dna(nbin_VN, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
library(pegas)
hap.div(nbin_VN, variance = TRUE)
nuc.div(nbin_VN, variance = TRUE, pairwise.deletion = FALSE)
library(picante)
x_VN <- as.matrix.DNAbin(nbin_VN)  #converting DNAbin to matrix

library(haplotypes)
h_VN <- pegas::haplotype(nbin_VN)
d_VN <- dist.dna(h_VN, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
nt_VN <- rmst(d_VN, quiet = TRUE)#constructs the haplotype network
sz_VN <- summary(h_VN)
nt.labs_VN <- attr(nt_VN, "labels")
sz_VN <- sz_VN[nt.labs_VN]

library(ape)
x_VN <- as.matrix.DNAbin(nbin_VN)
name_VN <- VN$name
R <- haploFreq(x_VN, fac = name_VN, haplo = h_VN)
R <- R[nt.labs_VN, ]
samp<-R %>% t()
library(ggmuller)
phylo_VN <- as.phylo(nt_VN)
dis <- cophenetic(phylo_VN)

dis<-dis[,phylo_VN$tip.label]
samp<-samp[,phylo_VN$tip.label]

library(picante) 
library(iCAMP)
pd(samp, phylo_VN, include.root = TRUE)
mpd(samp, dis, abundance.weighted=TRUE)
mpdn(R, cophenetic(phylo_VN), abundance.weighted = TRUE, time.output = FALSE)

tree_VN<-nj(dnbin_VN)
library(ggtree)
ggt_VN<-ggtree::ggtree(tree_VN, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_VN

## Cham

Cham <- VN %>% filter(ethnic=="Cham")
nbin_Cham <- nbin[labels(nbin) %in% Cham$name]
class(nbin_Cham)
dnbin_Cham <- dist.dna(nbin_Cham, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
hap.div(nbin_Cham, variance = TRUE)
nuc.div(nbin_Cham, variance = TRUE, pairwise.deletion = FALSE)
x_Cham <- as.matrix.DNAbin(nbin_Cham)  #converting DNAbin to matrix

h_Cham <- pegas::haplotype(nbin_Cham)
d_Cham <- dist.dna(h_Cham, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
nt_Cham <- rmst(d_Cham, quiet = TRUE)#constructs the haplotype network
sz_Cham <- summary(h_Cham)
nt.labs_Cham <- attr(nt_Cham, "labels")
sz_Cham <- sz_VN[nt.labs_Cham]

name_Cham <- Cham$name
R <- haploFreq(x_Cham, fac = name_Cham, haplo = h_Cham)
R <- R[nt.labs_Cham, ]
samp<-R %>% t()
phylo_Cham <- ape::as.phylo(nt_Cham)
dis <- cophenetic(phylo_Cham)

dis<-dis[,phylo_Cham$tip.label]
samp<-samp[,phylo_Cham$tip.label]

comdist(samp, dis, abundance.weighted=TRUE)
pd(samp, phylo_Cham, include.root = TRUE)
mpd(samp, dis)
mpdn(samp, dis, abundance.weighted = TRUE, time.output = FALSE)

tree_Cham<-nj(dnbin_Cham)
ggt_Cham<-ggtree::ggtree(tree_Cham, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(size=3)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_Cham

# Plotting Phylogenetic Distribution in Southeast Asia

library(sf)
library(readr)
library(proj4)
library(magick)
library(ggmap)
library(ggimage)
library(usethis)
library(scatterpie)
library(scales)
library(cowplot)

# library(sf)  # "plot"
# library(sp)  # "plot"
# # devtools::install_github("choisy/gadmSEA")
# # devtools::install_github("choisy/mcutils")
# library(mcutils)
# library(gadmSEA)
# library(magrittr)
# 
# ## defining colors:
# rgb2 <- function(...) rgb(..., max = 255)
# blue   <- rgb2(200, 237, 255)
# grey   <- rgb2(225, 225, 225)
# yellow <- rgb2(253, 252, 235)
# ## the map:
# plot(vietnam, xlab = "longitude (decimal degree)", border = NA,
#             ylab = "latitude (decimal degree)", bg = blue)
# mcutils::datasets("gadmSEA") %>%
# as.data.frame(stringsAsFactors = FALSE) %>%
# setNames("country") %>%
# mutate(col = ifelse(country == "vietnam", yellow, grey)) %$%
# invisible(purrr::map2(country, col, ~ plot(get(.x), col = .y, add = TRUE)))
# axis(1); axis(2); box(bty = "o")
# 
# # Plotting some countries in southeast Asia:
# library(sp)  # for "plot"
# library(gadmSEA) # for the countries
# library(magrittr) # for " %$% " and " %>% "
# 
# rgb2 <- function(...) rgb(..., max = 255)
# blue   <- rgb2(200, 237, 255)
# grey   <- rgb2(225, 225, 225)
# yellow <- rgb2(253, 252, 235)
# 
# ## the map:
# plot(vietnam, xlab = "longitude (decimal degree)", border = NA,
#               ylab = "latitude (decimal degree)", bg = blue,
#              xlim = c(68, 146), ylim = c(-11, 54))
# ctr <- c("vietnam", "philippines", "cambodia", "japan", "china", "thailand",
#        "singapore", "indonesia", "malaysia", "taiwan", "bangladesh", "laos",
#       "india", "nepal")
# mcutils::datasets("gadmSEA") %>%
# as.data.frame(stringsAsFactors = FALSE) %>%
# setNames("country") %>%
# dplyr::mutate(col = ifelse(country %in% ctr, yellow, grey)) %$%
# invisible(purrr::map2(country, col, ~ plot(get(.x), col = .y, add = TRUE)))
# axis(1); axis(2); box(bty = "o")
# 
# points(109.215971, 13.715266, col = "red", pch = 19)
# points(109.215971, 13.715266, col = "red", cex = 2)
# points(109.215971, 13.715266, col = "red", cex = 3)
# points(109.215971, 13.715266, col = "red", cex = 4)

# bangladesh_BGD0_sf <- readRDS("data/gadm36_BGD_0_sf.rds")
# bangladesh_BGD1_sf <- readRDS("data/gadm36_BGD_1_sf.rds")
# brunei_BRN0_sf <- readRDS("data/gadm36_BRN_0_sf.rds")
# brunei_BRN1_sf <- readRDS("data/gadm36_BRN_1_sf.rds")
# china_CHN0_sf <- readRDS("data/gadm36_CHN_0_sf.rds")
# china_CHN1_sf <- readRDS("data/gadm36_CHN_1_sf.rds")
# indonesia_IDN0_sf <- readRDS("data/gadm36_IDN_0_sf.rds")
# indonesia_IDN1_sf <- readRDS("data/gadm36_IDN_1_sf.rds")
# india_IND0_sf <- readRDS("data/gadm36_IND_0_sf.rds")
# india_IND1_sf <- readRDS("data/gadm36_IND_1_sf.rds")
# cambodia_KHM0_sf <- readRDS("data/gadm36_KHM_0_sf.rds")
# cambodia_KHM1_sf <- readRDS("data/gadm36_KHM_1_sf.rds")
# laos_LAO0_sf <- readRDS("data/gadm36_LAO_0_sf.rds")
# laos_LAO1_sf <- readRDS("data/gadm36_LAO_1_sf.rds")
# myanmar_MMR0_sf <- readRDS("data/gadm36_MMR_0_sf.rds")
# myanmar_MMR1_sf <- readRDS("data/gadm36_MMR_1_sf.rds")
# malaysia_MYS0_sf <- readRDS("data/gadm36_MYS_0_sf.rds")
# malaysia_MYS1_sf <- readRDS("data/gadm36_MYS_1_sf.rds")
# philippines_PHL0_sf <- readRDS("data/gadm36_PHL_0_sf.rds")
# philippines_PHL1_sf <- readRDS("data/gadm36_PHL_1_sf.rds")
# singapore_SGP0_sf <- readRDS("data/gadm36_SGP_0_sf.rds")
# singapore_SGP1_sf <- readRDS("data/gadm36_SGP_1_sf.rds")
# thailand_THA0_sf <- readRDS("data/gadm36_THA_0_sf.rds")
# thailand_THA1_sf <- readRDS("data/gadm36_THA_1_sf.rds")
# taiwan_TWN0_sf <- readRDS("data/gadm36_TWN_0_sf.rds")
# taiwan_TWN1_sf <- readRDS("data/gadm36_TWN_1_sf.rds")
# easttimor_TLS0_sf <- readRDS("data/gadm36_TLS_0_sf.rds")
# easttimor_TLS1_sf <- readRDS("data/gadm36_TLS_1_sf.rds")
# vietnam_VNM0_sf <- readRDS("data/gadm36_VNM_0_sf.rds")
# vietnam_VNM1_sf <- readRDS("data/gadm36_VNM_1_sf.rds")
# paracelislands_XPI0_sf <- readRDS("data/gadm36_XPI_0_sf.rds")
# spratlyislands_XSP0_sf <- readRDS("data/gadm36_XSP_0_sf.rds")

# SEA0_sf <- rbind(brunei_BRN0_sf, indonesia_IDN0_sf, cambodia_KHM0_sf, laos_LAO0_sf, myanmar_MMR0_sf, malaysia_MYS0_sf, philippines_PHL0_sf, singapore_SGP0_sf,
#                  thailand_THA0_sf, easttimor_TLS0_sf, vietnam_VNM0_sf, paracelislands_XPI0_sf, spratlyislands_XSP0_sf)
# 
# save(SEA0_sf, file = "data/SEA0_sf.RData")
# 
# SEA0p_sf <- rbind(bangladesh_BGD0_sf, brunei_BRN0_sf, china_CHN0_sf, indonesia_IDN0_sf, india_IND0_sf, cambodia_KHM0_sf, laos_LAO0_sf, myanmar_MMR0_sf, malaysia_MYS0_sf, philippines_PHL0_sf, singapore_SGP0_sf,
#                  thailand_THA0_sf, taiwan_TWN0_sf, easttimor_TLS0_sf, vietnam_VNM0_sf, paracelislands_XPI0_sf, spratlyislands_XSP0_sf)
# 
# save(SEA0p_sf, file = "data/SEA0p_sf.RData")

## Ancient mtDNA plus

load("data/SEA0p_sf.RData")
SEA0p_sf <- SEA0p_sf %>% dplyr::rename(country=NAME_0)

library(readxl)
ancient <- read_excel("all-ancient-dna-2-07-73-full.xlsx")
ancient_SEA <- ancient %>% filter(Country %in% c("Bangladesh", "Brunei", "Cambodia", "China", "Indonesia", "India", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Taiwan", "Timor-Leste", "Vietnam"))
dat_ancient_SEA <- ancient_SEA %>%
  rename(haplo="mtDNA-haplogroup",
         haploY="Y-Haplotree-Public",
         country="Country") %>%
  mutate(haplo1=ifelse(!(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189", "C4+152", "C4+152+16093", "D*", "E (95.07%)", "G1 (94.06%)", "M4″67", "n/a", "n/a (<2x)", "n/a (exome capture)", "R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)", "R2+13500", "U4'9", "W1+119")), str_extract(haplo, "^([A-Z])\\d\\w"), haplo),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo1),
         haplo1 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo1),
         haplo1 = ifelse(haplo %in% c("D*"), "D", haplo1),
         haplo1 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo1),
         haplo1 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo1),
         haplo1 = ifelse(haplo %in% c("M4″67"), "M4", haplo1),
         haplo1 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo1),
         haplo1 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo1),
         haplo1 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo1),
         haplo1 = ifelse(haplo %in% c("U4'9"), "U4", haplo1),
         haplo1 = ifelse(haplo %in% c("W1+119"), "W1", haplo1),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo3),
         haplo3 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo3),
         haplo3 = ifelse(haplo %in% c("D*"), "D", haplo3),
         haplo3 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo3),
         haplo3 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo3),
         haplo3 = ifelse(haplo %in% c("M4″67"), "M4", haplo3),
         haplo3 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo3),
         haplo3 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo3),
         haplo3 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo3),
         haplo3 = ifelse(haplo %in% c("U4'9"), "U4", haplo3),
         haplo3 = ifelse(haplo %in% c("W1+119"), "W1", haplo3),
         haplo3 = ifelse((is.na(haplo3) | haplo3==".."), haplo, haplo3)
  ) %>% setDT() %>% filter(haplo1!="Unspecified")

# library(writexl)
# write_xlsx(dat_ancient_SEA, "dat_ancient_SEA.xlsx")
dat_ancient_SEA <- read_excel("dat_ancient_SEA_extra.xlsx") %>%
  mutate(Time=paste(Time1, "(", Time1_LB, "-", Time1_UB, ")", "BP", sep = "")) %>% setDT() %>% filter(haplo1!="Unspecified")

hap_ancient_SEA <- dat_ancient_SEA[, .N, by = .(haplo1, country)] %>% arrange(desc(N))
hap_ancient_SEA1 <- dat_ancient_SEA[, .N, by = .(haplo1)] %>% arrange(desc(N))

### Haplo

country_ancient_SEA <- dat_ancient_SEA[, .N, by = .(country, haplo1)]
country_ancient_SEA <- country_ancient_SEA %>%
  group_by(country) %>% arrange(haplo1, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(viridis)
g5 <- ggplot(country_ancient_SEA) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplo1)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g5
ggsave(filename = file.path("figures", "country_ancient_haplo_plus.png"), width = 15, height = 10)

an_SEA <- dat_ancient_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1),
         count=1,
         haploY=ifelse((is.na(haploY) | haploY==".."), "Unspecified", haploY),
         countY=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>% ungroup() %>%
  group_by(country, haploY) %>% mutate(sumY=sum(countY), maxY=max(sumY)) %>%
  group_by(country) %>% arrange(desc(maxY)) %>% 
  mutate(orderY=order(maxY, decreasing = T), haploY_max=haploY[orderY==1]) %>% ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo1, haploY, haplo_max, haplo1_max, haploY_max, Age, Location, Label, Date, country)) %>% 
  filter(haplo1!="Unspecified")

an_SEA_sf <- merge(an_SEA, SEA0p_sf, by=c("country"))
an_SEA_plot <- an_SEA_sf %>% st_as_sf(crs = 4326)

countries <- SEA0p_sf
countries_coords <- st_coordinates(st_centroid(SEA0p_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_ancient_SEA %>%
  dplyr::rename(ID=country) %>%
  group_by(haplo1) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  dplyr::rename(key=haplo1, value=N) %>%
  select(-percent) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

an_SEA_plot_max <- an_SEA_plot %>% group_by(country) %>% slice(1)
an_SEA_max <- an_SEA_plot_max %>% st_drop_geometry() %>% select(country, haplo1_max) %>% rename(Country=country, `Max Haplogroup`=haplo1_max) %>% setDT()

library("ggpmisc")

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) +
  geom_sf(data=an_SEA_plot_max, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo1), data = an_SEA_plot, position=position_jitter(width=1,height=1), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo1, label = Date), data = an_SEA_plot, size = 2.5, hjust=-0.1, vjust=0.5, position=position_jitter(width=1,height=1), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:217], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(an_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=28), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15), 
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia (+)")
ggsave(filename = file.path("figures", "Ancient_SEA_plus_date.png"), width = 49, height = 33)

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) + 
  geom_sf(data=an_SEA_plot_max, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo1), data = an_SEA_plot, position=position_jitter(width=1,height=1), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo1, label = Time), data = an_SEA_plot, size = 2.5, hjust=-0.1, vjust=0.5, position=position_jitter(width=1,height=1), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:217], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(an_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=28), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15), 
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia (+)")
ggsave(filename = file.path("figures", "Ancient_SEA_plus_time.png"), width = 49, height = 33)

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) + 
  geom_sf(data=an_SEA_plot_max, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo1), data = an_SEA_plot, position=position_jitter(width=-0.5,height=0.5), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo1, label = haplo1), data = an_SEA_plot, size = 8, hjust=-0.1, vjust=0.5, position=position_jitter(width=-1.5,height=1.5), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:217], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(an_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=35), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia (+)")
ggsave(filename = file.path("figures", "Ancient_SEA_plus_label.png"), width = 49, height = 33)

## Ancient DNA

load("data/SEA0_sf.RData")
SEA0_sf <- SEA0_sf %>% rename(country=NAME_0)

library(readxl)
ancient <- read_excel("all-ancient-dna-2-07-73-full.xlsx")
ancient_SEA <- ancient %>% filter(Country %in% c("Brunei", "Cambodia", "Indonesia", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Timor-Leste", "Vietnam"))
dat_ancient_SEA <- ancient_SEA %>%
  rename(haplo="mtDNA-haplogroup",
         country="Country") %>%
  mutate(haplo1=str_extract(haplo, "^([A-Z])\\d\\w"),
         haplo1=ifelse(is.na(haplo1), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(is.na(haplo3), haplo, haplo3)
  ) %>% setDT()

hap_ancient_SEA <- dat_ancient_SEA[, .N, by = .(haplo1, country)] %>% arrange(desc(N))
hap_ancient_SEA1 <- dat_ancient_SEA[, .N, by = .(haplo1)] %>% arrange(desc(N))

### Haplo

country_ancient_SEA <- dat_ancient_SEA[, .N, by = .(country, haplo)]
country_ancient_SEA <- country_ancient_SEA %>%
  group_by(country) %>% arrange(haplo, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(viridis)
g5 <- ggplot(country_ancient_SEA) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplo)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=3, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g5
ggsave(filename = file.path("figures", "country_ancient_haplo.png"), width = 15, height = 10)

an_SEA <- dat_ancient_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1),
         count=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>%
  ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo1, haplo_max, haplo1_max, Age, Location, Label, Date, country))

an_SEA_sf <- merge(an_SEA, SEA0_sf, by=c("country"))
an_SEA_plot <- an_SEA_sf %>% st_as_sf(crs = 4326)

countries <- SEA0_sf
countries_coords <- st_coordinates(st_centroid(SEA0_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_ancient_SEA %>%
  mutate(haplo1=str_extract(haplo, "^([A-Z])\\d\\w"),
         haplo1=ifelse(is.na(haplo1), haplo, haplo1),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1)) %>%
  rename(ID=country) %>%
  group_by(haplo1) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haplo1, value=N) %>%
  select(-percent, -haplo) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

ggplot() + geom_sf() + geom_sf(data=an_SEA_plot, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo1), data = an_SEA_plot, size = 8) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo1, label = Location), data = an_SEA_plot, size = 7.5, hjust=-0.1, vjust=1, label.size = 1) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:33], color=NA, alpha=0.8) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Ancient_SEA.png"), width = 49, height = 33)

### Haplo 1

country_ancient_SEA1 <- dat_ancient_SEA[, .N, by = .(country, haplo1)]
country_ancient_SEA1 <- country_ancient_SEA1 %>%
  group_by(country) %>% arrange(haplo1, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup() %>%
  mutate(haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1))

g6 <- ggplot(country_ancient_SEA1) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplo1)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g6
ggsave(filename = file.path("figures", "country_ancient_haplo1.png"), width = 15, height = 10)

dat_f <- dat %>% mutate(count=1) %>% setDF()

an_SEA1 <- dat_ancient_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1),
         count=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>%
  ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo1, haplo_max, haplo1_max, Age, Location, Label, Date, country))

an_SEA1_sf <- merge(an_SEA1, SEA0_sf, by=c("country"))
an_SEA1_plot <- an_SEA1_sf %>% st_as_sf(crs = 4326)

countries <- SEA0_sf
countries_coords <- st_coordinates(st_centroid(SEA0_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_ancient_SEA1 %>%
  rename(ID=country) %>%
  group_by(haplo1) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haplo1, value=N) %>%
  select(-percent) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

ggplot() + geom_sf() + geom_sf(data=an_SEA1_plot, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo1), data = an_SEA1_plot, size = 8) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo1, label = Location), data = an_SEA1_plot, size = 7.5, hjust=-0.1, vjust=1, label.size = 1) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:33], color=NA, alpha=0.8) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Ancient_SEA1.png"), width = 49, height = 33)

### Haplo 3

country_ancient_SEA3 <- dat_ancient_SEA[, .N, by = .(country, haplo3)]
country_ancient_SEA3 <- country_ancient_SEA3 %>%
  group_by(country) %>% arrange(haplo3, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup() %>%
  mutate(haplo3=ifelse((is.na(haplo3) | haplo3==".."), "Unspecified", haplo3))

g7 <- ggplot(country_ancient_SEA3) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haplo3)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(1, "cm")) +
  coord_flip()
g7
ggsave(filename = file.path("figures", "country_ancient_haplo3.png"), width = 15, height = 10)

dat_f <- dat %>% mutate(count=1) %>% setDF()

an_SEA3 <- dat_ancient_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo3=ifelse((is.na(haplo3) | haplo3==".."), "Unspecified", haplo3),
         count=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>%
  group_by(country, haplo3) %>% mutate(sum3=sum(count), max3=max(sum3)) %>%
  group_by(country) %>% arrange(desc(max3)) %>%
  mutate(order3=order(max3, decreasing = T), haplo3_max=haplo3[order3==1]) %>%
  ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo3, haplo_max, haplo3_max, Age, Location, Label, Date, country))

an_SEA3_sf <- merge(an_SEA3, SEA0_sf, by=c("country"))
an_SEA3_plot <- an_SEA3_sf %>% st_as_sf(crs = 4326)

countries <- SEA0_sf
countries_coords <- st_coordinates(st_centroid(SEA0_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_ancient_SEA3 %>%
  rename(ID=country) %>%
  group_by(haplo3) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haplo3, value=N) %>%
  select(-percent) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

ggplot() + geom_sf() + geom_sf(data=an_SEA3_plot, aes(fill=haplo3_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haplo3), data = an_SEA3_plot, size = 8) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haplo3, label = Location), data = an_SEA3_plot, size = 7.5, hjust=-0.1, vjust=1, label.size = 1) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:28], color=NA, alpha=0.8) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Ancient_SEA3.png"), width = 49, height = 33)

## Ancient Y-DNA plus

load("data/SEA0p_sf.RData")
SEA0p_sf <- SEA0p_sf %>% rename(country=NAME_0)

library(readxl)
ancient <- read_excel("all-ancient-dna-2-07-73-full.xlsx")
ancient_SEA <- ancient %>% filter(Country %in% c("Bangladesh", "Brunei", "Cambodia", "China", "Indonesia", "India", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Taiwan", "Timor-Leste", "Vietnam"))
dat_ancientY_SEA <- ancient_SEA %>%
  rename(haplo="mtDNA-haplogroup",
         haploY="Y-Haplotree-Public",
         country="Country") %>%
  mutate(haploY1=ifelse(!is.na(`Y-New`), `Y-New`, `Y-Simple`),
         haploY1=ifelse(`Y-New` %in% c("CT low coverage", "D1(xD1a1a1a1b,xD1a1b1a,xD1a2) low coverage", "F low coverage", "F or J2a1a low coverage", "FTDNA: Forms a new branch D-Y65054 under D-PH344 with a Big Y tester from Kazakhstan", "ISOGG2015?", "ISOGG2015??", "ISOGG2016??", "ISOGG2018?", "J2b:FGC3945.1/FGC3945.2/FGC3945/Z526.1/Z526.2/Z526, J2b2a:AM01367/Z605", "K low coverage", "K2 low coverage", "M9, M526, M1221, L405, P295, F115", "NO1 low coverage", "O low coverage", "O1 low coverage", "O1a low coverage", "O1b1a1a1a(xO1b1a1a1a1a1)", "O2a1 low coverage", "O2a2b2a low coverage"), `Y-Simple`, haploY1),
         haploY2=ifelse(haploY1 %in% c("C", "C2", "CT", "D", "F", "K", "N", "N-L735", "NO", "O", "O2", "P*", "Q", "Q1", "R2"), haploY1, `Y-Simple`),
         haploY2=ifelse(!haploY1 %in% c("C", "C2", "CT", "D", "F", "K", "N", "N-L735", "NO", "O", "O2", "P*", "Q", "Q1", "R2"), str_extract(haploY1, "^([A-Z]+)\\d\\w"), haploY2),
         haplo1=ifelse(!(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189", "C4+152", "C4+152+16093", "D*", "E (95.07%)", "G1 (94.06%)", "M4″67", "n/a", "n/a (<2x)", "n/a (exome capture)", "R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)", "R2+13500", "U4'9", "W1+119")), str_extract(haplo, "^([A-Z])\\d\\w"), haplo),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haploY=ifelse((is.na(haploY) | haplo==".."), "Unspecified", haploY),
         haplo1 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo1),
         haplo1 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo1),
         haplo1 = ifelse(haplo %in% c("D*"), "D", haplo1),
         haplo1 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo1),
         haplo1 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo1),
         haplo1 = ifelse(haplo %in% c("M4″67"), "M4", haplo1),
         haplo1 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo1),
         haplo1 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo1),
         haplo1 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo1),
         haplo1 = ifelse(haplo %in% c("U4'9"), "U4", haplo1),
         haplo1 = ifelse(haplo %in% c("W1+119"), "W1", haplo1),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo3),
         haplo3 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo3),
         haplo3 = ifelse(haplo %in% c("D*"), "D", haplo3),
         haplo3 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo3),
         haplo3 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo3),
         haplo3 = ifelse(haplo %in% c("M4″67"), "M4", haplo3),
         haplo3 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo3),
         haplo3 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo3),
         haplo3 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo3),
         haplo3 = ifelse(haplo %in% c("U4'9"), "U4", haplo3),
         haplo3 = ifelse(haplo %in% c("W1+119"), "W1", haplo3),
         haplo3 = ifelse((is.na(haplo3) | haplo3==".."), haplo, haplo3)
  ) %>% setDT() %>% filter(!is.na(haploY) & haploY!="Unspecified")

# library(writexl)
# write_xlsx(dat_ancient_SEA, "dat_ancient_SEA.xlsx")
# dat_ancient_SEA <- read_excel("dat_ancient_SEA_extra.xlsx") %>%
#   mutate(Time=paste(Time1, "(", Time1_LB, "-", Time1_UB, ")", "BP", sep = "")) %>% setDT() %>% filter(haplo1!="Unspecified")

hapY_ancient_SEA <- dat_ancientY_SEA[, .N, by = .(haploY, country)] %>% arrange(desc(N))
hapY_ancient_SEA1 <- dat_ancientY_SEA[, .N, by = .(haploY)] %>% arrange(desc(N))

### Haplo

country_ancientY_SEA <- dat_ancientY_SEA[, .N, by = .(country, haploY)]
country_ancientY_SEA <- country_ancientY_SEA %>%
  group_by(country) %>% arrange(haploY, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

library(viridis)
g5Y <- ggplot(country_ancientY_SEA) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(country), y=percent, fill=factor(haploY)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=5, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Country") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "cm")) +
  coord_flip()
g5Y
ggsave(filename = file.path("figures", "country_ancient_haploY_plus.png"), width = 15, height = 10)

anY_SEA <- dat_ancientY_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1),
         count=1,
         haploY=ifelse((is.na(haploY) | haploY==".."), "Unspecified", haploY),
         countY=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>% ungroup() %>%
  group_by(country, haploY) %>% mutate(sumY=sum(countY), maxY=max(sumY)) %>%
  group_by(country) %>% arrange(desc(maxY)) %>% 
  mutate(orderY=order(maxY, decreasing = T), haploY_max=haploY[orderY==1]) %>% ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo1, haploY, haplo_max, haplo1_max, haploY_max, Age, Location, Label, Date, country)) %>% 
  filter(!is.na(haploY) & haploY!="Unspecified")

anY_SEA_sf <- merge(anY_SEA, SEA0p_sf, by=c("country"))
anY_SEA_plot <- anY_SEA_sf %>% st_as_sf(crs = 4326)

countries <- SEA0p_sf
countries_coords <- st_coordinates(st_centroid(SEA0p_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_ancientY_SEA %>%
  rename(ID=country) %>%
  group_by(haploY) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haploY, value=N) %>%
  select(-percent) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

anY_SEA_plot_max <- anY_SEA_plot %>% group_by(country) %>% slice(1)
anY_SEA_max <- anY_SEA_plot_max %>% st_drop_geometry() %>% select(country, haploY_max) %>% rename(Country=country, `Max Y-Haplogroup`=haploY_max) %>% setDT()

library("ggpmisc")

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) +
  geom_sf(data=anY_SEA_plot_max, aes(fill=haploY_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haploY), data = anY_SEA_plot, position=position_jitter(width=1,height=1), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haploY, label = Date), data = anY_SEA_plot, size = 2.5, hjust=-0.1, vjust=0.5, position=position_jitter(width=1,height=1), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:90], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(anY_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=28), 
        axis.text.x = element_text(size=15), 
        axis.text.y = element_text(size=15), 
        legend.text=element_text(size=20), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human Y-DNA (Y-DNA) Haplogroups in Southeast Asia (+)")
ggsave(filename = file.path("figures", "AncientY_SEA_plus_date.png"), width = 49, height = 33)

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) + 
  geom_sf(data=anY_SEA_plot_max, aes(fill=haploY_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haploY), data = anY_SEA_plot, position=position_jitter(width=-0.5,height=0.5), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haploY, label = haploY), data = anY_SEA_plot, size = 8, hjust=-0.1, vjust=0.5, position=position_jitter(width=-1.5,height=1.5), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:90], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(anY_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=40), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human Y-DNA (Y-DNA) Haplogroups in Southeast Asia (+)")
ggsave(filename = file.path("figures", "AncientY_SEA_plus_label.png"), width = 49, height = 33)

## O1b - Ancient Y-DNA plus

load("data/SEA0p_sf.RData")
SEA0p_sf <- SEA0p_sf %>% rename(country=NAME_0)

library(readxl)
ancient <- read_excel("all-ancient-dna-2-07-73-full.xlsx")
ancient_SEA <- ancient %>% filter(Country %in% c("Bangladesh", "Brunei", "Cambodia", "China", "Indonesia", "India", "Laos", "Malaysia", "Myanmar", "Philippines", "Singapore", "Thailand", "Taiwan", "Timor-Leste", "Vietnam"))
O_ancientY_SEA <- ancient_SEA %>%
  rename(haplo="mtDNA-haplogroup",
         haploY="Y-Haplotree-Public",
         country="Country") %>%
  mutate(haploY1=ifelse(!is.na(`Y-New`), `Y-New`, `Y-Simple`),
         haploY1=ifelse(`Y-New` %in% c("CT low coverage", "D1(xD1a1a1a1b,xD1a1b1a,xD1a2) low coverage", "F low coverage", "F or J2a1a low coverage", "FTDNA: Forms a new branch D-Y65054 under D-PH344 with a Big Y tester from Kazakhstan", "ISOGG2015?", "ISOGG2015??", "ISOGG2016??", "ISOGG2018?", "J2b:FGC3945.1/FGC3945.2/FGC3945/Z526.1/Z526.2/Z526, J2b2a:AM01367/Z605", "K low coverage", "K2 low coverage", "M9, M526, M1221, L405, P295, F115", "NO1 low coverage", "O low coverage", "O1 low coverage", "O1a low coverage", "O1b1a1a1a(xO1b1a1a1a1a1)", "O2a1 low coverage", "O2a2b2a low coverage"), `Y-Simple`, haploY1),
         haploY2=ifelse(haploY1 %in% c("C", "C2", "CT", "D", "F", "K", "N", "N-L735", "NO", "O", "O2", "P*", "Q", "Q1", "R2"), haploY1, `Y-Simple`),
         haploY2=ifelse(!haploY1 %in% c("C", "C2", "CT", "D", "F", "K", "N", "N-L735", "NO", "O", "O2", "P*", "Q", "Q1", "R2"), str_extract(haploY1, "^([A-Z]+)\\d\\w"), haploY2),
         haplo1=ifelse(!(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189", "C4+152", "C4+152+16093", "D*", "E (95.07%)", "G1 (94.06%)", "M4″67", "n/a", "n/a (<2x)", "n/a (exome capture)", "R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)", "R2+13500", "U4'9", "W1+119")), str_extract(haplo, "^([A-Z])\\d\\w"), haplo),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haploY=ifelse((is.na(haploY) | haplo==".."), "Unspecified", haploY),
         haplo1 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo1),
         haplo1 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo1),
         haplo1 = ifelse(haplo %in% c("D*"), "D", haplo1),
         haplo1 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo1),
         haplo1 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo1),
         haplo1 = ifelse(haplo %in% c("M4″67"), "M4", haplo1),
         haplo1 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo1),
         haplo1 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo1),
         haplo1 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo1),
         haplo1 = ifelse(haplo %in% c("U4'9"), "U4", haplo1),
         haplo1 = ifelse(haplo %in% c("W1+119"), "W1", haplo1),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(haplo %in% c("A+152", "A+152+16362", "A+152+16362+200", "A+152+16362+16189"), "A+", haplo3),
         haplo3 = ifelse(haplo %in% c("C4+152", "C4+152+16093"), "C4+", haplo3),
         haplo3 = ifelse(haplo %in% c("D*"), "D", haplo3),
         haplo3 = ifelse(haplo %in% c("E (95.07%)"), "E", haplo3),
         haplo3 = ifelse(haplo %in% c("G1 (94.06%)"), "G1", haplo3),
         haplo3 = ifelse(haplo %in% c("M4″67"), "M4", haplo3),
         haplo3 = ifelse(haplo %in% c("n/a", "n/a (<2x)", "n/a (exome capture)"), "Unspecified", haplo3),
         haplo3 = ifelse(haplo %in% c("R+16189", "R+16189C (76.45%)", "R+16189C (80.01%)", "R+16189C (81.67%)"), "R+", haplo3),
         haplo3 = ifelse(haplo %in% c("R2+13500"), "R2+", haplo3),
         haplo3 = ifelse(haplo %in% c("U4'9"), "U4", haplo3),
         haplo3 = ifelse(haplo %in% c("W1+119"), "W1", haplo3),
         haplo3 = ifelse((is.na(haplo3) | haplo3==".."), haplo, haplo3)
  ) %>% setDT() %>% filter(haploY2 %in% c("NO", "O", "O1a", "O1b", "O2", "O2a", "O3a"))

country_O_SEA <- O_ancientY_SEA[, .N, by = .(country, haploY1)]
country_O_SEA <- country_O_SEA %>%
  group_by(country) %>% arrange(haploY1, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

O_SEA <- O_ancientY_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1),
         count=1,
         haploY1=ifelse((is.na(haploY1) | haploY1==".."), "Unspecified", haploY1),
         countY1=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>% ungroup() %>%
  group_by(country, haploY1) %>% mutate(sumY1=sum(countY1), maxY1=max(sumY1)) %>%
  group_by(country) %>% arrange(desc(maxY1)) %>% 
  mutate(orderY1=order(maxY1, decreasing = T), haploY1_max=haploY1[orderY1==1]) %>% ungroup() %>%
  select(c(`Object-ID`, Latitude, Longitude, Sex, haplo, haplo1, haploY1, haplo_max, haplo1_max, haploY1_max, Age, Location, Label, Date, country)) %>% 
  filter(!is.na(haploY1) & haploY1!="Unspecified")

O_SEA_sf <- merge(O_SEA, SEA0p_sf, by=c("country"))
O_SEA_plot <- O_SEA_sf %>% st_as_sf(crs = 4326)

countries <- SEA0p_sf
countries_coords <- st_coordinates(st_centroid(SEA0p_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = countries$country)

res <- country_O_SEA %>%
  rename(ID=country) %>%
  group_by(haploY1) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haploY1, value=N) %>%
  select(-percent) %>%
  arrange(key)

res <- res %>% left_join(countries_coords)

dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

O_SEA_plot_max <- O_SEA_plot %>% group_by(country) %>% slice(1)
O_SEA_max <- O_SEA_plot_max %>% st_drop_geometry() %>% select(country, haploY1_max) %>% rename(Country=country, `Dominant Y-Haplogroup`=haploY1_max) %>% setDT()

library("ggpmisc")

ggplot() + 
  # geom_sf(data=SEA0p_sf, aes(fill="white"), alpha=0.1) + 
  geom_sf(data=O_SEA_plot_max, aes(fill=haploY1_max), lwd=0, alpha=0.6) +
  geom_point(aes(x = Longitude, y = Latitude,  colour = haploY1), data = O_SEA_plot, position=position_jitter(width=-0.5,height=0.5), size = 6) +
  geom_label(aes(x = Longitude, y = Latitude,  colour = haploY1, label = haploY1), data = O_SEA_plot, size = 8, hjust=-0.1, vjust=0.5, position=position_jitter(width=-1.5,height=1.5), label.size = 0.5) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:28], color=NA, alpha=0.8) +
  annotate(geom = "table", x = 80, y = -10, label = list(O_SEA_max), size = 6.5) +
  scale_fill_discrete(name="") +
  scale_color_discrete(name="") +
  guides(fill=guide_legend(nrow=4, byrow=TRUE)) +
  theme_bw() +
  theme(text = element_text(size=40), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(2, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Ancient Human Y-DNA Haplogroup O in Southeast Asia (+)")
ggsave(filename = file.path("figures", "AncientY_SEA_plus_haplo_O.png"), width = 49, height = 33)

### Make pie

# make_pie <- function(dt, title = NA, legend.position = 0){
#   if(is.na(title)){
#     title <- unique(dt$ID)
#   }
#   ggplot() +
#     geom_bar(data = dt,
#              aes(x = "", y = value, fill = key),
#              stat = "identity", width = 1) +
#     coord_polar("y") +
#     theme_void() +
#     theme(legend.position = legend.position) +
#     ggtitle(title)
# }
# 
# con1 <- make_pie(dplyr::filter(res, Country == 1))
# con2 <- make_pie(dplyr::filter(res, Country == 2))
# con3 <- make_pie(dplyr::filter(res, Country == 3))
# con4 <- make_pie(dplyr::filter(res, Country == 4))
# con5 <- make_pie(dplyr::filter(res, Country == 5))
# con6 <- make_pie(dplyr::filter(res, Country == 6))
# con7 <- make_pie(dplyr::filter(res, Country == 7))
# con8 <- make_pie(dplyr::filter(res, Country == 8))
# 
# (gg_countries <- ggplot(data = countries) +
#     geom_sf() +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0 )) +
#     theme(legend.position = 0,
#           plot.margin = unit(c(0,0,0,0), "cm"))
# )
# 
# leg <- get_legend(make_pie(res, "", legend.position = "left"))
# 
# draw_plot_loc <- function(plot, dt){
#   draw_plot(plot, x = dt$X[1], y = dt$Y[1],
#             height = 0.2)
# }
# 
# (all <-
# ggdraw(gg_countries) +
#     draw_plot_loc(con1, dplyr::filter(res, Country == 1)) +
#     draw_plot_loc(con2, dplyr::filter(res, Country == 2)) +
#     draw_plot_loc(con3, dplyr::filter(res, Country == 3)) +
#     draw_plot_loc(con4, dplyr::filter(res, Country == 4)) +
#     draw_plot_loc(con5, dplyr::filter(res, Country == 5)) +
#     draw_plot_loc(con6, dplyr::filter(res, Country == 6)) +
#     draw_plot_loc(con7, dplyr::filter(res, Country == 7)) +
#     draw_plot_loc(con8, dplyr::filter(res, Country == 8)) +
#     draw_plot_loc(con9, dplyr::filter(res, Country == 9))
#   )
# 
# cowplot::plot_grid(all, leg, rel_widths = c(1, 0.1))

## Present DNA

dat_f <- dat %>% mutate(count=1) %>% setDF()

pre_SEA <- dat %>% 
  mutate(haplo1=ifelse(!(haplogroup3 %in% c("A+152", "A+152+16362", " A+152+16362+200", "R+16189")), str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(is.na(haplo3), haplo, haplo3),
         count=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo1) %>% mutate(sum1=sum(count), max1=max(sum1)) %>%
  group_by(country) %>% arrange(desc(max1)) %>% 
  mutate(order1=order(max1, decreasing = T), haplo1_max=haplo1[order1==1]) %>%
  ungroup() %>% setDT()

pre_SEA_sf <- merge(pre_SEA, SEA0_sf, by=c("country"))
pre_SEA_plot <- pre_SEA_sf %>% st_as_sf(crs = 4326)

country_present_SEA <- pre_SEA[, .N, by = .(country, haplo1)]
country_present_SEA <- country_present_SEA %>%
  group_by(country) %>% arrange(haplo1, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

res_pre <- country_present_SEA %>%
  rename(ID=country) %>%
  group_by(haplo1) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haplo1, value=N) %>%
  select(-percent) %>%
  arrange(key)

res_pre <- res_pre %>% left_join(countries_coords)

dt_res <- spread(res_pre, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

pre_SEA_plot_max <- pre_SEA_plot %>% group_by(country) %>% slice(1)

ggplot() + geom_sf() + geom_sf(data=pre_SEA_plot_max, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_scatterpie(aes(x=x, y=y, r=0.6), data=dt_x, cols = colnames(dt_x)[1:180], color=NA, alpha=0.8) +
  guides(fill=guide_legend(nrow=18, byrow=TRUE)) +
  scale_fill_discrete(name="") +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Present Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Present_SEA.png"), width = 49, height = 33)

### Haplo 3

pre_SEA3 <- dat %>% 
  mutate(haplo1=ifelse(!(haplogroup3 %in% c("A+152", "A+152+16362", "A+152+16362+200", "R+16189")), str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), haplo, haplo1),
         haplo2 = substr(haplo, 1, 1),
         haplo3 = str_extract(haplo, "^([A-Z])\\d+"),
         haplo3 = ifelse(haplogroup3 %in% c("A+152", "A+152+16362", "A+152+16362+200"), "A+", haplo3),
         haplo3 = ifelse(haplogroup3 %in% c("R+16189"), "R+", haplo3),
         haplo3 = ifelse(is.na(haplo3), haplogroup3, haplo3),
         count=1) %>%
  group_by(country, haplo) %>%  mutate(sum=sum(count), max=max(sum)) %>%
  group_by(country) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo_max=haplo[order==1]) %>% ungroup %>%
  group_by(country, haplo3) %>% mutate(sum3=sum(count), max3=max(sum3)) %>%
  group_by(country) %>% arrange(desc(max3)) %>% 
  mutate(order3=order(max3, decreasing = T), haplo3_max=haplo3[order3==1]) %>%
  ungroup() %>% setDT()

pre_SEA3_sf <- merge(pre_SEA3, SEA0_sf, by=c("country"))
pre_SEA3_plot <- pre_SEA3_sf %>% st_as_sf(crs = 4326)

country_present_SEA3 <- pre_SEA3[, .N, by = .(country, haplo3)]
country_present_SEA3 <- country_present_SEA3 %>%
  group_by(country) %>% arrange(haplo3, .by_group = TRUE) %>% 
  mutate(percent=(N*100)/sum(N)) %>% ungroup()

res_pre <- country_present_SEA3 %>%
  rename(ID=country) %>%
  group_by(haplo3) %>%
  mutate(Country=order(ID)) %>%
  ungroup() %>%
  rename(key=haplo3, value=N) %>%
  select(-percent) %>%
  arrange(key)

res_pre <- res_pre %>% left_join(countries_coords)

dt_res <- spread(res_pre, key = key, value = value) %>% replace(is.na(.), 0)
DT <- dt_res %>% select(-c(ID, X, Y, Country))
m<-as.matrix(DT)
ID <- dt_res$ID
Country <- dt_res$Country
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(Country=order(ID)) %>% left_join(countries_coords) %>% rename(x=X, y=Y)
dt_x <- dt %>% select(-c(Country, ID))

pre_SEA3_plot_max <- pre_SEA3_plot %>% group_by(country) %>% slice(1)

ggplot() + geom_sf() + geom_sf(data=pre_SEA3_plot_max, aes(fill=haplo3_max), lwd=0, alpha=0.6) +
  geom_scatterpie(aes(x=x, y=y, r=1), data=dt_x, cols = colnames(dt_x)[1:139], color=NA, alpha=0.8) +
  guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
  scale_fill_discrete(name="") +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=30), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Present Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Present_SEA3.png"), width = 49, height = 33)

## Ethnicity

# SEA1_sf <- rbind(brunei_BRN1_sf, indonesia_IDN1_sf, cambodia_KHM1_sf, laos_LAO1_sf, myanmar_MMR1_sf, malaysia_MYS1_sf, philippines_PHL1_sf, singapore_SGP1_sf,
#                  thailand_THA1_sf, easttimor_TLS1_sf, vietnam_VNM1_sf)
# 
# save(SEA1_sf, file = "data/SEA1_sf.RData")

load("data/SEA1_sf.RData")
SEA1_sf <- SEA1_sf %>% dplyr::rename(country=NAME_0, location=NAME_1, type=ENGTYPE_1)
SEA1_sf$location <- stri_trans_general(SEA1_sf$location, "Latin-ASCII")
SEA1_sf$location <- trimws(gsub("\\s+", " ", SEA1_sf$location))
SEA1_sf <- SEA1_sf %>% dplyr::select(country, location, type, geometry)

library(readxl)
SEA <- read_excel("Changed_SEA_haplogroups.xlsx")
ethnic_SEA <- SEA[,-c(1,3,4,6,7,8)] %>% dplyr::rename(country=Country, ethnicity=Ethnicity, location=Location, sum=`Sample size`)
ethnic_SEA <- ethnic_SEA %>% 
  mutate(across(where(is.character), ~na_if(., "-"))) %>%
  mutate_at(c(5:612), as.numeric) %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% dplyr::select(-sum) %>% setDT()

ethnic_SEA <- gather(ethnic_SEA, haplo, N, -country, -ethnicity, -location, factor_key=TRUE) %>% setDT()

dat_ethnic_SEA <- ethnic_SEA %>%
  mutate(haplo1=ifelse(!(haplo %in% c("A+152", "A+152+16362", " A+152+16362+200", "R+16189")), str_extract(haplo, "^([A-Z])\\d\\w"), haplo),
         haplo1=ifelse(is.na(haplo1), haplo, haplo1),
         haplo1=case_when(haplo1 %in% c("134", "157", "161", "168", "171", "174", "241", "257", "259", "260", "262", "268", "279", "281", "295", "30", "304", 
                                        "313", "315", "32", "351", "375", "381", "425", "433", "439", "483", "489", "496", "50", "501", "51", "53", "530", 
                                        "563", "564", "567", "589", "607", "62", "73", "78", "90") ~ haplo, 
                          TRUE ~ haplo1),
         location=case_when(location=="Brunei (Borneo)" ~ "Brunei and Muara",
                            location=="Seim Riep (or Siem Riep)" ~ "Siemreab",
                            location=="Prey Veng" ~ "Prey Veng",
                            location=="Banteay Meanchey" ~ "Banteay Meanchey",
                            location=="Kampong Thom" ~ "Kampong Thum",
                            location=="Kampong Cham" ~ "Kampong Cham",
                            location=="Kratie" ~ "Kracheh",
                            location=="Takeo" ~ "Takev",
                            location=="Battanbang" ~ "Batdambang",
                            location=="Kampong Chahnang" ~ "Kampong Chhnang",
                            location=="Pursat" ~ "Pouthisat",
                            location=="Phnom Penh" ~ "Phnom Penh",
                            location=="Kampot" ~ "Kampot",
                            location=="Kandal" ~ "Kandal",
                            location=="Oddar Meanchey" ~ "Otdar Mean Chey",
                            location=="Koh Kong" ~ "Kaoh Kong",
                            location=="Svay Rieng" ~ "Svay Rieng",
                            location=="Alor" ~ "Nusa Tenggara Timur",
                            location=="Ambon" ~ "Maluku",
                            location=="Bali" ~ "Bali",
                            location=="South Kalimantan" ~ "Kalimantan Selatan",
                            location=="Padang" ~ "Sumatera Barat",
                            location=="Sumatra" ~ "Sumatera Barat",
                            location=="Java, Demak" ~ "Jawa Tengah",
                            location=="Sumatra, Kutaradja" ~ "Sumatera Selatan",
                            location=="Java" ~ "Jawa Tengah",
                            location=="West New Guinea" ~ "Papua Barat",
                            location=="Manado" ~ "Sulawesi Utara",
                            location=="Sulawesi" ~ "Sulawesi Tengah",
                            location=="Sulawesi, Manado" ~ "Sulawesi Utara",
                            location=="Sumatra, Padang" ~ "Sumatera Barat",
                            location=="Sumatra, Palembang" ~ "Sumatera Selatan",
                            location=="Palangkaraya" ~ "Kalimantan Tengah",
                            location=="South Borneo" ~ "Kalimantan Selatan",
                            location=="Pagaralam" ~ "Sumatera Selatan",
                            location=="Java" ~ "Jawa Tengah",
                            location=="Toraja" ~ "Sulawesi Selatan",
                            location=="Ujung Pandang" ~ "Sulawesi Selatan",
                            location=="Waigapu (Sumba)" ~ "Nusa Tenggara Timur",
                            location=="Sumba" ~ "Nusa Tenggara Timur",
                            location=="Banjarmasin (Borneo)" ~ "Kalimantan Selatan",
                            location=="Palangkaraya (Borneo)" ~ "Kalimantan Tengah",
                            location=="North Laos" ~ "Louangphrabang",
                            location=="Central Laos" ~ "Vientiane",
                            location=="Kedah" ~ "Kedah",
                            location=="Acheh-Kedah Yan" ~ "Kedah",
                            location=="Kelantan" ~ "Kelantan",
                            location=="Johor" ~ "Johor",
                            location=="Johor Pontian" ~ "Johor",
                            location=="Banjar Perak Kuala Kurau" ~ "Kedah",
                            location=="Perak, Banjar Malay" ~ "Kedah",
                            location=="Kelantan" ~ "Kelantan",
                            location=="Johor Semerah Jawa" ~ "Johor",
                            location=="Johor Muar Jawa" ~ "Johor",
                            location=="Sabah" ~ "Johor",
                            location=="Kelantan" ~ "Kelantan",
                            location=="Kelantan Kota Bahru" ~ "Kelantan",
                            location=="Negeri Sembilan" ~ "Negeri Sembilan",
                            location=="Minangkabau-Negeri Sembilan Lenggeng" ~ "Negeri Sembilan",
                            location=="Kelantan RantauPanjang" ~ "Kelantan",
                            location=="Perak, Rawa Malay" ~ "Perak",
                            location=="Kota Kinabalu (Borneo)" ~ "Sabah",
                            location=="Yangon" ~ "Yangon",
                            location=="Pakokku" ~ "Mon",
                            location=="Bataan" ~ "Bataan",
                            location=="Zambales" ~ "Zambales",
                            location=="Iriga" ~ "Camarines Sur",
                            location=="Batan Archipelago" ~ "Bataan",
                            location=="Cebu" ~ "Cebu",
                            location=="Philippines Bohol" ~ "Bohol",
                            location=="Luzon" ~ "Bulacan",
                            location=="North Thailand" ~ "Chiang Mai",
                            location=="Northeast Thailand" ~ "Samut Prakan",
                            location=="Northern Thailand" ~ "Chiang Mai",
                            location=="Central Thailand" ~ "Bangkok Metropolis",
                            location=="Mergui Archipelago" ~ "Tanintharyi",
                            location=="West Thailand" ~ "Surat Thani",
                            location=="Baucau" ~ "Baucau",
                            location=="Liquica" ~ "Liquica",
                            location=="Cova Lima" ~ "Covalima",
                            location=="Viqueque" ~ "Viqueque",
                            location=="Aileu" ~ "Aileu",
                            location=="Ermera" ~ "Ermera",
                            location=="Dili" ~ "Dili",
                            location=="Manufahi" ~ "Manufahi",
                            ethnicity=="Cham" ~ "Ninh Thuan",
                            ethnicity=="CoLao" ~ "Ha Giang",
                            ethnicity=="Dao" ~ "Ha Giang",
                            ethnicity=="Kinh" ~ "Ha Noi",
                            ethnicity=="Stieng" ~ "Binh Phuoc",
                            ethnicity=="Ede" ~ "Dak Lak",
                            ethnicity=="Giarai" ~ "Gia Lai",
                            ethnicity=="HaNhi" ~ "Lai Chau",
                            ethnicity=="Hmong" ~ "Dien Bien",
                            ethnicity=="Hui" ~ "Ha Giang",
                            ethnicity=="LaChi" ~ "Ha Giang",
                            ethnicity=="LaHu" ~ "Lai Chau",
                            ethnicity=="LoLo" ~ "Ha Giang",
                            ethnicity=="Mang" ~ "Lai Chau",
                            ethnicity=="Nung" ~ "Lang Son",
                            ethnicity=="PaThen" ~ "Ha Giang",
                            ethnicity=="PhuLa" ~ "Lao Cai",
                            ethnicity=="SiLa" ~ "Lai Chau",
                            ethnicity=="Tay" ~ "Lang Son",
                            ethnicity=="Thai" ~ "Son La",
                            ethnicity=="Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung" ~ "Ha Noi",
                            ethnicity=="Yao" ~ "Ha Giang",
                            ethnicity=="Kinh, Tay, Dao, Hmong, Muong, Hoa, Khmer, Nung" ~ "Ha Noi",
                            ethnicity=="Tay Nung" ~ "Lao Cai",
                            TRUE ~ location),) %>% setDT()

hap_ethnic_SEA <- dat_ethnic_SEA[, .N, by = .(haplo1, ethnicity, location)] %>% arrange(desc(N))
hap_ethnic_SEA1 <- dat_ethnic_SEA[, .N, by = .(haplo1)] %>% arrange(desc(N))

pre_ethnic_SEA <- dat_ethnic_SEA %>% 
  group_by(ethnicity, haplo1) %>% 
  mutate(N1=sum(N)) %>% ungroup() %>%
  arrange(desc(N1), .by_group = TRUE) %>%
  dplyr::select(ethnicity, haplo1, N1) %>%
  group_by(ethnicity, haplo1) %>% dplyr::slice(1) %>% ungroup() %>%
  group_by(ethnicity) %>% arrange(haplo1, .by_group = TRUE) %>% 
  mutate(percent=(N1*100)/sum(N1)) %>% ungroup()

g8 <- ggplot(pre_ethnic_SEA) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(ethnicity), y=percent, fill=factor(haplo1)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Ethnicity") +
  scale_y_continuous(name = "Percent") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  coord_flip()
g8
ggsave(filename = file.path("figures", "pre_ethnicity_haplo.png"), width = 20, height = 15)

g9 <- ggplot(pre_ethnic_SEA) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(ethnicity), y=percent, fill=factor(haplo1)), position = "stack", stat="identity", alpha=0.5) +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  scale_fill_viridis(discrete=TRUE) +
  scale_x_discrete(name = "Ethinicity") +
  scale_y_continuous(name = "Percent") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.3, "cm")) +
  coord_polar()
g9
ggsave(filename = file.path("figures", "pre_ethnicity_haplo_polar.png"), width = 20, height = 15)

# library(writexl)
# write_xlsx(dat_ethnic_SEA, "dat_ethnic_SEA.xlsx")
# write_xlsx(locations, "locations.xlsx")

ethnicity_SEA <- dat_ethnic_SEA %>% 
  mutate(haplo=ifelse((is.na(haplo) | haplo==".."), "Unspecified", haplo),
         haplo1=ifelse((is.na(haplo1) | haplo1==".."), "Unspecified", haplo1)) %>%
  group_by(country, ethnicity, location, haplo1) %>%  mutate(N1=sum(N)) %>% ungroup() %>%
  arrange(desc(N1), .by_group = TRUE) %>%
  dplyr::select(country, ethnicity, location, haplo1, N1) %>%
  group_by(country, ethnicity, location, haplo1) %>% dplyr::slice(1) %>% ungroup() %>%
  group_by(country, ethnicity, location) %>% arrange(haplo1, .by_group = TRUE) %>% 
  mutate(percent=(N1*100)/sum(N1)) %>% ungroup() %>%
  group_by(country, ethnicity, location, haplo1) %>%  mutate(sum=sum(N1), max=max(sum)) %>%
  group_by(country, ethnicity, location) %>% arrange(desc(max)) %>% mutate(order=order(max, decreasing = T), haplo1_max=haplo1[order==1]) %>% ungroup %>%
  dplyr::select(-country)

ethnicity_SEA_sf <- merge(ethnicity_SEA, SEA1_sf, by=c("location"))
ethnicity_SEA_plot <- ethnicity_SEA_sf %>% st_as_sf(crs = 4326)

locations <- SEA1_sf
locations_coords <- st_coordinates(st_centroid(SEA1_sf)) %>%
  data.frame(stringsAsFactors = FALSE) %>%
  mutate(ID = locations$location)

table(dat_ethnic_SEA$location %in% locations$location)

res <- ethnicity_SEA %>%
  dplyr::rename(ID=location) %>%
  dplyr::rename(key=haplo1, value=N1) %>%
  arrange(key)

res <- res %>% left_join(locations_coords) %>% filter(!is.na(ID))

# dt_res <- spread(res, key = key, value = value) %>% replace(is.na(.), 0)

dt_res <- res %>% 
  group_by(ethnicity, ID) %>% 
  mutate(Var = paste0("Val", row_number())) %>% 
  spread(key, value) %>% replace(is.na(.), 0) %>%
  ungroup()

DT <- dt_res %>% dplyr::select(-c(1:10))
m<-as.matrix(DT)
ID <- dt_res$ID
location <- dt_res$ID
dt <- aggregate(m, data.frame(ID),sum) %>% setDT()
# cbind(id = x[, 1], x[, -1]/rowSums(x[, -1]))
library(janitor)
dt <- dt %>% 
  adorn_percentages() %>% 
  dplyr::mutate_if(is.numeric, funs(. * 100)) %>%
  mutate(location=order(ID)) %>% left_join(locations_coords) %>% dplyr::rename(x=X, y=Y)
dt_x <- dt %>% dplyr::select(-c(location, ID))

ggplot() + geom_sf(data=SEA1_sf, aes(fill="white"), alpha=0.1) + 
  geom_sf(data=ethnicity_SEA_plot, aes(fill=haplo1_max), lwd=0, alpha=0.6) +
  geom_sf_text(data=ethnicity_SEA_plot, mapping=aes(label = ethnicity), stat = "sf_coordinates", position = "identity", check_overlap = T) +
  geom_scatterpie(aes(x=x, y=y, r=0.6), data=dt_x, cols = colnames(dt_x)[1:219], color=NA, alpha=0.8) +
  guides(fill=guide_legend(nrow=10, byrow=TRUE)) +
  scale_fill_discrete(name="") +
  theme_bw() +
  theme(text = element_text(size=36), 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        legend.text=element_text(size=25), 
        legend.key.size = unit(1, "cm"),
        legend.position = "bottom") +
  ggtitle("Geographic distribution of Present Human mitochondrial DNA (mtDNA) Haplogroups in Southeast Asia")
ggsave(filename = file.path("figures", "Ethnicity_SEA_edit2.png"), width = 49, height = 33)


################ SUBCLADES #######################

dat1 <- dat %>% mutate(haplogroup4 = str_extract(haplo, "^([A-Z])\\d\\w"),
                       haplogroup4 = ifelse(is.na(haplogroup4), haplogroup3, haplogroup4))
hap_present_SEA <- dat1[, .N, by = .(haplogroup4)] %>% arrange(desc(N))
hap_precountry_SEA <- dat1[, .N, by = .(country, haplogroup4)] %>% arrange(desc(N))

# Haplogroup F1f

hap_F1f <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "F1f")
nbin_F1f <- nbin[labels(nbin) %in% hap_F1f$name]
class(nbin_F1f)
dnbin_F1f<-dist.dna(nbin_F1f, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_F1f<-nj(dnbin_F1f)
library(ggtree)
ggt_F1f<-ggtree::ggtree(tree_F1f, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_F1f

library(treeio)
png("figures/Viet_Indo_F1f.png", width = 1200, height = 1600)
zoom(tree_F1f, grep("Vietnam|Indonesia", tree_F1f$tip.label, value = TRUE))
dev.off()
png("figures/Viet_Thai_F1f.png", width = 1200, height = 1600)
zoom(tree_F1f, grep("Vietnam|Thailand", tree_F1f$tip.label, value = TRUE))
dev.off()

# Haplogroup F1a
hap_F1a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "F1a")
nbin_F1a <- nbin[labels(nbin) %in% hap_F1a$name]
class(nbin_F1a)
dnbin_F1a<-dist.dna(nbin_F1a, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_F1a<-nj(dnbin_F1a)
library(ggtree)
ggt_F1a<-ggtree::ggtree(tree_F1a, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_F1a

library(treeio)
png("figures/Viet_Indo_F1a.png", width = 1200, height = 1600)
zoom(tree_F1a, grep("Vietnam|Indonesia", tree_F1a$tip.label, value = TRUE))
dev.off()
png("figures/Viet_Thai_F1a.png", width = 1200, height = 1600)
zoom(tree_F1a, grep("Vietnam|Thailand", tree_F1a$tip.label, value = TRUE))
dev.off()

F1a <- file[hap_F1a$name]
writeXStringSet(F1a, "data/F1a.fasta")

# cat(file="data/F1a.fasta", paste(paste0(">",names(nbin_F1a)),
#                                  sapply(nbin_F1a, paste, collapse=""), sep="\n"), sep="\n")

# Haplogroup B5a

hap_B5a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B5", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B5a")

B5a <- file[hap_B5a$name]
writeXStringSet(B5a, "data/B5a.fasta")

# Haplogroup M7b
hap_M7b <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "M7b")

M7b <- file[hap_M7b$name]
writeXStringSet(M7b, "data/M7b.fasta")

hap_M7b1a1 <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w\\d\\w\\d"), haplogroup3)) %>%
  filter(haplogroup4 == "M7b1a1")

M7b1a1 <- file[hap_M7b1a1$name]
writeXStringSet(M7b1a1, "data/M7b1a1.fasta")

# Haplogroup B4a

hap_B4a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B4a")

B4a <- file[hap_B4a$name]
writeXStringSet(B4a, "data/B4a.fasta")

# Haplogroup M7c
hap_M7c <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "M7c")

M7c <- file[hap_M7c$name]
writeXStringSet(M7c, "data/M7c.fasta")

# Haplogroup B4c

hap_B4c <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B4c")

B4c <- file[hap_B4c$name]
writeXStringSet(B4c, "data/B4c.fasta")

# Haplogroup F1f
hap_F1f <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "F1f")

F1f <- file[hap_F1f$name]
writeXStringSet(F1f, "data/F1f.fasta")

# Haplogroup N9a
hap_N9a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="N9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "N9a")

N9a <- file[hap_N9a$name]
writeXStringSet(N9a, "data/N9a.fasta")

# Haplogroup R9b
hap_R9b <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="R9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "R9b")

R9b <- file[hap_R9b$name]
writeXStringSet(R9b, "data/R9b.fasta")

# Haplogroup M74
hap_M74 <- dat %>%filter(haplogroup3 == "M74")

M74 <- file[hap_M74$name]
writeXStringSet(M74, "data/M74.fasta")

# Haplogroup C7a
hap_C7a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="C7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "C7a")

C7a <- file[hap_C7a$name]
writeXStringSet(C7a, "data/C7a.fasta")

################ SUBCLADE DISTRIBUTION ##########################

# Haplogroup F1a

F1a <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.F1a=(sum(subset(., haplogroup4=="F1a")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="F1a")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "F1a") %>% group_by(country) %>% slice(1)

F1a_sf <- merge(F1a, SEA0_sf, by=c("country"))
F1a_plot <- F1a_sf %>% select(prop.pop, prop.F1a, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=F1a_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(F1a_plot$prop.F1a[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup F1a")
ggsave(filename = file.path("figures", "F1a-PerPop.png"), width = 49, height = 33)

geom_point(aes(x = Lon, y = Lat,  colour = Facility), data = df, size = 0.5) + 
  theme(legend.position="bottom")

# Haplogroup B5a

B5a <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B5", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.B5a=(sum(subset(., haplogroup4=="B5a")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="B5a")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "B5a") %>% group_by(country) %>% slice(1)

B5a_sf <- merge(B5a, SEA0_sf, by=c("country"))
B5a_plot <- B5a_sf %>% select(prop.pop, prop.B5a, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=B5a_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(B5a_plot$prop.B5a[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup B5a")
ggsave(filename = file.path("figures", "B5a-PerPop.png"), width = 49, height = 33)

# Haplogroup F1f

F1f <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.F1f=(sum(subset(., haplogroup4=="F1f")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="F1f")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "F1f") %>% group_by(country) %>% slice(1)

F1f_sf <- merge(F1f, SEA0_sf, by=c("country"))
F1f_plot <- F1f_sf %>% select(prop.pop, prop.F1f, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=F1f_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(F1f_plot$prop.F1f[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup F1f")
ggsave(filename = file.path("figures", "F1f-PerPop.png"), width = 49, height = 33)

# Haplogroup F1c

F1c <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.F1c=(sum(subset(., haplogroup4=="F1c")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="F1c")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "F1c") %>% group_by(country) %>% slice(1)

F1c_sf <- merge(F1c, SEA0_sf, by=c("country"))
F1c_plot <- F1c_sf %>% select(prop.pop, prop.F1c, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=F1c_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=115, y=23, color="red", size=18, label= paste(round(F1c_plot$prop.F1c[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup F1c")
ggsave(filename = file.path("figures", "F1c-PerPop.png"), width = 49, height = 33)

# Haplogroup F1e

F1e <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="F1", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.F1e=(sum(subset(., haplogroup4=="F1e")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="F1e")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "F1e") %>% group_by(country) %>% slice(1)

F1e_sf <- merge(F1e, SEA0_sf, by=c("country"))
F1e_plot <- F1e_sf %>% select(prop.pop, prop.F1e, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=F1e_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=115, y=23, color="red", size=18, label= paste(round(F1e_plot$prop.F1e[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup F1e")
ggsave(filename = file.path("figures", "F1e-PerPop.png"), width = 49, height = 33)

# Haplogroup M7b

M7b <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.M7b=(sum(subset(., haplogroup4=="M7b")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="M7b")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "M7b") %>% group_by(country) %>% slice(1)

M7b_sf <- merge(M7b, SEA0_sf, by=c("country"))
M7b_plot <- M7b_sf %>% select(prop.pop, prop.M7b, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=M7b_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(M7b_plot$prop.M7b[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup M7b")
ggsave(filename = file.path("figures", "M7b-PerPop.png"), width = 49, height = 33)

# Haplogroup M7c

M7c <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.M7c=(sum(subset(., haplogroup4=="M7c")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="M7c")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "M7c") %>% group_by(country) %>% slice(1)

M7c_sf <- merge(M7c, SEA0_sf, by=c("country"))
M7c_plot <- M7c_sf %>% select(prop.pop, prop.M7c, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=M7c_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(M7c_plot$prop.M7c[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup M7c")
ggsave(filename = file.path("figures", "M7c-PerPop.png"), width = 49, height = 33)

# Haplogroup B4a

B4a <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.B4a=(sum(subset(., haplogroup4=="B4a")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="B4a")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "B4a") %>% group_by(country) %>% slice(1)

B4a_sf <- merge(B4a, SEA0_sf, by=c("country"))
B4a_plot <- B4a_sf %>% select(prop.pop, prop.B4a, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=B4a_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(B4a_plot$prop.B4a[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup B4a")
ggsave(filename = file.path("figures", "B4a-PerPop.png"), width = 49, height = 33)

# Haplogroup N9a

N9a <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="N9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.N9a=(sum(subset(., haplogroup4=="N9a")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="N9a")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "N9a") %>% group_by(country) %>% slice(1)

N9a_sf <- merge(N9a, SEA0_sf, by=c("country"))
N9a_plot <- N9a_sf %>% select(prop.pop, prop.N9a, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=N9a_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(N9a_plot$prop.N9a[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup N9a")
ggsave(filename = file.path("figures", "N9a-PerPop.png"), width = 49, height = 33)

# Haplogroup R9b

R9b <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="R9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.R9b=(sum(subset(., haplogroup4=="R9b")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="R9b")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "R9b") %>% group_by(country) %>% slice(1)

R9b_sf <- merge(R9b, SEA0_sf, by=c("country"))
R9b_plot <- R9b_sf %>% select(prop.pop, prop.R9b, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=R9b_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(R9b_plot$prop.R9b[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup R9b")
ggsave(filename = file.path("figures", "R9b-PerPop.png"), width = 49, height = 33)

# Haplogroup R9c

R9c <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="R9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.R9c=(sum(subset(., haplogroup4=="R9c")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="R9c")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "R9c") %>% group_by(country) %>% slice(1)

R9c_sf <- merge(R9c, SEA0_sf, by=c("country"))
R9c_plot <- R9c_sf %>% select(prop.pop, prop.R9c, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=R9c_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(R9c_plot$prop.R9c[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup R9c")
ggsave(filename = file.path("figures", "R9c-PerPop.png"), width = 49, height = 33)

# Haplogroup G2b

G2b <- dat_f %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="G2", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3),
         prop.G2b=(sum(subset(., haplogroup4=="G2b")[, "count"])/sum(.$count))*100) %>%
  group_by(country) %>% mutate(prop.pop=(sum(subset(dat_f, haplogroup4=="G2b")[, "count"])/sum(dat_f$count))*100) %>% ungroup() %>%
  filter(haplogroup4 == "G2b") %>% group_by(country) %>% slice(1)

G2b_sf <- merge(G2b, SEA0_sf, by=c("country"))
G2b_plot <- G2b_sf %>% select(prop.pop, prop.G2b, geometry) %>% st_as_sf(crs = 4326)
ggplot() + geom_sf() + geom_sf(data=G2b_plot, aes(fill=prop.pop), lwd=0) +
  annotate(geom="text", x=130, y=23, color="red", size=18, label= paste(round(G2b_plot$prop.G2b[1],1), "% dân số", sep = "")) +
  scale_fill_viridis("Tỉ lệ (%) trên dân số", direction = -1, option = "viridis") +
  theme(text = element_text(size=45), axis.text.x = element_text(size=30), axis.text.y = element_text(size=30), legend.text=element_text(size=30), legend.key.size = unit(2, "cm")) +
  ggtitle("Haplogroup G2b")
ggsave(filename = file.path("figures", "G2b-PerPop.png"), width = 49, height = 33)

# Haplogroup M

hap_M <- dat %>% filter(haplogroup1=="M")
nbin_M <- nbin[labels(nbin) %in% hap_M$name]
class(nbin_M)
dnbin_M<-dist.dna(nbin_M, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M<-nj(dnbin_M)
library(ggtree)
ggt_M<-ggtree::ggtree(tree_M, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_M

# Haplogroup M1
hap_M1 <- dat %>% filter(haplogroup2=="M1")
nbin_M1 <- nbin[labels(nbin) %in% hap_M1$name]
class(nbin_M1)
dnbin_M1<-dist.dna(nbin_M1, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M1<-nj(dnbin_M1)
library(ggtree)
ggt_M1<-ggtree::ggtree(tree_M1, layout = "circular", cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_M1

library(treeio)

png("figures/Viet_Indo_M1.png", width = 1200, height = 800)
zoom(tree_M1, grep("Vietnam|Indonesia", tree_M1$tip.label, value = TRUE))
dev.off()

zoom(tree_M1, grep("Vietnam|Thailand", tree_M1$tip.label, value = TRUE))

# Haplogroup M2
hap_M2 <- dat %>% filter(haplogroup2=="M2")
nbin_M2 <- nbin[labels(nbin) %in% hap_M2$name]
class(nbin_M2)
dnbin_M2<-dist.dna(nbin_M2, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M2<-nj(dnbin_M2)
library(ggtree)
ggt_M2<-ggtree::ggtree(tree_M2, layout = "circular", cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_M2

library(treeio)
png("figures/Viet_Indo_M2.png", width = 1200, height = 800)
zoom(tree_M2, grep("Vietnam|Indonesia", tree_M2$tip.label, value = TRUE))
dev.off()
zoom(tree_M2, grep("Vietnam|Thailand", tree_M2$tip.label, value = TRUE))

# Haplogroup M5
hap_M5 <- dat %>% filter(haplogroup2=="M5")
nbin_M5 <- nbin[labels(nbin) %in% hap_M5$name]
class(nbin_M5)
dnbin_M5<-dist.dna(nbin_M5, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M5<-nj(dnbin_M5)
library(ggtree)
ggt_M5<-ggtree::ggtree(tree_M5, layout = "circular", cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_M5

library(treeio)
png("figures/Viet_Indo_M5.png", width = 1200, height = 800)
zoom(tree_M5, grep("Vietnam|Indonesia", tree_M5$tip.label, value = TRUE))
dev.off()
zoom(tree_M5, grep("Vietnam|Thailand", tree_M5$tip.label, value = TRUE))

# Haplogroup M7
hap_M7 <- dat %>% filter(haplogroup2=="M7")
nbin_M7 <- nbin[labels(nbin) %in% hap_M7$name]
class(nbin_M7)
dnbin_M7<-dist.dna(nbin_M7, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M7<-nj(dnbin_M7)
library(ggtree)
ggt_M7<-ggtree::ggtree(tree_M7, layout = "circular", cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_M7

library(treeio)
png("figures/Viet_Indo_M7.png", width = 2400, height = 4800, res = 130)
zoom(tree_M7, grep("Vietnam|Indonesia", tree_M7$tip.label, value = TRUE))
dev.off()
zoom(tree_M7, grep("Vietnam|Thailand", tree_M7$tip.label, value = TRUE))

groupInfo <- split(tree_M7$tip.label, gsub("_\\w+", "", tree_M7$tip.label))
tree_M7 <- groupOTU(tree_M7, groupInfo)
p_M7 <- ggtree(tree_M7, aes(color=group)) + geom_tiplab() + xlim(NA, 23)
zoom(tree_M7, grep("Vietnam", tree_M7$tip.label), xmax_adjust=2)

tre_M7<-ladderize(tree_M7)
ggtree(tre_M7, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)+
  geom_highlight(node = 20, fill="purple", alpha=0.2)

tre<-ladderize(tree)
ggtree(tre, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)+
  geom_highlight(node = 20, fill="purple", alpha=0.2)

# Haplogroup M7b
hap_M7b <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "M7b")
nbin_M7b <- nbin[labels(nbin) %in% hap_M7b$name]
class(nbin_M7b)
dnbin_M7b<-dist.dna(nbin_M7b, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_M7b<-nj(dnbin_M7b)
library(ggtree)
ggt_M7b<-ggtree::ggtree(tree_M7b, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_M7b

M7b <- file[hap_M7b$name]
writeXStringSet(M7b, "data/M7b.fasta")

hap_M7b1a1 <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="M7", str_extract(haplo, "^([A-Z])\\d\\w\\d\\w\\d"), haplogroup3)) %>%
  filter(haplogroup4 == "M7b1a1")

M7b1a1 <- file[hap_M7b1a1$name]
writeXStringSet(M7b1a1, "data/M7b1a1.fasta")

library(treeio)
png("figures/Viet_Indo_M7b.png", width = 1200, height = 1600)
zoom(tree_M7b, grep("Vietnam|Indonesia", tree_M7b$tip.label, value = TRUE))
dev.off()
zoom(tree_M7b, grep("Vietnam|Thailand", tree_M7b$tip.label, value = TRUE))

# Haplogroup N9a
hap_N9a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="N9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "N9a")
nbin_N9a <- nbin[labels(nbin) %in% hap_N9a$name]
class(nbin_N9a)
dnbin_N9a<-dist.dna(nbin_N9a, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_N9a<-nj(dnbin_N9a)
library(ggtree)
ggt_N9a<-ggtree::ggtree(tree_N9a, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_N9a

library(treeio)
png("figures/Viet_Indo_N9a.png", width = 1200, height = 1600)
zoom(tree_N9a, grep("Vietnam|Indonesia", tree_N9a$tip.label, value = TRUE))
dev.off()
zoom(tree_N9a, grep("Vietnam|Thailand", tree_N9a$tip.label, value = TRUE))

# Haplogroup R9b
hap_R9b <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="R9", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "R9b")
nbin_R9b <- nbin[labels(nbin) %in% hap_R9b$name]
class(nbin_R9b)
dnbin_R9b<-dist.dna(nbin_R9b, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_R9b<-nj(dnbin_R9b)
library(ggtree)
ggt_R9b<-ggtree::ggtree(tree_R9b, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_R9b

library(treeio)
png("figures/Viet_Indo_R9b.png", width = 1200, height = 1600)
zoom(tree_R9b, grep("Vietnam|Indonesia", tree_R9b$tip.label, value = TRUE))
dev.off()
zoom(tree_R9b, grep("Vietnam|Thailand", tree_R9b$tip.label, value = TRUE))

# Haplogroup B4
hap_B4 <- dat %>% filter(haplogroup2=="B4")
nbin_B4 <- nbin[labels(nbin) %in% hap_B4$name]
class(nbin_B4)
dnbin_B4<-dist.dna(nbin_B4, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_B4<-nj(dnbin_B4)
library(ggtree)
ggt_B4<-ggtree::ggtree(tree_B4, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_B4

library(treeio)
zoom(tree_B4, grep("Vietnam|Indonesia", tree_B4$tip.label, value = TRUE))
zoom(tree_B4, grep("Vietnam|Thailand", tree_B4$tip.label, value = TRUE))

# Haplogroup B4a
hap_B4a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B4a")
nbin_B4a <- nbin[labels(nbin) %in% hap_B4a$name]
class(nbin_B4a)
dnbin_B4a<-dist.dna(nbin_B4a, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_B4a<-nj(dnbin_B4a)
library(ggtree)
ggt_B4a<-ggtree::ggtree(tree_B4a, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_B4a

library(treeio)
png("figures/Viet_Indo_B4a.png", width = 1200, height = 1600)
zoom(tree_B4a, grep("Vietnam|Indonesia", tree_B4a$tip.label, value = TRUE))
dev.off()
zoom(tree_B4a, grep("Vietnam|Thailand", tree_B4a$tip.label, value = TRUE))

# Haplogroup B4c
hap_B4c <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B4", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B4c")
nbin_B4c <- nbin[labels(nbin) %in% hap_B4c$name]
class(nbin_B4c)
dnbin_B4c<-dist.dna(nbin_B4c, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_B4c<-nj(dnbin_B4c)
library(ggtree)
ggt_B4c<-ggtree::ggtree(tree_B4c, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_B4c

library(treeio)
png("figures/Viet_Indo_B4c.png", width = 1200, height = 1600)
zoom(tree_B4c, grep("Vietnam|Indonesia", tree_B4a$tip.label, value = TRUE))
dev.off()
zoom(tree_B4c, grep("Vietnam|Thailand", tree_B4c$tip.label, value = TRUE))

# Haplogroup B5
hap_B5 <- dat %>% filter(haplogroup2=="B5")
nbin_B5 <- nbin[labels(nbin) %in% hap_B5$name]
class(nbin_B5)
dnbin_B5<-dist.dna(nbin_B5, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_B5<-nj(dnbin_B5)
library(ggtree)
ggt_B5<-ggtree::ggtree(tree_B5, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_B5

library(treeio)
zoom(tree_B5, grep("Vietnam|Indonesia", tree_B5$tip.label, value = TRUE))
zoom(tree_B5, grep("Vietnam|Thailand", tree_B5$tip.label, value = TRUE))

# Haplogroup B5a
hap_B5a <- dat %>% 
  mutate(haplogroup4 = ifelse(haplogroup3=="B5", str_extract(haplo, "^([A-Z])\\d\\w"), haplogroup3)) %>%
  filter(haplogroup4 == "B5a")
nbin_B5a <- nbin[labels(nbin) %in% hap_B5a$name]
class(nbin_B5a)
dnbin_B5a<-dist.dna(nbin_B5a, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_B5a<-nj(dnbin_B5a)
library(ggtree)
ggt_B5a<-ggtree::ggtree(tree_B5a, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_B5a

library(treeio)
png("figures/Viet_Indo_B5a.png", width = 1200, height = 1600)
zoom(tree_B5a, grep("Vietnam|Indonesia", tree_B5a$tip.label, value = TRUE))
dev.off()
zoom(tree_B5a, grep("Vietnam|Thailand", tree_B5a$tip.label, value = TRUE))

# Haplogroup D4
hap_D4 <- dat %>% filter(haplogroup2=="D4")
nbin_D4 <- nbin[labels(nbin) %in% hap_D4$name]
class(nbin_D4)
dnbin_D4<-dist.dna(nbin_D4, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_D4<-nj(dnbin_D4)
library(ggtree)
ggt_D4<-ggtree::ggtree(tree_D4, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_D4

library(treeio)
zoom(tree_D4, grep("Vietnam|Indonesia", tree_D4$tip.label, value = TRUE))
zoom(tree_D4, grep("Vietnam|Thailand", tree_D4$tip.label, value = TRUE))

# Haplogroup F1
hap_F1 <- dat %>% filter(haplogroup2=="F1")
nbin_F1 <- nbin[labels(nbin) %in% hap_F1$name]
class(nbin_F1)
dnbin_F1<-dist.dna(nbin_F1, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_F1<-nj(dnbin_F1)
library(ggtree)
ggt_F1<-ggtree::ggtree(tree_F1, cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 2)
ggt_F1

library(treeio)
zoom(tree_F1, grep("Vietnam|Indonesia", tree_F1$tip.label, value = TRUE))
zoom(tree_F1, grep("Vietnam|Thailand", tree_F1$tip.label, value = TRUE))

# Haplogroup R9
hap_R9 <- dat %>% filter(haplogroup2=="R9")
nbin_R9 <- nbin[labels(nbin) %in% hap_R9$name]
class(nbin_R9)
dnbin_R9<-dist.dna(nbin_R9, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree_R9<-nj(dnbin_R9)
library(ggtree)
ggt_R9<-ggtree::ggtree(tree_R9, layout = "circular", cex = 0.8, aes(color=branch.length))+
  scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_tiplab(align=TRUE, size=2)+
  geom_treescale(y = - 5, color = "coral4", fontsize = 4)
ggt_R9

library(treeio)
zoom(tree_R9, grep("Vietnam|Indonesia", tree_R9$tip.label, value = TRUE))
zoom(tree_R9, grep("Vietnam|Thailand", tree_R9$tip.label, value = TRUE))

#########   EXTRACTION SEQUENCE AND HAPLOTYPE INFORMATION    ###############

nrow(nm)#confirmation of number of samples

ncol(nm)#confirmation of sequences size

sat2 <- NULL
for (i in 1:nrow(nm)) {
  sat2[i] <- paste(nm[i, ], collapse="")
}

sat2 <- toupper(sat2) #converts all letters to uppercase
sat3 <- unique(sat2) #gives only unique sequences from all sequences
sat3#that is, it gives complete sequences of haplotypes (20x373).
hfreq <- NULL
for (i in 1:length(sat3)) {
  hcount = 0
  s3 <- sat3[i]
  for (j in 1:length(sat2)) {
    s2 <- sat2[j]
    if (s3 == s2) {
      hcount <- (hcount + 1) #counts the number of individuals with the same haplotype sequence. 
      #print(paste(i, "yes", hcount))
    }
    #print(s2)
  }
  hname<-(paste("H",i, sep =""))
  hfreq[i] <- hcount
  #print(paste(hname, hcount, collapse = ""))
}   #haplotype frequency in the all samples

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat4 = list()
for (j in 1:len) {
  same <- TRUE
  first <- substr(sat3[1], j, j)
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    if (first != ch1) {
      str <- paste(j, first, ch1)
      print(str)
      same <- FALSE
      break
    }
  }
  if (!same) {
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    }
    sat4[cnt] <- ss
    cnt <- cnt + 1
  }
}#it gives the mutation points and the nucleotide substitutions

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat5 = list() 
for (j in 1:len) { #scan all columnns and if all elements are the same do not copy
  same <- TRUE
  first <- substr(sat3[1], j, j)
  scol <- first
  for (i in 2:length(sat3)) {
    ch1 <- substr(sat3[i], j, j)
    scol <- paste(scol, ch1, sep="")
    if (first != ch1) {
      str <- paste(j, first, ch1)
      #print(str)
      same <- FALSE
      #break
    }
  }
  if (!same) {
    scol <- paste("V_", cnt, " ", scol, sep="")
    ss <- NULL
    for (i in 1:length(sat3)) {
      ss <- paste(ss, substr(sat3[i], j, j), sep="")
    } 
    sat5[cnt] <- ss
    cnt <- cnt + 1
  }
}

sat6 <- as.matrix(sat5)
mat6 = matrix(nrow=nrow(sat6), ncol=nchar(sat6[1]))
for (i in 1:nrow(mat6)) {
  s <- as.vector(strsplit(as.character(sat5[i]), ""))
  for (j in 1:ncol(mat6)) {
    mat6[i, j] <- as.character(s[[1]][j])
  }
}
mat7 <- t(mat6) #sequences of haplotypes and variable sites matrix (20x41)
write.table(mat7,file="mat7.txt", quote=FALSE, sep="\t")
hname<-paste("H", 1:nrow(mat7), sep = "")
rownames(mat7)=hname
write.table(mat7,file="mat7.txt", quote=FALSE, sep="\t") 

str4 <- NULL
str4[1] <- paste(mat7[1, ], collapse="")
for (i in 2:nrow(mat7)) {
  tmp <- NULL
  for (j in 1:ncol(mat7)) {
    chr = "."
    if(mat7[i, j] != mat7[1, j]) chr = mat7[i, j]
    tmp <- paste(tmp, chr, sep="")
  }
  str4[i] <- paste(tmp, collapse="")
}
nchar(str4[1]) #confirmation of number of variable sites
mstr4<-as.matrix(str4)
rownames(mstr4)<-hname
colnames(mstr4)<-paste("sequences length","(", ncol(mat7), "base pairs", ")")
pct<-round((as.matrix(hfreq)*100/colSums(as.matrix(hfreq))), 2)
colnames(pct)<-c("pct")
cmstr4<-as.data.frame(cbind(mstr4, hfreq, pct))
cmstr4
write.table(cmstr4,file="cmstr4.txt", quote=FALSE, sep="\t") 

#############  HAPLOTYPES FREQUENCY  ################
library(haplotypes)

kn<-as.dna(nbin)
kh<-haplotypes::haplotype(kn)

ncb <- as.matrix(labels(nbin))
n2 <- NULL
for (i in 1:nrow(ncb)) {
  n2[i] <- strsplit(ncb[i], '_')[[1]][1] #to get the names of the examples where the name and number are separated by an underscore
}
n2

for (i in 1:nrow(ncb)) {
  n2[i] <- strsplit(n2[i], ' ')[[1]][1] #to get the names of the examples where the name and number are separated by a space
}
n2

hf<-grouping(kh, factors=n2)
hf[["hapvec"]] <- NULL
dhf<-as.data.frame(hf$hapmat) #haplotype frequency matrix per populaion
rownames(dhf)<-paste("H", 1:nrow(mat7), sep = "")
dhf
write.table(dhf,file="dhf.txt", quote=FALSE, sep="\t")

############## DISTANCE MATRIX OF HAPLOTYPES  ############################
library(pegas)
dh<-dist.hamming(mat7) #from pegas package
dh
dhm<-as.matrix(dh)
write.table(dhm,file="dhm.txt", quote=FALSE, sep="\t")

############### HEAT MAP ###########################
sat2 <- NULL
for (i in 1:nrow(nm)) {
  sat2[i] <- paste(nm[i, ], collapse="")
}
sat2

sat2 <- toupper(sat2)
sat3 <- unique(sat2)
comat = matrix(nrow=length(sat3), ncol=length(sat3))
for (i in 1:length(sat3)) { 
  si <- sat3[i]
  for (j in 1:length(sat3)) { 
    sj <- sat3[j]
    difcnt = 0
    s1 = as.vector(strsplit(as.character(si), ""))
    s2 = as.vector(strsplit(as.character(sj), ""))
    for (k in 1:length(s1[[1]])) {
      if (s1[[1]][k] != s2[[1]][k]) {
        difcnt = difcnt + 1
      }
      comat[i, j] = difcnt
      #print(paste(i, " ", j, " ", difcnt))
    }
  }
}
comat	#is Hamming distance matrix
colnames(comat)<-paste("H", 1:nrow(comat), sep = "")
rownames(comat)<-paste("H", 1:nrow(comat), sep = "")
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE) #stats package

dev.new()
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)


pdf("heatmap.pdf", width = 7, height = 7)
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)
dev.off()

############## HAPLOTYPE NETWORK-INDIVIDUALS ##########################
library(pegas)
h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE)#constructs the haplotype network
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=rownames(nbin))
)

par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)

dev.new()
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=30, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)

pdf("hapind.pdf", width = 11, height = 5)
par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)
dev.off()

# M7b

library(pegas)
h_M7b<-pegas::haplotype(nbin_M7b, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname_M7b<-paste("H", 1:nrow(h_M7b), sep = "")
rownames(h_M7b)= paste(hname_M7b)
net_M7b<-haploNet(h_M7b, d = NULL, getProb = TRUE)#constructs the haplotype network
net_M7b
ind.hap<-with(
  utils::stack(setNames(attr(h_M7b, "index"), rownames(h_M7b))),
  table(hap=ind, individuals=rownames(nbin_M7b))
)

par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)

dev.new()
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=30, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)

pdf("hapind.pdf", width = 11, height = 5)
par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)
dev.off()


#For a better position of the name table, you can optimize x and y coordinates in the "legend()" argument.

############## HAPLOTYPE NETWORK-POPULATIONS ##########################

h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE) #constructs the haplotype network
net
ind.hap<-with(
  utils::stack(setNames(attr(h, "index"), rownames(h))),
  table(hap=ind, individuals=rownames(nbin)[values])
)

#We followed the name order in the ind.hap object for all coloring of samples in the haplotype network circles.
#You can see the order of sample names using the command "colnames(ind.hap)".
#The color of the samples has been assigned as the "bg" list for use in the plot () command, as below.
#The number of colors is determined by the number of samples in the same population (see "colnames (ind.hap)").
#Using the colors() command, you can choose different colors for your samples according to the populations.

bg<-c(rep("orange", 10), rep("red", 10))

#If you are not going to use a special order for the population names to be listed in the legend command (population names listed in the legend ("hapcol"), 
#you can automatically select and sort the population names using the following command ("sn2"). 
#This order is fully consistent with the order in the ind.hap object.

#un2<-unique(n2)
#sn2<-sort(un2)

#Then you must the use "sn2" list instead of "hapcol" list in the legend argument.

#But, we wanted to order names of the population names according to groups. 
#Therefore, we reordered the population names as "hapcol" list.
#So, you can see color transitions of the populations by the groups. 
#The population names were assigned using "hapcol" list in the legend () command.

hapcol<-c("Kinh", "Sumatra")

#At this point, the colors to be used in the list where the populations will be shown (legend) must match the colors you choose for the samples in the circles.
#Therefore, if you used the "sn2" list, you can directly use the "bg" list for the populations' color symbols. 
#Note, now the number of colors must be changed to "1" for each population.
#For example;
#bgp<-c(rep("dodgerblue4", 1), rep("olivedrab4",1),...)

#Then you must the use "bgp" list instead of "ubg" list in the legend argument (fill).

#But, we used different order for the population names (hapcol list).
#Thus, colors of the population were assigned as "ubg" list and used in the fill argument of the legend () command.
#The colors of the samples assigned in the bg list and the population colors (ubg) belonging samples must match each other.

ubg<-c(rep("orange", 1), rep("red", 1))

par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

########## USING "sn2" NAME LIST AND "bgp" COLOR LIST  ###############

#par(mar=c(0.001,0.001,0.001,0.001))
#plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
#legend(x=-36,y=53, sn2, fill=bgp, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

pdf("happop.pdf", width = 7, height = 4)
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)
dev.off()

############## HAPLOTYPE NETWORK-GROUPS ##########################

#To construct a haplotype network according to the groups, we assigned the group names of the samples according to the nbin object. 
#At this stage, we followed the name order in the nbin object using "labels(nbin)".

# nmname<-c(rep("Nature", 3), rep("Nature", 6), rep("Nature", 6), rep("Nature", 4), rep("Firm", 4), 
#           rep("Nature", 11), rep("Greenhouse", 3), rep("Nature", 3), rep("Greenhouse", 2), rep("Greenhouse", 12), 
#           rep("Nature", 6), rep("Firm", 4),  rep("Nature", 7), rep("Greenhouse", 13), rep("Nature", 6), rep("Greenhouse", 15), 
#           rep("Firm", 3), rep("Nature", 6), rep("Nature", 2), rep("Firm", 4))

nmname<-c(rep("Kinh", 1), rep("Sumatra", 2), rep("Kinh", 3), rep("Sumatra", 2), rep("Kinh", 5), rep("Sumatra", 6), rep("Kinh", 1))

ng<-nbin
rownames(ng)<-nmname
hg<-pegas::haplotype(ng, strict = FALSE, labels = hname, trailingGapsAsN = TRUE) #extracts haplotypes from DNAbin object
hg
netg<-haploNet(hg, d = NULL, getProb = TRUE) #constructs the haplotype network
netg
ind.hapg<-with(
  utils::stack(setNames(attr(hg, "index"), rownames(hg))),
  table(hap=ind, individuals=rownames(ng)[values])
)

#Group colors in the legend command which gbg argument below, were assigned according to the order of names in the ind.hapg object.
#To see name order, you can use the below command;
#as.data.frame(ind.hap)
#Color of samples was assigned using "gbg" argument in the plot () command.

gbg<-c(rep("red"), rep("blue"), rep("green"))

par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)


pdf("hapgrp.pdf", width = 8, height = 5) #save as pdf file
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)
dev.off()

############# CIRCULAR TREE-NJ-POPULATIONS #######################

library(pegas)
library(ggtree)
library(ggplot2)
library(phytools)
library(stats)

nname<-as.matrix(sort(labels(nbin)))#to construct population names 
krp<-list(Kinh = nname[1:10,],
          Sumatra = nname[11:20,])

nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin) #neighbor-joining tree estimation of Saitou and Nei (1987).

emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)

dev.new()
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)


pdf("njtreepop.pdf", width = 6.5, height = 5)
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)
dev.off()

##################### CIRCULAR TREE-NJ-DISTANCE ##############################

nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80")#computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin)#neighbor-joining tree estimation of Saitou and Nei (1987).

njdistree<-ggtree(tree,layout = 'circular', branch.length='branch.length', aes(color=branch.length), lwd = 0.5)+xlim(-0.1, NA)+
  geom_tiplab(names(nbin), size = 1.7, offset=0.002)+scale_color_continuous(high='lightskyblue1',low='coral4')+
  geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9) 

njdistree

dev.new()
njdistree


pdf("njdistree.pdf", width = 6.5, height = 5)
njdistree
dev.off()

############## HAPLOTYPE TREE ##########################

library(treeio)
D <- dist.hamming(mat7)#pegas package
class(D)
htre<-nj(D)#This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
bp <- boot.phylo(htre, mat7, B=100, function(x) nj(dist.hamming(x)))
bp2 <- data.frame(node=1:Nnode(htre) + Ntip(htre), bootstrap = bp)
htree <- full_join(htre, bp2, by="node")
boothap<-ggtree(htree, size=1, branch.length='branch.length')+geom_tiplab(size=4)+
  geom_nodepoint(aes(fill=cut(bootstrap, c(0,50,70,85,100))), shape=21, size=4)+
  theme_tree(legend.position=c(0.85, 0.2))+ 
  scale_fill_manual(values=c("black", "red", "pink1", "white"), guide='legend', name='Bootstrap Percentage (BP)',breaks=c('(85,100]', '(70,85]', '(50,70]', '(0,50]'), labels=expression(BP>=85, 70<=BP*"<85",50<=BP*"<70", BP<50))

boothap

dev.new()
boothap


pdf("boothap.pdf", width = 11, height = 5)
boothap
dev.off()

library("admixtools")

genotype_data = "/my/geno/prefix"
fit = qpgraph(genotype_data, example_graph)
fit$score

f2_blocks = f2_from_geno(genotype_data)
fit = qpgraph(f2_blocks, example_graph)

prefix = '/path/to/geno'
my_f2_dir = '/store/f2data/here/'

extract_f2(prefix, my_f2_dir)

admixtools::run_shiny_admixtools()

library(admixr)

snp_data <- eigenstrat(download_data())

result <- d(
  W = c("French", "Sardinian"), X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = snp_data
)

result

beast_file <- system.file("examples/MCC_FluA_H3.tree", 
                          package="ggtree")
beast_tree <- read.beast(beast_file)
ggtree(beast_tree, mrsd="2013-01-01") + theme_tree2()

samplelist <- read_tsv("analysis/samplelist.txt",
                       col_names = "sample")

read_delim("analysis/full_genome.filtered.numericChr.2.Q",
           col_names = paste0("Q",seq(1:2)),
           delim=" ")

### Coalescent MCMC

# Haplogroup F1
hap_F1 <- dat %>% filter(haplogroup2 == "F1")
nbin_F1 <- nbin[labels(nbin) %in% hap_F1$name]

data(woodmouse)
out <- coalescentMCMC(woodmouse)

#adjust plot margins
par(mar = c(1, 1, 1, 1))
plot(out)
getMCMCtrees()

res <- coalescentMCMC(woodmouse, 1e6, moves = c(1, 3)) # ~ 1 hr
plot(res) # surely hard to read
plot(subset(res, end = 1e3)) # plot only the first 1000 generations
acfplot(res)
acfplot(subset(res, 1e4, 100))

library(ape)
it = read.tree("example_newick.txt")  # autogenerated filename from the website.
getMRCA(it, c("Homo_sapiens", "Drosophila_melanogaster"))  # drosophila is one of only 2 arthropods in this tree, but arthropods all have an arthropod common ancestor
# prints "200", the node index of the MRCA.

library(phytools)
F1a_it = read.newick("F1a_newick.txt")
mrca(F1a_it, full = TRUE)
getMRCA(F1a_it, F1a_it$tip.label)

findMRCA(F1a_it, tips=NULL, type=c("node","height"))
anc<-findMRCA(F1a_it, F1a_it$tip.label)
plotTree(F1a_it,type="fan",fsize=0.7,lwd=1)
nodelabels(node=anc,frame="circle",pch=21,cex=1.5,
           bg="blue")
legend("topleft","most recent common ancestor\nof Puerto Rican TG anoles",
       pch=21,pt.cex=1.5,pt.bg="blue",cex=0.7,bty="n")
par(mar=c(5.1,4.1,4.1,2.1)) ## reset margin to default

library(devtools)
install_github("jhavsmith/startmrca")
library(startmrca)
help(run.startmrca)
run.startmrca(vcf.file = "examples/FIN_chr2_pos136608646.vcf.gz",
rec.file      = "examples/decode_recmap_sexaveraged.txt", 
sample.ids    = "examples/sample_ids.txt", 
refsample.ids = "examples/sample_ids.txt",
mut.rate      = 1.6e-8, 
nsel          = 50,
nanc          = 20,
chain.length  = 20,
nanc.post     = 10,
pos           = 136608646,
sel.allele    = 1)


