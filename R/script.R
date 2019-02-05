if (! requireNamespace("readr")) {
  install.packages("readr")
}

if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

if (! requireNamespace("igraph")) {
  install.packages("igraph")
}

myURL <- paste0("https://github.com/hyginn/",
                "BCB420-2019-resources/blob/master/HGNC.RData?raw=true")
load(url(myURL))  # loads HGNC data frame


tmp <- readr::read_tsv(file.path("../data", "homNOG.members.tsv"),
                         skip = 1, col_names = c("Species", "GeneID", "ID"))

head(tmp)

# A tibble: 6 x 3
#Species GeneID      ID
#<chr>   <chr>       <chr>
#1 homNOG  ENOG410RR6N 9606.ENSP00000473172
#2 homNOG  ENOG410S23Y 9606.ENSP00000473161
#3 homNOG  ENOG410RYSR 9606.ENSP00000473139
#4 homNOG  ENOG410S2EW 9606.ENSP00000473113
#5 homNOG  ENOG410S1H2 9606.ENSP00000473105
#6 homNOG  ENOG410RVFP 9606.ENSP00000473103


all(grepl("^9606\\.", tmp$ID)) #TRUE
tmp$ID <- gsub("^9606\\.", "", tmp$ID)
uENSP <- unique(tmp$ID) #Unique IDs to map

#==================== MAP ENSP TO HGNC =========================================

myMart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")

tmp <- biomaRt::getBM(filters = "ensembl_peptide_id",
                      attributes = c("ensembl_peptide_id",
                                     "hgnc_symbol"),
                      values = uENSP,
                      mart = myMart)
head(tmp)
#ensembl_peptide_id hgnc_symbol
#1    ENSP00000430236       ALG11
#2    ENSP00000430633      NPIPB5
#3    ENSP00000447378        MGAM
#4    ENSP00000446058         SHD
#5    ENSP00000468977      FLT3LG

nrow(tmp) #4872 symbols retrieved out of 5488

sum(duplicated(tmp$ensembl_peptide_id)) #1 duplicate
sum(! (uENSP) %in% tmp$ensembl_peptide_id) #616 - nothing returned

sum(is.na(tmp$ensembl_peptide_id)) #0
sum(tmp$ensembl_peptide_id == "") #0 None of the rows have ""

dupEnsp <- tmp$ensembl_peptide_id[duplicated(tmp$ensembl_peptide_id)]
tmp[tmp$ensembl_peptide_id %in% dupEnsp, ]
#ensembl_peptide_id hgnc_symbol
#2268    ENSP00000344961     PLEKHG7
#2269    ENSP00000344961    C12orf74

tmp[tmp$hgnc_symbol %in% c("C12orf74"), ]
tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("C12orf74")), ]
any(duplicated(tmp$ensembl_peptide_id)) #FALSE

ensp2sym <- tmp$hgnc_symbol
names(ensp2sym) <- tmp$ensembl_peptide_id
head(ensp2sym)
#ENSP00000430236 ENSP00000430633 ENSP00000447378 ENSP00000446058 ENSP00000468977
#    "ALG11"        "NPIPB5"          "MGAM"         "SHD"         "FLT3LG"

sel <- ! (uENSP %in% names(ensp2sym))
x <- rep(NA, sum( sel))
names(x) <- uENSP[ sel ]
any(duplicated(c(names(x), names(ensp2sym))))  # FALSE
ensp2sym <- c(ensp2sym, x)
all(uENSP %in% names(ensp2sym))  #TRUE

sel <- which(ensp2sym == "")
ensp2sym[head(sel)]
ensp2sym[sel] <- NA
ensp2sym[head(sel)]
all( uENSP %in% names(ensp2sym)) #TRUE
sum(is.na(ensp2sym)) #636


all(uENSP %in% names(ensp2sym)) #TRUE
sum(! is.na(ensp2sym)) #4835
sum(! is.na(ensp2sym)) * 100 / length(ensp2sym) #88.1%
all(ensp2sym[! is.na(ensp2sym)] %in% HGNC$sym) #TRUE

# Done.
# This concludes construction of our mapping tool.
# Save the map:

save(ensp2sym, file = file.path("inst", "extdata", "ensp2sym.RData"))

# From an RStudio project, the file can be loaded with
load(file = file.path("inst", "extdata", "ensp2sym.RData"))



#========================DATA STATISTICS========================================

tmp <- readr::read_tsv(file.path("../data", "homNOG.annotations.tsv"),
                       skip = 1, col_names = c("homNOG", "GeneID", "Protein", "Species", "Category", "ID"))

head(tmp)
# A tibble: 6 x 6
#homNOG GeneID      Protein Species Category ID
#<chr>  <chr>         <dbl>   <dbl> <chr>    <chr>
#1 homNOG ENOG410RR6N       4       4 L        9606.ENSP00000473172
#2 homNOG ENOG410S23Y       2       2 S        9606.ENSP00000473161
#3 homNOG ENOG410RYSR       4       3 S        9606.ENSP00000473139
#4 homNOG ENOG410S2EW       2       2 S        9606.ENSP00000473113
#5 homNOG ENOG410S1H2       3       2 S        9606.ENSP00000473105
#6 homNOG ENOG410RVFP       4       3 S        9606.ENSP00000473103

all(grepl("^9606\\.", tmp$ID)) #TRUE
tmp$ID <- gsub("^9606\\.", "", tmp$ID)

minProtein <- 0
maxProtein <- 20
hist(tmp$Protein[(tmp$Protein >= minProtein) & (tmp$Protein <= maxProtein)],
     xlim = c(minProtein, maxProtein),
     main = "EggNOG dataset Proteins",
     col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 2)(40),
     xlab = "Number of Proteins")

minSpecies <- 1
maxSpecies <- 4
hist(tmp$Species[(tmp$Species >= minSpecies) & (tmp$Species <= maxSpecies)],
     xlim = c(minSpecies, maxSpecies),
     main = "EggNOG dataset Species",
     col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 2)(40),
     xlab = "Number of Species")


tmp$ID <- ensp2sym[tmp$ID]
any(grepl("ENSP", tmp$ID)) #NONE
sum(is.na(tmp$ID)) #652
EggNOGedges <- tmp[( ! is.na(tmp$ID)), ]
save(EggNOGedges, file = file.path("..", "data", "STRINGedges.RData"))

