# BCB420.2019.EGGNOG

EggNOG is a database of orthologous groups and functional annotation

## 1 About this package:

This package describes the workflow to map the [the EggNOG database](http://eggnogdb.embl.de/#/app/home) IDs to [HGNC](https://www.genenames.org/) symbols and it provides examples of computing database statistics.The major part of the code is concerned with building an accurate ID map between ENSP identifiers and HGNC symbols.

The package serves dual duty, as an RStudio project, as well as an R package that can be installed. Package checks **pass without errors, warnings, or notes**.

&nbsp;

## 2 Data download and cleanup

To download the source data from EggNOG ... :

1. Navigate to the [**EggNOG** database](http://eggnogdb.embl.de/#/app/hom) and follow the link to the [download section](http://eggnogdb.embl.de/#/app/downloads).
2. Choose "Hominidae" as organism.
3. Download the data file: (Warning: large).
4. Place files into a sister directory of your working directory which is called `data`. (It should be reachable with `file.path("..", "data")`)

&nbsp;

## 3 Mapping ENSEMBL IDs to HGNC symbols

#### Preparations: packages, functions, files

To begin, we need to make sure the required packages are installed:

**`readr`** provides functions to read data which are particularly suitable for
large datasets. They are much faster than the built-in read.csv() etc. But caution: these functions return "tibbles", not data frames. ([Know the difference](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html).)
```R
if (! requireNamespace("readr")) {
  install.packages("readr")
}
```

**`biomaRt`** biomaRt is a Bioconductor package that implements the RESTful API of biomart,
the annotation framwork for model organism genomes at the EBI. It is a Bioconductor package, and as such it needs to be loaded via the **`BiocManager`**,
&nbsp;

```R
if (! requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (! requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}
```

**`igraph`** is THE go-to package for everything graph related. We use it here to
compute some statistics on the STRING- and example graphs and plot.
&nbsp;

```R
if (! requireNamespace("igraph")) {
  install.packages("igraph")
}
```

&nbsp;

#### 3.1 Step one: which IDs do we have to map?

&nbsp;

```R
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
```

&nbsp;

#### 3.2  Step two: mapping via biomaRt

To proceed with the mapping, we use biomaRt to fetch as many HGNC symbols as we can - first in bulk (mapping ENSP IDs to HGNC symbols), then individually for the remaining IDs we could not map.

&nbsp;

###### 3.2.1  Constructing an ID-mapping tool

```R

# Map ENSP to HGNC symbols: open a "Mart" object ..
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
```

&nbsp;

There are three possible problems with the data that biomart returns:

&nbsp;

**(1)** There might be more than one value returned. The ID appears more than
once in `tmp$ensembl_peptide_id`, with different mapped symbols.

```R
  sum(duplicated(tmp$ensembl_peptide_id)) #1 duplicate
```

&nbsp;

**(2)** There might be nothing returned for one ENSP ID. We have the ID in `uENSP`, but it does not appear in `tmp$ensembl_peptide_id`:

```R
  sum(! (uENSP) %in% tmp$ensembl_peptide_id) #616 
```
&nbsp;

**(3)** There might be no value returned: `NA`, or `""`. The ID appears in `tmp$ensembl_peptide_id`, but there is no symbol in `tmp$hgnc_symbol`.

```R
  sum(is.na(tmp$ensembl_peptide_id)) #0
  sum(tmp$ensembl_peptide_id == "") #0 None of the rows have ""
```

&nbsp;

Let's fix the "duplicates" problem.

&nbsp;

```R
dupEnsp <- tmp$ensembl_peptide_id[duplicated(tmp$ensembl_peptide_id)]
tmp[tmp$ensembl_peptide_id %in% dupEnsp, ]
#ensembl_peptide_id hgnc_symbol
#2268    ENSP00000344961     PLEKHG7
#2269    ENSP00000344961    C12orf74

tmp[tmp$hgnc_symbol %in% c("C12orf74"), ]
tmp <- tmp[ ! (tmp$hgnc_symbol %in% c("C12orf74")), ]
any(duplicated(tmp$ensembl_peptide_id)) #FALSE
```

&nbsp;

After this preliminary cleanup, defining the mapping tool is simple:

&nbsp;

```R
  ensp2sym <- tmp$hgnc_symbol
  names(ensp2sym) <- tmp$ensembl_peptide_id
  head(ensp2sym)
  #ENSP00000430236 ENSP00000430633 ENSP00000447378 ENSP00000446058 ENSP00000468977
  #    "ALG11"        "NPIPB5"          "MGAM"         "SHD"         "FLT3LG"

```

&nbsp;

###### 3.2.2  Cleanup and validation of `ensp2sym`

First, we add the symbols that were not returned by biomaRt to the map. They are present in uENSP, but not in ensp2sym$ensp:

&nbsp;

```R
  sel <- ! (uENSP %in% names(ensp2sym))
  x <- rep(NA, sum( sel))
  names(x) <- uENSP[ sel ]
  
  # confirm uniqueness
  any(duplicated(c(names(x), names(ensp2sym))))  # FALSE
  
  ensp2sym <- c(ensp2sym, x) # concatenate the two vectors
  
  # confirm
  all(uENSP %in% names(ensp2sym))  #TRUE
```

&nbsp;

Next, we set the symbols for which only an empty string was returned to `NA`:

&nbsp;

```R
  sel <- which(ensp2sym == "")
  ensp2sym[head(sel)] #before
  ensp2sym[sel] <- NA
  ensp2sym[head(sel)] #after
  all( uENSP %in% names(ensp2sym)) #TRUE
  sum(is.na(ensp2sym)) #636
```

&nbsp;

#### 3.3 Final validation

Validation and statistics of our mapping tool:

```R

# do we now have all ENSP IDs mapped?
all(uENSP %in% names(ensp2sym)) #TRUE

# how many symbols did we find?
sum(! is.na(ensp2sym)) #4835

# (in %)
sum(! is.na(ensp2sym)) * 100 / length(ensp2sym) #88.1%

# are all symbols current in our reference table?
all(ensp2sym[! is.na(ensp2sym)] %in% HGNC$sym) #TRUE

# Done.
# This concludes construction of our mapping tool.
# Save the map:

save(ensp2sym, file = file.path("inst", "extdata", "ensp2sym.RData"))

# From an RStudio project, the file can be loaded with
load(file = file.path("inst", "extdata", "ensp2sym.RData"))
```

&nbsp;

# 4 Annotating gene sets with EggNOG Data

Given our mapping tool, we can now annotate gene sets with EggNOG data. 

&nbsp;

```R
tmp <- readr::read_tsv(file.path("../data", "homNOG.annotations.tsv"),
                       skip = 1, 
                       col_names = c("homNOG", "GeneID", "Protein", "Species", "Category", "ID"))

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

# do they all have the right tax id?
all(grepl("^9606\\.", tmp$ID)) #TRUE
# remove "9606." prefix
tmp$ID <- gsub("^9606\\.", "", tmp$ID)

minProtein <- 0
maxProtein <- 20
hist(tmp$Protein[(tmp$Protein >= minProtein) & (tmp$Protein <= maxProtein)],
     xlim = c(minProtein, maxProtein),
     main = "EggNOG dataset Proteins",
     col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 2)(40),
     xlab = "Number of Proteins")

```

![](./inst/img/EggNOG_Proteins.jpg?sanitize=true "EggNOG dataset Proteins")

```R
minSpecies <- 1
maxSpecies <- 4
hist(tmp$Species[(tmp$Species >= minSpecies) & (tmp$Species <= maxSpecies)],
     xlim = c(minSpecies, maxSpecies),
     main = "EggNOG dataset Species",
     col = colorRampPalette(c("#FFFFFF","#8888A6","#FF6655"), bias = 2)(40),
     xlab = "Number of Species")

```

![](./inst/img/EggNOG_Species.jpg?sanitize=true "EggNOG dataset Species")

&nbsp;

Finally we map the ENSP IDs to HGNC symbols. Using our tool, this is a simple assignment:

&nbsp;

```R

tmp$ID <- ensp2sym[tmp$ID]

# Validate:
# how many rows could not be mapped
any(grepl("ENSP", tmp$ID)) #NONE
sum(is.na(tmp$ID)) #652

# we remove edges in which either one or the other node is NA to
# create our final data:
EggNOGedges <- tmp[( ! is.na(tmp$ID)), ]

# Done.
# Save result
save(EggNOGedges, file = file.path("..", "data", "STRINGedges.RData"))

```

&nbsp;

## 5 Annotation of the example gene set

To conclude, we annotate the example gene set, validate the annotation, and store the data.
&nbsp;

```R

# The specification of the sample set is copy-paste from the 
# BCB420 resources project.

xSet <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2",
          "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7",
          "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP",
          "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM",
          "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B",
          "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN",
          "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21",
          "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN",
          "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6",
          "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1",
          "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7",
          "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39",
          "VPS41", "VTI1B", "YKT6")

# which example genes are not among the known nodes?
x <- which( ! (xSet %in% c(EggNOGedges$ID)))
cat(sprintf("\t%s\t(%s)\n", HGNC[xSet[x], "sym"], HGNC[xSet[x], "name"]))

#ATP2A1	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 1)
#ATP2A2	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 2)
#ATP2A3	(ATPase sarcoplasmic/endoplasmic reticulum Ca2+ transporting 3)
#BECN2	(beclin 2)
#BLOC1S1	(biogenesis of lysosomal organelles complex 1 subunit 1)
#BORCS6	(BLOC-1 related complex subunit 6)
#BORCS7	(BLOC-1 related complex subunit 7)
#BORCS8	(BLOC-1 related complex subunit 8)
#CACNA1A	(calcium voltage-gated channel subunit alpha1 A)
#CALCOCO2	(calcium binding and coiled-coil domain 2)
#CTTN	(cortactin)
#DCTN1	(dynactin subunit 1)
#GABARAP	(GABA type A receptor-associated protein)
#GABARAPL1	(GABA type A receptor associated protein like 1)
#GABARAPL2	(GABA type A receptor associated protein like 2)
#HDAC6	(histone deacetylase 6)
#HSPB8	(heat shock protein family B (small) member 8)
#INPP5E	(inositol polyphosphate-5-phosphatase E)
#KXD1	(KxDL motif containing 1)
#LAMP2	(lysosomal associated membrane protein 2)
#MAP1LC3C	(microtubule associated protein 1 light chain 3 gamma)
#MGRN1	(mahogunin ring finger 1)
#MYO1C	(myosin IC)
#MYO6	(myosin VI)
#NAPA	(NSF attachment protein alpha)
#OPTN	(optineurin)
#OSBPL1A	(oxysterol binding protein like 1A)
#PI4K2A	(phosphatidylinositol 4-kinase type 2 alpha)
#PIK3C3	(phosphatidylinositol 3-kinase catalytic subunit type 3)
#PLEKHM1	(pleckstrin homology and RUN domain containing M1)
#PSEN1	(presenilin 1)
#RAB21	(RAB21, member RAS oncogene family)
#RAB34	(RAB34, member RAS oncogene family)
#RAB39A	(RAB39A, member RAS oncogene family)
#RAB7A	(RAB7A, member RAS oncogene family)
#RAB7B	(RAB7B, member RAS oncogene family)
#RPTOR	(regulatory associated protein of MTOR complex 1)
#RUBCN	(rubicon autophagy regulator)
#SNAP47	(synaptosome associated protein 47)
#SNAPIN	(SNAP associated protein)
#STX17	(syntaxin 17)
#SYT7	(synaptotagmin 7)
#TARDBP	(TAR DNA binding protein)
#TFEB	(transcription factor EB)
#TGM2	(transglutaminase 2)
#TIFA	(TRAF interacting protein with forkhead associated domain)
#TMEM175	(transmembrane protein 175)
#TPCN1	(two pore segment channel 1)
#TPCN2	(two pore segment channel 2)
#TPPP	(tubulin polymerization promoting protein)
#TXNIP	(thioredoxin interacting protein)
#VAMP7	(vesicle associated membrane protein 7)
#VAPA	(VAMP associated protein A)
#VPS11	(VPS11, CORVET/HOPS core subunit)
#VPS16	(VPS16, CORVET/HOPS core subunit)
#VPS18	(VPS18, CORVET/HOPS core subunit)
#VPS41	(VPS41, HOPS complex subunit)
#YKT6	(YKT6 v-SNARE homolog)



sel <- (EggNOGedges$ID %in% xSet)
xSetEdges <- EggNOGedges[sel, c("ID")]
# Statistics:
nrow(xSetEdges)   # 27

# Save the annotated set

# Save the annotated set
writeLines(c("ID",
             sprintf("%s", xSetEdges$ID)),
           con = "xSetEdges.tsv")
```

&nbsp;

## 6 References

&nbsp;

This package script was created based on and part of code to create the mapping tool was taken from the package [BCB420.2019.STRING](https://github.com/hyginn/BCB420.2019.STRING) 

The functionality and code of the toBrowser.R file was taken from Prof. Boris's [BCB course resources page](https://github.com/hyginn/BCB420-2019-resources)

&nbsp;

<!-- END -->
