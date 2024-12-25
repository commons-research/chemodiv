# Folowing https://cran.uib.no/web/packages/chemodiv/vignettes/chemodiv-vignette.html


# Install current version
install.packages("chemodiv")

# Install developmental version
install.packages("devtools") # Install devtools if not already installed
library("devtools")
install_github("hpetren/chemodiv")

# Load chemodiv
library(chemodiv)

data("alpinaSampData")

head(alpinaSampData)[,1:5]


data("alpinaCompData")

head(alpinaCompData)

chemoDivCheck(compoundData = alpinaCompData, sampleData = alpinaSampData)


data("alpinaPopData")

table(alpinaPopData)


alpinaNPC <- NPCTable(compoundData = alpinaCompData)


alpinaCompDis <- compDis(compoundData = alpinaCompData,
                         type = "PubChemFingerprint")

alpinaCompDis$fingerDisMat[1:4, 1:4] # Part of compound dissimilarity matrix


alpinaDiv <- calcDiv(sampleData = alpinaSampData, 
                     compDisMat = alpinaCompDis$fingerDisMat,
                     type = "FuncHillDiv",
                     q = 1)

alpinaDivProf <- calcDivProf(sampleData = alpinaSampData,
                             compDisMat = alpinaCompDis$fingerDisMat,
                             type = "FuncHillDiv")

alpinaBetaDiv <- calcBetaDiv(sampleData = alpinaSampData,
                             compDisMat = alpinaCompDis$fingerDisMat,
                             type = "FuncHillDiv")



### Sample dissimilarities


alpinaSampDis <- sampDis(sampleData = alpinaSampData,
                         compDisMat = alpinaCompDis$fingerDisMat,
                         type = "GenUniFrac")

alpinaSampDis$GenUniFrac[1:4, 1:4] # Part of sample dissimilarity matrix


# Molecular Network

alpinaNetwork <- molNet(compDisMat = alpinaCompDis$fingerDisMat,
                        npcTable = alpinaNPC,
                        cutOff = 0.75)

summary(alpinaNetwork)


molNetPlot(sampleData = alpinaSampData,
           networkObject = alpinaNetwork$networkObject,
           npcTable = alpinaNPC,
           plotNames = TRUE)

chemoDivPlot(compDisMat = alpinaCompDis$fingerDisMat,
             divData = alpinaDiv,
             divProfData = alpinaDivProf,
             sampDisMat = alpinaSampDis$GenUniFrac,
             groupData = alpinaPopData)

?chemodiv


### Example
## Here we will work on a real case dataset 

# Load chemodiv
library(chemodiv)
library(dplyr)
library(tidyverse)


# For this we run met_annot_enhancer which should conveniently fetch a GNPS job and proceed to annotations so we have the required tables
# https://github.com/mandelbrot-project/met_annot_enhancer



data("alpinaSampData")

head(alpinaSampData)[,1:5]

gnps_job_id <- 'd4aca8b0f5d342488086b66b755dfd21'

SampData <- paste('docs', gnps_job_id, quantification_table_reformatted)

# The samples composition data is fetched from the GNPS job folder
SampData <- read.csv('docs/d4aca8b0f5d342488086b66b755dfd21/quantification_table_reformatted/5f16bbdf89aa4ff687b85139f47b9615.csv')

# Formatting to fit the chemodiv needs

colnames(SampData)
SampData <- SampData %>%
    select(-row.m.z, -row.retention.time, -last_col())


SampData <- as.data.frame(t(SampData))

# 
names(SampData) <- SampData[1,]
SampData <- SampData[-1,]

myrmecoSampData <- SampData

### We now load the Compounds data

# The compounds data is fetched from the met_annot_enhancer results folder
CompData <- read.csv('docs/mapp_batch_00019/mapp_batch_00019_spectral_match_results_repond_flat.tsv', sep = '\t')


CompData <- CompData %>%
    select(feature_id, structure_smiles, structure_inchikey) %>%
    rename(
        compound = feature_id,
        smiles = structure_smiles,
        inchikey = structure_inchikey
    )  %>% 
    mutate_all(na_if,"")  %>% 
    arrange(compound)  %>% 
    distinct(compound, .keep_all= TRUE)


CompData <- na.omit(CompData)


Comp <- data.frame(compound = colnames(SampData))

Comp_merged <- merge(x = Comp, y = CompData, by = "compound", all = TRUE )

myrmecoCompData <- Comp_merged %>% 
    arrange(as.numeric(compound))

chemoDivCheck(compoundData = myrmecoCompData, sampleData = myrmecoSampData)



data("alpinaPopData")

table(alpinaPopData)


myrmecoMetaData <- read.csv('docs/d4aca8b0f5d342488086b66b755dfd21/metadata_table/metadata_table-00000.tsv', sep = '\t')


myrmecoMetaData$ATTRIBUTE_colony

myrmecoNPC <- NPCTable(compoundData = myrmecoCompData)


myrmecoNPC_nona <- myrmecoNPC %>% drop_na(smiles, class)
myrmecoCompData_nona <- myrmecoCompData %>% drop_na(smiles)



myrmecoCompDis <- compDis(compoundData = myrmecoCompData_nona,
                         type = "NPClassifier",
                         npcTable = myrmecoNPC_nona)


myrmecoCompDis$fingerDisMat[1:4, 1:4] # Part of compound dissimilarity matrix

myrmecoCompDis$npcDisMat[1:50, 1:50]

myrmecoSampData_filt <- myrmecoSampData[, which((names(myrmecoSampData) %in% myrmecoNPC_nona$compound)==TRUE)]


myrmecoDiv <- calcDiv(sampleData = myrmecoSampData_filt, 
                     compDisMat = myrmecoCompDis$npcDisMat,
                     type = "FuncHillDiv",
                     q = 1)



myrmecoDivProf <- calcDivProf(
    sampleData = myrmecoSampData_filt,
    compDisMat = myrmecoCompDis$npcDisMat,
    type = "FuncHillDiv",
    qMax = 2,
    step = 1
)


myrmecoBetaDiv <- calcBetaDiv(sampleData = myrmecoSampData_filt,
                             compDisMat = myrmecoCompDis$npcDisMat,
                             type = "FuncHillDiv")



### Sample dissimilarities


myrmecoSampDis <- sampDis(sampleData = myrmecoSampData_filt,
                         compDisMat = myrmecoCompDis$npcDisMat,
                         type = "GenUniFrac")

myrmecoSampDis$GenUniFrac[1:4, 1:4] # Part of sample dissimilarity matrix


# Molecular Network

myrmecoNetwork <- molNet(compDisMat = myrmecoCompDis$npcDisMat,
                        npcTable = myrmecoNPC_nona,
                        cutOff = 0.75)

summary(myrmecoNetwork)

library(igraph)

write_graph(
  myrmecoNetwork,
  'test.graphml',
  format = "graphml"
)




molNetPlot(sampleData = myrmecoSampData_filt,
           networkObject = myrmecoNetwork$networkObject,
           npcTable = myrmecoNPC_nona,
           plotNames = FALSE)

chemoDivPlot(compDisMat = myrmecoCompDis$npcDisMat
            #  divData = myrmecoDiv,
            #  divProfData = myrmecoDivProf,
            #  sampDisMat = myrmecoSampDis$GenUniFrac,
            #  groupData = myrmecoMetaData$ATTRIBUTE_colony
             )

?chemodiv


##### Working on the DBGI Sapindales datset


library(chemodiv)
library(dplyr)
library(tidyverse)




# The samples composition data is fetched from the GNPS job folder
SampData <- read.csv('/Users/pma/Dropbox/git_repos/COMMONS_Lab/DBGI/dbgi-sapindales/docs/dbgi_sapindales/dbgi_batch_000001/pos/results/met_annot_enhancer/86b3b9488fbb46bf9aad708106ff0a27/quantification_table_reformatted/364a3ecd1ce44fa59a0c10ff5d0dbd62.csv')

# Formatting to fit the chemodiv needs

colnames(SampData)
SampData <- SampData %>%
    select(-row.m.z, -row.retention.time, -last_col())


SampData <- as.data.frame(t(SampData))

# 
names(SampData) <- SampData[1,]
SampData <- SampData[-1,]

sapindalesSampData <- SampData


### We now load the Compounds data

# The compounds data is fetched from the met_annot_enhancer results folder
CompData <- read.csv('docs/mapp_batch_00019/mapp_batch_00019_spectral_match_results_repond_flat.tsv', sep = '\t')
CompData <- read.csv('/Users/pma/Dropbox/git_repos/COMMONS_Lab/DBGI/dbgi-sapindales/docs/dbgi_sapindales/dbgi_batch_000001/pos/results/met_annot_enhancer/dbgi_batch_000001/dbgi_batch_000001_spectral_match_results_repond_flat.tsv', sep = '\t')


CompData <- CompData %>%
    select(feature_id, structure_smiles, structure_inchikey, structure_taxonomy_npclassifier_01pathway, structure_taxonomy_npclassifier_02superclass, structure_taxonomy_npclassifier_03class) %>%
    rename(
        compound = feature_id,
        smiles = structure_smiles,
        inchikey = structure_inchikey,
        pathway = structure_taxonomy_npclassifier_01pathway,
        superclass = structure_taxonomy_npclassifier_02superclass,
        class = structure_taxonomy_npclassifier_03class
    )  %>%     
    # here we change the datatype of the feature_id to character
    mutate_at(vars(compound), ~as.character(.))  %>%
    mutate_all(na_if,"")  %>% 
    arrange(compound)  %>% 
    distinct(compound, .keep_all= TRUE)


CompData <- na.omit(CompData)

glimpse(CompData)

Comp <- data.frame(compound = colnames(SampData))

Comp_merged <- merge(x = Comp, y = CompData, by = "compound", all = TRUE )

sapindalesCompData <- Comp_merged %>% 
    arrange(as.numeric(compound))

chemoDivCheck(compoundData = sapindalesCompData, sampleData = sapindalesSampData)


# myrmecoMetaData <- read.csv('docs/d4aca8b0f5d342488086b66b755dfd21/metadata_table/metadata_table-00000.tsv', sep = '\t')
sapindalesMetaData <- read.csv('/Users/pma/Dropbox/git_repos/COMMONS_Lab/DBGI/dbgi-sapindales/docs/dbgi_sapindales/dbgi_batch_000001/pos/results/met_annot_enhancer/86b3b9488fbb46bf9aad708106ff0a27/metadata_table/metadata_table-00000.txt', sep = '\t')

# myrmecoMetaData$ATTRIBUTE_colony
# myrmecoNPC <- NPCTable(compoundData = myrmecoCompData)
# head(myrmecoNPC)
# glimpse(myrmecoNPC)


# myrmecoNPC_nona <- myrmecoNPC %>% drop_na(smiles, class)
# myrmecoCompData_nona <- myrmecoCompData %>% drop_na(smiles)

sapindalesNPC = sapindalesCompData

sapindalesNPC_nona <- sapindalesNPC %>% drop_na(smiles, class)

sapindalesCompData_nona <- sapindalesCompData %>% drop_na(smiles)



sapindalesCompDis <- compDis(compoundData = sapindalesCompData_nona,
                         type = "NPClassifier",
                         npcTable = sapindalesNPC_nona)


sapindalesCompDis$npcDisMat[1:50, 1:50]



sapindalesSampData_filt <- sapindalesSampData[, which((names(sapindalesSampData) %in% sapindalesNPC_nona$compound)==TRUE)]


sapindalesDiv <- calcDiv(sampleData = sapindalesSampData_filt, 
                     compDisMat = sapindalesCompDis$npcDisMat,
                     type = "FuncHillDiv",
                     q = 1)



myrmecoDivProf <- calcDivProf(
    sampleData = myrmecoSampData_filt,
    compDisMat = myrmecoCompDis$npcDisMat,
    type = "FuncHillDiv",
    qMax = 2,
    step = 1
)


myrmecoBetaDiv <- calcBetaDiv(sampleData = myrmecoSampData_filt,
                             compDisMat = myrmecoCompDis$npcDisMat,
                             type = "FuncHillDiv")



### Sample dissimilarities


myrmecoSampDis <- sampDis(sampleData = myrmecoSampData_filt,
                         compDisMat = myrmecoCompDis$npcDisMat,
                         type = "GenUniFrac")

myrmecoSampDis$GenUniFrac[1:4, 1:4] # Part of sample dissimilarity matrix


# Molecular Network

myrmecoNetwork <- molNet(compDisMat = myrmecoCompDis$npcDisMat,
                        npcTable = myrmecoNPC_nona,
                        cutOff = 0.75)

summary(myrmecoNetwork)

library(igraph)

write_graph(
  myrmecoNetwork,
  'test.graphml',
  format = "graphml"
)




molNetPlot(sampleData = myrmecoSampData_filt,
           networkObject = myrmecoNetwork$networkObject,
           npcTable = myrmecoNPC_nona,
           plotNames = FALSE)

chemoDivPlot(compDisMat = myrmecoCompDis$npcDisMat
            #  divData = myrmecoDiv,
            #  divProfData = myrmecoDivProf,
            #  sampDisMat = myrmecoSampDis$GenUniFrac,
            #  groupData = myrmecoMetaData$ATTRIBUTE_colony
             )

?chemodiv