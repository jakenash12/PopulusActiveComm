#WDs tests are done to confirm results of PERMANOVA
#PERMANOVA is sensitive to differences in beta dispersion
#WDs tests are less sensitive to differences in beta
#dispersion and are used as validation

#This script relies on the distance matrices and
#metadata files that were created in the PERMANOVA.R
#script. Make sure to run that first

#load Mic eco package for wds test function
library("MicEco")


#effect of site type - 16S DNA dataset
set.seed(72)
WdS.test(AUE2021_16S_rarefied_DNA_dm_bc, 
         as.factor(AUE2021_16S_rarefied_DNA_dm_metadata$SiteType), n=999)

#effect of season - 16S DNA dataset
set.seed(72)
WdS.test(AUE2021_16S_rarefied_DNA_dm_bc, 
         as.factor(AUE2021_16S_rarefied_DNA_dm_metadata$Time), n=999)

#effect of site type - 16S RNA dataset
set.seed(72)
WdS.test(AUE2021_16S_rarefied_RNA_dm_bc, 
         as.factor(AUE2021_16S_rarefied_RNA_dm_metadata$SiteType), n=999)

#effect of season - 16S RNA dataset
set.seed(72)
WdS.test(AUE2021_16S_rarefied_RNA_dm_bc, 
         as.factor(AUE2021_16S_rarefied_RNA_dm_metadata$Time), n=999)

#effect of site type - ITS DNA
set.seed(72)
WdS.test(AUE2021_ITS_rarefied_DNA_dm_bc, 
         as.factor(AUE2021_ITS_rarefied_DNA_dm_metadata$SiteType), n=9999)

#effect of season - ITS DNA
set.seed(72)
WdS.test(AUE2021_ITS_rarefied_DNA_dm_bc, 
         as.factor(AUE2021_ITS_rarefied_DNA_dm_metadata$Time), n=9999)

#effect of site type - ITS RNA
set.seed(72)
WdS.test(AUE2021_ITS_rarefied_RNA_dm_bc, 
         as.factor(AUE2021_ITS_rarefied_RNA_dm_metadata$SiteType), n=9999)

#effect of season - ITS RNA
set.seed(72)
WdS.test(AUE2021_ITS_rarefied_RNA_dm_bc, 
         as.factor(AUE2021_ITS_rarefied_RNA_dm_metadata$Time), n=9999)
