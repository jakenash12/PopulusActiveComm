#Applies the LinDA test for differential abundance based on Zhou et al 2022
#Does LINDA at Phylum/Class level for bacteria
#(summarises data to phylum level for everything except
#Proteobacteria, which is done at Class level)
#
#Then does LINDA at Class level for fungi

library("MicrobiomeStat")
library(dplyr)
library(magrittr)
library(ggplot2)
library(tibble)


lfc_cutoff_RvD=1
p_cutoff_RvD=0.05

#generates taxonomic summary tables from rarefied OTU tables

###########16S summary#################3#################################
AUE2021_16S_rarefied_phylum=
  AUE2021_16S_rarefied %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(tax_mat_16S_df) %>% 
  mutate(Phylum_edit=case_when(
    (Phylum=="Proteobacteria" & !(is.na(Class))) ~ Class,
    (Phylum=="Proteobacteria" & (is.na(Class))) ~ "OtherProteobacteria",
    is.na(Phylum) ~"Unassigned",
    .default=Phylum
  )) %>%
  group_by(Phylum_edit) %>%
  summarise_at(vars(matches('DNA')), sum) %>%  
  column_to_rownames("Phylum_edit")

Linda_16S_phylum_df_meta=
  filter(samples_df2_16S, sample %in% colnames(AUE2021_16S_rarefied_phylum))

Linda_16S_phylum_df_meta_DNA=
  samples_df2_16S %>%
  filter(source=="DNA",
         sample %in% colnames(AUE2021_16S_rarefied_phylum))

Linda_16S_phylum_df_DNA=
  AUE2021_16S_rarefied_phylum %>%
  dplyr::select(Linda_16S_phylum_df_meta_DNA$sample)

Linda_16S_phylum_df_meta_RNA=
  samples_df2_16S %>%
  filter(source=="cDNA",
         sample %in% colnames(AUE2021_16S_rarefied_phylum))

Linda_16S_phylum_df_RNA=
  AUE2021_16S_rarefied_phylum %>%
  dplyr::select(Linda_16S_phylum_df_meta_RNA$sample)


###########ITS Summary Tables###################################
AUE2021_ITS_rarefied_class=
  AUE2021_ITS_rarefied %>%
  otu_table %>%
  as.data.frame %>%
  mutate(otu=rownames(.)) %>%
  left_join(tax_mat_df) %>% 
  mutate(Class_edit=
           case_when(!(is.na(Class)) ~ Class,
                     .default="Unassigned")) %>%
  group_by(Class_edit) %>%
  summarise_at(vars(matches('DNA')), sum) %>%  
  column_to_rownames("Class_edit")

Linda_ITS_class_df_meta=
  filter(samples_df2, sample %in% colnames(AUE2021_ITS_rarefied_class))

Linda_ITS_class_df_meta_DNA=
  samples_df2 %>%
  filter(source=="DNA",
         sample %in% colnames(AUE2021_ITS_rarefied_class))

Linda_ITS_class_df_DNA=
  AUE2021_ITS_rarefied_class %>%
  dplyr::select(Linda_ITS_class_df_meta_DNA$sample)

Linda_ITS_class_df_meta_RNA=
  samples_df2 %>%
  filter(source=="cDNA",
         sample %in% colnames(AUE2021_ITS_rarefied_class))

Linda_ITS_class_df_RNA=
  AUE2021_ITS_rarefied_class %>%
  dplyr::select(Linda_ITS_class_df_meta_RNA$sample)

############################16S Data - RNA vs DNA ##################
Linda_16S_phylum_source =
  linda(
    feature.dat=AUE2021_16S_rarefied_phylum,
    meta.dat=Linda_16S_phylum_df_meta,
    formula='~SiteType+Time+source + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

DNAvsRNA_16S_Phylum_sig=
  Linda_16S_phylum_source$output$sourceDNA %>%
  mutate(Phylum=rownames(.)) %>%
  mutate(Sig=ifelse((padj<p_cutoff_RvD&log2FoldChange>lfc_cutoff_RvD),"DNA Enriched",ifelse((padj<p_cutoff_RvD&log2FoldChange< -lfc_cutoff_RvD),"RNA Enriched","Non Significant")))

############################ITS Data - RNA vs DNA ##################
Linda_ITS_class_source =
  linda(
    feature.dat=AUE2021_ITS_rarefied_class,
    meta.dat=Linda_ITS_class_df_meta,
    formula='~SiteType+Time+source + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

DNAvsRNA_ITS_class_sig=
  Linda_ITS_class_source$output$sourceDNA %>%
  mutate(Class=rownames(.)) %>%
  mutate(Sig=ifelse((padj<p_cutoff_RvD&log2FoldChange>lfc_cutoff_RvD),"DNA Enriched",ifelse((padj<p_cutoff_RvD&log2FoldChange< -lfc_cutoff_RvD),"RNA Enriched","Non Significant")))


######################Testing for effects of site type and Season###########################
#a log2fold change of 0 and an adjusted p-value of 0.05 are used to 
#determine significance for site type and season
lfc_cutoff_SS=0
p_cutoff_SS=0.05

##########################16S Data - DNA#########################

Linda_16S_Phylum_SiteType_DNA =
  linda(
    feature.dat=Linda_16S_phylum_df_DNA,
    meta.dat=mutate(Linda_16S_phylum_df_meta_DNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_16S_Phylum_SiteType_mid_June_DNA =
  linda(
    feature.dat=Linda_16S_phylum_df_DNA,
    meta.dat=mutate(Linda_16S_phylum_df_meta_DNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_DNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_mid_June_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_DNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_16S_Phylum_sig_DNA=
  Linda_16S_Phylum_SiteType_mid_June_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_16S_Phylum_variable_DNA=
  c(HvR_16S_Phylum_sig_DNA,HvM_16S_Phylum_sig_DNA,MvR_16S_Phylum_sig_DNA) %>%
  unique %>%
  data.frame(Phylum=.)

Time_16S_Phylum_variable_DNA=
  c(AvJ_16S_Phylum_sig_DNA,AvO_16S_Phylum_sig_DNA,JvO_16S_Phylum_sig_DNA) %>%
  unique %>%
  data.frame(Phylum=.)


##########################16S Data - RNA#########################

Linda_16S_Phylum_SiteType_RNA =
  linda(
    feature.dat=Linda_16S_phylum_df_RNA,
    meta.dat=mutate(Linda_16S_phylum_df_meta_RNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_16S_Phylum_SiteType_mid_June_RNA =
  linda(
    feature.dat=Linda_16S_phylum_df_RNA,
    meta.dat=mutate(Linda_16S_phylum_df_meta_RNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_RNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_mid_June_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_RNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_16S_Phylum_sig_RNA=
  Linda_16S_Phylum_SiteType_mid_June_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_16S_Phylum_variable_RNA=
  c(HvR_16S_Phylum_sig_RNA,HvM_16S_Phylum_sig_RNA,MvR_16S_Phylum_sig_RNA) %>%
  unique %>%
  data.frame(Phylum=.)

Time_16S_Phylum_variable_RNA=
  c(AvJ_16S_Phylum_sig_RNA,AvO_16S_Phylum_sig_RNA,JvO_16S_Phylum_sig_RNA) %>%
  unique %>%
  data.frame(Phylum=.)


##########################ITS Data - DNA#########################
Linda_ITS_class_SiteType_DNA =
  linda(
    feature.dat=Linda_ITS_class_df_DNA,
    meta.dat=mutate(Linda_ITS_class_df_meta_DNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_ITS_class_SiteType_mid_June_DNA =
  linda(
    feature.dat=Linda_ITS_class_df_DNA,
    meta.dat=mutate(Linda_ITS_class_df_meta_DNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_DNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_mid_June_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_DNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_ITS_class_sig_DNA=
  Linda_ITS_class_SiteType_mid_June_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_ITS_class_variable_DNA=
  c(HvR_ITS_class_sig_DNA,HvM_ITS_class_sig_DNA,MvR_ITS_class_sig_DNA) %>%
  unique %>%
  data.frame(class=.)

Time_ITS_class_variable_DNA=
  c(AvJ_ITS_class_sig_DNA,AvO_ITS_class_sig_DNA,JvO_ITS_class_sig_DNA) %>%
  unique %>%
  data.frame(class=.)

##########################ITS Data - RNA#########################
Linda_ITS_class_SiteType_RNA =
  linda(
    feature.dat=Linda_ITS_class_df_RNA,
    meta.dat=mutate(Linda_ITS_class_df_meta_RNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_ITS_class_SiteType_mid_June_RNA =
  linda(
    feature.dat=Linda_ITS_class_df_RNA,
    meta.dat=mutate(Linda_ITS_class_df_meta_RNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_RNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_mid_June_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_RNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_ITS_class_sig_RNA=
  Linda_ITS_class_SiteType_mid_June_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_ITS_class_variable_RNA=
  c(HvR_ITS_class_sig_RNA,HvM_ITS_class_sig_RNA,MvR_ITS_class_sig_RNA) %>%
  unique %>%
  data.frame(class=.)

Time_ITS_class_variable_RNA=
  c(AvJ_ITS_class_sig_RNA,AvO_ITS_class_sig_RNA,JvO_ITS_class_sig_RNA) %>%
  unique %>%
  data.frame(class=.)

