#Applies the LinDA test for differential abundance based on Zhou et al 2022
#This is used as a comparison for ANCOM because LinDA allows for 
#mixed model implementation

#Run the scripts ITS2_Analysis.R and 16S_Analysis.R prior to this script to generate the
#required data objects

library("MicrobiomeStat")
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggvenn)
library(cowplot)

#a log2fold change of 1 and an adjusted p-value of 0.05 are used to 
#determine significance for comparing RNA vs DNA datasets
lfc_cutoff_RvD=1
p_cutoff_RvD=0.05

#generates OTU tables and metadata formatted for Linda
Linda_ITS_df=as.data.frame(t(AUE2021_ITS_rarefied_df))
Linda_ITS_meta=filter(samples_df2, sample %in% rownames(AUE2021_ITS_rarefied_df))

Linda_ITS_df_DNA=as.data.frame(t(AUE2021_ITS_rarefied_df_DNA))
Linda_ITS_meta_DNA=filter(samples_df2, sample %in% rownames(AUE2021_ITS_rarefied_df_DNA))

Linda_ITS_df_RNA=as.data.frame(t(AUE2021_ITS_rarefied_df_RNA))
Linda_ITS_meta_RNA=filter(samples_df2, sample %in% rownames(AUE2021_ITS_rarefied_df_RNA))

Linda_16S_df=as.data.frame(t(AUE2021_16S_rarefied_df))
Linda_16S_meta=filter(samples_df2_16S, sample %in% rownames(AUE2021_16S_rarefied_df))

Linda_16S_df_DNA=as.data.frame(t(AUE2021_16S_rarefied_df_DNA))
Linda_16S_meta_DNA=filter(samples_df2_16S, sample %in% rownames(AUE2021_16S_rarefied_df_DNA))

Linda_16S_df_RNA=as.data.frame(t(AUE2021_16S_rarefied_df_RNA))
Linda_16S_meta_RNA=filter(samples_df2_16S, sample %in% rownames(AUE2021_16S_rarefied_df_RNA))

############################ITS Data - RNA vs DNA ##################
Linda_ITS_source =
  linda(
    feature.dat=Linda_ITS_df,
    meta.dat=Linda_ITS_meta,
    formula='~SiteType+Time+source + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

#creates df with a column titled "sig" indicating significance of test.
DNAvsRNA_ITS_sig=
  tax_mat_df %>%
  left_join(Linda_ITS_source$output$sourceDNA %>% mutate(otu=rownames(.))) %>%
  mutate(Sig=ifelse((padj<p_cutoff_RvD&log2FoldChange>lfc_cutoff_RvD),"DNA Enriched",ifelse((padj<p_cutoff_RvD&log2FoldChange< -(lfc_cutoff_RvD)),"RNA Enriched","Non Significant")))

############################16S Data - RNA vs DNA ##################
Linda_16S_source =
  linda(
    feature.dat=Linda_16S_df,
    meta.dat=Linda_16S_meta,
    formula='~SiteType+Time+source + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

#creates df with a column titled "sig" indicating significance of test.
#a log2fold change of 1 and an adjusted p-value of 0.01 are used to determine significance
DNAvsRNA_16S_sig=
  silva_tax_df %>%
  left_join(Linda_16S_source$output$sourceDNA %>% mutate(otu=rownames(.)))
  mutate(Sig=ifelse((padj<p_cutoff_RvD&log2FoldChange>lfc_cutoff_RvD),"DNA Enriched",ifelse((padj<p_cutoff_RvD&log2FoldChange< (-lfc_cutoff_RvD)),"RNA Enriched","Non Significant")))


######################Testing for effects of site type and Season###########################
#a log2fold change of 0 and an adjusted p-value of 0.05 are used to 
#determine significance for site type and season
lfc_cutoff_SS=0
p_cutoff_SS=0.05

############################ITS Data - DNA ##################

Linda_ITS_SiteType_DNA =
  linda(
    feature.dat=Linda_ITS_df_DNA,
    meta.dat=mutate(Linda_ITS_meta_DNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_ITS_SiteType_mid_June_DNA =
  linda(
    feature.dat=Linda_ITS_df_DNA,
    meta.dat=mutate(Linda_ITS_meta_DNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_ITS_sig_DNA=
  Linda_ITS_SiteType_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_ITS_sig_DNA=
  Linda_ITS_SiteType_DNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_ITS_sig_DNA=
  Linda_ITS_SiteType_mid_June_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_ITS_sig_DNA=
  Linda_ITS_SiteType_DNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_ITS_sig_DNA=
  Linda_ITS_SiteType_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_ITS_sig_DNA=
  Linda_ITS_SiteType_mid_June_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_ITS_variable_DNA=
  c(HvR_ITS_sig_DNA,HvM_ITS_sig_DNA,MvR_ITS_sig_DNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(tax_mat_df)
  
Time_ITS_variable_DNA=
  c(AvJ_ITS_sig_DNA,AvO_ITS_sig_DNA,JvO_ITS_sig_DNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(tax_mat_df)

############################ITS Data - RNA ##################

Linda_ITS_SiteType_RNA =
  linda(
    feature.dat=Linda_ITS_df_RNA,
    meta.dat=mutate(Linda_ITS_meta_RNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_ITS_SiteType_mid_June_RNA =
  linda(
    feature.dat=Linda_ITS_df_RNA,
    meta.dat=mutate(Linda_ITS_meta_RNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_ITS_sig_RNA=
  Linda_ITS_SiteType_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_ITS_sig_RNA=
  Linda_ITS_SiteType_RNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_ITS_sig_RNA=
  Linda_ITS_SiteType_mid_June_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_ITS_sig_RNA=
  Linda_ITS_SiteType_RNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_ITS_sig_RNA=
  Linda_ITS_SiteType_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_ITS_sig_RNA=
  Linda_ITS_SiteType_mid_June_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_ITS_variable_RNA=
  c(HvR_ITS_sig_RNA,HvM_ITS_sig_RNA,MvR_ITS_sig_RNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(tax_mat_df)

Time_ITS_variable_RNA=
  c(AvJ_ITS_sig_RNA,AvO_ITS_sig_RNA,JvO_ITS_sig_RNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(tax_mat_df)


##########################16S Data - DNA#########################

#creates a data matrix with DNA only
AUE2021_16S_rarefied_DNA=
  AUE2021_16S_rarefied %>%
  subset_samples(source=="DNA") %>%
  prune_taxa(taxa_sums(.) > 0, .) 

sample_data(AUE2021_16S_rarefied_DNA)$Time=factor(sample_data(AUE2021_16S_rarefied_DNA)$Time, levels=c("August", "June","October"))

Linda_16S_SiteType_DNA =
  linda(
    feature.dat=Linda_16S_df_DNA,
    meta.dat=mutate(Linda_16S_meta_DNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_16S_SiteType_mid_June_DNA =
  linda(
    feature.dat=Linda_16S_df_DNA,
    meta.dat=mutate(Linda_16S_meta_DNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_16S_sig=
  Linda_16S_SiteType_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_16S_sig=
  Linda_16S_SiteType_DNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_16S_sig=
  Linda_16S_SiteType_mid_June_DNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_16S_sig=
  Linda_16S_SiteType_DNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_16S_sig=
  Linda_16S_SiteType_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_16S_sig=
  Linda_16S_SiteType_mid_June_DNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_16S_variable_DNA=
  c(HvR_16S_sig,HvM_16S_sig,MvR_16S_sig) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(silva_tax_df)

Time_16S_variable_DNA=
  c(AvJ_16S_sig,AvO_16S_sig,JvO_16S_sig) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(silva_tax_df)


############################16S Data - RNA ##################

Linda_16S_SiteType_RNA =
  linda(
    feature.dat=Linda_16S_df_RNA,
    meta.dat=mutate(Linda_16S_meta_RNA,
                    SiteType=factor(SiteType, levels=c("High Elevation", "Mid Elevation", "Riparian")),
                    Time=factor(Time, levels=c("August","June", "October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

Linda_16S_SiteType_mid_June_RNA =
  linda(
    feature.dat=Linda_16S_df_RNA,
    meta.dat=mutate(Linda_16S_meta_RNA,
                    SiteType= factor(SiteType, levels=c("Mid Elevation", "Riparian","High Elevation")),
                    Time=factor(Time, levels=c("June", "August","October"))),
    formula='~SiteType+Time + (1|Site)',
    feature.dat.type="count",
    prev.filter=0.1,
    zero.handling="pseudo.cnt",
    pseudo.cnt = 0.5,
    p.adj.method="fdr",
    alpha=0.05)

HvR_16S_sig_RNA=
  Linda_16S_SiteType_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

HvM_16S_sig_RNA=
  Linda_16S_SiteType_RNA$output$`SiteTypeMid Elevation` %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

MvR_16S_sig_RNA=
  Linda_16S_SiteType_mid_June_RNA$output$SiteTypeRiparian %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvJ_16S_sig_RNA=
  Linda_16S_SiteType_RNA$output$TimeJune %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

AvO_16S_sig_RNA=
  Linda_16S_SiteType_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

JvO_16S_sig_RNA=
  Linda_16S_SiteType_mid_June_RNA$output$TimeOctober %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  rownames()

SiteType_16S_variable_RNA=
  c(HvR_16S_sig_RNA,HvM_16S_sig_RNA,MvR_16S_sig_RNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(silva_tax_df)

Time_16S_variable_RNA=
  c(AvJ_16S_sig_RNA,AvO_16S_sig_RNA,JvO_16S_sig_RNA) %>%
  unique %>%
  data.frame(otu=.) %>%
  left_join(silva_tax_df)

Linda_df16S_RNA=
  Linda_16S_SiteType_RNA$output$`SiteTypeMid Elevation` %>%
  mutate(otu=rownames(.)) %>%
  left_join(silva_tax_df) %>%
  filter(padj<p_cutoff_SS, (log2FoldChange>lfc_cutoff_SS | log2FoldChange<(-lfc_cutoff_SS))) %>%
  mutate(Direction=ifelse(log2FoldChange>0, "Positive","Negative")) %>%
  arrange(Direction,(log2FoldChange)) %>% 
  mutate(otu=factor(otu,levels=(.$otu)),
         TaxClass=ifelse(!is.na(Species),Species,ifelse(!is.na(Genus),Genus,ifelse(!is.na(Family),Family,ifelse(!is.na(Order),Order,ifelse(!is.na(Class),Class,ifelse(!is.na(Phylum),Phylum,Domain)))))))

ggplot(Linda_df16S_RNA,aes(x=otu, y=log2FoldChange,color=Direction)) +
  geom_point() +
  geom_errorbar(aes(ymax=log2FoldChange+lfcSE, ymin=log2FoldChange-lfcSE)) +
  geom_hline(yintercept = 0, color="grey", size=1) +
  scale_x_discrete(labels=Linda_df16S_RNA$TaxClass) +
  theme_minimal() +
  coord_flip()

####Generates dataframes of merged results to make into supplemental tables
format_linda=function(Linda_df, Treatment1, Treatment2){
  Out_df=
    Linda_df %>%
    as.data.frame %>%
    mutate(otu=rownames(.)) %>%
    select(otu, padj, log2FoldChange) %>%
    set_colnames(c("otu", paste(Treatment1, "_vs_",Treatment2,"_padj",sep=""), 
             paste(Treatment1, "_vs_",Treatment2,"_log2FC",sep="")))
  
  return(Out_df)
}

#########ITS DNA####################
HvM_ITS_DNA_supp=
  Linda_ITS_SiteType_DNA$output$`SiteTypeMid Elevation` %>% 
  format_linda("HighElevation","MidElevation")

HvR_ITS_DNA_supp=
  Linda_ITS_SiteType_DNA$output$SiteTypeRiparian %>% 
  format_linda("HighElevation","Riparian")

MvR_ITS_DNA_supp=
  Linda_ITS_SiteType_mid_June_DNA$output$SiteTypeRiparian %>% 
  format_linda("MidElevation","Riparian")

AvJ_ITS_DNA_supp=
  Linda_ITS_SiteType_DNA$output$TimeJune %>% 
  format_linda("August","June")

AvO_ITS_DNA_supp=
  Linda_ITS_SiteType_DNA$output$TimeOctober %>% 
  format_linda("August","October")

JvO_ITS_DNA_supp=
  Linda_ITS_SiteType_mid_June_DNA$output$TimeOctober %>% 
  format_linda("June","October")

LINDA_ITS_DNA_summary=
  tax_mat_df %>% 
  select(-TaxClass) %>%
  left_join(HvM_ITS_DNA_supp) %>%
  left_join(HvR_ITS_DNA_supp) %>% 
  left_join(MvR_ITS_DNA_supp) %>%
  left_join(AvJ_ITS_DNA_supp) %>%
  left_join(AvO_ITS_DNA_supp) %>%
  left_join(JvO_ITS_DNA_supp)

#########ITS RNA####################
HvM_ITS_RNA_supp=
  Linda_ITS_SiteType_RNA$output$`SiteTypeMid Elevation` %>% 
  format_linda("HighElevation","MidElevation")

HvR_ITS_RNA_supp=
  Linda_ITS_SiteType_RNA$output$SiteTypeRiparian %>% 
  format_linda("HighElevation","Riparian")

MvR_ITS_RNA_supp=
  Linda_ITS_SiteType_mid_June_RNA$output$SiteTypeRiparian %>% 
  format_linda("MidElevation","Riparian")

AvJ_ITS_RNA_supp=
  Linda_ITS_SiteType_RNA$output$TimeJune %>% 
  format_linda("August","June")

AvO_ITS_RNA_supp=
  Linda_ITS_SiteType_RNA$output$TimeOctober %>% 
  format_linda("August","October")

JvO_ITS_RNA_supp=
  Linda_ITS_SiteType_mid_June_RNA$output$TimeOctober %>% 
  format_linda("June","October")

LINDA_ITS_RNA_summary=
  tax_mat_df %>% 
  select(-TaxClass) %>%
  left_join(HvM_ITS_RNA_supp) %>%
  left_join(HvR_ITS_RNA_supp) %>% 
  left_join(MvR_ITS_RNA_supp) %>%
  left_join(AvJ_ITS_RNA_supp) %>%
  left_join(AvO_ITS_RNA_supp) %>%
  left_join(JvO_ITS_RNA_supp)

#########16S DNA####################
HvM_16S_DNA_supp=
  Linda_16S_SiteType_DNA$output$`SiteTypeMid Elevation` %>% 
  format_linda("HighElevation","MidElevation")

HvR_16S_DNA_supp=
  Linda_16S_SiteType_DNA$output$SiteTypeRiparian %>% 
  format_linda("HighElevation","Riparian")

MvR_16S_DNA_supp=
  Linda_16S_SiteType_mid_June_DNA$output$SiteTypeRiparian %>% 
  format_linda("MidElevation","Riparian")

AvJ_16S_DNA_supp=
  Linda_16S_SiteType_DNA$output$TimeJune %>% 
  format_linda("August","June")

AvO_16S_DNA_supp=
  Linda_16S_SiteType_DNA$output$TimeOctober %>% 
  format_linda("August","October")

JvO_16S_DNA_supp=
  Linda_16S_SiteType_mid_June_DNA$output$TimeOctober %>% 
  format_linda("June","October")

LINDA_16S_DNA_summary=
  silva_tax_df %>% 
  left_join(HvM_16S_DNA_supp) %>%
  left_join(HvR_16S_DNA_supp) %>% 
  left_join(MvR_16S_DNA_supp) %>%
  left_join(AvJ_16S_DNA_supp) %>%
  left_join(AvO_16S_DNA_supp) %>%
  left_join(JvO_16S_DNA_supp)

#########16S RNA####################
HvM_16S_RNA_supp=
  Linda_16S_SiteType_RNA$output$`SiteTypeMid Elevation` %>% 
  format_linda("HighElevation","MidElevation")

HvR_16S_RNA_supp=
  Linda_16S_SiteType_RNA$output$SiteTypeRiparian %>% 
  format_linda("HighElevation","Riparian")

MvR_16S_RNA_supp=
  Linda_16S_SiteType_mid_June_RNA$output$SiteTypeRiparian %>% 
  format_linda("MidElevation","Riparian")

AvJ_16S_RNA_supp=
  Linda_16S_SiteType_RNA$output$TimeJune %>% 
  format_linda("August","June")

AvO_16S_RNA_supp=
  Linda_16S_SiteType_RNA$output$TimeOctober %>% 
  format_linda("August","October")

JvO_16S_RNA_supp=
  Linda_16S_SiteType_mid_June_RNA$output$TimeOctober %>% 
  format_linda("June","October")

LINDA_16S_RNA_summary=
  silva_tax_df %>% 
  left_join(HvM_16S_RNA_supp) %>%
  left_join(HvR_16S_RNA_supp) %>% 
  left_join(MvR_16S_RNA_supp) %>%
  left_join(AvJ_16S_RNA_supp) %>%
  left_join(AvO_16S_RNA_supp) %>%
  left_join(JvO_16S_RNA_supp)

write.csv(LINDA_ITS_DNA_summary, "LINDA_ITS_DNA_summary.csv")
write.csv(LINDA_ITS_RNA_summary, "LINDA_ITS_RNA_summary.csv")
write.csv(LINDA_16S_DNA_summary, "LINDA_16S_DNA_summary.csv")
write.csv(LINDA_16S_RNA_summary, "LINDA_16S_RNA_summary.csv")

write.csv(DNAvsRNA_ITS_sig, "DNAvsRNA_ITS_sig.csv")
write.csv(DNAvsRNA_16S_sig, "DNAvsRNA_16S_sig.csv")

######################Venn diagrams############################
#venn diagrams are created to show # of OTUs detected as differentially abundant
#by site type and season in DNA and RNA datasets (and those detected in both datasets)
#done for both DNA and RNA

#creates list of site type-variable fungal OTUs
sitetype_ITS_venn=
  list(DNA=SiteType_ITS_variable_DNA$otu,
       RNA=SiteType_ITS_variable_RNA$otu)

#creates list of seasonally variable fungal OTUs
time_ITS_venn=
  list(DNA=Time_ITS_variable_DNA$otu,
       RNA=Time_ITS_variable_RNA$otu)

#creates list of site type-variable bacterial ASVs
sitetype_16S_venn=
  list(DNA=SiteType_16S_variable_DNA$otu,
       RNA=SiteType_16S_variable_RNA$otu)

#creates list of seasonally-variable bacterial ASVs
time_16S_venn=
  list(DNA=Time_16S_variable_DNA$otu,
       RNA=Time_16S_variable_RNA$otu)

#creates Venn Diagram of site type-variable fungal OTUs
sitetype_ITS_venn_plot=
  ggvenn(sitetype_ITS_venn, fill_color=c("#265DAB", "#8d4653"), auto_scale=TRUE, fill_alpha=0.55)

#creates Venn diagram of seasonally-variable fungal OTUs
time_ITS_venn_plot=
  ggvenn(time_ITS_venn, fill_color=c("#265DAB", "#8d4653"), auto_scale=TRUE, fill_alpha=0.55)

#creates Venn diagram of site type-variable bacterial ASVs
sitetype_16S_venn_plot=
  ggvenn(sitetype_16S_venn, fill_color=c("#265DAB", "#8d4653"), auto_scale=TRUE, fill_alpha=0.55)

#creates Venn diagram of seasonally-variable bacterial ASVs
time_16S_venn_plot=
  ggvenn(time_16S_venn, fill_color=c("#265DAB", "#8d4653"), auto_scale=TRUE, fill_alpha=0.55)

#creates multipanel figure with all Venn diagrams
Linda_venn=plot_grid(sitetype_ITS_venn_plot, sitetype_16S_venn_plot,
                     time_ITS_venn_plot, time_16S_venn_plot)

#saves multipanel Venn diagram as a pdf that will be incorporated into a figure
pdf(file="Linda_venn.pdf", 
    width = 3.5, height = 5.11811)
Linda_venn
dev.off()

