library(tidypopgen)
library(admixtools)
library(ggplot2)

dir("~/embo_popgen_2025_Betty/Betty/DAY5")
setwd("~/embo_popgen_2025_Betty/Betty/DAY5")


modern_gt <- tidypopgen::gen_tibble("./data/modern_samples.bed",
                                    valid_alleles = c("A","T","C","G"),
                                    missing_alleles = c("X"))
#To load again this file, use gt_load() to avoid doing the command again.

modern_gt

modern_gt %>% group_by(population) %>% tally()

#%>% is the same as the pipe | in bash

loci_report <- modern_gt %>%
  group_by(population) %>%
  qc_report_loci()

autoplot(loci_report, type = "missing")
  #The proportion of missing data is around 0.01

modern_gt <- modern_gt %>% 
  select_loci_if(loci_missingness(genotypes)<0.04)
#Remove loci with a genotype less than 0.04

ancient_gt <- tidypopgen::gen_tibble("./data/ancient_samples.vcf",
                                     valid_alleles = c("A","T","C","G","X"), 
                                     quiet = TRUE)

ancient_gt <- gt_pseudohaploid(ancient_gt)

ancient_gt$id

ancient_gt$id[ancient_gt$id == "GB20"] <- "Mota"

ancient_gt$population <- ancient_gt$id

merged_dry <- rbind_dry_run(modern_gt, ancient_gt, 
                            flip_strand = TRUE)
  #Chck if the SNP are the same in the modern and ancient dataset to remove the ones that aren't found in both

merged_gt <- rbind(modern_gt, ancient_gt, 
                   flip_strand = TRUE,
                   backingfile = "./data/merged_samples")
  #Merge the panels

merged_gt <- merged_gt %>% group_by(population)
  #Group data by populations

######################################################################################################

f2_dir <- "./data/f2_tidypopgen"


# compute f2
f2_tidypopgen <- gt_extract_f2(merged_gt, 
                               outdir = "./data/f2_tidypopgen", 
                               overwrite = TRUE)
    #Create the directory for the frequencies of the sample

f2_blocks = f2_from_precomp("./data/f2_tidypopgen")

######################################################################################################

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
                      "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

lbk_f3out <- f3(data = f2_blocks, 
                pop1 = "Mbuti", 
                pop2 = "LBK",
                pop3 = lbk_modern_panel)
lbk_f3out
  #Pop 1 = Outgroup
  #Pop 2 = Target
  #Pop 3 = Test 
    #If p wasn't significant, the ree topology would have been wrongly chosen
  
#Sardinians are the closest (also bc they are from an island that didn't experience as much migration so they look more like Neolithic populations than ohter Europeans populations)

lbk_f3out %>% arrange(desc(est))


ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


lbk_f3out$pop3<-factor(lbk_f3out$pop3, levels = lbk_f3out$pop3[order(lbk_f3out$est)])

ggplot(lbk_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with LBK", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Very close to European pop, this unknown sample is European.

######################################################################################################
MA1_f3out <- f3(data = f2_blocks, 
                pop1 = "Mbuti", 
                pop2 = "MA1",
                pop3 = lbk_modern_panel)
MA1_f3out
MA1_f3out %>% arrange(desc(est))

MA1_f3out$pop3<-factor(MA1_f3out$pop3, levels = MA1_f3out$pop3[order(MA1_f3out$est)])

ggplot(MA1_f3out, aes(pop3, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with MA1", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#MA1 is Central Asia : The closest populations are from America but they descended from Siberian people
    #They are also near to a set of European population : Contact with Caucasian populations

######################################################################################################


aa_f3admix <- f3(data = f2_blocks,
                 pop1 = "AA",
                 pop2 = "Yoruba",
                 pop3 = "French")
aa_f3admix
#F3 negative and different to 0 : Evidence for admixture

######################################################################################################

eurasian_sources <- c("French","Spanish","Sardinian","LBK")

somali_f3admix <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = eurasian_sources, 
                     pop3 = "Mota")
somali_f3admix
  #Very strong negative Z values : Admixture is present

somali_f3admix$pop2<-factor(somali_f3admix$pop2, levels = somali_f3admix$pop2[order(somali_f3admix$est)])

ggplot(somali_f3admix, aes(pop2, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with Eurasian sources", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Mota (from Ethiopia) : No ancestry from European people
  #They share more drift with the outgroup population 

######################################################################################################

somali_f3admix_2 <- f3(data = f2_blocks,
                     pop1 = "Somali",
                     pop2 = eurasian_sources, 
                     pop3 = "Dinka")
somali_f3admix_2
#Very strong negative Z values : Admixture is present

somali_f3admix_2$pop2<-factor(somali_f3admix_2$pop2, levels = somali_f3admix_2$pop2[order(somali_f3admix_2$est)])

ggplot(somali_f3admix_2, aes(pop2, est)) +
  geom_point() +
  geom_errorbar(aes(ymin = est - 2 * se, ymax = est + 2 * se)) +
  labs(y = "Shared drift with Eurasian sources", x = "populations") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#For Dinka (Sudan) : Here, the French and Spanish seems to "share" a bit more drift with Dinka
  #Back to Africa ?

######################################################################################################

neand_french_f4 <- f4(data = f2_blocks,
                      pop1 = "pan_troglodytes",
                      pop2 = "AltaiNea",
                      pop3 = "Mbuti",
                      pop4 = "French")
neand_french_f4

#Here, there is gene flow from AltaiNea into the French pop (Z positive and significant)

######################################################################################################

neand_han_f4 <- f4(data = f2_blocks,
                      pop1 = "pan_troglodytes",
                      pop2 = "AltaiNea",
                      pop3 = "Mbuti",
                      pop4 = "Han")
neand_han_f4

#This might be due to the fact that Nea and Den are closer together
  #Han populations have some Nea but this signal is mostly from Den

######################################################################################################

neand_AA_f4 <- f4(data = f2_blocks,
                   pop1 = "pan_troglodytes",
                   pop2 = "French",
                   pop3 = "Mbuti",
                   pop4 = "AA")
neand_AA_f4

#Gene flow from European into AA (more than African pop) : This results is the same for different European pop

######################################################################################################


pops <- c("Han", "pan_troglodytes", "AA","Yoruba","French")
qpf4ratio(data = f2_blocks, 
          pops = pops)

lbk_modern_panel <- c("Basque", "Bedouin2", "Druze", "Cypriot", "Tuscan",
                      "Sardinian", "French", "Spanish", "Onge", "Han", "Mayan", "Mixe", "Surui")

######################################################################################################

modern_panel_gt <- merged_gt %>% filter (population %in% lbk_modern_panel)
# remove monomorphic sites (ungrouping the tibble, as we want to use global frequencies
# rather than population frequencies)
modern_panel_gt <- modern_panel_gt %>% ungroup() %>% select_loci_if(loci_maf(genotypes)>0)

modern_panel_gt <- gt_update_backingfile(modern_panel_gt)
# reset the ploidy for this tibble, as it is now all diploids
attr(modern_panel_gt$genotypes, "ploidy") <- 2
modern_panel_gt <- modern_panel_gt %>% gt_impute_simple(method="mode")

modern_panel_gt <- modern_panel_gt %>% select_loci_if(loci_ld_clump(genotypes))

modern_pca <- modern_panel_gt %>% gt_pca_randomSVD()

autoplot(modern_pca, type = "scores") 

library(ggplot2)
autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")


lbk_gt <- merged_gt %>% filter(id == "LBK")
lbk_pca_scores <- predict(modern_pca, new_data = lbk_gt, project_method = "least_square")
lbk_pca_scores

autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population")+
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black")


autoplot(modern_pca, type = "scores") +
  aes(color = modern_panel_gt$population) +
  labs(color = "population") +
  geom_point(data=lbk_pca_scores, mapping=aes(x=.data$.PC1, y=.data$.PC2), col = "black") +
  lims(x=c(30, 70), y = c(-10, 15))

