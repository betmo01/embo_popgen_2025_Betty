library(admixtools)
library(tidypopgen)
install.packages("plotly")
library(plotly)

dir("~/embo_popgen_2025_Betty/Betty/DAY5")
setwd("~/embo_popgen_2025_Betty/Betty/DAY5")

f2_blocks = f2_from_precomp("./data/f2_tidypopgen", verbose = FALSE)

######################################################################################################
#Lower bound of the number of waves that have gone from the right pop to the left : minimum bound

neand_euras_wave <- qpwave(data = f2_blocks,
                           left = c("French","Spanish","Tuscan"),
                           right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave
#The rank 0 matrix model can't be rejected since p-value > 0.05 : Only one wave (but also the case for rank 1)
#Nested comparison (last column, p_nested) : if less than 0.05, the rank 1 (above) doesn't lead to a significant drop of explanatory power compared to rank 0
# No increase of explanatory power = we will choose the lowest rank (parsimony) for the minimum nb of waves

######################################################################################################

#Lower bound of the number of waves that have gone from the right pop to the left

neand_euras_wave_2 <- qpwave(data = f2_blocks,
                             left = c("French","Spanish","Tuscan", "Han", "Onge"),
                             right = c("AltaiNea","Mota", "Yoruba", "Denisova", "Mbuti")
)
neand_euras_wave_2
#The rank 0 matrix can be rejected since p-value < 0.05 but not the case for rank 1 + p_nested is significant so the minimum nb of wave is 2 (corresponding to rank 1)
#This second wave is most likely a Denisovian admixture

######################################################################################################
french_adm <- qpadm(data = f2_blocks,
                    left = c("Loschbour", "LBK", "Yamnaya"),
                    right = c("Mbuti", "Mota", "Dinka", "Yoruba", "Han"),
                    target= "French")

french_adm$popdrop

#p-value < 0.05 : rejected
#feasible = FALSE : rejected (because of negative component of a population)

#Here, two models could fit the data from these results but one, the 001, is very close to being rejected by the p-value (and is actually not true)
#Always choose the one with the less hypotheses in general (parsimony)

######################################################################################################
#R : root
#eAfr : East African last common ancestors bt Dinka and non-Africans
#outAfrica : MRCA bt Asians and Europeans

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

base_edges %>% edges_to_igraph() %>% plot_graph()

base_igraph %>% plotly_graph()


base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)
base_qpgraph$f3

#The graph is compatible to the data since the null hypothesis (H0 : It is compatible) is above 0.05

base_pqgaph$f3 %>% filter(abs(z)>2)
base_qpgraph$edges %>% plotly_graph()

######################################################################################################

base_swapped_edges <- matrix(
  c("R",	"Dinka",
    "R", "eAfr",
    "eAfr",	"Mbuti",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_swapped_igraph <- base_swapped_edges %>% edges_to_igraph()

is_valid(base_swapped_igraph)

base_swapped_edges %>% edges_to_igraph() %>% plot_graph()

base_swapped_igraph %>% plotly_graph()


base_swapped_qpgraph <- qpgraph(data = f2_blocks, graph = base_swapped_igraph)
base_swapped_qpgraph$f3

#The graph is compatible to the data since the null hypothesis (H0 : It is compatible) is above 0.05
#This is because the model will the Mbuti as just a very "drifted" pop from Dinka (which isn't the case)
#Here the root is "cosmetic", Dinka is still related to Han and Loschbour but since Mbuti isn't linked to another pop on this tree, it can be just seen as a result of drift
#So we need to have correct assumptions about our root before because the model won't be a good tool to find it


#African line that go all the way down. If we have Mbuti as a "unique" line, it will just be seen as a big proportion of drift by qpGraph

######################################################################################################

base_edges <- matrix(
  c("R",	"Mbuti",
    "R", "eAfr",
    "eAfr",	"Dinka",
    "eAfr",	"outAfrica",
    "outAfrica",	"Han",
    "outAfrica",	"Loschbour"),
  ncol=2,
  byrow = TRUE,
  dimnames=list(NULL, c("from","to")))

base_edges

base_igraph <- base_edges %>% edges_to_igraph()

is_valid(base_igraph)

base_igraph %>% plot_graph()

base_igraph %>% plotly_graph()

base_qpgraph <- qpgraph(data = f2_blocks, graph = base_igraph)

base_qpgraph$f3

base_qpgraph$f3 %>% filter(abs(z)>2)

base_qpgraph$edges %>% plot_graph()

fits = qpgraph_resample_multi(f2_blocks, 
                              graphlist = list(base_qpgraph[[1]], base_swapped_qpgraph[[1]]), 
                              nboot = 100)
compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

base_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

######################################################################################################

yamnaya_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "Dinka",
    "eAfr", "outAfrica",
    "outAfrica", "Han",
    "outAfrica", "wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "Loschbour"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))

yamnaya_igraph <- yamnaya_edges %>% edges_to_igraph()
yamnaya_igraph %>% plot_graph()


######################################################################################################

lbk_extra_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))

lbk_extra_igraph <- lbk_extra_edges %>% edges_to_igraph()
lbk_extra_igraph %>% plot_graph()

is_valid(lbk_extra_igraph)

lbk_extra_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

lbk_extra_qpgraph <- qpgraph(data = f2_blocks, graph = lbk_extra_igraph)
lbk_extra_qpgraph$edges %>% plot_graph()
#The number correspond to the drift

######################################################################################################

Sard_add_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "LBK", "Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))

Sard_add_igraph <- Sard_add_edges %>% edges_to_igraph()
Sard_add_igraph %>% plot_graph()

is_valid(Sard_add_igraph)

Sard_add_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

Sard_add_qpgraph <- qpgraph(data = f2_blocks, graph = Sard_add_igraph)
Sard_add_qpgraph$edges %>% plot_graph()


Yamn_add_edges <- matrix(
  c(
    "R",	"Mbuti",
    "R",	"eAfr",
    "eAfr", "pBasalEurasian",
    "eAfr", "Dinka",
    "pBasalEurasian", "BasalEurasian",
    "pBasalEurasian","outAfrica",
    "outAfrica", "Han",
    "outAfrica","wEurasian",
    "wEurasian", "Yamnaya",
    "wEurasian", "pLoschbour",
    "pLoschbour", "Loschbour",
    "pLoschbour","WHG",
    "BasalEurasian", "pLBK",
    "WHG", "pLBK",
    "pLBK","LBK",
    "LBK", "Sardinian",
    "Yamnaya", "Sardinian"),
  ncol = 2,
  byrow = TRUE,
  dimnames = list(NULL, c("from", "to")))

Yamn_add_igraph <- Yamn_add_edges %>% edges_to_igraph()
Yamn_add_igraph %>% plot_graph()

is_valid(Yamn_add_igraph)

Yamn_add_igraph %>% plot_graph(highlight_unidentifiable = TRUE)

Yamn_add_qpgraph <- qpgraph(data = f2_blocks, graph = Yamn_add_igraph)
Yamn_add_qpgraph$edges %>% plot_graph()

fits = qpgraph_resample_multi(f2_blocks,
                              graphlist = list(Sard_add_qpgraph[[1]], Yamn_add_qpgraph [[1]]),
                              nboot = 100)

compare_fits(fits[[1]]$score_test, fits[[2]]$score_test)

#No difference if we add Yamnaya, so the firs tmodel is best by parsimony (=we add less variables)
#Sardinians can be modeled with Yamnaya