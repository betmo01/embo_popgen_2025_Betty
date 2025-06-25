obs_data <- read.csv("~/embo_popgen_2025_Betty/Matthias_Steinrücken/day3/practicals/mosquito-observed.csv")
sim_data <- read.csv("~/mosquito-task2.csv")

library(abc)

abc (target = obs_data , #obs sum stat
     param = , #list of simulated param
     sumstat = , #simulated stat
     tol = , #epsilone
     method = ,#rejection)