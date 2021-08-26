library("tidyverse")

# calculate FPKM for our six species
Caus <- data.frame(samples= c("SRR2917198", "SRR2917152", "SRR2916699", 
                              "SRR2916026", "SRR2916025", "SRR2915407"),
                   FPKM=c(312325, 322874,326566, 301927, 303706, 331231),
                   phenotype=c("NR", "NR","NR", "R","R","R"),
                   stringsAsFactors = FALSE)

# total FPKM 1898629
sum(Caus$FPKM)

# mean FPKM per sample 316438.2
sum(Caus$FPKM)/nrow(Caus)

# mean FPKM per R sample 316438.2
sample_per_group <- 3
Caus %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE)/sample_per_group)
# NR                    106863.
# R                     104096 

###############################################################################


# calculate FPKM for our six species
Ccal <- data.frame(samples= c("SRR13281988","SRR13281982","SRR13281995", 
                              "SRR13282007","SRR13281993","SRR13282006"),
                   FPKM=c(500812, 449078, 441697, 510660, 517649, 539811),
                   phenotype=c("R","R","R","NR", "NR","NR"),
                   stringsAsFactors = FALSE)

# total FPKM 2959707
sum(Ccal$FPKM)

# mean FPKM per sample 493284.5
sum(Ccal$FPKM)/nrow(Ccal)

# mean FPKM per R sample 
sample_per_group <- 3
Ccal %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE)/sample_per_group)
# NR                    174236.
# R                     154621.

###############################################################################

# calculate FPKM for our six species
Lf <- data.frame(samples= c("SocF10","SocF11","SocF12","SocF14","SocF4","SocF5",
                              "SocF6","SocF7","SocF8","SocF9","SocR10","SocR1","SocR2",
                              "SocR3","SocR4","SocR5","SocR6","SocR7","SocR8"),
                   FPKM=c(436905,440441,449393,492183,445877,444259,440756,460104,
                          474888,456590,456233,438030,433900,423696,442500,429253,
                          425667,438659,441837),
                   phenotype=c("NR", "NR","NR","NR", "NR",
                               "NR","NR", "NR","NR", "NR",
                               "R","R","R",
                               "R","R","R",
                               "R","R","R"),
                   stringsAsFactors = FALSE)

# total FPKM 8471171
sum(Lf$FPKM)

# mean FPKM per sample 445851.1
sum(Lf$FPKM)/nrow(Lf)

# mean FPKM per R sample 
NR_sample_num <- 10
R_sample_num <- 9
Lf %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE))
# NR                    454140..
# R                     436642.
454140/NR_sample_num # 45414
436642/R_sample_num # 48515.78
###############################################################################




# calculate FPKM for our six species
Mgen <- data.frame(samples= c("SRR3948522","SRR3948580","SRR3948532",
                              "SRR3948545","SRR3948559","SRR3948571",
                              "SRR3948576","SRR3948549","SRR3948582",
                              "SRR3948523","SRR3948530","SRR3948535",
                              "SRR3948538","SRR3948546","SRR3948551",
                              "SRR3948569"),
                 FPKM=c(387728,384013,370058,389816,363871,351825,
                        370639,390251,424446,354827,402008,353445,
                        390396,346395,384567,398066),
                 phenotype=c("R","R","R",
                             "R","R","R",
                             "R",
                             "NR", "NR","NR",
                             "NR","NR", "NR",
                             "NR", "NR","NR"
                             ),
                 stringsAsFactors = FALSE)

# total FPKM 6062351
sum(Mgen$FPKM)

# mean FPKM per sample 378896.9
sum(Mgen$FPKM)/nrow(Mgen)

# mean FPKM per R sample 
NR_sample_num <- 9
R_sample_num <- 7
Mgen %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE))
# NR                    382711.
# R                     373993.
382711/NR_sample_num # 42523.44
373993/R_sample_num # 53427.57


###############################################################################

# calculate FPKM for our six species
Pcan <- data.frame(samples= c("SRR1519108","SRR1519109","SRR1519110",
                              "SRR1519111","SRR1519112","SRR1519113",
                              "SRR1519114","SRR1519115","SRR1519116",
                              "SRR1519117"),
                   FPKM=c(439317,449922,452254,442345,454696,
                          446554,445604,428060,470317,444633),
                   phenotype=c("R","R",
                               "R","R",
                               "NR", "NR","NR",
                               "NR","NR", "NR"),
                   stringsAsFactors = FALSE)

# total FPKM 4473702
sum(Pcan$FPKM)

# mean FPKM per sample 447370.2
sum(Pcan$FPKM)/nrow(Pcan)

# mean FPKM per R sample 
NR_sample_num <- 6
R_sample_num <- 4
Pcan %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE))
# NR                    448311.
# R                     445960.
448311/NR_sample_num # 74718.5
445960/R_sample_num # 111490

###############################################################################

# calculate FPKM for our six species
Pdom <- data.frame(samples= c("F03BL","F14GR","F38BL","F42OR","F45OR","F47RD",
                               "F48SL","F51RD","F65ORWH","F67YLBL","F69OR",
                               "F71BL","W01RDYL","W08RDOR","W09RDOR","W13YLSL",
                               "W28ORYL","W40GROR","W46ORYL","W49RDYL",
                               "W54ORYL","W56YLBL","W58BLWH","W70RDYL"),
                   FPKM=c(498481,446329,452346,444375,421681,424785,446986,
                          440807,463581,443491,448383,439222,449501,446851,
                          448616,470162,456376,438354,459054,478506,465422,
                          477532,507484,470962),
                   phenotype=c("R","R","R","R",
                               "R","R","R","R",
                               "R","R","R","R",
                               "NR", "NR","NR","NR",
                               "NR","NR", "NR","NR",
                               "NR","NR", "NR","NR"),
                   stringsAsFactors = FALSE)

# total FPKM 10939287
sum(Pdom$FPKM)

# mean FPKM per sample 455803.6
sum(Pdom$FPKM)/nrow(Pdom)

# mean FPKM per R sample 
NR_sample_num <- 12
R_sample_num <- 12
Pdom %>% 
  group_by(phenotype) %>%
  summarize(mean_per_sampleFPKM = mean(FPKM, na.rm = TRUE))
# NR                    464068.
# R                     447539.
464068/NR_sample_num # 38672.33
447539/R_sample_num # 37294.92


