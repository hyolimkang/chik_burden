setwd("D:/OneDrive - London School of Hygiene and Tropical Medicine/CHIK/1.Aim1/codes/CHIK")
setwd("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/CHIK/1.Aim1/codes/CHIK")
ssa_mcmc$mcmc9 <- mcmc9
ssa_mcmc$mcmc11 <- mcmc11
ssa_mcmc$mcmc12 <- mcmc12
ssa_mcmc$mcmc13 <- mcmc13
ssa_mcmc$mcmc14 <- mcmc14
ssa_mcmc$mcmc15 <- mcmc15
lac_mcmc$mcmc29 <- mcmc29
lac_mcmc$mcmc30 <- mcmc30
lac_mcmc$mcmc31 <- mcmc31
eap_mcmc$mcmc35 <- mcmc35
eap_mcmc$mcmc36 <- mcmc36
eap_mcmc$mcmc37 <- mcmc37
eap_mcmc$mcmc38 <- mcmc38
eap_mcmc$mcmc39 <- mcmc39
eap_mcmc$mcmc40 <- mcmc40
eap_mcmc$mcmc41 <- mcmc41
eap_mcmc$mcmc42 <- mcmc42
eap_mcmc$mcmc43 <- mcmc43
eap_mcmc$mcmc44 <- mcmc44
eap_mcmc$mcmc45 <- mcmc45
eap_mcmc$mcmc46 <- mcmc46
eap_mcmc$mcmc47 <- mcmc47
eap_mcmc$mcmc50 <- mcmc50
eap_mcmc$mcmc51 <- mcmc51
eap_mcmc$mcmc52 <- mcmc52
eap_mcmc$mcmc53 <- mcmc53

#manually add study_no
ssa_mcmc$mcmc1$ID <- 1
ssa_mcmc$mcmc2$ID <- 29
ssa_mcmc$mcmc3$ID <- 55
ssa_mcmc$mcmc4$ID <- 23
ssa_mcmc$mcmc5$ID <- 11
ssa_mcmc$mcmc6$ID <- 51
ssa_mcmc$mcmc7$ID <- 56
ssa_mcmc$mcmc8$ID <- 24
ssa_mcmc$mcmc9$ID <- 2
ssa_mcmc$mcmc10$ID <- 6
ssa_mcmc$mcmc11$ID <- 32
ssa_mcmc$mcmc12$ID <- 33
ssa_mcmc$mcmc13$ID <- 34
ssa_mcmc$mcmc14$ID <- 35
ssa_mcmc$mcmc15$ID <- 36
ssa_mcmc$mcmc16$ID <- 15
save(ssa_mcmc, file = "ssa_mcmc.RData")

lac_mcmc$mcmc17$ID <- 4
lac_mcmc$mcmc18$ID <- 7
lac_mcmc$mcmc19$ID <- 12
lac_mcmc$mcmc20$ID <- 14
lac_mcmc$mcmc21$ID <- 16
lac_mcmc$mcmc22$ID <- 20
lac_mcmc$mcmc23$ID <- 26
lac_mcmc$mcmc24$ID <- 27
lac_mcmc$mcmc25$ID <- 31
lac_mcmc$mcmc26$ID <- 43
lac_mcmc$mcmc27$ID <- 44
lac_mcmc$mcmc28$ID <- 52
lac_mcmc$mcmc29$ID <- 53
lac_mcmc$mcmc30$ID <- 54
lac_mcmc$mcmc31$ID <- 17
lac_mcmc$mcmc32$ID <- 25
lac_mcmc$mcmc33$ID <- 42

eap_mcmc$mcmc34$ID <- 28
eap_mcmc$mcmc35$ID <- 74
eap_mcmc$mcmc36$ID <- 75
eap_mcmc$mcmc37$ID <- 76
eap_mcmc$mcmc38$ID <- 45
eap_mcmc$mcmc39$ID <- 46
eap_mcmc$mcmc40$ID <- 47
eap_mcmc$mcmc41$ID <- 48
eap_mcmc$mcmc42$ID <- 8
eap_mcmc$mcmc43$ID <- 9
eap_mcmc$mcmc44$ID <- 10
eap_mcmc$mcmc45$ID <- 30
eap_mcmc$mcmc46$ID <- 21
eap_mcmc$mcmc47$ID <- 22
eap_mcmc$mcmc48$ID <- 49
eap_mcmc$mcmc49$ID <- 5
eap_mcmc$mcmc50$ID <- 38
eap_mcmc$mcmc51$ID <- 39
eap_mcmc$mcmc52$ID <- 40
eap_mcmc$mcmc53$ID <- 41
save(eap_mcmc, file = "eap_mcmc.RData")

sa_mcmc$mcmc54$ID <- 63
sa_mcmc$mcmc55$ID <- 18
sa_mcmc$mcmc56$ID <- 50
sa_mcmc$mcmc57$ID <- 57
sa_mcmc$mcmc58$ID <- 58
sa_mcmc$mcmc59$ID <- 59
sa_mcmc$mcmc60$ID <- 61
sa_mcmc$mcmc61$ID <- 62
sa_mcmc$mcmc62$ID <- 64
sa_mcmc$mcmc63$ID <- 65
sa_mcmc$mcmc64$ID <- 66
sa_mcmc$mcmc65$ID <- 67
sa_mcmc$mcmc66$ID <- 68
sa_mcmc$mcmc67$ID <- 69
sa_mcmc$mcmc68$ID <- 71
sa_mcmc$mcmc69$ID <- 72
sa_mcmc$mcmc70$ID <- 73
sa_mcmc$mcmc71$ID <- 70
sa_mcmc$mcmc72$ID <- 60

ecame_mcmc$mcmc73$ID <- 3
ecame_mcmc$mcmc74$ID <- 13
ecame_mcmc$mcmc75$ID <- 19
ecame_mcmc$mcmc76$ID <- 37


load("ssa_mcmc.RData")
load("lac_mcmc.RData")
load("eap_mcmc.RData")
load("sa_mcmc.RData")
load("ecame_mcmc.RData")
load("mcmc_tot.RData")

#-------------------------------------------------------------------------------
setwd("D:/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping") # window
num_samples <- 100 

mcmc_tot <- c(ssa_mcmc, lac_mcmc, eap_mcmc, sa_mcmc, ecame_mcmc)

mcmc_tot <- lapply(mcmc_tot, function(df) {
  # Adding a small value to probabilities that are 0 or 1
  small_value <- 0.000001
  df$adjusted_FOI <- ifelse(df[[2]] == 0, small_value, df[[2]])
  
  # Applying log transformation
  df$logFOI <- log(df$adjusted_FOI)
  df$logitFOI <- log(df$adjusted_FOI / (1 - df$adjusted_FOI))
  return(df)
})

new_colnames <- c("study_no", "FOI", "lat", "long", "ID", "adjusted_FOI", "logFOI", "logitFOI")

# Use lapply to iterate over the list and set the new names
mcmc_tot <- lapply(mcmc_tot, function(df) {
  colnames(df) <- new_colnames
  return(df)
})

mcmc_df <- do.call(rbind, mcmc_tot)
nan_rows <- is.nan(mcmc_df$logitFOI)
nan_df <- mcmc_df[nan_rows, ]

save(mcmc_tot, file = "mcmc_tot.RData")


