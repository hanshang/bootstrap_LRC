SST_data_Nino_1_2 = matrix(read.table("SST.txt", header = TRUE)[,3], 12, 71, byrow = FALSE)
SST_data_Nino_3   = matrix(read.table("SST.txt", header = TRUE)[,5], 12, 71, byrow = FALSE)
SST_data_Nino_4   = matrix(read.table("SST.txt", header = TRUE)[,7], 12, 71, byrow = FALSE)
SST_data_Nino_3_4 = matrix(read.table("SST.txt", header = TRUE)[,9], 12, 71, byrow = FALSE)

colnames(SST_data_Nino_1_2) = colnames(SST_data_Nino_3) =
colnames(SST_data_Nino_4) = colnames(SST_data_Nino_3_4) = 1950:2020

rownames(SST_data_Nino_1_2) = rownames(SST_data_Nino_3) =
rownames(SST_data_Nino_4) = rownames(SST_data_Nino_3_4) = 1:12

# sample long-run covariance

SST_lrc_sample_Nino_1_2 = long_run_covariance_estimation(dat = SST_data_Nino_1_2)
SST_lrc_sample_Nino_3   = long_run_covariance_estimation(dat = SST_data_Nino_3)
SST_lrc_sample_Nino_4   = long_run_covariance_estimation(dat = SST_data_Nino_4)
SST_lrc_sample_Nino_3_4 = long_run_covariance_estimation(dat = SST_data_Nino_3_4)

persp(1:12, 1:12, SST_lrc_sample_Nino_1_2, xlab = "\n Month", ylab = "\n Month", ticktype = "detailed",
      zlab = "\n Sample long-run covariance", theta = 45, phi = 45)

savepdf("Fig_3a", width = 12, height = 10, toplines = 0.8)
plot(head(ts(as.numeric(SST_data_Nino_1_2), start = 1950, end = 2021, frequency = 12), -1), xlab = "Year",
     ylab = "Monthly sea surface temperature in Celcius")
dev.off()

savepdf("Fig_3b", width = 12, height = 10, toplines = 0.8)
plot(fts(1:12, SST_data_Nino_1_2), xlab = "Month", ylab = "")
dev.off()

T_stationary(SST_data_Nino_1_2)


##################
# sieve bootstrap
# B = 400
##################

SST_sieve_Nino_1_2 = sieve_bootstrap(fun_dat = t(SST_data_Nino_1_2), grid = 1:12,
                                     ncomp_porder_selection = "eigenratio_AICC",
                                     B = 400, burn_in = 50)

SST_sieve_Nino_3 = sieve_bootstrap(fun_dat = t(SST_data_Nino_3), grid = 1:12,
                                   ncomp_porder_selection = "eigenratio_AICC",
                                   B = 400, burn_in = 50)

SST_sieve_Nino_4 = sieve_bootstrap(fun_dat = t(SST_data_Nino_4), grid = 1:12,
                                   ncomp_porder_selection = "eigenratio_AICC",
                                   B = 400, burn_in = 50)

SST_sieve_Nino_3_4 = sieve_bootstrap(fun_dat = t(SST_data_Nino_3_4), grid = 1:12,
                                     ncomp_porder_selection = "eigenratio_AICC",
                                     B = 400, burn_in = 50)

SST_sieve_lrc_boot_Nino_1_2 = SST_sieve_lrc_boot_Nino_3 =
SST_sieve_lrc_boot_Nino_4 = SST_sieve_lrc_boot_Nino_3_4 = array(NA, dim = c(12, 12, 400))
for(iwk in 1:400)
{
    SST_sieve_lrc_boot_Nino_1_2[,,iwk] = long_run_covariance_estimation(dat = SST_sieve_Nino_1_2$X_boot[,,iwk])
    SST_sieve_lrc_boot_Nino_3[,,iwk]   = long_run_covariance_estimation(dat = SST_sieve_Nino_3$X_boot[,,iwk])
    SST_sieve_lrc_boot_Nino_4[,,iwk]   = long_run_covariance_estimation(dat = SST_sieve_Nino_4$X_boot[,,iwk])
    SST_sieve_lrc_boot_Nino_3_4[,,iwk] = long_run_covariance_estimation(dat = SST_sieve_Nino_3_4$X_boot[,,iwk])
    print(iwk); rm(iwk)
}

############
# quantiles
############

# Nino 1+2

SST_sieve_lrc_boot_80_lb_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_Nino_1_2 = apply(SST_sieve_lrc_boot_Nino_1_2, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_Nino_1_2 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_80_ub_Nino_1_2 < SST_lrc_sample)))/(12 * 12), 4) # 0.7917

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_Nino_1_2 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_95_ub_Nino_1_2 < SST_lrc_sample)))/(12 * 12), 4) # 0.9028

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_Nino_1_2 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_99_ub_Nino_1_2 < SST_lrc_sample)))/(12 * 12), 4) # 0.9444

savepdf("Fig_4a", width = 12, height = 10, toplines = 0.8)
persp(1:12, 1:12, SST_lrc_sample_Nino_1_2, theta = 45, phi = 45, ticktype = "detailed",
      xlab = "Month", ylab = "Month", zlab = "")
dev.off()

savepdf("Fig_4b", width = 12, height = 10, toplines = 0.8)
persp(1:12, 1:12, SST_sieve_lrc_boot_80_lb_Nino_1_2, theta = 45, phi = 45, ticktype = "detailed",
      xlab = "Month", ylab = "Month", zlab = "")
dev.off()

savepdf("Fig_4c", width = 12, height = 10, toplines = 0.8)
persp(1:12, 1:12, SST_sieve_lrc_boot_80_ub_Nino_1_2, theta = 45, phi = 45, ticktype = "detailed",
      xlab = "Month", ylab = "Month", zlab = "")
dev.off()

# Nino_3

SST_sieve_lrc_boot_80_lb_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_Nino_3 = apply(SST_sieve_lrc_boot_Nino_3, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_Nino_3 > SST_lrc_sample_Nino_3)) +
               length(which(SST_sieve_lrc_boot_80_ub_Nino_3 < SST_lrc_sample_Nino_3)))/(12 * 12), 4) # 0.5

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_Nino_3 > SST_lrc_sample_Nino_3)) +
               length(which(SST_sieve_lrc_boot_95_ub_Nino_3 < SST_lrc_sample_Nino_3)))/(12 * 12), 4) # 0.7361

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_Nino_3 > SST_lrc_sample_Nino_3)) +
               length(which(SST_sieve_lrc_boot_99_ub_Nino_3 < SST_lrc_sample_Nino_3)))/(12 * 12), 4) # 0.9583

# Nino_4

SST_sieve_lrc_boot_80_lb_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_Nino_4 = apply(SST_sieve_lrc_boot_Nino_4, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_Nino_4 > SST_lrc_sample_Nino_4)) +
               length(which(SST_sieve_lrc_boot_80_ub_Nino_4 < SST_lrc_sample_Nino_4)))/(12 * 12), 4) # 0.9028

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_Nino_4 > SST_lrc_sample_Nino_4)) +
               length(which(SST_sieve_lrc_boot_95_ub_Nino_4 < SST_lrc_sample_Nino_4)))/(12 * 12), 4) # 1

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_Nino_4 > SST_lrc_sample_Nino_4)) +
               length(which(SST_sieve_lrc_boot_99_ub_Nino_4 < SST_lrc_sample_Nino_4)))/(12 * 12), 4) # 1

# Nino 3+4

SST_sieve_lrc_boot_80_lb_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_Nino_3_4 = apply(SST_sieve_lrc_boot_Nino_3_4, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_Nino_3_4 > SST_lrc_sample_Nino_3_4)) +
           length(which(SST_sieve_lrc_boot_80_ub_Nino_3_4 < SST_lrc_sample_Nino_3_4)))/(12 * 12), 4) # 1

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_Nino_3_4 > SST_lrc_sample_Nino_3_4)) +
               length(which(SST_sieve_lrc_boot_95_ub_Nino_3_4 < SST_lrc_sample_Nino_3_4)))/(12 * 12), 4) # 1

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_Nino_3_4 > SST_lrc_sample_Nino_3_4)) +
               length(which(SST_sieve_lrc_boot_99_ub_Nino_3_4 < SST_lrc_sample_Nino_3_4)))/(12 * 12), 4) # 1


##################################
# sensitivity analysis of no_boot
# B = 1000
##################################

SST_sieve_B_1000 = sieve_bootstrap(fun_dat = t(SST_data), grid = 1:12,
                            ncomp_porder_selection = "eigenratio_AICC",
                            B = 1000, burn_in = 50)

SST_sieve_lrc_boot_B_1000 = array(NA, dim = c(12, 12, 1000))
for(iwk in 1:1000)
{
    SST_sieve_lrc_boot_B_1000[,,iwk] = long_run_covariance_estimation(dat = SST_sieve_B_1000$X_boot[,,iwk])
    print(iwk); rm(iwk)
}

# quantiles

SST_sieve_lrc_boot_80_lb_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_B_1000 = apply(SST_sieve_lrc_boot_B_1000, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_B_1000 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_80_ub_B_1000 < SST_lrc_sample)))/(12 * 12), 4) # 0.7778

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_B_1000 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_95_ub_B_1000 < SST_lrc_sample)))/(12 * 12), 4) # 0.9028

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_B_1000 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_99_ub_B_1000 < SST_lrc_sample)))/(12 * 12), 4) # 0.9722

##################################
# sensitivity analysis of no_boot
# B = 1500
##################################

SST_sieve_B_1500 = sieve_bootstrap(fun_dat = t(SST_data), grid = 1:12,
                                   ncomp_porder_selection = "eigenratio_AICC",
                                   B = 1500, burn_in = 50)

SST_sieve_lrc_boot_B_1500 = array(NA, dim = c(12, 12, 1500))
for(iwk in 1:1500)
{
    SST_sieve_lrc_boot_B_1500[,,iwk] = long_run_covariance_estimation(dat = SST_sieve_B_1500$X_boot[,,iwk])
    print(iwk); rm(iwk)
}

# quantiles

SST_sieve_lrc_boot_80_lb_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.1, 0.9))[1,,]
SST_sieve_lrc_boot_80_ub_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.1, 0.9))[2,,]

SST_sieve_lrc_boot_95_lb_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.025, 0.975))[1,,]
SST_sieve_lrc_boot_95_ub_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.025, 0.975))[2,,]

SST_sieve_lrc_boot_99_lb_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.005, 0.995))[1,,]
SST_sieve_lrc_boot_99_ub_B_1500 = apply(SST_sieve_lrc_boot_B_1500, c(1, 2), quantile, c(0.005, 0.995))[2,,]

round(1 - (length(which(SST_sieve_lrc_boot_80_lb_B_1500 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_80_ub_B_1500 < SST_lrc_sample)))/(12 * 12), 4) # 0.7917

round(1 - (length(which(SST_sieve_lrc_boot_95_lb_B_1500 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_95_ub_B_1500 < SST_lrc_sample)))/(12 * 12), 4) # 0.8889

round(1 - (length(which(SST_sieve_lrc_boot_99_lb_B_1500 > SST_lrc_sample)) +
           length(which(SST_sieve_lrc_boot_99_ub_B_1500 < SST_lrc_sample)))/(12 * 12), 4) # 0.9444
