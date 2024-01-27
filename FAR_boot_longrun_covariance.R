###################
# FAR(1) bootstrap
# 200 replications
###################

#########
# MA_1_0
#########

registerDoMC(10)

# FAR(1)

FAR_MA_1_0_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR",
                                                            seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_MA_1_0_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR",
                                                            seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_MA_1_0_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR",
                                                            seed_number = iwk, bootnumber = 400, samplesize = 500)

# sieve

sieve_MA_1_0_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_MA_1_0_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_1_0_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

# FAR(1) + sieve

FAR_sieve_MA_1_0_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                  seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_MA_1_0_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                  seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_MA_1_0_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                  seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

FAR_MA_1_0_n100_dist = sieve_MA_1_0_n100_dist = FAR_sieve_MA_1_0_n100_dist =
FAR_MA_1_0_n300_dist = sieve_MA_1_0_n300_dist = FAR_sieve_MA_1_0_n300_dist =
FAR_MA_1_0_n500_dist = sieve_MA_1_0_n500_dist = FAR_sieve_MA_1_0_n500_dist = matrix(NA, 200, 10)
FAR_MA_1_0_n100_dist_true = FAR_MA_1_0_n300_dist_true = FAR_MA_1_0_n500_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    FAR_MA_1_0_n100_dist[ik,] = FAR_MA_1_0_n100[[ik]]$dist
    FAR_MA_1_0_n300_dist[ik,] = FAR_MA_1_0_n300[[ik]]$dist
    FAR_MA_1_0_n500_dist[ik,] = FAR_MA_1_0_n500[[ik]]$dist

    FAR_MA_1_0_n100_dist_true[ik] = FAR_MA_1_0_n100[[ik]]$disttrue
    FAR_MA_1_0_n300_dist_true[ik] = FAR_MA_1_0_n300[[ik]]$disttrue
    FAR_MA_1_0_n500_dist_true[ik] = FAR_MA_1_0_n500[[ik]]$disttrue
    print(ik); rm(ik)
}

for(ik in 1:200)
{
    sieve_MA_1_0_n100_dist[ik,] = sieve_MA_1_0_n100[[ik]]$dist
    sieve_MA_1_0_n300_dist[ik,] = sieve_MA_1_0_n300[[ik]]$dist
    sieve_MA_1_0_n500_dist[ik,] = sieve_MA_1_0_n500[[ik]]$dist
    print(ik); rm(ik)
}

sieve_MA_1_0_n500_dist[-146,]

for(ik in 113:200)
{
    FAR_sieve_MA_1_0_n100_dist[ik,] = FAR_sieve_MA_1_0_n100[[ik]]$dist
    FAR_sieve_MA_1_0_n300_dist[ik,] = FAR_sieve_MA_1_0_n300[[ik]]$dist
    FAR_sieve_MA_1_0_n500_dist[ik,] = FAR_sieve_MA_1_0_n500[[ik]]$dist
    print(ik); rm(ik)
}

index = c(38, 73, 112)

colnames(FAR_MA_1_0_n100_dist) = colnames(FAR_MA_1_0_n300_dist) = colnames(FAR_MA_1_0_n500_dist) =
colnames(sieve_MA_1_0_n100_dist) = colnames(sieve_MA_1_0_n300_dist) = colnames(sieve_MA_1_0_n500_dist) =
colnames(FAR_sieve_MA_1_0_n100_dist) = colnames(FAR_sieve_MA_1_0_n300_dist) = colnames(FAR_sieve_MA_1_0_n500_dist) = seq(0.5, 0.95, by = 0.05)

rownames(FAR_MA_1_0_n100_dist) = rownames(FAR_MA_1_0_n300_dist) = rownames(FAR_MA_1_0_n500_dist) =
rownames(sieve_MA_1_0_n100_dist) = rownames(sieve_MA_1_0_n300_dist) = rownames(sieve_MA_1_0_n500_dist) =
rownames(FAR_sieve_MA_1_0_n100_dist) = rownames(FAR_sieve_MA_1_0_n300_dist) = rownames(FAR_sieve_MA_1_0_n500_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 = FAR_count_measure_n500 =
sieve_count_measure_n100 = sieve_count_measure_n300 = sieve_count_measure_n500 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = FAR_sieve_count_measure_n500 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_MA_1_0_n100_dist_true[ik] <= FAR_MA_1_0_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_MA_1_0_n300_dist_true[ik] <= FAR_MA_1_0_n300_dist[ik,iw], 1, 0)
        FAR_count_measure_n500[ik,iw] = ifelse(FAR_MA_1_0_n500_dist_true[ik] <= FAR_MA_1_0_n500_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_1_0_n100_dist_true[ik] <= sieve_MA_1_0_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_1_0_n300_dist_true[ik] <= sieve_MA_1_0_n300_dist[ik,iw], 1, 0)
        sieve_count_measure_n500[ik,iw] = ifelse(FAR_MA_1_0_n500_dist_true[ik] <= sieve_MA_1_0_n500_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_1_0_n100_dist_true[ik] <= FAR_sieve_MA_1_0_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_1_0_n300_dist_true[ik] <= FAR_sieve_MA_1_0_n300_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n500[ik,iw] = ifelse(FAR_MA_1_0_n500_dist_true[ik] <= FAR_sieve_MA_1_0_n500_dist[ik,iw], 1, 0)
    }
}
FAR_MA_1_0_n100_EC = colSums(FAR_count_measure_n100)/200   # 0.630 0.670 0.750 0.795 0.800 0.850 0.895 0.930 0.960 0.980
FAR_MA_1_0_n300_EC = colSums(FAR_count_measure_n300)/200   # 0.705 0.760 0.785 0.850 0.890 0.915 0.945 0.980 0.985 0.995
FAR_MA_1_0_n500_EC = colSums(FAR_count_measure_n500)/200   # 0.745 0.795 0.830 0.875 0.905 0.940 0.955 0.975 0.995 1.000

sieve_MA_1_0_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_MA_1_0_n300_EC = colSums(sieve_count_measure_n300)/200
sieve_MA_1_0_n500_EC = colSums(sieve_count_measure_n500[-146,])/(200-1)

FAR_sieve_MA_1_0_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_MA_1_0_n300_EC = colSums(FAR_sieve_count_measure_n300)/200
FAR_sieve_MA_1_0_n500_EC = colSums(FAR_sieve_count_measure_n500[-c(38, 73, 112),])/(200-3)

MA_1_0_results_EC_mat = cbind(FAR_MA_1_0_n100_EC, sieve_MA_1_0_n100_EC, FAR_sieve_MA_1_0_n100_EC,
                              FAR_MA_1_0_n300_EC, sieve_MA_1_0_n300_EC, FAR_sieve_MA_1_0_n300_EC)
MA_1_0_results_EC = rbind(MA_1_0_results_EC_mat, apply(MA_1_0_results_EC_mat, 2, mean),
                          apply(MA_1_0_results_EC_mat, 2, median))
rownames(MA_1_0_results_EC) = c(grid_prob, "Mean", "Median")
colnames(MA_1_0_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)",
                                "FAR (n=500)", "sieve (n=500)", "FAR_sieve (n=500)")

# plot empirical coverage probability

savepdf("MA_1_0_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_1_0_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0.5, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[1](0)), col = 4)
points(grid_prob, MA_1_0_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
legend("topleft", c("FAR bootstrap", "Sieve bootstrap", "Target"),
       col = c(4, 2, 1), pch = c(18, 17, 16), cex = 0.8)
dev.off()

savepdf("MA_1_0_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_1_0_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0.5, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[1](0)), col = 4)
points(grid_prob, MA_1_0_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
legend("topleft", c("FAR bootstrap", "Sieve bootstrap", "Target"),
       col = c(4, 2, 1), pch = c(18, 17, 16), cex = 0.8)
dev.off()

# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_MA_1_0_n100_CPD = abs(FAR_MA_1_0_n100_EC - grid_prob)
FAR_MA_1_0_n300_CPD = abs(FAR_MA_1_0_n300_EC - grid_prob)

sieve_MA_1_0_n100_CPD = abs(sieve_MA_1_0_n100_EC - grid_prob)
sieve_MA_1_0_n300_CPD = abs(sieve_MA_1_0_n300_EC - grid_prob)

FAR_sieve_MA_1_0_n100_CPD = abs(FAR_sieve_MA_1_0_n100_EC - grid_prob)
FAR_sieve_MA_1_0_n300_CPD = abs(FAR_sieve_MA_1_0_n300_EC - grid_prob)

MA_1_0_results_mat = cbind(FAR_MA_1_0_n100_CPD, FAR_MA_1_0_n300_CPD,
                           sieve_MA_1_0_n100_CPD, sieve_MA_1_0_n300_CPD,
                           FAR_sieve_MA_1_0_n100_CPD, FAR_sieve_MA_1_0_n300_CPD)
MA_1_0_results = rbind(MA_1_0_results_mat, apply(MA_1_0_results_mat, 2, mean), apply(MA_1_0_results_mat, 2, median))
rownames(MA_1_0_results) = c(grid_prob, "Mean", "Median")
colnames(MA_1_0_results) = c("FAR (n=100)", "FAR (n=300)", "FAR (n=500)",
                             "sieve (n=100)", "sieve (n=300)", "sieve (n=500)",
                             "FAR_sieve (n=100)", "FAR_sieve (n=300)", "FAR_sieve (n=500)")


###########
# MA_0.5_1
###########

# FAR

FAR_MA_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_MA_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_MA_0.5_1_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 500)

# sieve

sieve_MA_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_MA_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_0.5_1_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 500)

# FAR sieve

FAR_sieve_MA_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_MA_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_MA_0.5_1_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 500)


FAR_MA_0.5_1_n100_dist = FAR_MA_0.5_1_n300_dist =
sieve_MA_0.5_1_n100_dist = sieve_MA_0.5_1_n300_dist =
FAR_sieve_MA_0.5_1_n100_dist = FAR_sieve_MA_0.5_1_n300_dist = matrix(NA, 200, 10)
FAR_MA_0.5_1_n100_dist_true = FAR_MA_0.5_1_n300_dist_true =
sieve_MA_0.5_1_n100_dist_true = sieve_MA_0.5_1_n300_dist_true =
FAR_sieve_MA_0.5_1_n100_dist_true = FAR_sieve_MA_0.5_1_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    # FAR boot

    FAR_MA_0.5_1_n100_dist[ik,] = FAR_MA_0.5_1_n100[[ik]]$dist
    FAR_MA_0.5_1_n100_dist_true[ik] = FAR_MA_0.5_1_n100[[ik]]$disttrue

    FAR_MA_0.5_1_n300_dist[ik,] = FAR_MA_0.5_1_n300[[ik]]$dist
    FAR_MA_0.5_1_n300_dist_true[ik] = FAR_MA_0.5_1_n300[[ik]]$disttrue

    # sieve boot

    sieve_MA_0.5_1_n100_dist[ik,] = sieve_MA_0.5_1_n100[[ik]]$dist
    sieve_MA_0.5_1_n100_dist_true[ik] = sieve_MA_0.5_1_n100[[ik]]$disttrue

    sieve_MA_0.5_1_n300_dist[ik,] = sieve_MA_0.5_1_n300[[ik]]$dist
    sieve_MA_0.5_1_n300_dist_true[ik] = sieve_MA_0.5_1_n300[[ik]]$disttrue

    # FAR_sieve boot

    FAR_sieve_MA_0.5_1_n100_dist[ik,] = FAR_sieve_MA_0.5_1_n100[[ik]]$dist
    FAR_sieve_MA_0.5_1_n100_dist_true[ik] = FAR_sieve_MA_0.5_1_n100[[ik]]$disttrue

    FAR_sieve_MA_0.5_1_n300_dist[ik,] = FAR_sieve_MA_0.5_1_n300[[ik]]$dist
    FAR_sieve_MA_0.5_1_n300_dist_true[ik] = FAR_sieve_MA_0.5_1_n300[[ik]]$disttrue
}
colnames(FAR_MA_0.5_1_n100_dist) = colnames(FAR_MA_0.5_1_n300_dist) =
colnames(sieve_MA_0.5_1_n100_dist) = colnames(sieve_MA_0.5_1_n300_dist) =
colnames(FAR_sieve_MA_0.5_1_n100_dist) = colnames(FAR_sieve_MA_0.5_1_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(FAR_MA_0.5_1_n100_dist) = rownames(FAR_MA_0.5_1_n300_dist) =
rownames(sieve_MA_0.5_1_n100_dist) = rownames(sieve_MA_0.5_1_n300_dist) =
rownames(FAR_sieve_MA_0.5_1_n100_dist) = rownames(FAR_sieve_MA_0.5_1_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        # FAR_boot

        FAR_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_1_n100_dist_true[ik] <= FAR_MA_0.5_1_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_1_n300_dist_true[ik] <= FAR_MA_0.5_1_n300_dist[ik,iw], 1, 0)

        # sieve_boot

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_1_n100_dist_true[ik] <= sieve_MA_0.5_1_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_1_n300_dist_true[ik] <= sieve_MA_0.5_1_n300_dist[ik,iw], 1, 0)

        # FAR_sieve_boot

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_1_n100_dist_true[ik] <= FAR_sieve_MA_0.5_1_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_1_n300_dist_true[ik] <= FAR_sieve_MA_0.5_1_n300_dist[ik,iw], 1, 0)
    }
}

# Empirical coverage probability

FAR_MA_0.5_1_n100_EC = colSums(FAR_count_measure_n100)/200
FAR_MA_0.5_1_n300_EC = colSums(FAR_count_measure_n300)/200

sieve_MA_0.5_1_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_MA_0.5_1_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_MA_0.5_1_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_MA_0.5_1_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

MA_0.5_1_EC_results_mat = cbind(FAR_MA_0.5_1_n100_EC, sieve_MA_0.5_1_n100_EC,
                                FAR_sieve_MA_0.5_1_n100_EC, FAR_MA_0.5_1_n300_EC,
                                sieve_MA_0.5_1_n300_EC, FAR_sieve_MA_0.5_1_n300_EC)
MA_0.5_1_results_EC = rbind(MA_0.5_1_EC_results_mat, apply(MA_0.5_1_EC_results_mat, 2, mean),
                            apply(MA_0.5_1_EC_results_mat, 2, median))
rownames(MA_0.5_1_results_EC) = c(grid_prob, "Mean", "Median")
colnames(MA_0.5_1_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

savepdf("MA_0.5_1_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_1_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0.45, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](1)), col = 4)
points(grid_prob, MA_0.5_1_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

savepdf("MA_0.5_1_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_1_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0.45, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](1)), col = 4)
points(grid_prob, MA_0.5_1_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()


# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_MA_0.5_1_n100_CPD = abs(FAR_MA_0.5_1_n100_EC - grid_prob)
FAR_MA_0.5_1_n300_CPD = abs(FAR_MA_0.5_1_n300_EC - grid_prob)

sieve_MA_0.5_1_n100_CPD = abs(sieve_MA_0.5_1_n100_EC - grid_prob)
sieve_MA_0.5_1_n300_CPD = abs(sieve_MA_0.5_1_n300_EC - grid_prob)

FAR_sieve_MA_0.5_1_n100_CPD = abs(FAR_sieve_MA_0.5_1_n100_EC - grid_prob)
FAR_sieve_MA_0.5_1_n300_CPD = abs(FAR_sieve_MA_0.5_1_n300_EC - grid_prob)

MA_0.5_1_CPD_results_mat = cbind(FAR_MA_0.5_1_n100_CPD, FAR_MA_0.5_1_n300_CPD,
                                 sieve_MA_0.5_1_n100_CPD, sieve_MA_0.5_1_n300_CPD,
                                 FAR_sieve_MA_0.5_1_n100_CPD, FAR_sieve_MA_0.5_1_n300_CPD)
MA_0.5_1_results = rbind(MA_0.5_1_CPD_results_mat, apply(MA_0.5_1_CPD_results_mat, 2, mean),
                         apply(MA_0.5_1_CPD_results_mat, 2, median))
rownames(MA_0.5_1_results) = c(grid_prob, "Mean", "Median")
colnames(MA_0.5_1_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)",
                               "FAR_sieve (n=100)", "FAR_sieve (n=300)")

###########
# MA_0.5_4
###########

registerDoMC(16)
theta = 0.5

# FAR

FAR_MA_0.5_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_MA_0.5_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_MA_0.5_4_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR",
                                                                  seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

# sieve

sieve_MA_0.5_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_MA_0.5_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_0.5_4_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

# FAR_sieve

FAR_sieve_MA_0.5_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_MA_0.5_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_MA_0.5_4_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                        seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)


FAR_MA_0.5_4_n100_dist = FAR_MA_0.5_4_n300_dist =
sieve_MA_0.5_4_n100_dist = sieve_MA_0.5_4_n300_dist =
FAR_sieve_MA_0.5_4_n100_dist = FAR_sieve_MA_0.5_4_n300_dist = matrix(NA, 200, 10)
FAR_MA_0.5_4_n100_dist_true = FAR_MA_0.5_4_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    # FAR

    FAR_MA_0.5_4_n100_dist[ik,] = FAR_MA_0.5_4_n100[[ik]]$dist
    FAR_MA_0.5_4_n100_dist_true[ik] = FAR_MA_0.5_4_n100[[ik]]$disttrue

    FAR_MA_0.5_4_n300_dist[ik,] = FAR_MA_0.5_4_n300[[ik]]$dist
    FAR_MA_0.5_4_n300_dist_true[ik] = FAR_MA_0.5_4_n300[[ik]]$disttrue

    # sieve

    sieve_MA_0.5_4_n100_dist[ik,] = sieve_MA_0.5_4_n100[[ik]]$dist
    sieve_MA_0.5_4_n300_dist[ik,] = sieve_MA_0.5_4_n300[[ik]]$dist

    # FAR_sieve

    FAR_sieve_MA_0.5_4_n100_dist[ik,] = FAR_sieve_MA_0.5_4_n100[[ik]]$dist
    FAR_sieve_MA_0.5_4_n300_dist[ik,] = FAR_sieve_MA_0.5_4_n300[[ik]]$dist
    print(ik); rm(ik)
}
colnames(FAR_MA_0.5_4_n100_dist) = colnames(FAR_MA_0.5_4_n300_dist) =
colnames(sieve_MA_0.5_4_n100_dist) = colnames(sieve_MA_0.5_4_n300_dist) =
colnames(FAR_sieve_MA_0.5_4_n100_dist) = colnames(FAR_sieve_MA_0.5_4_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(FAR_MA_0.5_4_n100_dist) = rownames(FAR_MA_0.5_4_n300_dist) =
rownames(sieve_MA_0.5_4_n100_dist) = rownames(sieve_MA_0.5_4_n300_dist) =
rownames(FAR_sieve_MA_0.5_4_n100_dist) = rownames(FAR_sieve_MA_0.5_4_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_4_n100_dist_true[ik] <= FAR_MA_0.5_4_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_4_n300_dist_true[ik] <= FAR_MA_0.5_4_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_4_n100_dist_true[ik] <= sieve_MA_0.5_4_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_4_n300_dist_true[ik] <= sieve_MA_0.5_4_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_4_n100_dist_true[ik] <= FAR_sieve_MA_0.5_4_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_4_n300_dist_true[ik] <= FAR_sieve_MA_0.5_4_n300_dist[ik,iw], 1, 0)
    }
}
FAR_MA_0.5_4_n100_EC = colSums(FAR_count_measure_n100)/200   # 0.185 0.190 0.200 0.200 0.215 0.225 0.250 0.265 0.295 0.320
FAR_MA_0.5_4_n300_EC = colSums(FAR_count_measure_n300)/200   # 0.415 0.430 0.455 0.460 0.470 0.500 0.515 0.545 0.555 0.595

sieve_MA_0.5_4_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_MA_0.5_4_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_MA_0.5_4_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_MA_0.5_4_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

MA_0.5_4_EC_results_mat = cbind(FAR_MA_0.5_4_n100_EC, sieve_MA_0.5_4_n100_EC,
                                FAR_sieve_MA_0.5_4_n100_EC, FAR_MA_0.5_4_n300_EC,
                                sieve_MA_0.5_4_n300_EC, FAR_sieve_MA_0.5_4_n300_EC)
MA_0.5_4_results_EC = rbind(MA_0.5_4_EC_results_mat, apply(MA_0.5_4_EC_results_mat, 2, mean),
                            apply(MA_0.5_4_EC_results_mat, 2, median))
rownames(MA_0.5_4_results_EC) = c(grid_prob, "Mean", "Median")
colnames(MA_0.5_4_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("MA_0.5_4_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_4_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0.1, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](4)), col = 4)
points(grid_prob, MA_0.5_4_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

savepdf("MA_0.5_4_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_4_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0.25, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](4)), col = 4)
points(grid_prob, MA_0.5_4_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_MA_0.5_4_n100_CPD = abs(FAR_MA_0.5_4_n100_EC - grid_prob)
FAR_MA_0.5_4_n300_CPD = abs(FAR_MA_0.5_4_n300_EC - grid_prob)

sieve_MA_0.5_4_n100_CPD = abs(sieve_MA_0.5_4_n100_EC - grid_prob)
sieve_MA_0.5_4_n300_CPD = abs(sieve_MA_0.5_4_n300_EC - grid_prob)

FAR_sieve_MA_0.5_4_n100_CPD = abs(FAR_sieve_MA_0.5_4_n100_EC - grid_prob)
FAR_sieve_MA_0.5_4_n300_CPD = abs(FAR_sieve_MA_0.5_4_n300_EC - grid_prob)

MA_0.5_4_results_mat = cbind(FAR_MA_0.5_4_n100_CPD, FAR_MA_0.5_4_n300_CPD,
                            sieve_MA_0.5_4_n100_CPD, sieve_MA_0.5_4_n300_CPD,
                            FAR_sieve_MA_0.5_4_n100_CPD, FAR_sieve_MA_0.5_4_n300_CPD)
MA_0.5_4_results = rbind(MA_0.5_4_results_mat, colMeans(MA_0.5_4_results_mat), apply(MA_0.5_4_results_mat, 2, median))

mean(FAR_MA_0.5_4_n100_CPD)  # 0.4905
mean(FAR_MA_0.5_4_n300_CPD)  # 0.231

###########
# MA_0.5_8
###########

# FAR

FAR_MA_0.5_8_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_MA_0.5_8_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "FAR",
                                                            seed_number = iwk, bootnumber = 400, samplesize = 300)

# Sieve

sieve_MA_0.5_8_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_MA_0.5_8_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

# FAR+Sieve

FAR_sieve_MA_0.5_8_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_MA_0.5_8_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)


FAR_MA_0.5_8_n100_dist = FAR_MA_0.5_8_n300_dist =
sieve_MA_0.5_8_n100_dist = sieve_MA_0.5_8_n300_dist =
FAR_sieve_MA_0.5_8_n100_dist = FAR_sieve_MA_0.5_8_n300_dist = matrix(NA, 200, 10)
FAR_MA_0.5_8_n100_dist_true = FAR_MA_0.5_8_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    # FAR

    FAR_MA_0.5_8_n100_dist[ik,] = FAR_MA_0.5_8_n100[[ik]]$dist
    FAR_MA_0.5_8_n100_dist_true[ik] = FAR_MA_0.5_8_n100[[ik]]$disttrue

    FAR_MA_0.5_8_n300_dist[ik,] = FAR_MA_0.5_8_n300[[ik]]$dist
    FAR_MA_0.5_8_n300_dist_true[ik] = FAR_MA_0.5_8_n300[[ik]]$disttrue

    # sieve

    sieve_MA_0.5_8_n100_dist[ik,] = sieve_MA_0.5_8_n100[[ik]]$dist
    sieve_MA_0.5_8_n300_dist[ik,] = sieve_MA_0.5_8_n300[[ik]]$dist

    # FAR_sieve

    FAR_sieve_MA_0.5_8_n100_dist[ik,] = FAR_sieve_MA_0.5_8_n100[[ik]]$dist
    FAR_sieve_MA_0.5_8_n300_dist[ik,] = FAR_sieve_MA_0.5_8_n300[[ik]]$dist
    print(ik); rm(ik)
}

colnames(FAR_MA_0.5_8_n100_dist) = colnames(FAR_MA_0.5_8_n300_dist) =
colnames(sieve_MA_0.5_8_n100_dist) = colnames(sieve_MA_0.5_8_n300_dist) =
colnames(FAR_sieve_MA_0.5_8_n100_dist) = colnames(FAR_sieve_MA_0.5_8_n300_dist) = seq(0.5, 0.95, by = 0.05)

rownames(FAR_MA_0.5_8_n100_dist) = rownames(FAR_MA_0.5_8_n300_dist) =
rownames(sieve_MA_0.5_8_n100_dist) = rownames(sieve_MA_0.5_8_n300_dist) =
rownames(FAR_sieve_MA_0.5_8_n100_dist) = rownames(FAR_sieve_MA_0.5_8_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_8_n100_dist_true[ik] <= FAR_MA_0.5_8_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_8_n300_dist_true[ik] <= FAR_MA_0.5_8_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_8_n100_dist_true[ik] <= sieve_MA_0.5_8_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_8_n300_dist_true[ik] <= sieve_MA_0.5_8_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_0.5_8_n100_dist_true[ik] <= FAR_sieve_MA_0.5_8_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_0.5_8_n300_dist_true[ik] <= FAR_sieve_MA_0.5_8_n300_dist[ik,iw], 1, 0)
    }
}
FAR_MA_0.5_8_n100_EC = colSums(FAR_count_measure_n100)/200
FAR_MA_0.5_8_n300_EC = colSums(FAR_count_measure_n300)/200

sieve_MA_0.5_8_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_MA_0.5_8_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_MA_0.5_8_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_MA_0.5_8_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

MA_0.5_8_EC_results_mat = cbind(FAR_MA_0.5_8_n100_EC, sieve_MA_0.5_8_n100_EC, FAR_sieve_MA_0.5_8_n100_EC,
                                FAR_MA_0.5_8_n300_EC, sieve_MA_0.5_8_n300_EC, FAR_sieve_MA_0.5_8_n300_EC)
MA_0.5_8_results_EC = rbind(MA_0.5_8_EC_results_mat, apply(MA_0.5_8_EC_results_mat, 2, mean),
                            apply(MA_0.5_8_EC_results_mat, 2, median))
rownames(MA_0.5_8_results_EC) = c(grid_prob, "Mean", "Median")
colnames(MA_0.5_8_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("MA_0.5_8_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_8_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](8)), col = 4)
points(grid_prob, MA_0.5_8_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

savepdf("MA_0.5_8_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_0.5_8_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FMA[0.5](8)), col = 4)
points(grid_prob, MA_0.5_8_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_MA_0.5_8_n100_CPD = abs(FAR_MA_0.5_8_n100_EC - grid_prob)
FAR_MA_0.5_8_n300_CPD = abs(FAR_MA_0.5_8_n300_EC - grid_prob)

sieve_MA_0.5_8_n100_CPD = abs(sieve_MA_0.5_8_n100_EC - grid_prob)
sieve_MA_0.5_8_n300_CPD = abs(sieve_MA_0.5_8_n300_EC - grid_prob)

FAR_sieve_MA_0.5_8_n100_CPD = abs(FAR_sieve_MA_0.5_8_n100_EC - grid_prob)
FAR_sieve_MA_0.5_8_n300_CPD = abs(FAR_sieve_MA_0.5_8_n300_EC - grid_prob)

MA_0.5_8_CPD_results_mat = cbind(FAR_MA_0.5_8_n100_CPD, FAR_MA_0.5_8_n300_CPD,
                                 sieve_MA_0.5_8_n100_CPD, sieve_MA_0.5_8_n300_CPD,
                                 FAR_sieve_MA_0.5_8_n100_CPD, FAR_sieve_MA_0.5_8_n300_CPD)
MA_0.5_8_results = rbind(MA_0.5_8_CPD_results_mat, apply(MA_0.5_8_CPD_results_mat, 2, mean),
                         apply(MA_0.5_8_CPD_results_mat, 2, median))
rownames(MA_0.5_8_results) = c(grid_prob, "Mean", "Median")
colnames(MA_0.5_8_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)",
                               "FAR_sieve (n=100)", "FAR_sieve (n=300)")

MA_process_results = rbind(MA_1_0_results[11,],
                            MA_0.5_1_results[11,],
                            MA_0.5_4_results[11,],
                            MA_0.5_8_results[11,])
colnames(MA_process_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)", "FAR+sieve (n=100)", "FAR+sieve (n=300)")
rownames(MA_process_results) = c("MA_0", "MA_1", "MA_4", "MA_8")

############
# FAR_0.5_1
############

# FAR

FAR_AR_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_AR_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_AR_0.5_1_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR",
                                                                  seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

# sieve

sieve_AR_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_AR_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_AR_0.5_1_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

# FAR_sieve

FAR_sieve_AR_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_AR_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_AR_0.5_1_n500 = foreach(iwk = 1:200) %dopar% try(bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 500), silent = TRUE)

FAR_AR_0.5_1_n100_dist = FAR_AR_0.5_1_n300_dist = FAR_AR_0.5_1_n500_dist =
sieve_AR_0.5_1_n100_dist = sieve_AR_0.5_1_n300_dist = sieve_AR_0.5_1_n500_dist =
FAR_sieve_AR_0.5_1_n100_dist = FAR_sieve_AR_0.5_1_n300_dist = FAR_sieve_AR_0.5_1_n500_dist = matrix(NA, 200, 10)
FAR_AR_0.5_1_n100_dist_true = FAR_AR_0.5_1_n300_dist_true = FAR_AR_0.5_1_n500_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    # FAR

    FAR_AR_0.5_1_n100_dist[ik,] = FAR_AR_0.5_1_n100[[ik]]$dist
    FAR_AR_0.5_1_n100_dist_true[ik] = FAR_AR_0.5_1_n100[[ik]]$disttrue

    FAR_AR_0.5_1_n300_dist[ik,] = FAR_AR_0.5_1_n300[[ik]]$dist
    FAR_AR_0.5_1_n300_dist_true[ik] = FAR_AR_0.5_1_n300[[ik]]$disttrue

    FAR_AR_0.5_1_n500_dist[ik,] = FAR_AR_0.5_1_n500[[ik]]$dist
    FAR_AR_0.5_1_n500_dist_true[ik] = FAR_AR_0.5_1_n500[[ik]]$disttrue

    # sieve

    sieve_AR_0.5_1_n100_dist[ik,] = sieve_AR_0.5_1_n100[[ik]]$dist
    sieve_AR_0.5_1_n300_dist[ik,] = sieve_AR_0.5_1_n300[[ik]]$dist
    sieve_AR_0.5_1_n500_dist[ik,] = sieve_AR_0.5_1_n500[[ik]]$dist

    # FAR_sieve

    FAR_sieve_AR_0.5_1_n100_dist[ik,] = FAR_sieve_AR_0.5_1_n100[[ik]]$dist
    FAR_sieve_AR_0.5_1_n300_dist[ik,] = FAR_sieve_AR_0.5_1_n300[[ik]]$dist
    FAR_sieve_AR_0.5_1_n500_dist[ik,] = FAR_sieve_AR_0.5_1_n500[[ik]]$dist
    print(ik); rm(ik)
}
colnames(FAR_AR_0.5_1_n100_dist) = colnames(FAR_AR_0.5_1_n300_dist) = colnames(FAR_AR_0.5_1_n500_dist) =
colnames(sieve_AR_0.5_1_n100_dist) = colnames(sieve_AR_0.5_1_n300_dist) = colnames(sieve_AR_0.5_1_n500_dist) =
colnames(FAR_sieve_AR_0.5_1_n100_dist) = colnames(FAR_sieve_AR_0.5_1_n300_dist) = colnames(FAR_sieve_AR_0.5_1_n500_dist) = seq(0.5, 0.95, by = 0.05)
rownames(FAR_AR_0.5_1_n100_dist) = rownames(FAR_AR_0.5_1_n300_dist) = rownames(FAR_AR_0.5_1_n500_dist) =
rownames(sieve_AR_0.5_1_n100_dist) = rownames(sieve_AR_0.5_1_n300_dist) = rownames(sieve_AR_0.5_1_n500_dist) =
rownames(FAR_sieve_AR_0.5_1_n100_dist) = rownames(FAR_sieve_AR_0.5_1_n300_dist) = rownames(FAR_sieve_AR_0.5_1_n500_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.5_1_n100_dist_true[ik] <= FAR_AR_0.5_1_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.5_1_n300_dist_true[ik] <= FAR_AR_0.5_1_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.5_1_n100_dist_true[ik] <= sieve_AR_0.5_1_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.5_1_n300_dist_true[ik] <= sieve_AR_0.5_1_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.5_1_n100_dist_true[ik] <= FAR_sieve_AR_0.5_1_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.5_1_n300_dist_true[ik] <= FAR_sieve_AR_0.5_1_n300_dist[ik,iw], 1, 0)
    }
    print(ik); rm(ik)
}
FAR_AR_0.5_1_n100_EC = colSums(FAR_count_measure_n100)/200   # 0.455 0.480 0.500 0.510 0.530 0.555 0.570 0.575 0.590 0.605
FAR_AR_0.5_1_n300_EC = colSums(FAR_count_measure_n300)/200   # 0.705 0.715 0.745 0.765 0.785 0.785 0.815 0.830 0.860 0.865

sieve_AR_0.5_1_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_AR_0.5_1_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_AR_0.5_1_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_AR_0.5_1_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

AR_0.5_1_EC_results_mat = cbind(FAR_AR_0.5_1_n100_EC, sieve_AR_0.5_1_n100_EC,
                                    FAR_sieve_AR_0.5_1_n100_EC, FAR_AR_0.5_1_n300_EC,
                                    sieve_AR_0.5_1_n300_EC, FAR_sieve_AR_0.5_1_n300_EC)
AR_0.5_1_results_EC = rbind(AR_0.5_1_EC_results_mat, apply(AR_0.5_1_EC_results_mat, 2, mean),
                            apply(AR_0.5_1_EC_results_mat, 2, median))
rownames(AR_0.5_1_results_EC) = c(grid_prob, "Mean", "Median")
colnames(AR_0.5_1_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("AR_0.5_1_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_0.5_1_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[0.5](1)), col = 4)
points(grid_prob, AR_0.5_1_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
legend("topleft", c("FAR bootstrap", "Sieve bootstrap", "Target"),
       col = c(4, 2, 1), pch = c(18, 17, 16), cex = 0.8)
dev.off()

savepdf("AR_0.5_1_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_0.5_1_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[0.5](1)), col = 4)
points(grid_prob, AR_0.5_1_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
legend("topleft", c("FAR bootstrap", "Sieve bootstrap", "Target"),
       col = c(4, 2, 1), pch = c(18, 17, 16), cex = 0.8)
dev.off()

# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_AR_0.5_1_n100_CPD = abs(FAR_AR_0.5_1_n100_EC - grid_prob)
FAR_AR_0.5_1_n300_CPD = abs(FAR_AR_0.5_1_n300_EC - grid_prob)

sieve_AR_0.5_1_n100_CPD = abs(sieve_AR_0.5_1_n100_EC - grid_prob)
sieve_AR_0.5_1_n300_CPD = abs(sieve_AR_0.5_1_n300_EC - grid_prob)

FAR_sieve_AR_0.5_1_n100_CPD = abs(FAR_sieve_AR_0.5_1_n100_EC - grid_prob)
FAR_sieve_AR_0.5_1_n300_CPD = abs(FAR_sieve_AR_0.5_1_n300_EC - grid_prob)

AR_0.5_1_CPD_results_mat = cbind(FAR_AR_0.5_1_n100_CPD, FAR_AR_0.5_1_n300_CPD,
                                 sieve_AR_0.5_1_n100_CPD, sieve_AR_0.5_1_n300_CPD,
                                 FAR_sieve_AR_0.5_1_n100_CPD, FAR_sieve_AR_0.5_1_n300_CPD)
AR_0.5_1_results = rbind(AR_0.5_1_CPD_results_mat, apply(AR_0.5_1_CPD_results_mat, 2, mean),
                         apply(AR_0.5_1_CPD_results_mat, 2, median))
rownames(AR_0.5_1_results) = c(grid_prob, "Mean", "Median")
colnames(AR_0.5_1_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)",
                               "FAR_sieve (n=100)", "FAR_sieve (n=300)")

mean(FAR_AR_0.5_1_n100_CPD)  # 0.188
mean(FAR_AR_0.5_1_n300_CPD)  # 0.091

############
# FAR_0.3_2
############

# FAR

FAR_AR_0.3_2_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "FAR",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_AR_0.3_2_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "FAR",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

# sieve

sieve_AR_0.3_2_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_AR_0.3_2_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

# FAR_sieve

FAR_sieve_AR_0.3_2_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_AR_0.3_2_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_AR_0.3_2_n100_dist = FAR_AR_0.3_2_n300_dist =
sieve_AR_0.3_2_n100_dist = sieve_AR_0.3_2_n300_dist =
FAR_sieve_AR_0.3_2_n100_dist = FAR_sieve_AR_0.3_2_n300_dist = matrix(NA, 200, 10)
FAR_AR_0.3_2_n100_dist_true = FAR_AR_0.3_2_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    FAR_AR_0.3_2_n100_dist[ik,] = FAR_AR_0.3_2_n100[[ik]]$dist
    FAR_AR_0.3_2_n100_dist_true[ik] = FAR_AR_0.3_2_n100[[ik]]$disttrue

    FAR_AR_0.3_2_n300_dist[ik,] = FAR_AR_0.3_2_n300[[ik]]$dist
    FAR_AR_0.3_2_n300_dist_true[ik] = FAR_AR_0.3_2_n300[[ik]]$disttrue

    sieve_AR_0.3_2_n100_dist[ik,] = sieve_AR_0.3_2_n100[[ik]]$dist
    sieve_AR_0.3_2_n300_dist[ik,] = sieve_AR_0.3_2_n300[[ik]]$dist

    FAR_sieve_AR_0.3_2_n100_dist[ik,] = FAR_sieve_AR_0.3_2_n100[[ik]]$dist
    FAR_sieve_AR_0.3_2_n300_dist[ik,] = FAR_sieve_AR_0.3_2_n300[[ik]]$dist
    print(ik); rm(ik)
}
colnames(FAR_AR_0.3_2_n100_dist) = colnames(FAR_AR_0.3_2_n300_dist) =
colnames(sieve_AR_0.3_2_n100_dist) = colnames(sieve_AR_0.3_2_n300_dist) =
colnames(FAR_sieve_AR_0.3_2_n100_dist) = colnames(FAR_sieve_AR_0.3_2_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(FAR_AR_0.3_2_n100_dist) = rownames(FAR_AR_0.3_2_n300_dist) =
rownames(sieve_AR_0.3_2_n100_dist) = rownames(sieve_AR_0.3_2_n300_dist) =
rownames(FAR_sieve_AR_0.3_2_n100_dist) = rownames(FAR_sieve_AR_0.3_2_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.3_2_n100_dist_true[ik] <= FAR_AR_0.3_2_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.3_2_n300_dist_true[ik] <= FAR_AR_0.3_2_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.3_2_n100_dist_true[ik] <= sieve_AR_0.3_2_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.3_2_n300_dist_true[ik] <= sieve_AR_0.3_2_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_0.3_2_n100_dist_true[ik] <= FAR_sieve_AR_0.3_2_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_0.3_2_n300_dist_true[ik] <= FAR_sieve_AR_0.3_2_n300_dist[ik,iw], 1, 0)
    }
}
FAR_AR_0.3_2_n100_EC = colSums(FAR_count_measure_n100)/200   # 0.520 0.530 0.545 0.570 0.580 0.590 0.605 0.605 0.620 0.650
FAR_AR_0.3_2_n300_EC = colSums(FAR_count_measure_n300)/200   # 0.755 0.765 0.780 0.805 0.825 0.845 0.855 0.865 0.875 0.890

sieve_AR_0.3_2_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_AR_0.3_2_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_AR_0.3_2_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_AR_0.3_2_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

AR_0.3_2_EC_results_mat = cbind(FAR_AR_0.3_2_n100_EC, sieve_AR_0.3_2_n100_EC,
                                FAR_sieve_AR_0.3_2_n100_EC, FAR_AR_0.3_2_n300_EC,
                                sieve_AR_0.3_2_n300_EC, FAR_sieve_AR_0.3_2_n300_EC)
AR_0.3_2_results_EC = rbind(AR_0.3_2_EC_results_mat, apply(AR_0.3_2_EC_results_mat, 2, mean),
                            apply(AR_0.3_2_EC_results_mat, 2, median))
rownames(AR_0.3_2_results_EC) = c(grid_prob, "Mean", "Median")
colnames(AR_0.3_2_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("AR_0.3_2_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_0.3_2_results_EC[1:10,1], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[0.3](2)), col = 4)
points(grid_prob, AR_0.3_2_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

savepdf("AR_0.3_2_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_0.3_2_results_EC[1:10,4], type = "p", pch = 18, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[0.3](2)), col = 4)
points(grid_prob, AR_0.3_2_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, grid_prob, pch = 16, col = 1)
dev.off()

# CPD

grid_prob = seq(0.5, 0.95, by = 0.05)
FAR_AR_0.3_2_n100_CPD = abs(FAR_AR_0.3_2_n100_EC - grid_prob)
FAR_AR_0.3_2_n300_CPD = abs(FAR_AR_0.3_2_n300_EC - grid_prob)

sieve_AR_0.3_2_n100_CPD = abs(sieve_AR_0.3_2_n100_EC - grid_prob)
sieve_AR_0.3_2_n300_CPD = abs(sieve_AR_0.3_2_n300_EC - grid_prob)

FAR_sieve_AR_0.3_2_n100_CPD = abs(FAR_sieve_AR_0.3_2_n100_EC - grid_prob)
FAR_sieve_AR_0.3_2_n300_CPD = abs(FAR_sieve_AR_0.3_2_n300_EC - grid_prob)

AR_0.3_2_CPD_results_mat = cbind(FAR_AR_0.3_2_n100_CPD, FAR_AR_0.3_2_n300_CPD,
                                 sieve_AR_0.3_2_n100_CPD, sieve_AR_0.3_2_n300_CPD,
                                 FAR_sieve_AR_0.3_2_n100_CPD, FAR_sieve_AR_0.3_2_n300_CPD)
AR_0.3_2_results = rbind(AR_0.3_2_CPD_results_mat, apply(AR_0.3_2_CPD_results_mat, 2, mean),
                         apply(AR_0.3_2_CPD_results_mat, 2, median))
rownames(AR_0.3_2_results) = c(grid_prob, "Mean", "Median")
colnames(AR_0.3_2_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)",
                               "FAR_sieve (n=100)", "FAR_sieve (n=300)")


mean(FAR_AR_0.3_2_n100_CPD)  # 0.1475
mean(FAR_AR_0.3_2_n300_CPD)  # 0.118

############
# FAR_psi_1
############

FAR_AR_psi_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_AR_psi_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_AR_psi_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_AR_psi_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_AR_psi_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_AR_psi_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_psi_1", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                    seed_number = iwk, bootnumber = 400, samplesize = 300)


FAR_AR_psi_1_n100_dist = FAR_AR_psi_1_n300_dist =
sieve_AR_psi_1_n100_dist = sieve_AR_psi_1_n300_dist =
FAR_sieve_AR_psi_1_n100_dist = FAR_sieve_AR_psi_1_n300_dist = matrix(NA, 200, 10)
FAR_AR_psi_1_n100_dist_true = FAR_AR_psi_1_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    FAR_AR_psi_1_n100_dist[ik,] = FAR_AR_psi_1_n100[[ik]]$dist
    FAR_AR_psi_1_n100_dist_true[ik] = FAR_AR_psi_1_n100[[ik]]$disttrue

    FAR_AR_psi_1_n300_dist[ik,] = FAR_AR_psi_1_n300[[ik]]$dist
    FAR_AR_psi_1_n300_dist_true[ik] = FAR_AR_psi_1_n300[[ik]]$disttrue

    sieve_AR_psi_1_n100_dist[ik,] = sieve_AR_psi_1_n100[[ik]]$dist
    sieve_AR_psi_1_n300_dist[ik,] = sieve_AR_psi_1_n300[[ik]]$dist

    FAR_sieve_AR_psi_1_n100_dist[ik,] = FAR_sieve_AR_psi_1_n100[[ik]]$dist
    FAR_sieve_AR_psi_1_n300_dist[ik,] = FAR_sieve_AR_psi_1_n300[[ik]]$dist
    print(ik); rm(ik)
}
colnames(FAR_AR_psi_1_n100_dist) = colnames(FAR_AR_psi_1_n300_dist) =
colnames(sieve_AR_psi_1_n100_dist) = colnames(sieve_AR_psi_1_n300_dist) =
colnames(FAR_sieve_AR_psi_1_n100_dist) = colnames(FAR_sieve_AR_psi_1_n300_dist) = seq(0.5, 0.95, by = 0.05)

rownames(FAR_AR_psi_1_n100_dist) = rownames(FAR_AR_psi_1_n300_dist) =
rownames(sieve_AR_psi_1_n100_dist) = rownames(sieve_AR_psi_1_n300_dist) =
rownames(FAR_sieve_AR_psi_1_n100_dist) = rownames(FAR_sieve_AR_psi_1_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_AR_psi_1_n100_dist_true[ik] <= FAR_AR_psi_1_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_AR_psi_1_n300_dist_true[ik] <= FAR_AR_psi_1_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_psi_1_n100_dist_true[ik] <= sieve_AR_psi_1_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_psi_1_n300_dist_true[ik] <= sieve_AR_psi_1_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_AR_psi_1_n100_dist_true[ik] <= FAR_sieve_AR_psi_1_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_AR_psi_1_n300_dist_true[ik] <= FAR_sieve_AR_psi_1_n300_dist[ik,iw], 1, 0)
    }
}

# empirical coverage probability

FAR_AR_psi_1_n100_EC = colSums(FAR_count_measure_n100)/200
FAR_AR_psi_1_n300_EC = colSums(FAR_count_measure_n300)/200

sieve_AR_psi_1_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_AR_psi_1_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_AR_psi_1_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_AR_psi_1_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

AR_psi_1_EC_results_mat = cbind(FAR_AR_psi_1_n100_EC, sieve_AR_psi_1_n100_EC,
                                FAR_sieve_AR_psi_1_n100_EC, FAR_AR_psi_1_n300_EC,
                                sieve_AR_psi_1_n300_EC, FAR_sieve_AR_psi_1_n300_EC)
AR_psi_1_results_EC = rbind(AR_psi_1_EC_results_mat, apply(AR_psi_1_EC_results_mat, 2, mean),
                            apply(AR_psi_1_EC_results_mat, 2, median))
grid_prob = seq(0.5, 0.95, by = 0.05)
rownames(AR_psi_1_results_EC) = c(grid_prob, "Mean", "Median")
colnames(AR_psi_1_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("AR_psi_1_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_psi_1_results_EC[1:10,1], type = "p", pch = 16, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[psi](1)))
points(grid_prob, AR_psi_1_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, AR_psi_1_results_EC[1:10,3], pch = 18, col = 3)
points(grid_prob, grid_prob, pch = 4, col = 4)
dev.off()

savepdf("AR_psi_1_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, AR_psi_1_results_EC[1:10,4], type = "p", pch = 16, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(FAR[psi](1)))
points(grid_prob, AR_psi_1_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, AR_psi_1_results_EC[1:10,6], pch = 18, col = 3)
points(grid_prob, grid_prob, pch = 4, col = 4)
dev.off()

FAR_AR_psi_1_n100_CPD = abs(FAR_AR_psi_1_n100_EC - grid_prob)
FAR_AR_psi_1_n300_CPD = abs(FAR_AR_psi_1_n300_EC - grid_prob)

sieve_AR_psi_1_n100_CPD = abs(sieve_AR_psi_1_n100_EC - grid_prob)
sieve_AR_psi_1_n300_CPD = abs(sieve_AR_psi_1_n300_EC - grid_prob)

FAR_sieve_AR_psi_1_n100_CPD = abs(FAR_sieve_AR_psi_1_n100_EC - grid_prob)
FAR_sieve_AR_psi_1_n300_CPD = abs(FAR_sieve_AR_psi_1_n300_EC - grid_prob)

AR_psi_1_CPD_results_mat = cbind(FAR_AR_psi_1_n100_CPD, FAR_AR_psi_1_n300_CPD,
                                 sieve_AR_psi_1_n100_CPD, sieve_AR_psi_1_n300_CPD,
                                 FAR_sieve_AR_psi_1_n100_CPD, FAR_sieve_AR_psi_1_n300_CPD)
AR_psi_1_results = rbind(AR_psi_1_CPD_results_mat, apply(AR_psi_1_CPD_results_mat, 2, mean),
                         apply(AR_psi_1_CPD_results_mat, 2, median))
rownames(AR_psi_1_results) = c(grid_prob, "Mean", "Median")
colnames(AR_psi_1_results) = c("FAR (n=100)", "FAR (n=300)", "sieve (n=100)", "sieve (n=300)",
                               "FAR_sieve (n=100)", "FAR_sieve (n=300)")

AR_process_results = rbind(AR_0.5_1_results[11,],
                            AR_0.3_2_results[11,],
                            AR_psi_1_results[11,])
rownames(AR_process_results) = c("AR_0.5_1", "AR_0.3_2", "AR_psi_1")

cbind(AR_process_results[,c(1,3,5)], AR_process_results[,c(2,4,6)])

###########
# MA_psi_4
###########

FAR_MA_psi_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_MA_psi_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "FAR",
                                                              seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_psi_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

sieve_MA_psi_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_sieve_MA_psi_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 100)

FAR_sieve_MA_psi_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", boot_method = "FAR_sieve",
                                                                seed_number = iwk, bootnumber = 400, samplesize = 300)

FAR_MA_psi_4_n100_dist = FAR_MA_psi_4_n300_dist =
sieve_MA_psi_4_n100_dist = sieve_MA_psi_4_n300_dist =
FAR_sieve_MA_psi_4_n100_dist = FAR_sieve_MA_psi_4_n300_dist = matrix(NA, 200, 10)
FAR_MA_psi_4_n100_dist_true = FAR_MA_psi_4_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    FAR_MA_psi_4_n100_dist[ik,] = FAR_MA_psi_4_n100[[ik]]$dist
    FAR_MA_psi_4_n100_dist_true[ik] = FAR_MA_psi_4_n100[[ik]]$disttrue

    FAR_MA_psi_4_n300_dist[ik,] = FAR_MA_psi_4_n300[[ik]]$dist
    FAR_MA_psi_4_n300_dist_true[ik] = FAR_MA_psi_4_n300[[ik]]$disttrue

    sieve_MA_psi_4_n100_dist[ik,] = sieve_MA_psi_4_n100[[ik]]$dist
    sieve_MA_psi_4_n300_dist[ik,] = sieve_MA_psi_4_n300[[ik]]$dist

    FAR_sieve_MA_psi_4_n100_dist[ik,] = FAR_sieve_MA_psi_4_n100[[ik]]$dist
    FAR_sieve_MA_psi_4_n300_dist[ik,] = FAR_sieve_MA_psi_4_n300[[ik]]$dist
    print(ik); rm(ik)
}
colnames(FAR_MA_psi_4_n100_dist) = colnames(FAR_MA_psi_4_n300_dist) =
colnames(sieve_MA_psi_4_n100_dist) = colnames(sieve_MA_psi_4_n300_dist) =
colnames(FAR_sieve_MA_psi_4_n100_dist) = colnames(FAR_sieve_MA_psi_4_n300_dist) = seq(0.5, 0.95, by = 0.05)

rownames(FAR_MA_psi_4_n100_dist) = rownames(FAR_MA_psi_4_n300_dist) =
rownames(sieve_MA_psi_4_n100_dist) = rownames(sieve_MA_psi_4_n300_dist) =
rownames(FAR_sieve_MA_psi_4_n100_dist) = rownames(FAR_sieve_MA_psi_4_n300_dist) = 1:200

# empirical coverage

FAR_count_measure_n100 = FAR_count_measure_n300 =
sieve_count_measure_n100 = sieve_count_measure_n300 =
FAR_sieve_count_measure_n100 = FAR_sieve_count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        FAR_count_measure_n100[ik,iw] = ifelse(FAR_MA_psi_4_n100_dist_true[ik] <= FAR_MA_psi_4_n100_dist[ik,iw], 1, 0)
        FAR_count_measure_n300[ik,iw] = ifelse(FAR_MA_psi_4_n300_dist_true[ik] <= FAR_MA_psi_4_n300_dist[ik,iw], 1, 0)

        sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_psi_4_n100_dist_true[ik] <= sieve_MA_psi_4_n100_dist[ik,iw], 1, 0)
        sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_psi_4_n300_dist_true[ik] <= sieve_MA_psi_4_n300_dist[ik,iw], 1, 0)

        FAR_sieve_count_measure_n100[ik,iw] = ifelse(FAR_MA_psi_4_n100_dist_true[ik] <= FAR_sieve_MA_psi_4_n100_dist[ik,iw], 1, 0)
        FAR_sieve_count_measure_n300[ik,iw] = ifelse(FAR_MA_psi_4_n300_dist_true[ik] <= FAR_sieve_MA_psi_4_n300_dist[ik,iw], 1, 0)
    }
}

# empirical coverage probability

FAR_MA_psi_4_n100_EC = colSums(FAR_count_measure_n100)/200
FAR_MA_psi_4_n300_EC = colSums(FAR_count_measure_n300)/200

sieve_MA_psi_4_n100_EC = colSums(sieve_count_measure_n100)/200
sieve_MA_psi_4_n300_EC = colSums(sieve_count_measure_n300)/200

FAR_sieve_MA_psi_4_n100_EC = colSums(FAR_sieve_count_measure_n100)/200
FAR_sieve_MA_psi_4_n300_EC = colSums(FAR_sieve_count_measure_n300)/200

MA_psi_4_EC_results_mat = cbind(FAR_MA_psi_4_n100_EC, sieve_MA_psi_4_n100_EC,
                                FAR_sieve_MA_psi_4_n100_EC, FAR_MA_psi_4_n300_EC,
                                sieve_MA_psi_4_n300_EC, FAR_sieve_MA_psi_4_n300_EC)
MA_psi_4_results_EC = rbind(MA_psi_4_EC_results_mat, apply(MA_psi_4_EC_results_mat, 2, mean),
                            apply(MA_psi_4_EC_results_mat, 2, median))
grid_prob = seq(0.5, 0.95, by = 0.05)
rownames(MA_psi_4_results_EC) = c(grid_prob, "Mean", "Median")
colnames(MA_psi_4_results_EC) = c("FAR (n=100)", "sieve (n=100)", "FAR_sieve (n=100)",
                                  "FAR (n=300)", "sieve (n=300)", "FAR_sieve (n=300)")

# plot empirical coverage probability

savepdf("MA_psi_4_n100", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_psi_4_results_EC[1:10,1], type = "p", pch = 16, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(MA[psi](4)))
points(grid_prob, MA_psi_4_results_EC[1:10,2], pch = 17, col = 2)
points(grid_prob, MA_psi_4_results_EC[1:10,3], pch = 18, col = 3)
points(grid_prob, grid_prob, pch = 4, col = 4)
dev.off()

savepdf("MA_psi_4_n300", width = 12, height = 10, toplines = 0.8)
plot(grid_prob, MA_psi_4_results_EC[1:10,4], type = "p", pch = 16, ylim = c(0, 1),
     xlab = "Nominal Coverage Probability", ylab = "Empirical Coverage Probability",
     main = expression(MA[psi](4)))
points(grid_prob, MA_psi_4_results_EC[1:10,5], pch = 17, col = 2)
points(grid_prob, MA_psi_4_results_EC[1:10,6], pch = 18, col = 3)
points(grid_prob, grid_prob, pch = 4, col = 4)
dev.off()
