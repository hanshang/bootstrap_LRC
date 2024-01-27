###################
# sieve bootstrap
# 200 replications
###################

# MA_1_0

registerDoMC(16)
single_MA_1_0_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve", 
                                                               seed_number = iwk, bootnumber = 400, samplesize = 100)

single_MA_1_0_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve", 
                                                               seed_number = iwk, bootnumber = 400, samplesize = 300)

single_MA_1_0_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "single", boot_method = "sieve", 
                                                               seed_number = iwk, bootnumber = 400, samplesize = 500)

double_MA_1_0_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "double", boot_method = "sieve", 
                                                               seed_number = iwk, bootnumber = 400, samplesize = 100)

double_MA_1_0_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_1_0", boot_procedure = "double", boot_method = "sieve", 
                                                               seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_1_0_n100_dist = sieve_MA_1_0_n300_dist = double_MA_1_0_n100_dist = double_MA_1_0_n300_dist = matrix(NA, 200, 10)
sieve_MA_1_0_n100_dist_true = sieve_MA_1_0_n300_dist_true = double_MA_1_0_n100_dist_true = double_MA_1_0_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_MA_1_0_n100_dist[ik,] = single_MA_1_0_n100[[ik]]$dist
    sieve_MA_1_0_n100_dist_true[ik] = single_MA_1_0_n100[[ik]]$disttrue
    
    sieve_MA_1_0_n300_dist[ik,] = single_MA_1_0_n300[[ik]]$dist
    sieve_MA_1_0_n300_dist_true[ik] = single_MA_1_0_n300[[ik]]$disttrue
    
    double_MA_1_0_n100_dist[ik,] = double_MA_1_0_n100[[ik]]$dist
    double_MA_1_0_n100_dist_true[ik] = double_MA_1_0_n100[[ik]]$disttrue
    
    double_MA_1_0_n300_dist[ik,] = double_MA_1_0_n300[[ik]]$dist
    double_MA_1_0_n300_dist_true[ik] = double_MA_1_0_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_MA_1_0_n100_dist) = colnames(sieve_MA_1_0_n300_dist) = 
colnames(double_MA_1_0_n100_dist) = colnames(double_MA_1_0_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_MA_1_0_n100_dist) = rownames(sieve_MA_1_0_n300_dist) = 
rownames(double_MA_1_0_n100_dist) = rownames(double_MA_1_0_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = double_measure_n100 = double_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_MA_1_0_n100_dist_true[ik] <= sieve_MA_1_0_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_MA_1_0_n300_dist_true[ik] <= sieve_MA_1_0_n300_dist[ik,iw], 1, 0)
        
        double_measure_n100[ik,iw] = ifelse(double_MA_1_0_n100_dist_true[ik] <= double_MA_1_0_n100_dist[ik,iw], 1, 0)
        double_measure_n300[ik,iw] = ifelse(double_MA_1_0_n300_dist_true[ik] <= double_MA_1_0_n300_dist[ik,iw], 1, 0)
    }
}
sieve_MA_1_0_n100_EC = colSums(count_measure_n100)/200   # 0.585 0.595 0.620 0.665 0.700 0.720 0.745 0.780 0.800 0.845
sieve_MA_1_0_n300_EC = colSums(count_measure_n300)/200   # 0.610 0.635 0.670 0.715 0.745 0.785 0.825 0.855 0.885 0.915
double_MA_1_0_n100_EC = colSums(double_measure_n100)/200 # 0.520 0.550 0.595 0.620 0.670 0.690 0.735 0.775 0.825 0.865
double_MA_1_0_n300_EC = colSums(double_measure_n300)/200 
rm(count_measure_n100); rm(count_measure_n300); rm(double_measure_n100); rm(double_measure_n300)

grid_prob = seq(0.5, 0.95, by = 0.05)
sieve_MA_1_0_n100_CPD = abs(sieve_MA_1_0_n100_EC - grid_prob)
sieve_MA_1_0_n300_CPD = abs(sieve_MA_1_0_n300_EC - grid_prob)
double_MA_1_0_n100_CPD = abs(double_MA_1_0_n100_EC - grid_prob)
double_MA_1_0_n300_CPD = abs(double_MA_1_0_n300_EC - grid_prob)

mean(sieve_MA_1_0_n100_CPD)  # 0.0525
mean(sieve_MA_1_0_n300_CPD)  # 0.049
mean(double_MA_1_0_n100_CPD) # 0.0445
mean(double_MA_1_0_n300_CPD) # 0.034

# MA_0.5_1

single_MA_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

single_MA_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

single_MA_0.5_1_n500 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 500)

double_MA_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "double",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

double_MA_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_1", boot_procedure = "double",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_0.5_1_n100_dist = sieve_MA_0.5_1_n300_dist = double_MA_0.5_1_n100_dist = double_MA_0.5_1_n300_dist = matrix(NA, 200, 10)
sieve_MA_0.5_1_n100_dist_true = sieve_MA_0.5_1_n300_dist_true = double_MA_0.5_1_n100_dist_true = double_MA_0.5_1_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_MA_0.5_1_n100_dist[ik,] = single_MA_0.5_1_n100[[ik]]$dist
    sieve_MA_0.5_1_n100_dist_true[ik] = single_MA_0.5_1_n100[[ik]]$disttrue
    
    sieve_MA_0.5_1_n300_dist[ik,] = single_MA_0.5_1_n300[[ik]]$dist
    sieve_MA_0.5_1_n300_dist_true[ik] = single_MA_0.5_1_n300[[ik]]$disttrue
    
    double_MA_0.5_1_n100_dist[ik,] = double_MA_0.5_1_n100[[ik]]$dist
    double_MA_0.5_1_n100_dist_true[ik] = double_MA_0.5_1_n100[[ik]]$disttrue
    
    print(ik); rm(ik)
}
colnames(sieve_MA_0.5_1_n100_dist) = colnames(sieve_MA_0.5_1_n300_dist) = colnames(double_MA_0.5_1_n100_dist) = colnames(double_MA_0.5_1_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_MA_0.5_1_n100_dist) = rownames(sieve_MA_0.5_1_n300_dist) = rownames(double_MA_0.5_1_n100_dist) = rownames(double_MA_0.5_1_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = double_measure_n100 = double_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_MA_0.5_1_n100_dist_true[ik] <= sieve_MA_0.5_1_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_MA_0.5_1_n300_dist_true[ik] <= sieve_MA_0.5_1_n300_dist[ik,iw], 1, 0)
        
        double_measure_n100[ik,iw] = ifelse(double_MA_0.5_1_n100_dist_true[ik] <= double_MA_0.5_1_n100_dist[ik,iw], 1, 0)
        double_measure_n300[ik,iw] = ifelse(double_MA_0.5_1_n300_dist_true[ik] <= double_MA_0.5_1_n300_dist[ik,iw], 1, 0)
    }
}
sieve_MA_0.5_1_n100_EC = colSums(count_measure_n100)/200   # 0.470 0.480 0.510 0.555 0.575 0.600 0.640 0.670 0.710 0.760
sieve_MA_0.5_1_n300_EC = colSums(count_measure_n300)/200   # 0.435 0.505 0.560 0.615 0.645 0.675 0.700 0.725 0.770 0.795
double_MA_0.5_1_n100_EC = colSums(double_measure_n100)/200 # 0.335 0.385 0.405 0.455 0.525 0.545 0.580 0.635 0.695 0.760
double_MA_0.5_1_n300_EC = colSums(double_measure_n300)/200 # 
rm(count_measure_n100); rm(count_measure_n300); rm(double_measure_n100); rm(double_measure_n300) 

sieve_MA_0.5_1_n100_CPD = abs(sieve_MA_0.5_1_n100_EC - grid_prob)
sieve_MA_0.5_1_n300_CPD = abs(sieve_MA_0.5_1_n300_EC - grid_prob)
double_MA_0.5_1_n100_CPD = abs(double_MA_0.5_1_n100_EC - grid_prob)


# MA_0.5_4

single_MA_0.5_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

single_MA_0.5_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

double_MA_0.5_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "double", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

double_MA_0.5_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_4", boot_procedure = "double", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_0.5_4_n100_dist = sieve_MA_0.5_4_n300_dist = double_MA_0.5_4_n100_dist = double_MA_0.5_4_n300_dist = matrix(NA, 200, 10)
sieve_MA_0.5_4_n100_dist_true = sieve_MA_0.5_4_n300_dist_true = double_MA_0.5_4_n100_dist_true = double_MA_0.5_4_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_MA_0.5_4_n100_dist[ik,] = single_MA_0.5_4_n100[[ik]]$dist
    sieve_MA_0.5_4_n100_dist_true[ik] = single_MA_0.5_4_n100[[ik]]$disttrue
    
    sieve_MA_0.5_4_n300_dist[ik,] = single_MA_0.5_4_n300[[ik]]$dist
    sieve_MA_0.5_4_n300_dist_true[ik] = single_MA_0.5_4_n300[[ik]]$disttrue
    
    double_MA_0.5_4_n100_dist[ik,] = double_MA_0.5_4_n100[[ik]]$dist
    double_MA_0.5_4_n100_dist_true[ik] = double_MA_0.5_4_n100[[ik]]$disttrue
    
    double_MA_0.5_4_n300_dist[ik,] = double_MA_0.5_4_n300[[ik]]$dist
    double_MA_0.5_4_n300_dist_true[ik] = double_MA_0.5_4_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_MA_0.5_4_n100_dist) = colnames(sieve_MA_0.5_4_n300_dist) = colnames(double_MA_0.5_4_n100_dist) = colnames(double_MA_0.5_4_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_MA_0.5_4_n100_dist) = rownames(sieve_MA_0.5_4_n300_dist) = rownames(double_MA_0.5_4_n100_dist) = rownames(double_MA_0.5_4_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = double_measure_n100 = double_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_MA_0.5_4_n100_dist_true[ik] <= sieve_MA_0.5_4_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_MA_0.5_4_n300_dist_true[ik] <= sieve_MA_0.5_4_n300_dist[ik,iw], 1, 0)
        
        double_measure_n100[ik,iw] = ifelse(double_MA_0.5_4_n100_dist_true[ik] <= double_MA_0.5_4_n100_dist[ik,iw], 1, 0)
        double_measure_n300[ik,iw] = ifelse(double_MA_0.5_4_n300_dist_true[ik] <= double_MA_0.5_4_n300_dist[ik,iw], 1, 0)
    }
}
sieve_MA_0.5_4_n100_EC = colSums(count_measure_n100)/200 # 0.215 0.250 0.275 0.290 0.320 0.330 0.345 0.370 0.425 0.495
sieve_MA_0.5_4_n300_EC = colSums(count_measure_n300)/200 # 0.245 0.280 0.310 0.340 0.385 0.420 0.460 0.505 0.565 0.635

double_MA_0.5_4_n100_EC = colSums(double_measure_n100)/200 # 0.110 0.150 0.175 0.195 0.245 0.270 0.325 0.355 0.420 0.485
double_MA_0.5_4_n300_EC = colSums(double_measure_n300)/200 # 0.180 0.210 0.250 0.260 0.310 0.360 0.400 0.460 0.535 0.650
rm(count_measure_n100); rm(count_measure_n300); rm(double_measure_n100); rm(double_measure_n300)

grid_prob = seq(0.5, 0.95, by = 0.05)
sieve_MA_0.5_4_n100_CPD = abs(sieve_MA_0.5_4_n100_EC - grid_prob)
sieve_MA_0.5_4_n300_CPD = abs(sieve_MA_0.5_4_n300_EC - grid_prob)
double_MA_0.5_4_n100_CPD = abs(double_MA_0.5_4_n100_EC - grid_prob)
double_MA_0.5_4_n300_CPD = abs(double_MA_0.5_4_n300_EC - grid_prob)

mean(sieve_MA_0.5_4_n100_CPD)  # 0.3935
mean(sieve_MA_0.5_4_n300_CPD)  # 0.3105
mean(double_MA_0.5_4_n100_CPD) # 0.452
mean(double_MA_0.5_4_n300_CPD) # 0.3635

# MA_0.5_8

theta = 0.5
single_MA_0.5_8_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)
single_MA_0.5_8_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

double_MA_0.5_8_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "double",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)
double_MA_0.5_8_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_0.5_8", boot_procedure = "double",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)


sieve_MA_0.5_8_n100_dist = sieve_MA_0.5_8_n300_dist = matrix(NA, 200, 10)
sieve_MA_0.5_8_n100_dist_true = sieve_MA_0.5_8_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_MA_0.5_8_n100_dist[ik,] = single_MA_0.5_8_n100[[ik]]$dist
    sieve_MA_0.5_8_n100_dist_true[ik] = single_MA_0.5_8_n100[[ik]]$disttrue
    
    sieve_MA_0.5_8_n300_dist[ik,] = single_MA_0.5_8_n300[[ik]]$dist
    sieve_MA_0.5_8_n300_dist_true[ik] = single_MA_0.5_8_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_MA_0.5_8_n100_dist) = colnames(sieve_MA_0.5_8_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_MA_0.5_8_n100_dist) = rownames(sieve_MA_0.5_8_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_MA_0.5_8_n100_dist_true[ik] < sieve_MA_0.5_8_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_MA_0.5_8_n300_dist_true[ik] < sieve_MA_0.5_8_n300_dist[ik,iw], 1, 0)
    }
}
sieve_MA_0.5_8_n100_EC = colSums(count_measure_n100)/200 # 0.055 0.055 0.070 0.080 0.080 0.090 0.100 0.120 0.160 0.195
sieve_MA_0.5_8_n300_EC = colSums(count_measure_n300)/200 # 0.075 0.090 0.120 0.130 0.150 0.170 0.180 0.215 0.240 0.290
rm(count_measure_n100); rm(count_measure_n300)

# FAR_0.5_1

single_AR_0.5_1_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)
single_AR_0.5_1_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.5_1", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_AR_0.5_1_n100_dist = sieve_AR_0.5_1_n300_dist = matrix(NA, 200, 10)
sieve_AR_0.5_1_n100_dist_true = sieve_AR_0.5_1_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_AR_0.5_1_n100_dist[ik,] = single_AR_0.5_1_n100[[ik]]$dist
    sieve_AR_0.5_1_n100_dist_true[ik] = single_AR_0.5_1_n100[[ik]]$disttrue
    
    sieve_AR_0.5_1_n300_dist[ik,] = single_AR_0.5_1_n300[[ik]]$dist
    sieve_AR_0.5_1_n300_dist_true[ik] = single_AR_0.5_1_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_AR_0.5_1_n100_dist) = colnames(sieve_AR_0.5_1_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_AR_0.5_1_n100_dist) = rownames(sieve_AR_0.5_1_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_AR_0.5_1_n100_dist_true[ik] < sieve_AR_0.5_1_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_AR_0.5_1_n300_dist_true[ik] < sieve_AR_0.5_1_n300_dist[ik,iw], 1, 0)
    }
}
sieve_AR_0.5_1_n100_EC = colSums(count_measure_n100)/200  # 0.235 0.255 0.285 0.320 0.360 0.395 0.435 0.475 0.525 0.575
sieve_AR_0.5_1_n300_EC = colSums(count_measure_n300)/200  # 0.315 0.370 0.395 0.425 0.450 0.455 0.505 0.540 0.585 0.675

# FAR_0.3_2

registerDoMC(15)
single_AR_0.3_2_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)

single_AR_0.3_2_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "FAR_0.3_2", boot_procedure = "single",
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_AR_0.3_2_n100_dist = sieve_AR_0.3_2_n300_dist = matrix(NA, 200, 10)
sieve_AR_0.3_2_n100_dist_true = sieve_AR_0.3_2_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_AR_0.3_2_n100_dist[ik,] = single_AR_0.3_2_n100[[ik]]$dist
    sieve_AR_0.3_2_n100_dist_true[ik] = single_AR_0.3_2_n100[[ik]]$disttrue
    
    sieve_AR_0.3_2_n300_dist[ik,] = single_AR_0.3_2_n300[[ik]]$dist
    sieve_AR_0.3_2_n300_dist_true[ik] = single_AR_0.3_2_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_AR_0.3_2_n100_dist) = colnames(sieve_AR_0.3_2_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_AR_0.3_2_n100_dist) = rownames(sieve_AR_0.3_2_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_AR_0.3_2_n100_dist_true[ik] <= sieve_AR_0.3_2_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_AR_0.3_2_n300_dist_true[ik] <= sieve_AR_0.3_2_n300_dist[ik,iw], 1, 0)
    }
}
sieve_AR_0.3_2_n100_EC = colSums(count_measure_n100)/200 # 0.305 0.340 0.355 0.410 0.445 0.490 0.510 0.550 0.580 0.615
sieve_AR_0.3_2_n300_EC = colSums(count_measure_n300)/200 # 0.365 0.420 0.445 0.460 0.495 0.525 0.535 0.595 0.655 0.740


# MA_psi_4

single_MA_psi_4_n100 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 100)
single_MA_psi_4_n300 = foreach(iwk = 1:200) %dopar% bootrep_eval(DGP_process = "MA_psi_4", boot_procedure = "single", 
                                                                 seed_number = iwk, bootnumber = 400, samplesize = 300)

sieve_MA_psi_4_n100_dist = sieve_MA_psi_4_n300_dist = matrix(NA, 200, 10)
sieve_MA_psi_4_n100_dist_true = sieve_MA_psi_4_n300_dist_true = vector("numeric", 200)
for(ik in 1:200)
{
    sieve_MA_psi_4_n100_dist[ik,] = single_MA_psi_4_n100[[ik]]$dist
    sieve_MA_psi_4_n100_dist_true[ik] = single_MA_psi_4_n100[[ik]]$disttrue
    
    sieve_MA_psi_4_n300_dist[ik,] = single_MA_psi_4_n300[[ik]]$dist
    sieve_MA_psi_4_n300_dist_true[ik] = single_MA_psi_4_n300[[ik]]$disttrue
    print(ik); rm(ik)
}
colnames(sieve_MA_psi_4_n100_dist) = colnames(sieve_MA_psi_4_n300_dist) = seq(0.5, 0.95, by = 0.05)
rownames(sieve_MA_psi_4_n100_dist) = rownames(sieve_MA_psi_4_n300_dist) = 1:200

# empirical coverage

count_measure_n100 = count_measure_n300 = matrix(NA, 200, 10)
for(ik in 1:200)
{
    for(iw in 1:10)
    {
        count_measure_n100[ik,iw] = ifelse(sieve_MA_psi_4_n100_dist_true[ik] <= sieve_MA_psi_4_n100_dist[ik,iw], 1, 0)
        count_measure_n300[ik,iw] = ifelse(sieve_MA_psi_4_n300_dist_true[ik] <= sieve_MA_psi_4_n300_dist[ik,iw], 1, 0)
    }
}
sieve_MA_psi_4_n100_EC = colSums(count_measure_n100)/200 # 0.015 0.025 0.025 0.045 0.070 0.070 0.090 0.120 0.155 0.320
sieve_MA_psi_4_n300_EC = colSums(count_measure_n300)/200 # something is wrong

