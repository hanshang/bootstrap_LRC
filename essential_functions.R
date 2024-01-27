##################
# load R packages
##################

packages <- c("sde", "far", "vars", "ftsa", "doMC", "meboot", "sandwich", "matrixcalc")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))
options(scipen=999)

# simulating Brownian motion

BM <- function (x = 0, t0 = 0, T = 1, N = 100)
{
  if (T <= t0)
    stop("wrong times")
  dt <- (T - t0)/N
  t <- seq(t0, T, length = N + 1)
  X <- ts(cumsum(c(x, rnorm(N) * sqrt(dt))), start = t0, deltat = dt)
  return(invisible(X))
}

# MA_1_0

BrownMat <- function(N,refinement)
{
  mat <- matrix(nrow=refinement,ncol=N)
  c <- 1
  while(c <= N)
  {
    vec <- BM(N=refinement-1)
    mat[,c] <- vec
    c <- c+1
  }
  return(mat)
}

# MA_0.5_1

BMlagMat <- function(N,refinement)
{
  firstmat <- BrownMat(N+1,refinement)
  mat <- matrix(nrow=refinement,ncol=N)
  c <- 1
  while(c <= N)
  {
    vec <- firstmat[,c]/2 + firstmat[,c+1]/2
    mat[,c] <- vec
    c <- c+1
  }
  return(mat)
}

# MA_0.5_4

MA5Mat <- function(N,refinement)
{
  firstmat <- BrownMat(N+4,refinement)
  mat <- matrix(nrow=refinement,ncol=N)
  c <- 1
  while(c <= N)
  {
    vec <- theta*firstmat[,c] + theta*firstmat[,c+1] + theta*firstmat[,c+2] + theta*firstmat[,c+3] + theta*firstmat[,c+4]
    mat[,c] <- vec
    c <- c+1
  }
  return(mat)
}

# MA_phi_4

MAphi = function(ref, k)
{
  Mat = matrix(nrow=ref, ncol=ref)
  for(i in 1:ref)
  {
    for(j in 1:ref)
    {
      Mat[i,j] = min(i,j)/ref
    }
  }
  return(k * Mat)
}

funMAMat = function(N, refinement)
{
  Mat = MAphi(refinement, 1.5)
  firstMat = BrownMat(N+4, refinement)
  finalMat = matrix(nrow=refinement, ncol=N)
  for(c in 1:N)
  {
    finalMat[,c] = funIntegral(refinement, Mat, firstMat[,c]) + funIntegral(refinement, Mat, firstMat[,(c+1)]) + funIntegral(refinement, Mat, firstMat[,(c+2)]) + funIntegral(refinement, Mat, firstMat[,(c+3)]) + firstMat[,(c+4)]
  }
  return(finalMat)
}

# MA_0.5_8

MA9Mat <- function(N,refinement)
{
  firstmat <- BrownMat(N+8,refinement)
  mat <- matrix(nrow=refinement,ncol=N)
  c <- 1
  while(c <= N)
  {
    vec <- theta*firstmat[,c] + theta*firstmat[,c+1] + theta*firstmat[,c+2] + theta*firstmat[,c+3] + theta*firstmat[,c+4] + theta*firstmat[,c+5] + theta*firstmat[,c+6] + theta*firstmat[,c+7] + theta*firstmat[,c+8]
    mat[,c] <- vec
    c <- c+1
  }
  return(mat)
}

# FAR_0.5_1

ZlagMat.5 <- function(N,refinement)
{
  firstmat <- BrownMat(N+50, refinement)
  mat <- matrix(nrow=refinement,ncol=N+50)
  mat[,1] = firstmat[,1]
  c <- 2
  M <- N+50
  while(c <= M)
  {
    vec <- firstmat[,c] + 0.5*mat[,c-1]
    mat[,c] <- vec
    c <- c+1
  }
  finalmat <- mat[1:refinement,51:M]
  return(finalmat)
}

# FAR_0.3_2

ZlagMat.3_2 <- function(N, refinement)
{
  firstmat <- BrownMat(N+50, refinement)
  mat <- matrix(nrow = refinement, ncol = N+50)
  mat[,1:2] <- firstmat[,1:2]
  c <- 3
  M <- N + 50
  while(c <= M)
  {
    mat[,c] <- firstmat[,c] + 0.6 * mat[,c-1] - 0.09 * mat[,c-2]
    c <- c+1
  }
  finalmat <- mat[1:refinement, 51:M]
  return(finalmat)
}

# FAR_psi_1

funKernel = function(ref, k)
{
  Mat = matrix(nrow=ref, ncol=ref)
  for(i in 1:ref)
  {
    for(j in 1:ref)
    {
      Mat[i,j] = k*exp(.5*((i/ref)^2+(j/ref)^2))
    }
  }
  return(Mat)
}

funIntegral = function(ref, Mat, X)
{
  Mat = Mat%*%X
  return(Mat/ref)
}

funARMat = function(N, refinement)
{
  Mat = funKernel(refinement, 0.34)
  firstmat <- matrix(nrow=refinement,ncol=N+50)
  firstmat[,1] = BM(N=refinement-1)
  for(i in 2:(N+50))
  {
    firstmat[,i] = (funIntegral(refinement, Mat, firstmat[,(i-1)]) + BM(N=refinement-1))
  }
  finalMat = firstmat[1:refinement,51:(N+50)]
  return(finalMat)
}

BMCov <- function(refinement)
{
  mat <- matrix(rep(0, refinement^2), ncol=refinement)
  t = s = 1:refinement/refinement
  c1 <- 1
  while(c1 <= refinement)
  {
    c2 <- 1
    while(c2 <= refinement)
    {
      mat[c1,c2] <- min(t[c1],s[c2])
      c2 <- c2+1
    }
    c1 <- c1+1
  }
  return(mat)
}

#################################
# seed: random seed
# N: sample size
# refinement: grid size
# dgp: choose among various DGPs
#################################

dgp_fun <- function(seed, N, refinement, dgp)
{
  set.seed(122 + seed)
  if(dgp == "MA_1_0")
  {
    data1 = BrownMat(N, refinement)
  }
  else if(dgp == "MA_0.5_1")
  {
    data1 = BMlagMat(N, refinement)
  }
  else if(dgp == "MA_0.5_4")
  {
    data1 = MA5Mat(N, refinement)
  }
  else if(dgp == "MA_psi_4")
  {
    data1 = funMAMat(N, refinement)
  }
  else if(dgp == "MA_0.5_8")
  {
    data1 = MA9Mat(N, refinement)
  }
  else if(dgp == "FAR_0.5_1")
  {
    data1 = ZlagMat.5(N, refinement)
  }
  else if(dgp == "FAR_0.3_2")
  {
    data1 = ZlagMat.3_2(N, refinement)
  }
  else if(dgp == "FAR_psi_1")
  {
    data1 = funARMat(N, refinement)
  }
  return(data1)
}

####################
# IID bootstrapping
####################

stscorebootstrapdata = function(dat, bootrep)
{
    dat = t(dat)
    n = nrow(dat)
    p = ncol(dat)
    boot = matrix(NA, p, bootrep)
    bootdataarray = array(NA, dim = c(p, n, bootrep))

    data = scale(dat, center = TRUE, scale = FALSE)
    meandata = matrix(rep(as.matrix(colMeans(dat)), n), p, n)

    dummy = svd(data)
    U = dummy$u
    D = diag(dummy$d)
    V = dummy$v
    for(b in 1:bootrep)
    {
      drawU = U[sample(1:n, n, replace = TRUE), ]
      bootdata = t(drawU%*%D%*%t(V)) + meandata
      bootdataarray[,,b] = bootdata
    }
    return(bootdataarray)
}

####################
# FAR bootstrapping
####################

# sim_dat: p by n
# bootrep: number of bootstraps

far_boot <- function(sim_dat, bootrep)
{
    dat = t(sim_dat)
    n = dim(dat)[1]
    p = dim(dat)[2]
    sim_fdata = as.fdata(as.numeric(sim_dat[,1:(n-1)]), col = 1, p = p, dates = 1:(n-1), name = "X")
    far_obj = far(data=sim_fdata, y="X", center=TRUE, na.rm = FALSE, joined = TRUE)

    pred1 = predict(far_obj, newdata=sim_fdata)
    far_resi = sim_dat[,2:n] - pred1[[1]]

    # stscorebootstrap decenters the data
    far_resi_boot = stscorebootstrapdata(dat = far_resi, bootrep = bootrep)

    far_boot = array(NA, dim = c(p, n, bootrep))
    boot = matrix(NA, p, bootrep)
    for(b in 1:bootrep)
    {
        far_boot[,,b] = cbind(sim_dat[,1], pred1[[1]] + far_resi_boot[,,b])
    }
    return(list(X_boot = far_boot))
}

############################
# FAR + sieve bootstrapping
############################

far_sieve <- function(sim_dat, bootrep)
{
    dat = t(sim_dat)
    n = dim(dat)[1]
    p = dim(dat)[2]
    sim_fdata = as.fdata(as.numeric(sim_dat[,1:(n-1)]), col = 1, p = p, dates = 1:(n-1), name = "X")
    far_obj = far(data=sim_fdata, y="X", center=TRUE, na.rm = FALSE, joined = TRUE)

    pred1 = predict(far_obj, newdata=sim_fdata)
    far_resi = sim_dat[,2:n] - pred1[[1]]

    # sieve bootstrap decenters the data
    far_resi_boot = sieve_bootstrap(fun_dat = t(far_resi), ncomp_porder_selection = "CPV_AICC", B = bootrep,
                                    burn_in = 50)$X_boot

    far_boot = array(NA, dim = c(p, n, bootrep))
    boot = matrix(NA, p, bootrep)
    for(b in 1:bootrep)
    {
        far_boot[,,b] = cbind(sim_dat[,1], pred1[[1]] + far_resi_boot[,,b])
    }
    return(list(X_boot = far_boot))
}


#########################################
# compute population long-run covariance
#########################################

# DGP: data generating process
# seed: random seed
# T: sample size
# no_grid: number of grids
# theta: moderate dependence

lrc_pop <- function(DGP = c("MA_1_0", "MA_0.5_1", "MA_0.5_4", "MA_psi_4", "MA_0.5_8", "FAR_0.5_1", "FAR_0.3_2", "FAR_psi_1"),
                    seed, T, no_grid = 101, theta = 0.5)
{
    ###############################################
    # calculate asympototic (theoretical) variance
    ###############################################

    set.seed(122 + seed)
    if(DGP == "MA_1_0")
    {
      realCov <- BMCov(no_grid)
    }
    else if(DGP == "MA_0.5_1")
    {
      realCov <- BMCov(no_grid)
    }
    else if(DGP == "MA_0.5_4")
    {
      realCov <- (theta^2)*25*BMCov(no_grid)
    }
    # theoretical is not available, so by simulation
    else if(DGP == "MA_psi_4")
    {
      dat = dgp_fun(seed = 0, N = 100000, refinement = 101, dgp = "MA_psi_4")
      center_dat_standardised_mean = as.matrix((ncol(dat)^(0.5) * rowMeans(dat)))
      decenter_dat = sweep(dat, 1, center_dat_standardised_mean, "-")
      realCov = cov(t(decenter_dat))
      rm(dat); rm(center_dat_standardised_mean); rm(decenter_dat)
    }
    else if(DGP == "MA_0.5_8")
    {
      realCov <- (theta^2)*81*BMCov(no_grid)
    }
    else if(DGP == "FAR_0.5_1")
    {
      realCov <- BMCov(no_grid)*4
    }
    else if(DGP == "FAR_0.3_2")
    {
      realCov <- BMCov(no_grid) * ((1/0.49)^2)
    }
    else if(DGP == "FAR_psi_1")
    {
      dat = dgp_fun(seed = 0, N = 100000, refinement = 101, dgp = "FAR_psi_1")
      center_dat_standardised_mean = as.matrix((ncol(dat)^(0.5) * rowMeans(dat)))
      decenter_dat = sweep(dat, 1, center_dat_standardised_mean, "-")
      realCov = cov(t(decenter_dat))
      rm(dat); rm(center_dat_standardised_mean); rm(decenter_dat)
    }
    else
    {
      warning("DGP is not on the list.")
    }

    dat_val = dgp_fun(seed = seed, N = T, refinement = no_grid, dgp = DGP)
    sampleCov = long_run_covariance_estimation(dat = dat_val)

    # fixed bandwidth and finite-order kernel

    err_FT_BT_proposed = hilbert.schmidt.norm(sampleCov - realCov)
    return(list(sample_data = dat_val, sampleCov = sampleCov, realCov = realCov, err = err_FT_BT_proposed))
}

#####################################
# compute sample long-run covariance
#####################################

# sample_dat: p by n matrix

lrc_sample <- function(sample_dat, sampleCov, no_boot, boot_method = c("sieve", "FAR", "FAR_sieve"))
{
    no_grid = nrow(sample_dat)
    # T = ncol(sample_dat)

    if(boot_method == "sieve")
    {
        dat = sieve_bootstrap(fun_dat = t(sample_dat), ncomp_porder_selection = "eigenratio_AICC",
                              B = no_boot, burn_in = 50)
    }
    else if(boot_method == "FAR")
    {
        dat = far_boot(sim_dat = sample_dat, bootrep = no_boot)
    }
    else if(boot_method == "FAR_sieve")
    {
        dat = far_sieve(sim_dat = sample_dat, bootrep = no_boot)
    }
    else
    {
        warning("boot_method is not on the list.")
    }

    err_FT_BT_boot = vector("numeric", no_boot)
    # bootCov_array = array(NA, dim = c(no_grid, no_grid, no_boot))
    for(iw in 1:no_boot)
    {
        bootCov = long_run_covariance_estimation(dat = dat$X_boot[,,iw])
        err_FT_BT_boot[iw] = hilbert.schmidt.norm(bootCov - sampleCov)
        rm(bootCov); rm(iw)
    }
    return(list(sample_boot_err = err_FT_BT_boot))
}


#######################
# bootstrap evaluation
#######################

bootrep_eval <- function(DGP_process, boot_procedure, boot_method, seed_number, bootnumber, samplesize,
                         cover = seq(0.5, 0.95, by = 0.05))
{
    pop_statistic = lrc_pop(DGP = DGP_process, seed = seed_number, T = samplesize)
    no_grid = nrow(pop_statistic$sample_data)

    sample_statistic = lrc_sample(sample_dat = pop_statistic$sample_data, sampleCov = pop_statistic$sampleCov,
                                  no_boot = bootnumber, boot_method = boot_method)

    if(boot_procedure == "single")
    {
        dist = vector("numeric", length(cover))
        for(ik in 1:length(cover))
        {
            dist[ik] = sort(sample_statistic$sample_boot_err)[cover[ik] * bootnumber]
            rm(ik)
        }
    }
    else if(boot_procedure == "double")
    {
        double_boot_sample = array(NA, dim = c(nrow(pop_statistic$sample_data), ncol(pop_statistic$sample_data), bootnumber))
        for(iw in 1:bootnumber)
        {
            double_boot_sample[,,iw] = sieve_bootstrap(fun_dat = t(sample_statistic$single_boot_sample[,,iw]),
                                                       ncomp_porder_selection = "eigenratio_AICC", B = 1, burn_in = 50)$X_boot
            rm(iw)
        }
        err_FT_BT_double_boot = vector("numeric", bootnumber)
        for(iw in 1:bootnumber)
        {
            bootCov = long_run_covariance_estimation(dat = double_boot_sample[,,iw])
            diff_C_FT_BT = bootCov - sample_statistic$single_bootCov[,,iw]
            err_FT_BT_double_boot[iw] = sum(diff_C_FT_BT^2)/(no_grid^2)
            rm(bootCov); rm(diff_C_FT_BT); rm(iw)
        }
        sample_double_boot_err = err_FT_BT_double_boot
        rm(err_FT_BT_double_boot)

        dist = vector("numeric", length(cover))
        for(ik in 1:length(cover))
        {
            dist[ik] = sort(sample_double_boot_err)[cover[ik] * bootnumber]
            rm(ik)
        }
    }
    else
    {
        warning("Bootstrap procedure is neither single nor double.")
    }
    rm(sample_statistic)
    return(list(dist = dist, disttrue = pop_statistic$err))
}
