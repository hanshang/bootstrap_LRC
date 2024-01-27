##################
# Sieve bootstrap
##################

# fun_dat: data matrix (n by p)
# grid: grid points in a curve
# p_m_order: give a pre-specified number of components and VAR order
# ncomp_porder_selection: CPV_AICC is the default choice
# CPV_percent: 0.85 (by default), for robustness check, we also consider 0.8 and 0.9
# VAR_type: default is none

sieve_bootstrap <- function(fun_dat, grid, p_m_order, ncomp_porder_selection,
                            CPV_percent = 0.85, VAR_type = "none", B, burn_in)
{
    #########
    # Step 1
    #########

    n_row = nrow(fun_dat)
    n_col = ncol(fun_dat)

    # de-centering functional data

    center_fun_dat = scale(fun_dat, center = TRUE, scale = FALSE)

    # computing the mean of functional data

    X_bar = colMeans(fun_dat)

    # determine p_order

    if(missing(grid))
    {
      if(is.null(colnames(fun_dat)))
      {
        grid = seq(1, n_col, length.out = n_col)
      }
      else
      {
        grid = seq(min(as.numeric(colnames(fun_dat))),
                   max(as.numeric(colnames(fun_dat))),
                   length.out = n_col)
      }
    }

    # determine m_component based on FPE of Aue et al. (2015)

    colnames(center_fun_dat) = grid
    rownames(center_fun_dat) = 1:n_row
    fun_dat_object = fts(grid, t(center_fun_dat))

    if(missing(p_m_order))
    {
        if(ncomp_porder_selection == "FPE")
        {
            ncomp_porder = ftsa:::method.FPE(object = fun_dat_object, D = 6, var_type = VAR_type, Pmax = 5)
        }
        else
        {
            if(ncomp_porder_selection == "CPV_AICC")
            {
                ncomp = head(which(cumsum(ftsm(y = fun_dat_object, order = min(nrow(center_fun_dat), ncol(center_fun_dat)))$varprop)>= CPV_percent), 1)
            }
            else if(ncomp_porder_selection == "eigenratio_AICC")
            {
                NIR_long_run_cov = long_run_covariance_estimation(dat = t(fun_dat))
                NIR_eigen_decomp = eigen(NIR_long_run_cov, symmetric = TRUE)

                lambda_val = NIR_eigen_decomp$values
                k_max = length(which(lambda_val >= mean(lambda_val)))

                tau = 1/log(max(lambda_val[1], n_row))
                eigen_val_ratio = vector("numeric", k_max)
                for(ik in 1:k_max)
                {
                    eigen_val_ratio[ik] = lambda_val[ik+1]/lambda_val[ik] * ifelse(lambda_val[ik]/lambda_val[1] >= tau, 1, 0) + ifelse(lambda_val[ik]/lambda_val[1] < tau, 1, 0)
                }
                ncomp = which.min(eigen_val_ratio)
                rm(NIR_long_run_cov); rm(NIR_eigen_decomp); rm(lambda_val); rm(k_max); rm(tau); rm(eigen_val_ratio)
            }
            else
            {
                warning("ncomp_porder_selection method is incorrectly specified.")
            }
            ftsm_object_scores = ftsm(y = fun_dat_object, order = ncomp, mean = FALSE)$coeff

            VAR_max_order = min(9, nrow(ftsm_object_scores) - 2)
            VAR_AICC = vector("numeric", VAR_max_order)
            if(ncomp == 1)
            {
                for(VAR_order in 1:VAR_max_order)
                {
                    ftsm_object_resi = na.omit(residuals(ar(x = ftsm_object_scores, aic = FALSE, order.max = VAR_order)))
                    sigma = det(crossprod(ftsm_object_resi)/n_row)
                    VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
                    rm(sigma); rm(ftsm_object_resi)
                }
            }
            else
            {
                for(VAR_order in 1:VAR_max_order)
                {
                    ftsm_object_resi = residuals(VAR(y = ftsm_object_scores, p = VAR_order, type = "none"))
                    sigma = det(crossprod(ftsm_object_resi)/n_row)
                    VAR_AICC[VAR_order] = n_row * log(sigma) + n_row * (n_row * ncomp + VAR_order * (ncomp^2))/(n_row - ncomp * (VAR_order + 1) - 1)
                    rm(sigma); rm(ftsm_object_resi)
                }
            }
            ncomp_porder = c(which.min(VAR_AICC[is.finite(VAR_AICC)][VAR_AICC[is.finite(VAR_AICC)]>0]), ncomp)
            rm(ncomp); rm(VAR_AICC); rm(ftsm_object_scores)
        }
    }
    else
    {
        ncomp_porder = p_m_order
    }

    # think about removing p_order = 0

    p_order = ncomp_porder[1]
    m_component = ncomp_porder[2]
    rm(ncomp_porder)

    # compute the matrix variance (X %*% t(X))

    variance_est = var(center_fun_dat)

    # compute eigen-decomposition

    eigen_decomp = eigen(variance_est)
    eigen_decomp_vector = as.matrix(eigen_decomp$vectors[,1:m_component])
    eigen_decomp_score  = center_fun_dat %*% eigen_decomp_vector
    colnames(eigen_decomp_score) = 1:ncol(eigen_decomp_score)

    # reconstruction of original curves (PC score %*% PCs + mean)

    reconstruction = eigen_decomp_score %*% t(eigen_decomp_vector[,1:m_component])
    reconstruction_mean = (sweep(reconstruction, 2, X_bar, FUN = "+"))
    reconstruction_err = fun_dat - reconstruction_mean

    # determine p_order

    rm(variance_est); rm(eigen_decomp); rm(reconstruction); rm(reconstruction_mean)

    #########
    # Step 2
    #########

    if(m_component == 1)
    {
      A_est = matrix(ar(eigen_decomp_score, aic = FALSE, order = p_order)$ar, p_order, 1)
      residual_e = matrix(ar(eigen_decomp_score, aic = FALSE, order = p_order)$resid[-(1:p_order)],,1)
    }
    if(m_component > 1)
    {
      VAR_forward = VAR(eigen_decomp_score, p = p_order, type = VAR_type)
      A_est = matrix(NA, p_order * m_component, m_component)
      for(iw in 1:m_component)
      {
        A_est[,iw] = VAR_forward$varresult[[iw]]$coef
      }
      residual_e = resid(VAR_forward)
      rownames(residual_e) = (p_order + 1):n_row
      rm(VAR_forward)
    }
    colnames(A_est) = colnames(eigen_decomp_score) = colnames(residual_e) = paste("score", 1:m_component, sep="_")
    A_est_list = list()
    for(j in 1:p_order)
    {
      A_est_list[[j]] = A_est[((j-1)*m_component+1):(j*m_component),]
    }

    # eigen_decomp_score[141,] %*% A_est_list[[2]] + eigen_decomp_score[142,] %*% A_est_list[[1]] + eigen_decomp_score[35,] %*% A_est_list[[1]](checked)

    # obtain centered residuals later for bootstrapping

    residual_e_centered = scale(residual_e, center = TRUE, scale = FALSE)

    eigen_decomp_boot = array(NA, dim = c(n_row + burn_in, m_component, B))
    for(iwk in 1:B)
    {
      eigen_decomp_boot[1:p_order,,iwk] = eigen_decomp_score[1:p_order,]
      for(iw in p_order:(n_row + (burn_in - 1)))
      {
        fitted_val = rep(0, m_component)
        for(ij in 1:p_order)
        {
          fitted_val = fitted_val + eigen_decomp_boot[rev(1:iw)[ij],,iwk] %*% A_est_list[[ij]]
        }
        eigen_decomp_boot[(iw+1),,iwk] = fitted_val + residual_e_centered[sample(1:nrow(residual_e_centered), size = 1),]
      }
    }

    # bootstrap samples (used to produce \widehat{X}_star and \widehat{X})

    U_centered = scale(reconstruction_err, center = TRUE, scale = FALSE)
    X_boot = array(NA, dim = c(n_col, (n_row + burn_in), B))
    for(ik in 1:B)
    {
      X_boot[,,ik] = (eigen_decomp_vector %*% t(eigen_decomp_boot[,,ik])) + t(U_centered[sample(1:n_row, size = (n_row + burn_in), replace = TRUE),]) + matrix(rep(X_bar, (n_row + burn_in)), n_col, (n_row + burn_in))
    }
    X_boot_sample = X_boot[,(burn_in+1):(n_row+burn_in),]
    return(list(X_boot = X_boot_sample, p_order = p_order, m_component = m_component))
}
