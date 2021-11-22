

bym2_logistic_allStrata_ixns_test_bias_simple_r1 <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf) 
  sigma_zip_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  # prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip) }
  for(i in 1:N_eth){ a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) }
  for(i in 1:2){ 
    a_inc[i] ~ dnorm(mean = 0, sd = sigma_inc)
    a_edu[i] ~ dnorm(mean = 0, sd = sigma_edu) 
  }
  for(i in 1:6){ a_hh[i] ~ dnorm(mean = 0, sd = sigma_hh) }
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(31*2)){ a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) }
  for(i in 1:(6*2)){ a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) }
  for(i in 1:(6*2)){ a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[1:N])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[1:N])))
  
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[1:N]) - b[3]*mean(x_zip[1:N]), scale = 1)
  
  # bym2
  for (i in 1:N_zip) {
    vc[i] ~ dnorm(0,1)
    Corr[i] <- sigma * uc[i] * sqrt(rho/scale)
    UCorr[i] <- sigma * vc[i] * sqrt((1-rho)) 
    bym2[i] <- Corr[i] + UCorr[i] 
  }
  
  uc[1:N_zip] ~ dcar_normal(adj=adj[1:L], weights=wei[1:L], num=num[1:N_zip], tau=1, zero_mean=1)  
  sigma ~ dunif(0,10)
  rho ~ dbeta(1,1)
  
  # test-bias --------
  # ortho
  y_sens_r1 ~ dbinom(prob = ilogit(logit_sens_r1), size = n_sens_r1)
  y_spec_r1 ~ dbinom(prob = ilogit(logit_spec_r1), size = n_spec_r1)
  
  logit_sens_r1 ~ dnorm(1.45, 0.5) 
  logit_spec_r1 ~ dnorm(5, 2.25) 
  
  # net se/sp 
  se_net_r1 <- ilogit(logit_sens_r1) 
  sp_net_r1 <- ilogit(logit_spec_r1) 
  
  for(i in 1:N){
    
    # likelihood
    p[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
                     a_age[age[i]] + a_eth[eth[i]] + a_hh[hh[i]] + 
                     a_inc[inc[i]+1] + 
                     a_edu[edu[i]] + 
                     bym2[zip[i]] + 
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                     a_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                     a_eth_inc[index_eth_inc[eth[i],round(inc[i])+1]] +
                     a_eth_age[index_eth_age[eth[i],age[i]]]
    )
    y[i] ~ dbern(prob = p_sample[i])
    p_sample[i] <- (p[i] * se_net_r1) + (1 - p[i]) * (1 - sp_net_r1) # prevalence corrected for imperfect tests
    
    # imputation for missing predictors --------
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # PS
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
          for(i_hh in 1:6) {
            for(i_inc in 1:2) {
              for (i_edu in 1:2) {
                p_pop[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]]
                  )
              }
            }
          }
        }
      }
    }
  }
  
  # pop avg
  p_avg <- inprod(N_pop[1:J], p_pop[1:J]) / sum(N_pop[1:J])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[1:J,i], p_pop[1:J]) / sum(zip_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[1:J,i], p_pop[1:J]) / sum(sex_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[1:J,i], p_pop[1:J]) / sum(age_mat[1:J,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[1:J,i], p_pop[1:J]) / sum(race_eth_mat[1:J,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[1:J,i], p_pop[1:J]) / sum(white_mat[1:J,i])
  }
  
  # hh avg
  for(i in 1:6) {
    p_hh[i] <- inprod(hh_mat[1:J,i], p_pop[1:J]) / sum(hh_mat[1:J,i])
  }
  
  # edu avg
  for(i in 1:2) {
    p_edu[i] <- inprod(edu_bin_mat[1:J,i], p_pop[1:J]) / sum(edu_bin_mat[1:J,i])
  }
  
  # income avg
  for(i in 1:2) {
    p_inc[i] <- inprod(income_bin_mat[1:J,i], p_pop[1:J]) / sum(income_bin_mat[1:J,i])
  }
  
  # income moa
  p_inc_rd <- p_inc[1] - p_inc[2] 
  p_inc_rr <- p_inc[1] / p_inc[2] 
  
  # education moa
  p_edu_rd <- p_edu[1] - p_edu[2] 
  p_edu_rr <- p_edu[1] / p_edu[2] 
  
  # white/non-white moa
  p_white_rd <- p_white[1] - p_white[2] 
  p_white_rr <- p_white[1] / p_white[2] 
  
})

bym2_logistic_allStrata_ixns_test_bias_simple_r12 <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf) 
  sigma_zip_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  # prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip) }
  for(i in 1:N_eth){ a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) }
  for(i in 1:2){ 
    a_inc[i] ~ dnorm(mean = 0, sd = sigma_inc)
    a_edu[i] ~ dnorm(mean = 0, sd = sigma_edu) 
  }
  for(i in 1:6){ a_hh[i] ~ dnorm(mean = 0, sd = sigma_hh) }
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(31*2)){ a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) }
  for(i in 1:(6*2)){ a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) }
  for(i in 1:(6*2)){ a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[1:N])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[1:N])))
  
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[1:N]) - b[3]*mean(x_zip[1:N]), scale = 1)
  
  # bym2
  for (i in 1:N_zip) {
    vc[i] ~ dnorm(0,1)
    Corr[i] <- sigma * uc[i] * sqrt(rho/scale)
    UCorr[i] <- sigma * vc[i] * sqrt((1-rho)) 
    bym2[i] <- Corr[i] + UCorr[i] 
  }
  
  uc[1:N_zip] ~ dcar_normal(adj=adj[1:L], weights=wei[1:L], 
                            num=num[1:N_zip], tau=1, zero_mean=1)  
  sigma ~ dunif(0,10)
  rho ~ dbeta(1,1)
  
  # test-bias --------
  
  # round 1
  # ortho
  y_sens_r1[1] ~ dbinom(prob = ilogit(logit_sens_r1[1]), size = n_sens_r1)
  y_spec_r1[1] ~ dbinom(prob = ilogit(logit_spec_r1[1]), size = n_spec_r1)
  
  logit_sens_r1[1] ~ dnorm(1.45, 0.6) # 2.5% 50% 97.5% 0.5695357 0.8107481 0.9316298 
  logit_spec_r1[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  se_net_r1 <- ilogit(logit_sens_r1[1])
  sp_net_r1 <- ilogit(logit_spec_r1[1])
  
  # round 2
  # ortho
  y_sens_r2[1] ~ dbinom(prob = ilogit(logit_sens_r2[1]), size = n_sens_r2[1])
  y_spec_r2[1] ~ dbinom(prob = ilogit(logit_spec_r2[1]), size = n_spec_r2[1])
  
  logit_sens_r2[1] ~ dnorm(2.1, 0.75) # 2.5% 50% 97.5% 0.6513360 0.8908752 0.9725546
  logit_spec_r2[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  # elisa
  y_sens_r2[2] ~ dbinom(prob = ilogit(logit_sens_r2[2]), size = n_sens_r2[2])
  y_spec_r2[2] ~ dbinom(prob = ilogit(logit_spec_r2[2]), size = n_spec_r2[2])
  
  logit_sens_r2[2] ~ dnorm(3.75, 1) # 2.5% 50% 97.5% 0.8554587 0.9769567 0.9966975 
  logit_spec_r2[2] ~ dnorm(5.25, sd = 1) # 2.5% 50% 97.5% 0.9640300 0.9947563 0.9992629 
  
  # net se/sp 
  se_net_r2 <- ilogit(logit_sens_r2[1]) * ilogit(logit_sens_r2[2])
  sp_net_r2 <- ilogit(logit_spec_r2[1]) + (1 - ilogit(logit_spec_r2[1])) * ilogit(logit_spec_r2[2])
  
  # round 1 & 2 cumulative
  se_net_r1_r2 <- se_net_r1 + se_net_r2 - (se_net_r1 * se_net_r2)
  sp_net_r1_r2 <- sp_net_r1 * sp_net_r2
  
  for(i in 1:N){
    
    # likelihood
    p[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
                     a_age[age[i]] + a_eth[eth[i]] + a_hh[hh[i]] + 
                     a_inc[inc[i]+1] + 
                     a_edu[edu[i]] + 
                     bym2[zip[i]] + 
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                     a_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                     a_eth_inc[index_eth_inc[eth[i],inc[i]+1]] +
                     a_eth_age[index_eth_age[eth[i],age[i]]]
    )
    
    y[i] ~ dbern(prob = p_sample[i])
    # prevalence corrected for imperfect tests
    p_sample[i] <- test_r1[i]*((p[i]*se_net_r1) + (1-p[i])*(1-sp_net_r1)) + 
                   test_r2[i]*((p[i]*se_net_r2) + (1-p[i])*(1-sp_net_r2)) +
                   test_r1_r2[i]*((p[i]*se_net_r1_r2) + (1-p[i])*(1-sp_net_r1_r2)) 
    
    # imputation for missing predictors --------
    
    # income
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # PS
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
          for(i_hh in 1:6) {
            for(i_inc in 1:2) {
              for (i_edu in 1:2) {
                p_pop[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]]
                  )
              }
            }
          }
        }
      }
    }
  }
  
  # pop avg
  p_avg <- inprod(N_pop[1:J], p_pop[1:J]) / sum(N_pop[1:J])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[1:J,i], p_pop[1:J]) / sum(zip_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[1:J,i], p_pop[1:J]) / sum(sex_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[1:J,i], p_pop[1:J]) / sum(age_mat[1:J,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[1:J,i], p_pop[1:J]) / sum(race_eth_mat[1:J,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[1:J,i], p_pop[1:J]) / sum(white_mat[1:J,i])
  }
  
  # hh avg
  for(i in 1:6) {
    p_hh[i] <- inprod(hh_mat[1:J,i], p_pop[1:J]) / sum(hh_mat[1:J,i])
  }
  
  # edu avg
  for(i in 1:2) {
    p_edu[i] <- inprod(edu_bin_mat[1:J,i], p_pop[1:J]) / sum(edu_bin_mat[1:J,i])
  }
  
  # income avg
  for(i in 1:2) {
    p_inc[i] <- inprod(income_bin_mat[1:J,i], p_pop[1:J]) / sum(income_bin_mat[1:J,i])
  }
  
  # income moa
  p_inc_rd <- p_inc[1] - p_inc[2] 
  p_inc_rr <- p_inc[1] / p_inc[2] 
  
  # education moa
  p_edu_rd <- p_edu[1] - p_edu[2] 
  p_edu_rr <- p_edu[1] / p_edu[2] 
  
  # white/non-white moa
  p_white_rd <- p_white[1] - p_white[2] 
  p_white_rr <- p_white[1] / p_white[2] 
  
})

bym2_logistic_allStrata_ixns_test_bias_simple_r123 <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf) 
  sigma_zip_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  # prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip) }
  for(i in 1:N_eth){ a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) }
  for(i in 1:2){ 
    a_inc[i] ~ dnorm(mean = 0, sd = sigma_inc)
    a_edu[i] ~ dnorm(mean = 0, sd = sigma_edu) 
  }
  for(i in 1:6){ a_hh[i] ~ dnorm(mean = 0, sd = sigma_hh) }
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(31*2)){ a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) }
  for(i in 1:(6*2)){ a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) }
  for(i in 1:(6*2)){ a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[1:N])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[1:N])))
  
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[1:N]) - b[3]*mean(x_zip[1:N]), scale = 1)
  
  # bym2
  for (i in 1:N_zip) {
    vc[i] ~ dnorm(0,1)
    Corr[i] <- sigma * uc[i] * sqrt(rho/scale)
    UCorr[i] <- sigma * vc[i] * sqrt((1-rho)) 
    bym2[i] <- Corr[i] + UCorr[i] 
  }
  
  uc[1:N_zip] ~ dcar_normal(adj=adj[1:L], weights=wei[1:L], num=num[1:N_zip], tau=1, zero_mean=1)  
  sigma ~ dunif(0,10)
  rho ~ dbeta(1,1)
  
  # test-bias --------
  
  # round 1
  # ortho
  y_sens_r1[1] ~ dbinom(prob = ilogit(logit_sens_r1[1]), size = n_sens_r1)
  y_spec_r1[1] ~ dbinom(prob = ilogit(logit_spec_r1[1]), size = n_spec_r1)
  
  logit_sens_r1[1] ~ dnorm(1.45, 0.6) # 2.5% 50% 97.5% 0.5695357 0.8107481 0.9316298 
  logit_spec_r1[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  se_net_r1 <- ilogit(logit_sens_r1[1])
  sp_net_r1 <- ilogit(logit_spec_r1[1])
  
  # round 2
  # ortho
  y_sens_r2[1] ~ dbinom(prob = ilogit(logit_sens_r2[1]), size = n_sens_r2[1])
  y_spec_r2[1] ~ dbinom(prob = ilogit(logit_spec_r2[1]), size = n_spec_r2[1])
  
  logit_sens_r2[1] ~ dnorm(2.1, 0.75) # 2.5% 50% 97.5% 0.6513360 0.8908752 0.9725546
  logit_spec_r2[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  # elisa
  y_sens_r2[2] ~ dbinom(prob = ilogit(logit_sens_r2[2]), size = n_sens_r2[2])
  y_spec_r2[2] ~ dbinom(prob = ilogit(logit_spec_r2[2]), size = n_spec_r2[2])
  
  logit_sens_r2[2] ~ dnorm(3.75, 1) # 2.5% 50% 97.5% 0.8554587 0.9769567 0.9966975 
  logit_spec_r2[2] ~ dnorm(5.25, 1) # 2.5% 50% 97.5% 0.9640300 0.9947563 0.9992629 
  
  # net se/sp 
  se_net_r2 <- ilogit(logit_sens_r2[1]) * ilogit(logit_sens_r2[2])
  sp_net_r2 <- ilogit(logit_spec_r2[1]) + (1 - ilogit(logit_spec_r2[1])) * ilogit(logit_spec_r2[2])
  # round 1 & 2 cumulative
  se_net_r1_r2 <- se_net_r1 + se_net_r2 - (se_net_r1 * se_net_r2)
  sp_net_r1_r2 <- sp_net_r1 * sp_net_r2
  
  # round 3
  # ortho
  y_sens_r3[1] ~ dbinom(prob = ilogit(logit_sens_r3[1]), size = n_sens_r3[1])
  y_spec_r3[1] ~ dbinom(prob = ilogit(logit_spec_r3[1]), size = n_spec_r3[1])
  
  logit_sens_r3[1] ~ dnorm(2.1, 0.75) # 2.5% 50% 97.5% 0.6513360 0.8908752 0.9725546
  logit_spec_r3[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  # elisa
  y_sens_r3[2] ~ dbinom(prob = ilogit(logit_sens_r3[2]), size = n_sens_r3[2])
  y_spec_r3[2] ~ dbinom(prob = ilogit(logit_spec_r3[2]), size = n_spec_r3[2])
  
  logit_sens_r3[2] ~ dnorm(1.9, 0.6) # 2.5% 50% 97.5% 0.8554587 0.9769567 0.9966975 
  logit_spec_r3[2] ~ dnorm(3.75, 0.6) # 2.5% 50% 97.5% 0.9640300 0.9947563 0.9992629 
  
  # net se/sp 
  se_net_r3 <- ilogit(logit_sens_r3[1]) * ilogit(logit_sens_r3[2])
  sp_net_r3 <- ilogit(logit_spec_r3[1]) + (1 - ilogit(logit_spec_r3[1])) * ilogit(logit_spec_r3[2])
  
  # round 1,2,3 cumulative
  se_net_r1_r2_r3 <- se_net_r1_r2 + se_net_r3 - (se_net_r1_r2 * se_net_r3)
  sp_net_r1_r2_r3 <- se_net_r1_r2 * sp_net_r3
  
  # round 1,3 cumulative
  se_net_r1_r3 <- se_net_r1 + se_net_r3 - (se_net_r1 * se_net_r3)
  sp_net_r1_r3 <- sp_net_r1 * sp_net_r3
  
  for(i in 1:N){
    
    # likelihood
    p[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
                     a_age[age[i]] + a_eth[eth[i]] + a_hh[hh[i]] + 
                     a_inc[inc[i]+1] + 
                     a_edu[edu[i]] + 
                     bym2[zip[i]] + 
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                     a_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                     a_eth_inc[index_eth_inc[eth[i],inc[i]+1]] +
                     a_eth_age[index_eth_age[eth[i],age[i]]]
    )
    y[i] ~ dbern(prob = p_sample[i])
    
    # prevalence corrected for imperfect tests
    p_sample[i] <- 
      test_r1[i]*((p[i]*se_net_r1) + (1-p[i])*(1-sp_net_r1)) + 
      test_r2[i]*((p[i]*se_net_r2) + (1-p[i])*(1-sp_net_r2)) +
      test_r3[i]*((p[i]*se_net_r3) + (1-p[i])*(1-sp_net_r3)) +
      test_r1_r2[i]*((p[i]*se_net_r1_r2) + (1-p[i])*(1-sp_net_r1_r2)) + 
      test_r2_r3[i]*((p[i]*se_net_r2_r3) + (1-p[i])*(1-sp_net_r2_r3)) +
      test_r1_r3[i]*((p[i]*se_net_r1_r3) + (1-p[i])*(1-sp_net_r1_r3)) +
      test_r1_r2_r3[i]*((p[i]*se_net_r1_r2_r3) + (1-p[i])*(1-sp_net_r1_r2_r3))
    
    
    # imputation for missing predictors --------
    
    # income
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # PS
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
          for(i_hh in 1:6) {
            for(i_inc in 1:2) {
              for (i_edu in 1:2) {
                p_pop[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]]
                  )
              }
            }
          }
        }
      }
    }
  }
  
  # pop avg
  p_avg <- inprod(N_pop[1:J], p_pop[1:J]) / sum(N_pop[1:J])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[1:J,i], p_pop[1:J]) / sum(zip_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[1:J,i], p_pop[1:J]) / sum(sex_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[1:J,i], p_pop[1:J]) / sum(age_mat[1:J,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[1:J,i], p_pop[1:J]) / sum(race_eth_mat[1:J,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[1:J,i], p_pop[1:J]) / sum(white_mat[1:J,i])
  }
  
  # hh avg
  for(i in 1:6) {
    p_hh[i] <- inprod(hh_mat[1:J,i], p_pop[1:J]) / sum(hh_mat[1:J,i])
  }
  
  # edu avg
  for(i in 1:2) {
    p_edu[i] <- inprod(edu_bin_mat[1:J,i], p_pop[1:J]) / sum(edu_bin_mat[1:J,i])
  }
  
  # income avg
  for(i in 1:2) {
    p_inc[i] <- inprod(income_bin_mat[1:J,i], p_pop[1:J]) / sum(income_bin_mat[1:J,i])
  }
  
  # income moa
  p_inc_rd <- p_inc[1] - p_inc[2] 
  p_inc_rr <- p_inc[1] / p_inc[2] 
  
  # education moa
  p_edu_rd <- p_edu[1] - p_edu[2] 
  p_edu_rr <- p_edu[1] / p_edu[2] 
  
  # white/non-white moa
  p_white_rd <- p_white[1] - p_white[2] 
  p_white_rr <- p_white[1] / p_white[2] 
  
})

bym2_logistic_allStrata_ixns_test_bias_simple_r123_spike <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf) 
  sigma_zip_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  l# prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip) }
  for(i in 1:N_eth){ a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) }
  for(i in 1:2){ 
    a_inc[i] ~ dnorm(mean = 0, sd = sigma_inc)
    a_edu[i] ~ dnorm(mean = 0, sd = sigma_edu) 
  }
  for(i in 1:6){ a_hh[i] ~ dnorm(mean = 0, sd = sigma_hh) }
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(31*2)){ a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) }
  for(i in 1:(6*2)){ a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) }
  for(i in 1:(6*2)){ a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[1:N]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[1:N])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[1:N])))
  
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[1:N]) - b[3]*mean(x_zip[1:N]), scale = 1)
  
  # bym2
  for (i in 1:N_zip) {
    vc[i] ~ dnorm(0,1)
    Corr[i] <- sigma * uc[i] * sqrt(rho/scale)
    UCorr[i] <- sigma * vc[i] * sqrt((1-rho)) 
    bym2[i] <- Corr[i] + UCorr[i] 
  }
  
  uc[1:N_zip] ~ dcar_normal(adj=adj[1:L], weights=wei[1:L], num=num[1:N_zip], tau=1, zero_mean=1)  
  sigma ~ dunif(0,10)
  rho ~ dbeta(1,1)
  
  # test-bias --------
  
  # round 1
  # ortho
  y_sens_r1[1] ~ dbinom(prob = ilogit(logit_sens_r1[1]), size = n_sens_r1)
  y_spec_r1[1] ~ dbinom(prob = ilogit(logit_spec_r1[1]), size = n_spec_r1)
  
  logit_sens_r1[1] ~ dnorm(1.45, 0.6) # 2.5% 50% 97.5% 0.5695357 0.8107481 0.9316298 
  logit_spec_r1[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  se_net_r1 <- ilogit(logit_sens_r1[1])
  sp_net_r1 <- ilogit(logit_spec_r1[1])
  
  # round 2
  # ortho
  y_sens_r2[1] ~ dbinom(prob = ilogit(logit_sens_r2[1]), size = n_sens_r2[1])
  y_spec_r2[1] ~ dbinom(prob = ilogit(logit_spec_r2[1]), size = n_spec_r2[1])
  
  logit_sens_r2[1] ~ dnorm(2.1, 0.75) # 2.5% 50% 97.5% 0.6513360 0.8908752 0.9725546
  logit_spec_r2[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  # elisa
  y_sens_r2[2] ~ dbinom(prob = ilogit(logit_sens_r2[2]), size = n_sens_r2[2])
  y_spec_r2[2] ~ dbinom(prob = ilogit(logit_spec_r2[2]), size = n_spec_r2[2])
  
  logit_sens_r2[2] ~ dnorm(3.75, 1) # 2.5% 50% 97.5% 0.8554587 0.9769567 0.9966975 
  logit_spec_r2[2] ~ dnorm(5.25, 1) # 2.5% 50% 97.5% 0.9640300 0.9947563 0.9992629 
  
  # net se/sp 
  se_net_r2 <- ilogit(logit_sens_r2[1]) * ilogit(logit_sens_r2[2])
  sp_net_r2 <- ilogit(logit_spec_r2[1]) + (1 - ilogit(logit_spec_r2[1])) * ilogit(logit_spec_r2[2])
  # round 1 & 2 cumulative
  se_net_r1_r2 <- se_net_r1 + se_net_r2 - (se_net_r1 * se_net_r2)
  sp_net_r1_r2 <- sp_net_r1 * sp_net_r2
  
  # round 3
  # ortho
  y_sens_r3[1] ~ dbinom(prob = ilogit(logit_sens_r3[1]), size = n_sens_r3[1])
  y_spec_r3[1] ~ dbinom(prob = ilogit(logit_spec_r3[1]), size = n_spec_r3[1])
  
  logit_sens_r3[1] ~ dnorm(2.1, 0.75) # 2.5% 50% 97.5% 0.6513360 0.8908752 0.9725546
  logit_spec_r3[1] ~ dnorm(5, 2.25) # 2.5% 50% 97.5% 0.6393143 0.9933017 0.9999159 
  
  # elisa
  y_sens_r3[2] ~ dbinom(prob = ilogit(logit_sens_r3[2]), size = n_sens_r3[2])
  y_spec_r3[2] ~ dbinom(prob = ilogit(logit_spec_r3[2]), size = n_spec_r3[2])
  
  logit_sens_r3[2] ~ dnorm(3.75, 1) # 2.5% 50% 97.5% 0.8554587 0.9769567 0.9966975 
  logit_spec_r3[2] ~ dnorm(5.25, 1) # 2.5% 50% 97.5% 0.9640300 0.9947563 0.9992629 
  
  # net se/sp 
  se_net_r3 <- ilogit(logit_sens_r3[1]) * ilogit(logit_sens_r3[2])
  sp_net_r3 <- ilogit(logit_spec_r3[1]) + (1 - ilogit(logit_spec_r3[1])) * ilogit(logit_spec_r3[2])
  
  # round 1,2,3 cumulative
  se_net_r1_r2_r3 <- se_net_r1_r2 + se_net_r3 - (se_net_r1_r2 * se_net_r3)
  sp_net_r1_r2_r3 <- se_net_r1_r2 * sp_net_r3
  
  # round 1,3 cumulative
  se_net_r1_r3 <- se_net_r1 + se_net_r3 - (se_net_r1 * se_net_r3)
  sp_net_r1_r3 <- sp_net_r1 * sp_net_r3
  
  for(i in 1:N){
    
    # likelihood
    p[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
                     a_age[age[i]] + a_eth[eth[i]] + a_hh[hh[i]] + 
                     a_inc[inc[i]+1] + 
                     a_edu[edu[i]] + 
                     bym2[zip[i]] + 
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                     a_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                     a_eth_inc[index_eth_inc[eth[i],inc[i]+1]] +
                     a_eth_age[index_eth_age[eth[i],age[i]]]
    )
    #y[i] ~ dbern(prob = p[i])
    y[i] ~ dbern(prob = p_sample[i])
    
    # prevalence corrected for imperfect tests
    p_sample[i] <- 
      test_r1[i]*((p[i]*se_net_r1) + (1-p[i])*(1-sp_net_r1)) + 
      test_r2[i]*((p[i]*se_net_r2) + (1-p[i])*(1-sp_net_r2)) +
      test_r3[i]*((p[i]*se_net_r3) + (1-p[i])*(1-sp_net_r3)) +
      test_r1_r2[i]*((p[i]*se_net_r1_r2) + (1-p[i])*(1-sp_net_r1_r2)) + 
      test_r2_r3[i]*((p[i]*se_net_r2_r3) + (1-p[i])*(1-sp_net_r2_r3)) +
      test_r1_r3[i]*((p[i]*se_net_r1_r3) + (1-p[i])*(1-sp_net_r1_r3)) +
      test_r1_r2_r3[i]*((p[i]*se_net_r1_r2_r3) + (1-p[i])*(1-sp_net_r1_r2_r3))
    
    
    # imputation for missing predictors --------
    
    # income
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # PS
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
          for(i_hh in 1:6) {
            for(i_inc in 1:2) {
              for (i_edu in 1:2) {
                p_pop[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]]
                  )
              }
            }
          }
        }
      }
    }
  }
  
  # pop avg
  p_avg <- inprod(N_pop[1:J], p_pop[1:J]) / sum(N_pop[1:J])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[1:J,i], p_pop[1:J]) / sum(zip_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[1:J,i], p_pop[1:J]) / sum(sex_mat[1:J,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[1:J,i], p_pop[1:J]) / sum(age_mat[1:J,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[1:J,i], p_pop[1:J]) / sum(race_eth_mat[1:J,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[1:J,i], p_pop[1:J]) / sum(white_mat[1:J,i])
  }
  
  # hh avg
  for(i in 1:6) {
    p_hh[i] <- inprod(hh_mat[1:J,i], p_pop[1:J]) / sum(hh_mat[1:J,i])
  }
  
  # edu avg
  for(i in 1:2) {
    p_edu[i] <- inprod(edu_bin_mat[1:J,i], p_pop[1:J]) / sum(edu_bin_mat[1:J,i])
  }
  
  # income avg
  for(i in 1:2) {
    p_inc[i] <- inprod(income_bin_mat[1:J,i], p_pop[1:J]) / sum(income_bin_mat[1:J,i])
  }
  
  # income moa
  p_inc_rd <- p_inc[1] - p_inc[2] 
  p_inc_rr <- p_inc[1] / p_inc[2] 
  
  # education moa
  p_edu_rd <- p_edu[1] - p_edu[2] 
  p_edu_rr <- p_edu[1] / p_edu[2] 
  
  # white/non-white moa
  p_white_rd <- p_white[1] - p_white[2] 
  p_white_rr <- p_white[1] / p_white[2] 
  
})




