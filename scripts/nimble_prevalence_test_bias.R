

# bym2 interactions, test-bias 

# Round 1
bym2_logistic_allStrata_ixns_test_bias_r1 <- nimbleCode({
  
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(mean = mu_logit_sens[1], sd = sigma_logit_sens[1]) # normal on log-odds scale
  logit_spec[1] ~ dnorm(mean = mu_logit_spec[1], sd = sigma_logit_spec[1]) # normal on log-odds scale
  
  mu_logit_sens[1] ~ dnorm(4,2)
  mu_logit_spec[1] ~ dnorm(2,1) 
  
  sigma_logit_sens[1] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[1] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) 
  spec_net <- ilogit(logit_spec[1]) 
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests
    
    # imputation for missing predictors --------
    
    # income
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # Poststratification
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
  
  # generate parameters of interest from PS strata predictions
  
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

})

# Round 1+2
bym2_logistic_allStrata_ixns_test_bias_r12 <- nimbleCode({
  
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(mean = mu_logit_sens[1], sd = sigma_logit_sens[1]) # normal on log-odds scale
  logit_spec[1] ~ dnorm(mean = mu_logit_spec[1], sd = sigma_logit_spec[1]) # normal on log-odds scale
  
  mu_logit_sens[1] ~ dnorm(4,2)
  mu_logit_spec[1] ~ dnorm(2,1) 
  
  sigma_logit_sens[1] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[1] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # elisa
  y_sens[2] ~ dbinom(prob = ilogit(logit_sens[2]), size = n_sens[2])
  y_spec[2] ~ dbinom(prob = ilogit(logit_spec[2]), size = n_spec[2])
  
  logit_sens[2] ~ dnorm(mean = mu_logit_sens[2], sd = sigma_logit_sens[2]) 
  logit_spec[2] ~ dnorm(mean = mu_logit_spec[2], sd = sigma_logit_spec[2]) 
  
  mu_logit_sens[2] ~ dnorm(4,2) 
  mu_logit_spec[2] ~ dnorm(4,2)

  sigma_logit_sens[2] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[2] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) * ilogit(logit_sens[2])
  spec_net <- ilogit(logit_spec[1]) + (1 - ilogit(logit_spec[1])) * ilogit(logit_spec[2])
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests

    # imputation for missing predictors --------
    
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # Poststratification
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
  
  # generate parameters of interest from PS strata predictions
  
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
  
})

# Round 1+2+3
bym2_logistic_allStrata_ixns_test_bias_r123 <- nimbleCode({
  
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(mean = mu_logit_sens[1], sd = sigma_logit_sens[1]) # normal on log-odds scale
  logit_spec[1] ~ dnorm(mean = mu_logit_spec[1], sd = sigma_logit_spec[1]) # normal on log-odds scale
  
  mu_logit_sens[1] ~ dnorm(4,2)
  mu_logit_spec[1] ~ dnorm(2,1) 
  
  sigma_logit_sens[1] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[1] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # roche
  y_sens[2] ~ dbinom(prob = ilogit(logit_sens[2]), size = n_sens[2])
  y_spec[2] ~ dbinom(prob = ilogit(logit_spec[2]), size = n_spec[2])
  
  logit_sens[2] ~ dnorm(mean = mu_logit_sens[2], sd = sigma_logit_sens[2]) 
  logit_spec[2] ~ dnorm(mean = mu_logit_spec[2], sd = sigma_logit_spec[2]) 
  
  mu_logit_sens[2] ~ dnorm(2,1) 
  mu_logit_spec[2] ~ dnorm(4,2)
  
  sigma_logit_sens[2] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[2] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) * ilogit(logit_sens[2])
  spec_net <- ilogit(logit_spec[1]) + (1 - ilogit(logit_spec[1])) * ilogit(logit_spec[2])
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests
     # imputation for missing predictors --------
    
    # income
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
  }
  
  # Poststratification
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
  
  # generate parameters of interest from PS strata predictions

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

# bym2 interactions, test-bias simple --------
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(1.45, 0.5) #     2.5%       50%     97.5% 0.6182797 0.8100915 0.9181854 
  logit_spec[1] ~ dnorm(5, 2.25) #   2.5%       50%     97.5% 0.6393143 0.9933017 0.9999159 
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) 
  spec_net <- ilogit(logit_spec[1]) 
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests
    
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

bym2_logistic_allStrata_ixns_test_bias_r12 <- nimbleCode({
  
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(mean = mu_logit_sens[1], sd = sigma_logit_sens[1]) # normal on log-odds scale
  logit_spec[1] ~ dnorm(mean = mu_logit_spec[1], sd = sigma_logit_spec[1]) # normal on log-odds scale
  
  mu_logit_sens[1] ~ dnorm(4,2)
  mu_logit_spec[1] ~ dnorm(2,1) 
  
  sigma_logit_sens[1] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[1] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # elisa
  y_sens[2] ~ dbinom(prob = ilogit(logit_sens[2]), size = n_sens[2])
  y_spec[2] ~ dbinom(prob = ilogit(logit_spec[2]), size = n_spec[2])
  
  logit_sens[2] ~ dnorm(mean = mu_logit_sens[2], sd = sigma_logit_sens[2]) 
  logit_spec[2] ~ dnorm(mean = mu_logit_spec[2], sd = sigma_logit_spec[2]) 
  
  mu_logit_sens[2] ~ dnorm(4,2) 
  mu_logit_spec[2] ~ dnorm(4,2)
  
  sigma_logit_sens[2] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[2] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) * ilogit(logit_sens[2])
  spec_net <- ilogit(logit_spec[1]) + (1 - ilogit(logit_spec[1])) * ilogit(logit_spec[2])
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests
    
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

bym2_logistic_allStrata_ixns_test_bias_r123 <- nimbleCode({
  
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
  y_sens[1] ~ dbinom(prob = ilogit(logit_sens[1]), size = n_sens[1])
  y_spec[1] ~ dbinom(prob = ilogit(logit_spec[1]), size = n_spec[1])
  
  logit_sens[1] ~ dnorm(mean = mu_logit_sens[1], sd = sigma_logit_sens[1]) # normal on log-odds scale
  logit_spec[1] ~ dnorm(mean = mu_logit_spec[1], sd = sigma_logit_spec[1]) # normal on log-odds scale
  
  mu_logit_sens[1] ~ dnorm(4,2)
  mu_logit_spec[1] ~ dnorm(2,1) 
  
  sigma_logit_sens[1] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[1] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # roche
  y_sens[2] ~ dbinom(prob = ilogit(logit_sens[2]), size = n_sens[2])
  y_spec[2] ~ dbinom(prob = ilogit(logit_spec[2]), size = n_spec[2])
  
  logit_sens[2] ~ dnorm(mean = mu_logit_sens[2], sd = sigma_logit_sens[2]) 
  logit_spec[2] ~ dnorm(mean = mu_logit_spec[2], sd = sigma_logit_spec[2]) 
  
  mu_logit_sens[2] ~ dnorm(2,1) 
  mu_logit_spec[2] ~ dnorm(4,2)
  
  sigma_logit_sens[2] ~ T(dnorm(0, logit_sens_prior_scale), min=0, max=Inf)
  sigma_logit_spec[2] ~ T(dnorm(0, logit_spec_prior_scale), min=0, max=Inf)
  
  # net se/sp 
  sens_net <- ilogit(logit_sens[1]) * ilogit(logit_sens[2])
  spec_net <- ilogit(logit_spec[1]) + (1 - ilogit(logit_spec[1])) * ilogit(logit_spec[2])
  
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
    p_sample[i] <- (p[i] * sens_net) + (1 - p[i]) * (1 - spec_net) # prevalence corrected for imperfect tests
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
# seleciton model

bym2_logistic_allStrata_ixn_selection_model <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  #sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # zip code
  sigma_zip_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
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
  
  sigma_age_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_eth_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_edu_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age_g ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  # prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ g_zip[i] ~ dnorm(mean = 0, sd = sigma_zip_g) }
  for(i in 1:N_eth){ 
    a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) 
    g_eth[i] ~ dnorm(mean = 0, sd = sigma_eth_g)
  }
  for(i in 1:N_age){ 
    g_age[i] ~ dnorm(mean = 0, sd = sigma_age_g)   
  }
  for(i in 1:2){ 
    a_inc[i] ~ dnorm(mean = 0, sd = sigma_inc)
    g_inc[i] ~ dnorm(mean = 0, sd = sigma_inc_g)
    a_edu[i] ~ dnorm(mean = 0, sd = sigma_edu) 
    g_edu[i] ~ dnorm(mean = 0, sd = sigma_edu_g) 
    
  }
  for(i in 1:6){ 
    a_hh[i] ~ dnorm(mean = 0, sd = sigma_hh) 
    g_hh[i] ~ dnorm(mean = 0, sd = sigma_hh_g) 
  }
  for(i in 1:(31*6)){ 
    a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) 
    g_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth_g) 
  }
  for(i in 1:(31*2)){ 
    a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) 
    g_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu_g) 
  }
  for(i in 1:(6*2)){ 
    a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) 
    g_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu_g) 
  }
  for(i in 1:(6*2)){ 
    a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) 
    g_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc_g)
  }
  for(i in 1:(6*5)){ 
    a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) 
    g_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age_g) 
  }
  
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  
  g[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  g[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  
 
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[]) - b[3]*mean(x_zip[]), scale = 1)
  g[1] ~ dlogis(location = 0 - g[2]*mean(male[]) - g[3]*mean(x_zip[]), scale = 1)
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[])))
  
  for (i in 1:N_zip) {
    vc[i] ~ dnorm(0,1)
    Corr[i] <- sigma * uc[i] * sqrt(rho/scale)
    UCorr[i] <- sigma * vc[i] * sqrt((1-rho)) 
    bym2[i] <- Corr[i] + UCorr[i]
  }
  
  uc[1:N_zip] ~ dcar_normal(adj=adj[1:L], weights=wei[1:L], num=num[1:N_zip], tau=1, zero_mean=1)  
  sigma ~ dunif(0,10)
  rho ~ dbeta(1,1)
  
  for(i in 1:N){
    
    
    # likelihood
    
    # full likelihood
    # P(S,Y) = P(S|Y)P(Y)
    
    # P(S|Y)
    s[i] ~ dbern(p_s[i])
    p_s[i] <- ilogit(g[1] + g[2]*(male[i]) + g[3]*(x_zip[i]) + 
                       g_age[age[i]] + g_eth[eth[i]] + g_hh[hh[i]] + 
                       g_inc[inc[i]+1] + g_edu[edu[i]] + g_zip[zip[i]] + 
                       g_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                       g_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                       g_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                       g_eth_inc[index_eth_inc[eth[i],inc[i]+1]] +
                       g_eth_age[index_eth_age[eth[i],age[i]]]
                     )
    
    # P(Y)
    y[i] ~ dbern(p_y[i])
    p_y[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
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
  
  # P(Y|S=s)
  p_sel1 <- inprod(sel_mat[1:N,1], p_y[1:N]) / sum(sel_mat[1:N,1])
  p_sel0 <- inprod(sel_mat[1:N,2], p_y[1:N]) / sum(sel_mat[1:N,2])
  
  
  # pop avg
  p_avg <- inprod(N_pop[], p_pop[]) / sum(N_pop[])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[,i], p_pop[]) / sum(zip_mat[,i])
  }
  
  # age avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[,i], p_pop[]) / sum(sex_mat[,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[,i], p_pop[]) / sum(age_mat[,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[,i], p_pop[]) / sum(race_eth_mat[,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[,i], p_pop[]) / sum(white_mat[,i])
  }
  
  # hh avg
  for(i in 1:6) {
    p_hh[i] <- inprod(hh_mat[,i], p_pop[]) / sum(hh_mat[,i])
  }
  
  # edu avg
  for(i in 1:2) {
    p_edu[i] <- inprod(edu_bin_mat[,i], p_pop[]) / sum(edu_bin_mat[,i])
  }
  
  # income avg
  for(i in 1:2) {
    p_inc[i] <- inprod(income_bin_mat[,i], p_pop[]) / sum(income_bin_mat[,i])
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


# va --------
bym2_logistic_allStrata_ixns_var4 <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # zip code
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # age category
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) # ethnicity
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) 
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0) 
  
  # prior for AR coefficient
  ar_rho ~ dbeta(0.5, 0.5) # same as uniform 0,1
  
  # AR coefficient:
  rho_transformed <- (ar_rho*2) - 1
  
  # random intercepts
  for(i in 1:N_zip){ a_zip[i] ~ dnorm(mean = 0, sd = sigma_zip) }
  for(i in 1:N_eth){ a_eth[i] ~ dnorm(mean = 0, sd = sigma_eth) }
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  
  # prior on centered intercept
  b[1] ~ dlogis(location = 0 - b[2]*mean(male[]) - b[3]*mean(x_zip[]), scale = 1)
  
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
  
  for(i in 1:N){
    
    # likelihood
    y[i] ~ dbern(p[i])
    p[i] <- ilogit(b[1] + b[2]*(male[i]) + b[3]*(x_zip[i]) + 
                     a_age[age[i]] + a_eth[eth[i]] + 
                     bym2[zip[i]] + 
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_eth_age[index_eth_age[eth[i],age[i]]]
                   )
                     
  }
  
  # PS
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
                p_pop[count_index[i_zip, i_age, i_eth, i_male]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + bym2[i_zip] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_eth_age[index_eth_age[i_eth,i_age]]
                         )
        }
      }
    }
  }
  
  # pop avg
  p_avg <- inprod(N_pop[], p_pop[]) / sum(N_pop[])
  
  # zip average
  for(i in 1:N_zip) {
    p_zip[i] <- inprod(zip_mat[,i], p_pop[]) / sum(zip_mat[,i])
  }
  
  # sex avg
  for(i in 1:2) {
    p_sex[i] <- inprod(sex_mat[,i], p_pop[]) / sum(sex_mat[,i])
  }
  
  # age avg
  for(i in 1:N_age) {
    p_age[i] <- inprod(age_mat[,i], p_pop[]) / sum(age_mat[,i])
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth[i] <- inprod(race_eth_mat[,i], p_pop[]) / sum(race_eth_mat[,i])
  }
  
  # white avg
  for(i in 1:2) {
    p_white[i] <- inprod(white_mat[,i], p_pop[]) / sum(white_mat[,i])
  }
  
  # white/non-white moa
  p_white_rd <- p_white[1] - p_white[2] 
  p_white_rr <- p_white[1] / p_white[2] 

})



