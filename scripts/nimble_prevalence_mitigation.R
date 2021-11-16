
bym2_logistic_allStrata_ixns_mit <- nimbleCode({
  
  # hyperpriors for sd of random intercepts
  sigma_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_hh ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_mit ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  sigma_zip_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_zip_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_eth_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
  sigma_mit_zip ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_mit_age ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_mit_eth ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_mit_inc ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  sigma_mit_edu ~ T(dnorm(mean = 0, sd = coef_prior_scale), min = 0, max=Inf)
  
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
  for(i in 1:2){ a_mit[i] ~ dnorm(mean = 0, sd = sigma_mit) }
  
  for(i in 1:(31*6)){ a_zip_eth[i] ~ dnorm(mean = 0, sd = sigma_zip_eth) }
  for(i in 1:(31*2)){ a_zip_edu[i] ~ dnorm(mean = 0, sd = sigma_zip_edu) }
  for(i in 1:(6*2)){ a_eth_edu[i] ~ dnorm(mean = 0, sd = sigma_eth_edu) }
  for(i in 1:(6*2)){ a_eth_inc[i] ~ dnorm(mean = 0, sd = sigma_eth_inc) }
  for(i in 1:(6*5)){ a_eth_age[i] ~ dnorm(mean = 0, sd = sigma_eth_age) }
  
  for(i in 1:(2*31)){ a_mit_zip[i] ~ dnorm(mean = 0, sd = sigma_mit_zip) }
  for(i in 1:(2*5)){ a_mit_age[i] ~ dnorm(mean = 0, sd = sigma_mit_age) }
  for(i in 1:(2*6)){ a_mit_eth[i] ~ dnorm(mean = 0, sd = sigma_mit_eth) }
  for(i in 1:(2*2)){ a_mit_inc[i] ~ dnorm(mean = 0, sd = sigma_mit_inc) }
  for(i in 1:(2*2)){ a_mit_edu[i] ~ dnorm(mean = 0, sd = sigma_mit_edu) }
  
  # incorporate structured AR(1) prior for age
  a_age[1] ~ dnorm(mean = 0, 1/sqrt(1 - rho_transformed^2))
  for(i in 2:N_age) { a_age[i] ~ dnorm(rho_transformed*a_age[i-1], 1) } 
  
  # coefficient priors
  b[2] ~ dnorm(mean = 0, sd = coef_prior_scale)
  b[3] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  
  w[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  w[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  w[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[])))
  w[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[])))
  
  v[1] ~ dnorm(mean = 0, sd = coef_prior_scale)
  v[2] ~ dnorm(mean = 0, sd = coef_prior_scale / sd(x_zip[]))
  v[3] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(male[])))
  v[4] ~ dnorm(mean = 0, sd = coef_prior_scale / (2*sd(edu[])))
  
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
                     a_age[age[i]] + a_eth[eth[i]] + a_hh[hh[i]] + 
                     a_inc[inc[i]+1] + 
                     a_edu[edu[i]] + 
                     bym2[zip[i]] + 
                     a_mit[mit[i]+1] +
                     a_zip_eth[index_zip_eth[zip[i],eth[i]]] + 
                     a_zip_edu[index_zip_edu[zip[i],edu[i]]] + 
                     a_eth_edu[index_eth_edu[eth[i],edu[i]]] +
                     a_eth_inc[index_eth_inc[eth[i],inc[i]+1]] +
                     a_eth_age[index_eth_age[eth[i],age[i]]] + 
                     a_mit_zip[index_mit_zip[mit[i]+1,zip[i]]] + 
                     a_mit_age[index_mit_age[mit[i]+1,age[i]]] + 
                     a_mit_eth[index_mit_eth[mit[i]+1,eth[i]]] + 
                     a_mit_inc[index_mit_inc[mit[i]+1,inc[i]+1]] + 
                     a_mit_edu[index_mit_edu[mit[i]+1,edu[i]]]
    )
    
    # imputation for missing predictors --------
    inc[i] ~ dbern(p_inc_imp[i])
    p_inc_imp[i] <- ilogit(w[1] + w[2]*(x_zip[i]) + w[3]*(male[i]) + w[4]*(edu[i]-1))
    
    mit[i] ~ dbern(p_mit_imp[i])
    p_mit_imp[i] <- ilogit(v[1] + v[2]*(x_zip[i]) + v[3]*(male[i]) + v[4]*(edu[i]-1))
    
  }
  
  # Poststratification
  for (i_zip in 1:N_zip) {
    for (i_age in 1:N_age) {
      for (i_eth in 1:N_eth) {
        for (i_male in 1:2) {
          for(i_hh in 1:6) {
            for(i_inc in 1:2) {
              for (i_edu in 1:2) {
                
                p_pop_mit1[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_mit[1] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]] + 
                           a_mit_zip[index_mit_zip[1,i_zip]] + 
                           a_mit_age[index_mit_age[1,i_age]] + 
                           a_mit_eth[index_mit_eth[1,i_eth]] + 
                           a_mit_inc[index_mit_inc[1,i_inc]] + 
                           a_mit_edu[index_mit_edu[1,i_edu]]
                  )
                
                p_pop_mit2[count_index[i_zip, i_age, i_eth, i_male, i_hh, i_inc, i_edu]] <- 
                  ilogit(b[1] + b[2]*(i_male-1) + b[3]*x_zip_uniq[i_zip] + 
                           a_eth[i_eth] + a_age[i_age] + 
                           a_inc[i_inc] + 
                           a_edu[i_edu] + a_hh[i_hh] + bym2[i_zip] + 
                           a_mit[2] + 
                           a_zip_eth[index_zip_eth[i_zip,i_eth]] + 
                           a_zip_edu[index_zip_edu[i_zip,i_edu]] + 
                           a_eth_edu[index_eth_edu[i_eth,i_edu]] + 
                           a_eth_inc[index_eth_inc[i_eth,i_inc]] +
                           a_eth_age[index_eth_age[i_eth,i_age]] + 
                           a_mit_zip[index_mit_zip[2,i_zip]] + 
                           a_mit_age[index_mit_age[2,i_age]] + 
                           a_mit_eth[index_mit_eth[2,i_eth]] + 
                           a_mit_inc[index_mit_inc[2,i_inc]] + 
                           a_mit_edu[index_mit_edu[2,i_edu]]
                  )
              }
            }
          }
        }
      }
    }
  }
  
  # generate parameters of interest from PS strata predictions
  
  # region
  p_avg_mit1 <- inprod(N_pop[1:J], p_pop_mit1[1:J]) / sum(N_pop[1:J])
  p_avg_mit2 <- inprod(N_pop[1:J], p_pop_mit2[1:J]) / sum(N_pop[1:J])
  p_avg_mit_rd <- p_avg_mit2 - p_avg_mit1
  p_avg_mit_rr <- p_avg_mit2 / p_avg_mit1
  
  
  # zip 
  for(i in 1:N_zip) {
    p_zip_mit1[i] <- inprod(zip_mat[,i], p_pop_mit1[1:J]) / sum(zip_mat[,i])
    p_zip_mit2[i] <- inprod(zip_mat[,i], p_pop_mit2[1:J]) / sum(zip_mat[,i])
    p_zip_mit_rd[i] <- p_zip_mit2[i] - p_zip_mit1[i]
    p_zip_mit_rr[i] <- p_zip_mit2[i] / p_zip_mit1[i]
    
  }
  
  # sex 
  for(i in 1:2) {
    p_sex_mit1[i] <- inprod(sex_mat[,i], p_pop_mit1[1:J]) / sum(sex_mat[,i])
    p_sex_mit2[i] <- inprod(sex_mat[,i], p_pop_mit2[1:J]) / sum(sex_mat[,i])
    p_sex_mit_rd[i] <- p_sex_mit2[i] - p_sex_mit1[i]
    p_sex_mit_rr[i] <- p_sex_mit2[i] - p_sex_mit1[i]
    
  }
  
  # age 
  for(i in 1:N_age) {
    p_age_mit1[i] <- inprod(age_mat[,i], p_pop_mit1[1:J]) / sum(age_mat[,i])
    p_age_mit2[i] <- inprod(age_mat[,i], p_pop_mit2[1:J]) / sum(age_mat[,i])
    p_age_mit_rd[i] <- p_age_mit2[i] - p_age_mit1[i]
    p_age_mit_rr[i] <- p_age_mit2[i] / p_age_mit1[i]
    
  }
  
  # race/eth avg
  for(i in 1:6) {
    p_race_eth_mit1[i] <- inprod(race_eth_mat[,i], p_pop_mit1[1:J]) / sum(race_eth_mat[,i])
    p_race_eth_mit2[i] <- inprod(race_eth_mat[,i], p_pop_mit2[1:J]) / sum(race_eth_mat[,i])
    p_race_eth_mit_rd[i] <- p_race_eth_mit2[i] - p_race_eth_mit1[i]
    p_race_eth_mit_rr[i] <- p_race_eth_mit2[i] / p_race_eth_mit1[i]
  }
  
  # white 
  for(i in 1:2) {
    p_white_mit1[i] <- inprod(white_mat[,i], p_pop_mit1[1:J]) / sum(white_mat[,i])
    p_white_mit2[i] <- inprod(white_mat[,i], p_pop_mit2[1:J]) / sum(white_mat[,i])
    p_white_mit_rd[i] <- p_white_mit2[i] - p_white_mit1[i]
    p_white_mit_rr[i] <- p_white_mit2[i] / p_white_mit1[i]
  }
  
  # hh 
  for(i in 1:6) {
    p_hh_mit1[i] <- inprod(hh_mat[,i], p_pop_mit1[1:J]) / sum(hh_mat[,i])
    p_hh_mit2[i] <- inprod(hh_mat[,i], p_pop_mit2[1:J]) / sum(hh_mat[,i])
    p_hh_mit_rd[i] <- p_hh_mit2[i] - p_hh_mit1[i]
    p_hh_mit_rr[i] <- p_hh_mit2[i] / p_hh_mit1[i]
  }
  
  # edu 
  for(i in 1:2) {
    p_edu_mit1[i] <- inprod(edu_bin_mat[,i], p_pop_mit1[1:J]) / sum(edu_bin_mat[,i])
    p_edu_mit2[i] <- inprod(edu_bin_mat[,i], p_pop_mit2[1:J]) / sum(edu_bin_mat[,i])
    p_edu_mit_rd[i] <- p_edu_mit2[i] - p_edu_mit1[i]
    p_edu_mit_rr[i] <- p_edu_mit2[i] / p_edu_mit1[i]
  }
  
  # income 
  for(i in 1:2) {
    p_inc_mit1[i] <- inprod(income_bin_mat[,i], p_pop_mit1[1:J]) / sum(income_bin_mat[,i])
    p_inc_mit2[i] <- inprod(income_bin_mat[,i], p_pop_mit2[1:J]) / sum(income_bin_mat[,i])
    p_inc_mit_rd[i] <- p_inc_mit2[i] - p_inc_mit1[i]
    p_inc_mit_rr[i] <- p_inc_mit2[i] / p_inc_mit1[i]
    } 

})  
  



