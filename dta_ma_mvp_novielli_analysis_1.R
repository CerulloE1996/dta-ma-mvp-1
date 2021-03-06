# Code Author: Enzo Cerullo
# Date: 06/03/2021
# Relating to manuscript: 
# Meta-analysis of dichotomous and polytomous diagnostic tests without a gold standard
# Authors: Enzo Cerullo,  Hayley E. Jones,  Terry Quinn, Nicola J. Cooper,  Alex J. Sutton
# For any queries contact via email at: enzo.cerullo.bath.edu

require(rstan)
require(bayesplot)
require(dplyr)
require(ggplot2)
require(devtools)
require(posterior)
require(patchwork)
require(loo)
require(bayesSurv)
require(scales)
require(ggrepel)

# INSTALL CmDStan - see https://mc-stan.org/cmdstanr/articles/cmdstanr.html
#devtools::install_github("stan-dev/cmdstanr", force = TRUE)
require(cmdstanr)
#install_cmdstan(overwrite = TRUE)

###################################################################################################
# SET WORKING DIRECTTORY
setwd("/media/enzo/A05C06A35C0673F6/Users/Enzo/Documents/latent variable modelling/latent trait/Wells")


rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

options(scipen = 999)
options(max.print = 1000000000)

################################################################################################
# Wells, D-Dimer and ref test(s) data (11 studies w/ complete data)
num <- 11
#T1 <- c(1,1,2,2,1,3,2,2,1,2,4,2,2) 
T1 <- rep(1, times = 11) # need to see what refs the studies used 
T2 <- rep(2, times = 11)
T <- matrix(c(T1, T2), ncol=2, nrow=11)
### GS=1
#DDIMER=1
r1 <- c(1, 3, 4, 15, 1, 4, 17, 18, 23, 6, 28) # low
r2 <- c(6, 9, 17, 31,22,6, 61, 16, 37, 15, 117) # moderate
r3 <- c(8, 30, 33, 54,35,13,79,21, 9,  10, 109) # high
# 2nd dichot (default ; Novielli et al)
r1_b <- r1 # wells = 0 
r2_b <- r2+r3 #  wells = 1 , GS = 1, DDIMER = 1
# 1st dichot 
#r1_b <- r1 + r2 # wells = 0 
#r2_b <- r3 #  wells = 1 , GS = 1, DDIMER = 1
#DDIMER=0
r4 <- c(0, 1, 1, 1, 0, 0, 3, 0, 5, 0, 1)
r5 <- c(0, 3, 7, 0, 0, 3, 15,1, 7, 3, 0)
r6 <- c(2, 0, 2, 1, 0, 2, 15, 0, 0, 1, 0)
r3_b <- r4
r4_b <- r5+r6 #  wells = 1 , GS = 1, DDIMER = 0
#r3_b <- r4 + r5 # wells = 0 
#r4_b <- r6 #  wells = 1 , GS = 1, DDIMER = 1
### GS=0
#DDIMER=1
r7 <- c(8, 8, 25, 49, 20, 17, 113, 85, 1, 42, 233)
r8 <- c(18, 12,51,51, 23, 9, 93,  83,  6, 59, 104)
r9 <- c(2, 8, 8, 36,  9,  2, 55,  30, 3, 16,  29)
r5_b <- r7
r6_b <- r8+r9 #  wells = 1 , GS = 0, DDIMER = 1
#r5_b <- r7 + r8 # wells = 0 
#r6_b <- r9 #  wells = 1 , GS = 1, DDIMER = 1
 #DDIMER=0
r10 <- c(32, 76, 176, 70, 17, 97, 313, 193, 3, 41, 243)
r11 <- c(20, 43, 113, 54, 19, 48, 23,  89, 5, 46, 16)
r12 <- c(5,  7,  6,  21,  12, 13, 50, 20,  2, 4, 3)
r7_b <- r10
r8_b <- r11+r12 #  wells = 1 , GS = 0, DDIMER = 0
#r7_b <- r10 + r11 # wells = 0 
#r8_b <- r12 #  wells = 1 , GS = 1, DDIMER = 1

#marginalise over d-dimer (ref vs wells) 
r1_t1 <- r1_b + r3_b ; r1_t1
r2_t1 <- r2_b + r4_b ; r2_t1
r3_t1 <- r5_b + r7_b ; r3_t1
r4_t1 <- r6_b + r8_b ; r4_t1

#marginalise over Wells (ref vs d-dimer) -
r1_t2 <- r1_b + r2_b
r2_t2 <- r3_b + r4_b
r3_t2 <- r5_b + r6_b
r4_t2 <- r7_b + r8_b

#marginalise over ref (wells vs d-dimer) -
r1_t3 <- r1_b + r5_b
r2_t3 <- r2_b + r6_b
r3_t3 <- r3_b + r7_b
r4_t3 <- r4_b + r8_b

prev <- (r1_t1+r2_t1)/(r2_t1+ r1_t1 + r3_t1 + r4_t1)
prevs <- round(prev,2)

ns <- c()
for (i in 1:num) {ns[i] <- r1[i] + r2[i] + r3[i] + r4[i] +
                           r5[i] + r6[i] + r7[i] + r8[i] +
                           r9[i] + r10[i] + r11[i] + r12[i] }
# order by test
data <- data.frame(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12, ns, t1 = T[,1], t2= T[,2]) #%>% arrange(t1)

r <- matrix(ncol = 12, nrow = num, c(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12)) ; r

t1 <- matrix(ncol = 4, nrow = num, c(r2_t1, r1_t1, r4_t1, r3_t1)) ;t1 ; sum(t1) # this is fine 

t2 <- matrix(ncol = 4, nrow = num, c(r1_t2, r2_t2, r3_t2, r4_t2)) # this is fine 

t3 <- matrix(ncol = 4, nrow = num, c(r2_t3, r1_t3, r4_t3, r3_t3)) ; sum(t3)
ns <- data$ns
data24 <-list()
data24 <- list( r = r, n = ns, NS= num , T=data$t1, num_ref=4, nt=2)
NS <- 11
sum(ns) # N= 4120

y_list <- list()
y1a <- list(length = max(ns))
y1b <- list(length = max(ns))
y1c <- list(length = max(ns))
pa <-  list(length = max(ns))

max <- max(ns)

for (i in 1:NS) {
  y1a[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(1, r[i,3]), rep(1, r[i,4]), rep(1, r[i,5]), rep(1, r[i,6]), rep(0, r[i,7]), rep(0, r[i,8]), rep(0, r[i,9]), rep(0, r[i,10]), rep(0, r[i,11]), rep(0, r[i,12]), rep(100,  max - ns[i] )) # ref test
  y1b[[i]] = c(rep(1, r[i,1]), rep(1, r[i,2]), rep(1, r[i,3]), rep(0, r[i,4]), rep(0, r[i,5]), rep(0, r[i,6]), rep(1, r[i,7]), rep(1, r[i,8]), rep(1, r[i,9]), rep(0, r[i,10]), rep(0, r[i,11]), rep(0, r[i,12]), rep(100,  max - ns[i] )) # D-Dimer
  y1c[[i]] = c(rep(1, r[i,1]), rep(2, r[i,2]), rep(3, r[i,3]), rep(1, r[i,4]), rep(2, r[i,5]), rep(3, r[i,6]), rep(1, r[i,7]), rep(2, r[i,8]), rep(3, r[i,9]), rep(1, r[i,10]), rep(2, r[i,11]), rep(3, r[i,12]), rep(100,  max - ns[i] ) ) # Wells
   y_list[[i]] =   matrix(ncol = 3, c(y1a[[i]] , y1b[[i]], y1c[[i]] )) 
}

y = array(data= unlist(y_list), dim = c( max(ns), 3, NS))
r = array(data = c(t2, t1, t3), dim = c( NS,4,3)) ; r

r_new <- array(dim = c(3, NS, 4))

for (s in 1:NS) { 
  for (cell in 1:4) { 
    for (i in 1:3) { 
      r_new[i,s,cell] = r[s,cell,i];
    }
  }
}

r_full <- array(data= c(r1,r2,r3,r4,r5, r6, r7, r8, r9, r10, r11, r12), dim = c(11, 12))
r_full
######################################################################################################

n = 11
m <- matrix(data = c(1, 0, 0, 0.80), nrow = 2, ncol = 2)
m2 <- array(data = rep(m,n), dim = c(2,2,n))
m3 <- array(dim = c(n,2,2))

for (s in 1:n) {
  for (i in 1:2) { 
    for (j in 1:2) {
      m3[s,i,j] = m2[i,j,s]
    }
  }
}

init = list(
            a1_m_raw = c(-2, 0.85), 
            a2_m_raw = c(-0.35, 1.60),
            a3_m_raw = c(-0.82, 0.80),
            mu = t(matrix(data = c(0.7,  2,  0.5, 
                                -2, -1, -1), ncol = 2))    ,
            sd1 = t(matrix(data = c(0.5,  0.2,  1, 
                                   0.5, 0.2, 0.2), ncol = 2))    ,
             C_d = t(array(dim = c(2,n), data = rep(c(-0.5, 0.5),2*n ) ))  ,
             C_nd = t(array(dim = c(2,n), data = rep(c(-0.5, 0.5),2*n ) ))  ,
            p = prevs[1:n],
            L_Omega_global_d = m3[1,,],
            L_Omega_global_nd = m3[1,,])

######################################################################################################
# select model to compile
################
# dichotomous Wells  ( imperfect GS )
#################
file <- file.path(file = "mvp_wells_CI_dichotomous.stan") # CI
file <- file.path(file = "mvp_wells_CD_dichotomous.stan") # CD

#################
# polytomous Wells 
#################

### perfect GS, CI
file <- file.path(file = "mvp_wells_CI_perfectGS.stan") 

### perfect GS, CD
file <- file.path(file = "mvp_wells_CD_perfectGS.stan") 

### imperfect GS, CI
file <- file.path(file = "mvp_wells_CI.stan") 

### imperfect GS, CD
file <- file.path(file = "mvp_wells_CD.stan") 

mod <- cmdstan_model(file) # compile model
#############################################################################################
# MCMC (HMC algorithm)
#############################################################################################
n <- 11
studies <- c(1:n)
num <- max(ns[studies])
data =  list( n_studies=n,  ns=ns[studies], y=y[1:num,,studies], 
              num_binary_tests = 2, Thr = c(1,1,2),         
              nt = 3, r = array(r_new[, studies , ], dim = c(3,n,4)), 
              numg = 5000, 
              n_patterns = 12,
              ns_cumsum = cumsum(ns[studies]),
              total_n = sum(ns[studies]),
              ind = c(0, rep(1, n-1)))
 
meta_model2 <- mod$sample(
  data = data,
  seed = 123,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000, 
  refresh = 20, 
  init = list(init,init,init,init), 
  adapt_delta = 0.98, 
  max_treedepth = 9)

meta_model2$cmdstan_diagnose() # view sampler diagnostics

meta_model2r <- rstan::read_stan_csv(meta_model2$output_files())  # convert to rstan CSV

print(meta_model2r, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                            "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "cov_global_d","cov_global_nd","p"
                         #   "alpha_d", "alpha_nd", "C_dm",
                            #"C_dm2",
                          #  "C_ndm", "C_ndm2",  
                        #    "sd1", "mu" , "log_diff",
                         #   "Se_pred", "Sp_pred", "Wells_DDimer_BTP_Se_pred", "Wells_DDimer_BTP_Sp_pred",
                      #    "Wells_DDimer_BTN_Se_pred", "Wells_DDimer_BTN_Sp_pred"
                         #   "L", "M", "H", "p_dm", "mu",
                       #   "mu_scale_log", "sd_scale"
),
probs = c(0.025,0.5, 0.975)) 

# save output 
saveRDS(meta_model2r, file = "Wells_CD_xu_2000_ad0.99_norm10.rds")

print(meta_model2r, pars= c("L_Omega_d", "L_Omega_nd", "alpha_d", "alpha_nd"),probs = c(0.025,0.5, 0.975))

#############################################################################################
# observed - expected correlation plots 
#############################################################################################

dc <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,5],3) ; dc
dc_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,4],3) ; dc_l
dc_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dc"))$summary[,6],3) ; dc_u

dc_data <- tibble(dc, dc_l, dc_u, `test-pair` = as.factor(c(rep("Ref vs D-Dimer", n), rep("Ref vs Wells",n),rep("D-Dimer vs Wells",n))), 
                  obs = seq(1, length(dc), by = 1))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("wells_ppc_corr_residuals_scale.tif",units = "in", width = 10, height=4, res=800, compression = "lzw")

ggplot(data = dc_data, aes(y = dc, x=obs, colour = `test-pair`)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=dc_l, ymax=dc_u), width= 0.75, position=position_dodge(.9)) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  ylim(-0.31, 0.31) + 
  ylab("Correlation Residuals") + 
  xlab("Study / test-pair") + 
  theme(text = element_text(size=14),
        axis.text.x = element_text()) 

#dev.off()

#############################################################################################
# observed - expected table count plots 
#############################################################################################
dt <-   round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,5],3) ; dt
dt_l <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,4],3) ; dt_l
dt_u <- round(summary(meta_model2r, probs = c(0.025,  0.5, 0.975), pars = c("dt"))$summary[,6],3) ; dt_u

dt_data <- tibble(dt, dt_l, dt_u, `test-pair` = as.factor(c(rep("Ref vs D-Dimer", n*4), rep("Ref vs Wells",n*4),
                                                           rep("D-Dimer vs Wells",n*4))), 
                                                           obs = seq(1, length(dt), by = 1))

cutoff <- data.frame( x = c(-Inf, Inf), y = 0, cutoff = factor(0) )

#tiff("wells_ppc_corr_tables_scale.tif",units = "in", width = 8, height=6, res=500, compression = "lzw")

    ggplot(data = dt_data, aes(y = dt, x=obs, colour = `test-pair`)) + 
      geom_point(size = 1) + 
      geom_errorbar(aes(ymin=dt_l, ymax=dt_u), width= 0.75, position=position_dodge(.9)) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      ylim(-50, 50) + 
      ylab("2x2 table residuals") + 
      xlab("Study / table cell") + 
      theme(text = element_text(size=14),
           # axis.text.x = element_text(),
            axis.text.x=element_blank()) +
      theme(legend.position="none") +
      facet_wrap( ~ `test-pair`, scales = "free", dir = "v")

#dev.off()
#############################################################################################
# Trace and posterior density plots
#############################################################################################
stan_trace(meta_model2r, pars = c("sd1"))
stan_trace(meta_model2r, pars = c("p"))
stan_trace(meta_model2r, pars = c("C_dm", "C_ndm"))
stan_trace(meta_model2r, pars = c("Se", "Sp"))
stan_trace(meta_model2r, pars = c("se", "sp"))
stan_trace(meta_model2r, pars = c("alpha"))
stan_trace(meta_model2r, pars = c("alpha_d", "alpha_nd"))

stan_trace(meta_model2r, pars = c("sd_scale", "nu_scales_log"))
stan_trace(meta_model2r, pars = c("nu"))

stan_trace(meta_model2r, pars = c("a1_m_raw", "a2_m_raw", "a3_m_raw"))
stan_trace(meta_model2r, pars = c("L_Omega_d"))
stan_trace(meta_model2r, pars = c("L_Omega_nd"))
stan_trace(meta_model2r, pars = c("mu"))

stan_dens(meta_model2r, pars =  c("p"), separate_chains = TRUE)
stan_dens(meta_model2r, pars =  c("p"))
stan_dens(meta_model2r, pars = c("L_Omega_d"))
stan_dens(meta_model2r, pars = c("L_Omega_nd"))
stan_dens(meta_model2r, pars =  c( "nu"))
stan_dens(meta_model2r, pars =  c( "sd", "z1"))
stan_dens(meta_model2r, pars =  c("Se", "Sp"))
stan_dens(meta_model2r, pars =  c("se", "sp"))
stan_dens(meta_model2r, pars =  c("L", "M", "H"))

#############################################################################################
# LOO
#############################################################################################

log_lik <- extract_log_lik(meta_model2r, parameter_name = "log_lik")
str(log_lik)

draws <-  dim(log_lik)[1]/4
length <- length(log_lik[1,][log_lik[1,] != "NaN" & !is.na(log_lik[1,]) ])
loglik2 <- array(data = log_lik[log_lik != "NaN" & !is.na(log_lik) ], dim = c(draws*4, length ))

r_eff <- relative_eff(exp(loglik2), cores = 4, chain_id = c(rep(1, draws), rep(2, draws), rep(3, draws), rep(4,draws)))

mod_loo <- loo(loglik2, cores = 4, r_eff = r_eff)

mod_loo


##############################################################################
## Figure for summary of results from 4 models for section 4.2 
##############################################################################

m1 <- readRDS(file = "Wells_CI_perfectGS_xu_1000_norm10.rds")   ###  M1 - perfect GS, CI
m2 <- readRDS(file = "Wells_CD_perfectGS_xu_1000_ad0.95_norm10.rds")   ###  M2 - perfect GS, CD
m3 <- readRDS(file = "Wells_CI_xu_1000_norm10.rds")   ###  M3 - IGS, CI
m4 <- readRDS(file = "Wells_CD_xu_1000_ad0.95_norm10.rds")   ###  M4 - IGS, CD

models <- list(m1, m2, m3, m4)

print(m4, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                            "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "cov_global_d","cov_global_nd","p",
                           # "alpha", 
                            "C_dm", "sd1", "mu", 
                          #  "mu_location", "sd_location", "mu_scale_log", "sd_scale",
                            "Se_pred", "Sp_pred", "Wells_DDimer_BTP_Se_pred", "Wells_DDimer_BTP_Sp_pred","Wells_DDimer_BTN_Se_pred", "Wells_DDimer_BTN_Sp_pred"
),
probs = c(0.025,0.5, 0.975))
#############################################################################################
# LOO TABLE 
#meta_model2r <- readRDS("Wells_fixed_thr_CI.rds") 

log_lik <- list()
loglik2 <- list()
draws <- list()
length <- list()
r_eff <- list()
mod_loo <- list()

require(loo)
for (i in 1:length(models)) { 
  log_lik[[i]] <- extract_log_lik(models[[i]], parameter_name = "log_lik")
  
  draws[[i]] <- dim(log_lik[[i]])[1]/4
  length[[i]] <- length(log_lik[[i]][1,][log_lik[[i]][1,] != "NaN" & !is.na(log_lik[[i]][1,]) ])
  loglik2[[i]] <- array(data = log_lik[[i]][ log_lik[[i]] != "NaN" & !is.na(log_lik[[i]] ) ] , dim = c(draws[[i]]*4, length[[i]] ))
  
  
  r_eff[[i]] <- relative_eff(exp(loglik2[[i]]), cores = 1, chain_id = c(rep(1, draws[[i]]), rep(2, draws[[i]]), rep(3, draws[[i]]), rep(4,draws[[i]])))
  
  mod_loo[[i]] <- loo(loglik2[[i]], cores = 1, r_eff = r_eff[[i]])
} 

mod_loo

loo_compare(mod_loo)

###############################################################

m_mod1 <- list(); m2_mod1 <- list()
l_mod1 <- list(); l2_mod1 <- list()
u_mod1 <- list(); u2_mod1 <- list()
m_mod1_cat <- list(); m2_mod1_cat <- list()
l_mod1_cat <- list(); l2_mod1_cat <- list()
u_mod1_cat <- list(); u2_mod1_cat <- list()
data_Se_mod <- list() 
data_Sp_mod <- list() 
data_d_mod1 <- list() 
data_nd_mod1 <- list() 

for (i in 1:length(models)) {
    m_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Se","Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,5],2)
    
    m2_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,5],2)
    
    
    l_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                            pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,4],2)  
    
    l2_mod1[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                             pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,4],2)
    
    
    u_mod1[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                             pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,6],2)  
    
    u2_mod1[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), 
                              pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,6],2)
    
    
    # wells categories
    m_mod1_cat[[i]] <- c(round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,5],2))
    m2_mod1_cat[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[2]", "M[2]", "H[2]"))$summary[,5],2)
    
    
    l_mod1_cat[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,4],2)  
    l2_mod1_cat[[i]] <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975),pars = c("L[2]", "M[2]", "H[2]"))$summary[,4],2)
    
    
    u_mod1_cat[[i]] <-  round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[1]", "M[1]", "H[1]"))$summary[,6],2)  
    u2_mod1_cat[[i]]  <- round(summary(models[[i]], probs = c(0.025,  0.5, 0.975), pars = c("L[2]", "M[2]", "H[2]"))$summary[,6],2)
    
    
    data_Se_mod[[i]] <- data.frame(m=m_mod1[[i]],l=l_mod1[[i]],u=u_mod1[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                            3.9+0.4*(i-1), 
                                                                                            6.9+0.4*(i-1), 
                                                                                            9.9+0.4*(i-1), 
                                                                                            12.9+0.4*(i-1)), 
                               label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                 "Wells & D-Dimer, BTN",
                                                 "Wells & D-Dimer, BTP"), 
                                               levels = c("Wells & D-Dimer, BTP",
                                                          "Wells & D-Dimer, BTN",
                                                          "Wells", 
                                                          "D-Dimer",
                                                          "Referemce")),
                               Model = factor(rep(i,  5)))
    
    data_Sp_mod[[i]] <- data.frame(m=m2_mod1[[i]],l=l2_mod1[[i]],u=u2_mod1[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                               3.9+0.4*(i-1), 
                                                                                               6.9+0.4*(i-1), 
                                                                                               9.9+0.4*(i-1), 
                                                                                               12.9+0.4*(i-1)),
                               label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                 "Wells & D-Dimer, BTN",
                                                 "Wells & D-Dimer, BTP"), 
                                               levels = c("Wells & D-Dimer, BTP",
                                                          "Wells & D-Dimer, BTN",
                                                          "Wells", 
                                                          "D-Dimer",
                                                          "Referemce")),
                               Model = factor(rep(i, 5)))
    
    ####### wells categories
    ## mod 1
    data_d_mod1[[i]] <- data.frame(m=m_mod1_cat[[i]], l=l_mod1_cat[[i]], u=u_mod1_cat[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                                          3.9+0.4*(i-1), 
                                                                                                          6.9+0.4*(i-1)), 
                              Category  = factor(c("Low", "Medium", "High"),
                                                 levels = c("Low", "Medium", "High")),
                              Model = factor(rep(i, 3)))
    
    data_nd_mod1[[i]] <- data.frame(m=m2_mod1_cat[[i]],l=l2_mod1_cat[[i]],u=u2_mod1_cat[[i]], location =  c(0.9+0.4*(i-1), 
                                                                                                            3.9+0.4*(i-1), 
                                                                                                            6.9+0.4*(i-1)),
                               Category  = factor(c("Low", "Medium", "High"),
                                                  levels = c("Low", "Medium", "High")),
                               Model = factor(rep(i, 3)))
    
    
}

require(data.table)
# use rbindlist to put all models in same data frame
data_Se <- rbindlist(data_Se_mod)
data_Sp <- rbindlist(data_Sp_mod)
data_d <- rbindlist(data_d_mod1)
data_nd <- rbindlist(data_nd_mod1)

data_Se$Model<- factor(data_Se$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))

data_Se$Mod = c(rep("Perfect GS, CI", 5), 
                rep("Perfect GS, CD", 5),
                rep("Imperfect GS, CI", 5),
                rep("Imperfect GS, CD", 5))
                
  
data_Sp$Model<- factor(data_Sp$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))

data_Sp$Mod = c(rep("Perfect GS, CI", 5),
                rep("Perfect GS, CD", 5),
                rep("Imperfect GS, CI", 5),
                rep("Imperfect GS, CD", 5))


data_d$Model<- factor(data_d$Model, labels =c("M1: Perfect GS, CI", 
                                              "M2: Perfect GS, CD",
                                              "M3: Imperfect GS, CI", 
                                              "M4: Imperfect GS, CD"))


data_nd$Model<- factor(data_nd$Model, labels =c("M1: Perfect GS, CI", 
                                                "M2: Perfect GS, CD",
                                                "M3: Imperfect GS, CI", 
                                                "M4: Imperfect GS, CD"))


################# plot 

cat_d_plot <- ggplot(data_d, aes(x=location, y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=Category)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Diseased") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=15), 
        axis.text.x = element_text(size = 10)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0, 1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = -0.2, vjust = -0.25), alpha = 0.35, size = 4.5) 

cat_d_plot

cat_nd_plot <- ggplot(data_nd, aes(x=location, y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=Category)) + 
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Non-diseased") + 
  theme(axis.ticks = element_blank()) + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16), 
        axis.text.x = element_text(size = 10))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0, 1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = -0.2, vjust = -0.25), alpha = 0.35, size = 4.5) +
  theme(legend.title = element_text(size = 0), 
        legend.text  = element_text(size = 12))

cat_nd_plot

require(patchwork)
cat_d_plot + cat_nd_plot

####################
se_plot <- ggplot(tibble(data_Se), aes(x=as.numeric(location), y = m, ymin= l,ymax= u, shape = Model)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Sensitivity") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16), 
        axis.text.x = element_text(size = 10))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0,1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 4) 

se_plot


sp_plot <- ggplot(data_Sp, aes(x=location, y = m, ymin= l,ymax= u, shape = Model )) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=(label))) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Specificity")  + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16), 
        axis.text.x = element_text(size = 10))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10), limits = c(0, 1)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35, size = 4) +
  theme(legend.title = element_text(size = 0), 
        legend.text  = element_text(size = 12))


sp_plot

tiff("wells_figure_summary.tif",units = "in", width = 9, height=7, res=500, compression = "lzw")
se_plot + sp_plot
dev.off()

tiff("wells_figure_cats.tif",units = "in", width = 9, height=5, res=500, compression = "lzw")
cat_d_plot + cat_nd_plot
dev.off()

stan_plot(models[[1]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[1]], pars = c("p"))

stan_plot(models[[2]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[2]], pars = c("p"))

stan_plot(models[[3]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[3]], pars = c("p"))

stan_plot(models[[4]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[4]], pars = c("p"))

stan_plot(models[[5]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.80, outer_level = 0.95)
stan_trace(models[[5]], pars = c("p"))


##############################
####### ROC plot
########################################
# load in the model
mod <- models[[4]] # M4


#tiff("figure_4.tif",units = "in", width = 11, height=7, res=500, compression = "lzw")

       
print(mod, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                  "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "rho_global_d","rho_global_nd",
                  "Se_pred", "Sp_pred"),
      probs = c(0.025,0.5, 0.975))

pnorm2 <- function(x) { plogis( 1.702*x ) }
qnorm2 <- function(x) { (1/1.702)*qlogis(x) }

#####
## credible region
cred_1 <- list()
cred_1p <- list()

for (t in 1:3) { 
  cred_1[[t]] <- tibble(y = qnorm2(extract(mod, pars = "Se")$Se[,t]), x = qnorm2(extract(mod, pars = "Sp")$Sp[,t]))
  cred_1p[[t]] <- tibble(y = (extract(mod, pars = "Se")$Se[,t]), x =  (extract(mod, pars = "Sp")$Sp[,t]))
} 

cred_1[[4]] <- tibble(y = qnorm2(extract(mod, pars = "Wells_DDimer_BTN_Se")$Wells_DDimer_BTN_Se), x = qnorm2(extract(mod, pars = "Wells_DDimer_BTN_Sp")$Wells_DDimer_BTN_Sp))
cred_1p[[4]] <- tibble(y = (extract(mod, pars = "Wells_DDimer_BTN_Se")$Wells_DDimer_BTN_Se), x =  (extract(mod, pars = "Wells_DDimer_BTN_Sp")$Wells_DDimer_BTN_Sp))

cred_1[[5]] <- tibble(y = qnorm2(extract(mod, pars = "Wells_DDimer_BTP_Se")$Wells_DDimer_BTP_Se), x = qnorm2(extract(mod, pars = "Wells_DDimer_BTP_Sp")$Wells_DDimer_BTP_Sp))
cred_1p[[5]] <- tibble(y = (extract(mod, pars = "Wells_DDimer_BTP_Se")$Wells_DDimer_BTP_Se), x =  (extract(mod, pars = "Wells_DDimer_BTP_Sp")$Wells_DDimer_BTP_Sp))


require(data.table)
cred <- rbindlist(cred_1, idcol = TRUE)
cred_p <- rbindlist(cred_1p, idcol = TRUE)

cred2 <- mutate(cred,  Test = factor(.id, label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells & D-Dimer, BTP"), 
                                                          levels = c("Wells & D-Dimer, BTP",
                                                                     "Wells & D-Dimer, BTN",
                                                                     "Wells", 
                                                                     "D-Dimer",
                                                                     "Referemce"))))

# in inv_probit space
g <- ggplot(data = cred2, aes(x = x, y = y, colour = Test))  + 
 # geom_point() + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
View(pb$data[[2]])

el = pb$data[[1]][c("x","y", "group")]

  
credible_region <- tibble(x = pnorm2(el$x), y = pnorm2(el$y), Test = factor(el$group, label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                                                                        "Wells & D-Dimer, BTN",
                                                                                                        "Wells & D-Dimer, BTP"), 
                                                                                                      levels = c("Wells & D-Dimer, BTP",
                                                                                                                 "Wells & D-Dimer, BTN",
                                                                                                                 "Wells", 
                                                                                                                 "D-Dimer",
                                                                                                                 "Referemce"))))
credible_region


####
## prediction region

pred_1 <- list()
pred_1p <- list()

for (t in 1:3) { 
  pred_1[[t]] <- tibble(y = qnorm2(extract(mod, pars = "Se_pred")$Se_pred[,t]), x = qnorm2(extract(mod, pars = "Sp_pred")$Sp_pred[,t]))
  pred_1p[[t]] <- tibble(y = (extract(mod, pars = "Se_pred")$Se_pred[,t]), x =  (extract(mod, pars = "Sp_pred")$Sp_pred[,t]))
} 

pred_1[[4]] <- tibble(y = qnorm2(extract(mod, pars = "Wells_DDimer_BTN_Se_pred")$Wells_DDimer_BTN_Se_pred), x = qnorm2(extract(mod, pars = "Wells_DDimer_BTN_Sp_pred")$Wells_DDimer_BTN_Sp_pred))
pred_1p[[4]] <- tibble(y = (extract(mod, pars = "Wells_DDimer_BTN_Se_pred")$Wells_DDimer_BTN_Se_pred), x =  (extract(mod, pars = "Wells_DDimer_BTN_Sp_pred")$Wells_DDimer_BTN_Sp_pred))

pred_1[[5]] <- tibble(y = qnorm2(extract(mod, pars = "Wells_DDimer_BTP_Se_pred")$Wells_DDimer_BTP_Se_pred), x = qnorm2(extract(mod, pars = "Wells_DDimer_BTP_Sp_pred")$Wells_DDimer_BTP_Sp_pred))
pred_1p[[5]] <- tibble(y = (extract(mod, pars = "Wells_DDimer_BTP_Se_pred")$Wells_DDimer_BTP_Se_pred), x =  (extract(mod, pars = "Wells_DDimer_BTP_Sp_pred")$Wells_DDimer_BTP_Sp_pred))


require(data.table)
pred <- rbindlist(pred_1, idcol = TRUE)
pred_p <- rbindlist(pred_1p, idcol = TRUE)

pred2 <- mutate(pred,  Test = factor(.id, label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells & D-Dimer, BTP"), 
                                                          levels = c("Wells & D-Dimer, BTP",
                                                                     "Wells & D-Dimer, BTN",
                                                                     "Wells", 
                                                                     "D-Dimer",
                                                                     "Referemce"))))


g <- ggplot(data = pred2, aes(x = x, y = y, colour = Test))  + 
  # geom_point() + 
  stat_ellipse()  

g

# Get ellipse coordinates from plot
pb <-  ggplot_build(g)
View(pb$data[[2]])

el = pb$data[[1]][c("x","y", "group")]


pred_region <- tibble(x = pnorm2(el$x), y = pnorm2(el$y), Test = factor(el$group, label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                                                                      "Wells & D-Dimer, BTN",
                                                                                                      "Wells & D-Dimer, BTP"), 
                                                                                                    levels = c("Wells & D-Dimer, BTP",
                                                                                                               "Wells & D-Dimer, BTN",
                                                                                                               "Wells", 
                                                                                                               "D-Dimer",
                                                                                                               "Referemce"))))
pred_region


## medians
print(mod, pars= c("Se", "Sp", "Wells_DDimer_BTN_Se", "Wells_DDimer_BTN_Sp",
                   "Wells_DDimer_BTP_Se", "Wells_DDimer_BTP_Sp", "rho_global_d","rho_global_nd",
                   "Se_pred", "Sp_pred"),
      probs = c(0.025,0.5, 0.975))




median_sens <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Se"))$summary[,5], 2), 
                    round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Wells_DDimer_BTN_Se", "Wells_DDimer_BTP_Se"))$summary[,5], 2))
median_spec <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Sp"))$summary[,5], 2), 
                 round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("Wells_DDimer_BTN_Sp", "Wells_DDimer_BTP_Sp"))$summary[,5], 2))

medians <- tibble(median_sens = median_sens, median_spec = median_spec, Test = factor( c(1:5), label  = factor(c("Referemce", "D-Dimer", "Wells", 
                                                                                                         "Wells & D-Dimer, BTN",
                                                                                                         "Wells & D-Dimer, BTP"), 
                                                                                                       levels = c("Wells & D-Dimer, BTP",
                                                                                                                  "Wells & D-Dimer, BTN",
                                                                                                                  "Wells", 
                                                                                                                  "D-Dimer",
                                                                                                                  "Referemce"))))

print(mod, pars= c("se"),   probs = c(0.025,0.5, 0.975))

s_sens <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("se"))$summary[,5], 5))
s_spec <- c(round(summary(mod, probs = c(0.025,  0.5, 0.975), pars = c("sp"))$summary[,5], 5))

ss <- tibble(s_sens = s_sens, s_spec = s_spec, Test = factor( rep(c(1:3), 11), label  = factor(c("Referemce", "D-Dimer", "Wells"),
                                                                                                               
                                                                                                               levels = c(
                                                                                                                          "Wells", 
                                                                                                                          "D-Dimer",
                                                                                                                                                                                                                                           "Reference"))))

#############################
## plot

g <- ggplot(data = medians, aes(y=median_sens, x = 1 - median_spec, colour = Test)) +    # summary points
  geom_point( size=2 ) + 
#  geom_point(data = ss, aes(y=s_sens, x = 1 - s_spec, colour = Test))  + 
 # xlim(0, 1) + 
  #ylim(0, 1) + 
  theme_bw() + 
  scale_x_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  + 
  scale_y_continuous(breaks = seq(0,1,0.1), limits = c(0,1))  +
  theme(legend.title=element_blank()) + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  geom_polygon(data = credible_region, aes(x = 1  - x, y = y, colour = Test), alpha=0.05, size=0.4)  + 
  geom_path(data = pred_region, aes(x = 1  - x, y = y, colour = Test), linetype = 2, size=0.4) + 
  theme(legend.position =  "bottom")
g


tiff("wells_sroc.tif",units = "in", width = 7, height=5, res=500, compression = "lzw")
g
dev.off()


##############################################################################
## Figure for summary of results from 4 models for section 4.1 (dichotomous Wells)
##############################################################################

## first dichot 
d1_ci <- readRDS(file = "Wells_dichotomous_CI_1st_dichot.rds")   ###  1st dichot, M1 - perfect GS, CI, random thresholds , diff SDs
d2_ci <- readRDS(file = "Wells_dichotomous_CI_2nd_dichot.rds")   ###  M2 - IGS, CI, fixed thresholds, diff SDs

d1_cd <- readRDS(file = "Wells_CD_dichot_1st.rds")   ###  M3 - IGS, CI, random thresholds, diff SDs
d2_cd <- readRDS(file = "Wells_CD_dichot_2nd.rds")   ###  M2 - IGS, CD, fixed thresholds, diff SDs

models2 <- list(d1_ci, d2_ci, d1_cd, d2_cd)

###############################################################

m_mod1 <- list(); m2_mod1 <- list()
l_mod1 <- list(); l2_mod1 <- list()
u_mod1 <- list(); u2_mod1 <- list()
m_mod1_cat <- list(); m2_mod1_cat <- list()
l_mod1_cat <- list(); l2_mod1_cat <- list()
u_mod1_cat <- list(); u2_mod1_cat <- list()
data_Se_mod <- list() 
data_Sp_mod <- list() 
data_d_mod1 <- list() 
data_nd_mod1 <- list() 

for (i in 1:length(models2)) {
  m_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("Se","Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,5],2)
  
  m2_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,5],2)
  
  
  l_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                               pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,4],2)  
  
  l2_mod1[[i]] <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,4],2)
  
  
  u_mod1[[i]] <-  round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                pars = c("Se", "Wells_DDimer_BTN_Se","Wells_DDimer_BTP_Se"))$summary[,6],2)  
  
  u2_mod1[[i]]  <- round(summary(models2[[i]], probs = c(0.025,  0.5, 0.975), 
                                 pars = c("Sp","Wells_DDimer_BTN_Sp","Wells_DDimer_BTP_Sp"))$summary[,6],2)
  
  data_Se_mod[[i]] <- data.frame(m=m_mod1[[i]],l=l_mod1[[i]],u=u_mod1[[i]], location =  c(0.9+0.3*(i-1), 
                                                                                          2.9+0.3*(i-1),
                                                                                          4.9+0.3*(i-1), 
                                                                                          6.9+0.3*(i-1), 
                                                                                          8.9+0.3*(i-1)), 
                                 label  = factor(c("Ultrasound", "D-Dimer", "Wells", 
                                                   "Wells & D-Dimer, BTN",
                                                   "Wells & D-Dimer, BTP"), 
                                                 levels = c("Wells & D-Dimer, BTP",
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells", 
                                                            "D-Dimer",
                                                            "Ultrasound")),
                                 Model = factor(rep(i,  5)))
  
  data_Sp_mod[[i]] <- data.frame(m=m2_mod1[[i]],l=l2_mod1[[i]],u=u2_mod1[[i]], location =  c(0.9+0.3*(i-1), 
                                                                                             2.9+0.3*(i-1),
                                                                                             4.9+0.3*(i-1),
                                                                                             6.9+0.3*(i-1),
                                                                                             8.9+0.3*(i-1)),
                                 label  = factor(c("Ultrasound", "D-Dimer", "Wells", 
                                                   "Wells & D-Dimer, BTN",
                                                   "Wells & D-Dimer, BTP"), 
                                                 levels = c("Wells & D-Dimer, BTP",
                                                            "Wells & D-Dimer, BTN",
                                                            "Wells", 
                                                            "D-Dimer",
                                                            "Ultrasound")),
                                 Model = factor(rep(i, 5)))
}

require(data.table)
# use rbindlist to put all models in same data frame
data_Se2 <- rbindlist(data_Se_mod)
data_Sp2 <- rbindlist(data_Sp_mod)

data_Se2$Dichot = c(rep("1st dichotomisation", 5),
                    rep("2nd dichotomisation", 5),
                    rep("1st dichotomisation", 5),
                    rep("2nd dichotomisation", 5))

data_Se2$Mod   =  c(rep("CI", 10),
                    rep("CD", 10))

data_Se2$Model<- factor(data_Se2$Model, labels =c("1st dichotomisation, CI", 
                                                "2nd dichotomisation, CI", 
                                                "1st dichotomisation, CD", 
                                                "2nd dichotomisation, CD"))

data_Sp2$Model<- factor(data_Sp2$Model, labels =c("1st dichotomisation, CI", 
                                                "2nd dichotomisation, CI", 
                                                "1st dichotomisation, CD", 
                                                "2nd dichotomisation, CD"))

data_Sp2$Dichot = c(rep("1st dichotomisation", 5),
                    rep("2nd dichotomisation", 5),
                    rep("1st dichotomisation", 5),
                    rep("2nd dichotomisation", 5))

data_Sp2$Mod   =  c(rep("CI", 10),
                    rep("CD", 10))

################# plot 

se_plot <- ggplot(tibble(data_Se2), aes(x=(location), y = m, ymin= l,ymax= u, shape = Mod)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label, linetype = Dichot)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Sensitivity") + 
  theme(legend.position = "none") + 
  theme(axis.ticks = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10)) + 
  scale_x_discrete( labels = rep("  ", length(data_Se2$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25),alpha = 0.35) 
#  scale_x_discrete( labels = c(" ", " ", " ", "Ref", 
#                               " ", " ", " ", "D-Dimer",
 #                              " ", " ", " ", "Wells", 
  #                             " ", " ", " ", "Wells & D-Dimer, BTN",
   #                            " ", " ", " ", "Wells & D-Dimer, BTP"))

se_plot


sp_plot <- ggplot(data_Sp2, aes(x=location, y = m, ymin= l,ymax= u, shape = Mod)) + 
  geom_point(size=4, alpha=0.2) + 
  geom_pointrange(aes(colour=label, linetype = Dichot)) +
  coord_flip() + 
  theme_bw()  + 
  xlab("") + 
  ylab("Specificity")  + 
  theme(legend.title = element_blank()) + 
  theme(text = element_text(size=16))+ 
  scale_y_continuous(breaks = seq(0, 1, 0.10))+ 
  scale_x_discrete( labels = rep("  ", length(data_Se2$location))) + 
  geom_text(aes( label = paste0( m*100," ", "[",l*100,",", u*100 , "]") ,hjust = 1.1, vjust = -0.25), alpha = 0.35) 


sp_plot


tiff("wells_figure_dichot_summary.tif",units = "in", width = 11, height=7, res=500, compression = "lzw")
se_plot + sp_plot
dev.off()


#####################
### plot of disease prevalences for 1st dichot vs 2nd dichot (CI on left, CD on right)



g1 <- stan_plot(models2[[1]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g1

g2<- stan_plot(models2[[2]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g2

params <- extract(models2[[1]])



g1 + g2

g1 <- stan_plot(models2[[3]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g1

g2<- stan_plot(models2[[4]], pars = c("p"), ncol = 1, show_density = T, ci_level = 0.95, outer_level = 0.95)
g2

g1 + g2





