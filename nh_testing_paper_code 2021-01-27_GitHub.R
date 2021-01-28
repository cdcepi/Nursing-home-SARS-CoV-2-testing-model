#Code for modeling of nursing home COVID-19 testing strategies last updated 12/17/20

options(stringsAsFactors=FALSE) # Personal preference; may not be necessary
setwd("./")

library(mgcv)
library(ggplot2)

# From https://github.com/cdcepi/nCov-2019/blob/master/Controllability%20analysis/Testing_travelers_function.R
# Best guess for antigen test sensitivity
gen_p_test_pos_default <- function(max_sensitivity = 0.99, shape = 7.2, scale = 0.80) {
  d_infectious <- function(t) dgamma(t, shape=shape, scale=scale)
  scale_sensitivity <- max_sensitivity / max(d_infectious(seq(0, 20, by=0.1))) 
  function(t) scale_sensitivity * d_infectious(t) # Sensitivity proportional to infectiousness
}

load("Clifford_RTPCR.RData") # Loads p_test_pos_Clifford (and other objects)

# Load and define infectiousness profiles
load("Clifford_infectiousness.RData")
load("Goyal_infectiousness.RData")
He_pars = readRDS("He_parameters.rds")
d_infectious_He = function(tt) with(He_pars, dgamma(tt, inf.par[1], inf.par[2]))
d_infectious_default = function(t) dgamma(t, shape = 7.2, scale = 0.80)
d_inf = list(Default=d_infectious_default,
             He=d_infectious_He,
             Goyal=function(tt) ifelse(tt>60, 0, d_infectious_ori_Goyal(tt)),
             Clifford=function(tt) ifelse(tt>60, 0, d_infectious_Clifford(tt)))
test_sensitivity = list(p_test_pos_Clifford, gen_p_test_pos_default(max_sensitivity=0.8))
saveRDS(d_inf, "Four_infectiousness_profiles.rds")

#################################################################

d_inf = readRDS("Four_infectiousness_profiles.rds")

# CDC parameters (5 scenarios, #5 is current best estimate [https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html]
CDCscenario = 5
f_presymp = c(0.3, 0.7, 0.3, 0.7, 0.5)[CDCscenario] # Proportion of infections from presymptomatic phase
f_asymp = c(0.1, 0.7, 0.1, 0.7, 0.4)[CDCscenario] # Proportion of infections that are asymptomatic
inf_asymp_over_symp = c(0.25, 1, 0.25, 1, 0.75)[CDCscenario]
exp_to_onset = 6 # days, from exposure to symptom onset

# Relative FOI (force of infection) from symptomatic infections (after symptom onset date)
f_symp = (1 - f_asymp)/(1 - f_asymp + f_asymp * inf_asymp_over_symp)

# Consistent estimate of other parameters
inf_max_vector = sapply(d_inf, function(ff) optimize(ff, c(0, 30), maximum=TRUE)$objective)
inf_to_onset_vector = sapply(d_inf, function(ff) uniroot(function(tt) integrate(ff, 0, tt)$value - f_presymp, c(0, 30))$root) # Formerly inf.par[3]


symp_eff = 1 # Proportion of symptomatic people you find without testing
eeee=.Machine$double.eps^0.25 #default rel.tol for integrate function
int_es=c(eeee,eeee,eeee,eeee*10) # increase tol for error for Clifford profile

# Pick profile
profile = 4
inf_profile = d_inf[[profile]]
inf_max = inf_max_vector[profile]
inf_to_onset = inf_to_onset_vector[profile]
int_error=int_es[profile]

# Sensitivity as function of result_lag
# r_sens is relative sensitivity compared to RT_PCR
sensitivity2 = function(tt, result_lag,r_sens=1) {
  sens_max = (result_lag <= 0.5) * 0.85 + (result_lag >= 4) * 0.95 + (result_lag>0.5 & result_lag<4) * (0.85 + 0.1/3.5 * (result_lag-0.5))
  r_sens * (sens_max/inf_max) * inf_profile(tt)
}

## If tested at t0 since start of infectiousness, and result/exclusion at t1=t0+result_lag, 
##   fraction of transmission prevented is int_t1^infinity sensitivity(t0) * inf_asymp(tt) dtt
inf_asymp = function(tt) inf_profile(tt) * ifelse(tt>inf_to_onset, 1-f_symp, 1)

transmission_symp = Vectorize(function(t_I) {
  f_symp * integrate(inf_profile, inf_to_onset, inf_to_onset+t_I)$value
})

#old, gave errors related to precision/rounding
#transmission_reduction2 = Vectorize(function(t0, result_lag, t_I,r_sens) {
#  sensitivity2(t0, result_lag,r_sens) * integrate(inf_asymp, lower=t0+result_lag, upper=t0+result_lag+t_I
#                                                  ,stop.on.error=FALSE
                                                  #,rel.tol=int_error
#                                                  )$value
#})

#Alternate hard-coded approximation of integrand above, avoids roundoff errors
transmission_reduction2 = Vectorize(function(t0, result_lag, t_I,r_sens) {
  sensitivity2(t0, result_lag,r_sens) * sum(inf_asymp(seq(t0+result_lag, min(100, t0+result_lag+t_I), 0.1)+0.05)) * 0.1
})


PI2 = Vectorize(function(tt, TT, tl, t_I, tmax=30,r_sens=1) {
  if (tt > tmax) {res = 0} else {
    res = transmission_reduction2(tt, tl, t_I,r_sens) + (1-sensitivity2(tt, tl,r_sens)) * PI2(tt+TT, TT, tl, t_I, tmax,r_sens)
  }
  res
})

mean_transmission_reduction2 = function(test_period, result_lag, t_isol,r_sens) {
  integrate(function(x) PI2(x, test_period, result_lag, t_isol,30,r_sens), lower=0, upper=test_period
            ,rel.tol=eeee*100 #Needed to avoid errors
            )$value/test_period
}

####
# code to estimate outbreak size and duration
#

# Branching process in finite populations

final.size.dist = function(init.state, R) {
  # init.state is vector (length population size + 1); (k+1)-th element = probability of k initial cases
  # R is basic reproduction (expected number of secondaries with one index and rest susceptible)
  N = length(init.state) - 1
  prob = R/(N-1)
  genmat = diag(init.state * 0)
  diag(genmat[1+0:N, 1+N-0:N]) = init.state
  for (k in 1:(N-1)) for (j in 1:k) {
    diag(genmat[1+0:(N-k),1+N-k-0:(N-k)]) = diag(genmat[1+0:(N-k),1+N-k-0:(N-k)]) + genmat[1+j, 1+N-k] * dbinom(0:(N-k), N-k, 1-(1-prob)^j)
  }
  for (j in 1:N) {
    genmat[1,1] = genmat[1,1] + genmat[1+j, 1] # Taking care f zero susceptibles left
  }
  rev(genmat[1,])
}

final.size = function(N, R) sum(0:N * final.size.dist(c(0, 1, rep(0, N-1)), R))

mean.generation.dist = function(init.state, R) {
  # init.state is vector (length population size + 1); (k+1)-th element = probability of k initial cases
  # R is basic reproduction (expected number of secondaries with one index and rest susceptible)
  N = length(init.state) - 1
  prob = R/(N-1)
  genmat = diag(init.state * 0)
  k.mat = genmat # mean number of generations needed to get to state k.mat[i,j]
  diag(genmat[1+0:N, 1+N-0:N]) = init.state
  for (k in 1:(N-1)) for (j in 1:k) {
    diag(genmat[1+0:(N-k),1+N-k-0:(N-k)]) = diag(genmat[1+0:(N-k),1+N-k-0:(N-k)]) + genmat[1+j, 1+N-k] * dbinom(0:(N-k), N-k, 1-(1-prob)^j)
  }
  for (j in 1:N) {
    genmat[1,1] = genmat[1,1] + genmat[1+j, 1] 
  }
  for (k in 1:N) for (i in 0:(N-k)) {
    k.mat[1+i, 1+N-k-i] = weighted.mean(1+k.mat[1+1:k, 1+N-k],  genmat[1+1:k, 1+N-k] * dbinom(i, N-k, 1-(1-prob)^(1:k)))
  }
  rev(k.mat[1,])
}

mean.generation = function(N, R) sum(mean.generation.dist(c(0, 1, rep(0, N-1)), R) * final.size.dist(c(0, 1, rep(0, N-1)), R))

# Maximum likelihood estimate of R from a single outbreak
estimate.R = function(N, n, max.R=20) { # outbreak of size n in a facility of size N
  optimize(function(rr) final.size.dist(c(0, 1, rep(0, N-1)), rr)[n+1], c(0, max.R), maximum=TRUE)$maximum
}


# typical size: 50 cases, 215 people
# unmitigated R0: 1.368889

mitigated_size_rf2 = function(population, unmit_R, test_p=7, res_d, symp_HCP_excl=TRUE, isol_eff=1,r_sens=1) {
  # unmit_size = unmitigated size; only mitigation is exclusion of symptomatics, if symp_HCP_excl=TRUE (default)
  # test_p = test periodicity (days); default weekly (7 days)
  # res_d = test result reporting lag
  # isol_eff = effectiveness of isolation once detected (default 1 i.e., 100%)
  # unmit_gen = number of generations of unmitigated (except exclusion of symptomatics if symp_HCP_excl=TRUE) transmission
  symp_trans = isol_eff * transmission_symp(Inf)
  asymp_trans = isol_eff * mean_transmission_reduction2(test_p, res_d, Inf,r_sens)
  r_mult = 1/(1 - symp_HCP_excl * transmission_symp(Inf)) # multiplier accounting for exclusion of symptomatics 
  est_mit_R0 = unmit_R * r_mult * (1 - symp_trans - asymp_trans)
  # obtain distribution after 1 distribution of unmitigated:
  dist_init=dbinom(0:(population-1),population-1,unmit_R/(population-2))
  # obtain size at this point
  extra_size_dist=final.size.dist(dist_init,est_mit_R0)
  final_size=sum(0:(population-1)*extra_size_dist)+1
  final_size
}

mitigated_length_rf2 = function(population, unmit_R, test_p=7, res_d, symp_HCP_excl=TRUE, isol_eff=1,r_sens=1) {
  # unmit_size = unmitigated size; only mitigation is exclusion of symptomatics, if symp_HCP_excl=TRUE (default)
  # test_p = test periodicity (days); default weekly (7 days)
  # res_d = test result reporting lag
  # isol_eff = effectiveness of isolation once detected (default 1 i.e., 100%)
  # unmit_gen = number of generations of unmitigated (except exclusion of symptomatics if symp_HCP_excl=TRUE) transmission
  symp_trans = isol_eff * transmission_symp(Inf)
  asymp_trans = isol_eff * mean_transmission_reduction2(test_p, res_d, Inf,r_sens)
  r_mult = 1/(1 - symp_HCP_excl * transmission_symp(Inf)) # multiplier accounting for exclusion of symptomatics 
  est_mit_R0 = unmit_R * r_mult * (1 - symp_trans - asymp_trans)
  # obtain distribution after 1 distribution of unmitigated:
  dist_init=dbinom(0:(population-1),population-1,unmit_R/(population-2))
  # obtain size at this point
  extra_size_dist=final.size.dist(dist_init,est_mit_R0)
  extra_length_dist=mean.generation.dist(dist_init,est_mit_R0)
  final_length=sum(extra_size_dist*extra_length_dist)+1
  final_length
}

# percent reduction in transmission from non-outbreak testing
non_outbreak2=function(p,t,t_i,r_sens=1) {
  mean_transmission_reduction2(p,t,t_i,r_sens)/(1-transmission_symp(Inf))
}

#set up values needed for table
sens<-c(1,1,1,1,0.85,0.85,0.5,0.5,0.5)
dels<-c(2,1,2,1,0,0,0,0,0)
pers<-c(7,7,3,3,7,3,7,3,1)
t_c<-rep(Inf,9)
true9<-rep(TRUE,9)
ipc1<-rep(1,9)
ipc2<-rep(0.9,9)
base_pop<-rep(215,9)
base_size<-rep(50,9)
base_r<-rep(1.368889,9)

#code to generate dataframe w results that can be exported to CSV
results_table=function()
{
  # 100% IPC effectiveness
  outbreak_size<-mapply(mitigated_size_rf2,base_pop,base_r,pers,dels,true9,ipc1,sens)
  outbreak_len<-mapply(mitigated_length_rf2,base_pop,base_r,pers,dels,true9,ipc1,sens)
  nonos<-mapply(non_outbreak2,pers,dels,t_c,sens)
  outbreak_tests=base_pop*7/pers*(outbreak_len*6.5/7+2)
  outbreak_effect=(base_size-outbreak_size)/base_size
  outbreak_efficiency=outbreak_tests/(base_size-outbreak_size)
  nonoutbreak_efficiency=129*7/pers/(base_size*nonos*0.1)
  combined_effect=1-(1-nonos)*(1-outbreak_effect)
  combined_efficiency=(129*7/pers+0.1*(1-nonos)*outbreak_tests)/(0.1*base_size*combined_effect)
  
  # 90% IPC
  outbreak_size2<-mapply(mitigated_size_rf2,base_pop,base_r,pers,dels,true9,ipc2,sens)
  outbreak_len2<-mapply(mitigated_length_rf2,base_pop,base_r,pers,dels,true9,ipc2,sens)
  nonos2<-mapply(non_outbreak2,pers,dels,t_c,sens)
  outbreak_tests2=base_pop*7/pers*(outbreak_len2*6.5/7+2)
  outbreak_effect2=(base_size-outbreak_size2)/base_size
  outbreak_efficiency2=outbreak_tests2/(base_size-outbreak_size2)
  nonoutbreak_efficiency2=129*7/pers/(base_size*nonos2*0.1)
  combined_effect2=1-(1-nonos2)*(1-outbreak_effect2)
  combined_efficiency2=(129*7/pers+0.1*(1-nonos2)*outbreak_tests2)/(0.1*base_size*combined_effect2)
  
  #make table
  sensitivities<-c("Baseline","Baseline","Baseline","Baseline","85%","85%","50%","50%","50%")
  turnaround_time<-c("48 hours","24 hours","48 hours","24 hours","POC","POC","POC","POC","POC")
  testing_frequency<-c("Weekly","Weekly","Every 3 days","Every 3 days","Weekly","Every 3 days","Weekly","Every 3 days","Daily")
  outbreak_pct_prevented=outbreak_effect
  outbreak_tests_per_case_prevented=outbreak_efficiency
  combined_pct_prevented=combined_effect
  combined_tests_per_case_prevented=combined_efficiency
  worse_IPC_pct_prevented=combined_effect2
  worse_IPC_tests_per_case_prevented=combined_efficiency2
  
  table_out<-data.frame(sensitivities,turnaround_time,
                        testing_frequency,outbreak_pct_prevented,
                        outbreak_tests_per_case_prevented,
                        combined_pct_prevented,
                        combined_tests_per_case_prevented,
                        worse_IPC_pct_prevented,
                        worse_IPC_tests_per_case_prevented)
  table_out
  
}



#new He results
profile = 2
inf_profile = d_inf[[profile]]
inf_max = inf_max_vector[profile]
inf_to_onset = inf_to_onset_vector[profile]
he_table<-results_table()
write.csv(he_table,"he_table2.csv",row.names=FALSE)


#"default" results, i.e., Johansson et al
profile = 1
inf_profile = d_inf[[profile]]
inf_max = inf_max_vector[profile]
inf_to_onset = inf_to_onset_vector[profile]
default_table<-results_table()
write.csv(default_table,"default_table2.csv",row.names=FALSE)


#Goyal results
profile = 3
inf_profile = d_inf[[profile]]
inf_max = inf_max_vector[profile]
inf_to_onset = inf_to_onset_vector[profile]
goyal_table<-results_table()
write.csv(goyal_table,"goyal_table2.csv",row.names=FALSE)


#Clifford results
profile = 4
inf_profile = d_inf[[profile]]
inf_max = inf_max_vector[profile]
inf_to_onset = inf_to_onset_vector[profile]
clifford_table<-results_table()
write.csv(clifford_table,"clifford_table2.csv",row.names=FALSE)



#plot infectiousness profiles using symptom onset as reference
#first translate x-axis to align with symptom onset date
default_g=function(x) { d_inf[[1]](x+inf_to_onset_vector[1])}
he_new_g=function(x) { d_inf[[2]](x+inf_to_onset_vector[2])}
goyal_g=function(x) { d_inf[[3]](x+inf_to_onset_vector[3])}
clifford_g=function(x) { d_inf[[4]](x+inf_to_onset_vector[4])}

#plot
tttt<-ggplot(data.frame(x=c(-10,20)),aes(x=x))+stat_function(fun=he_new_g,aes(colour="He"))+
    stat_function(fun=default_g,aes(colour="Johansson"))  +
    stat_function(fun=goyal_g,aes(colour="Goyal"))  +
    stat_function(fun=clifford_g,aes(colour="Clifford"))  +
    scale_x_continuous(name="Days from symptom onset",breaks=seq(-10,20,5),limits=c(-10,20))+scale_y_continuous(name="Infectiousness")+
    scale_colour_manual("Profiles",values=c("blue","red","purple","green")) +
    theme_bw()


#display                                                                                                                                                                                                                        
tttt