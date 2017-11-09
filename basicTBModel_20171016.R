###############
## This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
## <http://creativecommons.org/licenses/by-nc/4.0/> by James Trauer, Arathi Arakala
## This work is supported by -----

## Per the terms of this license, if you are making derivative use of this work, you must identify that 
## your work is a derivative work, give credit to the original work, provide a link to the license, 
## and indicate changes that were made.
###############



require(deSolve)
require(graphics)
require(ggplot2)

## Model structure and Equations
TBmodel_basic<-function(t,n,parameters){
  with(as.list(parameters),{
    
    S_v=n[1] # number of people vaccinated
    S_u=n[2] #number of people unvaccinated
    L_a=n[3] #number in early latency
    L_b=n[4] #number in late latency
    I=n[5] # number infected
    T_r=n[6] #number of people given treatment
    S_r=n[7] # number of people recovered after treatment
    
    N=S_v+S_u+L_a+L_b+I+T_r+S_r #total population
    lamda=rho*(beta* (I+(o*T_r) ) )/N # force of infection
    
    dS_vdt=(i*pi*N) - (mu*S_v) - (alpha*lamda*S_v) 
    dS_udt=((1-i)*pi*N) - (mu*S_u) - lamda*S_u
    dL_adt=(alpha*lamda*S_v) + lamda*S_u + (alpha*lamda*S_r) + (chi*lamda*L_b) - (mu+epsilon+kappa)*L_a
    dL_bdt= kappa*L_a + gamma*I -(mu+nu+(chi*lamda) )*L_b
    dIdt= nu*L_b + epsilon*L_a + omega*T_r - (mu_I+gamma+delta) * I
    dT_rdt= delta*I - (omega+phi+mu_T) * T_r
    dS_rdt= phi*T_r - ( (alpha*lamda) + mu) * S_r
    
    return(list(c(dS_vdt, dS_udt, dL_adt,dL_bdt, dIdt, dT_rdt,dS_rdt )))
  })
  
}



## Parameters of the model, all taken from Trauer et al 2014, table 1. All per year rates
parameters<-c(
  ### demography parameters based on geography ####
  #################################################
  pi=0.0193/365, #birth rate into the population, per year, for India in 2014
  i=0.97, #coverage of vaccination, BCG rate, for India average
  mu = (1/70)/365, # standard death rate, in India in 2014, if average life expectancy is 70
  
  ### disease parameters that vary with geography ###
  ####################################################
  mu_I = 0.5/(3*365), # daily death rate due to disease in India, WHO world TB report
  mu_T = (0.15*2)/365, # death rate during treatment, half of mu_I
  #delta = 0.59*73/(0.41*365*70*3), # based on CDR 2015, glbal TB report for India, probability of treatment=0.59 , will be typically time dependent
  delta =0,
  ### disease based parameters ####################
  #################################################
  epsilon = 0.0011, #rate of progression to active disease from early latency per day, 
  kappa = 0.01, # rate of progression from early to late latency
  nu = 5.5e-6, #rate of progression to active disease from late latency per day.
  gamma = 0.5/(3*365), # rate of spontaneous cure, Dowdy et al; 
  o=0.21, # fraction of treated individuals contributing to the transmission rate
  phi= (0.74*2)/365, # rate of recovery, per year. is 74% from golabl TB report for India
  omega = (0.11*2)/365, # probability of defaulting treatment
  beta=36.5/365, #contact rate of infection
  chi = 0.49, # fractional reduction in the force of infection corresponding to return from late to early latency
  alpha = 0.5, # fractional reduction in force of infection due to vaccination
  rho = 0.3 # proportion of infected people who are infectious. avg over years, see Notifications table , WHO TB report.
)

# #Burn in time, 100 years
# time<-seq(from=0, to=365*100, by=1)
# nstart=c(S_v=1000000,S_u=1000000, L_a=0, L_b=0, I=1, T_r=0, S_r=0)
# output<-as.data.frame(ode(nstart,time,TBmodel_basic,parameters))
# eqbm<-output[dim(output)[1], ]
# eqbm


## Test the model output
time<-seq(from=0, to=365*500, by=1)
nstart=c(S_v=1000000,S_u=1000000, L_a=0, L_b=0, I=100, T_r=0, S_r=0)
output<-as.data.frame(ode(nstart,time,TBmodel_basic,parameters))
eqbm<-output[dim(output)[1], ]
eqbm
# # plot of cumulative incidence
# I_cum<-cumsum(output$I)
# quartz()
# plot(round(output$time/365), I_cum, type="l", lwd=2, col=1, ylab="Cumulative Incidence", xlab="Time from infection (years)")

N=output$S_v+output$S_u+output$L_a+output$L_b+output$I+output$T_r+output$S_r #total population
lamda=as.numeric(parameters["rho"])*( as.numeric(parameters["beta"]) * (output$I+(as.numeric(parameters["o"])*output$T_r) ) )/N # force of infection

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))

plot(output$time/365, output$S_v/N, type='l', lwd=2, col=1, ylim=c(0, max(c(output$S_v/N,output$S_u/N, output$S_r/N)) ), lty=1, xlab="years", ylab="Susceptibles", cex.axis=1.5, cex.lab=1.5)
lines(output$time/365, output$S_u/N, lwd=2, col=2)
lines(output$time/365, output$S_r/N, lwd=2, col=3)
legend(x=40, y=max(c(output$S_v/N,output$S_u/N, output$S_r/N)), legend=c("S_v", "S_u", "S_r"), lwd=2, lty=1, col=c(1,2,3), cex=1.5)


plot(output$time/365, output$L_a/N, type='l', lwd=2, col=4, lty=2, ylim=c(0, max( c(output$L_a/N,output$L_b/N) )), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
lines(output$time/365, output$L_b/N, lwd=2, col=5, lty=2)
legend(x=3, y=max( c(output$L_a/N,output$L_b/N) ), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,5), cex=1.5)


plot(output$time/365, output$I/N, type='l', lwd=2, col=6, lty=2, ylim=c(0, max(output$I/N)) , xlab="years", ylab="Infected", cex.axis=1.5, cex.lab=1.5)
legend(x=3, y=max(output$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)

plot(output$time/365, output$T_r/N, type='l', lwd=2, col=7, lty=2, ylim=c(0, max(output$T_r/N)), xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 )
legend(x=3, y=max(output$T_r/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)
title(main="Compartmental TB model output", outer=TRUE)

Prev<-((output$I+output$T_r)/N) * 100000
Inc<-( ((as.numeric(parameters["nu"])*output$L_b) + ( parameters["epsilon"]*output$L_a) )/N )*365 * 100000
quartz()
par(mfrow=c(1,2), oma=c(0,0,2,0))
plot(output$time/365, Prev, type="l", lwd=2, col=8, xlab="years", ylab="Prevelance", cex.axis=1.5, cex.lab=1.5  )
plot(output$time/365, Inc, type="l", lwd=2, col=9, xlab="years", ylab="Incidence", cex.axis=1.5, cex.lab=1.5  )

# quartz()
# par(mfrow=c(2,2), oma=c(0,0,2,0))
# plot(output$time/365, lamda, type="l", lwd=2, main="lamda")
# plot(output$time/365, N, type="l", lwd=2, main="total population")
# plot(output$time/365, output$L_a, type='l', lwd=2, col=4, lty=2, ylim=c(0, max(output$L_a+output$L_b)), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
# lines(output$time/365, output$L_b, lwd=2, col=2, lty=2)
# plot(output$time/365, output$I, type='l', lwd=2, col=5, lty=2, ylim=c(0, max(output$I)) , xlab="years", ylab="Infecteds", cex.axis=1.5, cex.lab=1.5)


