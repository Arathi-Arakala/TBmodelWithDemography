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
    lamda=rho*(beta*I)/N # force of infection
    
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



## Parameters of the model, all taken from Trauer et al 2014, table 1.
parameters<-c(
  pi=0.025, #birth rate into the population, per year
  i=0.65, #rate of vaccination, BCG rate
  beta=24, #transmission rate of infection
  o=0.21, # fraction of treated individuals contributing to the transmission rate
  epsilon = 0.129, #rate of progression to active disease from early latency
  nu = 0.075, #rate of progression to active disease from late latency
  kappa = 0.821, # rate of progression from early to late latency
  gamma = 0.63, # rate of spontaneous cure
  mu = 0.016, # standard death rate
  mu_I = 0.37, # death rate due to disease
  mu_T = 0.5*0.37, # death rate during treatment, half of mu_I
  delta = 0.72, # probability of treatment , will be typically time dependent
  phi= 2, # rate of recovery, per year
  omega = 0.25, # probability of defaulting treatment
  chi = 0.49, # fractional reduction in the force of infection corresponding to return from late to early latency
  alpha = 10, # fractional reduction in force of infection due to vaccination
  rho = 0.35 #infectious proportion
)

## Test the model output
time<-seq(from=0, to=2*365, by=1)
nstart=c(S_v=1000000,S_u=1000000, L_a=0, L_b=0, I=1, T_r=0, S_r=0)
output<-as.data.frame(ode(nstart,time,TBmodel_basic,parameters))
eqbm<-output[dim(output)[1], ]


quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(output$time, output$S_v, type='l', lwd=2, col=1, ylim=c(0, max(output$S_v+output$S_u)), lty=1)
lines(output$time, output$S_u, lwd=2, col=3, lty=2)
legend(x=1000, y=max(output$S_v+output$S_u), legend=c("S_v", "S_u"), lwd=2, lty=c(1,1), col=c(1,3))

plot(output$time, output$L_a, type='l', lwd=2, col=4, lty=2, ylim=c(0, max(output$L_a+output$L_b)) )
lines(output$time, output$L_b, lwd=2, col=2, lty=2)
legend(x=1000, y=max(output$L_a+output$L_b), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,2))

plot(output$time, output$I, type='l', lwd=2, col=5, lty=2, ylim=c(0, max(output$I)) )
legend(x=1000, y=max(output$I), legend=c("I"), lwd=2, lty=1, col=c(5))

plot(output$time, output$T_r, type='l', lwd=2, col=6, lty=2, ylim=c(0, max(output$T_r+output$S_r)) )
lines(output$time, output$S_r, lwd=2, col=7, lty=2)
legend(x=1000, y=max(output$T_r+output$S_r), legend=c("T_r", "S_r"), lwd=2, lty=c(1,1), col=c(6,7))



