###############
## This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
## <http://creativecommons.org/licenses/by-nc/4.0/> by James Trauer, Arathi Arakala
## This work is supported by -----

## Per the terms of this license, if you are making derivative use of this work, you must identify that 
## your work is a derivative work, give credit to the original work, provide a link to the license, 
## and indicate changes that were made.
###############
## This code will introduce variable CDR over 20 years, specifically based on what was seen in India ##
######################


require(deSolve)
require(graphics)
require(ggplot2)
require(xlsx)
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
  #delta =0,
  delta=0*73/(0.41*365*70*3),
  ### disease based parameters ####################
  #################################################
  epsilon = 0.0011, #rate of progression to active disease from early latency per day, 
  kappa = 0.01, # rate of progression from early to late latency
  nu = 5.5e-6, #rate of progression to active disease from late latency per day.
  gamma = 0.5/(3*365), # rate of spontaneous cure, Dowdy et al; 
  o=0, # fraction of treated individuals contributing to the transmission rate
  phi= (0.74*2)/365, # rate of recovery, per year. is 74% from golabl TB report for India
  omega = (0.11*2)/365, # probability of defaulting treatment
  beta=24/365, #contact rate of infection
  chi = 0.21, # fractional reduction in the force of infection corresponding to return from late to early latency
  alpha = 0.5, # fractional reduction in force of infection due to vaccination
  rho = 0.3 # proportion of infected people who are infectious. avg over years, see Notifications table , WHO TB report.
)
###########



##########
## CDR taken from global TB report and can vary. So create a function to alter that parameter
#####
changeCDR_parameter<-function(parameters, year, cdr_data_India){
  rowindex<-which(cdr_data_India[,1]==year)
  cdr_new<-cdr_data_India[rowindex,2]/100
  parameters_new<-parameters
  delta_new<-(cdr_new/(1-cdr_new))*(73/(365*70*3))
  parameters_new["delta"]<-delta_new
  parameters_new
}

## Model structure and Equations, CDR varying with time
TBmodel_varCDR<-function(t,n,parameters){
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
    
    year<-trunc(t) # change delta based on CDR
    params<-changeCDR_parameter(parameters, year, cdr_data_India)
    delta<-params["delta"]
    
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



## Read data in from Global TB report file on India

#inputfile<-"/Users/Arathi/Documents/2017/Monash/TB/WHO Global TB report 2016 data/TB_burden_countries_2017-10-10.csv"
inputfile<-"/Users/Arathi/Documents/2017/Monash/TB/WHO Global TB report 2017/TB_burden_countries_2017-11-09.csv"
MyData<- read.csv(file=inputfile, header=TRUE, sep=",")
MyData_India<-MyData[which(MyData$country=="India"), ]
cdr_India<-MyData_India$c_cdr
year_India<-MyData_India$year
cdr_data_India<-cbind(year_India, cdr_India)
#pre GTB data
years<-1950:1999
cdr<-c(rep(0, times=10), rep(10, times=33), rep(20, times=7))
preGTR_India<-cbind(years,cdr)
cdr_data_India<-rbind(preGTR_India, cdr_data_India)
colnames(cdr_data_India)<-c("year_India", "cdr_India")

# quartz()
# par(mfrow=c(1,1), oma=c(0,0,2,0))
# plot(year_India, cdr_India, type="b")
# year<-2015
# params<-changeCDR_parameter(parameters, year, cdr_data_India)
# params["delta"]


#Burn in time,> 500 years
time<-seq(from=1, to=365*200, by=1 )
nstart=c(S_v=1000000,S_u=1000000, L_a=0, L_b=0, I=3, T_r=0, S_r=0)
output<-as.data.frame(ode(nstart,time,TBmodel_basic,parameters))
eqbm<-output[dim(output)[1], ]

#### Calibrate the model beta so that model output in 2000 is approx. that in year 2000 from global TB report
N=output$S_v+output$S_u+output$L_a+output$L_b+output$I+output$T_r+output$S_r #total population

newI<-((output$L_b)*parameters["nu"])+(output$L_a*parameters["epsilon"]) #new I at every timestep
# timestep = 1 day = 1/365 years.
Incidence<- (365)*(newI/N)*100000

quartz()
par(mfrow=c(1,1), oma=c(0,0,2,0))
plot(output$time, Incidence, type='l', lwd=2, col=1, pch=1, ylim=c(0, max(MyData_India$e_inc_100k)+200), xlab="years", ylab="Incidence per 100k", cex.axis=1.5, cex.lab=1.5)
abline(h=MyData_India$e_inc_100k[1]+75, col=2, lwd=2)
title(paste("beta =", parameters["beta"], sep=""), outer=TRUE)

###########################################################################

N=output$S_v+output$S_u+output$L_a+output$L_b+output$I+output$T_r+output$S_r #total population
print(c(eqbm$S_v/N[dim(output)[1]], eqbm$S_u/N[dim(output)[1]], eqbm$L_a/N[dim(output)[1]], eqbm$L_b/N[dim(output)[1]], eqbm$I/N[dim(output)[1]], eqbm$T_r/N[dim(output)[1]], eqbm$S_r/N[dim(output)[1]] ))
print(eqbm$I/N[dim(output)[1]]+eqbm$T_r/N[dim(output)[1]])
quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(output$time/365, output$S_v/N, type='l', lwd=2, col=1, ylim=c(0, max(c(output$S_v/N,output$S_u/N, output$S_r/N)) ), lty=1, xlab="years", ylab="Susceptibles", cex.axis=1.5, cex.lab=1.5)
lines(output$time/365, output$S_u/N, lwd=2, col=2)
lines(output$time/365, output$S_r/N, lwd=2, col=3)
legend(x=head(output$time/365, n=1), y=max(c(output$S_v/N,output$S_u/N, output$S_r/N)), legend=c("S_v", "S_u", "S_r"), lwd=2, lty=1, col=c(1,2,3), cex=1.5)

plot(output$time/365, output$L_a/N, type='l', lwd=2, col=4, lty=2, ylim=c(0, max( c(output$L_a/N,output$L_b/N) )), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
lines(output$time/365, output$L_b/N, lwd=2, col=5, lty=2)
legend(x=head(output$time/365, n=1), y=max( c(output$L_a/N,output$L_b/N) ), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,5), cex=1.5)

plot(output$time/365, output$I/N, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Infectious", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
legend(x=head(output$time/365, n=1), y=max(output$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)

plot(output$time/365, output$T_r/N, type='l', lwd=2, col=7, lty=2, xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 ) #, ylim=c(0, max(output_all$T_r/N))
legend(x=head(output$time/365, n=1), y=max(output$T_r/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)
title(main="Compartmental TB model run-in output", outer=TRUE)


output_all<-numeric()
delta_all<-numeric()

output_all<-rbind(output_all, output[(dim(output)[1]-365):dim(output)[1], ])
output_all$time<-seq(from=1949, to=1950, by=1/365)

for(year in 1950:2016){ #2016
  nstart=c(S_v=eqbm$S_v, S_u=eqbm$S_u, L_a=eqbm$L_a, L_b=eqbm$L_b, I=eqbm$I, T_r=eqbm$T_r, S_r=eqbm$S_r)
  time<-seq(from=1,to=(366), by=1)
  year_seq<-seq(from=year, to=(year+1), by=1/365)
  params<-changeCDR_parameter(parameters, year, cdr_data_India)
  #params["delta"]<-0
  delta_all<-c(delta_all, params["delta"])
  output<-as.data.frame(ode(nstart,time,TBmodel_basic,params))
  eqbm<-output[dim(output)[1], ]

  output<-output[-dim(output)[1],] # remove last as it will be starting values for next year
  output$time<-year_seq[1:365]
  output_all<-rbind(output_all, output)
  
  eqbm
  
  if(year==2016){
    quartz()
    N=output_all$S_v+output_all$S_u+output_all$L_a+output_all$L_b+output_all$I+output_all$T_r+output_all$S_r #total population
    
    plot(output_all$time, ((output_all$I+output_all$T_r)/N)*100000, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Prev", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
    legend(x=2000, y=max(output_all$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)
    
  }
  
  
    
}
N=output_all$S_v+output_all$S_u+output_all$L_a+output_all$L_b+output_all$I+output_all$T_r+output_all$S_r #total population

quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(1950:2016, delta_all, type='b', ylab="delta")

plot(cdr_data_India, type='b', ylab="cdr")


plot(output_all$time, output_all$I/N, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Infected", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
legend(x=2000, y=max(output_all$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)

plot(output_all$time, output_all$T_r/N, type='l', lwd=2, col=7, lty=2, xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 ) #, ylim=c(0, max(output_all$T_r/N))
legend(x=2000, y=max(output_all$T_r/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)



quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))

plot(output_all$time, output_all$S_v/N, type='l', lwd=2, col=1, ylim=c(0, max(c(output_all$S_v/N,output_all$S_u/N, output_all$S_r/N)) ), lty=1, xlab="years", ylab="Susceptibles", cex.axis=1.5, cex.lab=1.5)
lines(output_all$time, output_all$S_u/N, lwd=2, col=2)
lines(output_all$time, output_all$S_r/N, lwd=2, col=3)
legend(x=1970, y=max(c(output_all$S_v/N,output_all$S_u/N, output_all$S_r/N)), legend=c("S_v", "S_u", "S_r"), lwd=2, lty=1, col=c(1,2,3), cex=1.5)


plot(output_all$time, output_all$L_a/N, type='l', lwd=2, col=4, lty=2, ylim=c(0, max( c(output_all$L_a/N,output_all$L_b/N) )), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
lines(output_all$time, output_all$L_b/N, lwd=2, col=5, lty=2)
legend(x=1970, y=max( c(output_all$L_a/N,output_all$L_b/N) ), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,5), cex=1.5)


plot(output_all$time, output_all$I/N, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Infected", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
legend(x=1970, y=max(output_all$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)

plot(output_all$time, output_all$T_r/N, type='l', lwd=2, col=7, lty=2, xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 ) #, ylim=c(0, max(output_all$T_r/N))
legend(x=1970, y=max(output_all$T_r/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)
title(main="Compartmental TB model output", outer=TRUE)


###############################################
#### Plot incidence and prevelance over time
#### compare with global tb report

newI<-((output_all$L_b)*params["nu"])+(output_all$L_a*params["epsilon"]) #new I at every timestep
# timestep = 1 day = 1/365 years.
Incidence<- (365)*(newI/N)*100000

quartz()
plot(output_all$time, Incidence, type='l', lwd=2, col=1, pch=1, ylim=c(0, max(MyData_India$e_inc_100k)), xlab="years", ylab="Incidence per 100k")
lines(MyData_India$year, MyData_India$e_inc_100k, type='l', lty=2, col=2, pch=2, lwd=2)
legend(1960,max(MyData_India$e_inc_100k), legend = c("model output","global TB report"), lty=c(1,2), col=c(1,2) )
abline(v=2000, lty=3, col=3, lwd=2)

prev<-(output_all$I/N+output_all$T_r/N)*100000
quartz()
plot(output_all$time, prev, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Prevalence", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 


quartz()
par(mfrow=c(2,2), oma=c(0,0,2,0))

plot(output_all$time, Incidence, type='l', lwd=2, col=1, pch=1, ylim=c(0, max(MyData_India$e_inc_100k)), xlab="years", ylab="Incidence per 100k", cex.axis=1.5, cex.lab=1.5)
lines(MyData_India$year, MyData_India$e_inc_100k, type='l', lty=2, col=2, pch=2, lwd=2)
legend(1960,max(MyData_India$e_inc_100k), legend = c("model output","global TB report"), lty=c(1,2), col=c(1,2) )
abline(v=2000, lty=3, col=3, lwd=2)

plot(output_all$time, prev, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Prevalence", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
plot(cdr_data_India[,1], cdr_data_India[,2], type='l', ylab="cdr", xlab="years", cex.axis=1.5, cex.lab=1.5)



#Inc<-(1/timestep)*(I_new/N)*100000
#Prev<-(1/timestep)*(I+T)/N *100000
