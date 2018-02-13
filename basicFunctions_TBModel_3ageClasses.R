## Parameters of the model, all taken from Trauer et al 2014, table 1. All per year rates
parameters<-list(
  ### demography parameters based on geography ####
  #################################################
  pi=0.0193/365, #birth rate into the population, per year, for India in 2014
  i=0.97, #coverage of vaccination, BCG rate, for India average
  #mu = c( 1000000 , 100000,100000 ),
  #mu = c( (1/70)/365 , (1/70)/365,(1/70)/365 ), # standard death rate, in India in 2014, if average life expectancy is 70
  mu = c((1/68)/365, (1/63)/365, (1/53)/365 ) ,# age structured death rate, exp lifespan model, if life expectancy is 68
  tau_in = c(0, 1/(5*365), 1/(10*365)), # ageing parameters
  tau_out = c(1/(5*365), 1/(10*365), 0), # ageing parameters
  
  ### disease parameters that vary with geography ###
  ####################################################
  mu_I = 0.5/(3*365), # daily death rate due to disease in India, WHO world TB report
  mu_T = (0.15*2)/365, # death rate during treatment, half of mu_I
  #delta = 0.59*73/(0.41*365*70*3), # based on CDR 2015, glbal TB report for India, probability of treatment=0.59 , will be typically time dependent
  #delta =0,
  delta=0*73/(0.41*365*70*3),# prob of treatment=0, no treatment case
  ### disease based parameters ####################
  #################################################
  #epsilon = rep(0.0011, times=3), #rate of progression to active" disease from early latency per day, 
  epsilon = c(6.6E-03, 2.7E-03, 2.7E-04), #from Romain's paper "optimally......."
  #kappa = rep(0.01, times=3), # rate of progression from early to late latency
  kappa = c(1.2E-02 ,1.2E-02 , 5.4E-03), # rate of progression from early to late latency
  #nu = rep(5.5e-6, times=3) , #rate of progression to active disease from late latency per day.
  nu = c( 1.9E-11,6.4E-06,3.3E-06 ), #rate of progression to active disease from late latency per day.
  #iota = c(1,1,1) , # same for all age classes
  iota = c(0.1, 0.5, 1), # change in infectiousness compared to an adult for the 3 age classes
  
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
  parameters_new$delta<-delta_new
  parameters_new
}

TBmodel_3AgeClasses<-function(t,n,parameters){
  with(as.list(parameters),{
    
    #read in parameters
    pi<-parameters$pi
    i<-parameters$i
    mu<-parameters$mu
    tau_in<-parameters$tau_in
    tau_out<-parameters$tau_out
    mu_I<-parameters$mu_I
    mu_T<-parameters$mu_T
    delta<-parameters$delta
    epsilon<-parameters$epsilon
    kappa<-parameters$kappa
    nu<-parameters$nu
    gamma<-parameters$gamma
    o<-parameters$o
    phi<-parameters$phi
    omega<-parameters$omega
    beta<-parameters$beta
    chi<-parameters$chi
    alpha<-parameters$alpha
    iota<-parameters$iota
    rho<-parameters$rho

    
    # read in intial values of state variables
    S_v<-n[1:3] # number of people vaccinated
    S_u<-n[4:6] #number of people unvaccinated
    L_a<-n[7:9] #number in early latency
    L_b<-n[10:12] #number in late latency
    I<-n[13:15] # number infected
    T_r<-n[16:18] #number of people given treatment
    S_r<-n[19:21] # number of people recovered after treatment
    
    N<-rep(0, times=3) #total population in each age class
    N_tot<-sum(S_v)+sum(S_u)+sum(L_a)+sum(L_b)+sum(I)+sum(T_r)+sum(S_r)
    dS_vdt<-rep(0, times=3)
    dS_udt<-rep(0, times=3)
    dL_adt<-rep(0, times=3)
    dL_bdt<-rep(0, times=3)
    dIdt<-rep(0, times=3)
    dT_rdt<-rep(0, times=3)
    dS_rdt<-rep(0, times=3)
    for(a in 1:3){
      N[a]=S_v[a]+S_u[a]+L_a[a]+L_b[a]+I[a]+T_r[a]+S_r[a] #total population in age class a
      lamda=iota[a]*rho*(beta* (sum(I)+(o*sum(T_r)))/ N_tot ) # force of infection, homogenous mixing
      #lamda=rho*(beta* (sum(I)+(o*sum(T_r)))/ N_tot ) # force of infection, homogenous mixing
      
      # 
      # dS_vdt[a]=(i*pi*N_tot) - (mu*S_v[a]) - (alpha*lamda*S_v[a]) - (tau_out[a]*S_v[a])+(tau_in[a]*S_v[a-1])
      # dS_udt[a]=((1-i)*pi*N_tot) - (mu*S_u[a]) - lamda*S_u[a] - (tau_out[a]*S_u[a])+(tau_in[a]*S_u[a-1])
      # 
      # dL_adt[a]=(alpha*lamda*S_v[a]) + lamda*S_u[a] + (alpha*lamda*S_r[a]) + (chi*lamda*L_b[a]) - (mu+epsilon+kappa)*L_a[a] - (tau_out[a]*L_a[a])
      # dL_bdt[a]= kappa*L_a[a] + gamma*I[a] -(mu+nu+(chi*lamda) )*L_b[a] - (tau_out[a]*L_b[a])
      # dIdt[a]= nu*L_b[a] + epsilon*L_a[a] + omega*T_r[a] - (mu_I+gamma+delta) * I[a] - (tau_out[a]*I[a])
      # dT_rdt[a]= delta*I[a] - (omega+phi+mu_T) * T_r[a] - (tau_out[a]*T_r[a])
      # dS_rdt[a]= phi*T_r[a] - ( (alpha*lamda) + mu) * S_r[a] - (tau_out[a]*S_r[a])
      # 
      
      if(a==1){
        dS_vdt[a]=(i*pi*N_tot) - (mu[a]*S_v[a]) - (alpha*lamda*S_v[a]) - (tau_out[a]*S_v[a])
        dS_udt[a]=((1-i)*pi*N_tot) - (mu[a]*S_u[a]) - lamda*S_u[a] - ((tau_out[a])*S_u[a])
        
        dL_adt[a]=(alpha*lamda*S_v[a]) + lamda*S_u[a] + (alpha*lamda*S_r[a]) + (chi*lamda*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] - ((tau_out[a])*L_a[a])
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda) )*L_b[a] - ((tau_out[a])*L_b[a])
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] + omega*T_r[a] - (mu_I+gamma+delta) * I[a] - ((tau_out[a])*I[a])
        dT_rdt[a]= delta*I[a] - (omega+phi+mu_T) * T_r[a] - ((tau_out[a])*T_r[a])
        dS_rdt[a]= phi*T_r[a] - ( (alpha*lamda) + mu[a]) * S_r[a] - ((tau_out[a])*S_r[a])
      }
      if(a==2){
        dS_vdt[a]= - (mu[a]*S_v[a]) - (alpha*lamda*S_v[a]) +(tau_in[a])*S_v[a-1] - (tau_out[a])*S_v[a]
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda*S_u[a] +(tau_in[a])*S_u[a-1] - (tau_out[a])*S_u[a]
        
        dL_adt[a]=(alpha*lamda*S_v[a]) + lamda*S_u[a] + (alpha*lamda*S_r[a]) + (chi*lamda*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1] - (tau_out[a])*L_a[a]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda) )*L_b[a] + (tau_in[a])*L_b[a-1] - (tau_out[a])*L_b[a]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] + omega*T_r[a] - (mu_I+gamma+delta) * I[a] + (tau_in[a])*I[a-1] - (tau_out[a])*I[a]
        dT_rdt[a]= delta*I[a] - (omega+phi+mu_T) * T_r[a] + (tau_in[a])*T_r[a-1] - (tau_out[a])*T_r[a]
        dS_rdt[a]= phi*T_r[a] - ( (alpha*lamda) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1] - (tau_out[a])*S_r[a]
      }
      
      if(a==3){
        dS_vdt[a]= - (mu[a]*S_v[a]) - (alpha*lamda*S_v[a]) + (tau_in[a])*S_v[a-1]
        dS_udt[a]= - (mu[a]*S_u[a]) - lamda*S_u[a] + (tau_in[a])*S_u[a-1]
        
        dL_adt[a]=(alpha*lamda*S_v[a]) + lamda*S_u[a] + (alpha*lamda*S_r[a]) + (chi*lamda*L_b[a]) - (mu[a]+epsilon[a]+kappa[a])*L_a[a] + (tau_in[a])*L_a[a-1]
        dL_bdt[a]= kappa[a]*L_a[a] + gamma*I[a] -(mu[a]+nu[a]+(chi*lamda) )*L_b[a] + (tau_in[a])*L_b[a-1]
        dIdt[a]= nu[a]*L_b[a] + epsilon[a]*L_a[a] + omega*T_r[a] - (mu_I+gamma+delta) * I[a] + (tau_in[a])*I[a-1]
        dT_rdt[a]= delta*I[a] - (omega+phi+mu_T) * T_r[a] + (tau_in[a])*T_r[a-1]
        dS_rdt[a]= phi*T_r[a] - ( (alpha*lamda) + mu[a]) * S_r[a] + (tau_in[a])*S_r[a-1]
        
      }
      
    }
    
    return(list(c(dS_vdt, dS_udt, dL_adt,dL_bdt, dIdt, dT_rdt,dS_rdt )))
  })
}


plot_output_3ageClasses<-function(output){
  eqbm<-output[dim(output)[1], ]
  for(a in 2:4){
    N<-output[,a]+output[,(a+3)]+output[,(a+6)]+output[,(a+9)]+output[,(a+12)]+output[,(a+15)]+output[,(a+18)] #total population
    print(c(eqbm[,a]/N[dim(output)[1]], eqbm[,(a+3)]/N[dim(output)[1]], eqbm[,(a+6)]/N[dim(output)[1]], eqbm[,(a+9)]/N[dim(output)[1]], eqbm[,(a+12)]/N[dim(output)[1]], eqbm[,(a+15)]/N[dim(output)[1]], eqbm[,(a+18)]/N[dim(output)[1]] ))
    print(eqbm[,(a+12)]/N[dim(output)[1]]+eqbm[,(a+15)]/N[dim(output)[1]])
    quartz()
    par(mfrow=c(2,2), oma=c(0,0,2,0))
    plot(output$time/365, output[,(a)]/N, type='l', lwd=2, col=1, ylim=c(0, max(c(output[,(a)]/N,output[,(a+3)]/N, output[,(a+18)]/N)) ), lty=1, xlab="years", ylab="Susceptibles", cex.axis=1.5, cex.lab=1.5)
    lines(output$time/365, output[,(a+3)]/N, lwd=2, col=2)
    lines(output$time/365, output[,(a+18)]/N, lwd=2, col=3)
    legend(x=head(output$time/365, n=1), y=max(c(output[,(a+3)]/N,output[,(a+6)]/N, output[,(a+18)]/N)), legend=c("S_v", "S_u", "S_r"), lwd=2, lty=1, col=c(1,2,3), cex=1.5)
    
    plot(output$time/365, output[,(a+6)]/N, type='l', lwd=2, col=4, lty=2, ylim=c(0, max( c(output[,(a+6)]/N,output[,(a+9)]/N) )), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
    lines(output$time/365, output[,(a+9)]/N, lwd=2, col=5, lty=2)
    legend(x=head(output$time/365, n=1), y=max( c(output[,(a+6)]/N,output[,(a+9)]/N) ), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,5), cex=1.5)
    
    plot(output$time/365, output[,(a+12)]/N, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Infectious", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
    legend(x=head(output$time/365, n=1), y=max(output[,(a+12)]/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)
    
    plot(output$time/365, output[,(a+15)]/N, type='l', lwd=2, col=7, lty=2, xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 ) #, ylim=c(0, max(output_all$T_r/N))
    legend(x=head(output$time/365, n=1), y=max(output[,(a+15)]/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)
    title(main=paste("Compartmental TB model run-in output for age class ", a-1, sep=" "), outer=TRUE)
    
  }
}


plot_aggregatedoutput_3ageClasses<-function(output){
  S_v_agg<-rowSums(output[,2:4])
  S_u_agg<-rowSums(output[,5:7])
  L_a_agg<-rowSums(output[,8:10])
  L_b_agg<-rowSums(output[,11:13])
  I_agg<-rowSums(output[,14:16])
  T_r_agg<-rowSums(output[,17:19])
  S_r_agg<-rowSums(output[,20:22])
  
  N<-rowSums(output[,2:22])
  quartz()
  par(mfrow=c(2,2), oma=c(0,0,2,0))
  
  plot(output$time/365, S_v_agg/N, type='l', lwd=2, col=1, ylim=c(0, max(c(S_v_agg/N,S_u_agg/N, S_r_agg/N)) ), lty=1, xlab="years", ylab="Susceptibles", cex.axis=1.5, cex.lab=1.5)
  lines(output$time/365, S_u_agg/N, lwd=2, col=2)
  lines(output$time/365, S_r_agg/N, lwd=2, col=3)
  legend(x=1970, y=max(c(S_v_agg/N,S_u_agg/N, S_r_agg/N)), legend=c("S_v", "S_u", "S_r"), lwd=2, lty=1, col=c(1,2,3), cex=1.5)
  
  
  plot(output$time/365, L_a_agg/N, type='l', lwd=2, col=4, lty=2, ylim=c(0, max( c(L_a_agg/N,L_b_agg/N) )), xlab="years", ylab="Latents", cex.axis=1.5, cex.lab=1.5 )
  lines(output$time/365, L_b_agg/N, lwd=2, col=5, lty=2)
  legend(x=1970, y=max( c(L_a_agg/N,L_b_agg/N) ), legend=c("L_a", "L_b"), lwd=2, lty=c(1,1), col=c(4,5), cex=1.5)
  
  
  plot(output$time/365, I_agg/N, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Infected", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
  legend(x=1970, y=max(I_agg/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)
  
  plot(output$time/365, T_r_agg/N, type='l', lwd=2, col=7, lty=2, xlab="years", ylab="Treated ", cex.axis=1.5, cex.lab=1.5 ) #, ylim=c(0, max(output_all$T_r/N))
  legend(x=1970, y=max(T_r_agg/N), legend=c("T_r"), lwd=2, lty=c(1,1), col=c(7), cex=1.5)
  title(main="Compartmental TB model aggregated output from 3 age classes", outer=TRUE)
  
  
}