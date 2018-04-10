###### Create an age stratified model ###########
###### Split into 3 age classes - preschool being less than 5 (p), 
###### schoolage being 5-15 years (s) and adults being >15 years (a) ####
###### Input varying CDR data from India, plot the prevelance or incidence of infection
#####################################################

require(deSolve)
require(graphics)
require(ggplot2)
require(xlsx)

source("basicFunctions_TBModel_3AgeClasses.R")

nstart=c(c(1000000/3,1000000/3,1000000/3),c(1000000/3,1000000/3,1000000/3), rep(0, times=3), rep(0, times=3), rep(1/3, times=3), rep(0, times=3), rep(0, times=3) )
time<-seq(from=1, to=365*1000, by=1 )
output_pre<-as.data.frame(ode(nstart,time,TBmodel_3AgeClasses,parameters))
eqbm<-output_pre[dim(output_pre)[1], ]
eqbm



inputfile<-"/Users/Arathi/Documents/2018/Monash/TBmodelWithDemography/WHO Global TB report 2017/TB_burden_countries_2017-11-09.csv"
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

output_all<-numeric()
delta_all<-numeric()

output_all<-rbind(output_all, output[(dim(output_pre)[1]-365):dim(output_pre)[1], ])
output_all$time<-seq(from=1949, to=1950, by=1/365)

output_all<-output_all[-dim(output_all)[1],]

for(year in 1950:2016){ #2016
  nstart<-as.numeric(eqbm[2:22])
  #time<-seq(from=1,to=(365), by=1)
  year_seq<-seq(from=year, to=(year+1), by=1/365)
  params<-changeCDR_parameter(parameters, year, cdr_data_India)
  #params["delta"]<-0
  delta_all<-c(delta_all, params["delta"])
  output<-as.data.frame(ode(nstart,year_seq,TBmodel_3AgeClasses,params))
  eqbm<-output[dim(output)[1], ]
  
  output<-output[-dim(output)[1],] # remove last as it will be starting values for next year
  output_all<-rbind(output_all, output)
  
  eqbm
  
  if(year==2016){
    quartz()
    par(mfrow=c(2,1), oma=c(0,0,2,0))
    N<-rowSums(output_all[,2:22])
    tot_infected<-rowSums(output_all[,14:19]) #includes infected and treated
    plot(output_all$time, (tot_infected/N)*100000, type='l', lwd=2, col=6, lty=2, xlab="years", ylab="Prev", cex.axis=1.5, cex.lab=1.5) #ylim=c(0, max(output_all$I/N)) , 
   
    plot(cdr_data_India[,1], cdr_data_India[,2],type='b', pch=15, col="brown", cex.axis=1.5, cex.lab=1.5, xlab="years", ylab="cdr")
    
    # legend(x=2000, y=max(output_all$I/N), legend=c("I"), lwd=2, lty=1, col=c(6), cex=1.5)
    
  }
  
  
  
}
