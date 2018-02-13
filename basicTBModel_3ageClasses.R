###### Create an age stratified model ###########
###### Split into 3 age classes - preschool being less than 5 (p), 
###### schoolage being 5-15 years (s) and adults being >15 years (a) ####
#####################################################

require(deSolve)
require(graphics)
require(ggplot2)
require(xlsx)

source("basicFunctions_TBModel_3AgeClasses.R")

nstart=c(c(1000000/3,1000000/3,1000000/3),c(1000000/3,1000000/3,1000000/3), rep(0, times=3), rep(0, times=3), rep(1/3, times=3), rep(0, times=3), rep(0, times=3) )
time<-seq(from=1, to=365*1000, by=1 )
output<-as.data.frame(ode(nstart,time,TBmodel_3AgeClasses,parameters))
eqbm<-output[dim(output)[1], ]
eqbm

#plot the run in output from each of the age classes
plot_output_3ageClasses(output)

#plot the total/aggregate output from all compartments
plot_aggregatedoutput_3ageClasses(output)

nstart<-as.numeric(eqbm[2:22])
parameters$delta<-0.59*73/(0.41*365*70*3) # 59% treatment
time<-seq(from=1, to=365*100, by=1 )
output<-as.data.frame(ode(nstart,time,TBmodel_3AgeClasses,parameters))

#plot the run in output from each of the age classes
plot_output_3ageClasses(output)

#plot the total/aggregate output from all compartments
plot_aggregatedoutput_3ageClasses(output)

