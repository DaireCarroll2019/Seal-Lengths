#Daire Carroll, Gothenburg University, 2021
#script to develop some potential demographic scenarios for an observed populaiton growht rate under the assumption of assymptotic growth
#for harbour seals aquired from https://sharkweb.smhi.se/hamta-data/

graphics.off() 
setwd("C:/Users/daire/Desktop/Drones")

my.data = read.csv("Shark_2_Max_Means_Pups.csv", sep = ",", header = TRUE, fileEncoding = 'UTF-8-BOM')
attach(my.data)
head(my.data)

my.data$Urs_Rams_78 = Urs_Rams_78/0.65
my.data$Urs_Rams_88 = Urs_Rams_88/0.65
my.data$Urs_Rams_02 = Urs_Rams_02/0.65

attach(my.data)

nlc = nls.control(maxiter = 100000,minFactor = 1/100024)

nls1 = nls(Urs_Rams_88  ~ a*exp(log(r)*Year), 
           start = list(a = 7e-82, r = exp(0.0968)),control = nlc)
nls4 = nls(Pups_88  ~ a*exp(log(r)*Year), 
           start = list(a = 7e-27, r = exp(0.0332)),control = nlc)
nls2 = nls(Urs_Rams_02  ~ a*exp(log(r)*Year), 
           start = list(a = 1e-22, r = exp(0.03)),control = nlc)
nls3 = nls(Urs_Rams_78  ~ a*exp(log(r)*Year), 
           start = list(a = 7e-27, r = exp(0.0332)),control = nlc)
nls5 = nls(Pups_02~ a*exp(log(r)*Year), 
          start = list(a = 2e17, r = exp(-0.017)),control = nlc )

ratio1 = (Pups_88)[12:25] / (Urs_Rams_88)[12:25]
ratio2 = (Pups_02)[26:length(Year)]/ (Urs_Rams_02)[26:length(Year)]
ratio3 = ((summary(nls4)$par[1]*exp(summary(nls4)$par[2]*Year[12:25])))/((summary(nls1)$par[1]*exp(summary(nls1)$par[2]*Year[12:25])))
ratio4 = ((summary(nls5)$par[1]*exp(summary(nls5)$par[2]*Year[36:length(Year)])))/((summary(nls2)$par[1]*exp(summary(nls2)$par[2]*Year[36:length(Year)])))

ratio5= (Pups_88)[12:25] / ((summary(nls1)$par[1]*exp(summary(nls1)$par[2]*Year[12:25])))
ratio6 = (Pups_02)[36:length(Year)]/ ((summary(nls2)$par[1]*exp(summary(nls2)$par[2]*Year[36:length(Year)])))

m1 = (summary(nls1)$par[1]*exp(summary(nls1)$par[2]*Year[12:25]))
m2 = (summary(nls2)$par[1]*exp(summary(nls2)$par[2]*Year[27:length(Year)]))
m3 = (summary(nls3)$par[1]*exp(summary(nls3)$par[2]*Year[1:11]))
m4 = (summary(nls4)$par[1]*exp(summary(nls4)$par[2]*Year[12:25]))
m5 = (summary(nls5)$par[1]*exp(summary(nls5)$par[2]*Year[36:length(Year)]))
m1_b = data.frame(cbind(m1,Year[12:25]))
m2_b = data.frame(cbind(m2,Year[27:length(Year)]))
m3_b = data.frame(cbind(m3,Year[1:11]))
m4_b = data.frame(cbind(m4,Year[12:25]))
m5_b = data.frame(cbind(m5,Year[36:length(Year)]))

ratio1 = data.frame(cbind(ratio1,Year[12:25]))
ratio2 = data.frame(cbind(ratio2,Year[26:length(Year)]))
ratio3 = data.frame(cbind(ratio3,Year[12:25]))
ratio4 = data.frame(cbind(ratio4,Year[36:length(Year)]))
ratio5 = data.frame(cbind(ratio5,Year[12:25]))
ratio6 = data.frame(cbind(ratio6,Year[36:length(Year)]))

N = 1 #pop size
A = 38 #maximum age (https://doi.org/10.1002/ecs2.3343)
class = as.factor(c(1:A))

B = rep(0,A)
B[4] = 0.17
B[5] = 0.33
B[6:27] = 0.47
B[27:(A)] = 0.35 #birthrates (https://doi.org/10.1002/ecs2.3343)

S = rep(0,A)
S[1] = 0.75
S[2:4] = 0.89
S[5:(A)] = 0.95 #survival (https://doi.org/10.1002/ecs2.3343)

B = 2*B/S #age specific fertility is back calculated

ev_goal = coef(nls2)[2] #the target eiganvalue
ev_1 = round(ev_goal*1000)/1000  #target based on exponential fit

structure = function(B,S,N){
  L = matrix((B*S/2),ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev = as.numeric(ev$vectors[,1])
  P = c()
  for(i in 1:A){
    P[i] = ev[i]/sum(ev)*N
  } 
  return(P)
} #function to give eigan vector given B/S and a structure if N != 1 

structure_ev = function(B,S){
  L = matrix((B*S/2),ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev = as.numeric(ev$value[1])
  return(ev)
} #function to give eigan value

structure_ev(B,S)
structure(B,S,1)

eig_val = matrix(ncol = 4)

for(i in 250:1000){
  s0 = i*0.001 
  for(j in 500:1000){ 
    s_sub_A = j*0.001 
    for(k in 500){
      b_A = k*0.001
      B.1 = B/2
      S.1 = S
      S.1[1] = s0
      S.1[2:4] = s_sub_A
      vals = c(structure_ev(B.1,S.1),s0,s_sub_A,b_A)
      eig_val = rbind(eig_val,vals)
    }
  }
} #eigan values for a range S & B which give target lambda

eig_val = eig_val[2:nrow(eig_val),]
colnames(eig_val) = c("Eiganvalue","Pup","Sub-adult","Fertility_Reduction")
eig_val = data.frame(eig_val)
eig_val$Eiganvalue = round(eig_val$Eiganvalue*1000)/1000

eig_val_pup = eig_val[eig_val$Eiganvalue == ev_1 & 
                      eig_val$Pup < eig_val$Sub.adult &
                      eig_val$Fertility_Reduction == 0.5
                      ,] #pup and subadult survival values which will give you a lambda of approximatly 1.028

ev_goal = round(ev_goal*100)/100

parm1 = function(par){
  
  S1 = S
  
  S1[1] = par[1]
  
  S1[2:4] = par[2]
  
  B1 = B
  
  B1 = B1*S1/2
  
  B1 = B1/2 #year skipping
  
  L = matrix(B1,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S1[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev2 = as.numeric(ev$value[1])
  err = sqrt((ev_goal-ev2)^2)
  return(err)
} #function given vector par containing pup and subadult survval, returns error based on distance from ev_goal

o1 = optim(c(0.75, 0.89),
           parm1,
           method = "L-BFGS-B",
           lower = c(0,0), 
           upper = c(095,0.95)) #can't be higher than adult survival

o1$par #one possible option, gives the minimum error but given variance of the data used to establish expected lambda many options for pup and subadult survival are valid

S1 = S
S1[1] = o1$par[1]
S1[2:4] = o1$par[2]

B1 = B/2

structure_ev(B1,S1)
  
Sen0 = structure(B1,S1,1)

Sen0 = data.frame(cbind(class,Sen0))#age structure

plot(Sen0$Sen0~Sen0$class,
     xlab = "Age Class",
     ylab = "Proportion")
