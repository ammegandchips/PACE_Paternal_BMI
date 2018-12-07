#weighted mean and standard deviation of total maternal and paternal BMI and meta-analysed correlation coefficient
mat <- data.frame(es = c(22.63,26.99,24.05,23.88,24.55,24.45),
sd= c(33.35,4.43,3.87,4.44,5.14,4.73),
n= c(570,108,335,516,177,276))
pat <- data.frame(es = c(25.01,27.63,25.01,26.58,25.99,26.40),
sd = c(3.11,3.43,3.28,3.71,3.64,3.73),
n= c(570,108,335,516,177,276))
mat$se <- mat$sd/sqrt(mat$n)
pat$se <- pat$sd/sqrt(pat$n)
require(metafor)
mat.estimate<-rma.uni(yi=mat$es,sei=mat$se,method="FE",weighted=TRUE)
mat.estimate$se*sqrt(sum(mat$n))#standard deviation
pat.estimate<-rma.uni(yi=pat$es,sei=pat$se,method="FE",weighted=TRUE)
pat.estimate$se*sqrt(sum(pat$n))#standard deviation

correlation <- c(0.19,0.21,0.22,0.28,0.29,0.42)
Ps <- c(5.1e-6,2.9e-2,1.5e-5,1.1e-10,9.6e-5,4.2e-13)
T <- qt(Ps,mat$n-1,lower=F)
SEs <- correlation/T
rma.uni(yi=correlation,sei=SEs,method="FE",weighted=TRUE)
