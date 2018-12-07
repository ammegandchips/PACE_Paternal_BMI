#weighted mean and standard deviation of total maternal and paternal BMI and meta-analysed correlation coefficient
mat <- data.frame(es = c(22.60, 26.17, 26.66, 26.48, 23.24, 23.84, 23.95, 24.29, 24.02, 22.56, 24.33, 25.09),
sd= c(3.34, 5.76, 6.10, 4.51, 3.86, 4.51, 4.08, 4.56, 3.94, 3.88, 4.93, 5.47),
n= c(531,70, 115, 158, 947, 352, 982, 621, 212, 98, 324, 94))
pat <- data.frame(es = c(34.98, 26.57, 27.45, 28.04, 25.17, 25.82, 25.61, 25.74, 25.98, 24.93, 26.33, 27.23),
sd = c(3.01, 5.25, 4.73, 4.18, 3.21, 3.54, 3.09, 3.06, 3.10, 3.03, 3.58, 3.95),
n= c(531,70, 115, 158, 947, 352, 982, 621, 212, 98, 324, 94))
mat$se <- mat$sd/sqrt(mat$n)
pat$se <- pat$sd/sqrt(pat$n)
require(metafor)
mat.estimate<-rma.uni(yi=mat$es,sei=mat$se,method="FE",weighted=TRUE)#23.70
mat.estimate$se*sqrt(sum(mat$n))#standard deviation = 4.11
pat.estimate<-rma.uni(yi=pat$es,sei=pat$se,method="FE",weighted=TRUE)#26.98
pat.estimate$se*sqrt(sum(pat$n))#standard deviation = 3.24

correlation <- c(0.21,0.19, 0.35, 0.11,0.17,0.22,0.13,0.22,0.25,0.30, 0.25, 0.26)
Ps <- c(7.6e-7,1.1e-1,9.1e-4,1.1e-1,8.8e-8,4.7e-5,6.2e-6,2.2e-8,2.0e-4,2.6e-3,5.5e-6,1.1e-2)
T <- qt(Ps,mat$n-1,lower=F)
SEs <- correlation/T
rma.uni(yi=correlation,sei=SEs,method="FE",weighted=TRUE)
