source("myread.ENVI.R")
source("imagematrix.R")
require(ggplot2)
require(reshape2)
require(ggthemes)
require(maxLik)
imagepath <- "../Statistics-SAR-Intensity-master/Data/Images/ESAR/"
HH_Complex <- myread.ENVI(paste(imagepath, "ESAR97HH.DAT", sep = ""),
                          paste(imagepath, "ESAR97HH.hdr", sep = ""))
HH_Intensity <- (Mod(HH_Complex))^2
example <- HH_Intensity[1300:1400,2280:2480]
vec_example <- data.frame(HH=as.vector(example))

##select forest region
plot(imagematrix(equalize(example)))
imagematrixPNG(name = "./forest.png", imagematrix(equalize(example)))

vec_example <- data.frame(HH=as.vector(example))
summary(vec_example)

binwidth_complete <- 2*IQR(vec_example$HH)*length(vec_example$HH)^(-1/3)

ggplot(data=vec_example, aes(x=HH)) + 
  geom_histogram(aes(y=..density..), 
  binwidth = binwidth_complete) + 
  xlab("Intensities") +
  ylab("Proportions") +
  ggtitle("Complete Histogram") +
  theme_few()
ggsave(filename = "./HistogramExample.pdf")

## Estimation

GI0.Estimator.m1m2 <- function(z, L) {
  m1 <- mean(z)
  m2 <- mean(z^2)
  m212 <- m2/m1^2
  
  a <- -2 - (L+1) / (L * m212)
  g <- m1 * (2 + (L+1) / (L * m212))
  
  return(list("alpha"=a, "gamma"=g))
}

estim.example <- GI0.Estimator.m1m2(example, 1)

LogLikelihoodLknown <- function(params) {
  
  p_alpha <- -abs(params[1])
  p_gamma <- abs(params[2])
  p_L <- abs(params[3])
  
  n <- length(z)
  
  return(
    n*(lgamma(p_L-p_alpha) - p_alpha*log(p_gamma) - lgamma(-p_alpha)) + 
      (p_alpha-p_L)*sum(log(p_gamma + z*p_L)) 
  )
}

z <- vec_example$HH

estim.exampleML <- maxNR(LogLikelihoodLknown, 
                       start=c(estim.example$alpha, estim.example$gamma,1), 
                       activePar=c(TRUE,TRUE,FALSE))$estimate[1:2]

