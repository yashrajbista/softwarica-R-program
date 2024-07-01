
#Yash Raj Bista : College ID - 240224

#Libraries used in the assignment
install.packages('matlib')#library for Linear Algebra and Multivariate Statistics
library(matlib)

install.packages("ggplot2") #Library for creating time series plot
library(ggplot2)

install.packages("rsample")
library(rsample)

# Three CSV Files are provided where X.csvfile contains 4 input EEG signals
# X = {x1, x2, x3, x4}, y.csv file contains the output EEG signal y
# and time.csv contains the sampling time of all EEG data in seconds
#--------------------------------------------------------------------------------
#Import data (X.csv)
Xdata = (read.csv(file = "C:/Users/Yash Raj/Downloads/ds.csv", header = T))
#Convert XData to Matrix
X = as.matrix(Xdata)
# Set Column names of X to X1, X2, X3 and X4
colnames(X) = c("X1", "X2", "X3", "X4", "X5")
#--------------------------------------------------------------------------------
#Import data (time.csv)
TimeData = (read.csv(file = "C:/Users/Yash Raj/Downloads/ds.csv", header = T))
#Convert TimeData to Matrix
time = as.matrix(TimeData)

#Set Column names to Time in Time
#colnames(Time) = c("time")
#--------------------------------------------------------------------------------
#Task 1: Preliminary data analysis
#--------------------------------------------------------------------------------

#Time series plots (of input and output EEG signals)
XSignals = ts(X) #ts function create time series object
plot(XSignals, main = "Time Series Plot: All input Signals", xlab = "time", ylab =" Signals")

#YSignals = ts(Y) #ts function create time series object
#plot(YSignals, main = "Time Series Plot of Output Signals", xlab = "time", ylab =" Output Signals")

#--------------------------------------------------------------------------------
#Task 1: Distribution of Each EEG Signal
#--------------------------------------------------------------------------------

#All Input gene
par(mfrow=c(1,2)) #plots 2x2=4 plots in same display
#Density plot of all input genes
plot(density(X), type = "l", xlab = "Input Genes", main = "Density curve of all input Genes")
#Creating histogram of X Signal
hist(X, freq = FALSE, xlab = "Input Genes", main = "Histogram and density curve of all input Genes") #plot Histogram of X
#Plot density in histogram
lines(density(X), col = "red", lwd = 1) #Plot density

par(mfrow=c(2,2)) #plots 2x2=4 plots in same display
#Creating histogram of X1 Signal
hist(X[,"X1"], freq = FALSE, xlab = "X1", main = "Histogram and density curve of X1 Genes")
#Plot density in histogram
lines(density(X[,"X1"]), lwd = 1, xlab = "X1")

#Creating histogram of X2 Signal
hist(X[,"X2"], freq = FALSE, xlab = "X2", main = "Histogram and density curve of X2 Genes")
#Plot density in histogram
lines(density(X[,"X2"]), lwd = 1, xlab = "X2")

#Creating histogram of X3 Signal
hist(X[,"X3"], freq = FALSE, xlab = "X3", main = "Histogram and density curve of X3 Genes")
#Plot density in histogram
lines(density(X[,"X3"]), lwd = 1, xlab = "X3")

#Creating histogram of X4 Signal
hist(X[,"X4"], freq = FALSE, xlab = "X4", main = "Histogram and density curve of X4 Genes")
#Plot density in histogram
lines(density(X[,"X4"]), lwd = 1, xlab = "X4")

#Creating histogram of X4 Signal
hist(X[,"X5"], freq = FALSE, xlab = "X5", main = "Histogram and density curve of X5 Genes")
#Plot density in histogram
lines(density(X[,"X5"]), lwd = 1, xlab = "X5")

#Output Signals
par(mfrow=c(1,2)) #plots 1x2 plots in same display
#Density plot of Y signals
#plot(density(Y), type = "l", xlab = "Y", main = "Density curve of Y Signal")
#Creating histogram of Y Signal
#hist(Y, freq = FALSE, xlab = "Y", main = "Histogram and density curve of Y signal") #plot Histogram of X
#Plot density in histogram
lines(density(Y), col = "black", lwd = 1) #Plot density


#--------------------------------------------------------------------------------
#Task 1: Correlation and scatter plots to examine their dependencies 
#--------------------------------------------------------------------------------
par(mfrow=c(2,2)) #plots 2x2=4 plots in same display
#Creating Plot
plot(X[,"X1"], X[,"X2"], pch = 1, col="blue", main = "Correlation between X1 and X2", xlab = "X1")
#regression line
abline(lm(X[,"X2"]~X[,"X1"]), col = "red", lwd = 1)

#Creating Plot
#plot(X[,"X2"], Y, pch = 1, col="blue", main = "Correlation between X2 and Y", xlab = "X2")
#regression line
#abline(lm(Y ~ X[,"X2"]), col = "red", lwd = 1)

#Creating Plot
plot(X[,"X3"], X[,"X2"], pch = 1, col="blue", main = "Correlation between X3 and X2", xlab = "X3")
#regression line
abline(lm(X[,"X2"] ~ X[,"X3"]), col = "red", lwd = 1)

#Creating Plot
plot(X[,"X4"], X[,"X2"], pch = 1, col="blue", main = "Correlation between X4 and X2", xlab = "X4")
#regression line
abline(lm(X[,"X2"] ~ X[,"X4"]), col = "red", lwd = 1)

#Creating Plot
plot(X[,"X5"], X[,"X2"], pch = 1, col="blue", main = "Correlation between X5 and X2", xlab = "X5")
#regression line
abline(lm(X[,"X2"] ~ X[,"X5"]), col = "red", lwd = 1)

#--------------------------------------------------------------------------------
#Task 2:  Regression – modeling the relationship between EEG signals
#--------------------------------------------------------------------------------
#onesMatrix = matrix(1 , length(X)/4,1) # Creating matrix of ones
onesMatrix = matrix(1, nrow(X), 1)



# Creating Data for model 1 from given equation
XData_model1 = cbind(onesMatrix, X[,"X4"], X[,"X3"]^2)
XData_model1
# Creating Data for model 2 from given equation
XData_model2 = cbind(onesMatrix, X[,"X4"], X[,"X3"]^2, X[,"X5"])
XData_model2
# Creating Data for model 3 from given equation
XData_model3 = cbind(onesMatrix, (X[,"X3"]), (X[,"X4"]), X[,"X5"]^3)
XData_model3
# Creating Data for model 4 from given equation
XData_model4 = cbind(onesMatrix, X[,"X4"], (X[,"X3"])^2, (X[,"X5"])^3)
XData_model4
# Creating Data for model 5 from given equation
XData_model5 = cbind(onesMatrix, (X[,"X4"]), (X[,"X1"])^2,  (X[,"X3"])^2)
XData_model5
#--------------------------------------------------------------------------------
#Task 2.1: 
#--------------------------------------------------------------------------------
#Calculation the least square (thetaHat)
#for model1
model1_Thetahat = solve(t(XData_model1) %*% XData_model1) %*% t(XData_model1) %*% X[,"X2"]
model1_Thetahat

# Model 2
model2_Thetahat = solve(t(XData_model2) %*% XData_model2) %*% t(XData_model2) %*% X[,"X2"]
model2_Thetahat

# Model 3
model3_Thetahat = solve(t(XData_model3) %*% XData_model3) %*% t(XData_model3) %*% X[,"X2"]
model3_Thetahat

# Model 4
model4_Thetahat = solve(t(XData_model4) %*% XData_model4) %*% t(XData_model4) %*% X[,"X2"]
model4_Thetahat

# Model 5
model5_Thetahat = solve(t(XData_model5) %*% XData_model5) %*% t(XData_model5) %*% X[,"X2"]
model5_Thetahat


#-------------------------------------------------------------------------------- 
#Task 2.2: Model Residual Error
#--------------------------------------------------------------------------------
# Calculating y-hat
# Model 1
model1_YHat = XData_model1 %*% model1_Thetahat
model1_YHat

# Model 2
model2_YHat = XData_model2 %*% model2_Thetahat
model2_YHat

# Model 3
model3_YHat = XData_model3 %*% model3_Thetahat
model3_YHat

# Model 4
model4_YHat = XData_model4 %*% model4_Thetahat
model4_YHat

# Model 5
model5_YHat = XData_model5 %*% model5_Thetahat
model5_YHat


# Calculating RSS
# Model 1
RSS_model1 = sum((X[,"X2"]-model1_YHat)^2)
RSS_model1

# Model 2
RSS_model2 = sum((X[,"X2"]-model2_YHat)^2)
RSS_model2

# Model 3
RSS_model3 = sum((X[,"X2"]-model3_YHat)^2)
RSS_model3

# Model 4
RSS_model4 = sum((X[,"X2"]-model4_YHat)^2)
RSS_model4

# Model 5
RSS_model5 = sum((X[,"X2"]-model5_YHat)^2)
RSS_model5

#-------------------------------------------------------------------------------- 
#Task 2.3: Calculating Likelihood and Variance for all models
#--------------------------------------------------------------------------------
n = length(X[,"X2"]) #Calculating length of Y
# Calculating Variance for
# Model 1
VAR_model1 = RSS_model1/(n-1)
VAR_model1

# Model 2
VAR_model2 = RSS_model2/(n-1)
VAR_model2

# Model 3
VAR_model3 = RSS_model3/(n-1)
VAR_model3

# Model 4
VAR_model4 = RSS_model4/(n-1)
VAR_model4

# Model 5
VAR_model5 = RSS_model5/(n-1)
VAR_model5

# Calculating likelihood for
# Model 1
Likelihood_model1 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model1))-(1/(2*VAR_model1))*RSS_model1
Likelihood_model1

# Model 2
Likelihood_model2 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model2))-(1/(2*VAR_model2))*RSS_model2
Likelihood_model2

# Model 3
Likelihood_model3 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model3))-(1/(2*VAR_model3))*RSS_model3
Likelihood_model3

# Model 4
Likelihood_model4 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model4))-(1/(2*VAR_model4))*RSS_model4
Likelihood_model4

# Model 5
Likelihood_model5 = -(n/2)*(log(2*pi))-(n/2)*(log(VAR_model5))-(1/(2*VAR_model5))*RSS_model5
Likelihood_model5


#-------------------------------------------------------------------------------- 
#Task 2.4: Calculating AIC and BIC for all models
#--------------------------------------------------------------------------------
# Calculating AIC
# Model 1
AIC_model1 = 2*(length(model1_Thetahat))-2*Likelihood_model1
AIC_model1

# Model 2
AIC_model2 = 2*(length(model2_Thetahat))-2*Likelihood_model2
AIC_model2

# Model 3
AIC_model3 = 2*(length(model3_Thetahat))-2*Likelihood_model3
AIC_model3

# Model 4
AIC_model4 = 2*(length(model4_Thetahat))-2*Likelihood_model4
AIC_model4

# Model 5
AIC_model5 = 2*(length(model1_Thetahat))-2*Likelihood_model5
AIC_model5

# Calculating BIC
# Model 1
BIC_model1 = length(model1_Thetahat)*log(n)-2*Likelihood_model1
BIC_model1

# Model 2
BIC_model2 = length(model2_Thetahat)*log(n)-2*Likelihood_model2
BIC_model2

# Model 3
BIC_model3 = length(model3_Thetahat)*log(n)-2*Likelihood_model3
BIC_model3

# Model 4
BIC_model4 = length(model4_Thetahat)*log(n)-2*Likelihood_model4
BIC_model4

# Model 5
BIC_model5 = length(model5_Thetahat)*log(n)-2*Likelihood_model5
BIC_model5


#-------------------------------------------------------------------------------- 
#Task 2.5: Calculating Error for all models and Plotting Q-Q plot with Q-Q line for them
#--------------------------------------------------------------------------------

par(mfrow = c(3, 2))
# Model 1
Error_model1 = X[,"X2"] - model1_YHat #Error
qqnorm(Error_model1, col = "#336600", main = "Q-Q plot of Model 1") # Plots Graph
qqline(Error_model1, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 2
Error_model2 = X[,"X2"] - model2_YHat #Error
qqnorm(Error_model2, col = "#336600", main = "Q-Q plot of Model 2") # Plots Graph
qqline(Error_model2, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 3
Error_model3 = X[,"X2"] - model3_YHat #Error
qqnorm(Error_model3, col = "#336600", main = "Q-Q plot of Model 3") # Plots Graph
qqline(Error_model3, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 4
Error_model4 = X[,"X2"] - model4_YHat #Error
qqnorm(Error_model4, col = "#336600", main = "Q-Q plot of Model 4") # Plots Graph
qqline(Error_model4, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

# Model 5
Error_model5 = X[,"X2"] - model5_YHat #Error
qqnorm(Error_model5, col = "#336600", main = "Q-Q plot of Model 5") # Plots Graph
qqline(Error_model5, col = "#e60000", lwd = 1) # Adds Q-Q line on graph

#-------------------------------------------------------------------------------- 
#Task 2.6: Selecting a model
#--------------------------------------------------------------------------------
# Calculating Mean Squared Error (MSE)
# Model 1
MSE_model1 = sum(Error_model1^2)/length(Error_model1)
MSE_model1

# Model 2
MSE_model2 = sum(Error_model2^2)/length(Error_model2)
MSE_model2

# Model 3
MSE_model3 = sum(Error_model3^2)/length(Error_model3)
MSE_model3

# Model 4
MSE_model4 = sum(Error_model4^2)/length(Error_model4)
MSE_model4

# Model 5
MSE_model5 = sum(Error_model5^2)/length(Error_model5)
MSE_model5



#-------------------------------------------------------------------------------- 
#Task 2.7: Splitting input and output data sets into training and testing sets in ratio of 7:3
#--------------------------------------------------------------------------------
Y = X[,"X2"]
# Splitting the data (Training Data set)
XSplit = initial_split(data = as.data.frame(X), prop = .7)
YSplit = initial_split(data = as.data.frame(Y), prop = .7)
# Training Data set
# Y Data
Y_Training_Set = training(YSplit) # Y Training data set
Y_Training_Data = as.matrix(Y_Training_Set) # Y Training data
# X Data
X_Training_Set = training(XSplit) # X Training data set
X_Training_Data = as.matrix(X_Training_Set) # X Training data

# Testing Data set
# Y Data
Y_Testing_Set = testing(YSplit) # Y Testing data set
Y_Testing_Data = as.matrix(Y_Testing_Set) # Y Testing data
# X Data
X_Testing_Set = testing(XSplit) # X Testing data set
X_Testing_Data = as.matrix(X_Testing_Set) # X Testing data

# Selecting Model 5, estimating model parameters using training data
TrainingOneXMatrix = matrix(1, length(X_Training_Set$X1), 1) # ones matrix for training set
TrainingOneXMatrix = matrix(1, length(X_Training_Set$X1), 1)
n = length(X_Training_Set$X1)
TrainingOneXMatrix = matrix(1, n, 1)
TrainingXModel = cbind(TrainingOneXMatrix, X_Training_Set[, "X4"], (X_Training_Set[, "X1"])^2, (X_Training_Set[, "X3"])^2) #Training Model
TrainingThetaHat = solve(t(TrainingXModel) %*% TrainingXModel) %*% t(TrainingXModel) %*% Y_Training_Data
TrainingThetaHat

# Computing output/prediction of Model2 using testing data set
TestingYHat = X_Testing_Data %*% TrainingThetaHat
TestingYHat
RSStesting = sum((Y_Testing_Set - TestingYHat)^2)
RSStesting 

#used t-test as sample size is less than 30 i.e 4 and varriance is know (calculated)
t.test(Y_Training_Data, mu = 500, alternative = "two.sided", conf.level = 0.95)
C_I1 = -0.06331197
C_I2 = 0.51992124
meu = 0.2283046

# With 95% of confidence interval, predicting the model and plotting them with testing data and error bars
par(mfrow = c(1, 1))
TrainingDensity = density(Y_Training_Data) # Density of training data of output signal
TrainingDensity
plot(TrainingDensity, col="#336600", lwd = 2, main="Distribution of Output Signal Training Data")
abline(v = C_I1, col = "#e60000", lty=2)
abline(v = C_I2, col = "#e60000", lty=2)
abline(v = meu, col = "#1a1a1a", lty=2)

residual = ((Y_Testing_Set - TestingYHat)) # Calculating Error
residual

# plotting Error Bars
# Calculating Standard Deviation (Sigma)
Sigma = sqrt(VAR_model2) #Variance of model2 from task 2.3
Sigma
XData_model2 #Data model 2 from task 2.1

dataFrame = data.frame(
  xAxis = XData_model2,
  yAxis = Y
)
dataFrame

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.1, y=Y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.1, ymin=Y1-Sigma, ymax=Y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 2 - X1)", x="Model 2 - X1", y = "Output Signal Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.2, y=Y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.2, ymin=Y1-Sigma, ymax=Y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 2 - X2)", x="Model 2 - X2", y = "Output Signal Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.3, y=Y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.3, ymin=Y1-Sigma, ymax=Y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 2 - X3)", x="Model 2 - X3", y = "Output Signal Data")

ggplot(dataFrame) +
  geom_bar( aes(x=xAxis.4, y=Y1), stat="identity", fill="#336600", alpha=0.7) +
  geom_errorbar( aes(x=xAxis.4, ymin=Y1-Sigma, ymax=Y1+Sigma), width=0.4, colour="#e60000", alpha=0.9, linewidth=1) +
  labs(title="Error Bar (Model 2 - X4)", x="Model 2 - X4", y = "Output Signal Data")


#-------------------------------------------------------------------------------- 
#Task 3: Approximate Bayesian Computation (ABC)
#--------------------------------------------------------------------------------
array1 = 0
array2 = 0
f_value = 0
s_value = 0

# Model 5 thetahat values from Task 2.1
ThetaBias = 1.2951518 # chosen parameter
ThetaA = 0.5385828 # chosen parameter
ThetaB = 0.8312983 # set constant
ThetaC = 0.1096679 # set constant
Epsilon = RSS_model5 * 2 ## fixing value of epsilon, RSS_model5 from task 2.2
num = 100 # number of iteration
##Calculating Y-hat for performing rejection ABC
counter <- 0
for (i in 1:num) {
  range1 = runif(1, -1.2951518, 1.2951518) # calculating the range
  range1
  range2 = runif(1, -0.5385828, 0.5385828)
  range2
  NewThetahat = matrix(c(range1, range2, ThetaB, ThetaC))
  NewYHat = XData_model5 %*% NewThetahat ## New Y hat and model2_Thetahat from task2.1
  NewRSS = sum((Y - NewYHat)^2)
  NewRSS
  if (NewRSS > Epsilon){ #Performing rejection ABC
    array1[i] = range1
    array2[i] = range2
    counter = counter + 1
    Fvalue = matrix(array1)
    Svalue = matrix(array2)
  }
}
# Plotting the graph
plot(Fvalue, Svalue, col = c("#ff1a1a", "#3366ff"), main = "Joint and Marginal Posterior Distribution Model 2")

