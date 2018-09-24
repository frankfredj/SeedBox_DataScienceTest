#-------- WARNING --------

#Custom packages are used in this file, namely:
#ggplot2, foreach, DoParallel, reshape2, cowplot, grid, ggfortify

#-------- WARNING --------



#Load files

library(readr)

testSamples <- read_csv("C:/Users/Francis/Downloads/testSamples.csv")
transData <- read_csv("C:/Users/Francis/Downloads/transData.csv")

#Sort the data

Id <- testSamples$sample_id

Transaction_group <- unique(intersect(testSamples$sample_id, transData$sample_id))
No_Transaction_group <- testSamples$sample_id[-which(!is.na(match(testSamples$sample_id,Transaction_group)))]

Transaction_test <- Transaction_group[which(testSamples$test_group[which(!is.na(match(testSamples$sample_id,Transaction_group)))] == 1)]
Transaction_control <- Transaction_group[which(testSamples$test_group[which(!is.na(match(testSamples$sample_id,Transaction_group)))] == 0)]

#Build analysis functions to be run in loops

#This function pulls out transactions associated with an ID

pull_transactions <- function(id){

if(length(which(Transaction_group == id)) != 0){

tS.index <- which(testSamples$sample_id == id)
tD.index <- which(transData$sample_id == id)

output <- transData[tD.index,]

return(output)}

else {print("This client did not perform any transaction", quote = FALSE)}

}

#This function gives basic descriptive statistics based off of a pull_transactions output

analyse_transaction <- function(frame){

output <- matrix(nrow = 1, ncol = 9)
colnames(output) <- c("TrNum","TrTotal" , "TrMean", "TrMedian", "TrSD", "REBILL", "CHARGEBACK", "REFUND", "RF + ChB")

output[1] <- nrow(frame)
output[2] <- sum(frame$transaction_amount)
output[3] <- mean(frame$transaction_amount)
output[4] <- median(frame$transaction_amount)
output[5] <- sd(frame$transaction_amount)
output[6] <- length(which(frame$transaction_type == "REBILL"))
output[7] <- length(which(frame$transaction_type == "CHARGEBACK"))
output[8] <- length(which(frame$transaction_type == "REFUND"))
output[9] <- output[7] + output[8]

if(is.na(output[5])){output[5] <- 0}

return(output)

}

#This function is a combination of the previous 2

pull_n_analyse <- function(id){return(analyse_transaction(pull_transactions(id)))}

#Run parallel loops to extract data from both groups (test, control)

library(foreach)
library(doParallel)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

Transaction_test_analysis <- foreach(i = 1:length(Transaction_test), .combine = "rbind") %dopar% pull_n_analyse(Transaction_test[i])

stopCluster(cl)


cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

Transaction_control_analysis <- foreach(i = 1:length(Transaction_control), .combine = "rbind") %dopar% pull_n_analyse(Transaction_control[i])

stopCluster(cl)



library(ggplot2)

#This function extract the frequency (in %) of the different values of a pull_n_analyse output
#This will be used to graph numerical pdfs

Extract_Frequencies <- function(frame){

output <- list()

for(i in 1:ncol(frame)){

x <- sort(unique(frame[,i]))
y <- x
for(j in 1:length(y)){y[j] <- length(which(frame[,i] == x[j]))}
y <- y / nrow(frame)
y <- as.matrix(y)
colnames(y) <- colnames(frame)[i]
rownames(y) <- c(x)

output[[i]] <- y

}

return(output)

}

Transaction_control_frequencies <- Extract_Frequencies(Transaction_control_analysis)
Transaction_test_frequencies <- Extract_Frequencies(Transaction_test_analysis)

library(reshape2)
library(grid)
library(cowplot)

#This function takes a frequency-based output (from Transaction_control_frequencies and Transaction_test_frequencies)
#It spits out a data frame with percentile-based linear interpolations so to be able to compare both the control and test respective outputs of the Extract_Frequencies function

Graph_bind <- function(k){

x <- Transaction_control_frequencies[[k]]
y <- Transaction_test_frequencies[[k]]

x.p <- as.numeric(rownames(x))
y.p <- as.numeric(rownames(y))

p <- seq(from = min(c(x.p,y.p)), to = max(c(x.p,y.p)), length.out = 100)

x.i <- approx(x.p, x, p, method = "linear")$y
y.i <- approx(y.p, y, p, method = "linear")$y

if(length(which(is.na(x.i))) != 0){x.i[which(is.na(x.i))] <- 0}
if(length(which(is.na(y.i))) != 0){y.i[which(is.na(y.i))] <- 0}

plot.frame <- as.data.frame(cbind(p,x.i,y.i))

colnames(plot.frame) <- c("Value", paste("control", colnames(Transaction_control_frequencies[[k]])), paste("test", colnames(Transaction_control_frequencies[[k]])))

return(plot.frame)

}

stack_freq_plot <- function(index){

plots <- list()

for(i in 1:length(index)){

graph_data <- melt(Graph_bind(index[i]), id = "Value")

plots[[i]] <- ggplot(data=graph_data,
       aes(x=Value, y=value, colour=variable)) +
       geom_line() + ylab("%")


}


return(plots)

}

#As a result, we can easily draw comparative numerical pdfs graph with 4 lines of code, i.e:

plots.1 <- stack_freq_plot(c(1,2,3))
plots.2 <- stack_freq_plot(c(6,7,8,9))

plot_grid(plots.1[[1]], plots.1[[2]], plots.1[[3]], labels = c("Number of bills", "Total Revenue" ,"Mean Transaction Value"))
plot_grid(plots.2[[1]], plots.2[[2]], plots.2[[3]], plots.2[[4]], labels = c("REBILL", "CHARGEBACK", "REFUND", "ChB + RF"))


#The refund and chargeback distributions are nearly identical
#High rebill rates appear to be much more frequent in the control group

#Low number of bills seem to be more frequent in the test group
#Negative transaction values appear to be much more frequent in the control group


total.frame <- rbind(Transaction_control_analysis, Transaction_test_analysis)

identifier <- matrix(nrow = nrow(total.frame), ncol = 1)
identifier[1:nrow(Transaction_control_analysis)] <-  "Control"
identifier[-c(1:nrow(Transaction_control_analysis))] <- "Test"
identifier <- as.data.frame(identifier)
colnames(identifier) <- c("Group")

identifier$Group <- as.factor(identifier$Group)
levels(identifier$Group) <- make.names(levels(factor(identifier$Group)))

library(ggfortify)

autoplot(prcomp(scale(total.frame)), data = identifier, colour = 'Group')

#The data isn't linearly separable, and there aren't enough predictor to use kernels or other fancy separation technique
#This also mean that linear models will be unreliable

#-----------------------------------------------------

#What is the approximate probability distribution between the two group?

total.matrix <- total.frame[,c(2,6,9)]

autoplot(prcomp(total.matrix), data = identifier, colour = 'Group')

total.Tmat <- as.matrix(eigen(var(total.matrix))$vectors)

total.PC <- total.matrix %*% total.Tmat

colnames(total.PC) <- colnames(total.matrix)

library(caret)

TrainingParameters.Prob <-  trainControl(method = "cv", number = 7, verboseIter = TRUE, allowParallel = TRUE, savePredictions = TRUE, summaryFunction = twoClassSummary, classProbs = TRUE)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

multinom.Model <- train(total.frame, identifier$Group,
                  method = "multinom",
                  trControl= TrainingParameters.Prob,
                  tuneLength = 10,
                  metric = "ROC"

)

multinomPredictions <- predict(object = multinom.Model , total.frame)
# Create confusion matrix
CF <-confusionMatrix(multinomPredictions, identifier$Group)

print(CF)


plot(varImp(multinom.Model))


#There doesn't seem to be a clear way to distinguish between both groups based on data alone
#I will use Statistical Simulation techniques to simulate from the two distributions {test, control}

#I want to simulate data that satisfy three conditions:

# 1: means of columns are preserved
# 2: sd of columns are preserved
# 3: the variance-covariance is preserved

#To do this, I will apply a linear transformation (PCA) so to obtain non-correlated columns
#I will simulate data from the empirical pdf function of these columns
#I will apply an inverse transformation (to go back from PCA to the original data structure)

#Using such a tool, I will be able to run statistical tests with regards to questions 2,3 and 4

#I will only consider three variables: {total revenue, REBILL and (REFUND + CHARGEBACK)}


test.matrix <- Transaction_test_analysis[,c(2,6,9)]
control.matrix <- Transaction_control_analysis[,c(2,6,9)]

test.Tmat <- as.matrix(eigen(var(test.matrix))$vectors)
control.Tmat <- as.matrix(eigen(var(control.matrix))$vectors)

test.InvMat <- solve(test.Tmat)
control.InvMat <- solve(control.Tmat)

test.PC <- test.matrix %*% test.Tmat
control.PC <- control.matrix %*% control.Tmat

colnames(test.PC) <- colnames(test.matrix)
colnames(control.PC) <- colnames(control.matrix)

test.PC.freq <- Extract_Frequencies(test.PC)
control.PC.freq <- Extract_Frequencies(control.PC)

#This function takes an output from the Extract_Frequencies function, and spits out a variable in accordance to the empirical cdf (based off the frequencies)

Simulate_from_frequency <- function(freq,k,p){

y <- freq
x <- as.numeric(rownames(y))

n <- length(y)

output <- matrix(nrow = k, ncol = 1)


for(i in 1:k){


repeat{

u <- c(ceiling(runif(min = 1, max = n, n = 1)), runif(min = 0, max = 1, n = 1))

if(u[2] <= y[u[1]]){break}

}

output[i] <- x[u[1]]


}

colnames(output) <- colnames(freq)
return(output)

}

#To show that we are indeed simulating RVs from the exact same distribution as our sample, let's do a test run:

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

col.1 <- foreach(i = 1:10000, .combine = "c") %dopar% Simulate_from_frequency(test.PC.freq[[1]],10)

stopCluster(cl)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

col.2 <- foreach(i = 1:10000, .combine = "c") %dopar% Simulate_from_frequency(test.PC.freq[[2]],10)

stopCluster(cl)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

col.3 <- foreach(i = 1:10000, .combine = "c") %dopar% Simulate_from_frequency(test.PC.freq[[3]],10)

stopCluster(cl)



Simulated.Data <- (as.matrix(cbind(col.1, col.2, col.3)) %*% test.InvMat)
colnames(Simulated.Data) <- colnames(test.matrix)

print(cor(Simulated.Data))
print(cor(test.matrix))

print("Respective means:")
for(i in 1:ncol(Simulated.Data)){print(c(mean(Simulated.Data[,i]),mean(test.matrix[,i])))}

print("Respective sd's:")
for(i in 1:ncol(Simulated.Data)){print(c(sd(Simulated.Data[,i]),sd(test.matrix[,i])))}

#So we are effectively able to sample from any numerical multivariate distribution with mean vector, sd vector and vcov matrix {u, o, M} using the above tool
#This will be usefull to run statistical simulations rather than drawing conclusion from a sample which size is relatively modest

#This will also enable us to generate RVs {PR1, PR2, PR3} which then can be fed to a matrix that will generate {Revenue, REBILL, CHARGEBACK + REFUND}
#From the increment's artifical cdf, we will fit distributions for {test, control}

#There is, however, a probability of not performing any transactions at all.

p.transaction <- matrix(nrow = 2, ncol = 1)
colnames(p.transaction) <- c("P(Transaction)")
rownames(p.transaction) <- c("Test", "Control")

p.transaction[1] <- length(Transaction_test) / length(which(testSamples$test_group == 1))
p.transaction[2] <- length(Transaction_control) / length(which(testSamples$test_group == 0))

print(p.transaction)



#Additionnaly, we'll build a function that generates two whole vectors from two numerical distributions, then compare their difference
#Arguments p.1 and p.2 represent the probability of generating a transaction

Simulate_Diff <- function(freq.list.1, freq.list.2, index, mat.1, mat.2, p.1, p.2){

u <- matrix(nrow = length(freq.list.1), ncol = 1)
v <- u

p.vec <- runif(min = 0, max = 1, n = 2)

if(p.vec[1] <= p.1){for(i in 1:length(u)){u[i] <- Simulate_from_frequency(freq.list.1[[i]],1)}} else {u[] <- 0}
if(p.vec[2] <= p.2){for(i in 1:length(v)){v[i] <- Simulate_from_frequency(freq.list.2[[i]],1)}} else {v[] <- 0}

u <- mat.1 %*% u
v <- mat.2 %*% v

return(u[index] - v[index])

}



#---------------------------------------
#Code for computing p-values
#---------------------------------------



cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

ReBill.Differences <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(test.PC.freq, control.PC.freq, 2, test.InvMat, control.InvMat, p.transaction[1], p.transaction[2])

stopCluster(cl)


library(DescTools)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

ReBill.Differences.control <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(control.PC.freq, control.PC.freq, 2, control.InvMat, control.InvMat, p.transaction[2], p.transaction[2])

stopCluster(cl)


GrtThenOneReBill.Ps <- matrix(nrow = 2, ncol = 1)
colnames(GrtThenOneReBill.Ps) <- "p"
rownames(GrtThenOneReBill.Ps) <- c("test vs control", "control vs control")
GrtThenOneReBill.Ps[1] <- length(which(ReBill.Differences >= 1)) / length(ReBill.Differences)
GrtThenOneReBill.Ps[2] <- length(which(ReBill.Differences.control  >= 1)) / length(ReBill.Differences.control)


print(GrtThenOneReBill.Ps)

p.val.test <- matrix(nrow = 1, ncol = 1)
colnames(p.val.test) <- c("Binomial p-value")

x <- GrtThenOneReBill.Ps[1] - GrtThenOneReBill.Ps[2]
sigma <- GrtThenOneReBill.Ps[1]*(1-GrtThenOneReBill.Ps[1])/length(ReBill.Differences.control) + GrtThenOneReBill.Ps[2]*(1-GrtThenOneReBill.Ps[2])/length(ReBill.Differences.control)

p.val.test[1] <- pnorm(x / sqrt(sigma))

print(p.val.test)




cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

Revenue.Differences <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(test.PC.freq, control.PC.freq, 1, test.InvMat, control.InvMat, p.transaction[1], p.transaction[2])

stopCluster(cl)

cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

Revenue.Differences.control <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(control.PC.freq, control.PC.freq, 1, control.InvMat, control.InvMat, p.transaction[2], p.transaction[2])

stopCluster(cl)


MoreRev.Ps <- matrix(nrow = 2, ncol = 1)
colnames(MoreRev.Ps) <- "p"
rownames(MoreRev.Ps) <- c("test vs control", "test vs test")
MoreRev.Ps[1] <- length(which(Revenue.Differences > 0)) / length(Revenue.Differences)
MoreRev.Ps[2] <- length(which(Revenue.Differences.control  > 0)) / length(Revenue.Differences.control)


print(MoreRev.Ps)

MoreRev.p.val.test <- matrix(nrow = 1, ncol = 1)
colnames(MoreRev.p.val.test) <- c("Binomial p-value")

x <- MoreRev.Ps[1] - MoreRev.Ps[2]
sigma <- MoreRev.Ps[1]*(1-MoreRev.Ps[1])/length(ReBill.Differences.control) + MoreRev.Ps[2]*(1-MoreRev.Ps[2])/length(ReBill.Differences.control)

MoreRev.p.val.test[1] <- pnorm(x / sqrt(sigma))

print(MoreRev.p.val.test)




test.matrix.2 <- Transaction_test_analysis[,c(2,6,7,8)] 
control.matrix.2 <- Transaction_control_analysis[,c(2,6,7,8)]

test.matrix.2[,2] <- test.matrix.2[,3] / test.matrix.2[,2]
control.matrix.2[,2] <- control.matrix.2[,3] / control.matrix.2[,2]

test.matrix.2 <- test.matrix.2[,-3]
control.matrix.2 <- control.matrix.2[,-3]

test.matrix.2[which(is.infinite(test.matrix.2))] <- 1
control.matrix.2[which(is.infinite(control.matrix.2))] <- 1

test.matrix.2[which(is.na(test.matrix.2))] <- 0
control.matrix.2[which(is.na(control.matrix.2))] <- 0

test.Tmat.2 <- as.matrix(eigen(var(test.matrix.2))$vectors)
control.Tmat.2 <- as.matrix(eigen(var(control.matrix.2))$vectors)

test.InvMat.2 <- solve(test.Tmat.2)
control.InvMat.2 <- solve(control.Tmat.2)

test.PC.2 <- test.matrix.2 %*% test.Tmat.2
control.PC.2 <- control.matrix.2 %*% control.Tmat.2

colnames(test.PC.2) <- colnames(test.matrix.2)
colnames(control.PC.2) <- colnames(control.matrix.2)

test.PC.freq.2 <- Extract_Frequencies(test.PC.2)
control.PC.freq.2 <- Extract_Frequencies(control.PC.2)



cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

ChargeBackRates.Differences <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(test.PC.freq.2, control.PC.freq.2, 2, test.InvMat.2, control.InvMat.2, p.transaction[1], p.transaction[2])

stopCluster(cl)


cl <- makeCluster(detectCores() - 2) 
clusterEvalQ(cl, library(foreach))
registerDoParallel(cl)

ChargeBackRates.Differences.control <- foreach(i = 1:100000, .combine = "c") %dopar% Simulate_Diff(control.PC.freq.2, control.PC.freq.2, 2, control.InvMat.2, control.InvMat.2, p.transaction[2], p.transaction[2])

stopCluster(cl)


MoreCBr.Ps <- matrix(nrow = 2, ncol = 1)
colnames(MoreCBr.Ps) <- "p"
rownames(MoreCBr.Ps) <- c("test vs control", "test vs test")
MoreCBr.Ps[1] <- length(which(ChargeBackRates.Differences >= 0)) / length(ChargeBackRates.Differences)
MoreCBr.Ps[2] <- length(which(ChargeBackRates.Differences.control  >= 0)) / length(ChargeBackRates.Differences.control)


print(MoreCBr.Ps)

MoreCBr.p.val.test <- matrix(nrow = 1, ncol = 1)
colnames(MoreCBr.p.val.test) <- c("Binomial p-value")

x <- MoreCBr.Ps[1] - MoreCBr.Ps[2]
sigma <- MoreCBr.Ps[1]*(1-MoreCBr.Ps[1])/length(ReBill.Differences.control) + MoreCBr.Ps[2]*(1-MoreCBr.Ps[2])/length(ReBill.Differences.control)

MoreCBr.p.val.test[1] <- pnorm(x / sqrt(sigma))

print(MoreCBr.p.val.test)







