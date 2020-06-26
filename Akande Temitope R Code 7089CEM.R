
library(ggplot2)
library(ggthemes)

gene <- read.csv("gene_data.csv")

str(gene)

################################################################
# Plot of Time Series
################################################################
df = data.matrix(gene)
#select each columns
xt = df[, 1]
x1 = df[, 2]
x2 = df[, 3]
x3 = df[, 4]
x4 = df[, 5]
x5 = df[, 6]


plot.ts(df, main = "Time series plot of Gene Dataset", las = 1,)

###################################################################
#Gene Distributions
##################################################################
summary(gene)
#the function to create the histogram of each gene with the normal distribution superimposed over it
distribution.plot <- function(g) {
  m <- mean(g) #mean
  std <- sqrt(var(g)) #Standard deviation
  h <- hist(
    g,
    breaks = 10,
    density = 10,
    col = "blue",
    xlab = "gene",
    main = "Histogram with Normal Curve",
    las = 1
  )
  xfit <- seq(min(g), max(g), length = 20)
  yfit <- dnorm(xfit, mean = mean(g), sd = sd(g))
  yfit <- yfit * diff(h$mids[1:2]) * length(g)
  
  lines(xfit, yfit, col = "black", lwd = 2)
}


#histograms Gene 1 -5
distribution.plot(gene$x1)
distribution.plot(gene$x2)
distribution.plot(gene$x3)
distribution.plot(gene$x4)
distribution.plot(gene$x5)

##################################################################
#correlation matrix
#################################################################
cor(gene, method = "pearson")
cor(gene, method = "spearman")
cor.test(gene$x3, gene$x5, method = "pearson")
cor.test(gene$x3, gene$x5, method = "spearman")

#scatterplot matrix between all the gene
pairs( ~ x1 + x2 + x3 + x4 + x5, data = gene, main = "scatterplot matrix")

plot(
  gene$x3,
  gene$x5,
  xlabel = "Gene3",
  ylabel = "Gene5",
  las = 1,
  col = 'blue',
  main = "Scatterplot of Gene4 vs Gene5"
)
#abline(lm(gene$x3~gene$x5), col = 'red')
lines(smooth.spline(gene$x3, gene$x5),
      lty = 2,
      lwd = 5)

###############################################################
#Eigan Decomposition PCA
###############################################################

pca <- prcomp(gene[c(2:6)], center = TRUE, scale = TRUE)

#standard deviation, proportion of variance and cumulative proportion
summary(pca)
expl.var <- round(pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100)

#function to create colors
Cols = function (vec) {
  cols = rainbow (length (unique (vec)))
  return (cols[as.numeric (as.factor (vec))])
}


#PCA Screenplot
screeplot(pca, type = "l", main = "Screeplot of the 5 PCs")
abline(h = 1, col = "red", lty = 5)
legend(
  "topright",
  legend = c("Eigenvalue = 1"),
  col = c("red"),
  lty = 5,
  cex = 0.6
)

#cumulative variance plot
cumpro <- cumsum(pca$sdev ^ 2 / sum(pca$sdev ^ 2))
plot(cumpro[0:5],
     xlab = "PC #",
     ylab = "Amount of explained variance",
     main = "Cumulative variance plot")
abline(v = 2, col = "blue", lty = 5)
abline(h = 0.975, col = "blue", lty = 5)
legend(
  "topleft",
  legend = c("Cut-off @ PC2"),
  col = c("blue"),
  lty = 5,
  cex = 0.6
)

#biplot of PCA
biplot(pca, cex = 0.5, cex.axis = 0.5)

#Scatterplot of PC1 and PC2
plot(
  pca$x[, 1],
  pca$x[, 2],
  las = 1,
  col = Cols(gene),
  pch = 19,
  xlab = "PC1 (70%)",
  ylab = "PC2 (27%)",
  main = "PC1 / PC2 - plot"
)
########################################################
# (MSE and residual analysis for Linear model)
DF = gene[c(2:6)]
#DF = data.frame(data, check.name = FALSE)

x1 <- DF$x4
x2 <- DF$x5
y <- DF$x3
ones = matrix(1, length(y), 1)

Y = matrix(y)
X_line = cbind(ones, x2)
theta_hat = solve(t(X_line) %*% X_line) %*% t(X_line) %*% y
print(theta_hat)

y_hat = X_line %*% theta_hat

par(mfrow = c(1, 1))
plot(
  x2,
  y,
  type = "p",
  col = "green",
  main = "Linear Regression Model",
  xlab = "X value",
  ylab = "Y value"
)
lines(x2, y_hat, col = "red", lw = 3)

MSE = mean((y_hat - y) ^ 2)
print(MSE)

# residual analysis
## residual scatter plot
res = y - y_hat
par(mfrow = c(1, 2))
plot(res, type = "p", main = "Residual Scatter Plot")

qqnorm(res, main = "Residual Q-Q Normal Plot")
qqline(res, col = "red", lw = 3)

##########################################################
# Non-linear polynomial regression 
###########################################################
X = cbind(ones, x1, x1^2, x1^3, x1^4,x2, x2^2, x2^3, x2^4)

colnames(X) = c("intercept", "a1", "a2", "a3", "a4", "b1", "b2", "b3", "b4")

set.seed(1000)
train = sample(nrow(X), 0.8 * nrow(X))
train
X_train = X[train, ]
X_test = X[-train, ]

Y_train = Y[train, ]
Y_test = Y[-train, ]

select_eval = function(X_tr, Y_tr = Y_train, X_ts, Y_ts = Y_test){
  theta = solve(t(X_tr) %*% X_tr) %*% t(X_tr) %*% Y_tr
  y_hat = X_ts %*% theta
  #calculate MSE
  MSE = mean((Y_ts - y_hat)^2)
  return (MSE)
}
selection_1 = data.frame(term = colnames(X), MSE = rep(0, 9))
for (i in 1:9){
  selection_1[i, 2] = select_eval(X_tr = X_train[,i], X_ts = X_test[,i])
}
index = which.min(selection_1$MSE)
best1 = (as.character(selection_1$term[which.min(selection_1$MSE)]))
cat("Best model selection 1:", paste0("y~", best1, "+error"),"(MSE=", selection_1$MSE[index],")")

final_X_tr = X_train[,index]
final_X_ts = X_test[,index]
X2_train = X_train[,-index]
X2_test = X_test[,-index]

selection_2 = data.frame(term = paste0(best1, "+", colnames(X2_train)), MSE = rep(0, 8))
for (i in 1:8){
  selection_2[i, 2] = select_eval(X_tr = cbind(X2_train[,i], final_X_tr), X_ts = cbind(X2_test[,i], final_X_ts))
}
index = which.min(selection_2$MSE)
best2 = (as.character(selection_2$term[index]))
cat("\nBest model selection 2:", paste0("y~", best2, "+ error"),"(MSE=", selection_2$MSE[index],")")

final_X_tr = cbind(final_X_tr, X2_train[,index])
final_X_ts = cbind(final_X_ts, X2_test[,index])
X3_train = X2_train[,-index]
X3_test = X2_test[,-index]

selection_3 = data.frame(term = paste0(best2, " + ", colnames(X3_train)), MSE = rep(0,7))

for (i in 1:7) {
  selection_3[i,2] = select_eval(X_tr = cbind(X3_train[,i],final_X_tr), X_ts = cbind(X3_test[,i], final_X_ts))
}

index = which.min(selection_3$MSE)
best3 = (as.character(selection_3$term[which.min(selection_3$MSE)]))
cat("\nBest model selection 3:", paste0("y ~ ",best3, " + error"), "(MSE=", selection_3$MSE[index], ")")

final_X_tr = cbind(final_X_tr, X3_train[,index])
final_X_ts = cbind(final_X_ts, X3_test[,index])
X4_train = X3_train[,-index]
X4_test = X3_test[,-index]

selection_4 = data.frame(term = paste0(best3, " + ", colnames(X4_train)), MSE = rep(0,6))

for (i in 1:6) {
  selection_4[i,2] = select_eval(X_tr = cbind(X4_train[,i],final_X_tr), X_ts = cbind(X4_test[,i], final_X_ts))
}

index = which.min(`selection_4`$MSE)
best4 = (as.character(selection_4$term[which.min(selection_4$MSE)]))
cat("\nBest model selection 4:", paste0("y ~ ",best4, " + error"), "(MSE=", selection_4$MSE[index], ")")

final_X_tr = cbind(final_X_tr, X4_train[,index])
final_X_ts = cbind(final_X_ts, X4_test[,index])
X5_train = X4_train[,-index]
X5_test = X4_test[,-index]

selection_5 = data.frame(term = paste0(best4, " + ", colnames(X5_train)), MSE = rep(0,5))

for (i in 1:5) {
  selection_4[i,2] = select_eval(X_tr = cbind(X5_train[,i],final_X_tr), X_ts = cbind(X5_test[,i], final_X_ts))
}

index = which.min(selection_5$MSE)
best5 = (as.character(selection_5$term[which.min(selection_5$MSE)]))
cat("\nBest model selection 5:", paste0("y ~ ",best5, " + error"), "(MSE=", selection_5$MSE[index], ")")


model <- lm(y ~ I(x2^3) + I(x1^4) + I(x2^4) + I(x1^3), data = DF)
summary(model)
plot(model)
#Covariance Matrix
vcov(model)


#AIC selection
select_eval = function(X, Y){
  X = as.matrix(X)
  Y = as.matrix(Y)
  theta = solve(t(X) %*% X) %*% t(X) %*% Y
  y_hat2 = X %*% theta
  #error
  residuals = Y - y_hat2
  sigma_sq = sum(residuals^2)/(nrow(Y) - 1)
  loglik= sum(log(dnorm(DF$x3, mean = y_hat2, sd = sqrt(sigma_sq))))
  
  k = ncol(X)
  AIC = 2*k - 2*loglik
  return(AIC)  
}  
parameters = 1:9
param_1 = t(as.matrix(parameters))
param_2 = combn(parameters, m = 2)
param_3= combn(parameters, m = 3)

# 1 parameter
AIC_1 = data.frame(combination = rep(0, ncol(param_1)),
                   AIC = rep(0, ncol(param_1)))
for (i in 1:ncol(param_1)) {
  x = X[,param_1[,i]]
  AIC_1[i, 1] = colnames(X)[i]
  AIC_1[i, 2] = select_eval(X, Y)
}

# 2 parameters
AIC_2 = data.frame(combination = rep(0, ncol(param_2)),
                   AIC = rep(0, ncol(param_2)))
for (i in 1:ncol(param_2)) {
  x = X[,param_2[,i]]
  AIC_2[i, 1] = paste(colnames(x), collapse = " + ")
  AIC_2[i, 2] = select_eval(X, Y)
}

# 3 parameters
AIC_3 = data.frame(combination = rep(0, ncol(param_3)),
                   AIC = rep(0, ncol(param_3)))
for (i in 1:ncol(param_3)) {
  x = X[,param_3[,i]]
  AIC_3[i, 1] = paste(colnames(x), collapse = " + ")
  AIC_3[i, 2] = select_eval(X, Y)
}
AIC_results = rbind(AIC_1, AIC_2, AIC_3)

best_AIC = AIC_results$combination[which.min(AIC_results$AIC)]
min_AIC = AIC_results$AIC[which.min(AIC_results$AIC)]
min_AIC



#model validation
X = cbind(x2^3,
          x1^4,
          x2^4,
          x1^3
)

set.seed(1000)
train = sample(nrow(X), 0.7 * nrow(X))

X_train = X[train, ]
X_test = X[-train, ]

Y_train = Y[train, ]
Y_test = Y[-train, ]
theta = solve(crossprod(X_train), crossprod(X_train, Y_train))
y_hat3 = X_train %*% theta

residuals = (Y_train - y_hat3)
SSE = sum(residuals^2)
sigma_sq = SSE/(nrow(X_train) - 1)
R2_adj = function(x, y, p) {
  n = length(x)
  R2 = cor(x, y) ^ 2
  R2_adj = 1 - (1 - R2)*((n - 1)/(n - p - 1))
  return(R2_adj)
}

r_adjusted = R2_adj(x = Y_train, y = y_hat3, p = 3)
residuals_pl = as.data.frame(residuals)


cov = sigma_sq * (solve(t(X_train) %*% X_train))
print(round(cov, 6))
#Models prediction
conf = {}
for (i in 1:nrow(X_train)) {
  Xi = X_train[i,]
  Xi = matrix(Xi, 1, 4)
  v = Xi %*% cov %*% t(Xi)
  conf = rbind(conf, v)
}
conf_95 = 1.96*sqrt(conf)

#Plot the predictions + Confidence interval
pred_data = data.frame(y_hat = y_hat3, lower = y_hat3 - conf_95, upper = y_hat3 + conf_95)
pred_data = cbind(X_train[,1],Y_train, pred_data)
colnames(pred_data)[1:2] = c("x", "y")
(error_pred = ggplot(pred_data, aes(x = x, y = y_hat3))+
    geom_point(color = "blue")+
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)+
    ylab("y")+
    theme_few())
(y_y_pred = ggplot(pred_data, aes(x = y, y = y_hat3))+
    geom_point()+
    theme_few())
#model validation
#split into 70:30
train = DF[1:211,]
test = DF[212:301,]
X_train = cbind(train$x5^3,
                train$x4^4,
                train$x5^4,
                train$x4^3
)
X_test = cbind(test$x5^3,
               test$x4^4,
               test$x5^4,
               test$x4^3
)
fit_evaluate = function(X_tr, y_tr, X_ts, y_ts){
  theta = solve(t(X_tr) %*% X_tr) %*% t(X_tr) %*% y_tr
  y_hat4 = X_ts %*% theta
  #calculate error
  MSE_test = mean((y_ts - y_hat4)^2)
  R2_test = R2_adj(y_ts, y_hat4, 4)
  y_hat4 = X_tr %*% theta
  #calculate error
  MSE_train = mean((y_tr - y_hat4)^2)
  R2_train = R2_adj(y_tr, y_hat4, 4)
  return(list(MSE_train, MSE_test, R2_train, R2_test))
}

fit = fit_evaluate(X_tr = X_train, y_tr = train$x3, X_ts = X_test, y_ts = test$x3)

cat("MSE on training data = ", fit[[1]],
    "\nMSE on testing data = ", fit[[2]],
    "\nR2 on training data = ", fit[[3]],
    "\nR2 on testing data = ", fit[[4]])

k = 10
folds = cut(seq_len(nrow(DF)),breaks=k,labels=FALSE)
#Perform 10 fold cross validation
result = data.frame(fold = 1:k,
                    MSE_train = rep(0, k),
                    MSE_test = rep(0, k),
                    R2_train = rep(0, k),
                    R2_test = rep(0, k)
)


for(i in 1:k){
  #Segment your data by fold using the which() function
  test_idx = which(folds==i,arr.ind=TRUE)
  train = DF[-test_idx,]
  test = DF[test_idx,]
  X_train = cbind(train$x5^3,
                  train$x4^4,
                  train$x5^4,
                  train$x4^3
  )
  X_test = cbind(test$x5^3,
                 test$x4^4,
                 test$x5^4,
                 test$x4^3
  )
  fit = fit_evaluate(X_tr = X_train, y_tr = train$x3,
                     X_ts = X_test, y_ts = test$x3)
  result[i, 2:5] = fit
}

plot_MSE = data.frame(mean = c(mean(result$MSE_train),
                               mean(result$MSE_test)),
                      sd = c(sd(result$MSE_train),
                             sd(result$MSE_test)),
                      type = c("training", "testing"))
R2_plot = data.frame(mean = c(mean(result$R2_train),
                              mean(result$R2_test)),
                     sd = c(sd(result$R2_train),
                            sd(result$R2_test)),
                     type = c("training", "testing"))

plot_MSE$type = factor(plot_MSE$type, levels = c("training", "testing"))
R2_plot$type = factor(R2_plot$type, levels = c("training", "testing"))
(plot_MSE = ggplot(result, aes(MSE_train, MSE_test))+
    geom_line()+
    theme_few())
# 
#rm(list = ls())

###########################################

