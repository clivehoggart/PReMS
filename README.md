# PRReMS
# Installation
#Required packages parallel, glmnet and pROC can be installed in R from CRAN mirror, eg
install.packages('parallel',repos='http://cran.ma.imperial.ac.uk')

#Save BayesLogit_quiet.tar.gz and prrems.R to directory of choice and install from R
install.packages('path/BayesLogit_quiet.tar.gz', repos = NULL, type="source")

#Install PRReMS
source('path/prrems.R')

# Example analysis of the SPECTF Heart Data Set
#Download SPECTF.train and SPECTF.test from https://archive.ics.uci.edu/ml/machine-learning-databases/spect/

# Read and format data
test <- read.csv('SPECTF.test')

train <- read.csv('SPECTF.train')

x.test <- as.matrix(test[,-1])

x.train <- as.matrix(train[,-1])

y.test <- test[,1]

y.train <- train[,1]

colnames(x.train) <- paste("X",1:44,sep='')

colnames(x.test) <- paste("X",1:44,sep='')

# Estimate the Guassian precison paramter -- the penalty
tau <-  TauEst(  y=y.train, x=x.train, nfolds=50, lasso.penalty="lambda.min" )

# Set PRReMS parameters
no.cores <- 12 # Set according to the number of available processors
k.max <- 10 # Maximum model size to explore
max.s <- 100 # Number of models to take forward at each increase in model size

# PRReMS analysis
prrems.fit <- prrems( y=y.train, x=x.train, family='binomial', tau=tau$tau, k.max=k.max, max.s=max.s, standardize=TRUE, no.cores=no.cores )
prrems.fit <- fill.ICs( fitted.models=prrems.fit, y=y.train, x=x.train, n.waic=10000, model.sizes=1:k.max, no.cores=no.cores )

# Evaluate AUC for best model of each size
r.test <- my.auc( prrems.fit, 1:k.max, x.test, y.test )
t(sapply(r.test,getElement,'ci'))

# Croos-validation to determine optimal model size
cv.prrems.fit <- cv.prrems( y=y.train, x=x.train, k.min=1, k.max=k.max, no.cores=no.cores, nfolds=20, max.s=100 )

# Predict in test set using model size determined by cross-validation
pred.test <- predict.prrems( prrems.fit, x.test, size=cv.prrems.fit$best, no.cores=no.cores )
roc( y.test, pred.test[,1], ci=TRUE )
