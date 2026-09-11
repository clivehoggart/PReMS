# PReMS

PReMS (Parallel Regularised Model Search) searches for sparse predictive models by combining regularised regression with a structured search over model space. The current implementation is supplied as an R script (`prems.R`).

For each model size, PReMS ranks candidate models using an approximate marginal-likelihood score and retains the highest-ranking models for expansion to the next model size. Model size can then be chosen by cross-validation.

## Current scope

The current script supports Gaussian, binomial and Cox proportional-hazards outcomes. Candidate models are ranked by the PReMS model-search criterion (`criteria = "ML"`), and model size can be selected using `cv.prems()`.

For binary outcomes, `y` must be coded `0/1`. For Cox models, `y` should be a `survival::Surv` object; both ordinary right-censored `Surv(time, status)` responses and counting-process `Surv(start, stop, status)` responses are supported. Predictor matrices should be numeric, contain no missing values and have unique column names. Prediction data must use the same predictor names as the training data.

## Installation

PReMS is currently used by sourcing the R script rather than installing an R package.

The only external package required for the core Gaussian and binomial workflow is `glmnet`. Cox models additionally require `survival`.

```r
install.packages("glmnet")

# Required only for Cox models if survival is not already installed:
install.packages("survival")

source("prems.R")
```

`parallel` is a recommended R package distributed with R and does not normally need to be installed separately. PReMS calls `parallel::mclapply()` explicitly. On Windows use `no.cores = 1`; multiple cores can be used with `mclapply()` on macOS and Linux.

The SPECTF example uses `pROC` to calculate AUC:

```r
install.packages("pROC")
```

The diabetes and NKI70 data are downloaded directly in the examples below, so the `lars` and `penalized` packages do not need to be installed. `pROC` is otherwise required only if the convenience functions `cv.auc()` or `my.auc()` are used. `gplots`, `igraph`, `MASS` and `exvatools` are not required by the current script.

## Main functions

| Function | Purpose |
| --- | --- |
| `prems()` | Search and fit candidate models up to a specified model size. |
| `ModelSearchIncrease()` | Extend an existing PReMS search by one additional model size without repeating the smaller-model searches. |
| `cv.prems()` | Cross-validation over model size. |
| `TauEst()` | Estimate the ridge-prior precision used by PReMS from a Lasso fit. |
| `getModelFit()` | Extract a ranked fitted model of a specified size. |
| `predict.prems()` | Generate predictions from a fitted PReMS model. |
| `plot.cv.prems()` | Plot the cross-validation criterion against model size. |

The most important arguments are:

| Argument | Meaning |
| --- | --- |
| `x` | Matrix of predictors considered for selection. |
| `x.fixed` | Optional covariates included in every model and not subject to selection. The default is `NULL`. |
| `family` | Outcome family: `"binomial"`, `"gaussian"` or `"cox"`. |
| `tau` | Ridge-prior precision. Usually estimated with `TauEst()`. |
| `k.max` | Largest candidate model size. |
| `max.s` | Number of highest-ranking models retained at each model size for subsequent expansion. |
| `max2way` | `"all"` evaluates all two-predictor models; an integer limits the two-predictor search by expanding only that many top one-predictor models. |
| `standardize` | Standardise predictors internally before model fitting. |
| `no.cores` | Number of cores used by the model search. |

Within the returned fit object, the quantity named `ML` is the negative approximate log marginal likelihood, so smaller values indicate better-ranked models. Cross-validation uses predictive performance to select model size and returns both the size with the best mean performance (`best`) and a smaller one-standard-error alternative where available (`one.se`).

## Example: UCI SPECTF Heart data

The SPECTF Heart data set contains 267 individuals measured for 44 continuous SPECT-derived predictors. UCI supplies predefined training (`SPECTF.train`, 80 individuals) and test (`SPECTF.test`, 187 individuals) sets. The outcome is binary, with `0` and `1` denoting the two diagnostic classes.

Download and extract the SPECTF Heart data from the UCI Machine Learning Repository:

https://archive.ics.uci.edu/dataset/96/spectf+heart

Place `SPECTF.train` and `SPECTF.test` in the working directory.

### 1. Read the data

```r
source("prems.R")

feature.names <- as.vector(rbind(
  paste0("F", 1:22, "R"),
  paste0("F", 1:22, "S")
))

train <- read.csv("SPECTF.train", header = FALSE)
test  <- read.csv("SPECTF.test",  header = FALSE)

colnames(train) <- c("diagnosis", feature.names)
colnames(test)  <- c("diagnosis", feature.names)

y.train <- train$diagnosis
x.train <- as.matrix(train[, feature.names])

y.test <- test$diagnosis
x.test <- as.matrix(test[, feature.names])

table(y.train)
table(y.test)
```

### 2. Select model size by cross-validation

For a direct comparison with Lasso, first create a single set of stratified folds and use these folds for both methods.

```r
foldid <- make.folds2(y.train, folds = 10, seed = 1)

cv.fit <- cv.prems(
  y = y.train,
  x = x.train,
  family = "binomial",
  k.min = 1,
  k.max = 5,
  foldid = foldid,
  max2way = "all",
  max.s = 50,
  no.cores = 4,
  criteria = "ML"
)

cv.fit$best
cv.fit$one.se

plot.cv.prems(cv.fit)
```

For a faster exploratory run, reduce `max.s` or supply an integer such as `max2way = 20`. For higher-dimensional applications, evaluating all two-predictor models can become expensive and a restricted two-predictor search is generally more practical.

### 3. Fit a cross-validated Lasso for comparison

The same fold assignments can be supplied to `glmnet`. The example retains both the minimum-deviance (`Lasso.min`) and one-standard-error (`Lasso.1se`) models.

```r
lasso.cv <- glmnet::cv.glmnet(
  x = x.train,
  y = y.train,
  family = "binomial",
  alpha = 1,
  foldid = foldid,
  type.measure = "deviance",
  standardize = TRUE
)

lasso.coef.min <- as.matrix(coef(lasso.cv, s = "lambda.min"))
lasso.coef.1se <- as.matrix(coef(lasso.cv, s = "lambda.1se"))

lasso.size.min <- sum(lasso.coef.min[-1, 1] != 0)
lasso.size.1se <- sum(lasso.coef.1se[-1, 1] != 0)

c(
  Lasso.min = lasso.size.min,
  Lasso.1se = lasso.size.1se
)
```

Equivalently, coefficients can be obtained directly from the cross-validation object with `coef(lasso.cv, s = "lambda.min")` and `coef(lasso.cv, s = "lambda.1se")`.

### 4. Fit PReMS to the complete training set

Estimate `tau` on the full training data, then repeat the PReMS search using the same maximum model size considered during cross-validation.

```r
set.seed(1)

tau.fit <- TauEst(
  y = y.train,
  x = x.train,
  family = "binomial",
  nfolds = 10
)

tau <- tau.fit$tau.opt

fit <- prems(
  y = y.train,
  x = x.train,
  family = "binomial",
  tau = tau,
  k.max = 5,
  max2way = "all",
  max.s = 50,
  no.cores = 4,
  standardize = TRUE
)
```

### 5. Inspect the selected PReMS model

Here the model size with the best mean cross-validation performance is used. `cv.fit$one.se` can instead be used when a smaller model within the one-standard-error criterion is preferred.

```r
selected.size <- cv.fit$best

selected.model <- getModelFit(
  fit,
  size = selected.size,
  rank = 1,
  criteria = "ML"
)

selected.model$beta
```

The returned coefficient vector contains the intercept (`I`) followed by the selected predictors. Alternative high-ranking models of the same size can be inspected by increasing `rank`.

```r
second.model <- getModelFit(
  fit,
  size = selected.size,
  rank = 2,
  criteria = "ML"
)

second.model$beta
```

### 6. Extend an existing search to a larger model size

`ModelSearchIncrease()` can be used when a fitted search should be extended by one additional model size without repeating the searches already completed. For example, the fit above searched models of sizes 1--5. The following extends it to size 6:

```r
fit <- ModelSearchIncrease(
  fitted.models = fit,
  y = y.train,
  x = x.train,
  no.cores = 4
)

size6.model <- getModelFit(
  fit,
  size = 6,
  rank = 1,
  criteria = "ML"
)

size6.model$beta
```

By default, `ModelSearchIncrease()` uses the number of models retained at the current largest model size to determine the breadth of the expansion. Supply `max.s` explicitly to use a different number. Each call adds one additional model size, so the function can be called repeatedly if a search needs to be extended further.

### 7. Compare PReMS and Lasso in the independent test set

Generate held-out predictions for PReMS, `Lasso.min` and `Lasso.1se` and compare both AUC and model size.

```r
prems.pred <- predict.prems(
  fit,
  newx = x.test,
  size = selected.size,
  rank = 1,
  criteria = "ML"
)

lasso.pred.min <- as.numeric(
  predict(
    lasso.cv,
    newx = x.test,
    s = "lambda.min",
    type = "response"
  )
)

lasso.pred.1se <- as.numeric(
  predict(
    lasso.cv,
    newx = x.test,
    s = "lambda.1se",
    type = "response"
  )
)

prems.auc <- as.numeric(
  pROC::auc(pROC::roc(y.test, as.numeric(prems.pred), quiet = TRUE))
)

lasso.auc.min <- as.numeric(
  pROC::auc(pROC::roc(y.test, lasso.pred.min, quiet = TRUE))
)

lasso.auc.1se <- as.numeric(
  pROC::auc(pROC::roc(y.test, lasso.pred.1se, quiet = TRUE))
)

comparison <- data.frame(
  method = c("PReMS", "Lasso.min", "Lasso.1se"),
  predictors = c(
    selected.size,
    lasso.size.min,
    lasso.size.1se
  ),
  AUC = c(
    prems.auc,
    lasso.auc.min,
    lasso.auc.1se
  )
)

comparison
```

The predefined SPECTF test set is highly imbalanced, so AUC is more informative than raw classification accuracy for this example. The comparison is intended to illustrate the PReMS workflow and its sparsity relative to a standard Lasso analysis rather than to provide a definitive benchmark based on a single small data set.


## Example: Efron et al. diabetes data

The diabetes data used by Efron et al. (2004) provide a compact example with a continuous outcome. The `lars` package contains 442 individuals, a continuous measure of diabetes progression (`y`) and two predictor matrices: `x` contains the 10 original baseline measurements, while `x2` contains 64 main-effect, quadratic and interaction terms. The 64-predictor `x2` matrix is used here because it provides a more informative sparse-model search example.

Unlike SPECTF Heart, these data do not have a predefined training and test split. The example therefore creates an 80% training set and reserves the remaining 20% for independent evaluation.

### 1. Download, load and split the data

The exact `diabetes` object distributed with the `lars` package can be downloaded directly from the read-only CRAN package mirror on GitHub:

```r
diabetes.url <- paste0(
  "https://github.com/cran/lars/raw/refs/heads/master/",
  "data/diabetes.RData"
)

download.file(diabetes.url, "diabetes.RData", mode = "wb")
load("diabetes.RData")

x <- as.matrix(diabetes$x2)
y <- diabetes$y

set.seed(1)
train.id <- sample(seq_len(nrow(x)), size = floor(0.8 * nrow(x)))
test.id <- setdiff(seq_len(nrow(x)), train.id)

x.train <- x[train.id, , drop = FALSE]
y.train <- y[train.id]

x.test <- x[test.id, , drop = FALSE]
y.test <- y[test.id]
```

### 2. Select PReMS model size by cross-validation

Use the same fold assignments for PReMS and Lasso.

```r
foldid <- make.folds.continuous(length(y.train), folds = 10)

cv.fit <- cv.prems(
  y = y.train,
  x = x.train,
  family = "gaussian",
  k.min = 1,
  k.max = 10,
  foldid = foldid,
  max2way = "all",
  max.s = 50,
  no.cores = 4,
  criteria = "ML"
)

cv.fit$best
cv.fit$one.se

plot.cv.prems(cv.fit)
```

With 64 predictors, an exhaustive two-predictor search remains manageable. For substantially larger predictor sets, use an integer value for `max2way` to restrict the initial search.

### 3. Fit Lasso for comparison

```r
lasso.cv <- glmnet::cv.glmnet(
  x = x.train,
  y = y.train,
  family = "gaussian",
  alpha = 1,
  foldid = foldid,
  type.measure = "deviance",
  standardize = TRUE
)

lasso.coef.min <- as.matrix(coef(lasso.cv, s = "lambda.min"))
lasso.coef.1se <- as.matrix(coef(lasso.cv, s = "lambda.1se"))

lasso.size.min <- sum(lasso.coef.min[-1, 1] != 0)
lasso.size.1se <- sum(lasso.coef.1se[-1, 1] != 0)
```

### 4. Refit PReMS to the complete training set

```r
set.seed(1)

tau.fit <- TauEst(
  y = y.train,
  x = x.train,
  family = "gaussian",
  nfolds = 10
)

fit <- prems(
  y = y.train,
  x = x.train,
  family = "gaussian",
  tau = tau.fit$tau.opt,
  k.max = 10,
  max2way = "all",
  max.s = 50,
  no.cores = 4,
  standardize = TRUE
)

selected.size <- cv.fit$best
selected.model <- getModelFit(
  fit,
  size = selected.size,
  rank = 1,
  criteria = "ML"
)

selected.model$beta
```

### 5. Compare prediction in the independent test set

For a continuous outcome, predictive performance can be summarised using root mean squared error (RMSE) and squared Pearson correlation (`R2`).

```r
prems.pred <- as.numeric(
  predict.prems(
    fit,
    newx = x.test,
    size = selected.size,
    rank = 1,
    criteria = "ML"
  )
)

lasso.pred.min <- as.numeric(
  predict(lasso.cv, newx = x.test, s = "lambda.min")
)

lasso.pred.1se <- as.numeric(
  predict(lasso.cv, newx = x.test, s = "lambda.1se")
)

rmse <- function(y, pred) {
  sqrt(mean((y - pred)^2))
}

r2 <- function(y, pred) {
  cor(y, pred)^2
}

comparison <- data.frame(
  method = c("PReMS", "Lasso.min", "Lasso.1se"),
  predictors = c(
    selected.size,
    lasso.size.min,
    lasso.size.1se
  ),
  RMSE = c(
    rmse(y.test, prems.pred),
    rmse(y.test, lasso.pred.min),
    rmse(y.test, lasso.pred.1se)
  ),
  R2 = c(
    r2(y.test, prems.pred),
    r2(y.test, lasso.pred.min),
    r2(y.test, lasso.pred.1se)
  )
)

comparison
```

Because the train/test split is generated randomly, exact results will depend on the split. Setting the random seed makes the example reproducible.

## Example: NKI70 breast-cancer survival data

The `nki70` data supplied with the `penalized` package provide a compact survival example with an omics-like predictor structure. The data contain 144 lymph-node-positive breast-cancer patients, metastasis-free follow-up time, an event indicator, five clinical risk factors and expression measurements for 70 genes. The 70-gene signature was identified in an earlier study, so this example is useful for demonstrating the software workflow but should not be interpreted as an independent biomarker-discovery benchmark.

This data set is particularly useful for illustrating `x.fixed`: the 70 gene-expression measurements can be searched by PReMS while the clinical predictors are included in every model.

### 1. Download and prepare the data

The exact `nki70` object distributed with the `penalized` package can be downloaded directly from the read-only CRAN package mirror on GitHub:

```r
nki70.url <- paste0(
  "https://github.com/cran/penalized/raw/refs/heads/master/",
  "data/nki70.RData"
)

download.file(nki70.url, "nki70.RData", mode = "wb")
load("nki70.RData")

x <- as.matrix(nki70[, 8:77])

x.fixed <- model.matrix(
  ~ ER + Age + Diam + N + Grade,
  data = nki70
)[, -1, drop = FALSE]

y <- survival::Surv(
  time = nki70$time,
  event = nki70$event
)
```

Here `x` contains the 70 candidate molecular predictors and `x.fixed` contains the encoded clinical covariates. The intercept column produced by `model.matrix()` is removed because the Cox model does not require an intercept.

### 2. Select PReMS model size by cross-validation

The intended interface mirrors the binary and Gaussian examples:

```r
foldid <- make.folds2(nki70$event, folds = 10, seed = 1)

cv.fit <- cv.prems(
  y = y,
  x = x,
  x.fixed = x.fixed,
  family = "cox",
  k.min = 1,
  k.max = 5,
  foldid = foldid,
  max2way = "all",
  max.s = 50,
  no.cores = 4,
  criteria = "ML"
)
```

### 3. Refit PReMS to the complete data

After selecting model size, estimate `tau` on the complete data and repeat the search.

```r
tau.fit <- TauEst(
  y = y,
  x = x,
  x.fixed = x.fixed,
  family = "cox",
  nfolds = 10
)

fit <- prems(
  y = y,
  x = x,
  x.fixed = x.fixed,
  family = "cox",
  tau = tau.fit$tau.opt,
  k.max = 5,
  max2way = "all",
  max.s = 50,
  no.cores = 4
)

selected.size <- cv.fit$best

selected.model <- getModelFit(
  fit,
  size = selected.size,
  rank = 1,
  criteria = "ML"
)

selected.model$beta
```

The returned coefficients contain the unpenalised clinical covariates followed by the selected gene-expression predictors.

### 4. Cox Lasso comparison

A corresponding Lasso can be fit with the clinical covariates unpenalised by assigning them a penalty factor of zero:

```r
x.lasso <- cbind(x.fixed, x)
penalty.factor <- c(
  rep(0, ncol(x.fixed)),
  rep(1, ncol(x))
)

y.glmnet <- survival::Surv(nki70$time, nki70$event)

lasso.cv <- glmnet::cv.glmnet(
  x = x.lasso,
  y = y.glmnet,
  family = "cox",
  alpha = 1,
  foldid = foldid,
  penalty.factor = penalty.factor,
  standardize = TRUE
)

coef.min <- as.matrix(coef(lasso.cv, s = "lambda.min"))
coef.1se <- as.matrix(coef(lasso.cv, s = "lambda.1se"))

# Count selected genes only; clinical covariates are unpenalised.
lasso.size.min <- sum(coef.min[(ncol(x.fixed) + 1):nrow(coef.min), 1] != 0)
lasso.size.1se <- sum(coef.1se[(ncol(x.fixed) + 1):nrow(coef.1se), 1] != 0)
```

## Fixed covariates

Covariates supplied through `x.fixed` are included in every model but are not subject to predictor selection. If there are no fixed covariates, leave `x.fixed = NULL` (the default).

For example, if age and sex are to be included in all models:

```r
fit <- prems(
  y = y,
  x = biomarkers,
  x.fixed = cbind(age = age, sex = sex),
  family = "binomial",
  tau = tau,
  k.max = 5,
  max.s = 50,
  no.cores = 4
)
```

Use the same fixed-covariate columns, in the same order, when calling `predict.prems()`.

## Practical notes

PReMS standardises selected and fixed predictors internally by default. Zero-variance candidate predictors are excluded from the search with a warning. Missing values should be handled before calling PReMS. Because prediction identifies selected predictors by column name, training and prediction matrices should have consistent, unique column names.

The search can become computationally expensive when the number of predictors is large. `max2way` controls the breadth of the initial two-predictor search and `max.s` controls the number of candidate models propagated to larger model sizes. These parameters therefore trade off search breadth against computation.

If a completed search later needs to be extended, `ModelSearchIncrease()` adds one model size at a time without recomputing the smaller models. The same training `y`, `x` and, where applicable, `x.fixed` used for the original fit should be supplied.

## Citation

A manuscript citation should be added here before public release.

The example data are from:

Cios K, Kurgan L and Goodenday L. SPECTF Heart. UCI Machine Learning Repository, 2001. DOI: 10.24432/C5N015.

Efron B, Hastie T, Johnstone I and Tibshirani R. Least angle regression. *Annals of Statistics*. 2004;32:407-499. DOI: 10.1214/009053604000000067.

van de Vijver MJ, He YD, van 't Veer LJ, et al. A gene-expression signature as a predictor of survival in breast cancer. *New England Journal of Medicine*. 2002;347:1999-2009. DOI: 10.1056/NEJMoa021967.
