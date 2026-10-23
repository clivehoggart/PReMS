# PReMS

**Documentation:** [Main README](README.md) · [Parallel CV and model-size continuation](PARALLEL_CV.md)

PReMS (Parallel Regularised Model Search) searches for sparse predictive models by combining regularised regression with a structured search over model space. The current implementation is supplied as an R script (`prems.R`).

For each model size, PReMS ranks candidate models using an approximate marginal-likelihood score and retains the highest-ranking models for expansion to the next model size. Model size is chosen by cross-validation. When the appropriate search range is not known in advance, PReMS can extend the search sequentially in blocks and use a paired cross-validation stopping rule to determine whether larger models need to be explored.

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

`parallel` is a recommended R package distributed with R and does not normally need to be installed separately. PReMS calls `parallel::mclapply()` explicitly. On Windows use `no.cores = 1`; multiple cores can be used with `mclapply()` on macOS and Linux. The optional `prems_parallel.R` command-line wrapper uses only base R in addition to the packages required by `prems.R`.

The SPECTF example uses `pROC` to calculate AUC:

```r
install.packages("pROC")
```

The Tecator and NKI70 data are downloaded directly in the examples below, so no additional R packages are needed to obtain these data. `pROC` is otherwise required only if the convenience functions `cv.auc()` or `my.auc()` are used. `gplots`, `igraph`, `MASS` and `exvatools` are not required by the current script.

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
  k.max = 15,
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

Searching through size 15 places the observed optimum inside the evaluated range. The SPECTF CV curve nevertheless has a broad plateau around its optimum. With a small data set, the selected size can depend on the fold allocation, so repeated CV replicates are preferable when computation permits. Repeated CV was used for model-size selection in the benchmark analyses in the paper.

For a faster exploratory run, reduce `max.s` or use a restricted two-predictor search such as `max2way = 50`.

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
  k.max = 15,
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

The returned coefficient vector contains the intercept (`I`) followed by the selected predictors. Increase `rank` to inspect alternative high-ranking models of the same size.

### 6. Extend an existing search to a larger model size

`ModelSearchIncrease()` can be used when a fitted search should be extended by one additional model size without repeating the searches already completed. For example, the fit above searched models of sizes 1--15. The following extends it to size 16:

```r
fit <- ModelSearchIncrease(
  fitted.models = fit,
  y = y.train,
  x = x.train,
  no.cores = 4
)

size16.model <- getModelFit(
  fit,
  size = 16,
  rank = 1,
  criteria = "ML"
)

size16.model$beta
```

Each call to `ModelSearchIncrease()` adds one model size. By default it retains the same search breadth as the current largest size; supply `max.s` to override this.

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

The predefined SPECTF test set is highly imbalanced, so AUC is more informative than raw accuracy here. This is a workflow example rather than a definitive benchmark.


## Model-size selection and sequential stopping

The largest useful model size is usually not known in advance. When computation is modest, the simplest approach is to search a broad range in one `cv.prems()` call, as in the SPECTF example above. For larger analyses, model sizes can instead be explored sequentially in blocks.

`prems_cv_stopping_rule.R` compares the best model in each new block with the current best model using paired differences in held-out CV score across folds. The search continues only if the new block improves mean predictive performance by more than one paired standard error. The defaults are blocks of 10 and a one-SE threshold.

### SPECTF example in blocks of five

The same SPECTF analysis can be rerun sequentially to illustrate the stopping rule. Use the same `foldid` throughout.

```r
run.spectf.cv <- function(k.max) {
  cv.prems(
    y = y.train,
    x = x.train,
    family = "binomial",
    k.min = 1,
    k.max = k.max,
    foldid = foldid,
    max2way = "all",
    max.s = 50,
    no.cores = 4,
    criteria = "ML"
  )
}

cv.fit <- run.spectf.cv(5)
saveRDS(cv.fit, "spectf_cv.rds")
```

Apply the stopping rule from the command line:

```sh
Rscript prems_cv_stopping_rule.R \
  --cv spectf_cv.rds \
  --block-size 5 \
  --out-prefix spectf_cv
```

With only sizes 1--5 evaluated, another complete block is required. Rerun through size 10 and apply the same command again:

```r
cv.fit <- run.spectf.cv(10)
saveRDS(cv.fit, "spectf_cv.rds")
```

If the 6--10 block improves by more than one paired SE, continue through size 15:

```r
cv.fit <- run.spectf.cv(15)
saveRDS(cv.fit, "spectf_cv.rds")
```

After each extension, rerun `prems_cv_stopping_rule.R`. The output reports the retained model size and whether the explored range is adequate.

For large analyses, rerunning `cv.prems()` would refit all smaller model sizes. The parallel/checkpointed workflow below avoids this by extending saved fold-specific searches with `ModelSearchIncrease()`.

## Parallel cross-validation and model-size continuation

For larger analyses, the supplied `prems_parallel.R` wrapper can distribute outer cross-validation folds across scheduler jobs and checkpoint searches so that completed model sizes do not need to be refitted. It also supports extending an existing full-data or fold-specific fit with `ModelSearchIncrease()`.

See **[Parallel CV and model-size continuation](PARALLEL_CV.md)** for a complete SPECTF example, including data preparation, an LSF array job, CV aggregation, application of the stopping rule, blockwise extension of the search and final model fitting.

## Example: Tecator meat spectroscopy data

The Tecator data provide a compact continuous-outcome example with a predictor structure closer to high-dimensional molecular data. The data contain 215 meat samples, each measured for absorbance at 100 wavelengths between 850 and 1050 nm. Fat, water and protein percentages were measured by analytic chemistry. Here the 100 raw absorbance measurements are candidate predictors and fat percentage is the continuous outcome.

The predictors are measurements of the same underlying spectrum and are therefore strongly correlated, making the data useful for illustrating sparse selection among correlated predictors. The first 129 samples were originally used as a training set; this example retains that historical split and uses the remaining 86 samples for independent evaluation.

### 1. Download, load and prepare the data

The `tecator` object distributed with the `fda.usc` package can be downloaded directly from the package's GitHub repository; `fda.usc` itself does not need to be installed.

```r
tecator.url <- paste0(
  "https://github.com/moviedo5/fda.usc/raw/refs/heads/master/",
  "data/tecator.rda"
)

download.file(tecator.url, "tecator.rda", mode = "wb")
load("tecator.rda")

x <- as.matrix(tecator$absorp.fdata$data)
y <- tecator$y$Fat

# Give the spectral predictors informative, unique names.
colnames(x) <- paste0("wl_", tecator$absorp.fdata$argvals)

train.id <- 1:129
test.id <- 130:nrow(x)

x.train <- x[train.id, , drop = FALSE]
y.train <- y[train.id]

x.test <- x[test.id, , drop = FALSE]
y.test <- y[test.id]
```

Only the absorbance measurements are used as predictors. Water and protein percentage are not included because, like fat percentage, they are laboratory-measured composition outcomes rather than spectral predictors.

### 2. Select PReMS model size by cross-validation

Use the same fold assignments for PReMS and Lasso. With 100 candidate predictors, an exhaustive search of all two-predictor models is possible but makes a README example unnecessarily slow. Here `max2way = 20` expands the 20 highest-ranking one-predictor models when constructing the initial two-predictor search.

```r
set.seed(1)
foldid <- make.folds.continuous(length(y.train), folds = 10)

cv.fit <- cv.prems(
  y = y.train,
  x = x.train,
  family = "gaussian",
  k.min = 1,
  k.max = 10,
  foldid = foldid,
  max2way = 20,
  max.s = 50,
  no.cores = 4,
  criteria = "ML"
)

cv.fit$best
cv.fit$one.se

plot.cv.prems(cv.fit)
```

Set `max2way = "all"` if an exhaustive search of all two-predictor models is preferred.

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
  max2way = 20,
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

The selected predictor names correspond to wavelengths, which makes it straightforward to inspect whether alternative high-ranking models select nearby, highly correlated parts of the spectrum.

### 5. Compare prediction in the independent test set

For a continuous outcome, predictive performance is summarised here using root mean squared error (RMSE) and squared Pearson correlation (`R2`).

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

Because the historical test set is not used for either model search or tuning, it provides a direct held-out comparison of PReMS and Lasso. The purpose of this example is to demonstrate continuous-outcome analysis and sparse selection in the presence of correlated predictors, rather than to treat individual wavelengths as independent biological signals.

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
