# PReMS

**Documentation:** [Main README](README.md) · [Parallel CV and model-size continuation](PARALLEL_CV.md)

PReMS (Parallel Regularised Model Search) searches for sparse predictive models by combining regularised regression with a structured search over model space. The current implementation is supplied as an R script (`prems.R`).

For each model size, PReMS ranks candidate models by the marginal-likelihood and retains the highest-ranking models for expansion to the next model size. Model size is chosen by cross-validation. When the appropriate search range is not known in advance, PReMS can extend the search sequentially in blocks and use a paired cross-validation stopping rule to determine whether larger models need to be explored.

## Contents

- [Current scope](#current-scope)
- [Installation](#installation)
- [Function reference](#function-reference)
  - [Main functions](#main-functions)
  - [Helper functions](#helper-functions)
  - [Auxiliary functions](#auxiliary-functions)
- [Example: Binary UCI SPECTF Heart data](#example-uci-spectf-heart-data)
- [Model-size selection and sequential stopping](#model-size-selection-and-sequential-stopping)
- [Parallel cross-validation on HPC](#parallel-cross-validation-on-hpc)
- [Example: Continuous tecator meat spectroscopy data](#example-tecator-meat-spectroscopy-data)
- [Example: Survival NKI70 breast-cancer data](#example-nki70-breast-cancer-survival-data)
- [Fixed covariates](#fixed-covariates)
- [Practical notes](#practical-notes)
- [Citation](#citation)

## Current scope

The current script supports Gaussian, binomial and Cox
proportional-hazards outcomes. Models are searched using `prems()` and
model size is selected using `cv.prems()`.

For binary outcomes, `y` is coded `0/1`. For Cox models, `y` should be a `survival::Surv` object; both ordinary right-censored `Surv(time, status)` responses and counting-process `Surv(start, stop, status)` responses are supported. Predictor matrices should be numeric, contain no missing values and have unique column names. Prediction data must use the same predictor names as the training data.

## Installation

PReMS is currently used by sourcing the R script rather than installing an R package.

The core Gaussian and binomial workflow requires `glmnet`. Cox models
additionally require `survival`. For multicore implementation of
`glmnet` used to estimate ridge prior precision parameter `tau`, install `doParallel`.

```r
install.packages("glmnet")

# Required for multicore automatic estimation of tau:
install.packages("doParallel")

# Required only for Cox models:
install.packages("survival")

source("prems.R")
```

Multiple cores can be used on macOS and Linux, on Windows use `no.cores = 1`. Parallel cross-validation across HPC jobs is documented in **[Parallel CV and model-size continuation](PARALLEL_CV.md)**.

The SPECTF example uses `pROC` to calculate AUC:

```r
install.packages("pROC")
```

The Tecator and NKI70 data are downloaded directly in the examples below.

## Function reference

### Main functions

The two principal user-facing functions are `prems()` and `cv.prems()`. PReMS ranks candidate models by the negative log marginal likelihood (`ML`), so smaller values indicate higher-ranked models.

#### `prems()`

```r
prems(
  y, x, family, k.max, x.fixed = NULL, max2way = "all",
  tau = NULL, max.s = 50, no.cores = 1, standardize = TRUE,
  verbose = TRUE, tau.nfolds = 10, tau.attempts = 20
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `y` | — | Outcome. Use a numeric vector for Gaussian or binomial models and a `survival::Surv` object for Cox models. Binary outcomes should be coded `0/1`. |
| `x` | — | Numeric matrix of candidate predictors subject to model selection. Columns should have unique names. |
| `family` | — | Outcome family: `"gaussian"`, `"binomial"` or `"cox"`. |
| `k.max` | — | Largest model size to search. |
| `x.fixed` | `NULL` | Optional matrix of covariates included in every model and not subject to selection. |
| `max2way` | `"all"` | Controls the initial two-predictor search. `"all"` evaluates all two-predictor models; an integer constructs the two-predictor search from that many highest-ranking one-predictor models. |
| `tau` | `NULL` | Ridge-prior precision. When `NULL` estimated automatically. Supply a positive numeric value to override automatic estimation. |
| `max.s` | `50` | Number of highest-ranking models retained at each model size for expansion to the next size. |
| `no.cores` | `1` | Number of cores used for automatic `tau` estimation and the PReMS model search. |
| `standardize` | `TRUE` | Standardise predictors internally before fitting. |
| `verbose` | `TRUE` | Print progress information. |
| `tau.nfolds` | `10` | Number of folds used by `TauEst()` when `tau = NULL`. |
| `tau.attempts` | `20` | Maximum number of attempts to obtain a valid automatic estimate of `tau`. |

Both `family` and `k.max` must be supplied explicitly.

**Returns**

`prems()` returns a PReMS fit object with the following components:

| Component | Meaning |
| --- | --- |
| `null` | Fit of the model containing no selected candidate predictors; fixed covariates are included when supplied. |
| `fitted.models` | Retained fitted models for each model size. Each fitted model contains `ML`, `aic` and `beta`. |
| `model.indicator` | Predictor indices defining each retained model. |
| `cnames` | Names of the candidate predictors supplied in `x`. |
| `cnames.fixed` | Names of the fixed covariates supplied in `x.fixed`, when present. |
| `m`, `sd` | Predictor means and standard deviations used to transform the candidate predictors. |
| `m.fixed`, `sd.fixed` | Means and standard deviations used to transform the fixed covariates. |
| `tau` | Ridge-prior precision used for fitting, including an automatically estimated value when `tau = NULL`. |
| `standardize` | Whether internal predictor standardisation was used. |
| `family` | Fitted outcome family. |

Within each retained model, `ML` is the negative approximate log marginal likelihood, so smaller values indicate better-ranked models. `beta` contains coefficients on the original predictor scale.

#### `cv.prems()`

```r
cv.prems(
  y, x, family, k.max, x.fixed = NULL, no.cores = 1, k.min = 1,
  tau = NULL, max.s = 50, max2way = "all", standardize = TRUE,
  nfolds = 10, foldid = NULL, verbose = TRUE, cv.fit = NULL,
  return.fits = FALSE, tau.nfolds = 10, tau.attempts = 20
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `y` | — | Outcome, in the same form as for `prems()`. |
| `x` | — | Numeric matrix of candidate predictors. |
| `family` | — | Outcome family: `"gaussian"`, `"binomial"` or `"cox"`. |
| `k.max` | — | Largest model size to search and evaluate. When continuing an existing CV fit, this is the new target maximum size. |
| `x.fixed` | `NULL` | Optional covariates included in every model and not subject to selection. |
| `no.cores` | `1` | Number of cores used within each fold for automatic `tau` estimation and the PReMS search. CV folds themselves are fitted sequentially by `cv.prems()`. |
| `k.min` | `1` | Smallest model size for which CV performance is evaluated. |
| `tau` | `NULL` | Optional fixed value of `tau` used in every CV fold. When `NULL`, `tau` is estimated separately within each training fold. |
| `max.s` | `50` | Number of highest-ranking models retained at each model size. When continuing an existing CV search and `max.s` is omitted, the retained search breadth is inherited from the saved fold fits. |
| `max2way` | `"all"` | Controls the two-predictor search as in `prems()`. |
| `standardize` | `TRUE` | Standardise predictors internally within each training fold. |
| `nfolds` | `10` | Number of CV folds. If `foldid` is supplied, the number of folds is inferred from `foldid` when `nfolds` is omitted. |
| `foldid` | `NULL` | Optional vector assigning each observation to a CV fold. Supplying this is useful for reproducibility and for using identical folds across methods. |
| `verbose` | `TRUE` | Print fold and model-search progress. |
| `cv.fit` | `NULL` | Existing `cv.prems()` object to continue to a larger `k.max`. The object must contain retained fold fits. |
| `return.fits` | `FALSE` | Retain the complete fold-specific PReMS fits needed for later continuation. |
| `tau.nfolds` | `10` | Number of folds used for automatic `tau` estimation within each training fold. |
| `tau.attempts` | `20` | Maximum number of attempts to obtain a valid automatic estimate of `tau` within each training fold. |

`family` and `k.max` must be supplied explicitly. By default
`cv.prems()` uses 10-fold CV; supplying `foldid` overrides the
generated fold allocation and, when `nfolds` is omitted, determines
the number of folds. `max.s = 50` is shared with `prems()` so that,
unless overridden, model-size selection and the final full-data fit
use the same search breadth.

**Returns**

`cv.prems()` returns a `cv.prems` object with the following components:

| Component | Meaning |
| --- | --- |
| `best` | Model size with the highest mean cross-validated predictive score. |
| `one.se` | Smaller model-size alternative selected by the one-standard-error rule. |
| `cvm` | Mean predictive score across folds for each evaluated model size. |
| `cvsd` | Standard error of the fold-level predictive scores for each model size. |
| `selected.coef` | Candidate predictors selected in the highest-ranking model at each model size and in each fold. |
| `pwll` | Matrix of fold-level predictive scores, with folds in rows and model sizes in columns. |
| `foldid` | Fold assignment for each observation. |
| `family` | Outcome family used for the CV analysis. |
| `standardize` | Whether internal predictor standardisation was used. |
| `x.names` | Names of the candidate predictors. |
| `x.fixed.names` | Names of the fixed covariates, when present. |
| `has.fixed` | Indicates whether fixed covariates were included. |
| `fold.fits` | Complete fold-specific PReMS fits, included only when `return.fits = TRUE`. These allow the search to be continued without refitting smaller model sizes. |

**RAM usage:** Set `return.fits = TRUE` when you want to continue a CV search later. This retains the complete PReMS search for every CV fold and can therefore substantially increase RAM use. On the final continuation call, set `return.fits = FALSE` if you want the returned object to contain only the CV summaries.

### Helper functions

These functions extract models, make predictions, display CV results and apply the sequential model-size stopping rule.

#### `getModelFit()`

```r
getModelFit(fitted.models, size = 1, rank = 1)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `fitted.models` | — | Object returned by `prems()` or `ModelSearchIncrease()`. |
| `size` | `1` | Model size to extract. |
| `rank` | `1` | Rank of the model to extract within that size according to `ML`; `1` returns the highest-ranking model. |

**Returns**

A fitted model with three components:

| Component | Meaning |
| --- | --- |
| `ML` | Negative approximate log marginal likelihood used to rank models; smaller values indicate a higher-ranked model. |
| `aic` | AIC for the fitted model. |
| `beta` | Named fitted coefficients on the original predictor scale. For Gaussian and binomial models this includes the intercept `I`; fixed covariates are included before the selected candidate predictors. |

#### `predict.prems()`

```r
predict.prems(
  fitted.models, newx, newx.fixed = NULL, size = 1, rank = 1
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `fitted.models` | — | Fitted PReMS object. |
| `newx` | — | Matrix of candidate predictors for the observations to be predicted. Predictor names must match those used for fitting. |
| `newx.fixed` | `NULL` | Fixed covariates for the new observations if fixed covariates were used when fitting. |
| `size` | `1` | Model size used for prediction. |
| `rank` | `1` | `ML` rank of the model used within the selected model size. |

**Returns**

Predictions for the rows of `newx`. For binomial models these are probabilities, for Gaussian models they are predicted outcome values and for Cox models they are linear predictors.

#### `plot.cv.prems()`

```r
plot.cv.prems(cv.fit, ylim = NULL, cex = 1)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `cv.fit` | — | Object returned by `cv.prems()`. |
| `ylim` | `NULL` | Optional y-axis limits. When `NULL`, limits are chosen from the CV means and standard errors. |
| `cex` | `1` | Expansion factor used for axis and label text. |

**Returns**

Produces a plot in the current graphics device. The plot shows mean cross-validated score with standard-error bars against model size and marks the `best` and `one.se` selections.

#### `prems.cv.stop()`

```r
prems.cv.stop(
  cv.fit, block.size = 10, se.multiplier = 1,
  require.complete.block = TRUE
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `cv.fit` | — | `cv.prems()` object containing fold-level scores in `pwll`. |
| `block.size` | `10` | Number of consecutive model sizes in each block. |
| `se.multiplier` | `1` | Number of paired standard errors by which a new block must improve on the currently retained model for the search to continue. |
| `require.complete.block` | `TRUE` | Require a complete final block before the search can be declared adequate. |

**Returns**

| Component | Meaning |
| --- | --- |
| `selected.size` | Model size retained by the sequential stopping rule. |
| `selected.mean` | Mean held-out CV score for `selected.size`. |
| `adequate` | Whether the evaluated model-size range is sufficient for the stopping rule to stop. |
| `reason` | Explanation of the stopping decision. |
| `comparisons` | Block-by-block paired comparisons used by the stopping rule. |
| `mean.score` | Mean held-out CV score for each evaluated model size. |

The fields within `comparisons` and the distinction between this sequential stopping rule and selecting the maximum CV score over a fixed range are described in [Model-size selection and sequential stopping](#model-size-selection-and-sequential-stopping).

### Auxiliary functions

These auxiliary functions support model-search continuation and separate estimation of `tau`.

#### `ModelSearchIncrease()`

```r
ModelSearchIncrease(
  fitted.models, y, x, x.fixed = NULL, k.max = NULL,
  no.cores = 1, max.s = NULL, verbose = TRUE,
  checkpoint = NULL
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `fitted.models` | — | Existing PReMS fit to extend. |
| `y` | — | Original training outcome. |
| `x` | — | Original matrix of candidate predictors. |
| `x.fixed` | `NULL` | Original fixed-covariate matrix, if used. |
| `k.max` | `NULL` | Target maximum model size. When `NULL`, one additional model size is fitted. |
| `no.cores` | `1` | Number of cores used for the added model searches. |
| `max.s` | `NULL` | Number of models retained at each newly fitted size. When `NULL`, the search breadth retained at the current largest size is preserved. |
| `verbose` | `TRUE` | Print progress as each new model size is completed. |
| `checkpoint` | `NULL` | Optional function called after each newly completed model size. See `PARALLEL_CV.md` for checkpointed HPC use. |

The existing `tau`, scaling parameters and fitted smaller model sizes are inherited from `fitted.models`.

**Returns**

An updated PReMS fit object in the same format as `prems()`, extended through the requested `k.max`. Previously fitted smaller model sizes are retained.

#### `TauEst()`

```r
TauEst(
  y, x, family, x.fixed = NULL, standardize = TRUE,
  n.coef = 1, fit = NULL, nfolds = 10, parallel = FALSE,
  no.cores = 1
)
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `y` | — | Outcome. |
| `x` | — | Matrix of candidate predictors. |
| `family` | — | Outcome family passed to `glmnet::cv.glmnet()`: `"gaussian"`, `"binomial"` or `"cox"`. |
| `x.fixed` | `NULL` | Optional fixed covariates, assigned zero Lasso penalty. |
| `standardize` | `TRUE` | Whether predictors are standardised by `glmnet` and when converting the selected coefficient magnitude to `tau`. |
| `n.coef` | `1` | Number of the largest absolute Lasso coefficients used in the calculation of `tau`. |
| `fit` | `NULL` | Optional existing `cv.glmnet` fit. When supplied, the Lasso CV step is not repeated. |
| `nfolds` | `10` | Number of folds used by `cv.glmnet()`. |
| `parallel` | `FALSE` | Use parallel cross-validation when estimating `tau`. |
| `no.cores` | `1` | Number of cores used when `parallel = TRUE`. |

`prems()` invokes `TauEst()` automatically when `tau = NULL`. Call `TauEst()` directly when you want to estimate or inspect `tau` separately.

**Returns**

| Component | Meaning |
| --- | --- |
| `tau.opt` | Estimate of `tau` based on the Lasso model selected at `lambda.min`; this is the value used by `prems()` for automatic estimation. |
| `tau.1se` | Corresponding estimate based on the Lasso `lambda.1se` model; `NA` if that model contains no penalised predictors. |
| `fit.lasso` | Cross-validated `glmnet` fit used to estimate `tau`. |

## Example: Binary UCI SPECTF Heart data

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
  no.cores = 4
)

cv.fit$best
cv.fit$one.se

plot.cv.prems(cv.fit)
```

Searching through size 15 places the observed optimum inside the evaluated range. The SPECTF CV curve nevertheless has a broad plateau around its optimum. With a small data set, the selected size can depend on the fold allocation, so repeated CV replicates are preferable when computation permits. Repeated CV was used for model-size selection in the benchmark analyses in the paper.

For a faster exploratory run, reduce `max.s` or use a restricted two-predictor search such as `max2way = 20`.

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

Repeat the PReMS search using the complete training data and the same maximum model size considered during cross-validation. Because `tau` is omitted, `prems()` estimates it automatically using `TauEst()`.

```r
set.seed(1)

fit <- prems(
  y = y.train,
  x = x.train,
  family = "binomial",
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
  rank = 1
)

selected.model$beta
```

The returned coefficient vector contains the intercept (`I`) followed by the selected predictors. Increase `rank` to inspect alternative high-ranking models of the same size.

### 6. Extend an existing search to a larger model size

`ModelSearchIncrease()` can be used when a fitted search should be extended to larger model sizes without repeating the searches already completed. For example, the fit above searched models of sizes 1--15. The following extends it through size 20:

```r
fit <- ModelSearchIncrease(
  fitted.models = fit,
  y = y.train,
  x = x.train,
  k.max = 20,
  no.cores = 4
)

size20.model <- getModelFit(
  fit,
  size = 20,
  rank = 1
)

size20.model$beta
```

`k.max` is the largest model size to be fitted. If it is omitted, `ModelSearchIncrease()` adds one model size. By default the search breadth is inherited from the current largest model size; supply `max.s` to override this.

### 7. Compare PReMS and Lasso in the independent test set

Generate held-out predictions for PReMS, `Lasso.min` and `Lasso.1se` and compare both AUC and model size.

```r
prems.pred <- predict.prems(
  fit,
  newx = x.test,
  size = selected.size,
  rank = 1
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

When computation is modest, a broad range of model sizes can be searched in one `cv.prems()` call, as in the SPECTF example above. For more expensive analyses, model sizes can instead be explored sequentially in blocks.

`prems.cv.stop()` applies a different model-size rule from simply taking the model size with the highest mean CV score over a fixed range. It starts with the best model in the first block, then compares this retained model with the best model in each subsequent block using paired differences in held-out CV score across folds. The search is extended only when the new block improves predictive performance by more than `se.multiplier` paired standard errors; the default is one paired standard error.

The object returned by `prems.cv.stop()` contains:

| Component | Meaning |
| --- | --- |
| `selected.size` | Model size retained by the sequential stopping rule. This is not necessarily the largest model size evaluated. |
| `selected.mean` | Mean held-out CV score for `selected.size`. |
| `adequate` | `TRUE` when the evaluated range is sufficient for the stopping rule to stop; `FALSE` when another complete block should be explored. |
| `reason` | Text explaining why the search stopped or why another block is required. |
| `comparisons` | Data frame containing the block-by-block paired comparisons. |
| `mean.score` | Mean held-out CV score for every model size evaluated. |

Within `comparisons`, `previous_size` is the currently retained model size and `candidate_size` is the best model in the new block. `mean_difference` is the mean fold-wise difference in predictive score between them, `se_difference` is its standard error and `threshold` is `se.multiplier * se_difference`. `improves` is `TRUE` only when `mean_difference > threshold`. `block_complete` indicates whether the full requested block has been evaluated.

### SPECTF example in blocks of five

The initial search evaluates model sizes 1--5:

```r
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
  return.fits = TRUE
)

stop.fit <- prems.cv.stop(cv.fit, block.size = 5)

stop.fit$selected.size
stop.fit$adequate
stop.fit$reason
```

With only the first block evaluated, `adequate` is `FALSE`: the stopping rule requires another block before deciding whether the explored range is sufficient.

`return.fits = TRUE` retains the fold-specific PReMS searches needed for continuation. If another block is required, pass the returned object back to `cv.prems()`. The following extends each fold from size 5 to size 10; sizes 1--5 are not refitted:

**RAM usage:** retaining fold fits can make a `cv.prems` object large because it stores the complete PReMS search for every CV fold. The default `return.fits = FALSE` keeps the returned object compact. Set `return.fits = TRUE` when continuation is required, particularly for analyses with large `max.s` or searches extending to many model sizes.

```r
cv.fit <- cv.prems(
  y = y.train,
  x = x.train,
  k.max = 10,
  cv.fit = cv.fit,
  no.cores = 4
)

stop.fit <- prems.cv.stop(cv.fit, block.size = 5)

stop.fit$selected.size
stop.fit$adequate
stop.fit$reason
stop.fit$comparisons
```

For the SPECTF data, the search is adequate after model sizes 1--10 have been evaluated, so `stop.fit$adequate` is `TRUE` and the search ends at this point rather than extending through size 15. This means that the **search stops after evaluating model size 10**; it does not imply that `selected.size` must equal 10.

This sequential result need not agree with the model size selected in the first SPECTF example. There, `cv.fit$best` is simply the model size with the highest mean CV predictive score among all sizes 1--15. In contrast, `prems.cv.stop()` uses a blockwise stopping criterion based on paired fold-level improvements and may stop before the full range is evaluated. The two procedures therefore answer slightly different questions and can select different model sizes, particularly when the CV curve is fairly flat around its maximum.

When `cv.fit` is supplied, the original folds, outcome family and fold-specific shrinkage estimates are retained. On a continuation call, omitting `return.fits` retains the fold fits so that the search can be extended again. On the final continuation call, set `return.fits = FALSE` if you want the returned object to contain only the CV summaries. The same `y`, `x` and, where applicable, `x.fixed` used in the original CV analysis should be supplied.

The SPECTF predictive log-likelihood is fairly flat around its optimum. For small or moderate datasets, repeated CV replicates can therefore give a more stable estimate of model size; this was the approach used for the benchmark datasets in the paper. Sequential stopping is most useful when repeatedly searching a wide range of model sizes would be computationally expensive.

## Parallel cross-validation on HPC

For large analyses, cross-validation folds can be fitted independently and model searches checkpointed and continued on an HPC system. See **[Parallel CV and model-size continuation](PARALLEL_CV.md)** for the batch/HPC workflow.

## Example: Continuous tecator meat spectroscopy data

The Tecator data provide a compact continuous-outcome example with a predictor structure closer to high-dimensional molecular data. The data contain 215 meat samples, each measured for absorbance at 100 wavelengths between 850 and 1050 nm. Fat, water and protein percentages were measured by analytic chemistry. Here the 100 raw absorbance measurements are candidate predictors and fat percentage is the continuous outcome.

The predictors are measurements of the same underlying spectrum and are therefore strongly correlated, making the data useful for illustrating sparse selection among correlated predictors. The first 129 samples were originally used as a training set; this example retains that historical split and uses the remaining 86 samples for independent evaluation.

### 1. Download, load and prepare the data

The `tecator` object can be downloaded directly from the `fda.usc` package's GitHub repository, so the example can be run without installing the package.

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

The absorbance measurements are used as predictors. Water and protein percentage are laboratory-measured composition outcomes, like fat percentage, so they are excluded from the predictor matrix.

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
  no.cores = 4
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

With `tau = NULL` (the default), `prems()` estimates the shrinkage precision automatically from the training data.

```r
set.seed(1)

fit <- prems(
  y = y.train,
  x = x.train,
  family = "gaussian",
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
  rank = 1
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
    rank = 1
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

## Example: Survival NKI70 breast-cancer data

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
  no.cores = 4
)
```

### 3. Refit PReMS to the complete data

After selecting model size, repeat the search using the complete data. `prems()` estimates `tau` automatically because no explicit value is supplied.

```r
fit <- prems(
  y = y,
  x = x,
  x.fixed = x.fixed,
  family = "cox",
  k.max = 5,
  max2way = "all",
  max.s = 50,
  no.cores = 4
)

selected.size <- cv.fit$best

selected.model <- getModelFit(
  fit,
  size = selected.size,
  rank = 1
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

Covariates supplied through `x.fixed` are included in every model and excluded from predictor selection. When no fixed covariates are required, use the default `x.fixed = NULL`.

For example, if age and sex are to be included in all models:

```r
fit <- prems(
  y = y,
  x = biomarkers,
  x.fixed = cbind(age = age, sex = sex),
  family = "binomial",
  k.max = 5,
  max.s = 50,
  no.cores = 4
)
```

Use the same fixed-covariate columns, in the same order, when calling `predict.prems()`.

## Practical notes

PReMS standardises selected and fixed predictors internally by default. Zero-variance candidate predictors are excluded from the search with a warning. Missing values should be handled before calling PReMS. Because prediction identifies selected predictors by column name, training and prediction matrices should have consistent, unique column names.

If `tau = NULL` (the default), `prems()` estimates the ridge-prior precision automatically using `TauEst()`. The estimated value is stored as `fit$tau`. Supply a positive numeric `tau` to override automatic estimation. `tau.nfolds` and `tau.attempts` control this estimation step.

The search can become computationally expensive when the number of predictors is large. `max2way` controls the breadth of the initial two-predictor search and `max.s` controls the number of candidate models propagated to larger model sizes. These parameters therefore trade off search breadth against computation.

By default, `cv.prems()` returns the CV summaries without retaining the fold-specific PReMS fits. Set `return.fits = TRUE` when model-size continuation is required; this stores a complete search object for every CV fold and can substantially increase RAM use.

If a completed search later needs to be extended, `ModelSearchIncrease()` can extend it directly to a specified `k.max` without recomputing the smaller models. The same training `y`, `x` and, where applicable, `x.fixed` used for the original fit should be supplied.

## Citation

A manuscript citation should be added here before public release.

The example data are from:

Cios K, Kurgan L and Goodenday L. SPECTF Heart. UCI Machine Learning Repository, 2001. DOI: 10.24432/C5N015.

Efron B, Hastie T, Johnstone I and Tibshirani R. Least angle regression. *Annals of Statistics*. 2004;32:407-499. DOI: 10.1214/009053604000000067.

van de Vijver MJ, He YD, van 't Veer LJ, et al. A gene-expression signature as a predictor of survival in breast cancer. *New England Journal of Medicine*. 2002;347:1999-2009. DOI: 10.1056/NEJMoa021967.
