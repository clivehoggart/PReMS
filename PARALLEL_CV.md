# Parallel CV and model-size continuation

**Documentation:** [Main README](README.md) · **Parallel CV and model-size continuation**

`prems_parallel.R` provides a command-line wrapper for running PReMS cross-validation folds as independent jobs and for extending saved searches without repeating model sizes that have already been fitted. This is intended for analyses where `cv.prems()` would be inefficient because the outer folds are run sequentially.

Place `prems_parallel.R` in the same directory as `prems.R`. The wrapper locates and sources `prems.R` automatically, so the path to `prems.R` does not need to be supplied on the command line. The wrapper can therefore be called from any working directory.

The wrapper has two commands:

- `fit` fits either the complete data or one CV training fold. If the output RDS already exists, it loads the saved fit and extends it only through model sizes that have not yet been searched.
- `summarise` loads the independently fitted CV folds, evaluates each saved model size on the corresponding held-out fold and returns a `cv.prems()`-style summary containing `best`, `one.se`, `cvm`, `cvsd` and `selected.coef`.

There are two distinct levels of parallelism. Scheduler array jobs can run the outer CV folds in parallel, while `--cores` controls candidate-model evaluation within each fold.

## Example: SPECTF Heart

This example uses the UCI SPECTF Heart training data for 10-fold cross-validation, followed by evaluation of the selected model in the supplied independent test set.

### 1. Prepare the analysis RDS

Download and extract `SPECTF.train` and `SPECTF.test` from the UCI SPECTF Heart dataset and place them in the working directory.

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

# Generate the folds once. Every parallel job uses these assignments.
foldid <- make.folds2(y.train, folds = 10, seed = 1)

parallel.data <- list(
  y = y.train,
  x = x.train,
  x.fixed = NULL,
  family = "binomial",
  foldid = foldid
)

saveRDS(parallel.data, "spectf_parallel.rds")
```

The analysis RDS must contain `y`, `x` and `family`. For parallel CV it must also contain `foldid`; `x.fixed` is optional. For Cox analyses, `y` can be either `Surv(time, status)` or `Surv(start, stop, status)`.

### 2. Run one fold

A single fold can be fitted directly:

```sh
mkdir -p cv

Rscript prems_parallel.R fit \
  --data spectf_parallel.rds \
  --fold 1 \
  --out cv/spectf_fold1.rds \
  --initial-k 5 \
  --target-k 5 \
  --max-s 50 \
  --max2way all \
  --cores 4 \
  --seed 1
```

The specified fold is left out. `tau` is estimated from the remaining observations, so its estimation is contained within the outer training fold. A fixed value can instead be supplied with `--tau`.

For a new fit, `--initial-k` is the largest model size searched by the initial call to `prems()`. `--target-k` is the final model size requested. If the target is larger than the initial search, the additional sizes are fitted sequentially using `ModelSearchIncrease()`.

### 3. Run all ten folds in parallel

On an LSF cluster, the ten folds can be submitted as an array job:

```sh
#!/bin/sh

#BSUB -L /bin/sh
#BSUB -J spectf_prems[1-10]
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 02:00
#BSUB -M 8000
#BSUB -o logs/spectf.o%J.%I
#BSUB -eo logs/spectf.e%J.%I

module load R/4.5.1

fold=${LSB_JOBINDEX}

Rscript prems_parallel.R fit \
  --data spectf_parallel.rds \
  --fold ${fold} \
  --out cv/spectf_fold${fold}.rds \
  --initial-k 5 \
  --target-k 5 \
  --max-s 50 \
  --max2way all \
  --cores 4 \
  --seed 1
```

Submit it with, for example:

```sh
mkdir -p logs cv
bsub < spectf_cv.sh
```

Each array task runs one outer fold. Within a task, `--cores 4` parallelises candidate-model evaluation. The number supplied to `--cores` should not exceed the number of cores allocated to that scheduler task.

### 4. Combine the completed folds

After all ten fold jobs have completed:

```sh
Rscript prems_parallel.R summarise \
  --data spectf_parallel.rds \
  --fit-template 'cv/spectf_fold{fold}.rds' \
  --out spectf_cv.rds \
  --criteria ML
```

The summarisation step does not refit the models. Each saved fold-specific fit is applied to that fold's held-out observations. For binomial outcomes the score is predictive log-likelihood, for Gaussian outcomes it is negative squared prediction error and for Cox outcomes it is held-out Cox partial log-likelihood.

Inspect the result in R:

```r
source("prems.R")

cv.fit <- readRDS("spectf_cv.rds")

cv.fit$best
cv.fit$one.se
cv.fit$cvm
cv.fit$cvsd

plot.cv.prems(cv.fit)
```

The selected predictors in each fold can also be inspected. For example:

```r
k <- cv.fit$one.se

cv.fit$selected.coef[[k]]

sort(
  table(as.vector(cv.fit$selected.coef[[k]])),
  decreasing = TRUE
)
```

These selection frequencies are useful diagnostically but are not themselves the PReMS optimisation criterion.

### 5. Apply the blockwise CV stopping rule

For larger searches, model size should not be chosen simply because the current CV optimum is below the largest size evaluated. The supplied `prems_cv_stopping_rule.R` implements the paired blockwise stopping rule used for the large PReMS analyses. It compares the best model in each new block with the previously retained model using the fold-level difference in predictive performance. The search is stopped only when the best model in a complete new block fails to improve mean held-out performance by more than the specified multiple of its paired standard error.

`prems_parallel.R summarise` stores the fold-by-model-size scores in `cv.fit$pwll`, so the stopping rule can be applied directly to its output. For this small SPECTF example, use blocks of five model sizes:

```sh
Rscript prems_cv_stopping_rule.R \
  --cv spectf_cv.rds \
  --block-size 5 \
  --out-prefix spectf_cv
```

The command writes:

- `spectf_cv_stopping_summary.csv`, containing the currently selected model size, whether the explored search range is adequate and the reason for the decision;
- `spectf_cv_stopping_details.rds`, containing the fold-paired comparisons between successive blocks.

After an initial search over sizes 1--5, the stopping rule will report that the search is not yet adequate because only the initial block has been evaluated. At least one further complete block is required. After extending the folds through size 10 and rerunning `summarise`, run the stopping rule again. It then compares the best model among sizes 6--10 with the retained best model from sizes 1--5 using the within-fold score differences.

For example, inspect the result directly in R with:

```r
stop.result <- readRDS("spectf_cv_stopping_details.rds")
stop.result$selected.size
stop.result$adequate
stop.result$reason
stop.result$comparisons
```

If `adequate` is `FALSE` because the final block improved by more than the required paired SE, extend the fold fits by another block and repeat the `summarise` and stopping-rule steps. If `adequate` is `TRUE`, the search can stop and `selected.size` is the model size selected by the sequential rule.

The SPECTF example uses a block size of five only to keep the demonstration small. The UK Biobank analyses described in the manuscript used blocks of 10, corresponding to:

```sh
Rscript prems_cv_stopping_rule.R \
  --cv my_cv.rds \
  --block-size 10 \
  --out-prefix my_cv
```

The default is one paired standard error (`--se-multiplier 1`). By default the final block must be complete before lack of improvement is accepted as evidence to stop.

### 6. Extend the model-size search

If the stopping rule indicates that additional model sizes are required, extend the existing fold fits rather than starting again. For example, if the original search stopped at size 5:

```sh
#!/bin/sh

#BSUB -L /bin/sh
#BSUB -J spectf_prems_extend[1-10]
#BSUB -n 4
#BSUB -R "span[hosts=1]"
#BSUB -W 02:00
#BSUB -M 8000
#BSUB -o logs/spectf_extend.o%J.%I
#BSUB -eo logs/spectf_extend.e%J.%I

module load R/4.5.1

fold=${LSB_JOBINDEX}

Rscript prems_parallel.R fit \
  --data spectf_parallel.rds \
  --fold ${fold} \
  --out cv/spectf_fold${fold}.rds \
  --target-k 10 \
  --max-s 50 \
  --cores 4
```

Because each output file already contains the models through size 5, the wrapper loads the existing fit and calls `ModelSearchIncrease()` for sizes 6 through 10. The RDS is saved after each newly completed model size, so an interrupted job can be resubmitted without repeating completed sizes.

After all folds have been extended, rerun the same `summarise` command. The summary automatically uses the largest model size available in all fold fits unless `--k-max` is supplied explicitly. Then rerun `prems_cv_stopping_rule.R`. Continue adding complete blocks while the stopping-rule output has `adequate = FALSE` because the final available block provides an improvement greater than the paired-SE threshold.

### 7. Fit the complete training data

Once the stopping rule reports that the explored range is adequate, use its `selected.size` for the final PReMS fit. For example:

```r
stop.result <- readRDS("spectf_cv_stopping_details.rds")
k <- stop.result$selected.size
k
```

Then omit `--fold` to fit all training observations through at least that model size. If, for example, `k = 7`:

```sh
Rscript prems_parallel.R fit \
  --data spectf_parallel.rds \
  --out spectf_fit.rds \
  --initial-k 5 \
  --target-k 7 \
  --max-s 50 \
  --max2way all \
  --cores 8 \
  --seed 1
```

If a full-data fit already exists through a smaller size, the same command can be rerun with the larger `--target-k`; only the additional sizes are fitted. Continuation uses the `tau`, scaling information, predictor names and family stored in the existing fit. If the candidate predictors or fixed covariates change, start a new fit rather than continuing an existing RDS.

If the stopping rule is not being used, `cv.fit$best` or `cv.fit$one.se` can instead be used according to the pre-specified model-size selection rule.

### 8. Examine the final model and independent test performance

Using the model size selected by the stopping rule:

```r
library(pROC)
source("prems.R")

stop.result <- readRDS("spectf_cv_stopping_details.rds")
fit <- readRDS("spectf_fit.rds")

k <- stop.result$selected.size

model <- getModelFit(
  fit,
  size = k,
  criteria = "ML"
)

model$beta
```

Generate predictions in the independent SPECTF test set:

```r
pred <- predict.prems(
  fit,
  newx = x.test,
  size = k,
  criteria = "ML"
)

roc.test <- pROC::roc(
  response = y.test,
  predictor = as.numeric(pred),
  quiet = TRUE
)

pROC::auc(roc.test)
```

The selected predictors can be extracted with:

```r
names(model$beta)[-1]
```

## Command reference

Run:

```sh
Rscript prems_parallel.R --help
```

for the complete option list. The principal `fit` options are:

| Option | Meaning |
| --- | --- |
| `--data` | Analysis RDS containing `y`, `x`, `family`, optional `x.fixed` and, for CV, `foldid`. |
| `--fold` | Fold to leave out. Omit for a full-data fit. |
| `--out` | Saved PReMS fit. Existing files are continued rather than overwritten from scratch. |
| `--initial-k` | Largest size searched in the initial `prems()` call for a new fit. |
| `--target-k` | Largest model size required. |
| `--max-s` | Number of highest-ranking models retained at each size. |
| `--max2way` | `all` or the number of top one-predictor models used to construct the initial two-predictor search. |
| `--cores` | Number of cores used within a PReMS search. |
| `--tau` | Optional fixed prior precision. If omitted, `tau` is estimated from the training data. |
| `--tau-folds` | Number of `glmnet` folds used when estimating `tau`. |
| `--seed` | Random seed used for a new fit and `tau` estimation. |

The principal `summarise` options are:

| Option | Meaning |
| --- | --- |
| `--fit-template` | Fold-fit path containing `{fold}`, for example `cv/fold{fold}.rds`. |
| `--out` | Output CV-summary RDS. |
| `--k-min` | Smallest model size to evaluate. |
| `--k-max` | Largest model size to evaluate. If omitted, use the largest size common to all folds. |
| `--criteria` | `ML` or `aic` for ranking models within each size. |

## Workflow summary

The intended workflow is:

1. Generate the fold assignments once and save them with the analysis data.
2. Fit the outer folds independently, preferably as scheduler array jobs.
3. Summarise held-out predictive performance across folds.
4. Apply `prems_cv_stopping_rule.R` to the fold-level CV scores.
5. If the search is not adequate, extend the saved fold fits by another block using a larger `--target-k`, rerun `summarise` and reapply the stopping rule.
6. Fit or extend the full training-data model to the model size selected by the stopping rule.
7. Evaluate the chosen model once on any independent test data.
