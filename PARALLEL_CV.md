# Parallel CV and model-size continuation

**Documentation:** [Main README](README.md) · **Parallel CV and model-size continuation**

`prems_parallel.R` provides R functions for fitting CV folds independently, checkpointing model searches and extending saved fits without repeating completed model sizes. The same functions are available through command-line wrappers for batch or HPC use.

The primary R interface is:

```r
source("prems.R")
source("prems_parallel.R")
```

- `prems_parallel_fit()` fits one CV training fold or the complete data and saves a checkpointed PReMS fit.
- `prems_parallel_summarise()` combines fold-specific fits into a `cv.prems()`-style object.
- `prems.cv.stop()` applies the paired blockwise stopping rule to the resulting fold-level CV scores.

## Example: SPECTF Heart

Assume `y.train`, `x.train` and the 10-fold `foldid` have been created as in the main README.

### 1. Prepare the analysis object

```r
parallel.data <- list(
  y = y.train,
  x = x.train,
  x.fixed = NULL,
  family = "binomial",
  foldid = foldid
)

dir.create("cv", showWarnings = FALSE)
```

### 2. Fit the CV folds

A single fold can be fitted with:

```r
prems_parallel_fit(
  data = parallel.data,
  fold = 1,
  out.file = "cv/spectf_fold1.rds",
  initial.k = 5,
  target.k = 5,
  max.s = 50,
  max2way = "all",
  no.cores = 4
)
```

`tau` is estimated using the training observations for that fold. Each completed model size is saved to `out.file`.

On macOS or Linux, the outer folds can also be run in parallel from R. In this case it is usually simplest to allocate one core within each PReMS fit and parallelise across folds:

```r
folds <- sort(unique(foldid))

parallel::mclapply(
  folds,
  function(f) {
    prems_parallel_fit(
      data = parallel.data,
      fold = f,
      out.file = sprintf("cv/spectf_fold%d.rds", f),
      initial.k = 5,
      target.k = 5,
      max.s = 50,
      max2way = "all",
      no.cores = 1
    )
  },
  mc.cores = 4
)
```

Alternatively, run folds sequentially and use multiple cores within each fold. Avoid allocating the full core count at both levels.

### 3. Summarise CV performance

```r
cv.fit <- prems_parallel_summarise(
  data = parallel.data,
  fit.template = "cv/spectf_fold{fold}.rds"
)

cv.fit$best
cv.fit$one.se
plot.cv.prems(cv.fit)
```

The returned object contains `pwll`, the fold-by-model-size predictive score matrix used by the stopping rule.

### 4. Apply the stopping rule

```r
stop.fit <- prems.cv.stop(
  cv.fit,
  block.size = 5
)

stop.fit$selected.size
stop.fit$adequate
stop.fit$comparisons
```

With only sizes 1--5 fitted, another complete block is required.

### 5. Extend the existing fold fits

The same fold files can be extended to size 10 without repeating sizes 1--5:

```r
parallel::mclapply(
  folds,
  function(f) {
    prems_parallel_fit(
      data = parallel.data,
      fold = f,
      out.file = sprintf("cv/spectf_fold%d.rds", f),
      target.k = 10,
      max.s = 50,
      no.cores = 1
    )
  },
  mc.cores = 4
)
```

Then resummarise and reapply the stopping rule:

```r
cv.fit <- prems_parallel_summarise(
  data = parallel.data,
  fit.template = "cv/spectf_fold{fold}.rds"
)

stop.fit <- prems.cv.stop(cv.fit, block.size = 5)
```

Continue in complete blocks while the final block improves predictive performance by more than the paired-SE threshold.

### 6. Fit the complete training data

Once the stopping rule is satisfied, fit the complete training set through the selected size:

```r
k <- stop.fit$selected.size

fit <- prems_parallel_fit(
  data = parallel.data,
  out.file = "spectf_fit.rds",
  initial.k = min(5, k),
  target.k = k,
  max.s = 50,
  max2way = "all",
  no.cores = 4
)
```

If `spectf_fit.rds` already exists at a smaller size, only the additional sizes are fitted.

## Optional Bash/HPC wrappers

For cluster execution, save the analysis object once:

```r
saveRDS(parallel.data, "spectf_parallel.rds")
```

One fold can then be fitted from Bash with:

```sh
Rscript prems_parallel.R fit \
  --data spectf_parallel.rds \
  --fold 1 \
  --out cv/spectf_fold1.rds \
  --initial-k 5 \
  --target-k 5 \
  --max-s 50 \
  --max2way all \
  --cores 4
```

An LSF array job can substitute `${LSB_JOBINDEX}` for the fold number:

```sh
#!/bin/sh
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
  --cores 4
```

After the array completes, either return to R and use `prems_parallel_summarise()`, or use the wrapper:

```sh
Rscript prems_parallel.R summarise \
  --data spectf_parallel.rds \
  --fit-template 'cv/spectf_fold{fold}.rds' \
  --out spectf_cv.rds
```

The stopping rule likewise has an optional wrapper:

```sh
Rscript prems_cv_stopping_rule.R \
  --cv spectf_cv.rds \
  --block-size 5 \
  --out-prefix spectf_cv
```

These scripts call the same R functionality described above; they are intended for scheduling and batch execution rather than as a separate analysis interface.

## Command-line help

```sh
Rscript prems_parallel.R --help
Rscript prems_cv_stopping_rule.R --help
```
