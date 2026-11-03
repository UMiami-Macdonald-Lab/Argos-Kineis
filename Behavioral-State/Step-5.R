library(aniMotum)

# Uses aligned_run_dir from the previous chunk.
stopifnot(
  exists("aligned_run_dir"),
  dir.exists(aligned_run_dir)
)

# Inspect all eight fits, including the two with warnings.
for (shark_id in c("34885", "34890", "45270", "45273")) {
  
  for (system_name in c("Argos", "Kineis")) {
    
    fit_path <- file.path(
      aligned_run_dir,
      paste0(
        system_name,
        "_",
        shark_id,
        "_fit.rds"
      )
    )
    
    cat(
      "\n\n",
      system_name,
      "— shark",
      shark_id,
      "\n"
    )
    
    if (!file.exists(fit_path)) {
      cat("Saved fit not found:", fit_path, "\n")
      next
    }
    
    saved_fit <- readRDS(fit_path)
    
    print(
      summary(saved_fit)
    )
  }
}