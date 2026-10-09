# Run after benchmark-runtime.R to compare the preserved baseline with new results.
baseline <- read.csv("scripts/benchmark-runtime-baseline.csv")
updated <- read.csv("scripts/benchmark-runtime-results.csv")
keys <- c("target","countries","variables","observations","SV","draws","burnin","repetitions")
comparison <- merge(baseline,updated,by=keys,suffixes=c("_old","_new"))
stopifnot(nrow(comparison)==nrow(baseline),nrow(comparison)==nrow(updated))
comparison$speedup <- comparison$median_seconds_old/comparison$median_seconds_new
comparison$runtime_reduction_pct <- 100*(1-comparison$median_seconds_new/comparison$median_seconds_old)
comparison$throughput_increase_pct <- 100*(comparison$speedup-1)
write.csv(comparison,"scripts/benchmark-runtime-comparison.csv",row.names=FALSE)
print(comparison[,c(keys,"median_seconds_old","median_seconds_new","speedup","runtime_reduction_pct")],row.names=FALSE)
