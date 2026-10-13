# Prepares the speed_acc example data from rtdists::speed_acc (Wagenmakers, Ratcliff,
# Gomez, & McKoon, 2008, Experiment 1), as distributed with the rtdists package
# (Singmann et al.). Run from the package root to recreate data/speed_acc.rda.

d <- rtdists::speed_acc

# Remove the trials excluded in the original analyses (uninterpretable responses, RTs
# below 180 ms or above 3 s), as flagged by rtdists in `censor`
d <- d[!d$censor & d$response != "error", ]

# Stimulus type: word frequency for words; nonwords form one category
stimulus <- as.character(d$frequency)
stimulus[d$stim_cat == "nonword"] <- "nonword"

speed_acc <- data.frame(
  id = factor(as.integer(as.character(d$id))),
  condition = factor(as.character(d$condition), levels = c("accuracy", "speed")),
  stimulus = factor(stimulus, levels = c("high", "low", "very_low", "nonword")),
  response = as.integer(d$response == "word"),
  rt = d$rt
)
rownames(speed_acc) <- NULL

usethis::use_data(speed_acc, overwrite = TRUE)
