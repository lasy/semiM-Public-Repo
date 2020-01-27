par = list()
par$data_type = "full" # "synthetic" "subset" "full"
par$local_user = Sys.getenv("LOGNAME")
par$n_cores = detectCores() - 1
par$reset = FALSE
par$reset_steps = c("","","")

source("Scripts/00_variables_IO.R")





mucus.dict = data.frame(
  names = c("none", 
            "sticky_little",
            "sticky_medium",
            "sticky_lots",
            "creamy_little",
            "creamy_medium",
            "creamy_lots",
            "eggwhite_little",
            "eggwhite_medium",
            "eggwhite_lots",
            "watery_little",
            "watery_medium",
            "watery_lots"
  ),
  score = c(0, # none
            rep(0,3), # sticky
            0.1,
            0.2,
            0.3,
            0.5, #  "eggwhite_little",
            0.9, # "eggwhite_medium",
            1,   # "eggwhite_lots",
            0.5, # "watery_little",
            0.9,
            1
  )
)
