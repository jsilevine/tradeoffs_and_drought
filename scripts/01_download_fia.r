##---------------------------------------------------------------
## 01_download_fia.r: Download all FIA data
##
## Author: Jacob Levine; jacob.levine@utah.edu
##---------------------------------------------------------------

##--------------------------------------------------------------
## 00. Load libraries                                                       
##--------------------------------------------------------------

library(rFIA)
library(dplyr)
library(here)

## increase timeout to 1 hour
options(timeout = 3600)

##--------------------------------------------------------------
## 01. Download FIA data                                                       
##--------------------------------------------------------------

states <- read.csv(here("data", "state_list.csv"))

for (state in states[,1]) {
  if (!file.exists(here("data", "fia", state))) {
    dir.create(here("data", "fia", state), showWarnings = FALSE)
  }
  getFIA(states = state, dir = file.path(here("data", "fia", state)), common = TRUE)
}
