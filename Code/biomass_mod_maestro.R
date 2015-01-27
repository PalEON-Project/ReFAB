rm(list=ls())

SAVE = TRUE
DRAW = TRUE

data.dir = "/Users/paleolab/babySTEPPS/Data/"
dump.dir = "/Users/paleolab/babySTEPPS/Dump/"
source("/Users/paleolab/babySTEPPS/Code/Format_Data.R")

model.dir = "/Users/paleolab/babySTEPPS/Code/"

n.adapt = 1
n.update = 1
n.iter = 2

source("/Users/paleolab/babySTEPPS/Code/Run_Model.R")

source("/Users/paleolab/babySTEPPS/Code/Run_Model_Pick_Splines.R")


