rm(list=ls())

SAVE = FALSE
DRAW = FALSE

data.dir = "/Users/paleolab/babySTEPPS/Data/"
dump.dir = "/Users/paleolab/babySTEPPS/Dump/"
source("/Users/paleolab/babySTEPPS/Code/Format_Data.R")

model.dir = "/Users/paleolab/babySTEPPS/Code/"

n.adapt = 400
n.update = 400
n.iter = 2000

source("/Users/paleolab/babySTEPPS/Code/Run_Model.R")

source("/Users/paleolab/babySTEPPS/Code/Run_Model_Pick_Splines.R")


