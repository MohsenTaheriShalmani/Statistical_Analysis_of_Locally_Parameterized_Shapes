# Author: Mohsen Taheri, (2022)
# Email: mohsen.taherishalmani@uis.no
# Title: Simulation and hypothesis testing of LP-ds-rep and GP-ds-rep
# Place: University of Stavanger, Norway

#####################################################################################################
#####################################################################################################
# Libraries and data

library(shapes)
library(rgl)
library(matlib)
library(RiemBase)
library(data.table)
library(Directional)
library(truncnorm)
library(pracma)
library(numDeriv)
library(ggplot2)
library(dplyr)
library(rotations)
library(ks)
library(fields)
library(data.table)
library(numDeriv)

#clear the environment
remove(list=ls())

#set working directory to file location
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# load functions
if(TRUE){
  source("subFunctions/MathFunctions.R")
  source("subFunctions/euclideanization.R")
  source("subFunctions/normalsOfSkeletalSheetByTriangles.R")
  source("subFunctions/rotateFrameForwardAndBackward.R")
  source("subFunctions/ds_rep_Analysis.R")
}


# load initial ds-rep data for simulation
load("simulationData.RData")

#####################################################################################################
#####################################################################################################
# Main function

ds_rep_Analysis(simulationData=simulationData,
                typeOfAnalysis="LP_ds_rep", # choose "LP_ds_rep", "GP_ds_rep", or "EDMA"
                typeOfStudy="shapeAnalysis", # removing or preserving scale 
                typeOfStudy4directions="tangent space", # choose euclideanization method
                typeOfMeanDirection="Frechet", # choose type of mean direction
                typeOfTest="Parametric", # choose typeOfTest as "Parametric" or "Permutation"
                nSamplesG1=50, # number of samples for group 1
                nSamplesG2=50, # number of samples for group 2
                rotatingSpinalFramesId=c(19,22,25), # Three chosen spinal frames from 62,33,32,31,28,25,22,19,13,10,7,4,1,2,3,52
                rotateAboutWhichAxis=c(2,2,2), # choose the axis of rotation of the selected spinal frame. 1 for X-axis, 2 for Y-axis, and 3 or Z-axis.
                thetaG1=c(0,0,0), # thetaG1's elements are bending parameter for each group 1 (G1) spinal frame.
                thetaG2=c(-pi/12,-pi/12,-pi/12), # thetaG2's elements are bending parameter for each group 2 (G2) spinal frame.
                kappa1=100,  # kappa1 is for theta variation
                kappa2=600,  # kappa2 is to add noise to frame directions
                kappa3=250,  # kappa3 is to add noise to spoke directions
                kappa4=5000, # kappa4 is to add noise to connection directions. NB! adding large variation to connection direction makes ds-reps messy, be careful !!!
                noiseSDradii=0.3, #add noise to spoke lengths
                noiseSDConnectionLength=0.2, # add noise to connection lengths
                plotting=TRUE)

