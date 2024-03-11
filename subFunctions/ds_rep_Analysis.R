
ds_rep_Analysis <- function(simulationData=simulationData,
                            typeOfAnalysis="LP_ds_rep", # choose "LP_ds_rep", "GP_ds_rep", or "EDMA"
                            typeOfStudy="shapeAnalysis", # removing or preserving scale 
                            typeOfStudy4directions="tangent space", # choose euclideanization method
                            typeOfMeanDirection="Frechet", # choose type of mean direction
                            typeOfTest="Parametric", # choose typeOfTest as "Parametric" or "Permutation"
                            nSamplesG1=100, # number of samples for group 1
                            nSamplesG2=100, # number of samples for group 2
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
                            plotting=TRUE) {
  
  #read data
  upSpoeksNumber<-simulationData$upSpoeksNumber
  downSpoeksNumber<-simulationData$downSpoeksNumber
  crestSpoksNumber<-simulationData$crestSpoksNumber
  nTotalRadii<-simulationData$nTotalRadii
  skelPointNo<-simulationData$skelPointNo
  skelRange<-simulationData$skelRange
  framesCenters<-simulationData$framesCenters
  framesParents<-simulationData$framesParents
  framesBackPoints<-simulationData$framesBackPoints
  framesFronts<-simulationData$framesFronts
  referenceFramesBasedOnParents<-simulationData$referenceFramesBasedOnParents
  spokesBasedOnFrames4deformation<-simulationData$spokesBasedOnFrames4deformation
  radii4deformation<-simulationData$radii4deformation
  connectionLength4deformation<-simulationData$connectionLength4deformation
  connectionsBasedOnParent4deformation<-simulationData$connectionsBasedOnParent4deformation
  numberOfFrames<-length(framesCenters)
  
  #simulation
  # simulate frames' directions
  # rotate selected spinal frames and add noise to all frames directions
  referenceFramesG1_temp<-referenceFramesBasedOnParents
  framesBasedOnParents_G1<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG1))
  I_tilde<-matrix(c(0,0,1,1,0,0,0,1,0),nrow = 3,byrow = T)
  for (j in 1:nSamplesG1) {
    for (k in 1:length(rotatingSpinalFramesId)) {
      frameIndexTemp<-rotatingSpinalFramesId[k]
      if(rotateAboutWhichAxis[k]==1){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG1[k]),sin(thetaG1[k])) , k=kappa1)
        vonMisesS1ToS2<-c(0,vonMises4S1[1],vonMises4S1[2])
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(0,1,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG1_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else if(rotateAboutWhichAxis[k]==2){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG1[k]),sin(thetaG1[k])) , k=kappa1)
        vonMisesS1ToS2<-c(vonMises4S1[1],0,vonMises4S1[2])
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(1,0,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG1_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else if(rotateAboutWhichAxis[k]==3){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG1[k]),sin(thetaG1[k])) , k=kappa1)
        vonMisesS1ToS2<-c(vonMises4S1[1],vonMises4S1[2],0)
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(1,0,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG1_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else{
        stop("Please define rotation axis of each frame as 1 or 2 or 3 !")
      }
    }
    for (i in 1:numberOfFrames) {
      tempFrame<-referenceFramesG1_temp[,,i]
      tempDirection1<-as.vector(rvmf(n = 1, mu = tempFrame[1,], k = kappa2))
      R1<-rotMat(tempFrame[1,],tempDirection1)
      tempFrame<-tempFrame%*%t(R1)
      tempDirection2<-as.vector(rvmf(n = 1, mu = tempFrame[2,], k = kappa2))
      R2<-rotMat(tempFrame[2,],tempDirection2)
      tempFrame<-tempFrame%*%t(R2)
      framesBasedOnParents_G1[,,i,j]<-tempFrame%*%t(R1)
    }
  }
  referenceFramesG2_temp<-referenceFramesBasedOnParents
  framesBasedOnParents_G2<-array(NA,dim = c(3,3,numberOfFrames,nSamplesG2))
  I_tilde<-matrix(c(0,0,1,1,0,0,0,1,0),nrow = 3,byrow = T)
  for (j in 1:nSamplesG2) {
    for (k in 1:length(rotatingSpinalFramesId)) {
      frameIndexTemp<-rotatingSpinalFramesId[k]
      if(rotateAboutWhichAxis[k]==1){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG2[k]),sin(thetaG2[k])) , k=kappa1)
        vonMisesS1ToS2<-c(0,vonMises4S1[1],vonMises4S1[2])
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(0,1,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG2_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else if(rotateAboutWhichAxis[k]==2){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG2[k]),sin(thetaG2[k])) , k=kappa1)
        vonMisesS1ToS2<-c(vonMises4S1[1],0,vonMises4S1[2])
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(1,0,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG2_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else if(rotateAboutWhichAxis[k]==3){
        vonMises4S1<-rvmf(n=1, mu = c(cos(thetaG2[k]),sin(thetaG2[k])) , k=kappa1)
        vonMisesS1ToS2<-c(vonMises4S1[1],vonMises4S1[2],0)
        myFrame<-referenceFramesBasedOnParents[,,frameIndexTemp]
        R1<-rotMat(myFrame[1,],c(0,0,1))
        rotatedFrame1<-myFrame%*%t(R1)
        R2<-rotMat(rotatedFrame1[2,],c(1,0,0))
        R3<-rotMat(c(1,0,0),vonMisesS1ToS2)
        rotatedMyFrame<-(I_tilde%*%t(R3))%*%solve(t(R2%*%R1))
        referenceFramesG2_temp[,,frameIndexTemp]<-rotatedMyFrame 
      }else{
        stop("Please define rotation axis of each frame as 1 or 2 or 3 !")
      }
    }
    for (i in 1:numberOfFrames) {
      tempFrame<-referenceFramesG2_temp[,,i]
      tempDirection1<-as.vector(rvmf(n = 1, mu = tempFrame[1,], k = kappa2))
      R1<-rotMat(tempFrame[1,],tempDirection1)
      tempFrame<-tempFrame%*%t(R1)
      tempDirection2<-as.vector(rvmf(n = 1, mu = tempFrame[2,], k = kappa2))
      R2<-rotMat(tempFrame[2,],tempDirection2)
      tempFrame<-tempFrame%*%t(R2)
      framesBasedOnParents_G2[,,i,j]<-tempFrame%*%t(R1)
    }
  }
  
  # simulate spokes' lengths
  lowerBound <- min(radii4deformation)/10
  upperBound <- mean(radii4deformation)*10
  radii_G1<-array(NA,dim = c(length(radii4deformation),nSamplesG1))
  for (j in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      radii_G1[i,j]<-rtruncnorm(n = 1, mean = radii4deformation[i],
                                a =lowerBound,b = upperBound ,sd = noiseSDradii)
    } 
  }
  radii_G2<-array(NA,dim = c(length(radii4deformation),nSamplesG2))
  for (j in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      radii_G2[i,j]<-rtruncnorm(n = 1, mean = radii4deformation[i],
                                a =lowerBound,b = upperBound ,sd = noiseSDradii)
      
    } 
  }
  
  # simulate connections' lengths
  lowerBound <- min(connectionLength4deformation)/10
  upperBound <- mean(connectionLength4deformation)*10
  connectionsLengths_G1<-array(NA,dim = c(length(connectionLength4deformation),nSamplesG1))
  for (j in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      connectionsLengths_G1[i,j]<-rtruncnorm(n = 1, mean = connectionLength4deformation[i],
                                             a =lowerBound,b = upperBound ,sd = noiseSDConnectionLength)
    } 
  }
  connectionsLengths_G2<-array(NA,dim = c(length(connectionLength4deformation),nSamplesG2))
  for (j in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      connectionsLengths_G2[i,j]<-rtruncnorm(n = 1, mean = connectionLength4deformation[i],
                                             a =lowerBound,b = upperBound ,sd = noiseSDConnectionLength)
      
    } 
  }
  
  # simulate spokes' directions
  spokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (i in 1:nSamplesG1) {
    for (j in 1:nTotalRadii) {
      spokesDirectionsBasedOnFrames_G1[j,,i]<-rvmf(n = 1, mu = spokesBasedOnFrames4deformation[j,], k = kappa3)
    }
  }
  spokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (i in 1:nSamplesG2) {
    for (j in 1:nTotalRadii) {
      spokesDirectionsBasedOnFrames_G2[j,,i]<-rvmf(n = 1, mu = spokesBasedOnFrames4deformation[j,], k = kappa3)
    }
  }
  
  # simulate connections' directions
  connectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  for (i in 1:nSamplesG1) {
    for (j in 1:numberOfFrames) {
      connectionsBasedOnParentFrames_G1[j,,i]<-rvmf(n = 1, mu = connectionsBasedOnParent4deformation[j,], k = kappa4)
    }
  }
  connectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  for (i in 1:nSamplesG2) {
    for (j in 1:numberOfFrames) {
      connectionsBasedOnParentFrames_G2[j,,i]<-rvmf(n = 1, mu = connectionsBasedOnParent4deformation[j,], k = kappa4)
    }
  }
  
  # convert LP-ds-reps to GP-ds-rep
  
  framesGlobalCoordinate_G1<-array(NA,dim = dim(framesBasedOnParents_G1))
  framesGlobalCoordinate_G1[,,framesCenters[1],]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (j in 1:nSamplesG1) {
    for (k in 2:numberOfFrames) {
      parent_Index<-framesParents[k]
      child_Index<-framesCenters[k]
      updatedParent<-framesGlobalCoordinate_G1[,,parent_Index,j]
      framesGlobalCoordinate_G1[,,child_Index,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                           vectorsInMainAxes = framesBasedOnParents_G1[,,child_Index,j])
    }
  }
  framesGlobalCoordinate_G2<-array(NA,dim = dim(framesBasedOnParents_G2))
  framesGlobalCoordinate_G2[,,framesCenters[1],]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
  for (j in 1:nSamplesG2) {
    for (k in 2:numberOfFrames) {
      parent_Index<-framesParents[k]
      child_Index<-framesCenters[k]
      updatedParent<-framesGlobalCoordinate_G2[,,parent_Index,j]
      framesGlobalCoordinate_G2[,,child_Index,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                           vectorsInMainAxes = framesBasedOnParents_G2[,,child_Index,j])
    }
  }
  
  spokesGlobalCoordinate_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (j in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesGlobalCoordinate_G1[i,,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G1[,,frameOfSpokeNo,j],
                                           vectorsInMainAxes = spokesDirectionsBasedOnFrames_G1[i,,j])
    }
  }
  spokesGlobalCoordinate_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (j in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      spokeNo<-i # 1<=spokeNo<=nTotalRadii
      frameOfSpokeNo<-NA
      if(spokeNo<=upSpoeksNumber){
        frameOfSpokeNo<-spokeNo
      }else if(spokeNo<=2*upSpoeksNumber){
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
      }
      
      spokesGlobalCoordinate_G2[i,,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G2[,,frameOfSpokeNo,j],
                                           vectorsInMainAxes = spokesDirectionsBasedOnFrames_G2[i,,j])
    }
  }
  
  
  connectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  for (j in 1:nSamplesG1) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      connectionsGlobalCoordinate_G1[k1,,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G1[,,k2,j],
                                           vectorsInMainAxes = connectionsBasedOnParentFrames_G1[k1,,j])
    }
  }
  connectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  for (j in 1:nSamplesG2) {
    for (i in 1:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      connectionsGlobalCoordinate_G2[k1,,j]<-
        rotateFrameToMainAxesAndRotateBack(myFrame = framesGlobalCoordinate_G2[,,k2,j],
                                           vectorsInMainAxes = connectionsBasedOnParentFrames_G2[k1,,j])
    }
  }
  
  
  positions_G1<-array(NA,dim = c(numberOfFrames,3,nSamplesG1))
  positions_G1[16,,]<-c(0,0,0)
  for (j in 1:nSamplesG1) {
    for (i in 2:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      positions_G1[k1,,j]<-positions_G1[k2,,j]+
        connectionsLengths_G1[k1,j]*connectionsGlobalCoordinate_G1[k1,,j]
    }
  }
  positions_G2<-array(NA,dim = c(numberOfFrames,3,nSamplesG2))
  positions_G2[16,,]<-c(0,0,0)
  for (j in 1:nSamplesG2) {
    for (i in 2:numberOfFrames) {
      k1<-framesCenters[i]
      k2<-framesParents[i]
      positions_G2[k1,,j]<-positions_G2[k2,,j]+
        connectionsLengths_G2[k1,j]*connectionsGlobalCoordinate_G2[k1,,j]
    }
  }
  
  spokesTails_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  spokesTips_G1<-array(NA,dim = c(nTotalRadii,3,nSamplesG1))
  for (j in 1:nSamplesG1) {
    for (i in 1:nTotalRadii) {
      frameOfSpokeNo<-NA
      if(i<=upSpoeksNumber){
        frameOfSpokeNo<-i
      }else if(i<=2*upSpoeksNumber){
        frameOfSpokeNo<-(i-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(i-upSpoeksNumber)
      }
      
      spokesTails_G1[i,,j]<-positions_G1[frameOfSpokeNo,,j]
      spokesTips_G1[i,,j]<-positions_G1[frameOfSpokeNo,,j]+
        spokesGlobalCoordinate_G1[i,,j]*radii_G1[i,j]
    }
  }
  spokesTails_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  spokesTips_G2<-array(NA,dim = c(nTotalRadii,3,nSamplesG2))
  for (j in 1:nSamplesG2) {
    for (i in 1:nTotalRadii) {
      frameOfSpokeNo<-NA
      if(i<=upSpoeksNumber){
        frameOfSpokeNo<-i
      }else if(i<=2*upSpoeksNumber){
        frameOfSpokeNo<-(i-upSpoeksNumber)
      }else{ #crest
        frameOfSpokeNo<-(i-upSpoeksNumber)
      }
      
      spokesTails_G2[i,,j]<-positions_G2[frameOfSpokeNo,,j]
      spokesTips_G2[i,,j]<-positions_G2[frameOfSpokeNo,,j]+
        spokesGlobalCoordinate_G2[i,,j]*radii_G2[i,j]
    }
  }
  
  
  #plot some samples of the simulated ds-reps
  if(plotting==TRUE){
    open3d()
    numberOfsamples4plot<-5 #choose numberOfsamples4plot from 1 to min(nSampleG1,nSampleG2)
    for (j in 1:numberOfsamples4plot) {
      srep1<-rbind(spokesTails_G1[,,j],spokesTips_G1[,,j])
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.5,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
    }
    for (j in 1:numberOfsamples4plot) {
      srep1<-rbind(spokesTails_G2[,,j],spokesTips_G2[,,j])
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.5,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
    }
    decorate3d(xlab = "x", ylab = "y", zlab = "z",
               # xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
               box = TRUE, axes = TRUE, main = "Five samples of each group", sub = NULL,
               top = TRUE, aspect = FALSE, expand = 1.1)
  }
  
  if(typeOfAnalysis=="LP_ds_rep"){
    
    # calculate LP sizes
    LP_sizes_G1<-rep(NA,nSamplesG1)
    for (i in 1:nSamplesG1) {
      LP_sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[,i]),log(radii_G1[,i]))))
    }
    LP_sizes_G2<-rep(NA,nSamplesG2)
    for (i in 1:nSamplesG2) {
      LP_sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[,i]),log(radii_G2[,i]))))
    }
    
    #Removing or preserving the scale by LP-size
    if(typeOfStudy=="sizeAndShapeAnalysis"){
      
      # sizes_G1 and sizes_G2 are LP size but we use them here for plot and scaling
      
      sizes_G1<-rep(1,nSamplesG1)
      sizes_G2<-rep(1,nSamplesG2)
      
      radiiScaled_G1<-radii_G1 #we don't have scaling in size-and-shape analysis
      radiiScaled_G2<-radii_G2
      
      connectionsLengthsScaled_G1<-connectionsLengths_G1
      connectionsLengthsScaled_G2<-connectionsLengths_G2
      
    }else if(typeOfStudy=="shapeAnalysis"){
      sizes_G1<-rep(NA,nSamplesG1)
      for (i in 1:nSamplesG1) {
        sizes_G1[i]<-exp(mean(c(log(connectionsLengths_G1[,i]),log(radii_G1[,i]))))
      }
      sizes_G2<-rep(NA,nSamplesG2)
      for (i in 1:nSamplesG2) {
        sizes_G2[i]<-exp(mean(c(log(connectionsLengths_G2[,i]),log(radii_G2[,i]))))
      }
      
      radiiScaled_G1<-sweep(radii_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
      radiiScaled_G2<-sweep(radii_G2, 2, sizes_G2, "/") 
      
      connectionsLengthsScaled_G1<-sweep(connectionsLengths_G1, 2, sizes_G1, "/") #2 indicate operation "/" on columns
      connectionsLengthsScaled_G2<-sweep(connectionsLengths_G2, 2, sizes_G2, "/") 
      
    }
    
    #calculate mean frames in local and global coordinate systems
    if(TRUE){
      
      framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
      for (i in 1:nSamplesG1) {
        for (k in 1:numberOfFrames) {
          framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
        }
      }
      framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
      for (i in 1:nSamplesG2) {
        for (k in 1:numberOfFrames) {
          framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
        }
      }
      
      meanFramesBasedOnParents_G1<-array(NA, dim = c(3,3,numberOfFrames))
      for (k in framesCenters) {
        if(k==16){
          meanFramesBasedOnParents_G1[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
        }else{
          tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G1[k,,])),type = 'geometric')
          meanFramesBasedOnParents_G1[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
        }
      }
      meanFramesBasedOnParents_G2<-array(NA, dim = c(3,3,numberOfFrames))
      for (k in framesCenters) {
        if(k==16){
          meanFramesBasedOnParents_G2[,,k]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
        }else{
          tempVec<-mean(as.SO3(t(framesBasedOnParentsVectorized_G2[k,,])),type = 'geometric')
          meanFramesBasedOnParents_G2[,,k]<-matrix(tempVec,nrow = 3,byrow = TRUE)
        }
      }
      
      meanFramesGlobalCoordinate_G1<-array(NA,dim = dim(meanFramesBasedOnParents_G1))
      meanFramesGlobalCoordinate_G1[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
      for (k in 2:numberOfFrames) {
        parent_Index<-framesParents[k]
        child_Index<-framesCenters[k]
        updatedParent<-meanFramesGlobalCoordinate_G1[,,parent_Index]
        meanFramesGlobalCoordinate_G1[,,child_Index]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                             vectorsInMainAxes = meanFramesBasedOnParents_G1[,,child_Index])
      }
      meanFramesGlobalCoordinate_G2<-array(NA,dim = dim(meanFramesBasedOnParents_G2))
      meanFramesGlobalCoordinate_G2[,,framesCenters[1]]<-rbind(c(0,0,1),c(1,0,0),c(0,1,0))
      for (k in 2:numberOfFrames) {
        parent_Index<-framesParents[k]
        child_Index<-framesCenters[k]
        updatedParent<-meanFramesGlobalCoordinate_G2[,,parent_Index]
        meanFramesGlobalCoordinate_G2[,,child_Index]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = updatedParent,
                                             vectorsInMainAxes = meanFramesBasedOnParents_G2[,,child_Index])
      }
      print("Calculation of the mean frames in GCS is done!")
    }
    
    # Calculate mean spokes' directions based on frames
    if(TRUE){
      meanSpokesDirectionsBasedOnFrames_G1<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        # For extremely concentrated data we use Mardia mean direction 
        pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G1[i,,]))
        if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
          meanSpokesDirectionsBasedOnFrames_G1[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G1[i,,])))
        }else if(typeOfMeanDirection=="Frechet"){
          meanSpokesDirectionsBasedOnFrames_G1[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G1[i,,]) 
        }else if(typeOfMeanDirection=="PNS"){
          sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G1[i,,])
          meanSpokesDirectionsBasedOnFrames_G1[i,]<-pns(spokesDirectionsBasedOnFrames_G1[i,,],sphere.type = sphereType)$PNS$mean 
        }else{
          stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
        }
      }
      meanSpokesDirectionsBasedOnFrames_G2<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        # For extremely concentrated data we use Mardia mean direction 
        pcaTemp<-prcomp(t(spokesDirectionsBasedOnFrames_G2[i,,]))
        if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
          meanSpokesDirectionsBasedOnFrames_G2[i,]<-convertVec2unitVec(colMeans(t(spokesDirectionsBasedOnFrames_G2[i,,])))
        }else if(typeOfMeanDirection=="Frechet"){
          meanSpokesDirectionsBasedOnFrames_G2[i,]<-frechetMean(spokesDirectionsBasedOnFrames_G2[i,,])
        }else if(typeOfMeanDirection=="PNS"){
          sphereType<-kurtosisTestFunction(spokesDirectionsBasedOnFrames_G2[i,,])
          meanSpokesDirectionsBasedOnFrames_G2[i,]<-pns(spokesDirectionsBasedOnFrames_G2[i,,],sphere.type = sphereType)$PNS$mean 
        }else{
          stop("Please specify the typeOfMeanDirection by PNS or Frechet!")
        }
      }
      print("Calculation of the mean directions is done!")
    }
    
    # Calculate mean spokes' directions based on global coordinate (using mean frames in global coordinate)
    if(TRUE){
      meanSpokesDirectionsGlobalCoordinate_G1<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        spokeNo<-i # 1<=spokeNo<=nTotalRadii
        frameOfSpokeNo<-NA
        if(spokeNo<=upSpoeksNumber){
          frameOfSpokeNo<-spokeNo
        }else if(spokeNo<=2*upSpoeksNumber){
          frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
        }else{ #crest
          frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
        }
        meanSpokesDirectionsGlobalCoordinate_G1[i,]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,frameOfSpokeNo],
                                             vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G1[i,])
        
      }
      meanSpokesDirectionsGlobalCoordinate_G2<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        spokeNo<-i # 1<=spokeNo<=nTotalRadii
        frameOfSpokeNo<-NA
        if(spokeNo<=upSpoeksNumber){
          frameOfSpokeNo<-spokeNo
        }else if(spokeNo<=2*upSpoeksNumber){
          frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
        }else{ #crest
          frameOfSpokeNo<-(spokeNo-upSpoeksNumber)
        }
        meanSpokesDirectionsGlobalCoordinate_G2[i,]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,frameOfSpokeNo],
                                             vectorsInMainAxes = meanSpokesDirectionsBasedOnFrames_G2[i,])
        
      }
    }
    
    # Calculate geometric mean of spokes' lengths
    if(TRUE){
      radiiMean_G1<-exp(rowMeans(log(radiiScaled_G1)))
      radiiMean_G2<-exp(rowMeans(log(radiiScaled_G2)))
    }
    
    # Calculate geometric mean of connections' lengths
    if(TRUE){
      meanConnectionsLengths_G1<-exp(rowMeans(log(connectionsLengthsScaled_G1)))
      meanConnectionsLengths_G2<-exp(rowMeans(log(connectionsLengthsScaled_G2)))
    }
    
    # Calculate mean connection directions based on frames
    if(TRUE){
      meanConnectionsBasedOnParentFrames_G1<-array(NA,dim = c(numberOfFrames,3))
      for (i in 1:numberOfFrames) {
        if(i==framesCenters[1]){
          meanConnectionsBasedOnParentFrames_G1[i,]<-c(0,0,0)
        }else{
          # For extremely concentrated data we use Mardia mean direction 
          pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G1[i,,]))
          if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
            meanConnectionsBasedOnParentFrames_G1[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G1[i,,])))
          }else if(typeOfMeanDirection=="Frechet"){
            meanConnectionsBasedOnParentFrames_G1[i,]<-frechetMean(connectionsBasedOnParentFrames_G1[i,,]) 
          }else if(typeOfMeanDirection=="PNS"){
            sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G1[i,,])
            meanConnectionsBasedOnParentFrames_G1[i,]<-pns(connectionsBasedOnParentFrames_G1[i,,],sphere.type = sphereType)$PNS$mean  
          }else{
            stop("Please specify the typeOfMeanDirection by PNS or Frechet")
          }
        }
      }
      meanConnectionsBasedOnParentFrames_G2<-array(NA,dim = c(numberOfFrames,3))
      for (i in 1:numberOfFrames) {
        if(i==framesCenters[1]){
          meanConnectionsBasedOnParentFrames_G2[i,]<-c(0,0,0)
        }else{
          # For extremely concentrated data we use Mardia mean direction 
          pcaTemp<-prcomp(t(connectionsBasedOnParentFrames_G2[i,,]))
          if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
            meanConnectionsBasedOnParentFrames_G2[i,]<-convertVec2unitVec(colMeans(t(connectionsBasedOnParentFrames_G2[i,,])))
          }else if(typeOfMeanDirection=="Frechet"){
            meanConnectionsBasedOnParentFrames_G2[i,]<-frechetMean(connectionsBasedOnParentFrames_G2[i,,]) 
          }else if(typeOfMeanDirection=="PNS"){
            sphereType<-kurtosisTestFunction(connectionsBasedOnParentFrames_G2[i,,])
            meanConnectionsBasedOnParentFrames_G2[i,]<-pns(connectionsBasedOnParentFrames_G2[i,,],sphere.type = sphereType)$PNS$mean  
          }else{
            stop("Please specify the typeOfMeanDirection by PNS or Frechet")
          }
        }
      }
      print("Calculation of the mean connections is done!")
    }
    
    # Calculate mean connection based on global coordinate (using mean frames in global coordinate)
    if(TRUE){
      meanConnectionsGlobalCoordinate_G1<-array(NA,dim = c(numberOfFrames,3))
      for (i in 1:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        meanConnectionsGlobalCoordinate_G1[k1,]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G1[,,k2],
                                             vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G1[k1,])
      }
      meanConnectionsGlobalCoordinate_G2<-array(NA,dim = c(numberOfFrames,3))
      for (i in 1:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        meanConnectionsGlobalCoordinate_G2[k1,]<-
          rotateFrameToMainAxesAndRotateBack(myFrame = meanFramesGlobalCoordinate_G2[,,k2],
                                             vectorsInMainAxes = meanConnectionsBasedOnParentFrames_G2[k1,])
      }
    }
    
    # Convert mean LP-ds-rep to GP-ds-rep
    if(TRUE){
      meanPositions_G1<-array(NA,dim = c(numberOfFrames,3))
      meanPositions_G1[16,]<-c(0,0,0)
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        meanPositions_G1[k1,]<-meanPositions_G1[k2,]+
          meanConnectionsLengths_G1[k1]*meanConnectionsGlobalCoordinate_G1[k1,]
        
      }
      meanPositions_G2<-array(NA,dim = c(numberOfFrames,3))
      meanPositions_G2[16,]<-c(0,0,0)
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        meanPositions_G2[k1,]<-meanPositions_G2[k2,]+
          meanConnectionsLengths_G2[k1]*meanConnectionsGlobalCoordinate_G2[k1,]
        
      }
      meanSpokesTails_G1<-array(NA,dim = c(nTotalRadii,3))
      meanSpokesTips_G1<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        frameOfSpokeNo<-NA
        if(i<=upSpoeksNumber){
          frameOfSpokeNo<-i
        }else if(i<=2*upSpoeksNumber){
          frameOfSpokeNo<-(i-upSpoeksNumber)
        }else{ #crest
          frameOfSpokeNo<-(i-upSpoeksNumber)
        }
        
        meanSpokesTails_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]
        meanSpokesTips_G1[i,]<-meanPositions_G1[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G1[i,]*radiiMean_G1[i]
      }
      meanSpokesTails_G2<-array(NA,dim = c(nTotalRadii,3))
      meanSpokesTips_G2<-array(NA,dim = c(nTotalRadii,3))
      for (i in 1:nTotalRadii) {
        frameOfSpokeNo<-NA
        if(i<=upSpoeksNumber){
          frameOfSpokeNo<-i
        }else if(i<=2*upSpoeksNumber){
          frameOfSpokeNo<-(i-upSpoeksNumber)
        }else{ #crest
          frameOfSpokeNo<-(i-upSpoeksNumber)
        }
        
        meanSpokesTails_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]
        meanSpokesTips_G2[i,]<-meanPositions_G2[frameOfSpokeNo,]+meanSpokesDirectionsGlobalCoordinate_G2[i,]*radiiMean_G2[i]
      }
      print("Calculation of the mean GP_ds_rep is done!")
    }
    
    # Plot overlaid LP-ds-rep means of G1 and G2
    if(plotting==TRUE){
      open3d()
      srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1) #we scale back to the original size by *mean(sizes_G1)
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      srep2<-rbind(meanSpokesTails_G2,meanSpokesTips_G2)*mean(sizes_G2)
      for (i in 1:nTotalRadii) {
        plot3d(srep2[c(i,(i+nTotalRadii)),],type="l",lwd = 1,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      # legend3d("topright", legend = paste(c('Mean G1', 'Mean G2')), pch = 16, col = c("blue","red"), cex=1, inset=c(0.02))
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 # xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
                 box = TRUE, axes = TRUE, main = "overlaid mean LP-ds-reps", sub = NULL,
                 top = TRUE, aspect = FALSE, expand = 1.1)
    }
    
    
    ##############################
    # LP-ds-rep Hypothesis testing 
    
    # hypothesis test on LP size
    pValues_LP_sizes<-meanDifferenceTest1D(log(LP_sizes_G1),log(LP_sizes_G2),type = typeOfTest) 
    cat("pValue of LP sizes is:",pValues_LP_sizes,"\n")
    boxplot(LP_sizes_G1, LP_sizes_G2, names = c("G1","G2"),main="LP size")
    cat("sd LP size G1:",sd(LP_sizes_G1),"sd LP size G2:",sd(LP_sizes_G2),"\n")
    cat("Mean LP size G1:",mean(LP_sizes_G1),"mean LP size G2:",mean(LP_sizes_G2),"\n")
    
    # hypothesis test on spokes' lengths
    pValues_TtestRadii<-rep(NA,nTotalRadii)
    pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
    for (i in 1:nTotalRadii) {
      setTxtProgressBar(pb, i) #create progress bar
      
      pValues_TtestRadii[i]<-meanDifferenceTest1D(log(radiiScaled_G1[i,]),
                                                  log(radiiScaled_G2[i,]),
                                                  type = typeOfTest)
    }
    # which(pValues_TtestRadii<=0.05)
    
    
    # hypothesis test on connections' length
    pValues_TtestConnectionsLengths<-rep(NA,numberOfFrames)
    pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
    for (i in 1:numberOfFrames) {
      setTxtProgressBar(pb, i) #create progress bar
      if(i==16){
        pValues_TtestConnectionsLengths[i]<-1
      }else{
        pValues_TtestConnectionsLengths[i]<-meanDifferenceTest1D(log(connectionsLengthsScaled_G1[i,]),
                                                                 log(connectionsLengthsScaled_G2[i,]),
                                                                 type = typeOfTest)
      }
    }
    # which(pValues_TtestConnectionsLengths<=0.05)
    
    # hypothesis test on spokes' directions based on local frames
    euclideanizedSpokesDirBasedOnFramesG1<-array(NA,dim = c(nSamplesG1,2,nTotalRadii))
    euclideanizedSpokesDirBasedOnFramesG2<-array(NA,dim = c(nSamplesG2,2,nTotalRadii))
    pValspokesDirectionsBasedOnFrames<-rep(NA,nTotalRadii)
    pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
    for(i in 1:nTotalRadii){
      setTxtProgressBar(pb, i) #create progress bar
      #NB! euclideanization must contain two groups because it uses the pooled mean 
      euclideanizedTemp<-euclideanization(spokesDirectionsBasedOnFrames_G1[i,,],
                                          spokesDirectionsBasedOnFrames_G2[i,,],
                                          type = typeOfStudy4directions)
      
      pValspokesDirectionsBasedOnFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                           euclideanizedTemp$euclideanG2,
                                                                           type=typeOfTest)
      
      euclideanizedSpokesDirBasedOnFramesG1[,,i]<-euclideanizedTemp$euclideanG1
      euclideanizedSpokesDirBasedOnFramesG2[,,i]<-euclideanizedTemp$euclideanG2
      
    }
    # which(pValspokesDirectionsBasedOnFrames<=0.05)
    
    
    # hypothesis test on connections' directions based on local frames
    euclideanizedConnectionsBasedOnParentFramesG1<-array(0,dim = c(nSamplesG1,2,numberOfFrames))
    euclideanizedConnectionsBasedOnParentFramesG2<-array(0,dim = c(nSamplesG2,2,numberOfFrames))
    pValConnectionsBasedOnParentFrames<-rep(NA,numberOfFrames)
    pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
    for(i in 1:numberOfFrames){
      setTxtProgressBar(pb, i) #create progress bar
      if(i==16){
        pValConnectionsBasedOnParentFrames[i]<-1
        next
      }
      
      euclideanizedTemp<-euclideanization(connectionsBasedOnParentFrames_G1[i,,],
                                          connectionsBasedOnParentFrames_G2[i,,],
                                          type = typeOfStudy4directions)
      
      pValConnectionsBasedOnParentFrames[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                            euclideanizedTemp$euclideanG2,
                                                                            type=typeOfTest)
      
      euclideanizedConnectionsBasedOnParentFramesG1[,,i]<-euclideanizedTemp$euclideanG1
      euclideanizedConnectionsBasedOnParentFramesG2[,,i]<-euclideanizedTemp$euclideanG2
      
    }
    # which(pValConnectionsBasedOnParentFrames<=0.05)
    
    framesBasedOnParentsVectorized_G1<-array(NA,dim = c(numberOfFrames,9,nSamplesG1))
    for (i in 1:nSamplesG1) {
      for (k in 1:numberOfFrames) {
        framesBasedOnParentsVectorized_G1[k,,i]<-as.vector(t(framesBasedOnParents_G1[,,k,i]))
      }
    }
    framesBasedOnParentsVectorized_G2<-array(NA,dim = c(numberOfFrames,9,nSamplesG2))
    for (i in 1:nSamplesG2) {
      for (k in 1:numberOfFrames) {
        framesBasedOnParentsVectorized_G2[k,,i]<-as.vector(t(framesBasedOnParents_G2[,,k,i]))
      }
    }
    
    # hypothesis test on frames' normal directions based on parent frames
    euclideanizedFrameBasedOnParentG1<-array(0,dim = c(nSamplesG1,3,numberOfFrames))
    euclideanizedFrameBasedOnParentG2<-array(0,dim = c(nSamplesG2,3,numberOfFrames))
    pValFramesBasedOnParent<-rep(NA,numberOfFrames)
    pb <- txtProgressBar(min = 0, max = numberOfFrames, style = 3) #progress bar
    for(i in 1:numberOfFrames){
      setTxtProgressBar(pb, i) #create progress bar
      if(i==16){
        pValFramesBasedOnParent[i]<-1
        next
      }
      
      # conversion to unit quaternions
      Q4Temp<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G1[i,,])))
      Q4Temp2<-matrix(as.numeric(t(Q4Temp)),ncol = 4,byrow = TRUE)
      for (j in 1:dim(Q4Temp2)[1]) {
        Q4Temp2[j,]<-Q4Temp2[j,]/norm(Q4Temp2[j,],type = '2')
      }
      Q4Temp3<-as.Q4(as.SO3(t(framesBasedOnParentsVectorized_G2[i,,])))
      Q4Temp4<-matrix(as.numeric(t(Q4Temp3)),ncol = 4,byrow = TRUE)
      for (j in 1:dim(Q4Temp4)[1]) {
        Q4Temp4[j,]<-Q4Temp4[j,]/norm(Q4Temp4[j,],type = '2')
      }
      
      euclideanizedTemp<-euclideanization(t(Q4Temp2),t(Q4Temp4),type = typeOfStudy4directions)
      
      
      pValFramesBasedOnParent[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                                 euclideanizedTemp$euclideanG2,
                                                                 type=typeOfTest)
      
      euclideanizedFrameBasedOnParentG1[,,i]<-euclideanizedTemp$euclideanG1
      euclideanizedFrameBasedOnParentG2[,,i]<-euclideanizedTemp$euclideanG2
      
    }
    # which(pValFramesBasedOnParent<=0.05)
    
    
    #################################
    # plot significant LP-ds-rep GOPs
    
    pvalues_LP_ds_rep <- c(pValues_TtestRadii,                  #n_s = 122
                           pValues_TtestConnectionsLengths,     #n_f = 71
                           pValspokesDirectionsBasedOnFrames,   #n_s = 122
                           pValConnectionsBasedOnParentFrames,  #n_f = 71
                           pValFramesBasedOnParent,            #n_f = 71
                           pValues_LP_sizes)
    
    n_s<-nTotalRadii
    n_f<-numberOfFrames
    
    alpha<-0.05
    significantPvalues<-which(pvalues_LP_ds_rep<=alpha)
    significantPvalues
    
    #adjust p-values by Benjamini-Hochberg
    FDR<-0.05
    pvalues_LP_ds_rep_BH<-p.adjust(pvalues_LP_ds_rep,method = "BH")
    pvalues_LP_ds_rep_Bonferroni<-p.adjust(pvalues_LP_ds_rep,method = "bonferroni")
    significantPvalues_BH<-which(pvalues_LP_ds_rep_BH<=FDR)
    significantPvalues_BH
    significantPvalues_Bonferroni<-which(pvalues_LP_ds_rep_Bonferroni<=FDR)
    significantPvalues_Bonferroni
    
    
    # plot by ggplot
    if(plotting==TRUE){
      df_LP <- data.frame(Type=c(rep("Raw p-value",length(pvalues_LP_ds_rep)),
                                 rep("Bonferroni",length(pvalues_LP_ds_rep)),
                                 rep("BH",length(pvalues_LP_ds_rep))),
                          ordereOfPvalues=1:length(pvalues_LP_ds_rep),
                          Values=c(sort(pvalues_LP_ds_rep),sort(pvalues_LP_ds_rep_Bonferroni),sort(pvalues_LP_ds_rep_BH)))
      p<-ggplot(df_LP, aes(x=ordereOfPvalues, y=Values, group=Type))
      p + geom_line(aes(linetype=Type),size=1)+
        # geom_line(aes(linetype=Type, col=Type),size=1)+
        # geom_point(aes(shape=Type),alpha=0.7)+
        geom_hline(yintercept=0.05,linetype="solid", color = "red")+
        scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
        theme_bw()+
        theme(plot.title = element_text(size = 17, hjust = 0.5),
              legend.text=element_text(size=17),
              legend.title=element_blank(),
              # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              # legend.title=element_text(size=17),
              legend.position="bottom",
              axis.title=element_text(size=17),
              legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
        guides(colour = guide_legend(title.hjust = 0.5))+
        xlab("Ranking of p-values") + ylab("p-values")
      
      }
    
    #1
    significantRadii<-significantPvalues[which(significantPvalues<=n_s)]
    significantRadii
    significantRadii_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_s)]
    significantRadii_BH
    #2
    significantConnectionsLengths<-significantPvalues[which(n_s+1<=significantPvalues
                                                            & significantPvalues<=(n_s+n_f))]-n_s
    significantConnectionsLengths
    significantConnectionsLengths_BH<-significantPvalues_BH[which((n_s+1)<=significantPvalues_BH
                                                                  & significantPvalues_BH<=(n_s+n_f))]-n_s
    significantConnectionsLengths_BH
    #3
    significantspokesDirections<-significantPvalues[which((n_s+n_f+1)<=significantPvalues
                                                          & significantPvalues<=(2*n_s+n_f))]-(n_s+n_f)
    significantspokesDirections
    significantspokesDirections_BH<-significantPvalues_BH[which((n_s+n_f+1)<=significantPvalues_BH
                                                                & significantPvalues_BH<=(2*n_s+n_f))]-(n_s+n_f)
    significantspokesDirections_BH
    #4
    significantConnectionsDirections<-significantPvalues[which((2*n_s+n_f+1)<=significantPvalues
                                                               & significantPvalues<=(2*n_s+2*n_f))]-(2*n_s+n_f)
    significantConnectionsDirections
    significantConnectionsDirections_BH<-significantPvalues_BH[which((2*n_s+n_f+1)<=significantPvalues_BH
                                                                     & significantPvalues_BH<=(2*n_s+2*n_f))]-(2*n_s+n_f)
    significantConnectionsDirections_BH
    #5
    significantFrame<-significantPvalues[which((2*n_s+2*n_f+1)<=significantPvalues
                                               & significantPvalues<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
    significantFrame
    significantFrame_BH<-significantPvalues_BH[which((2*n_s+2*n_f+1)<=significantPvalues_BH
                                                     & significantPvalues_BH<=(2*n_s+3*n_f))]-(2*n_s+2*n_f)
    significantFrame_BH
    
    
    #plot significant GOPs before and after the BH adjustment
    if(plotting==TRUE){
      #1 plot
      srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' lengths",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantRadii_BH) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' lengths after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #2 plot
      open3d()
      skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
      }
      for (i in significantConnectionsLengths) {
        k1<-i
        k2<-framesParents[which(framesCenters==i)]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant connections' lengths",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
      }
      for (i in significantConnectionsLengths_BH) {
        k1<-i
        k2<-framesParents[which(framesCenters==i)]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant connections' lengths after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #3 plot
      srep1<-rbind(meanSpokesTails_G1,meanSpokesTips_G1)*mean(sizes_G1)
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantspokesDirections) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spoks' directions",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantspokesDirections_BH) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' directions after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #4 plot
      open3d()
      skelG1_1<-meanSpokesTails_G1[skelRange,]*mean(sizes_G1)
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
      }
      for (i in significantConnectionsDirections) {
        k1<-i
        k2<-framesParents[which(framesCenters==i)]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant connections' directions",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 2:numberOfFrames) {
        k1<-framesCenters[i]
        k2<-framesParents[i]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="blue", lwd=1)
      }
      for (i in significantConnectionsDirections_BH) {
        k1<-i
        k2<-framesParents[which(framesCenters==i)]
        vectors3d(skelG1_1[k1,],origin = skelG1_1[k2,],headlength = 0.1,radius = 1/6, col="red", lwd=3)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant connections' directions after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #5 plot
      open3d()
      skeletalSheet<-meanPositions_G1*mean(sizes_G1)
      for (i in 2:numberOfFrames) {
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=0.8)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=0.8)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=0.8)
      }
      for (i in significantFrame) {
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant frames",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 2:numberOfFrames) {
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="darkblue", lwd=0.8)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="blue", lwd=0.8)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="lightblue", lwd=0.8)
      }
      for (i in significantFrame_BH) {
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[1,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[2,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
        vectors3d(skeletalSheet[i,]+meanFramesGlobalCoordinate_G1[3,,i],origin = skeletalSheet[i,],headlength = 0.1,radius = 1/10, col="red", lwd=5)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant frames after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
    }
    
    
  }else if(typeOfAnalysis=="GP_ds_rep"){
    
    # combine spokes' information
    nTotalRadii<-dim(spokesTails_G1)[1]
    boundaryPlusSkeletal_G1<-array(NA,dim=c(2*nTotalRadii,3,nSamplesG1))
    boundaryPlusSkeletal_G2<-array(NA,dim=c(2*nTotalRadii,3,nSamplesG2))
    for (i in 1:nSamplesG1) {
      boundaryPlusSkeletal_G1[,,i]<-rbind(spokesTails_G1[,,i],spokesTips_G1[,,i])
    }
    for (i in 1:nSamplesG2) {
      boundaryPlusSkeletal_G2[,,i]<-rbind(spokesTails_G2[,,i],spokesTips_G2[,,i])
    }
    
    # Alignment by GPA based on spokes' tips and tails
    if(typeOfStudy=="shapeAnalysis"){
      # sreps alignment by spokes' tips and tails
      all_G1_G2<-abind(boundaryPlusSkeletal_G1,boundaryPlusSkeletal_G2)
      procAllSkeletalPlusBoundary<-procGPA(all_G1_G2, scale = TRUE)   # scale = TRUE for shape analysis
      meanOfAll_G1_G2<-procAllSkeletalPlusBoundary$mshape
    }else if(typeOfStudy=="sizeAndShapeAnalysis"){
      # sreps alignment by spokes' tips and tails
      all_G1_G2<-abind(boundaryPlusSkeletal_G1,boundaryPlusSkeletal_G2)
      procAllSkeletalPlusBoundary<-procGPA(all_G1_G2, scale = FALSE)   # scale = TRUE for shape analysis
      meanOfAll_G1_G2<-procAllSkeletalPlusBoundary$mshape
    }else{
      stop("Please choose the type of study as shapeAnalysis or sizeAndShapeAnalysis!")
    }
    
    #pooled group
    alignedAllSkeletalPlusBoundary<-procAllSkeletalPlusBoundary$rotated  
    
    # plot all spokes' tips and tails of CG ds-reps after the alignment
    if(plotting==TRUE){
      open3d()
      for (i in 1:nSamplesG1) {
        plot3d(alignedAllSkeletalPlusBoundary[,,i],type="p",col = "red" ,expand = 10,box=TRUE,add = TRUE)
      }
      for (i in (nSamplesG1+1):(nSamplesG1+nSamplesG2)) {
        plot3d(alignedAllSkeletalPlusBoundary[,,i],type="p",col = "blue" ,expand = 10,box=TRUE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 # xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
                 box = TRUE, axes = TRUE, main = "Spokes' tips and tailes after alignment", sub = NULL,
                 top = TRUE, aspect = FALSE, expand = 1.1)
    }
    
    # extract skeletal and boundary after the alignment
    alignedSkeletalPlusBoundaryG1<-alignedAllSkeletalPlusBoundary[,,1:nSamplesG1]  #G1 group
    alignedSkeletalPlusBoundaryG2<-alignedAllSkeletalPlusBoundary[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)] #G2 group
    
    skeletalAlignedG1<-alignedSkeletalPlusBoundaryG1[1:nTotalRadii,,]
    boundaryAlignedG1<-alignedSkeletalPlusBoundaryG1[(nTotalRadii+1):(2*nTotalRadii),,]
    skeletalAlignedG2<-alignedSkeletalPlusBoundaryG2[1:nTotalRadii,,]
    boundaryAlignedG2<-alignedSkeletalPlusBoundaryG2[(nTotalRadii+1):(2*nTotalRadii),,]
    
    # Plot GP-ds-rep
    # Note that after the alignment centroid of the objects would be at 
    # the origin of the coordinate system. Thus GP-ds-rep skeletal positions (i.e., p_i) are 
    # vectors with tail at (0,0,0)^T
    
    #plot spokes
    sampleNo<-1 #choose sampleNo between 1 to nSamplesG1=108 to see other ds-reps
    if(plotting==TRUE){
      #plot spokes
      open3d()
      for (i in 1:nTotalRadii) {
        vectors3d(alignedSkeletalPlusBoundaryG1[i,,sampleNo],origin = c(0,0,0),headlength = 0.05,radius = 1/10, col="blue", lwd=1)
        vectors3d(alignedSkeletalPlusBoundaryG1[i+nTotalRadii,,sampleNo],origin = alignedSkeletalPlusBoundaryG1[i,,sampleNo],headlength = 0.1,radius = 1/10, col="red", lwd=1)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "First sample of the first group",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
    }
    
    # Calculate spokes' lengths (i.e., radii) after the alignment (or removing scale)
    radiiG1<-array(NA, dim=c(nTotalRadii,nSamplesG1))
    for (k in 1:nSamplesG1) {
      for (i in 1:nTotalRadii) {
        radiiG1[i,k]<-norm(skeletalAlignedG1[i,,k]-boundaryAlignedG1[i,,k],type = "2")
      }
    }
    radiiG2<-array(NA, dim=c(nTotalRadii,nSamplesG2))
    for (k in 1:nSamplesG2) {
      for (i in 1:nTotalRadii) {
        radiiG2[i,k]<-norm(skeletalAlignedG2[i,,k]-boundaryAlignedG2[i,,k],type = "2")
      }
    }
    
    # Calculate spokes' directions after the alignment (or removing scale)
    DirectionsG1<-array(NA, dim=c(nTotalRadii,3,nSamplesG1))
    for (k in 1:nSamplesG1) {
      for (i in 1:nTotalRadii) {
        DirectionsG1[i,,k]<-convertVec2unitVec(boundaryAlignedG1[i,,k]-skeletalAlignedG1[i,,k])
      }
    }
    DirectionsG2<-array(NA, dim=c(nTotalRadii,3,nSamplesG2))
    for (k in 1:nSamplesG2) {
      for (i in 1:nTotalRadii) {
        DirectionsG2[i,,k]<-convertVec2unitVec(boundaryAlignedG2[i,,k]-skeletalAlignedG2[i,,k])
      }
    }
    
    #remove duplicated up and down skeltal positions
    skelG1<-skeletalAlignedG1[skelRange,,]
    skelG2<-skeletalAlignedG2[skelRange,,]
    
    #############################
    # Mean GP-ds-rep of CG and PD
    
    # mean of skeletal PDM by GPA
    procSkelG1<-procGPA(skelG1, scale = F)  #scale F because we managed the scaling previously
    skelMeanG1<-procSkelG1$mshape           #mean of skeletal of G1
    procSkelG2<-procGPA(skelG2, scale = F)  #scale F because we managed the scaling previously
    skelMeanG2<-procSkelG2$mshape           #mean of skeletal of G2
    
    
    directionMeanG1<-array(NA,dim = c(nTotalRadii,3))             
    for (i in 1:nTotalRadii) {
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(DirectionsG1[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        directionMeanG1<-convertVec2unitVec(colMeans(t(DirectionsG1[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        directionMeanG1[i,]<-frechetMean(DirectionsG1[i,,])
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(DirectionsG1[i,,])
        directionMeanG1[i,]<-pns(DirectionsG1[i,,],sphere.type = sphereType)$PNS$mean
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
    directionMeanG2<-array(NA,dim = c(nTotalRadii,3))               
    for (i in 1:nTotalRadii) {
      # For extremely concentrated data we use Mardia mean direction 
      pcaTemp<-prcomp(t(DirectionsG2[i,,]))
      if(pcaTemp$sdev[1]<1e-02 | pcaTemp$sdev[2]<1e-02){
        directionMeanG2<-convertVec2unitVec(colMeans(t(DirectionsG2[i,,])))
      }else if(typeOfMeanDirection=="Frechet"){
        directionMeanG2[i,]<-frechetMean(DirectionsG2[i,,])
      }else if(typeOfMeanDirection=="PNS"){
        sphereType<-kurtosisTestFunction(DirectionsG2[i,,])
        directionMeanG2[i,]<-pns(DirectionsG2[i,,],sphere.type = sphereType)$PNS$mean
      }else{
        stop("Please specify the typeOfMeanDirection by PNS or Frechet")
      }
    }
    
    #geometric mean of the spokes' lengths
    meanLogR_G1 <-rowMeans(log(radiiG1))
    meanRadii_G1 <- exp( meanLogR_G1 )  
    meanLogR_G2 <-rowMeans(log(radiiG2))
    meanRadii_G2 <- exp( meanLogR_G2 )  
    
    
    #plot mean GP-ds-reps
    skelMeanG12<-c()
    skelMeanG12<-rbind(skelMeanG1[1:downSpoeksNumber,],
                       skelMeanG1[1:downSpoeksNumber,],skelMeanG1[(downSpoeksNumber+1):skelPointNo,])
    skelMeanG22<-c()
    skelMeanG22<-rbind(skelMeanG2[1:downSpoeksNumber,],
                       skelMeanG2[1:downSpoeksNumber,],skelMeanG2[(downSpoeksNumber+1):skelPointNo,])
    tipSpokesMeanG1<-skelMeanG12+directionMeanG1*meanRadii_G1
    tipSpokesMeanG2<-skelMeanG22+directionMeanG2*meanRadii_G2
    
    # overlaid plot of CG and PD means
    if(plotting==TRUE){
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(rbind(skelMeanG12[i,],tipSpokesMeanG1[i,]),type="l",lwd = 1,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      # plot3d(skelMeanG2,type="p",radius = 0.2,col = "red",expand = 10,box=FALSE,add = TRUE)
      for (i in 1:nTotalRadii) {
        plot3d(rbind(skelMeanG22[i,],tipSpokesMeanG2[i,]),type="l",lwd = 1,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Overlaid mean GP-ds-reps",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
    }
    
    ##############################
    # Hypothesis test on GP-ds-rep
    
    # hypothesis test on GP size of spokes' tips (i.e., centroid size)
    centroidSize_G1<-rep(NA,nSamplesG1)
    for (i in 1:nSamplesG1) {
      centroidSize_G1[i]<-centroid.size(spokesTips_G1[,,i])
    }
    centroidSize_G2<-rep(NA,nSamplesG1)
    for (i in 1:nSamplesG2) {
      centroidSize_G2[i]<-centroid.size(spokesTips_G2[,,i])
    }
    pValues_GP_sizes<-meanDifferenceTest1D(log(centroidSize_G1),log(centroidSize_G2),type = typeOfTest) 
    cat("pValue of GP sizes is:",pValues_GP_sizes,"\n")
    boxplot(centroidSize_G1, centroidSize_G2, names = c("G1","G2"),main="Centroid size")
    cat("sd GP size G1:",sd(centroidSize_G1),"sd GP size G2:",sd(centroidSize_G2),"\n")
    cat("Mean GP size G1:",mean(centroidSize_G1),"mean GP size G2:",mean(centroidSize_G2),"\n")
    
    
    # hypothesis test on spokes' lengths
    pValues_TtestRadii<-rep(NA,nTotalRadii)
    pb <- txtProgressBar(min = 0, max = nTotalRadii, style = 3) #progress bar
    for (i in 1:nTotalRadii) {
      setTxtProgressBar(pb, i) #create progress bar
      
      pValues_TtestRadii[i]<-meanDifferenceTest1D(log(radiiG1[i,]),
                                                  log(radiiG2[i,]),
                                                  type = typeOfTest)
    }
    
    # hypothesis test on  skeletal positions
    pValues_Skeletal<-rep(NA,skelPointNo)
    for (i in 1:skelPointNo) {
      pValues_Skeletal[i]<-meanDifferenceTestMultivariate(t(skelG1[i,,]),
                                                          t(skelG2[i,,]),
                                                          type = typeOfTest)
    }
    
    # hypothesis test on  directions
    pValues_Directions<-rep(NA,nTotalRadii)
    for (i in 1:nTotalRadii) {
      euclideanizedTemp<-euclideanization(DirectionsG1[i,,],
                                          DirectionsG2[i,,],
                                          type = typeOfStudy4directions)
      
      pValues_Directions[i]<-meanDifferenceTestMultivariate(euclideanizedTemp$euclideanG1,
                                                            euclideanizedTemp$euclideanG2,
                                                            type=typeOfTest)
    }
    
    #################################
    # plot significant GP-ds-rep GOPs
    
    # all GP-ds-rep p-values of spokes' positions, directions and radii
    pvalues_GP_ds_rep<-c(pValues_Skeletal,
                         pValues_Directions,
                         pValues_TtestRadii,
                         pValues_GP_sizes)
    
    alpha<-0.05
    significantPvalues<-which(pvalues_GP_ds_rep<=alpha)
    significantPvalues
    #adjust p-values
    FDR<-0.05
    pvalues_GP_ds_rep_BH<-p.adjust(pvalues_GP_ds_rep, method = "BH")
    pvalues_GP_ds_rep_Bonferroni<-p.adjust(pvalues_GP_ds_rep,method = "bonferroni")
    significantPvalues_BH<-which(pvalues_GP_ds_rep_BH<=FDR)
    significantPvalues_BH
    significantPvalues_Bonferroni<-which(pvalues_GP_ds_rep_Bonferroni<=FDR)
    significantPvalues_Bonferroni
    
    if(plotting==TRUE){
      # plot p-values and adjusted p-values
      df_GP <- data.frame(Type=c(rep("Raw p-value",length(pvalues_GP_ds_rep)),
                                 rep("Bonferroni",length(pvalues_GP_ds_rep)),
                                 rep("BH",length(pvalues_GP_ds_rep))),
                          ordereOfPvalues=1:length(pvalues_GP_ds_rep),
                          Values=c(sort(pvalues_GP_ds_rep),sort(pvalues_GP_ds_rep_Bonferroni),sort(pvalues_GP_ds_rep_BH)))
      p<-ggplot(df_GP, aes(x=ordereOfPvalues, y=Values, group=Type))
      p + geom_line(aes(linetype=Type),size=1)+
        # geom_line(aes(linetype=Type, col=Type),size=1)+
        # geom_point(aes(shape=Type),alpha=0.7)+
        geom_hline(yintercept=0.05,linetype="solid", color = "red")+
        scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
        theme_bw()+
        theme(plot.title = element_text(size = 17, hjust = 0.5),
              legend.text=element_text(size=17),
              legend.title=element_blank(),
              # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              # legend.title=element_text(size=17),
              legend.position="bottom",
              axis.title=element_text(size=17),
              legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
        guides(colour = guide_legend(title.hjust = 0.5))+
        xlab("Ranking of p-values") + ylab("p-values") 
    }
    
    n_s<-nTotalRadii
    n_f<-skelPointNo
    
    # plot significant Pos
    significantSkeletalPos<-significantPvalues[which(significantPvalues<=n_f)]
    significantSkeletalPos
    significantSkeletalPos_BH<-significantPvalues_BH[which(significantPvalues_BH<=n_f)]
    significantSkeletalPos_BH
    
    # plot significant directions
    significantDir<-significantPvalues[which(significantPvalues>n_f &
                                               significantPvalues<=(n_f+n_s))]-n_f
    significantDir
    significantDir_BH<-significantPvalues_BH[which(significantPvalues_BH>n_f &
                                                     significantPvalues_BH<=(n_f+n_s))]-n_f
    significantDir_BH
    
    # plot significant radii
    significantRadii<-significantPvalues[which(significantPvalues>(n_f+n_s) & significantPvalues<=(n_f+2*n_s))]-(n_f+n_s)
    significantRadii
    significantRadii_BH<-significantPvalues_BH[which(significantPvalues_BH>(n_f+n_s) & significantPvalues_BH<=(n_f+2*n_s))]-(n_f+n_s)
    significantRadii_BH
    
    #plot significant GOPs before and after the BH adjustment
    if(plotting==TRUE){
      
      #shape for plot
      srep1<-meanOfAll_G1_G2
      
      #without correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' lengths",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantRadii_BH) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' lengths after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #without correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantDir) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' directions",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      for (i in 1:nTotalRadii) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 0.2,col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in significantDir_BH) {
        plot3d(srep1[c(i,(i+nTotalRadii)),],type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
      }
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant spokes' directions after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #shape 4 plot 
      skelShape<-meanOfAll_G1_G2[skelRange,]
      
      #without correction
      open3d()
      plot3d(skelShape,type="s", radius = 0.4,col = "#90a0ff",expand = 10,box=FALSE,add = TRUE)
      plot3d(skelShape[significantSkeletalPos,],type="s", radius = 0.45,col = "red",expand = 10,box=FALSE,add = TRUE)
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant skeletal points",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
      #with correction
      open3d()
      plot3d(skelShape,type="s", radius = 0.4,col = "#90a0ff",expand = 10,box=FALSE,add = TRUE)
      plot3d(skelShape[significantSkeletalPos_BH,],type="s", radius = 0.45,col = "red",expand = 10,box=FALSE,add = TRUE)
      decorate3d(xlab = "x", ylab = "y", zlab = "z",
                 box = F, axes = TRUE, main = "Significant skeletal points after BH adjustment",
                 sub = NULL, top = T, aspect = FALSE, expand = 1.1)
      
    }
    
    
  }else if(typeOfAnalysis=="EDMA"){
    
    
    group1<-spokesTails_G1[skelRange,,]
    group2<-spokesTails_G2[skelRange,,]
    
    #plot
    if(plotting==TRUE){
      open3d()
      for (i in 1:nSamplesG1) {
        plot3d(group1[,,i],type="p",col = "blue",expand = 10,box=FALSE,add = TRUE)
      }
      for (i in 1:nSamplesG2) {
        plot3d(group2[,,i],type="p",col = "red",expand = 10,box=FALSE,add = TRUE)
      }
    }
    
    # calculate 3d tensors of distances (All EDMs)
    q<-dim(group1)[1]
    tensor1<-array(NA, dim=c(q,q,dim(group1)[3]))
    for (k in 1:dim(group1)[3]) {
      tensor1[,,k]<-as.matrix(dist(group1[,,k],diag = TRUE,upper = TRUE))
    }
    tensor2<-array(NA, dim=c(q,q,dim(group2)[3]))
    for (k in 1:dim(group2)[3]) {
      tensor2[,,k]<-as.matrix(dist(group2[,,k],diag = TRUE,upper = TRUE))
    }
    
    #EDM size 
    EDM_sizes_G1<-rep(NA,dim(group1)[3])
    for(k in 1:dim(group1)[3]){
      tempVec<-as.vector(tensor1[,,k][upper.tri(tensor1[,,k])])
      EDM_sizes_G1[k]<-exp(mean(log(tempVec)))
    }
    EDM_sizes_G2<-rep(NA,dim(group2)[3])
    for(k in 1:dim(group2)[3]){
      tempVec<-as.vector(tensor2[,,k][upper.tri(tensor2[,,k])])
      EDM_sizes_G2[k]<-exp(mean(log(tempVec)))
    }
    
    if(typeOfStudy=="shapeAnalysis"){
      for (k in 1:dim(group1)[3]) {
        tensor1[,,k]<-tensor1[,,k]/EDM_sizes_G1[k]
      }
      for (k in 1:dim(group2)[3]) {
        tensor2[,,k]<-tensor2[,,k]/EDM_sizes_G2[k]
      }
    }else if(typeOfStudy=="sizeAndShapeAnalysis"){
      tensor1<-tensor1
      tensor2<-tensor2
    }else{
      stop("Please choose the type of study as shapeAnalysis or sizeAndShapeAnalysis!")
    }
    
    # EDMA hypothesis testing
    
    pValues_EDM_sizes<-meanDifferenceTest1D(log(EDM_sizes_G1),log(EDM_sizes_G2),type = typeOfTest) 
    cat("pValue of EDM sizes is:",pValues_EDM_sizes,"\n")
    boxplot(EDM_sizes_G1, EDM_sizes_G2, names = c("G1","G2"),main="EDM size")
    cat("sd EDM size G1:",sd(EDM_sizes_G1),"sd EDM size G2:",sd(EDM_sizes_G2),"\n")
    cat("Mean EDM size G1:",mean(EDM_sizes_G1),"mean EDM size G2:",mean(EDM_sizes_G2),"\n")
    
    
    # p-value matrix
    pvalMatrix<-array(NA, dim=c(q,q))
    for (i in 1:q) {
      for (j in 1:q) {
        if(i!=j){
          pvalMatrix[i,j]<-t.test(tensor1[i,j,],tensor2[i,j,])$p.value
        }else{
          pvalMatrix[i,j]<-1
        }
      }
    }
    
    pvalMatrix

    #adjust p-values by Benjamini-Hochberg
    pvalues<-as.vector(pvalMatrix[upper.tri(pvalMatrix)])
    FDR<-0.05
    pvalues_BH<-p.adjust(pvalues,method = "BH")
    pvalues_Bonferroni<-p.adjust(pvalues,method = "bonferroni")
    significantPvalues_BH<-which(pvalues_BH<=FDR)
    significantPvalues_BH
    significantPvalues_Bonferroni<-which(pvalues_Bonferroni<=FDR)
    significantPvalues_Bonferroni
    
    
    # plot by ggplot
    if(plotting==TRUE){
      df_pval <- data.frame(Type=c(rep("Raw p-value",length(pvalues)),
                                   rep("Bonferroni",length(pvalues)),
                                   rep("BH",length(pvalues))),
                            ordereOfPvalues=1:length(pvalues),
                            Values=c(sort(pvalues),sort(pvalues_Bonferroni),sort(pvalues_BH)))
      p<-ggplot(df_pval, aes(x=ordereOfPvalues, y=Values, group=Type))
      p + geom_line(aes(linetype=Type),size=1)+
        # geom_line(aes(linetype=Type, col=Type),size=1)+
        # geom_point(aes(shape=Type),alpha=0.7)+
        geom_hline(yintercept=0.05,linetype="solid", color = "red")+
        scale_linetype_manual(values=c("solid","dotdash", "dotted")) +
        theme_bw()+
        theme(plot.title = element_text(size = 17, hjust = 0.5),
              legend.text=element_text(size=17),
              legend.title=element_blank(),
              # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              # legend.title=element_text(size=17),
              legend.position="bottom",
              axis.title=element_text(size=17),
              legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) +
        guides(colour = guide_legend(title.hjust = 0.5))+
        xlab("Ranking of p-values") + ylab("p-values") 
    }
    
  }else{
    stop('Please specify the type of analysis!')
  }
  
}
