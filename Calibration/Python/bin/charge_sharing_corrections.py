import numpy as np
from datetime import datetime

def CSA(frameIn,neigh_range=1):
    """
    Function to conduct charge sharing correction on frame of data
    """

    xPix = frameIn.shape[0]
    yPix = frameIn.shape[1]
    
    #defining variables
    eventThresh = 50
    binMax = 8000
    nCSA = 0 
    
    #setting up raw frame array and inputting frame data into array
    frameRaw = np.zeros((xPix+4,yPix+4))
    frameRaw[2:2+xPix, 2:2+yPix] = frameIn
    #setting up empty array for CSA frame
    CSAframe = np.zeros((xPix,yPix))
    
    [x,y] = np.where(frameRaw>eventThresh)
    
    for (n,j) in zip(x,y):
        if frameRaw[n,j] > eventThresh:
            neighbours = frameRaw[n-neigh_range:n+1+neigh_range,\
                                           j-neigh_range:j+1+neigh_range][frameRaw[n-neigh_range:n+1+neigh_range,\
                                           j-neigh_range:j+1+neigh_range]>eventThresh]
            #will return array of length one if only hit in central event pixel detected
            if len(neighbours) > 1:
                #summing total energy across charge sharing event
                FullE = (frameRaw[n-neigh_range:n+1+neigh_range,\
                                         j-neigh_range:j+1+neigh_range]).sum()
                #if FullE is less than Emax then place all event energy in central pixel
                if FullE < binMax:
                    [maxX,maxY] = np.unravel_index(frameRaw[n-neigh_range:n+1+neigh_range,\
                                                    j-neigh_range:j+1+neigh_range].argmax(),frameRaw[n-neigh_range:n+1+neigh_range,\
                                                    j-neigh_range:j+1+neigh_range].shape)
                    frameRaw[n-neigh_range:n+1+neigh_range,j-neigh_range:j+1+neigh_range] = 0
                    frameRaw[n-neigh_range+maxX,j-neigh_range+maxY] = FullE
                    nCSA += 1
    
    CSAframe = frameRaw[2:2+xPix,2:2+yPix]

    return CSAframe, nCSA

def CSD(frameIn,neigh_range=1):
    """
    Function to conduct charge sharing correction of frame of data
    """
    
    xPix = frameIn.shape[0]
    yPix = frameIn.shape[1]
    
    #defining variables
    eventThresh = 50
    binMax = 8000
    nCSD = 0 
    
    #setting up raw frame array and inputting frame data into array
    frameRaw = np.zeros((xPix+4,yPix+4))
    frameRaw[2:2+xPix, 2:2+yPix] = frameIn
    #setting up empty array for CSA frame
    CSDframe = np.zeros((xPix,yPix))
      
    [x,y] = np.where(frameRaw>eventThresh)
    
    for (n,j) in zip(x,y):
        if frameRaw[n,j] > eventThresh:
# =============================================================================
#             neighbours = np.where(frameRaw[n-neigh_range:n+1+neigh_range,\
#                                            j-neigh_range:j+1+neigh_range]>eventThresh)[0]
# =============================================================================
            neighbours = frameRaw[n-neigh_range:n+1+neigh_range,\
                                           j-neigh_range:j+1+neigh_range][frameRaw[n-neigh_range:n+1+neigh_range,\
                                           j-neigh_range:j+1+neigh_range]>eventThresh]
            #will return array of length one if only hit in central event pixel detected
            if len(neighbours) > 1:
                #summing total energy across charge sharing event
                FullE = (frameRaw[n-neigh_range:n+1+neigh_range,\
                                         j-neigh_range:j+1+neigh_range]).sum()
                #if FullE is less than Emax then remove all data from charge-sharing event
                if FullE < binMax:
                    frameRaw[n-neigh_range:n+1+neigh_range,j-neigh_range:j+1+neigh_range] = 0
                    nCSD += 1

    CSDframe = frameRaw[2:2+xPix,2:2+yPix]

    return CSDframe, nCSD

