# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:08:54 2020

@author: rhian
"""
def read80x80frame(binaryFile):
    import numpy as np
    frame = np.zeros([80,80])
    n=0
    dt = np.dtype(np.uint16)     
    dataVector = np.fromfile(binaryFile, dtype=dt, count=6400)
    yaxis = np.arange(0,80,1)
    block = np.arange(0,20,1)
    nBlocks = np.arange(1,5,1)    
    for row in yaxis:
        for j in block:
            for k in nBlocks:
                col = j+((k-1)*20)
                pixVal = dataVector[n]
                frame[row,col] = pixVal
                n=n+1   
                
    return frame

#%%
import numpy as np, matplotlib.pyplot as plt, os    

filePath = r'C:\Users\dml36958\Desktop\DesktopDocs\BeamTime\Calibration\Am241_28C_600V_noCSD_50ADUThresh_300s_190925_155355.bin'
binaryFile = open(filePath)
binaryFile.seek(0)
nFrames = (os.stat(filePath).st_size/(80*80))/2

threshold=80
for slice in np.arange(0,nFrames,1):
    frame=read80x80frame(binaryFile)
    frame[frame<threshold] = 0
    
    plt.figure(1)
    plt.imshow(frame)
    plt.pause(1)
