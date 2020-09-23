# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:36:57 2020

@author: dml36958
"""
def read80x80frame(binaryFile):
    import numpy as np
    frame = np.zeros([80,80])
    n=0
#    binaryFile.seek(0,0)
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
from tkinter.filedialog import askopenfilename
from scipy.ndimage import measurements

import numpy as np, matplotlib.pyplot as plt, os, docx
import cProfile, pstats, io

pr = cProfile.Profile()
pr.enable()

threshold = 80
bins = np.arange(1,8000,1)
s = np.ones((3,3))

edgePxX = 80
edgePxY = 80

pixelHistRaw = np.zeros((edgePxX,edgePxY,len(bins)))
pixelHistCSD = np.zeros((edgePxX,edgePxY,len(bins)))
pixelHistCSA = np.zeros((edgePxX,edgePxY,len(bins)))

globalHistRaw = np.zeros((len(bins)))
globalHistCSD = np.zeros((len(bins)))
globalHistCSA = np.zeros((len(bins)))

filePath = askopenfilename()
print(filePath)

binaryFile = open(filePath)
binaryFile.seek(0)
nFrames = 1000#(os.stat(filePath).st_size/(80*80))/2
for slice in np.arange(0,nFrames,1):
    frame=read80x80frame(binaryFile)
    frame[frame<threshold] = 0
    FrameCSD = np.copy(frame)
    FrameCSA = np.copy(frame)
    frameLabelled, numEvents = measurements.label(frame,structure=s)
    frameLabelled = frameLabelled
    for hitNum in range(1,numEvents):
        Xhits, Yhits = np.where(frameLabelled==hitNum)
        if len(Xhits)>1:
            FrameCSD[Xhits,Yhits]=0
            FrameCSA[Xhits,Yhits]=0
            totalHit = sum(frame[Xhits,Yhits])
            maxHit = np.where(frame[Xhits,Yhits]==max(frame[Xhits,Yhits]))
            FrameCSA[Xhits[maxHit[0][0]],Yhits[maxHit[0][0]]] = totalHit
    
    RawBins = np.digitize(frame, bins)
    [x,y] = np.where(RawBins>0)
    for z in range(0, len(x)):
        if RawBins[x[z],y[z]]>7998:
            continue
        pixelHistRaw[x[z],y[z],RawBins[x[z],y[z]]] = pixelHistRaw[x[z],y[z],RawBins[x[z],y[z]]] + 1
        
    CSDBins = np.digitize(FrameCSD, bins)
    [x,y] = np.where(CSDBins>0)
    for z in range(0, len(x)):
        if CSDBins[x[z],y[z]]>7998:
            continue
        pixelHistCSD[x[z],y[z],CSDBins[x[z],y[z]]] = pixelHistCSD[x[z],y[z],CSDBins[x[z],y[z]]] + 1
        
    CSABins = np.digitize(FrameCSA, bins)
    [x,y] = np.where(CSABins>0)
    for z in range(0, len(x)):
        if CSABins[x[z],y[z]]>7998:
            continue        
        pixelHistCSA[x[z],y[z],CSABins[x[z],y[z]]] = pixelHistCSA[x[z],y[z],CSABins[x[z],y[z]]] + 1
        
        
    RawBins = RawBins[RawBins>0]
    for z in range(0, len(RawBins)):
        if RawBins[z]>7998:
            continue        
        globalHistRaw[RawBins[z]] = globalHistRaw[RawBins[z]] + 1
        
    CSDBins = CSDBins[CSDBins>0]
    for z in range(0, len(CSDBins)):
        if CSDBins[z]>7998:
            continue               
        globalHistCSD[CSDBins[z]] = globalHistCSD[CSDBins[z]] + 1
        
    CSABins = CSABins[CSABins>0]
    for z in range(0, len(CSABins)):
        if CSABins[z]>7998:
            continue            
        globalHistCSA[CSABins[z]] = globalHistCSA[CSABins[z]] + 1
    
#plt.figure()
#plt.plot(bins,globalHistRaw)    
#plt.plot(bins,globalHistCSD) 
#plt.plot(bins,globalHistCSA)    

pr.disable()   
s = io.StringIO()
sortby = SortKey.CUMULATIVE
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())
##%% Calibration on CSD
#energyPeaks = [59.54, 13.94, 17.75]
#searchLow = [1500, 300, 400]
#searchHigh = [2500, 400, 650]
#
#peakAmplitudes = np.zeros((edgePxX, edgePxY, len(energyPeaks)))
#peakCentroids = np.zeros((edgePxX, edgePxY, len(energyPeaks)))
#peakFWHM = np.zeros((edgePxX, edgePxY))
#
#m = np.zeros((edgePxX, edgePxY))
#c = np.zeros((edgePxX, edgePxY))
#
#searchLowLoc = np.zeros((len(searchLow)))
#searchHighLoc = np.zeros((len(searchHigh)))
#
#for energy in range(0, len(energyPeaks)):
#    searchLowLoc[energy] = np.where(bins==searchLow[energy])[0]
#    searchHighLoc[energy] = np.where(bins==searchHigh[energy])[0]
#
#searchLowLoc=searchLowLoc.astype(int)
#searchHighLoc=searchHighLoc.astype(int)
#
#x=0
#y=0
#for z in range(0, (edgePxX*edgePxY)-1):
#    for energy in range(0, len(energyPeaks)):
#        currentPixel = pixelHistCSD[x,y,searchLowLoc[energy]:searchHighLoc[energy]]
#        peakAmplitudes[x,y,energy] = np.max(currentPixel)
#        peakCentroids[x,y,energy] = np.where(currentPixel == peakAmplitudes[x,y,energy])[0][0] + searchLowLoc[energy]
#        if energy == 0:
#            FWHMLow = np.where(currentPixel>(peakAmplitudes[x,y,energy])/2)[0][0]
#            FWHMHigh = np.where(currentPixel>(peakAmplitudes[x,y,energy])/2)[0][-1]
#            peakFWHM[x,y] = FWHMHigh-FWHMLow
#    
#    fit = np.polyfit(peakCentroids[x,y,:],energyPeaks,1)
#    m[x,y] = fit[0]
#    c[x,y] = fit[1]
#    x=x+1
#    if x==edgePxX:
#        x=0
#        y=y+1
#    
##%% Creating Figures
##figurePath = os.getcwd() + '\Figures'
##if os.path.exists(figurePath) == False:
##    os.mkdir(figurePath)
##
##os.chdir (figurePath)
#
#cmap=plt.get_cmap('jet') 
#plt.figure()
#plt.imshow(m, vmin=np.median(m)*0.5, vmax=np.median(m)*1.5, cmap=cmap)
#cbar = plt.colorbar()
#cbar.set_label('Gradient (keV/ADU)',fontsize=18)
#plt.xlabel('X (Pix)',fontsize=16)  
#plt.ylabel('Y (Pix)',fontsize=16)   
##plt.savefig('gradients.png')
#plt.close()
#
#cmap=plt.get_cmap('jet') 
#plt.figure()
#plt.imshow(c, vmin=np.median(c)*0.5, vmax=np.median(c)*1.5, cmap=cmap)
#cbar = plt.colorbar()
#cbar.set_label('Intercept (keV)',fontsize=18)
#plt.xlabel('X (Pix)',fontsize=16)  
#plt.ylabel('Y (Pix)',fontsize=16)   
##plt.savefig('intercepts.png')
#plt.close()
#
#cmap=plt.get_cmap('jet') 
#plt.figure()
#plt.imshow(peakFWHM, vmin=np.median(peakFWHM)*0.5, vmax=np.median(peakFWHM)*1.5, cmap=cmap)
#cbar = plt.colorbar()
#cbar.set_label('FWHM about 60keV Peak (keV)',fontsize=18)
#plt.xlabel('X (Pix)',fontsize=16)  
#plt.ylabel('Y (Pix)',fontsize=16)   
##plt.savefig('FWHM.png')
#plt.close()
#
#
#
#
#       
#    
#    
##%% Publishing Results
#
#SensorName = 'D212___Blah'
#doc = docx.Document()
#doc.add_heading('HEXITEC: Detector Test Results', 0)
#doc.add_heading(SensorName, 1)
#doc.save('test.docx')
#    
##    hitsMask = np.zeros((edgePxX,edgePxY))
##    hitsMask[frameLabelled>1] = 1
##    Xhits, Yhits = np.where(hitsMask==1)
#    