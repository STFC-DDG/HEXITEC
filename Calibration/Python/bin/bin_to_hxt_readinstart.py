from read80x80Frame_readinstart import *
from charge_sharing_corrections import *
from hxtV3Write import hxtV3Write
from tqdm import tqdm
import numpy as np
from datetime import datetime
import cv2

def bin_to_hxt(filePath,gradFile=None,intFile=None,numX=80,numY=80,binWidth=10,numBins=800):
    
    #getting the file name for later saving of .hxt files
    filename = filePath.split('/')[-1].split('.')[0]
    directory = filePath.split('/')[:-1]
    
    frame = np.zeros((numX,numY))
    
    bins = np.arange(0,(numBins*binWidth),binWidth)
    
    #setting up empty aarrays and variables
    globalHist = np.zeros(numBins)
    globalHistCSD= np.zeros(numBins)
    globalHistCSA= np.zeros(numBins)
    pixelHist = np.zeros((numBins,80,80))
    pixelHistCSA = np.zeros((numBins,80,80))
    pixelHistCSD = np.zeros((numBins,80,80))
    
    nCSDTot=0
    nCSATot=0
   
    lowThresh = 50
    highThresh = 7000
    flag = 1
    
    #opening bin file
    fid = open(filePath,'rb')
    
    #calculating the number of frames in file
    end = fid.seek(0,2)
    nFrames = np.floor(end/(80*80)).astype(int)
    print(nFrames)
    
    #resetting position to beginning of file
    fid.seek(0,0)
    fNum = 0
    
    data = np.fromfile(fid,dtype=np.dtype(np.uint16),count=int((nFrames/10)*6400))
    
    for fNum in tqdm(range(np.ceil((nFrames/2)-1).astype(int))):
        #reading in frame
        frame = read80x80frame(fNum,data)
        noise_frame = ((frame>lowThresh)&(frame<highThresh))*frame
        #[XFrameRaw,YFrameRaw] = np.where(frame>lowThresh)
        RawNonZero = cv2.findNonZero(noise_frame)
        
        if RawNonZero is not None :
            [YFrameRaw,XFrameRaw] = [RawNonZero[:,0,i] for i in range(2)]
            bin_ind = bins.searchsorted(noise_frame[XFrameRaw,YFrameRaw])
            pixelHist[bin_ind,XFrameRaw,YFrameRaw] += 1
            globalHist[bin_ind] += 1
            
            #processing CSD and CSA
            [CSDFrame,nCSD] = CSD(noise_frame)
            [CSAFrame,nCSA] = CSA(noise_frame)
            
            nCSDTot += nCSD
            nCSATot += nCSA
            
            CSDNonZero = cv2.findNonZero(CSDFrame)    
            if CSDNonZero is not None:
                [YFrameCSD,XFrameCSD] = [CSDNonZero[:,0,i] for i in range(2)]
                bin_ind = bins.searchsorted(CSDFrame[XFrameCSD,YFrameCSD])
                pixelHistCSD[bin_ind,XFrameCSD,YFrameCSD] += 1
                globalHistCSD[bin_ind] += 1
            
            #[XFrameCSA,YFrameCSA] = np.where((CSAFrame>lowThresh) & (CSAFrame<highThresh))
            CSANonZero = cv2.findNonZero(CSAFrame) 
            if CSANonZero is not None:
                [YFrameCSA,XFrameCSA] = [CSANonZero[:,0,i] for i in range(2)]
                bin_ind = bins.searchsorted(CSAFrame[XFrameCSA,YFrameCSA])
                pixelHistCSA[bin_ind,XFrameCSA,YFrameCSA] += 1
                globalHistCSA[bin_ind] += 1
    fid.close()
    
    #writing spectra and bins to .hxt files
    hxtV3Write(pixelHist,bins,directory+'/'+filename+'_Raw.hxt')
    hxtV3Write(pixelHistCSD,bins,directory+'/'+filename+'_CSD.hxt')
    hxtV3Write(pixelHistCSA,bins,directory+'/'+filename+'_CSA.hxt')
        
        