import numpy as np
from datetime import datetime
from tqdm import tqdm

def read80x80frame(frameNum,data,flag=1):
    """
    Function to return frames of 80x80 together in 80x80 frame of data.
    Used for analysis of HEXITEC .bin file data.
    """
    
    frame = np.zeros((80,80))
        
    if flag==1:
       #Sorting data into 240x80 pixels with correct orientation
       dataIn = data[frameNum*6400:((frameNum+1)*6400)]
       #n = 0
       row=[]
       [row.extend([i]*80) for i in range(80)]
       col = []
       [col.extend([0+i,20+i,40+i,60+i]) for i in range(20)]
       col = col*80
       
       frame[row,col] = dataIn
       
    elif flag == 0:
        #Read data into 80x80 pxiels with interleaved format for speed
        frame = dataIn
    else:
        print('Invalid Flag. Use Flag=1 for sorting frame into correct orientation or Flag=0 for interleaved frame which will be faster.')
    return frame


def threading_read_80x80_frame(coreID,numCores,fid,end,flag=1):
    
    fNums = np.arange(coreID*6400,end,numCores*6400)
    
    for i in tqdm(range(len(fNums))):
        read80x80frame(fNums[i], fid, end)
        
        