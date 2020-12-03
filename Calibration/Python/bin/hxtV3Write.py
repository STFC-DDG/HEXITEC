import numpy as np
import tkinter as tk
from tkinter import filedialog
from datetime import datetime
import os

def hxtV3Write(M, commonX=None, filePath=None):
    
    nRows = M.shape[0]
    nCols = M.shape[1]
    nBins = M.shape[2]
    
    
    hxtVersion = 3
    hxtLabel = 'HEXITECH'
    dummyMotor = 2222222
    
    #generating x axis information if not provided within input
    if not commonX.shape:
        commonX = np.arange(1,nBins)
    
    #generating file path if not provided within input
    if filePath == None:
        root = tk.Tk()
        #ensures file dialog windows open up on top
        root.attributes("-topmost", True)
        root.withdraw()
        #obtaining directory and filepath
        filePath = filedialog.askopenfilenames(title='Save Data as',filetypes=[("hxt files",'.hxt')],parent=root)
        PathName = '/'.join(filePath[0].split('/')[:-1])
        FileName = filePath[0].split('/')[-1].split('.')[0]
    else:
        FileName = 'Dummy'
    
    #getting time to calculate duration of function
    time_init = datetime.now()
    
    #generating dummy PreFix based on length of FileName
    nCharFPrefix = len(FileName)
    dummyPreFix = 'x' * (100-nCharFPrefix)
        
    dummyTimeStamp = '00000000_0000000'
    
    #opening file
    fid = open(filePath,'ab')
        
    #writing data packet
    fid.write(hxtLabel.encode('utf-8'))
    #uint64 is an 8 bit integer
    fid.write(hxtVersion.to_bytes(8,'little',signed=True))
    
    #dummy motor positions and time stamps
    #uint32 is a 4 bit integer
    mssX = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    mssY = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    mssZ = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    mssRot = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    GalX = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    GalY = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    GalZ = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    GalRot = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    GalRot2 = fid.write(dummyMotor.to_bytes(4,'little',signed=False))
    nCharFPrefix = fid.write(nCharFPrefix.to_bytes(4,'little',signed=True))
    filePreFix = fid.write(FileName.encode('utf-8'))
    dummy = fid.write(dummyPreFix.encode('utf-8'))
    timeStamp = fid.write(dummyTimeStamp.encode('utf-8'))
    
    #writing number of rows,bins,columns
    fid.write(nRows.to_bytes(4,'little', signed=False))
    fid.write(nCols.to_bytes(4,'little', signed=False))
    fid.write(nBins.to_bytes(4,'little', signed=False))
    
    #write X axis
    fid.write(commonX.astype('double'))
    
    #write data
    fid.write(M.astype('double'))
    
    fid.close
    time_final = datetime.now()
    print(('.hxt file produced in {0}s').format((time_final-time_init).total_seconds()))