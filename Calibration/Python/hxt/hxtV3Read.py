import numpy as np

def hxtV3Read(filePath):
    """
    """
    fid = open(filePath,'rb')
    
    #reading first 8 charachters to distinguish file type
    label = fid.read(8).decode('ISO-8859-1')
    
    if label == 'HEXITECH':
        version = np.fromfile(fid,dtype='uint64',count=1)
        if version == 3:
            mssX = np.fromfile(fid,dtype='uint32',count=1)[0]
            mssY = np.fromfile(fid,dtype='uint32',count=1)[0]
            mssZ = np.fromfile(fid,dtype='uint32',count=1)[0]
            mssRot = np.fromfile(fid,dtype='uint32',count=1)[0]
            GalX = np.fromfile(fid,dtype='uint32',count=1)[0]
            GalY = np.fromfile(fid,dtype='uint32',count=1)[0]
            GalZ = np.fromfile(fid,dtype='uint32',count=1)[0]
            GalRot = np.fromfile(fid,dtype='uint32',count=1)[0]
            GalRot2 = np.fromfile(fid,dtype='uint32',count=1)[0]
            nCharFPreFix = np.fromfile(fid,dtype='int32',count=1)[0]
            filePrefix = fid.read(nCharFPreFix).decode('ISO-8859-1')
            dummy = fid.read(100-nCharFPreFix).decode('ISO-8859-1')
            timestamp = fid.read(16).decode('ISO-8859-1')
            nRows = np.fromfile(fid,dtype='uint32',count=1)[0]
            nCols = np.fromfile(fid,dtype='uint32',count=1)[0]
            nBins = np.fromfile(fid,dtype='uint32',count=1)[0]
            bins = np.fromfile(fid,dtype='double', count=nBins)
            d = np.fromfile(fid,dtype=np.float64,count=nBins*nRows*nCols)
            M = np.reshape(d,((nCols,nRows,nBins)))

    else:
        print('Not Version 3 of HXT File- Zeros Returned')
        M = 0
        bins = 0
    
    fid.close()
    
    return M, bins

    
