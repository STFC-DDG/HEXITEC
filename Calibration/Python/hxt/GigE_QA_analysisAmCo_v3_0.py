import numpy as np
import matplotlib.pyplot as plt
from hxtV3Read import hxtV3Read
import tkinter as tk
from tkinter import filedialog
import os
import time
import pandas as pd


#setting directory into which results are saved; not neccessarily
#same directory as where data is saved
data_dir = 'C:/Users/uhf41485/Documents/HEXITEC_analysis'

#getting screen width and height for images
root = tk.Tk()
width = root.winfo_screenwidth()/100
height = root.winfo_screenheight()/100
root.withdraw()

#Assuming both Co-57 and Am241 files are in the same folder
#selecting Co-57 file for analysis
root = tk.Tk()
#esnures file dialog windows open up on top
root.attributes("-topmost", True)
root.withdraw()
#obtaining Co57 and Am241 file names
CoFileName = filedialog.askopenfilenames(title='Select Co57 .hxt file',filetypes=[("hxt files",'.hxt')],parent=root)
FileName = filedialog.askopenfilenames(title='Select Am241 .hxt file',filetypes=[("hxt files",'.hxt')],parent=root)
#obtaining folder name and experiment name
folderLoc = '/'.join(FileName[0].split('/')[:-1])
CoExpName = CoFileName[0].split('/')[-1].split('.')[0]
AmExpName = (FileName[0].split('/')[-1]).split('.')[0]

#making directory in which to put results based on AmExpName
results_dir = data_dir + '/' + AmExpName + '_' + time.strftime("%Y%m%d-%H%M%S")
os.mkdir(results_dir)

#creating empty pd dataframe detailing issues with non-valid pixels
pixel_array = np.unravel_index(np.arange(6400),(80,80))
non_valid_dataframe = pd.DataFrame(index=pd.MultiIndex.from_arrays(pixel_array,names=('x','y')),columns=('60keV Height','60keV Position','60keV FWHM'))

# =============================================================================
# folderLoc = 'D:\Redlen_CZT'
# CoFileName = folderLoc + '/BenTestLG_RedlenCZT_ASIC45_Co57_201023_115519.hxt'
# FileName = folderLoc + '/BenTestLG_RedlenCZT_ASIC45_Am241_201023_143415.hxt'
# =============================================================================

#Detector parameters
nCols = 80
nRows = 80

#Gain parameters 0=High,1=Low. This effects where peaks are searched for
gain = 1

#Gain ratio set by gain
if gain == 0:
    gain_ratio = 1
elif gain == 1:
    gain_ratio = 1/4

#loading in data from .hxt file
#Am241 files
for i in range(len(FileName)):
    if i == 0:
        specPix, bins = hxtV3Read(FileName[i])
    else:
        specPix += hxtV3Read(FileName[i])[0]
        binsNew = hxtV3Read(FileName[i])[1]
        if not (bins==binsNew).all():
            print('Bins in individual Am241 files are not the same. Rebin data or select different files')
            break
#Co57 files
for i in range(len(CoFileName)):
    if i == 0:
        specPixCo, binsCo = hxtV3Read(CoFileName[i])
    else:
        specPixCo += hxtV3Read(CoFileName[i])[0]
        binsNew = hxtV3Read(CoFileName[i])[1]
        if not (binsCo==binsNew).all():
            print('Bins in individual Co57 files are not the same. Rebin data or select different files')
            break

#loading in bins from .hxt file. Removes requirement to change values when using
#different bin setup
binWidth = (bins[-1]-bins[-2]).astype(int)
nBins = len(bins)
if not (bins == binsCo).all():
    print('Bin values of Co and Am files are different. Current script does not support this functionality')

#Am241 info
AmADU = (np.ceil(2000/binWidth)*gain_ratio).astype(int)
AmSearchADU = (np.ceil(400/binWidth)*gain_ratio).astype(int)
AmSearchADU2 = (np.ceil(100/binWidth)*gain_ratio).astype(int)
AmKeV = [59.54, 17.75, 13.94, 122]

heightThreshold = 10
PeakMin = AmADU - AmSearchADU
PeakMax = AmADU + AmSearchADU

#setting up empty arrays
#peak position of Am and Co peaks
PeakPos = np.zeros((80,80))
PeakPos2 = np.zeros((80,80))
PeakPos3 = np.zeros((80,80))
PeakPos4 = np.zeros((80,80))

#peak height of Co peak
peakMax4 = np.zeros((80,80))

valid = np.ones((80,80))
validFWHM = np.zeros((80,80))
threshold = np.ones((80,80))*50
gradient = np.zeros((80,80))
intercept = np.zeros((80,80))
FWHM_ADU = np.zeros((80,80))
FWHM_keV = np.zeros((80,80))
LHWHM = np.zeros((80,80))
UHWHM = np.zeros((80,80))
FWHM_ADU_Co = np.zeros((80,80))
FWHM_keV_Co = np.zeros((80,80))
LHWHM_Co = np.zeros((80,80))
UHWHM_Co = np.zeros((80,80))

#find the Am241 60keV peak heights
PeakHeights = np.amax(specPix[:,:,AmADU-AmSearchADU:AmADU+AmSearchADU],axis=2)


for iRow in range(nRows):
    for iCol in range(nCols):
        #if peak height is too low then set calibration to zero and threshold high
        if PeakHeights[iRow,iCol] < heightThreshold:
            valid[iRow,iCol] = 0
            threshold[iRow,iCol] = 2000
            #other values already zero
            #adding peak_height value to dataframe
            non_valid_dataframe['60keV Height'][iRow,iCol] = PeakHeights[iRow,iCol]
            print('Peak less than threshold {0},{1}'.format(iRow,iCol))
        else:
            #find peak position
            PeakPos[iRow,iCol] = (binWidth*(PeakMin + np.where(specPix[iRow,iCol,AmADU-AmSearchADU:AmADU+AmSearchADU] == PeakHeights[iRow,iCol])[0][0])).astype(int)
            #if peak position is out of range set to zero then invalid
            if (PeakPos[iRow,iCol] < (binWidth*PeakMin) + 1) or (PeakPos[iRow,iCol] > (binWidth*PeakMax)-1):
                valid[iRow,iCol] = 0
                threshold[iRow,iCol] = 0
                #other values already zero
                #adding peak_pos value to datframe
                non_valid_dataframe['60keV Position'][iRow,iCol]=PeakPos[iRow,iCol]
            else:
                #go on to complete calibration and analysis
                roughCal = PeakPos[iRow,iCol]/AmKeV[0]
                roughPos2 = np.ceil((roughCal*AmKeV[1])-(10*binWidth)).astype(int)
                peakMax2 = np.amax(specPix[iRow,iCol,(int(np.floor((roughPos2/binWidth)-AmSearchADU2))):int(np.ceil(((roughPos2/binWidth)+AmSearchADU2)))])
                roughPos3 = np.ceil((roughCal*AmKeV[2])-(10*binWidth)).astype(int)
                peakMax3 = np.amax(specPix[iRow,iCol,(int(np.floor((roughPos3/binWidth)-AmSearchADU2))):int(np.ceil(((roughPos3/binWidth)+AmSearchADU2)))])
                roughPos4 = np.ceil((roughCal*AmKeV[3])/binWidth).astype(int)
                peakMax4[iRow,iCol] = np.amax(specPixCo[iRow,iCol,roughPos4-AmSearchADU:roughPos4+AmSearchADU])
                
                PeakPos2[iRow,iCol] = (np.floor((roughPos2/binWidth))- AmSearchADU2 + np.where(specPix[iRow,iCol,(int(np.floor((roughPos2/binWidth)))-AmSearchADU2):(int(np.ceil((roughPos2/binWidth)))+AmSearchADU2)] == peakMax2)[0][0]).astype(int)
                PeakPos3[iRow,iCol] = (np.floor((roughPos3/binWidth))- AmSearchADU2 + np.where(specPix[iRow,iCol,(int(np.floor((roughPos3/binWidth)))-AmSearchADU2):(int(np.ceil((roughPos3/binWidth)))+AmSearchADU2)] == peakMax3)[0][0]).astype(int)
                PeakPos4[iRow,iCol] = (roughPos4 - AmSearchADU + np.where(specPixCo[iRow,iCol,(roughPos4-AmSearchADU):(roughPos4+AmSearchADU)] == peakMax4[iRow,iCol])[0][0]).astype(int)
            
                #set threshold
                #find point where level is 30% of peak height
                if specPix[iRow,iCol,11] > PeakHeights[iRow,iCol]*0.3:
                    step = 12
                    while specPix[iRow,iCol,step] > PeakHeights[iRow,iCol]*0.3:
                        step = step + 1
                    threshold[iRow,iCol] = step*binWidth
                
# =============================================================================
#                 if iCol == 5:
#                     fig101 = plt.figure()
#                     plt.plot(bins,np.squeeze(specPix[iRow,iCol,:]))
#                     plt.plot(PeakPos2[iRow,iCol]*binWidth,peakMax2,color='red',marker='o')
#                     plt.plot(PeakPos3[iRow,iCol]*binWidth,peakMax3,color='green',marker='o')
#                     plt.bar(threshold[iRow,iCol],PeakHeights[iRow,iCol])
#                     plt.title('('+str(iRow)+','+ str(iCol)+')')
#                     plt.show()
# =============================================================================
                        
                if (threshold[iRow,iCol] < 500) and (specPix[iRow,iCol,int((PeakPos[iRow,iCol]+500)/binWidth)] < 10):
                    ADU = [PeakPos3[iRow,iCol]*binWidth, PeakPos2[iRow,iCol]*binWidth, PeakPos[iRow,iCol], PeakPos4[iRow,iCol]*binWidth]
                    E = sorted(AmKeV)
                    coeffs = np.polyfit(ADU,E,1)
                    gradient[iRow,iCol] = coeffs[0]
                    intercept[iRow,iCol] = coeffs[1]
                    
                    #interpolate spectrum to get resolution
                    interpSpec = np.interp(np.transpose(np.squeeze(np.arange(1,nBins*binWidth))),np.transpose(bins),np.squeeze(specPix[iRow,iCol,:]))
                    #measure the FWHM
                    try:
                        LHWHM[iRow,iCol] = PeakPos[iRow,iCol] - (AmSearchADU*binWidth) + np.where(interpSpec[int(PeakPos[iRow,iCol]-(AmSearchADU*binWidth)):int(PeakPos[iRow,iCol])] > PeakHeights[iRow,iCol]/2)[0][0]
                    except:
                        LHWHM[iRow,iCol] = PeakPos[iRow,iCol] + np.where(interpSpec[:int(PeakPos[iRow,iCol])] < PeakHeights[iRow,iCol]/2)[0][-1]
                    try:
                        UHWHM[iRow,iCol] = PeakPos[iRow,iCol] + np.where(interpSpec[int(PeakPos[iRow,iCol]):int(PeakPos[iRow,iCol] + (AmSearchADU*binWidth))] < PeakHeights[iRow,iCol]/2)[0][0]
                    except:
                        UHWHM[iRow,iCol] = PeakPos[iRow,iCol] + np.where(interpSpec[int(PeakPos[iRow,iCol]):] < PeakHeights[iRow,iCol]/2)[0][0]
                    FWHM_ADU[iRow,iCol] = UHWHM[iRow,iCol] - LHWHM[iRow,iCol]
                    FWHM_keV[iRow,iCol] = FWHM_ADU[iRow,iCol]*gradient[iRow,iCol]
                    
                    if peakMax4[iRow,iCol] > 30:
                        #interpolate spectrum to get resolution of Co peak
                        interpSpec_Co = np.interp(np.transpose(np.squeeze(np.arange(1,nBins*binWidth))),np.transpose(bins),np.squeeze(specPixCo[iRow,iCol,:]))
                        #measure the FWHM
                        LHWHM_Co[iRow,iCol] = (PeakPos4[iRow,iCol]-(AmSearchADU))*binWidth + np.where(interpSpec_Co[int((PeakPos4[iRow,iCol]-AmSearchADU))*binWidth:int(PeakPos4[iRow,iCol]*binWidth)] > peakMax4[iRow,iCol]/2)[0][0]
                        UHWHM_Co[iRow,iCol] = (PeakPos4[iRow,iCol]*binWidth) + np.where(interpSpec_Co[int((PeakPos4[iRow,iCol]*binWidth)):int((PeakPos4[iRow,iCol]+AmSearchADU)*binWidth)] < peakMax4[iRow,iCol]/2)[0][0]
                        FWHM_ADU_Co[iRow,iCol] = UHWHM[iRow,iCol]-LHWHM[iRow,iCol]
                        FWHM_keV_Co[iRow,iCol] = FWHM_ADU_Co[iRow,iCol]*gradient[iRow,iCol]
# =============================================================================
#                         if (iRow==30) and (iCol==30):
#                             fig20 = plt.figure()
#                             plt.plot(interpSpec_Co)
#                             plt.show()
# =============================================================================
                    
                    #make valid FWHM if <2keV
                    if FWHM_keV[iRow,iCol] < (2.5/gain_ratio):
                        validFWHM[iRow,iCol] = 1
                    #adding non_valid FWHM to dataframe
                    else:
                        non_valid_dataframe['60keV FWHM'][iRow,iCol] = FWHM_keV[iRow,iCol]
                        

#removing good pixel entries from non_valid datframe 
non_valid_dataframe.dropna('index',how='all',inplace=True)
#saving dataframe to file
non_valid_dataframe.to_csv(results_dir+'/invalid_pixel_info.csv',sep=',')

#adding ability to scroll through spectra of non_valid pixels
#titles are based on what makes them invalid
nvalid_x = non_valid_dataframe.index.get_level_values('x')
nvalid_y = non_valid_dataframe.index.get_level_values('y')
specPix_nvalid = specPix[nvalid_x,nvalid_y,:] 

class IndexTracker(object):
    def __init__(self, ax, bins,specPix, x_ind, y_ind,dataframe):
        #loading in parameters
        self.ax = ax
        self.bins = bins
        self.specPix = specPix
        self.x_ind = x_ind
        self.y_ind = y_ind
        self.dataframe = dataframe
        #calculating number of indexes
        self.slices = len(self.x_ind)
        #setting initial index to 0 and setting plot title accordingly
        self.ind = 0
        #setting x and y axis titles
        self.ax.set_xlabel('ADU')
        self.ax.set_ylabel('Counts')
        #update title based on issues with indexed pixel
        self.title_update()
        #ploting spectrum of indexed pixel
        self.im, = ax.plot(bins,np.squeeze(self.specPix[self.ind,:]))
        self.update()

    def onscroll(self, event):
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        self.im.set_data(bins,np.squeeze(self.specPix[self.ind,:]))
        self.title_update()
        self.im.axes.figure.canvas.draw()
    def title_update(self):
        self.title_str = 'Pixel ({0},{1}): Issues with '.format(self.x_ind[self.ind],self.y_ind[self.ind])
        num_issues = 0
        for ind,col in enumerate(self.dataframe.columns):
            if pd.isnull(self.dataframe[col][self.x_ind[self.ind],self.y_ind[self.ind]])==False:
                if num_issues==0:
                    self.title_str += col
                else:
                    self.title_str += ', ' + col
                num_issues += 1
        self.ax.set_title(self.title_str)


fig,ax = plt.subplots(1,1)
tracker_FWHM = IndexTracker(ax,bins,specPix_nvalid,nvalid_x,nvalid_y,non_valid_dataframe)
fig.canvas.mpl_connect('scroll_event', tracker_FWHM.onscroll)
plt.show()

nValid = len(np.where(valid==1)[0])
pcValid = 100*nValid/6400

fig1 = plt.figure(figsize=(width,height))
fig1.add_subplot(1,3,1)
plt.imshow(PeakHeights,cmap='jet')
plt.colorbar(label='Counts',fraction=0.046, pad=0.04)
plt.title('Peak Heights')
fig1.add_subplot(1,3,2)
plt.imshow(PeakPos,cmap='jet')
plt.colorbar(label='ADU',fraction=0.046, pad=0.04)
plt.title('Peak Positions')
fig1.add_subplot(1,3,3)
plt.imshow(valid,cmap='gray')
plt.colorbar(label='Valid',fraction=0.046, pad=0.04)
plt.clim(0,1)
plt.title(str("{:.2f}".format(pcValid))+'% Valid Based on Peak Height and Position')
plt.subplots_adjust(wspace=0.5)
plt.savefig(results_dir + '\Am60keV_' + AmExpName + '_' + CoExpName + '.png')

fig2 = plt.figure(figsize=(width,height))
fig2.add_subplot(1,3,1)
plt.imshow(threshold)
plt.colorbar(label='ADU',fraction=0.046, pad=0.04)
plt.clim(0,300)
plt.title('Thresholds')
fig2.add_subplot(1,3,2)
plt.imshow(gradient)
plt.colorbar(label='keV ADU$^{-1}$',fraction=0.046, pad=0.04)
plt.clim(0.025,0.028)
plt.title('Gradients')
fig2.add_subplot(1,3,3)
plt.imshow(intercept)
plt.colorbar(label='keV',fraction=0.046, pad=0.04)
plt.title('Intercepts')
plt.clim(0,5)
plt.subplots_adjust(wspace=0.5)
plt.savefig(results_dir + '\GradientIntercepts_' + AmExpName + '_' + CoExpName + '.png')

fig3 = plt.figure(figsize=(width,height))
fig3.add_subplot(1,3,1)
plt.imshow(PeakPos)
plt.colorbar(label='ADU',fraction=0.046, pad=0.04)
plt.title('59.54keV Peak Positions')
fig3.add_subplot(1,3,2)
plt.imshow(PeakPos2)
plt.colorbar(label='ADU',fraction=0.046, pad=0.04)
plt.title('17.75keV Peak Positions')
fig3.add_subplot(1,3,3)
plt.imshow(PeakPos3)
plt.title('13.98keV Peak Positions')
plt.colorbar(label='ADU',fraction=0.046, pad=0.04)
plt.subplots_adjust(wspace=0.5)
plt.savefig(results_dir + '\PeakPosPerPixel_' + AmExpName + '_' + CoExpName + '.png')

FWHMbins = np.arange(0,(2.5/gain_ratio),0.1)
FWHMdist = np.histogram(np.reshape(FWHM_keV,(6400,1)), bins = FWHMbins)[0]
nValidFWHM = len(np.where(validFWHM==1)[0])
pcValidFWHM = 100*nValidFWHM/6400

FWHMbins_Co = np.arange(0,(2.5/gain_ratio),0.1)
FWHMdist_Co = np.histogram(np.reshape(FWHM_keV_Co,(6400,1)),FWHMbins_Co)[0]

fig4 = plt.figure(figsize=(width,height))
fig4.add_subplot(1,3,1)
plt.imshow(FWHM_keV)
plt.colorbar(label='keV',fraction=0.046, pad=0.04)
plt.clim(0.5,2.5)
plt.title('FWHM')
plt.clim(0,2)
fig4.add_subplot(1,3,2)
plt.bar(FWHMbins[:-1],FWHMdist)
plt.title('Distribution of FWHM')
plt.xlim(0,(2.5/gain_ratio))
plt.ylim(0,3000)
plt.xlabel('FWHM (keV)')
plt.ylabel('Number of Pixels')
fig4.add_subplot(1,3,3)
plt.imshow(validFWHM,cmap='gray')
plt.colorbar(label='Valid',fraction=0.046, pad=0.04)
plt.title(str("{:.2f}".format(pcValidFWHM)) + '% Valid Based on FWHM after peak rejection')
plt.subplots_adjust(wspace=0.5)
plt.savefig(results_dir + '\Am241FWHM_'  + AmExpName + '_' + CoExpName + '.png')

fig5 = plt.figure(figsize=(width,height))
fig5.add_subplot(1,3,1)
plt.imshow(FWHM_keV_Co)
plt.colorbar(label='keV',fraction=0.046, pad=0.04)
plt.clim(0.5,2.5)
plt.title('Co-57 FWHM')
plt.clim(0,2)
fig5.add_subplot(1,3,2)
plt.bar(FWHMbins_Co[:-1],FWHMdist_Co)
plt.title('Distribution of Co-57 FWHM')
plt.xlim(0,(2.5/gain_ratio))
plt.ylim(0,3000)
plt.xlabel('FWHM (keV)')
plt.ylabel('Number of Pixels')
fig5.add_subplot(1,3,3)
plt.imshow(peakMax4)
plt.colorbar(label='Counts',fraction=0.046, pad=0.04)
plt.title('Co-57 Peak Height')
plt.subplots_adjust(wspace=0.5)
plt.savefig(results_dir + '\Co57FWHM_' + AmExpName + '_' + CoExpName + '.png')

#printing values
print(('The 60keV FWHM was {0}+/-{1}keV').format(np.mean(np.mean(FWHM_keV)),np.std(np.std(FWHM_keV))))
print(('The 122keV FWHM was {0}+/-{1}keV').format(np.mean(np.mean(FWHM_keV_Co)),np.std(np.std(FWHM_keV_Co))))
print(('The number of valid pixels was {0} which is {1}%').format(nValidFWHM,pcValidFWHM))

#saving the data
np.savetxt(results_dir + '/Gradients_' + AmExpName + '_' + CoExpName + '.txt', np.transpose(gradient),delimiter=' ')
np.savetxt(results_dir + '/Intercepts_' + AmExpName + '_' + CoExpName +  '.txt', np.transpose(intercept),delimiter= ' ')
np.savetxt(results_dir + '/Thresholds_' + AmExpName + '_' + CoExpName + '.txt' , np.transpose(threshold),delimiter = ' ')
   
font = 15
ticks=13
pixX = pixY = 40
fig10 = plt.figure()
plt.plot(bins,specPix[pixX,pixY],color='darkblue',linewidth=2)
plt.xlabel('ADU',fontsize=font)
plt.ylabel('Counts',fontsize=font)
plt.xlim(left=0,right=4000)
plt.ylim(bottom=0)
plt.xticks(fontsize=ticks)
plt.yticks(fontsize=ticks)
plt.title('Am241 binned spectrum for pixel ({0},{1})'.format(str(pixX),str(pixY)),fontsize=font)
plt.show()

font = 24
ticks=18
pixX = pixY = 40
fig10 = plt.figure()
plt.plot(bins,specPix[pixX,pixY],color='darkblue',linewidth=2)
plt.plot(PeakPos[pixX,pixY],PeakHeights[pixX,pixY],marker='o',markersize=10,linestyle='None',color='red',label='Am-241 59.54 keV peak position = {0} ADU'.format(int(PeakPos[pixX,pixY])))
plt.plot([LHWHM[pixX,pixY],UHWHM[pixX,pixY]],[PeakHeights[pixX,pixY]/2,PeakHeights[pixX,pixY]/2],label='FWHM = {0} ADU = {1} keV'.format(int(FWHM_ADU[pixX,pixY]),'{:.2f}'.format(FWHM_keV[pixX,pixY])),linestyle='--',color='green',marker='o',markersize=10)
plt.xlabel('ADU',fontsize=font)
plt.ylabel('Counts',fontsize=font)
plt.xlim(left=0,right=750)
plt.ylim(bottom=0)
plt.xticks(fontsize=ticks)
plt.yticks(fontsize=ticks)
plt.legend(fontsize=font-3)
#plt.title('Am241 binned spectrum for pixel ({0},{1})'.format(str(pixX),str(pixY)),fontsize=font)
plt.show()

fig11 = plt.figure()
plt.bar(FWHMbins[:-1],FWHMdist,color='darkblue')
plt.title('Distribution of FWHM in High-Gain mode',fontsize=font)
plt.xlim(0,2.5)
plt.ylim(0,3000)
plt.xticks(fontsize=ticks)
plt.yticks(fontsize=ticks)
plt.xlabel('FWHM (keV)',fontsize=font)
plt.ylabel('Number of Pixels',fontsize=font)
plt.show()

np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/FWHM_ADU_CO',FWHM_ADU_Co)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHM_keV'),FWHM_keV)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHM_keV_Co'),FWHM_keV_Co)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHMbins'),FWHMbins)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHMbins_Co'),FWHMbins_Co)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHMdist'),FWHMdist)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('FWHMdist_Co'),FWHMdist_Co)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('bins'),bins)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('binsCo'),binsCo)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('specPix'),specPix)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('specPixCo'),specPixCo)
np.save('C:/Users/uhf41485/Documents/Saved arrays/' + str(gain) + '/' + str('peakPos'),PeakPos)











                                                                           