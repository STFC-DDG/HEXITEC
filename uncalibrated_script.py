# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 11:25:38 2020

@author: vqs78369
"""

import numpy
from matplotlib import pyplot as plt

import time
import skimage.measure as ski
 
path = 'C:\\Users\\vqs78369\\Desktop\\Am241_600s_500V\\Am241_600s_500V.bin'
f = open(path, "r")

threshold = 100
nSeconds = 60   # for 600s maybe nSeconds = 60 and nLoops = 10
frame_rate = 1589
nFrames = nSeconds*frame_rate
#nLoops = 10


absolute_total = 0
#throw_away = numpy.fromfile(f, dtype=numpy.uint16, count=24)

overall_start_time = time.time()

import os
newpath = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

for nLoops in range(0,10):
    data_cube = numpy.empty((80,80,nFrames))
    timestamp = str(nLoops).zfill(4)
    for frame_index in range(0,nFrames):
        #frame_throw_away = numpy.fromfile(f, dtype=numpy.uint16, count=6)
        frame_values = numpy.fromfile(f, dtype=numpy.uint16, count=80*80)
        frame = numpy.zeros((80,80))
        
        coords_over_thresh = numpy.nonzero(frame_values>threshold)[0]
        for c in range(0,len(coords_over_thresh)):
            index = coords_over_thresh[c]
            row = int(numpy.floor(index/80))
            column = int(20*(index%4) + numpy.floor((index-(80*numpy.floor(index/80)))/4))
            frame[row,column] = frame_values[index]
            
        #frame = numpy.reshape(a, (80,80),'F')
        original_frame = frame
        data_cube[:,:,frame_index] = frame
        
    #mean_frame = numpy.mean(data_cube,axis=2)
    #for frame_index in range(0,nFrames):
    #    data_cube[:,:,frame_index] = data_cube[:,:,frame_index] - mean_frame
        
    D = data_cube > threshold
    E = numpy.zeros((80,80,nFrames))
    E[D] = 1
    percent_events = numpy.sum(E,2)
    percent_events = 100*percent_events/nFrames
    del D 
    del E
    
    hot_pixels = percent_events > 2   #change to adjust definition of a hot pixel
    data_cube[hot_pixels] = 0

    D = data_cube > threshold
    E = numpy.zeros((80,80,nFrames))
    E[D] = 1
    total_number_of_events = sum(sum(sum(E)))
    absolute_total += total_number_of_events
    spatial_number_of_events = numpy.sum(E,axis=2)

    bins = numpy.arange(0,8000,10)
    spatial_EOF = numpy.zeros((80,80))
    EOF_kept = numpy.zeros((len(bins),1))
    EOF_thrown = numpy.zeros((len(bins),1))
    
    start_time = time.time()
    for N in range(0,nFrames-1):
        frame = data_cube[:,:,N:N+2]
        events = frame > threshold
        
        eof_correction_frame = numpy.logical_and(
            events[:, :, 0],
            events[:, :, 1]
        )
        
        eof_correction_events = numpy.where(eof_correction_frame)
        
        for (x, y) in zip(*eof_correction_events):
            
            spatial_EOF[x,y] += 1
            
            ind_1 = int(numpy.floor(frame[x,y,0]/10))
            ind_2 = int(numpy.floor(frame[x,y,1]/10))
            
            if 0.95 <= frame[x,y,0]/frame[x,y,1] <= 1.05:
                data_cube[x,y,N] = 0
                data_cube[x,y,N+1] = 0
            elif frame[x,y,0] > frame[x,y,1]:
                data_cube[x,y,N+1] = 0
                EOF_kept[ind_1,0] += 1
                EOF_thrown[ind_2,0] += 1
            elif frame[x,y,1] > frame[x,y,0]:
                data_cube[x,y,N] = 0
                EOF_kept[ind_2,0] += 1
                EOF_thrown[ind_1,0] += 1
    
    #save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\EOF_kept_' + timestamp
    #numpy.save(save_path,EOF_kept)
    #save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\EOF_thrown_' + timestamp
    #numpy.save(save_path,EOF_thrown)
    #save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\EOF_spatial_' + timestamp
    #numpy.save(save_path,spatial_EOF)
    
    end_time = time.time()
    print(end_time - start_time)
    
    start_time = time.time()
    bins = numpy.arange(0,8000,10)
    
    spectra_per_pixel_raw = numpy.zeros((80,80,len(bins)))
    spectra_per_pixel_CSD = numpy.zeros((80,80,len(bins)))
    spectra_per_pixel_CSA = numpy.zeros((80,80,len(bins)))
    
    #multiplicity = numpy.zeros((162,162,25))
    
    for N in range(0,nFrames):
        frame = numpy.copy(data_cube[:,:,N])
        events = frame > threshold
        event_coords = list(zip(*numpy.where(events)))
        
        # raw spectra
        for (x, y) in event_coords:
            event_val = frame[x, y]
            if event_val < 8000:
                ind = int(numpy.floor(event_val/10))
                spectra_per_pixel_raw[x, y, ind] += 1
                
        binary = numpy.zeros((80,80))
        binary[events] = 1
        
        labels = ski.label(binary)
        regions = ski.regionprops(labels, frame)
        for region in regions:
            coord_pairs = region.coords
            cluster_size = len(coord_pairs)
            
            cluster_values = [frame[x, y] for (x, y) in coord_pairs]
            
            S = numpy.sum(cluster_values)
            
            max_index = coord_pairs[numpy.argmax(cluster_values)]
            
            # CSA (add up charge share events - more events but wider FWHM)
            if S < 8000:
                ind = int(numpy.floor(S/10))
                spectra_per_pixel_CSA[max_index[0],max_index[1],ind] += 1
                #multiplicity[max_index[0],max_index[1],cluster_size] +=1
            
            # CSD (delete any charge share events - fewer events but smaller FWHM)
            if cluster_size == 1:
                ind = int(numpy.floor(S/10))
                spectra_per_pixel_CSD[max_index[0],max_index[1],ind] += 1
    
    end_time = time.time()
    print(end_time - start_time)
    
    total_raw_spec = numpy.sum(numpy.sum(spectra_per_pixel_raw,0),0)
    total_CSD_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSD,0),0)
    total_CSA_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSA,0),0)
    
    plt.figure(nLoops)
    plt.plot(bins[0:500],total_raw_spec[0:500],label='Raw')
    plt.plot(bins[0:500],total_CSA_spec[0:500],label='CSA')
    plt.plot(bins[0:500],total_CSD_spec[0:500],label='CSD',linestyle='dashed')
    plt.legend()
    plt.show()
    
    save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\raw_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_raw)
    save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\CSD_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_CSD)
    save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\CSA_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_CSA)

overall_end_time = time.time()
print(overall_end_time-overall_start_time)

#%%
#plt.figure(2)
#plt.plot(bins[0:500],EOF_kept[0:500])
#plt.plot(bins[0:500],EOF_thrown[0:500])
#plt.show()