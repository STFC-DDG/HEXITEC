# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:56:28 2020

@author: vqs78369
"""

import numpy
from matplotlib import pyplot as plt

import time
import skimage.measure as ski

from scipy.stats import norm
 
path = 'C:\\Users\\vqs78369\\Desktop\\Am241_600s_500V\\Am241_600s_500V.bin'
f = open(path, "r")

path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\gradient.npy'
gradient = numpy.load(path)
path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\intercept.npy'
intercept = numpy.load(path)

threshold = 100
nSeconds = 60   # for 600s maybe nSeconds = 60 and nLoops = 10
frame_rate = 1589
nFrames = nSeconds*frame_rate

overall_start_time = time.time()

import os
newpath = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\calibrated' 
if not os.path.exists(newpath):
    os.makedirs(newpath)

for nLoops in range(0,10):
    data_cube = numpy.empty((80,80,nFrames))
    timestamp = str(nLoops).zfill(4)
    for frame_index in range(0,nFrames):
        #frame_throw_away = numpy.fromfile(f, dtype=numpy.uint16, count=6)
        frame_values = numpy.fromfile(f, dtype=numpy.uint16, count=80*80)
        frame = numpy.zeros((80,80))
        for index in range(0,len(frame_values)):
            row = int(numpy.floor(index/80))
            column = int(20*(index%4) + numpy.floor((index-(80*numpy.floor(index/80)))/4))
            frame[row,column] = frame_values[index]
         
        #coords_over_thresh = numpy.nonzero(frame_values>threshold)[0]
        #for c in range(0,len(coords_over_thresh)):
        #    index = coords_over_thresh[c]
        #    row = int(numpy.floor(index/80))
        #    column = int(20*(index%4) + numpy.floor((index-(80*numpy.floor(index/80)))/4))
        #    frame[row,column] = frame_values[index]
            
        #frame = numpy.reshape(a, (80,80),'F')
        original_frame = frame
        data_cube[:,:,frame_index] = frame
        
    #mean_frame = numpy.mean(data_cube<100,axis=2)
    #for frame_index in range(0,nFrames):
    #    data_cube[:,:,frame_index] = data_cube[:,:,frame_index] - mean_frame
    

    noise_array = numpy.copy(data_cube)
    start_time = time.time()
    
    above_threshold = noise_array > threshold
    noise_array[above_threshold] = numpy.nan
    
    spatial_mean = numpy.zeros((80,80))
    standard_dev = numpy.zeros((80,80))                
    for x in range(0,80):
        for y in range(0,80):
            elements = noise_array[x,y,:]
            mu, std = norm.fit(elements[numpy.isfinite(elements)])
            spatial_mean[x,y] = mu
            standard_dev[x,y] = std

    del noise_array 
    
    end_time = time.time()
    print(end_time - start_time)
    
    threshold_per_pixel = spatial_mean+8*standard_dev

    save_path = newfile + '\\threshold_' + timestamp
    numpy.save(save_path,threshold_per_pixel)
    
    D = numpy.zeros((80,80,nFrames))
    for frame_index in range(0,nFrames):
        D[:,:,frame_index] = data_cube[:,:,frame_index] > threshold_per_pixel

    D = numpy.array(D, dtype=bool)    
    #D = data_cube > threshold
    E = numpy.zeros((80,80,nFrames))
    E[D] = 1
    percent_events = numpy.sum(E,2)
    percent_events = 100*percent_events/nFrames
    del D 
    del E
    
    hot_pixels = percent_events > 2
    data_cube[hot_pixels] = 0

    bins = numpy.arange(0,8000,10)
    spatial_EOF = numpy.zeros((80,80))
    EOF_kept = numpy.zeros((len(bins),1))
    EOF_thrown = numpy.zeros((len(bins),1))
    
    start_time = time.time()
    for N in range(0,nFrames-1):
        frame = data_cube[:,:,N:N+2]
        events = numpy.zeros((80,80,2))
        for k in range(0,2):
            events[:,:,k] = frame[:,:,k] > threshold_per_pixel
        events = numpy.array(events, dtype=bool)
        #events = frame > threshold
        
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
    
    save_path = newpath + '\\EOF_kept_' + timestamp
    numpy.save(save_path,EOF_kept)
    save_path = newpath + '\\EOF_thrown_' + timestamp
    numpy.save(save_path,EOF_thrown)
    save_path = newpath + '\\EOF_spatial_' + timestamp
    numpy.save(save_path,spatial_EOF)
    
    end_time = time.time()
    print(end_time - start_time)
    
    start_time = time.time()
    bins = numpy.arange(0,150,0.2)
    
    spectra_per_pixel_raw = numpy.zeros((80,80,len(bins)))
    spectra_per_pixel_CSD = numpy.zeros((80,80,len(bins)))
    spectra_per_pixel_CSA = numpy.zeros((80,80,len(bins)))
    
    multiplicity = numpy.zeros((80,80,25))
    
    for N in range(0,nFrames):
        frame = numpy.copy(data_cube[:,:,N])
        events = frame > threshold_per_pixel
        #events = frame > threshold
        
        frame = frame*gradient + intercept
        
        event_coords = list(zip(*numpy.where(events)))
        
        # raw spectra
        for (x, y) in event_coords:
            event_val = frame[x, y]
            if event_val < 150:
                ind = int(numpy.floor(event_val/0.2))
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
            if S < 150 and cluster_size < 25:
                ind = int(numpy.floor(S/0.2))
                spectra_per_pixel_CSA[max_index[0],max_index[1],ind] += 1
                multiplicity[max_index[0],max_index[1],cluster_size] +=1
            
            # CSD (delete any charge share events - fewer events but smaller FWHM)
            if cluster_size == 1 and S < 150:
                ind = int(numpy.floor(S/0.2))
                spectra_per_pixel_CSD[max_index[0],max_index[1],ind] += 1
    
    end_time = time.time()
    print(end_time - start_time)
    
    #total_raw_spec = numpy.sum(numpy.sum(spectra_per_pixel_raw,0),0)
    #total_CSD_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSD,0),0)
    #total_CSA_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSA,0),0)
    
    #plt.plot(bins[0:500],total_raw_spec[0:500],label='Raw')
    #plt.plot(bins[0:500],total_CSA_spec[0:500],label='CSA')
    #plt.plot(bins[0:500],total_CSD_spec[0:500],label='CSD',linestyle='dashed')
    #plt.legend()
    #plt.show()
    
    save_path = newpath + '\\raw_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_raw)
    save_path = newpath + '\\CSD_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_CSD)
    save_path = newpath + '\\CSA_spectrum_ADU_' + timestamp
    numpy.save(save_path,spectra_per_pixel_CSA)
    
    save_path = newpath + '\\multiplicity_' + timestamp
    numpy.save(save_path,multiplicity)   
    
    
    energies = numpy.zeros((3,1))    
    energies[0] = 13.94
    energies[1] = 17.75
    energies[2] = 59.54
    
    search_indices = numpy.zeros((3,2))
    search_indices[0,0] = 62
    search_indices[1,0] = 77
    search_indices[2,0] = 274
    
    search_indices[0,1] = 77
    search_indices[1,1] = 98
    search_indices[2,1] = 315
    
    s = int(search_indices[2,0])
    e = int(search_indices[2,1])
    
    full_width = numpy.zeros((80,80))
    dif_FWHM = numpy.zeros((80,80))
    funny = 0
    for x in range(0,80):
        for y in range(0,80):
            peak_index = numpy.argmax(spectra_per_pixel_CSD[x,y,s:e]) + s
            peak_value = numpy.max(spectra_per_pixel_CSD[x,y,s:e])

            test = spectra_per_pixel_CSD[x,y,:] - peak_value/2
            test2 = abs(test)
            if peak_index > s and peak_index < e:
                upper_FWHM_index = numpy.argmin(test2[peak_index:e]) + peak_index
                lower_FWHM_index = numpy.argmin(test2[s:peak_index]) + s
                upper_FWHM_value = spectra_per_pixel_CSD[x,y,upper_FWHM_index]
                lower_FWHM_value = spectra_per_pixel_CSD[x,y,lower_FWHM_index]
                full_width[x,y] = 0.2*(upper_FWHM_index-lower_FWHM_index)
                dif_FWHM[x,y] = upper_FWHM_value-lower_FWHM_value
            else: 
                funny += 1
                    
    #plt.figure(0)    
    #plt.imshow(full_width)
    #plt.title('FWHM per pixel [ADU]')
    #cbar = plt.colorbar()
    #plt.show()
    #plt.figure(1)
    #plt.imshow(dif_FWHM)
    #plt.title('Difference between upper and lower sides of FWHM [ADU]')
    #cbar = plt.colorbar()
    #plt.show()  
    
    save_path = newpath + '\\FWHM_' + timestamp
    numpy.save(save_path,full_width)      
    
overall_end_time = time.time()
print(overall_end_time-overall_start_time)