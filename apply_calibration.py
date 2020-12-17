# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:49:36 2020

@author: vqs78369
"""

import numpy
from matplotlib import pyplot as plt

bins = numpy.arange(0,8000,10)
total_raw_spectrum = numpy.zeros((80,80,len(bins)))
total_CSA_spectrum = numpy.zeros((80,80,len(bins)))
total_CSD_spectrum = numpy.zeros((80,80,len(bins)))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\raw_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_raw = numpy.load(path)
    total_raw_spectrum += spectra_per_pixel_raw
    path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\CSA_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_CSA = numpy.load(path)
    total_CSA_spectrum += spectra_per_pixel_CSA
    path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\uncalibrated\\CSD_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_CSD = numpy.load(path)
    total_CSD_spectrum += spectra_per_pixel_CSD
    
#%%
total_raw = sum(sum(total_raw_spectrum,0),0)
plt.plot(bins[0:250],total_raw[0:250],label='Raw')
total_CSA = sum(sum(total_CSA_spectrum,0),0)
plt.plot(bins[0:250],total_CSA[0:250],label='CSA')
total_CSD = sum(sum(total_CSD_spectrum,0),0)
plt.plot(bins[0:250],total_CSD[0:250],label='CSD',linestyle='dashed')
plt.legend()
plt.show()

#%%
from scipy import stats
energies = numpy.zeros(3)
energies[0] = 13.94
energies[1] = 17.75
energies[2] = 59.54

# Need to change search indices for Co-57 source
search_indices = numpy.zeros((3,2))
search_indices[0,0] = 35
search_indices[1,0] = 45
search_indices[2,0] = 188

search_indices[0,1] = 45
search_indices[1,1] = 60
search_indices[2,1] = 210

ADU_indices = numpy.zeros((80,80,3))
ADU_max_values = numpy.zeros((80,80,3))
gradient = numpy.zeros((80,80))
intercept = numpy.zeros((80,80))

r_value = numpy.zeros((80,80))

for x in range(0,80):
    for y in range(0,80):
        for i in range(0,3):
            s = int(search_indices[i,0])
            e = int(search_indices[i,1])
            ADU_indices[x,y,i] = 10*(numpy.argmax(total_CSD_spectrum[x,y,s:e]) + s)
            ADU_max_values[x,y,i] = numpy.max(total_CSD_spectrum[x,y,s:e])
            #fit = numpy.polyfit(numpy.squeeze(ADU_indices[x,y,:]), energies, 1, full=True)
            gradient_single_pixel, intercept_single_pixel, r_value_single_pixel, p_value, std_err = stats.linregress(numpy.squeeze(ADU_indices[x,y,:]), energies)
            gradient[x,y] = gradient_single_pixel
            intercept[x,y] = intercept_single_pixel
            r_value[x,y] = r_value_single_pixel
    
    
plt.figure(0)
plt.imshow(gradient)
plt.title('Gradient')
plt.show()
plt.figure(1)
plt.imshow(intercept)
plt.title('Intercept')
plt.show()
plt.figure(2)
plt.imshow(r_value)
plt.title('R value')
plt.show()

save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\gradient.npy'
numpy.save(save_path,gradient)
save_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\intercept.npy'
numpy.save(save_path,intercept)