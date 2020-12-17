# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:06:44 2020

@author: vqs78369
"""

import numpy
from matplotlib import pyplot as plt

saved_path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\calibrated'

bins = numpy.arange(0,150,0.2)
total_raw_spectrum = numpy.zeros((80,80,len(bins)))
total_CSA_spectrum = numpy.zeros((80,80,len(bins)))
total_CSD_spectrum = numpy.zeros((80,80,len(bins)))
total_ECC_spectrum = numpy.zeros((80,80,len(bins)))
for nLoops in range(0,10):#360):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\raw_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_raw = numpy.load(path)
    total_raw_spectrum += spectra_per_pixel_raw
    path = saved_path + '\\CSA_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_CSA = numpy.load(path)
    total_CSA_spectrum += spectra_per_pixel_CSA
    path = saved_path + '\\CSD_spectrum_ADU_' + timestamp + '.npy'
    spectra_per_pixel_CSD = numpy.load(path)
    total_CSD_spectrum += spectra_per_pixel_CSD
    #path = saved_path + '\\ECC_spectrum_' + timestamp + '.npy'
    #spectra_per_pixel_ECC = numpy.load(path)
    #total_ECC_spectrum += spectra_per_pixel_ECC
    
#%%
total_raw = sum(sum(total_raw_spectrum,0),0)
plt.plot(bins[0:350],total_raw[0:350],label='Raw')
#total_ECC = sum(sum(total_ECC_spectrum,0),0)
#plt.plot(bins[0:350],total_ECC[0:350],label='ECC')
total_CSA = sum(sum(total_CSA_spectrum,0),0)
plt.plot(bins[0:350],total_CSA[0:350],label='CSA')
total_CSD = sum(sum(total_CSD_spectrum,0),0)
plt.plot(bins[0:350],total_CSD[0:350],label='CSD',linestyle='dashed')
plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
plt.legend()
plt.xlabel('Energy [keV]')
plt.ylabel('Number of events')
plt.title('Global Spectra')
plt.show()

#%%
total_ECC = sum(sum(total_ECC_spectrum,0),0)
plt.plot(bins[0:350],total_ECC[0:350],label='ECC')
total_CSA = sum(sum(total_CSA_spectrum,0),0)
plt.plot(bins[0:350],total_CSA[0:350],label='CSA')
plt.legend()
plt.xlabel('Energy [keV]')
plt.ylabel('Number of events')
plt.title('CdTe Am-241 Source')
plt.show()
    
#%% FWHM from CSD
search_indices = numpy.zeros((3,2))
#search_indices[0,0] = 35
#search_indices[1,0] = 45
#search_indices[2,0] = 188

#search_indices[0,1] = 45
#search_indices[1,1] = 60
#search_indices[2,1] = 210

search_indices[0,0] = 62
search_indices[1,0] = 77
search_indices[2,0] = 274

search_indices[0,1] = 77
search_indices[1,1] = 98
search_indices[2,1] = 315

s = int(search_indices[2,0])
e = int(search_indices[2,1])
from scipy.signal import find_peaks, peak_widths

FW = numpy.zeros((80,80))
funny = 0 

for x in range(0,80):
    for y in range(0,80):
        peaks, _ = find_peaks(total_CSD_spectrum[x,y,s:e])
        results_half = peak_widths(total_CSD_spectrum[x,y,s:e], peaks, rel_height=0.5)
        
        if len(results_half[0]) !=0:
            FW[x,y] = numpy.max(results_half[0])*0.2
#%%
plt.figure(0)
plt.imshow(FW)
plt.colorbar(label='FWHM [keV]')
plt.title('FWHM [kev]')
plt.show()
#%%

FW_min = numpy.min(FW)
FW_max = numpy.max(FW)
FW_bin_width = (FW_max-FW_min)/100
FW_bins = numpy.arange(FW_min,FW_max,FW_bin_width)
FW_counts = numpy.zeros((len(FW_bins),1))
for x in range(0,80):
    for y in range(0,80):
        for z in range(0,len(FW_bins)-1):
            if FW[x,y]>=FW_bins[z] and FW[x,y]<FW_bins[z+1]:
                FW_counts[z] += 1

FW_cumulative_counts = numpy.zeros((len(FW_bins),1))
for z in range(0,len(FW_bins)):
    FW_cumulative_counts[z] = numpy.sum(FW_counts[0:z])

plt.figure(1)
plt.plot(FW_bins,FW_cumulative_counts/6400)
plt.xlabel('FWHM [kev]')
plt.ylabel('Normalised Counts')
plt.title('Cumulative distribution of FWHM')
    

mean_FW= numpy.mean(FW[FW>0])
print('Mean FWHM = ' + str(mean_FW))

percent_bad_pixels = 100*numpy.sum(FW>2)/(80*80*4)
print('Percent bad pixels = ' + str(percent_bad_pixels))

#%%
total_multiplicity = numpy.zeros((80,80,25))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\multiplicity_' + timestamp + '.npy'   
    multiplicity = numpy.load(path)
    total_multiplicity += multiplicity

#%%
all_multiplicity = numpy.sum(numpy.sum(total_multiplicity,0),0)
plt.bar(numpy.arange(1,10),all_multiplicity[1:10],tick_label=numpy.arange(1,10))
plt.yscale('log')
plt.xlabel('Cluster Size [pixels]')
plt.ylabel('Number of events [log scale]')
plt.title('Multiplicity Histogram')

#%%
s = numpy.zeros((80,80))
n = numpy.zeros((80,80))
for m in range(0,25):
    s += total_multiplicity[:,:,m]*m
    n += total_multiplicity[:,:,m]

weighted_mean_multiplicity = s/n
plt.imshow(weighted_mean_multiplicity)
plt.colorbar()
plt.title('Mean Multiplicity')
plt.show()

mean_multiplicity = numpy.mean(weighted_mean_multiplicity[weighted_mean_multiplicity>0])
print('Mean multplicity = ' + str(mean_multiplicity))

#%%
one_and_two_pixel_percent = 100*numpy.sum(all_multiplicity[1:3])/numpy.sum(all_multiplicity)
print('Percent of clusters that are one and two pixels = ' + str(one_and_two_pixel_percent))
plt.pie(all_multiplicity[1:6],labels = ['n=1','n=2','n=3','n=4','n=5'])
plt.title('Multiplicity')

#%%
total_events_per_pixel = numpy.sum(total_multiplicity,2)
plt.imshow(total_events_per_pixel)
plt.colorbar()
plt.title('Total number of events per pixel')
plt.show()

mean_number_of_events = numpy.mean(total_events_per_pixel[total_events_per_pixel>0])
print('Mean events per pixel = ' + str(mean_number_of_events))
#%%
percent_CS_per_pixel = 100*(1-total_multiplicity[:,:,1]/total_events_per_pixel)
plt.imshow(percent_CS_per_pixel)
plt.colorbar(label='[%]')
plt.title('Percentage charge sharing per pixel')
plt.show()
mean_percent_CS_per_pixel = numpy.mean(percent_CS_per_pixel[percent_CS_per_pixel>0])
print('Mean charge share per pixel [%] = ' + str(mean_percent_CS_per_pixel))

#%%
path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing' + '\\gradient.npy'
gradient = numpy.load(path)
plt.figure(0)
plt.imshow(gradient,vmin=0.028, vmax=0.031)
plt.colorbar(label='Gradient [keV/ADU]')
plt.title('Gradient')
plt.show()

grad_min = numpy.min(gradient)
grad_max = numpy.max(gradient)
grad_bin_width = (grad_max-grad_min)/10
grad_bins = numpy.arange(grad_min,grad_max,grad_bin_width)
grad_counts = numpy.zeros((len(grad_bins),1))
for x in range(0,80):
    for y in range(0,80):
        for z in range(0,len(grad_bins)-1):
            if gradient[x,y]>=grad_bins[z] and gradient[x,y]<grad_bins[z+1]:
                grad_counts[z] += 1

grad_cumulative_counts = numpy.zeros((len(grad_bins),1))
for z in range(0,len(grad_bins)):
    grad_cumulative_counts[z] = numpy.sum(grad_counts[0:z])

plt.figure(1)
plt.plot(grad_bins,grad_cumulative_counts/6400)
plt.xlabel('Gradient')
plt.ylabel('Normalised Counts')
plt.title('Cumulative distribution of gradient')

mean_gradient = numpy.mean(gradient[gradient>0])
print('Mean gradient = ' + str(mean_gradient))
#%%
path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing' + '\\intercept.npy'
intercept = numpy.load(path)
plt.figure(0)
plt.imshow(intercept)
plt.colorbar(label='Intercept [keV]')
plt.title('Intercept')
plt.show()

intercept_min = numpy.min(intercept)
intercept_max = numpy.max(intercept)
intercept_bin_width = (intercept_max-intercept_min)/30
intercept_bins = numpy.arange(intercept_min,intercept_max,intercept_bin_width)
intercept_counts = numpy.zeros((len(intercept_bins),1))
for x in range(0,80):
    for y in range(0,80):
        for z in range(0,len(intercept_bins)-1):
            if intercept[x,y]>=intercept_bins[z] and intercept[x,y]<intercept_bins[z+1]:
                intercept_counts[z] += 1

intercept_cumulative_counts = numpy.zeros((len(intercept_bins),1))
for z in range(0,len(intercept_bins)):
    intercept_cumulative_counts[z] = numpy.sum(intercept_counts[0:z])

plt.figure(1)
plt.plot(intercept_bins,intercept_cumulative_counts/6400)
plt.xlabel('Intercept')
plt.ylabel('Normalised Counts')
plt.title('Cumulative distribution of intercept')


mean_intercept = numpy.mean(intercept[intercept>0])
print('Mean intercept = ' + str(mean_intercept))

#%%
total_threshold = numpy.zeros((80,80))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\threshold_' + timestamp + '.npy'
    threshold = numpy.load(path)
    total_threshold += threshold

total_threshold = total_threshold/10
plt.figure(0)
plt.imshow(total_threshold)
plt.colorbar(label='Threshold [ADU]')
plt.title('Threshold [ADU]')
plt.show()

threshold_min = 0
threshold_max = 100
threshold_bin_width = (threshold_max-threshold_min)/100
threshold_bins = numpy.arange(threshold_min,threshold_max,threshold_bin_width)
threshold_counts = numpy.zeros((len(threshold_bins),1))
for x in range(0,80):
    for y in range(0,80):
        for z in range(0,len(threshold_bins)-1):
            if total_threshold[x,y]>=threshold_bins[z] and total_threshold[x,y]<threshold_bins[z+1]:
                threshold_counts[z] += 1

threshold_cumulative_counts = numpy.zeros((len(threshold_bins),1))
for z in range(0,len(threshold_bins)):
    threshold_cumulative_counts[z] = numpy.sum(threshold_counts[0:z])

plt.figure(1)
plt.plot(threshold_bins,threshold_cumulative_counts/6400)
plt.xlabel('Threshold [ADU]')
plt.ylabel('Normalised Counts')
plt.title('Cumulative distribution of threshold [ADU]')
#%%
mean_total_threshold = numpy.mean(total_threshold[total_threshold>0])
print('Mean total_threshold [ADU] = ' + str(mean_total_threshold))

total_threshold = total_threshold*gradient + intercept
plt.figure(2)
plt.imshow(total_threshold)
plt.colorbar(label='Threshold [keV]')
plt.title('Threshold [keV]')
plt.show()

threshold_min = 0
threshold_max = 10
threshold_bin_width = (threshold_max-threshold_min)/100
threshold_bins = numpy.arange(threshold_min,threshold_max,threshold_bin_width)
threshold_counts = numpy.zeros((len(threshold_bins),1))
for x in range(0,80):
    for y in range(0,80):
        for z in range(0,len(threshold_bins)-1):
            if total_threshold[x,y]>=threshold_bins[z] and total_threshold[x,y]<threshold_bins[z+1]:
                threshold_counts[z] += 1

threshold_cumulative_counts = numpy.zeros((len(threshold_bins),1))
for z in range(0,len(threshold_bins)):
    threshold_cumulative_counts[z] = numpy.sum(threshold_counts[0:z])

plt.figure(3)
plt.plot(threshold_bins,threshold_cumulative_counts/6400)
plt.xlabel('Threshold [keV]')
plt.ylabel('Normalised Counts')
plt.title('Cumulative distribution of threshold [keV]')

mean_total_threshold = numpy.mean(total_threshold[total_threshold>0])
print('Mean total_threshold [keV] = ' + str(mean_total_threshold))

#%%
EOF_bins = numpy.arange(0,8000,10)
total_EOF_kept = numpy.zeros((800,1))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\EOF_kept_' + timestamp + '.npy'
    EOF_kept = numpy.load(path)
    total_EOF_kept += EOF_kept

total_EOF_thrown = numpy.zeros((800,1))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\EOF_thrown_' + timestamp + '.npy'
    EOF_thrown = numpy.load(path)
    total_EOF_thrown += EOF_thrown

plt.plot(EOF_bins[0:250],total_EOF_kept[0:250],label='Events kept')
plt.plot(EOF_bins[0:250],total_EOF_thrown[0:250],label='Events removed')
plt.legend()
plt.ticklabel_format(style='sci', axis='y',scilimits=(0,0))
plt.xlabel('Energy [ADU]')
plt.ylabel('Number of events')
plt.title('EOF corrections')
plt.show()
#%%
total_EOF_spatial = numpy.zeros((80,80))
for nLoops in range(0,10):
    if nLoops%10 == 0:
        print(nLoops)
    timestamp = str(nLoops).zfill(4)
    path = saved_path + '\\EOF_spatial_' + timestamp + '.npy'
    EOF_spatial = numpy.load(path)
    total_EOF_spatial += EOF_spatial

EOF_percent = numpy.zeros((80,80))
for x in range(0,80):
    for y in range(0,80):
        EOF_percent[x,y] = 100*2*total_EOF_spatial[x,y]/(numpy.sum(total_raw_spectrum[x,y,:])+total_EOF_spatial[x,y])

plt.imshow(EOF_percent)
plt.colorbar(label='[%]')
plt.title('Percentage EOF events')

mean_EOF_percent = numpy.mean(EOF_percent[EOF_percent>0])
print('Mean EOF percent = ' + str(mean_EOF_percent))
#%%
for i in range(0,10):
    x = random.randint(0,80)
    y = random.randint(0,80)
    peaks, _ = find_peaks(total_CSD_spectrum[x,y,s:e])
    results_half = peak_widths(total_CSD_spectrum[x,y,s:e], peaks, rel_height=0.5)
    plt.figure(i)
    plt.plot(total_CSD_spectrum[x,y,s:e])
    plt.plot(peaks, total_CSD_spectrum[x,y,s:e][peaks], "x")
    plt.hlines(*results_half[1:], color="C2")
    plt.title('x = ' + str(x) + '  y = ' + str(y))
    
