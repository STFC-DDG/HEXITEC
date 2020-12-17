# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 12:25:39 2020

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
nSeconds = 5   # for 600s maybe nSeconds = 60 and nLoops = 10
frame_rate = 1589
nFrames = nSeconds*frame_rate

overall_start_time = time.time()

#fit = numpy.polyfit(curve[:,0],curve[:,1],2)

#import os
#newpath = 'I:\\20201015_145751\\calibrated' 
#if not os.path.exists(newpath):
#   os.makedirs(newpath)

N_pairs = 0
R_E = numpy.zeros((1000000,2))   #391068

#N_pairs_over_peak = 0
#R_E_over_peak = numpy.zeros((68384,3))

XRF_energies = numpy.array([27.4723,
27.2017,
26.0955,
23.1736,
22.9841
])
    
XRF_energies_with_tolerance = numpy.zeros((len(XRF_energies),2))
for i in range(0,len(XRF_energies)):
    XRF_energies_with_tolerance[i,0] = XRF_energies[i]*0.9
    XRF_energies_with_tolerance[i,1] = XRF_energies[i]*1.1

#N_pairs_1 = 0
#N_pairs_2 = 0
#RE_1 = numpy.zeros((6550,2)) # 10% tolerance 6550    1% tolerance 714
#RE_2 = numpy.zeros((6355,2)) # 10% tolerance 6355    1% tolerance 700
    
#N_fl = 0
#RE_fl = numpy.zeros((4143,2))

#N_fl_1 = 0
#N_fl_2 = 0
#RE_fl_1 = numpy.zeros((9671,2))
#RE_fl_2 = numpy.zeros((6994,2))
    
for nLoops in range(0,1):
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
            
    
    data_cube[:,:,frame_index] = frame

    #noise_array = numpy.copy(data_cube)
    #start_time = time.time()
    
    #above_threshold = noise_array > threshold
    #noise_array[above_threshold] = numpy.nan
    
    #spatial_mean = numpy.zeros((162,162))
    #standard_dev = numpy.zeros((162,162))                
    #for x in range(0,162):
    #        # Ignore blank columns
    #        if x in [80, 81]:
    #            continue
    #        for y in range(0,162):
    #            # Ignore blank rows
    #            if y in [80, 81]:
    #                continue
    #            elements = noise_array[x,y,:]
    #            mu, std = norm.fit(elements[numpy.isfinite(elements)])
    #            spatial_mean[x,y] = mu
    #            standard_dev[x,y] = std
    
#    del noise_array 
#    
#    end_time = time.time()
#    print(end_time - start_time)
    
#    threshold_per_pixel = spatial_mean+8*standard_dev
    
    path = 'C:\\Users\\vqs78369\\Desktop\\Github_Testing\\calibrated\\threshold_' + timestamp + '.npy'
    threshold_per_pixel = numpy.load(path)
    
    D = numpy.zeros((80,80,nFrames))
    for frame_index in range(0,nFrames):
        D[:,:,frame_index] = data_cube[:,:,frame_index] > threshold_per_pixel
    
    D = numpy.array(D, dtype=bool)    
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
    
    end_time = time.time()
    print(end_time - start_time)
    
    start_time = time.time()
    bins = numpy.arange(0,150,0.2)
    
    #spectra_per_pixel_raw = numpy.zeros((162,162,len(bins)))
    #spectra_per_pixel_CSD = numpy.zeros((162,162,len(bins)))
    #spectra_per_pixel_CSA = numpy.zeros((162,162,len(bins)))
    spectra_per_pixel_CSA_adj = numpy.zeros((80,80,len(bins)))
    
    multiplicity = numpy.zeros((80,80,25))
    
    for N in range(0,nFrames):
        frame = numpy.copy(data_cube[:,:,N])
        events = frame > threshold_per_pixel
        
        frame = frame*gradient + intercept
        
#        event_coords = list(zip(*numpy.where(events)))
        
#        # raw spectra
#        for (x, y) in event_coords:
#            event_val = frame[x, y]
#            if event_val < 150:
#                ind = int(numpy.floor(event_val/0.2))
#                spectra_per_pixel_raw[x, y, ind] += 1
                
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
            
#            # CSA (add up charge share events - more events but wider FWHM)
#            if S < 150 and cluster_size < 25:
#                ind = int(numpy.floor(S/0.2))
#                spectra_per_pixel_CSA[max_index[0],max_index[1],ind] += 1
#                multiplicity[max_index[0],max_index[1],cluster_size] +=1
                
            # corrected CSA (with fudge factor)
            if S < 150:
                if cluster_size == 2 and 2/3<S/59.5<1:
                    R = (cluster_values[0]-cluster_values[1])/59.5
                    R_E[N_pairs,0] = R
                    R_E[N_pairs,1] = S
                    #R_E[N_pairs,2] = S_adj
                    N_pairs += 1
                else:
                    ind = int(numpy.floor(S/0.2))
                    spectra_per_pixel_CSA_adj[max_index[0],max_index[1],ind] += 1
            
#            # CSD (delete any charge share events - fewer events but smaller FWHM)
#            if cluster_size == 1 and S < 150:
#                ind = int(numpy.floor(S/0.2))
#                spectra_per_pixel_CSD[max_index[0],max_index[1],ind] += 1
    

overall_end_time = time.time()
print(overall_end_time-overall_start_time)
#%%
#total_raw_spec = numpy.sum(numpy.sum(spectra_per_pixel_raw,0),0)
#total_CSD_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSD,0),0)
#total_CSA_spec = numpy.sum(numpy.sum(spectra_per_pixel_CSA,0),0)

#total_CSA_spec_adj = numpy.sum(numpy.sum(spectra_per_pixel_CSA_adj,0),0)

R_E = R_E[0:N_pairs,:]
import numpy as np

from sklearn.cluster import DBSCAN
#from sklearn import metrics
#from sklearn.datasets import make_blobs
#from sklearn.preprocessing import StandardScaler


# #############################################################################
# Generate sample data
#centers = [[1, 1], [-1, -1], [1, -1]]
#X, labels_true = make_blobs(n_samples=750, centers=centers, cluster_std=0.4,
                     #       random_state=0)

#X = StandardScaler().fit_transform(X)

# #############################################################################
# Compute DBSCAN db = DBSCAN(eps=0.075, min_samples=85).fit(R_E)     (eps=0.05, min_samples=200)
db = DBSCAN(eps=0.05, min_samples=20).fit(R_E)
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)

print('Estimated number of clusters: %d' % n_clusters_)
print('Estimated number of noise points: %d' % n_noise_)
#print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
#print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
#print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
#print("Adjusted Rand Index: %0.3f"
 #     % metrics.adjusted_rand_score(labels_true, labels))
#print("Adjusted Mutual Information: %0.3f"
#      % metrics.adjusted_mutual_info_score(labels_true, labels))
#print("Silhouette Coefficient: %0.3f"
 #     % metrics.silhouette_score(X, labels))
#%%
R = R_E[labels==0][:,0]
E = R_E[labels==0][:,1]
bins_0 = numpy.arange(numpy.min(R),numpy.max(R),(numpy.max(R)-numpy.min(R))/15)
energy_vals_0 = numpy.zeros((len(bins_0),1))
number_of_values_0 = numpy.zeros((len(bins_0),1))
for i in range(0,len(R)):
    for b in range(0,len(bins_0)-1):
        if R[i]>=bins_0[b] and R[i]<bins_0[b+1]:
            energy_vals_0[b] += E[i]
            number_of_values_0[b] += 1

mean_energy_0 = energy_vals_0/number_of_values_0
plt.plot(bins_0,mean_energy_0)
#fit_0 = numpy.polyfit(bins_0,mean_energy_0,2)

R = R_E[labels==1][:,0]
E = R_E[labels==1][:,1]
bins_1 = numpy.arange(numpy.min(R),numpy.max(R),(numpy.max(R)-numpy.min(R))/15)
energy_vals_1 = numpy.zeros((len(bins_1),1))
number_of_values_1 = numpy.zeros((len(bins_1),1))
for i in range(0,len(R)):
    for b in range(0,len(bins_1)-1):
        if R[i]>=bins_1[b] and R[i]<bins_1[b+1]:
            energy_vals_1[b] += E[i]
            number_of_values_1[b] += 1

mean_energy_1 = energy_vals_1/number_of_values_1
plt.plot(bins_1,mean_energy_1)
plt.show()
#fit_1 = numpy.polyfit(bins_1,mean_energy_1,2)

bins_0 = bins_0[0:14]
bins_1 = bins_1[0:14]
mean_energy_0 = mean_energy_0[0:14]
mean_energy_1 = mean_energy_1[0:14]
bins = numpy.concatenate((bins_0,bins_1))
mean_energies = numpy.concatenate((mean_energy_0,mean_energy_1))
fit = numpy.polyfit(bins,mean_energies,2)
X = numpy.arange(-1,1,1/100)
Y = fit[0]*X*X + fit[1]*X + fit[2]
plt.plot(X,Y)
plt.plot(bins,mean_energies,marker='x', linestyle='')
plt.show()

ECC = 59.54 - (fit[0]*X*X + fit[1]*X + fit[2])

#%%
# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = ['deepskyblue','deepskyblue','limegreen','k','k','k','k','k','k','k','k','k','k']
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = (labels == k)

    xy = R_E[class_member_mask & core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=4)
#tuple(col)
    xy = R_E[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=1)

#plt.title('Estimated number of clusters: %d' % n_clusters_)
plt.plot(X,Y,linewidth=4,color='crimson')
plt.xlabel('R = ($E_1$-$E_2$)/$E_0$')
plt.ylabel('$E_1$ + $E_2$ [keV]')
plt.title('CdTe Am-241 Source')  
plt.show()
#%%

bins = numpy.concatenate((bins_0,bins_1))
mean_energies_corrected = 59.5-numpy.concatenate((mean_energy_0,mean_energy_1))
plt.figure(1)
plt.plot(X,ECC,'k')
plt.xlabel('R = ($E_1$-$E_2$)/$E_0$')
plt.ylabel('ECC(R) [keV]')
plt.title('CdTe Am-241 Source')  
plt.show()

#%%
plt.plot(R_E[:,0],R_E[:,1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=1)
plt.xlabel('R = ($E_1$-$E_2$)/$E_0$')
plt.ylabel('$E_1$ + $E_2$ [keV]')
plt.title('CZT Am-241 Source')  
plt.show()