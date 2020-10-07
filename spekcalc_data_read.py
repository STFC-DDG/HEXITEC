# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 14:49:14 2020

@author: uhf41485
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.table import Table

os.chdir(r'C:\Users\uhf41485\Documents\Dan Ryan simulation\hexitec_simulation-master\hexitec_simulation-master')
from hexitec_simulation import *

file_full = open(r'C:\Users\uhf41485\Documents\Dan Ryan simulation\100kVp 30deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Ti 0C 0Wa.spec', 'r')
content = file_full.readlines()

for i in range(len(content)):
    if 'Energy[keV]' in content[i]:
        print('True')
        start_index = i+1
        break

#max energy of spekcalc file
max_energy = 100
#setting the number of photons to be incident upon the detector
num_photons = 4000
#setting the threshold energy for an event to be registered
threshold = 4*u.keV

#loading in the spectrum from the spekcalc file and converting into an appropriate astropy table
energies_counts = np.loadtxt(r'C:\Users\uhf41485\Documents\Dan Ryan simulation\100kVp 30deg 1000Air 0Be 1Al 0Cu 0Sn 0W 0Ta 0Ti 0C 0Wa.spec',skiprows=start_index,delimiter='  ')
diff = energies_counts[1,0] - energies_counts[0,0]
low = (energies_counts[:,0]+(energies_counts[:,0]-1))/2
up =  (energies_counts[:,0]+(energies_counts[:,0]+1))/2
energy_counts_adjusted = (energies_counts[:,1]/energies_counts[:,1].sum())*num_photons
incident_spectrum = Table([low*u.keV,up*u.keV,energy_counts_adjusted],names=('lower_bin_edges','upper_bin_edges','counts'))

#running the HEXITEC simulation - see notation in Dan Ryan's code for futher details on options chosen
hexitec_sim = simulate_hexitec_on_spectrum(incident_spectrum,5000/u.s,num_photons,incident_xpixel_range=(4,9),incident_ypixel_range=(4,9),readout_xpixel_range=(5,8),readout_ypixel_range=(5,8),threshold=threshold)

#binning incident and measured photon data
hist_hex = np.histogram(hexitec_sim.measured_photons['energy'],np.arange(0,max_energy))
hist_in = np.histogram(hexitec_sim.incident_photons['energy'],np.arange(0,max_energy))

#calculating the percentage of events which display charge sharing
#estimating charge sharing event if more than one measured photon has the same subframe time
measured_photon_times, num_pixels_in_event = np.unique(hexitec_sim.measured_photons['subframe time'],return_counts=True)
hist_charge_share = np.histogram(num_pixels_in_event-1,np.arange(0,8))
fraction_charge_share = float(hist_charge_share[0].sum()-hist_charge_share[0][0])/hist_charge_share[0].sum()
fig10 = plt.figure()
plt.title('Number of charge sharing events for measured incident photons: Threshold = {0}'.format(threshold))
plt.hist(num_pixels_in_event-1,np.arange(0,8))
plt.xlabel('Number of charge sharing events')
plt.ylabel('Counts')
plt.show()

#plotting the comparison between the incident and the measured spectra
fig0 = plt.figure()
plt.title('Comparison between incident and HEXITEC readout spectra')
plt.plot(hist_hex[1][1:],hist_hex[0],label='HEXITEC readout',color='blue')
plt.plot(hist_in[1][1:],hist_in[0], label='Incident',color='red')
plt.xlabel('Energy /keV')
plt.ylabel('Counts')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.show()


fig1 = plt.figure()
plt.plot(incident_spectrum['lower_bin_edges'],incident_spectrum['counts'])
plt.title('Incident spectrum')
plt.xlabel('Energy /keV')
plt.ylabel('Counts')
plt.show()

fig2 = plt.figure()
plt.title('Comparison between incident and HEXITEC readout spectra')
plt.hist(hexitec_sim.incident_photons['energy'],energies_counts[:,0],label='Incident',color='red')
plt.hist(hexitec_sim.measured_photons['energy'],energies_counts[:,0],label='HEXITEC readout',color='blue')
plt.xlabel('Energy /keV')
plt.ylabel('Counts')
plt.legend()
plt.show()

