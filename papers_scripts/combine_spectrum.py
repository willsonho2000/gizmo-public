import numpy as np
import matplotlib.pyplot as plt
import h5py
import os

font = {'family': 'serif',
        'weight': 'normal',
        'size': 15,
}

snap_folder = './'
snap_list = [int(snap[9:12]) for snap in os.listdir(snap_folder) if (snap.startswith("snapshot_") & snap.endswith(".hdf5"))]

N_i = min(snap_list)
N_f = max(snap_list) + 1
# N_s = 54
N_s = N_f - 1

# Calculate the mach number of the specific snapshots

f = h5py.File( "snapshot_%03d.hdf5" % N_s )

scale_factor        = f['Header'].attrs['Time']
h_factor            = f['Header'].attrs['HubbleParam']
unit_length_cgs     = f['Header'].attrs['UnitLength_In_CGS']
unit_mass_cgs       = f['Header'].attrs['UnitMass_In_CGS']
unit_vel_cgs        = f['Header'].attrs['UnitVelocity_In_CGS']
time = 1/scale_factor - 1

# Find all spectrum files
spectrum_file = [f for f in os.listdir() if f.startswith("spectrum")]
spectrum_scale = [int(f.split("_")[1].split(".")[0]) for f in spectrum_file]
# sort by scale from small to large
spectrum_scale.sort(reverse=True)

# load data
data = np.loadtxt("spectrum_%03d.txt" % spectrum_scale[0])
k = data[:, 0]
E_spectrum = data[:, 1]
k /= (scale_factor * unit_length_cgs / 3.08e21)

offset = k[0]
scale = spectrum_scale[0]
for i in spectrum_scale[1:]:
    data = np.loadtxt("spectrum_%03d.txt" % i)
    k_new = data[:, 0] / (scale_factor * unit_length_cgs / 3.08e21)
    E_new = data[:, 1]
    
    old_offset = np.where(k > 50)[0][0]
    new_offset = np.where(k_new > 50)[0][0]

    # if np.shape(find_offset)[0] == 0:
    #     continue
    # else:
    #     find_offset = find_offset[-1]

    amplitude_ratio = E_new[new_offset] / E_spectrum[old_offset]

    k = np.concatenate( (k[:old_offset], k_new[new_offset:]) )
    E_spectrum = np.concatenate( (E_spectrum[:old_offset], E_new[new_offset:] / amplitude_ratio)) 
    offset = k[0]

k_off = np.where(k > 110)[0][0]
k = k[:k_off]
E_spectrum = E_spectrum[:k_off]

kmax = k[np.argmax(E_spectrum[2:]) + 2]
Emax = E_spectrum[np.argmax(E_spectrum[2:]) + 2]

plt.figure( figsize=(7, 6), dpi= 100 )

plt.axvline(kmax, color="r", linestyle="--", label=r"k$_1$ = %.2f pc" % (1e3/kmax))
plt.loglog(k[:k_off], E_spectrum[:k_off], c='b')
plt.loglog(k[:k_off], Emax * (k[:k_off] / kmax) ** (-5.0 / 3.0), ls=":", color="0.5", label=r"k$^{-5/3}$")
# plt.xlabel(r"k = 1/L [kpc$^{-1}$]", fontdict=font)
# plt.ylabel(r"E(k)", fontdict=font)
plt.xticks( fontsize=15 )
plt.yticks( fontsize=15 )

plt.legend(prop=font)