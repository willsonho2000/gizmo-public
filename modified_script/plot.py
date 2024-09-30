import numpy as np
import matplotlib.pyplot as plt
import h5py

font = {'family': 'serif',
#         'color':  'dark',
        'weight': 'normal',
        'size': 10,
        }

for i in range(102):
    f1 = h5py.File("kh_mcnally_ultra/snapshot_%03d.hdf5" % (i),'r')
    
    x, y = f1['PartType0']['Coordinates'][:,0], f1['PartType0']['Coordinates'][:,1]
    z = np.log10(np.array(f1['PartType0']['Density']))

    plt.figure(figsize=(7,5))
    ax = plt.axes()
    ax.set_facecolor((0.0, 0.0, 0.5, 1.0))
    plt.scatter(x, y, c=z, s=0.06, vmin=-0.05, vmax=0.4, cmap="jet")
    cbar = plt.colorbar()
    font['size'] = 10
    cbar.set_label('log density (g/cm$^3$)', rotation=270, labelpad=14, fontdict=font)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(0.1, 0.9, "%.2f sec." % (i / 10), bbox=props)

    # plt.xlim(-5, 5)
    plt.xlabel("cm", fontdict=font)

    # plt.ylim(-5, 5)
    plt.ylabel("cm", fontdict=font)

    font['size'] = 16
    plt.title("Density", fontdict=font)

#     plt.show()
    plt.savefig("kh_mcnally_ultra/%03d.png" % (i))
    plt.close()
