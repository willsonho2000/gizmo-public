import h5py
import numpy as np

f = h5py.File( "snapshot_007.hdf5", "r+" )

metalliticy = np.zeros( len(f['PartType0/Velocities']) ) + 1e-5
f['PartType0'].create_dataset('Metallicity', data=metalliticy)

f.close()
