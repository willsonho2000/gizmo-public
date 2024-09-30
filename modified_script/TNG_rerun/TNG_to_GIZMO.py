import h5py
import numpy as np
import os

file_name = os.listdir()

for i in file_name:
    if i.endswith('.hdf5'):
        file_name = i.split('.hdf5')[0]
        break

f  = h5py.File( file_name + '.hdf5' )
f2 = h5py.File( "TNG_modified.hdf5", 'w' )
f1 = h5py.File( "/tiara/home/myho/gizmo-public/music/level9.hdf5" )

tng_header      = list(  f['Header'].attrs.keys() )
target_header   = list( f1['Header'].attrs.keys() )

# Select the attributions wants to copy and transform

field_copy = []
for i in tng_header:
    if i in target_header: field_copy.append(i)
field_transform_f = ['UnitLength_in_cm', 'UnitMass_in_g', 'UnitVelocity_in_cm_per_s', 'NumPart_ThisFile', 'Omega0', 'OmegaBaryon', 'OmegaLambda']
field_transform_t = ['UnitLength_In_CGS', 'UnitMass_In_CGS', 'UnitVelocity_In_CGS', 'NumPart_Total', 'Omega_Matter', 'Omega_Baryon', 'Omega_Lambda']

# field_copy
f1.copy( 'Header', f2)
f2.create_group( name='PartType0' )
f2.create_group( name='PartType1' )

# Function of copy attributes
def attr_copy( file_to, file_1, ind_to, ind_1, hubble_scale=False ):
    file_to['Header'].attrs[ ind_to ] = file_1['Header'].attrs[ ind_1 ]

    if hubble_scale:
        file_to['Header'].attrs[ ind_to ] /=  file_1['Header'].attrs['HubbleParam']

# Copy the value of attributes from tng file
for i, ind in enumerate( field_copy ):
    if ind == 'MassTable':
        attr_copy( f2, f1, ind, ind)
    elif ind == "NumFilesPerSnapshot":
        f2['Header'].attrs.__delitem__(ind)
        f2['Header'].attrs.create(ind, 1, dtype='int32')
    else:
        attr_copy( f2, f, ind, ind)

# Transfrom the specific attributes
for i, ind in enumerate( field_transform_t ):
    f_ind = field_transform_f[i]   # Correspondings from the tng

    # set the case
    match ind:
        case 'UnitLength_In_CGS':
            attr_copy( f2, f, ind, f_ind, hubble_scale=True )
        case 'UnitMass_In_CGS':
            attr_copy( f2, f, ind, f_ind, hubble_scale=True )
        case 'NumPart_Total':
            attr_copy( f2, f, ind, f_ind )
        case 'Omega_Matter':
            attr_copy( f2, f, ind, f_ind )
        case 'Omega_Baryon':
            attr_copy( f2, f, ind, f_ind )
        case 'Omega_Lambda':
            attr_copy( f2, f, ind, f_ind )

# First calculate the range of x, y, z
min_set = []
max_set = []
pos_com = []
vel_com = []

for i in range(3):
    min_set.append( min( np.min( f['PartType0/Coordinates'][:, i] ) , np.min( f['PartType1/Coordinates'][:, i] ) ) )
    max_set.append( max( np.max( f['PartType0/Coordinates'][:, i] ) , np.max( f['PartType1/Coordinates'][:, i] ) ) )
    pos_com.append(np.average( f['PartType0/Coordinates'][:, i] ))
    vel_com.append( np.average( f['PartType0/Velocities'][:, i] ) )

# print( "center:", pos_com )
# print( "velocity:", vel_com )
print( "particle range: ", np.max( np.array( max_set ) - np.array( min_set ) ) )

if np.max( np.array( max_set ) - np.array( min_set ) ) < 30:
    f2['Header'].attrs['BoxSize'] = 30.0
else:
    f2['Header'].attrs['BoxSize'] = np.max( np.array( max_set ) - np.array( min_set ) ) * 1.1

## Copy the data from TNG

# For PartType0
field_copy = list(  f['PartType0'].keys() )
field_gizm = list( f1['PartType0'].keys() )

for vi, ind in enumerate( field_copy ):
    if ind in field_gizm:
        match ind:
            case 'ParticleIDs':
                gas_IDs = len( f['PartType0'][ind] )
                f2['PartType0'].create_dataset( name=ind, data=np.arange(0, gas_IDs, 1), dtype='u4' )
            case 'Coordinates':
                f_data = f['PartType0'][ind][:]
                f2['PartType0'].create_dataset( name=ind, data=f_data, dtype='f8' )
            case _:
                f2['PartType0'].create_dataset( name=ind, data=f['PartType0'][ind] )

field_more = ['ParticleChildIDsNumber', 'ParticleIDGenerationNumber']

# Create additional datasets
gas_len  = len( f['PartType0/Coordinates'] )
gas_data = np.zeros( gas_len )
for ind in field_more:
    f2['PartType0'].create_dataset( name=ind, data=gas_data, dtype='u4' )

# For PartType1
field_copy = list(  f['PartType1'].keys() )
field_gizm = list( f1['PartType1'].keys() )

for vi, ind in enumerate( field_copy ):
    if ind in field_gizm:
        match ind:
            case 'ParticleIDs':
                dm_IDs = len( f['PartType1'][ind] )
                f2['PartType1'].create_dataset( name=ind, data=np.arange(gas_IDs, gas_IDs+dm_IDs, 1), dtype='u4' )
            case 'Coordinates':
                f_data = f['PartType1'][ind][:]
                f2['PartType1'].create_dataset( name=ind, data=f_data, dtype='f8' )
            case _:
                f2['PartType1'].create_dataset( name=ind, data=f['PartType1'][ind] )

# Create DM masses' dataset
dm_len  = len( f['PartType1/Coordinates'] )
dm_data = np.ones( dm_len ) * 4.53e5 / 1e10
f2['PartType1'].create_dataset( name='Masses', data=dm_data, dtype='f4' )

# Create additional datasets
dm_data = np.zeros( dm_len )
for ind in field_more:
    f2['PartType1'].create_dataset( name=ind, data=dm_data, dtype='u4' )

# Modify the velocity and the position of the original data
part = [0, 1]

for i in part:
    f2['PartType' + str(i) ]['Coordinates'][:]   -= pos_com
    f2['PartType' + str(i) ]['Velocities'][:]    -= vel_com

unit_mass_cgs = f2['Header'].attrs['UnitMass_In_CGS']
unit_vel_cgs = f2['Header'].attrs['UnitVelocity_In_CGS']

gas_intern = f2['PartType0']['InternalEnergy'][:]
ele_abund  = f2['PartType0']['ElectronAbundance'][:]
f2['PartType0']['StarFormationRate'][:] = 0.0

# Calculate the temperature of gas particles
gamma = 5/3
k = 1.3807e-16
X_H = 0.76
m_p = 1.67e-24
mu = 4 / (1 + 3*X_H + 4*X_H*ele_abund) * m_p

T = (gamma-1) * gas_intern/k * pow(unit_vel_cgs, 2) * mu

f2['PartType0'].create_dataset('Temperature', data=T)

metalliticy = np.zeros_like( T ) + 1e-6
f2['PartType0'].create_dataset('Metallicity', data=metalliticy, dtype='f8')

h2_fraction = np.zeros_like( T )
f2['PartType0'].create_dataset('MolecularMassFraction', data=h2_fraction, dtype='f8')

f.close()
f1.close()
f2.close()