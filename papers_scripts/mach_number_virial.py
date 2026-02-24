import numpy as np
import h5py
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.constants import G
import os

# Constants
G = 6.67430e-11         # Gravitational constant in m^3 kg^-1 s^-2
G *= (100)**3 / 1000 # to cgs
Number_Density_Fraction = 1.16e-24

# Functions

def Virial_Radius( halo_mass, redshift):
    # Example values
    halo_mass = halo_mass * u.g  # Halo mass
    overdensity = 200  # Overdensity factor

    # Calculate critical density of the universe at given redshift
    critical_density = cosmo.critical_density(redshift)

    # Calculate the virial radius
    # M = (4/3) * pi * R_vir^3 * Delta * rho_c
    # Rearranging for R_vir:
    virial_radius = ((3 * halo_mass) / (4 * np.pi * overdensity * critical_density))**(1/3)

    return virial_radius.value

def Mass_Shells( coord, mass, center, num_bins ):
    # Calculate the distance from the center of mass
    dist = np.sqrt( np.sum( (coord - center)**2, axis=1 ) )

    # Calculate the virial radius
    max_radius = max( dist )

    # Calculate the bin edges
    bin_edges = np.linspace( 0, max_radius, num_bins + 1 )

    # Calculate the mass in each shell
    mass_in_shells, bin_edges = np.histogram( dist, bins=bin_edges, weights=mass )

    # Calculate the mass in each shell
    mass_in_shells = mass_in_shells * u.g

    # Calculate the volume of each shell
    shell_volume = (4/3) * np.pi * (bin_edges[1:]**3 - bin_edges[:-1]**3)

    # Calculate the density in each shell
    density_in_shells = mass_in_shells.value / shell_volume

    return density_in_shells, bin_edges

def Center_of_Quan( quan, mass ):
    # Calculate the center of mass
    coq = np.zeros(3)
    for i in range(3):
        coq[i] = np.sum( quan[:,i] * mass ) / np.sum( mass )

    return coq

def Potential( coord, mass ):
    # Calculate the potential energy of the system
    G = 6.67430e-11         # Gravitational constant in m^3 kg^-1 s^-2
    G = G * (100)**3 / 1000 # to cgs

    potential = 0
    for i in range(len(mass)):
        for j in range(i+1, len(mass)):
            r = np.sqrt( np.sum( (coord[i] - coord[j])**2 ) )
            potential += - G * mass[i] * mass[j] / r
    
    return potential

def Density_In_Beams( mass, density, num_bins, to_number_density=False ):

    if to_number_density:
        density_cal = density / Number_Density_Fraction
        bin_edges = np.logspace( -6, 9, num_bins + 1 )
    else:
        density_cal = density
        bin_edges = np.logspace( -30, -13, num_bins + 1 )

    # Calculate the mass in each beams
    density_in_beams, bin_edges = np.histogram( density_cal, bins=bin_edges, weights=mass )
    density_in_beams_fraction = density_in_beams / np.sum(mass)

    return density_in_beams_fraction, bin_edges

def Read_File( file ):
    # Open the file
    # file = h5py.File(file, 'r')
    
    # Read the file
    scale_factor        = file['Header'].attrs['Time']
    unit_length_cgs     = file['Header'].attrs['UnitLength_In_CGS']
    unit_mass_cgs       = file['Header'].attrs['UnitMass_In_CGS']
    unit_vel_cgs        = file['Header'].attrs['UnitVelocity_In_CGS']

    # Calculate the position and the velocity of the center of mass

    gas_coord   = file['PartType0']['Coordinates'][:]
    gas_vel     = file['PartType0']['Velocities'][:]
    gas_mass    = file['PartType0']['Masses'][:]
    gas_density = file['PartType0']['Density'][:]
    gas_temp    = file['PartType0']['Temperature'][:]

    coord   = gas_coord * scale_factor * unit_length_cgs
    vel     = gas_vel * np.sqrt(scale_factor) * unit_vel_cgs
    mass    = gas_mass * unit_mass_cgs
    density = gas_density * unit_mass_cgs / (unit_length_cgs * scale_factor)**3

    # Close the file

    return coord, mass, density, vel, gas_temp

def Speed_of_Sound( temperature, ele_abund ):
    # Constants for air
    gamma = 5/3  # Adiabatic index for idea gas
    R = 8.31  # Specific gas constant for ISM in J/mol·K
    m = 1.67e-3  # Molar mass of ISM in kg/mol
    X_H = 0.76  # Hydrogen mass fraction
    mu = 4 / (1 + 3*X_H + 4*X_H*ele_abund)  # Mean molecular weight of ISM

    # Calculate the speed of sound in the medium
    speed_of_sound = np.sqrt(gamma * R * temperature / (mu * m)) * 100 # in cm/s

    return speed_of_sound

def Mach_Number( velocities, masses, temperature, elec_abund ):
    # Calculate the average velocity
    # If you have masses and want to do a mass-weighted average:
    total_mass = np.sum(masses)
    average_velocity = np.sum(velocities.T * masses, axis=1) / total_mass

    # Calculate the velocity fluctuations
    velocity_fluctuations = velocities - average_velocity

    # If needed, calculate the magnitude of these fluctuations
    velocity_fluctuation_magnitudes = np.linalg.norm(velocity_fluctuations, axis=1)

    # Calculate the speed of sound in the medium
    speed_of_sound = Speed_of_Sound( temperature, elec_abund )

    # Calculate the Mach number for each particle
    mach_numbers = velocity_fluctuation_magnitudes / speed_of_sound

    return mach_numbers, velocity_fluctuation_magnitudes
    

def fft_comp(ds, irho, iu, nindex_rho, level, low, delta):
    cube = ds.covering_grid(level, left_edge=low, dims=delta, fields=[irho, iu])

    rho = cube[irho].d
    u = cube[iu].d

    nx, ny, nz = rho.shape

    # do the FFTs -- note that since our data is real, there will be
    # too much information here.  fftn puts the positive freq terms in
    # the first half of the axes -- that's what we keep.  Our
    # normalization has an '8' to account for this clipping to one
    # octant.
    ru = np.fft.fftn(rho**nindex_rho * u)[
        0 : nx // 2 + 1, 0 : ny // 2 + 1, 0 : nz // 2 + 1
    ]
    ru = 8.0 * ru / (nx * ny * nz)

    return np.abs(ru) ** 2

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

# Particle Type 0
gas_coord   = f['PartType0']['Coordinates'][:] 
gas_vel     = f['PartType0']['Velocities'][:]
gas_mass    = f['PartType0']['Masses'][:]
gas_density = f['PartType0']['Density'][:]
gas_potent  = f['PartType0']['Potential'][:]
gas_temp    = f['PartType0']['Temperature'][:]
gas_intern  = f['PartType0']['InternalEnergy'][:]
gas_electron= f['PartType0']['ElectronAbundance'][:]
gas_smooth  = f['PartType0']['SmoothingLength'][:]

coord       = gas_coord * scale_factor * unit_length_cgs
mass        = gas_mass * unit_mass_cgs
velocity    = gas_vel * np.sqrt(scale_factor) * unit_vel_cgs
density     = gas_density * unit_mass_cgs / (unit_length_cgs * scale_factor)**3
potential   = gas_potent * unit_vel_cgs**2 / scale_factor
internal    = gas_intern * unit_vel_cgs**2
smoothing   = gas_smooth * scale_factor * unit_length_cgs

# Particle Type 1
dm_coord    = f['PartType1']['Coordinates'][:] * scale_factor * unit_length_cgs
dm_mass     = f['PartType1']['Masses'][:] * unit_mass_cgs

time = 1/scale_factor - 1

# Alternitive way to calculate the center of mass (The shortest softening length)
COM_2 = coord[np.argmin( gas_smooth )]

# Find the final position by finding the maximum density
rel_x = coord[:, 0] - COM_2[0]
rel_y = coord[:, 1] - COM_2[1]
rel_z = coord[:, 2] - COM_2[2]

distance = np.sqrt( rel_x**2 + rel_y**2 + rel_z**2)
r_offset = 300 * 3.08e18
ind = np.where(distance < r_offset)[0]
COM = coord[ind][np.argmax(density[ind])]
# COM = Center_of_Quan( coord_over[ind], mass_over[ind] )

# Position and velocity relative to center of mass
rel_x = coord[:, 0] - COM[0]
rel_y = coord[:, 1] - COM[1]
rel_z = coord[:, 2] - COM[2]

# Calculate the total energy respect to the center of mass
distance = np.logspace(-2, 2, 200)

gas_distance = np.sqrt( np.sum( (coord - COM)**2, axis=1 ) )
dm_distance  = np.sqrt( np.sum( (dm_coord - COM)**2, axis=1 ) )

total_gas   = []
total_dm    = []

# Do at the very first
gas_ind = np.where( gas_distance < distance[0] * scale_factor * unit_length_cgs )[0]
dm_ind  = np.where( dm_distance < distance[0] * scale_factor * unit_length_cgs )[0]
total_gas.append(np.sum( mass[gas_ind] ))
total_dm.append(np.sum( dm_mass[dm_ind] ))

for d in distance[1:]:
    d *= (scale_factor * unit_length_cgs)
    # Find the gas and dm particles within the distance
    gas_ind = np.where( gas_distance < d )[0]
    dm_ind  = np.where( dm_distance < d )[0]

    total_gas.append(np.sum( mass[gas_ind] ))
    total_dm.append(np.sum( dm_mass[dm_ind] ))


# Compute the radius of the halo

total_mass = np.array(total_gas) + np.array(total_dm)
density_distance = total_mass / (4/3 * np.pi * (distance*scale_factor*unit_length_cgs)**3)

# Find the point where the density drops below 200 times the critical density
ind = np.where( density_distance > 200*cosmo.critical_density(time).value )[0]

if len(ind) == 0:
    R_vir = distance[ np.argmax(density_distance) ] * scale_factor * unit_length_cgs
    M_vir = total_mass[ np.argmax(density_distance) ] / unit_mass_cgs
    M_vir_gas = total_gas[ np.argmax(density_distance) ] / unit_mass_cgs
    M_vir_dm = total_dm[ np.argmax(density_distance) ] / unit_mass_cgs
else:
    R_vir = distance[ind[-1]] * scale_factor * unit_length_cgs
    M_vir = total_mass[ind[-1]] / unit_mass_cgs
    M_vir_gas = total_gas[ind[-1]] / unit_mass_cgs
    M_vir_dm = total_dm[ind[-1]] / unit_mass_cgs

# Calculate the Mach number within the virial radius

ind = np.where( gas_distance < R_vir )[0]

coord_in_1  = coord[ind]
mass_in_1   = mass[ind]
vel_in_1    = velocity[ind]
temp_in_1   = gas_temp[ind]
ele_abund_in_1 = gas_electron[ind]
density_in_1 = density[ind]

mach_numbers, velocity_fluctuation_magnitudes = Mach_Number( vel_in_1, mass_in_1, temp_in_1, ele_abund_in_1 )
mean_mach = np.sum(mach_numbers * mass_in_1) / np.sum(mass_in_1)
std_mach = np.sqrt( np.sum( (mach_numbers - mean_mach)**2 * mass_in_1 ) / np.sum(mass_in_1) )

# Compute velocity dispersion (simple finite difference approximation)
f = open("mach_number_%03d_radius.txt" % N_s, "w")
f.write("Mean,\tStd,\tVirial Radius [pc],\tVirial Mass [Msun],\tGas Mass[Msun],\tDM Mass[Msun],\tz\n")
f.write("%.4f,\t%.4f,\t%2.6e,\t\t%2.6e,\t\t%2.6e,\t\t%2.6e,\t%.2f\n" % 
        (mean_mach, std_mach, R_vir/3.08e18,M_vir * 1e10/h_factor, M_vir_gas * 1e10/h_factor, M_vir_dm * 1e10/h_factor, time))
f.close()

# Save the mach number curve
data = np.histogram( mach_numbers, bins=100, weights=mass_in_1/unit_mass_cgs*1e10/h_factor )
f = open("mach_number_%03d_curve.txt" % N_s, "w")
f.write("Mach Number,\tMass [Msun]\n")

for i in range(len(data[0])):
    f.write("%.4f,\t%.4f\n" % (data[1][i], data[0][i]))

f.close()

# Save the mach number curve for density
data = np.histogram( density_in_1, bins=np.logspace(-27,-12, 100), weights=mass_in_1/unit_mass_cgs*1e10/h_factor )
f = open("mach_number_%03d_curve_density.txt" % N_s, "w")
f.write("Density,\tMass [Msun]\n")

for i in range(len(data[0])):
    f.write("%.4e,\t%.4f\n" % (data[1][i], data[0][i]))

# f = open("record.txt", "w")
# Center = COM / (scale_factor * unit_length_cgs)
# f.write("Center [code length] \n")
# f.write("%2.6e,\t%2.6e,\t%2.6e\n" % (Center[0], Center[1], Center[2]) )
# f.write("Virial Radius [code length] \n")
# f.write("%2.6e\n" % (R_vir/(scale_factor * unit_length_cgs)))
f.close()
