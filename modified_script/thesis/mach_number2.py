import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo
from astropy.constants import G
import yt

font = {'family': 'serif',
        'weight': 'normal',
        'size': 10,
}

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

def Mach_Number( velocities, masses, temperature ):
    # Calculate the average velocity
    # If you have masses and want to do a mass-weighted average:
    total_mass = np.sum(masses)
    average_velocity = np.sum(velocities.T * masses, axis=1) / total_mass

    # Calculate the velocity fluctuations
    velocity_fluctuations = velocities - average_velocity

    # If needed, calculate the magnitude of these fluctuations
    velocity_fluctuation_magnitudes = np.linalg.norm(velocity_fluctuations, axis=1)

    # Constants for air
    gamma = 5/3  # Adiabatic index for idea gas
    R = 8.31  # Specific gas constant for ISM in J/mol·K
    m = 1e-3  # Molar mass of ISM in g/mol

    # Calculate the speed of sound in the medium
    speed_of_sound = np.sqrt(gamma * R * temperature / m) * 100

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

def Speed_of_Sound( temperature, ele_abund ):
    gamma = 5 / 3
    k_B = 1.380649e-23  # Boltzmann constant in J/K
    m_p = 1.6726219e-27  # Proton mass in kg
    X_H = 0.76  # Hydrogen mass fraction
    mu = 4 / (1 + 3*X_H + 4*X_H*ele_abund)  # Mean molecular weight of ISM

    # Compute speed of sound
    speed_of_sound = np.sqrt(gamma * k_B * temperature / (mu * m_p)) * 100

    return speed_of_sound

N_s = 30

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

# Compute velocity dispersion (simple finite difference approximation)
from scipy.spatial import KDTree
tree = KDTree(coord)    # Create KDTree for fast neighbor search

# Create density beams from 1e-22 to 1e-13
beams = np.logspace(-26, -15, 100)

f = open("mach_number_%03d_2.txt" % N_s, "w")
f.write("Density,Mean,Std\n")

for j in range(len(beams) - 1):
    ind = np.where( (density > beams[j]) & (density < beams[j+1]) )[0]
    if len(ind) == 0:
        continue

    velocities = velocity[ind]
    positions = coord[ind]
    masses = mass[ind]
    temperature = gas_temp[ind]
    ele_abund = gas_electron[ind]

    speed_of_sound = Speed_of_Sound( temperature, ele_abund )
    velocity_disper = np.zeros(np.shape(velocities)[0])

    for i in range(np.shape(velocities)[0]):
        distances, indices = tree.query(positions[i], 32)  # Find neighbors
        neighbors = velocity[indices[1:]]  # Exclude the point itself
        local_velocities = velocities[i] - neighbors
        squared_deviations = np.sum(local_velocities**2, axis=1)
        velocity_disper[i] = np.sqrt(np.mean(squared_deviations))

    mach_number_of_disperse = velocity_disper / speed_of_sound
    mean_value = np.mean( mach_number_of_disperse )
    std_value = np.std( mach_number_of_disperse )

    f.write("%.3e,%.3f,%.3f\n" % (beams[j], mean_value, std_value))

f.close()