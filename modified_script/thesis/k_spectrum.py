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

coord       = gas_coord * scale_factor * unit_length_cgs
mass        = gas_mass * unit_mass_cgs
velocity    = gas_vel * np.sqrt(scale_factor) * unit_vel_cgs
density     = gas_density * unit_mass_cgs / (unit_length_cgs * scale_factor)**3
potential   = gas_potent * unit_vel_cgs**2 / scale_factor
internal    = gas_intern * unit_vel_cgs**2

# Particle Type 1
dm_coord    = f['PartType1']['Coordinates'][:] * scale_factor * unit_length_cgs
dm_mass     = f['PartType1']['Masses'][:] * unit_mass_cgs

time = 1/scale_factor - 1

ind = np.where( density > 200*cosmo.critical_density(time).value )[0]

coord_over  = coord[ind]
mass_over   = mass[ind]
density_over = density[ind]

COM_1 = Center_of_Quan( coord_over, mass_over )

# Repeat it again
rel_x = coord_over[:, 0] - COM_1[0]
rel_y = coord_over[:, 1] - COM_1[1]
rel_z = coord_over[:, 2] - COM_1[2]

distance = np.sqrt( rel_x**2 + rel_y**2 + rel_z**2)
r_offset = 5 * scale_factor * unit_length_cgs
ind = np.where(distance < r_offset)[0]

# if ind.size == 0:
#     mass_gt_0 = mass_over[np.where(coord_over[:, 1] > 0)[0]]
#     mass_lt_0 = mass_over[np.where(coord_over[:, 1] < 0)[0]]

#     if np.sum(mass_gt_0) > np.sum(mass_lt_0):
#         ind = np.where(coord_over[:, 1] > 0)[0]
#     else:
#         ind = np.where(coord_over[:, 1] < 0)[0]

COM_2 = Center_of_Quan( coord_over[ind], mass_over[ind] )

# Find the final position by finding the maximum density
rel_x = coord_over[:, 0] - COM_2[0]
rel_y = coord_over[:, 1] - COM_2[1]
rel_z = coord_over[:, 2] - COM_2[2]

distance = np.sqrt( rel_x**2 + rel_y**2 + rel_z**2)
r_offset = 2.5 * scale_factor * unit_length_cgs
ind = np.where(distance < r_offset)[0]
# COM = coord_over[np.argmax(density_over[ind])]
COM = Center_of_Quan( coord_over[ind], mass_over[ind] )

print("Center of Mass: ", COM/(scale_factor * unit_length_cgs))

# Position and velocity relative to center of mass
rel_x = coord[:, 0] - COM[0]
rel_y = coord[:, 1] - COM[1]
rel_z = coord[:, 2] - COM[2]

# radius within 1 ckpc/h
radius = 1.5
distance = np.sqrt( rel_x**2 + rel_y**2 + rel_z**2)
r_offset = radius * scale_factor * unit_length_cgs
ind_1 = np.where(distance < r_offset)[0]
print("Search radius: %.3f pc" % (r_offset/3.08e18))

coord_in_1  = coord[ind_1]
mass_in_1   = mass[ind_1]
vel_in_1    = velocity[ind_1]

# Calculate the total energy respect to the center of mass
distance = np.logspace(-2, 2, 200)

gas_distance = np.sqrt( np.sum( (coord - COM)**2, axis=1 ) )
dm_distance  = np.sqrt( np.sum( (dm_coord - COM)**2, axis=1 ) )

f = open("record.txt", "w")
r_offset = radius * scale_factor * unit_length_cgs / 3.08e18
f.write("Plot area: \n")
for i in range(3):
    com_i = COM[i] / 3.08e18
    f.write("%2.6e %2.6e\n" % (com_i - r_offset, com_i + r_offset))

f.write("Center of Mass: \n")
for i in range(3):
    f.write("%2.6e\n" % (COM[i] / 3.08e18))

# Plot the Kolmogorov Spectrum

box_len = [32, 8]
for i in box_len:

    len = i
    print("Box_length: %d pc" % (i*scale_factor*unit_length_cgs/3.08e18))
    CM = COM / (scale_factor * unit_length_cgs)

    bbox = [[-len + CM[0], len + CM[0]], [-len + CM[1], len + CM[1]], [-len + CM[2], len + CM[2]]]
    unit_base = {
        "length": (scale_factor * unit_length_cgs / 3.08e21, "kpc"),
        "velocity": (1.0, "km/s"),
        "mass": (unit_mass_cgs / 1.933e33, "Msun"),
    }

    ds = yt.load("./snapshot_%03d.hdf5" % (N_s), unit_base=unit_base, bounding_box=bbox)
    ds.force_periodicity()

    # a FFT operates on uniformly gridded data.  We'll use the yt
    # covering grid for this.

    max_level = 9

    low = ds.domain_left_edge
    dims = 2**max_level

    nx = ny = nz = dims

    nindex_rho = 1.0 / 3.0

    Kk = np.zeros((nx // 2 + 1, ny // 2 + 1, nz // 2 + 1))

    for vel in [("gas", "velocity_x"), ("gas", "velocity_y"), ("gas", "velocity_z")]:
        Kk += 0.5 * fft_comp(
            ds, ("gas", "density"), vel, nindex_rho, max_level, low, dims
        )

    # wavenumbers
    L = (ds.domain_right_edge - ds.domain_left_edge).d

    kx = np.fft.rfftfreq(nx) * nx / L[0]
    ky = np.fft.rfftfreq(ny) * ny / L[1]
    kz = np.fft.rfftfreq(nz) * nz / L[2]

    # physical limits to the wavenumbers
    kmin = np.min(1.0 / L)
    kmax = np.min(0.5 * dims / L)

    kbins = np.arange(kmin, kmax, kmin)
    N = np.shape(kbins)[0]

    # bin the Fourier KE into radial kbins
    kx3d, ky3d, kz3d = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx3d**2 + ky3d**2 + kz3d**2)

    whichbin = np.digitize(k.flat, kbins)
    ncount = np.bincount(whichbin)

    E_spectrum = np.zeros(np.shape(ncount)[0] - 1)

    for n in range(1, np.shape(ncount)[0]):
        E_spectrum[n - 1] = np.sum(Kk.flat[whichbin == n])

    k = 0.5 * (kbins[0 : N - 1] + kbins[1:N])
    E_spectrum = E_spectrum[1:N] / (kmin)**3

    index = np.argmax(E_spectrum[:-5])
    kmax = k[int(np.shape(k)[0] / 2)]
    Emax = E_spectrum[int(np.shape(k)[0] / 2)]

    # output the result to txt file
    np.savetxt("spectrum_%03d.txt" % len, np.column_stack((k, E_spectrum)))