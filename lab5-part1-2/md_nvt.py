#!/usr/bin/env python

import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *
import tsase
import ase

##### Setting up the atoms object with ASE and TSASE #########
##############################################################

#### Import the starting structure from a file ###############
p = tsase.io.read_con('cluster_38.con')

#### Define the PES ##########################################
lj = tsase.calculators.lj(cutoff=35.0)
p.center(50.0)
p.set_calculator(lj)

###############################################################
#### MOLECULAR DYNAMICS #######################################
###############################################################

#### a few important functions for MD  #####
#### VELOCITY-VERLET ALGORITHM ##############
def step(p,dt,f,fixcm=True):
    m = p.get_momenta()
    m += 0.5 * dt * f    
    if fixcm:
        msum = m.sum(axis=0) / float(len(m))
        m = m - msum
    p.set_positions(p.get_positions() + dt * m / p.get_masses()[:,numpy.newaxis])
    p.set_momenta(m)
    f = p.get_forces()
    p.set_momenta(p.get_momenta() + 0.5 * dt * f)
    return f

################################################################


def run_md(numdt, dt): 
    f = p.get_forces()
    tsase.io.write_con('movie_nvt.con',p,w='w')  # initiates a movie for the MD trajectory 
    for i in range(int(time)):
        therm.apply_thermostat()                 # apply the thermostat at each MD step
        f = step(p,dt,f)                         # take an MD step 
        if i%100 == 0:                           ### get snap shot for movie every 100 timesteps 
            tsase.io.write_con('movie_nvt.con',p,w='a')   # Append to movie of the MD trajectory

##### specify MD parameters (Note: In order to convert time into 
##### fs you need to multiple by ase.units.fs
T = 300.                     # Set the initial temperature
dt_fs = 0.5                  # Set the initial time step for the MD simulations
dt = dt_fs * ase.units.fs    # Convert time step to correct units
kT = T * ase.units.kB        # Convert temperature to units of energy 


#### this is how to call the thermostat ########################
alpha = 0.8                  # Set Alpha for MD simulation
tcol = 100. * ase.units.fs   # Set tcol for MD simulations
therm = tsase.md.nvtandersen(p,dt,kT,alpha,tcol,fixcm=False)  # Call the thermostat from MD simulation

##### time of MD trajecotory in fs #############################
time = 50000.
numdt = time/dt_fs          # Set total time of the MD trajectory (Note: This trajectory is for 50,000 femtoseconds and will work for any timestep)

run_md(numdt,dt)


sys.exit()








