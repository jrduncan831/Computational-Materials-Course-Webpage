#!/usr/bin/env python

### importing important libraries 
import numpy
import sys
import matplotlib
matplotlib.use("agg") # Change to "agg" to run on FRI
from pylab import *

###### Lennard-Jones Potential ################
###### The function below returns the potential energy for a given bond length (r), sigma and epsilon
def LJ(r,sigma,epsilon):
    Vlj = 4* epsilon * ( numpy.power(sigma/r,12) - numpy.power(sigma/r,6) )
    return Vlj

##### This function creates an array of various bond lengths (rarray) 
##### and its corresponding energy (energyarray)
##### this is required to plot the potential
def get_energy_array(rmin,rmax,numbersamples,sigma,epsilon):
    rarray = numpy.linspace(rmin,rmax,numbersamples)
    engarray = numpy.zeros(int(numbersamples))
    for i in range(len(rarray)):
        engarray[i] = LJ(rarray[i],sigma,epsilon)
    return engarray,rarray ## Note: python can return as many variables as you would like from a functions
                           ## Just seperate the output by a comma 

    
# The function below plots how the potential energy changes with bondlength using matplotlib 
# Try commenting out the plot and scatter functions below to see how they are different! 
def plot_pes(energyar,rar,filename):
    figure()
    plot(rar,energyar)     #  This will plot a line connecting all points
    scatter(rar,energyar)  #  This will plot the points sampled
    ylim([-5,5])                 #  This specifies the range of the y-axis in your plot; adjust as needed
    xlabel('Bond Length')
    ylabel('Potential Energy')
    title('Lennard-Jones potential for diatomic')
    savefig(filename)      # This is how you generate an output file which stores the image you are creating.
#  filename should be a string with an extension like '.png' or '.eps' (image file extensions) 

# Below will be the main part of the function where you run the functions required to do your computations
# The script is currently setup to create your energy and bondlengths arrays and 
# then plot the PES with an output file 'LJ_s1_e1.png

# reminder of variables: rmin,rmax,numbersamples,sigma,epsilon
energyarray,rarray = get_energy_array(0.5, 10., 300.0, 1., 1.)  ## Note: since this function returns two values 
                                                                ## there are two variable names for each returned value.
plot_pes(energyarray,rarray,'LJ_s1_e1.png')



