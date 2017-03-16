import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import sys

h = 0.72
H0inv = 9.778e9*h
O_m0 = 0.3
w = 1.0 - O_m0

def find_age(h):
    """Method which returns t_0 and the value of H_0"""
    year = 60*60*24*365
    Mpc = 3.08567758e22
    km = 1e3    
    H_0 = 100*km*h/Mpc #[s-1]
    a_ml = (O_m0/w)**(1./3.)
    t_0 = 2./(3*H_0*np.sqrt(w))*np.arcsinh((1./a_ml)**(3./2.))
    print "The age of the universe with LambdaCDM parameters is %.2f billion years" %(t_0/year/1e9)
    return H_0, t_0
H0, t0 = find_age(h)

def find_particle_horizon_today(t0):
    lyTom = 9.4605284e15
    c = 3e8
    print "The proper distance to the particle horizon today with LambdaCDM parameters is %.2f Gly" %(c*t0/lyTom/1e9)
    return c*t0
find_particle_horizon_today(t0)

def find_event_horizon_today(t0):
    
