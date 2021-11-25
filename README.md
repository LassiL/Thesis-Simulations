# Thesis-Simulations
Python scripts to simulate the time-resolved x-ray diffraction intensity of nearly forbidden Bragg reflections of 
the InSb zincblende geometry from THz excited phonons.

1. pre_xrd.py script
    - calculates the atomic form factors needed for the xrd intensity
2. ho_solve_odeint.py 
    - solves driven damped harmonic oscillator equations for user defined THz laser pulse and anharmonic lifetime
3. txrd.py 
    - takes the harmonic oscillator solutions and calculates the txrd intensity of chosen Bragg reflections
