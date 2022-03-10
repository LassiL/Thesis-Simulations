"""I should be able to just define the input efield here as a function and then run all the scripts. I.e., just write all the other scripts as functions and solve them here. Do also the plotting here s.t. I can reproduce all the files in my thesis"""

# See the updated files in path: /Users/lassi/Desktop/all-thesis/scripts-for-thz-txrd/trxd.py

# just import ho_solve_odeint
# 1. Have to remove the fact that txrd.py reads the .txt files: 
#   - I can just call Uys[:,0], Uys[:,0], Uzs[:,0] which contain the displacements from ho_solve_odeint.py to lines 33-35 in txrd.py
#        - Easy fix as this is just unnecessary step
"""This is way harder"""
# 2. I need to write the harmonic oscillator equations as a function of the electric field components
# 3. My electric field function needs to return the electric field components, which are then taken as arguments in the harmonic oscillator
#   - How do I do this
# 4. Should define the time globally here as well (see the note on line 93 in ho_solve_odeint.py)

# I would like to be able to loop over different field parameters and plot similar graphs as the final graphs in my thesis
# def E_field(arguments ...):
#     User needs to input the polarization, total field strength, FWHM and the carrier frequency
#           - Maybe I just define some of them globally.
#     plot the time and fourier tf of the pulse.
#     return [array of the field components? and the harmonic oscillator equations take the elements of the arrays?]
#       - I.e., 
# the constant a and b in the harmonic oscillator functions need to call this pulse function, lines 100-101 etc


# Maybe should write as a class, where I have different pulses for different polarizations?
# class Person:
#   def __init__(self, name, age):
#     self.name = name
#     self.age = age

#   def myfunc(self):
#     print("Hello my name is " + self.name)

# p1 = Person("John", 36)
# p1.myfunc()

"""Order of business to be done"""
# 1. Remove that txrd reads the txt files 
#   - test if it works
# 2. Rewrite the harmonic oscillators functions as explained in lines 110 and 111 in ho_solve_odeint.py
#    - test that I can run as before when I define the electric field components globally
# 3. Write the electric field function in this script (I could just define the pulse here without writing it as a function)
#   - (do I really need to write it as a function? or what do I actually gain by this?)
#       - I could just do define E_0tot = [0,5,100]E8 and then for i in E_tot: calculate polarization and solve the Uxs
#   - how would I run the program? E_field(100,10,10) returns an array
#   - i.e., in my harmonic oscillator I would write  E_field(100,10,10)[0] for E_0x component
# 4. Incorporate the if name main testing in each script
# 5. Could I define the double gaussian pulse field? I would need to modify the force terms
#   - Maybe I should just modify the force terms to start with?