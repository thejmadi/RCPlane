# -*- coding: utf-8 -*-
"""
Created on Fri Aug  5 12:45:35 2022

@author: tarun
"""

# Input desired coefficient of lift, pitch moment coefficient, scaling constant,
# chord length, maximum thickness. Program will calculate camber line using cubic
# polynomial. Program will then calculate upper and lower edges using y_t. Program
# will then plot all three lines, and print the coefficient of lift and
# pitch moment coefficient. Program will then output upper and lower edge coordinates
# to csv file. 

# Uses NACA y_t = 5*t*((0.2969*np.sqrt(x)) - (0.126*x) - (0.3516*(x**2)) + (0.2843*(x**3)) - (0.1015*(x**4)))
# for upper and lower edges of airfoil

import numpy as np
import matplotlib.pyplot as plt


############## VARIABLES ###############


# Variables for camber line adjustment

a = 1 # DO NOT CHANGE
c = 1 # DO NOT CHANGE


angle_of_attack = 9 # Used in final sectional lift coefficient, NOT in airfoil shape calculation
c_l = .7 # desired coefficient of lift
c_m_quarter = 0.0 # desired pitch moment coefficient


A = np.array([[0.0, 0.0, 0.0, 1.0],
             [1.0, 1.0, 1.0, 1.0],
             [-3*np.pi/2, -np.pi, 0.0, 0.0],
             [3*np.pi/8, 0.0, 0.0, 0.0]])


B = np.array([0.0 / a, 0.0 / a, c_l / a, ((4 * c_m_quarter) + c_l) / a])

# Variables for thickness adjustment

t = .15 # max thickness in tenths of chord length


################### FUNCTIONS ###############


# Function to calc theta to be used in calc of coords of upper and lower edges
def thetaCalc(x, coeff):
    theta = []
    # theta = arctan(dy_camber/dx)
    theta = np.arctan(a*((3*coeff[0]*(x[:]**2)) + (2*coeff[1]*x[:]) + coeff[2])) # took out a*c
    return theta
    
# Function to calc coords of camber line, upper and lower edges. 
def coordCalc(x, coeff, theta):
    y_camber = []
    y_upper = []
    x_upper = []
    y_lower = []
    x_lower = []
    
    # y = ac[B_0(x/c)^3 + B_1(x/c)^2 + B_2(x/c) + B_3]
    y_camber = a*(coeff[0]*(x[:]**3) + coeff[1]*(x[:]**2) + coeff[2]*(x[:]) + coeff[3]) #took out a*c
    
    y_t = []
    y_t = 5*t*((0.2969*np.sqrt(x)) - (0.126*x) - (0.3516*(x**2)) + (0.2843*(x**3)) - (0.1015*(x**4)))
    
    # For loop to calculate upper and lower edges
    for i in range(0, len(y_camber)):
         y_upper.append(y_camber[i] + (y_t[i] * np.cos(theta[i])))
         y_lower.append(y_camber[i] - (y_t[i] * np.cos(theta[i])))
         x_upper.append(x[i] - (y_t[i] * np.sin(theta[i])))
         x_lower.append(x[i] + (y_t[i] * np.sin(theta[i])))
    
    return y_camber, x_upper, y_upper, x_lower, y_lower

# Calculates the sectional lift at a given angle of attack
def coeff_lift(AoA_list, coeff, angle_of_attack):
    A_0 = 0
    A_1 = a * (-(3 * coeff[0] / 2) - coeff[1])
    
    sectional_c_l = []
    AoA_list_radians = np.radians(AoA_list)
    
    for aoa in AoA_list_radians:
        A_0 = aoa - (a * ((9 * coeff[0] / 8) + coeff[1] + coeff[2]))
        sectional_c_l.append(2*np.pi*A_0 + np.pi*A_1) # c_l = 2*pi*A_0 + pi*A_1
        
    return sectional_c_l
    
################## MAIN ################


# Solves for coefficients of camber line cubic polynomial
coeff = np.linalg.solve(A, B)

# Creates x,y coord lists for plotting, and theta list for y coord calc
x_camber = np.linspace(0, c, 400)
theta_vals = thetaCalc(x_camber, coeff)
y_camber, x_upper, y_upper, x_lower,y_lower = coordCalc(x_camber, coeff, theta_vals)

# Generates lists of angles of attack and c_l to be plotted using coeff_lift() function
angle_of_attack_list = np.linspace(-15, 15, 400)
sectional_c_l_list = coeff_lift(angle_of_attack_list, coeff, angle_of_attack)

################## Plotting ##################

# Plots airfoil
plt.plot(x_camber, y_camber)
plt.plot(x_upper, y_upper)
plt.plot(x_lower, y_lower)
plt.ylim(-1.0, 1.0)
plt.show()

# Plots c_l vs. alpha graph
plt.plot(angle_of_attack_list, sectional_c_l_list)
plt.xticks(np.arange(angle_of_attack_list[0], angle_of_attack_list[-1], 2))
plt.yticks(np.arange(sectional_c_l_list[0], sectional_c_l_list[-1], 0.2))
plt.grid()
plt.show()

# Necessary code to write airfoil coords into text file
x_lower.reverse()
y_lower.reverse()

x_airfoil = x_upper + x_lower
y_airfoil = y_upper + y_lower

# Writes airfoil coords into text file
with open('wing_coords_proto_5.txt', 'w') as file:
    for i in range(0, len(x_airfoil)):
        file.writelines('%f\t%f\t%f\n' % (x_airfoil[i], y_airfoil[i], 0.0))

with open('wing_camberline_proto_5.txt', 'w') as file:
    for i in range(0, len(x_camber)):
        file.writelines('%f\t%f\t%f\n' % (x_camber[i], y_camber[i], 0.0))
