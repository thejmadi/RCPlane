import numpy as np
import matplotlib.pyplot as plt

#######     VARIABLES    ###########

N = 100                                                                          # Number of points over both sides of wing to be evaluated; 
root_chord = .0889
b = .7112                                                                       # Wingspan

c_zero_over_b = root_chord / b                                                           # Root chord over wingspan
wing_tip_length = 0.00000000001                                                           # Difference between root and tip chord
alpha_L_equals_zero = -5 * np.pi / 180                                           # AoA when Lift = 0
angle_of_attack = 0 * np.pi / 180 

total_alpha = angle_of_attack - alpha_L_equals_zero

AR = b / root_chord                                                             # Wing aspect ratio
S = b * root_chord                                                              # Wing area 


q = 0.5 * (19**2) * (1.144)                                                      # Dynamic Pressure = 0.5 * (freestream velocity ** 2) * (freestream air density)

B_matrix = np.zeros(shape = (N,N))
alpha_coeff_matrix = np.ones(N)
alpha_coeff_matrix.transpose()

aircraft_mass = 0.55                                                             # Total mass of plane kg
aircraft_weight = aircraft_mass * 9.81                                          # Total weight of plane N


#########     FUNCTIONS     #########

def thetaM(N):                                                                  # Calc theta_m array
    terms = []
    div = np.pi / (2 * N)
    for m in range(1, N + 1):
        terms.append(div * ((2 * m) - 1))
    return terms

def ymOverB(theta_m):
    return 0.5 * np.cos(theta_m[:])

def cmOverB(c_zero_over_b, theta_m, wing_tip_length):                           # Calc cm/b for linear taper
    terms = 1 - (abs(np.cos(theta_m[:])) * wing_tip_length)                     # Taper of 1/2
    return terms * c_zero_over_b

def B(N, c_zero_over_b, wing_tip_length, theta_m):                                       # Build B Matrix
    theta_m_arr = np.array(theta_m)
    cm_over_b = np.array(cmOverB(c_zero_over_b, theta_m_arr, wing_tip_length))
    B_list = []
    n = np.arange(1, N + 1, 1)
    theta = 0
    for i in range(0, len(theta_m_arr)):
        theta = theta_m_arr[i]
        B_list.append(np.sin(n[:] * theta) * ((2 / (np.pi * cm_over_b[i])) + (n[:] / np.sin(theta))))
    return B_list
    
def CL_calc(A_1, AR):                                               # C_L  = pi * AR * A_1, where A_1 = a coefficient * angle
    return np.pi * AR * A_1

def CD_inducedCalc(A_coeff, AR):                                                # C_D,i = pi * AR * SUM(n * (A_n ** 2))
    CDi = 0
    for i in range(0, len(A_coeff), 1):
        if abs(A_coeff[i]) >= 1e-12:
           CDi += (i+1) * (A_coeff[i] ** 2)
    return CDi * AR * np.pi
            
def eCalc(C_L, AR, C_Di):                                                       # e = (C_L ** 2) / (AR * pi * C_D,i)
    return (C_L ** 2) / (np.pi * AR * C_Di)
    
def liftDistributionCalc(A_coeff, theta_m):
    lift_distr = [0.0]
    for theta in theta_m:
        if abs(A_coeff[0]) >= 1e-12:
            lift_distr.append((A_coeff[0]) * np.sin(1 * theta))
        else:                                                                   # For non-symmetric Coefficients about midspan
            lift_distr.append(0.0)
        for i in range(1, len(A_coeff), 1):
            if abs(A_coeff[i]) >= 1e-12:
                lift_distr[-1] += (A_coeff[i]) * np.sin((i+1) * theta)
        lift_distr[-1] = 2 * lift_distr[-1]
    lift_distr.append(0.0)
    return lift_distr       

def ellipticDistr(A_1, theta_m):
    elliptic_lift = [0.0]
    for theta in theta_m:
        elliptic_lift.append(2*(A_1) * np.sin(theta))
    elliptic_lift.append(0.0)
    return elliptic_lift

def lift(S, C_L, total_alpha, q):
    lift = q * S * (C_L * total_alpha)
    return lift

def induced_drag(S, C_Di, total_alpha, q):
    induced_drag = q * S * (C_Di * (total_alpha**2))
    return induced_drag
           
############   MAIN   ###########

theta_m = thetaM(N)


B_matrix = np.array(B(N, c_zero_over_b, wing_tip_length, theta_m))

B_inv = np.linalg.inv(B_matrix)

A_coeff = np.dot(B_inv, alpha_coeff_matrix)


C_L = CL_calc(A_coeff[0], AR)


C_Di = CD_inducedCalc(A_coeff, AR)

e = eCalc(C_L, AR, C_Di)

L = lift(S, C_L, total_alpha, q)

D_induced = induced_drag(S, C_Di, total_alpha, q)

#############     Printing       #############

print('A coefficients:\n')
for i in range(0, len(A_coeff), 1):
    print(A_coeff[i], '\u03B1')

print('\n\nC_L value:\n')
print(C_L, '\u03B1, where \u03B1 includes \u03B1 and \u03B1 where lift = 0')

print('\n\nC_D,i value:\n')
print(C_Di, '\u03B1\u00B2, where \u03B1 includes \u03B1 and \u03B1 where lift = 0')

print('\n\nSpan Efficiency Factor e:\n')
print(e)

print('\n\nTotal Wing Lift L:\n')
print(L, 'N')

print('\n\nTotal Wing Induced Drag D_i:\n')
print(D_induced, 'N')

print('\n\nTotal aircraft weight W:\n')
print(aircraft_weight, 'N')


#############      Plotting        ###############

ym_over_b = ymOverB(theta_m)

lift_distr = liftDistributionCalc(A_coeff, theta_m)

elliptic_distr = ellipticDistr(A_coeff[0], theta_m)

ym_over_b_list = ym_over_b.tolist()
ym_over_b_list.insert(0, 0.5)
ym_over_b_list.append(-0.5)

plt.plot(ym_over_b_list, lift_distr, label = 'Total \u0393 Distribution')
plt.plot(ym_over_b_list, elliptic_distr, label = 'Elliptic Distribution')
plt.title('\u0393 Distribution Along Wing')
plt.xlabel('y/b')
plt.ylabel('\u0393/bU')
plt.legend(loc = 'best')
plt.show()