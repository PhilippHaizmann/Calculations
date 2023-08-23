import math
import numpy as np
import matplotlib.pyplot as plt

# Constants
ELECTRON_MASS = 5.68563e-12  # in eVm^-1s^-1
REDUCED_PLANC = 6.5821995e-16
BOLTZMANN_CONSTANT = 1.38064825e-23
ELEMENTARY_CHARGE = 1.602176634e-19
VACUUM_PERMITTIVITY = 8.854187817e-12
AVOGADRO_NUMBER = 6.02214086e23

# Calculations related to electrons and their energy
def calc_k_paralell(E_ex, ac_ang):
    k_para = (math.sin(math.radians(ac_ang)) * math.sqrt(2 * ELECTRON_MASS * E_ex) / REDUCED_PLANC**2) / 1e10
    return k_para

def calc_RT_energy(temp_room):
    temp_kelvin = 273.15 + temp_room
    E = (3/2) * (BOLTZMANN_CONSTANT / ELEMENTARY_CHARGE * temp_kelvin) * 1000
    return E

def plot_e_v_t(t_start, t_end, points):
    x = np.linspace(t_start, t_end, points)
    y = [calc_RT_energy(i) for i in x]

    plt.plot(x, y, linewidth=2.0)
    plt.xlabel("Temperature [C]")
    plt.ylabel("Energy [meV]")
    plt.show()

# Calculations for mean free path
def calc_SeahDench_organic_fmp(E_ex, dens):
    fmp = 49/E_ex**2 + 0.11 * math.sqrt(E_ex)
    return fmp, fmp/dens

def fmp_seah(M, rho, N, E):
    a = (M * 1e21 / (rho * N * AVOGADRO_NUMBER))**(1/3)
    l_m = (2170 / E**2) + (0.72 * math.sqrt(a * E))
    return a * l_m

def fmp_seah_S1(g, h, Z_g, Z_h, M, rho, Eg, E):
    a = (M * 1e21 / (rho * (g + h) * AVOGADRO_NUMBER))**(1/3)
    Z = (g * Z_g + h * Z_h) / (g + h)
    return (4 + 0.44 * Z**0.5 + 0.104 * E**0.872) * a**1.7 / (Z**0.3 * (1 - 0.02 * Eg))

def tpp_2m(N_v, rho, M, E_g, E):
    U = N_v * rho / M
    E_p = 28.8 * math.sqrt(U)
    gamma = 0.191 * rho**-0.50
    beta = -0.1 + 0.944 / math.sqrt(E_p**2 + E_g**2) + 0.069 * rho**0.1
    C = 1.97 - 0.91 * U
    D = 53.4 - 20.8 * U
    l = E / (E_p**2 * (beta * math.log(gamma * E) - C / E + D / E**2))
    return l * 0.1

# Calculations related to semiconductors
def debey_semiconductor(temperature, electron_density_power, epsilon_r):
    epsilon = epsilon_r * VACUUM_PERMITTIVITY
    electron_density = 10**electron_density_power
    lambda_D = math.sqrt((epsilon * BOLTZMANN_CONSTANT * temperature) / (ELEMENTARY_CHARGE**2 * electron_density))
    return lambda_D

# Uncomment the below line to test the function
print(fmp_seah_S1(1, 2, 42, 16, 160.07, 5.06, 1.23, 1252.1))
