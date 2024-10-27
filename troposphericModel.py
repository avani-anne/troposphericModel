import numpy as np
import math as m
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp
matplotlib.use('Agg')

# Follows pep-8 style guidelines


class troposphericModel():

    # Definition of rate constants
    k4 = 6.9E-8
    k5 = 1.9E-13
    k6 = 1.5E-15
    k7 = 7.9E-15
    k8 = 3.9E-3
    k9 = 9.6E-12
    k10 = 1.3E-11

    def __init__(self, FCO, FNO, FO3, t0, tf, X0):

        # time range
        self.t0 = t0
        self.tf = tf

        # emissions/static sources (forcings)
        self.FCO = FCO
        self.FNO = FNO
        self.FO3 = FO3
        self.t0 = t0
        self.tf = tf

        # Definition of inital conditions
        self.X0 = X0

    def dXdt(self, t, X):
        """ Returns differential equations for different compounds

        Args:
            t: float, time
            X: list, initial conditions
        """
        CO, NO, O3, NO2, OH, HO2 = X

        k4 = troposphericModel.k4
        k5 = troposphericModel.k5
        k6 = troposphericModel.k6
        k7 = troposphericModel.k7
        k8 = troposphericModel.k8
        k9 = troposphericModel.k9
        k10 = troposphericModel.k10

        dCOdt = + self.FCO - k5*CO*OH
        dNOdt = + self.FNO + k8*NO2 - k7*NO*O3 - k9*HO2*NO
        dO3dt = + self.FO3 + k8*NO2 - k4*O3-k6*HO2*O3 - k7*NO*O3
        dNO2dt = + k7*NO*O3+k9*HO2*NO - k8*NO2 - k10*OH*NO2
        dOHdt = + 2*k4*O3 + k6*HO2*O3 + k9*HO2*NO - k5*OH*CO - k10*OH*NO2
        dHO2dt = + k5*OH*CO - k6*HO2*O3 - k9*HO2*NO

        return dCOdt, dNOdt, dO3dt, dNO2dt, dOHdt, dHO2dt

    def solve(self, spin_up):
        """ Returns array of concentrations for each compound and array of time

        Args:
            spin_up: float, decimal fraction of transient data
        """
        # defining integration time
        X = solve_ivp(self.dXdt, (self.t0, self.tf), self.X0, method='LSODA')
        X.t = X.t/(24*60*60)  # convert time to days

        CO, NO, O3, NO2, OH, HO2 = X.y

        idx = round(len(CO)*spin_up)

        return CO[idx:], NO[idx:], O3[idx:], NO2[idx:], \
            OH[idx:], HO2[idx:], X.t[idx:]


def run():
    """
    This function creates a gif of a [CO] vs. [OH] phase plot
    """
    FCO = 5E5
    FNO = 6.3E4
    FO3 = 6E4
    t0, tf = 0, 10E7

    X0 = 5E11, 5E9, 3E11, 5E9, 5E6, 2.5E8

    spin_up = .5  # amount of spin up, arbitrarily chosen
    model = troposphericModel(FCO, FNO, FO3, t0, tf, X0)
    CO, NO, O3, NO2, OH, HO2, time = model.solve(spin_up)

    # Plotting and creating animation
    # Collaborated with Ana Beatriz Studart on this portion
    fig, ax = plt.subplots()
    line, = ax.plot([], [], 'r-')  # Initialize an empty line
    ax.grid()
    ax.set_xlabel("[CO] molec/$cm^3$", fontsize=30)
    ax.set_ylabel("[OH] molec/$cm^3$", fontsize=30)
    ax.set_xlim(min(CO), max(CO))
    ax.set_ylim(min(OH), max(OH))
    ax.set_title("[CO] molec/$cm^3$ vs. [OH] molec/$cm^3$", fontsize=35)

    def animate(i):
        line.set_data(CO[:i], OH[:i])
        print("Animating frame", i)
        return line,

    ani = FuncAnimation(fig, animate, frames=len(CO), interval=50)
    ani.save('ani.gif', writer='pillow')
