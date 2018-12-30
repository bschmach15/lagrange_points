import numpy as np
from scipy import integrate

def rk45(dxdt, x_0, t0, tf, max_step):
    """
    Explicitly use RK45 scheme to solve an ODE. This is a lower level function that gives more control than the one
    liner solver that is also included in scipy.
    :param dxdt: Function to integrate
    :param x_0: initial conditions
    :param t0: initial time
    :param tf: final time
    :param step: time steps
    :return: times, solution array [x, y, z, u, v, w]
    """
    integrator = integrate.RK45(dxdt, t0, x_0, tf, max_step)
    times = np.array([[t0]])
    soln = np.array([x_0])
    while integrator.status != 'finished':
        integrator.step()
        soln = np.append(soln, np.array([integrator.y]), axis=0)
        times = np.append(times, np.array([[integrator.t]]), axis=0)
    return times, soln


def V(x, y, mu):
    """
    Calcualte the potential energy of a 3rd mass in the CR3PB system
    :param x: x-coordinate of the mass
    :param y: y-coordinate of the mass
    :param mu: mass fraction of the system
    :return: the calculated potential energy of the 3rd mass at any point in the plane.
    """
    a = 0.5 * (x ** 2 + y ** 2)
    b = (1 - mu) / (np.sqrt((x + mu) ** 2 + y ** 2)) + mu / (np.sqrt((x - 1 + mu) ** 2 + y ** 2))
    return a + b


def Vx(x, y, mu):
    """
    x-derivative of the potential energy V(x,y,mu)
    :param x:
    :param y:
    :param mu:
    :return:
    """
    a = mu * (-1 + x + mu)
    b = (1 - mu) * (x + mu)
    a_denom = (y ** 2 + (x - 1 + mu) ** 2) ** (3 / 2)
    b_denom = (y ** 2 + (x + mu) ** 2) ** (3 / 2)
    return x - a / a_denom - b / b_denom


def Vy(x, y, mu):
    """
    y-derivative of the potential energy V(x,y,mu)
    :param x:
    :param y:
    :param mu:
    :return:
    """
    a = y * mu
    b = y * (1 - mu)
    a_demon = (y ** 2 + (x - 1 + mu) ** 2) ** (3 / 2)
    b_denom = (y ** 2 + (x + mu) ** 2) ** (3 / 2)
    return y - a / a_demon - b / b_denom


def jacobi(x, y, x_dot, y_dot, mu):
    """
    Calculate Jacobi's constant for a 3rd mass. This defines the forbidden regions and zero-velocity curves for an object
    :param x: x-position of the mass
    :param y: y-position of the mass
    :param x_dot: x-velocity of the mass
    :param y_dot: y-velocity of the mass
    :param mu: mass fraction of the system
    :return:
    """
    a = 0.5 * (x_dot ** 2 + y_dot ** 2)
    return a - V(x, y, mu)


def rotating_CR3BP_DEs(mu, t, x_0):
    """
    Differential equations for the CR3PB in the rotatin reference frame. For use in RK45 scheme to solve for motion
    through time.
    :param mu: Mass fraction of the system - must be bound to function with something like functools.partial for the
    RK45 scheme defined above to work
    :param t: time of evaluation
    :param x_0: initial parameters
    :return:
    """
    x = x_0[0]
    y = x_0[1]

    u = x_0[2]  # x dot
    v = x_0[3]  # y dot

    du = 2 * v + Vx(x, y, mu)
    dv = - 2 * u + Vy(x, y, mu)

    return [u, v, du, dv]