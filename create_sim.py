import CR3BP
import numpy as np
from scipy import optimize
from functools import partial

m_Earth = 5.972E24  # kg
m_Moon = 7.34E22  # kg
m_Sun = 1.98E30  # kg

mu_Earth = m_Earth / (m_Earth + m_Sun)
mu_Moon = m_Moon / (m_Earth + m_Moon)


class Simulation:

    def __init__(self):
        """
        Containter for simulation solution and necessary plotting variables and items
        :param mu:
        :param type: Earth-Sun or Earth-Moon system
        :param direction: Forwards or backwards time integration. This defines unstable/stable manifolds for trajectory
        simulations
        """
        self.mu = None
        self.type = None
        self.direction = None
        self.eq_points = None
        self.init_conds = None
        self.time = None
        self.trajectory = None
        self.initial_point = None
        self.contour_levels = []


    def find_lagrange_points(self, mu=None):

        if mu is None:
            mu = self.mu

        bound_func = partial(CR3BP.Vx, y=0,
                             mu=mu)  # Bind the potential to be on the line y=0 with the value of mu so optimize.fsolve only has to do a 1D solve for the x-point
        L1 = optimize.fsolve(bound_func, np.array([0.5]))
        L2 = optimize.fsolve(bound_func, np.array([1]))
        L3 = optimize.fsolve(bound_func, np.array([-1]))

        return {"L1": np.array([L1[0], 0]), "L2": np.array([L2[0], 0]), "L3": np.array([L3[0], 0]),
                "L4": np.array([0.5 - mu, np.sqrt(3) / 2]), "L5": np.array([0.5 - mu, -np.sqrt(3) / 2])}



    def pick_random_initial_conditions(self, determine_dir = False, determine_mass = False, pos_epsilon=5E-3, vel_espilon=1E-7):
        """
        There are a few random conditions to pick
        Type: Earth-Sun or Earth-Moon system
        Lagrange Point: Which point to start near (L1, L2, L3, L4 or L5)
        Direction: integrate forwards or backwards in time
        Initial conditions: x, y, xdot and ydot all need to be defined appropriately
        """

        # Pick Sun or Moon System
        if determine_mass:
            if np.random.random() > 0.5:
                self.type = 'Sun'
                self.mu = mu_Earth
            else:
                self.type = 'Moon'
                self.mu = mu_Moon
        else:
            self.type = 'Moon'
            self.mu = mu_Moon

        self.eq_points = self.find_lagrange_points(self.mu)

        for k, point, in self.eq_points.items():
            self.contour_levels.append(CR3BP.V(point[0], point[1], self.mu))

        self.contour_levels.sort()

        # pick lagrange point
        lagrange_point = np.random.randint(1,6)
        init_xy = self.eq_points["L{}".format(lagrange_point)]
        self.initial_point = self.eq_points["L{}".format(lagrange_point)]
        self.initial_point_str = "L{}".format(lagrange_point)

        # pick random offsets from lagrange point and a random velocity in x-y direction
        self.init_conds = np.concatenate([np.random.uniform(init_xy[0] - pos_epsilon, init_xy[0] + pos_epsilon, [1,1]),
                                     np.random.uniform(init_xy[1] - pos_epsilon, init_xy[1] + pos_epsilon,[1,1]),
                                     np.random.uniform(-1 * vel_espilon, vel_espilon, [1, 2])], axis=1).reshape(-1)

        if determine_dir:
        # Pick to either simulate forwards or backwards in time
        # Forwards in time will show unstable manifolds while backwards in time will show stable manifolds
            if np.random.random() > 0.5:
                # integrate fowards in time
                self.direction = True
            else:
                # integrate backwards in time
                self.direction = False
        else:
            self.direction = True


    def simulate_trajectory(self, time_length = 6 * np.pi):

        if not self.direction:
            time_length = -1 * time_length

        bound_differential_equation = partial(CR3BP.rotating_CR3BP_DEs, self.mu)
        self.time, self.trajectory = CR3BP.rk45(bound_differential_equation, self.init_conds, 0, time_length, 0.01)

        if not self.direction:
            # If you integrate backwards in time, the last data point becomes the first data point
            self.time = self.time + np.min(self.time)


    def make_simulation(self, time_length=6*np.pi):
        self.pick_random_initial_conditions()
        self.simulate_trajectory(time_length)


if __name__ == '__main__':
    Sims = Simulation()
    Sims.make_simulation()