import matplotlib
import matplotlib.pyplot as plot
from matplotlib import animation
import os
from datetime import datetime
from create_sim import Simulation
from CR3BP import V
import numpy as np
# matplotlib.use('Agg')

def plot_traj(simulation):

    x = np.arange(-2, 2, 0.01)
    X, Y = np.meshgrid(x, x)
    Z = V(X, Y, simulation.mu)

    contour_levels = V(simulation.initial_point[0], simulation.initial_point[1], simulation.mu)
    contours = plot.contour(X, Y, Z, levels=contour_levels, linestyles='dashed')
    reuslts = simulation.trajectory
    plot.plot(reuslts[:, 0], reuslts[:, 1], label='Position')
    plot.plot(reuslts[0, 0], reuslts[0, 1], 'go')
    plot.plot(reuslts[-1, 0], reuslts[-1, 1], 'bo')


    plot.plot(-1 * simulation.mu, 0, 'yo', markersize=12)
    plot.plot(1 - simulation.mu, 0, 'ro', markersize=12)
    # plot.plot(point[0], point[1], 'rx')
    plot.xlabel("x (mu)")
    plot.ylabel("y (mu)")
    # plot.title(title)

    for key, point in simulation.eq_points.items():
        plot.plot(point[0], point[1], 'rx')

    # legend_elememnts = [Line2D([0], [0], marker='o', color='w', label='Start', markerfacecolor='g', markersize=5),
    #                     Line2D([0], [0], marker='o', color='w', label='End', markerfacecolor='b', markersize=5)]
    # plot.legend(handles=legend_elememnts)

    plot.show()

def create_animation(simulation, save=False, hide_axes=True, filename=None):
    fig = plot.figure(figsize=(15,10))
    ax = plot.axes()



    # ax.plot(0,0,'o', ms=8, c='C0')

    if simulation.type == 'Moon':
        color_1 = 'b'
        color_2 = 'gray'
    else:
        color_1 = 'yellow'
        color_2 = 'b'
    size1 = 20
    for key, point in simulation.eq_points.items():
        plot.plot(point[0], point[1], 'rx')

    # line, = ax.plot([],[],'-',lw=3,c='C0', animated=True)
    m3, = ax.plot([],[],'o',ms = 3, c='g')
    path3, = ax.plot([],[],'-', c='k')
    m1 = ax.plot(-1 * simulation.mu, 0, 'o', c=color_1, ms=size1)
    m2 = ax.plot(1 - simulation.mu, 0, 'o', c=color_2, ms=size1/2)

    lims = 1.5
    ax.set_xlim(-1 * lims,lims)
    ax.set_ylim(-1 * lims,lims)

    if hide_axes:
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
    ax.set_aspect('equal')

    # xlims, ylims = ax.get_xlim(), ax.get_ylim()
    x, y = np.arange(-2,2, 0.01), np.arange(-2,2, 0.01)
    X, Y = np.meshgrid(x, y)
    Z = V(X, Y, simulation.mu)

    print(simulation.contour_levels)
    contours = plot.contour(X, Y, Z, levels=simulation.contour_levels, linestyles='dashed')


    fig.tight_layout()

    skip = int(0.02/(simulation.time[-1]/len(simulation.time)))  # Should be the time step?
    # skip = int()
    def animate(i):
        # i *= skip
        # line.set_data([0,simulation.trajectory[i,0]])
        # m3.set_data(simulation.trajectory[0:i:skip,0], simulation.trajectory[0:i:skip,1])
        path3.set_data(simulation.trajectory[0:i,0], simulation.trajectory[0:i,1])
        m3.set_data(simulation.trajectory[i,0], simulation.trajectory[i,1])
        return m3, path3

    anim = animation.FuncAnimation(fig=fig, func=animate, frames=len(simulation.time)//skip,
                                   interval=0.2, blit=True)

    plot.title("Orbit Starting at {0} in the Earth-{1} System".format(simulation.initial_point_str, simulation.type))

    if save:
        filename = datetime.now().strftime('%Y-%m-%d_%H-%M') if filename is None else filename

        if not os.path.isdir(os.path.join(os.getcwd(), 'animations')):
            os.mkdir('animations')

        try:
            writer = animation.AVConvWriter(fps=50, bitrate=-1)
            anim.save('./animations/{}.mp4'.format(filename), writer=writer)
        except FileNotFoundError:
            writer = animation.FFMpegWriter(fps=50, bitrate=-1)
            anim.save('./animations/{}.mp4'.format(filename), writer=writer)
    else:
        plot.show()

if __name__ == '__main__':
    sim = Simulation()
    sim.make_simulation(36*np.pi)
    print(len(sim.time), (len(sim.time)/30)/60)
    # plot_traj(sim)
    create_animation(sim, save=False)