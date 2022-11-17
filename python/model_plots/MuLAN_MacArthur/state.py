"""MuLAN_MacArthur-model specific plot function"""

from math import log
from posixpath import dirname
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib import ticker
from matplotlib.collections import LineCollection
from scipy.signal import find_peaks
import itertools
import networkx as nx

from utopya import DataManager, UniverseGroup, MultiverseGroup
from utopya.plotting import (is_plot_func, PlotHelper, UniversePlotCreator, MultiversePlotCreator)

import os.path

import xarray as xr

from .helper import *


# -----------------------------------------------------------------------------


def save_to_file(dm: DataManager, *,
                 uni: UniverseGroup,
                 out_path: str,
                 save_kwargs: dict = None,
                 **plot_kwargs):

    keywords = get_argument_or_default(plot_kwargs, "save", ["_types", "_traits", "_masses"])

    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']
    config = uni["cfg"]

    time_list = []
    for t in grp:
        time_list.append(int(t))

    time_list = sorted(time_list)
    frame = grp[str(time_list[-1])]

    save_list = []
    for keyword in keywords:
        save_list.append(frame[keyword])

    out_path = out_path.split(".pdf")[0] + ".graph"

    with open(out_path, "w") as f:

        for v in zip(*save_list):
            f.write(','.join([str(np.array([x])[0]) for x in v]) + "\n")
    return


# graph in 3d plot
@is_plot_func(creator_type=UniversePlotCreator)
def plot_graph(dm: DataManager, *,
                uni: UniverseGroup,
                hlpr: PlotHelper,
                xlim: tuple = None,
                ylim: tuple = None,
                zlim: tuple = None,
                time: int = -1,
                save_kwargs: dict = None,
                **plot_kwargs):

    ''' Plot the Graph in a 3D trait_space embedding

    Args:
        dm: datamanager
        uni: universe group
        hlpr: plot helper
        xlim: x axis limit
        ylim: y axis limit
        zlim: z axis limit
        time: The time index to plot. timestep = time * write_every
        save_kwargs:
        plot_kwargs:

    '''

    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "x")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "y")
    zlabel = get_argument_or_default(plot_kwargs, "zlabel", "mass")

    azim = get_argument_or_default(plot_kwargs, "azim", -60)
    elev = get_argument_or_default(plot_kwargs, "elev", None)
    rot = get_argument_or_default(plot_kwargs, "rotation", 0)

    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "_niche_w")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_types")
    s_name = get_argument_or_default(plot_kwargs, "s_name", "_masses")
    w_name = get_argument_or_default(plot_kwargs, "w_name", "_graph_edge_v")

    s_scale = get_argument_or_default(plot_kwargs, "s_scale", 200)
    w_scale = get_argument_or_default(plot_kwargs, "w_scale", 3)

    timestamp = get_argument_or_default(plot_kwargs, "timestamp", True)

    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']
    config = uni["cfg"]


    fig = hlpr.fig
    out_path = hlpr.out_path

    time_list = get_time_grps(grp)
    frame = grp[str(time_list[time])]

    vertices = frame["_graph_vertex"]
    edges = frame["_graph_edges"]
    weights = frame[w_name]
    types = frame[z_name]
    traits = frame[x_name]
    masses = frame[s_name]
    nw = frame[y_name]


    edge_l = []
    weight_l = []
    for i in zip(*edges):
       edge_l.append(np.array(i))
    for i in weights:
       weight_l.append(float(i.data))


    pos = {}
    node_l = []
    for i, j, k, l in zip(vertices, traits, nw, types):
        v_i = int(i.data)
        node_l.append(v_i)
        pos[v_i] = np.array([j, k, l])

    mass_l = []
    for i in masses:
        mass_l.append(float(i.data))


    G = nx.DiGraph(edge_l)
    G.add_nodes_from(node_l)

    fig.clf()
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(azim=azim)
    fig.tight_layout()

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if zlim is not None:
        ax.set_ylim(zlim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel, rotation=270)

    #nx.draw(G, ax=hlpr.ax)

    if (timestamp):
        dt2 = float(config["MuLAN_MA"]["dt2"])
        # t_max = int(config["num_steps"])
        # write_every = int(config["write_every"])
        time_s = time_list[time]

        fig.text(0.1, 0.9, "t=" + str(time_s * dt2) + " (" + str(time_s) + "dt)")

    edge_xyz = np.array([(pos[u], pos[v]) for u, v in edge_l])
    ax.scatter(traits, nw, types, s=(np.array(mass_l)/max(mass_l)*s_scale))

    m_w = max(weight_l)

    for (vizedge, w) in zip(edge_xyz, weight_l):
        ax.plot(*vizedge.T, color="tab:gray", linewidth=w/m_w*w_scale)

    fig.tight_layout()
    return


@is_plot_func(creator_type=UniversePlotCreator)
def plot_graph_properties(dm: DataManager, *,
               uni: UniverseGroup,
               hlpr: PlotHelper,
               xlim: tuple = None,
               ylim: tuple = None,
               zlim: tuple = None,
               time: int = -1,
               save_kwargs: dict = None,
               **plot_kwargs):

    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "x")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "y")

    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "_niche_w")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_types")
    m_name = get_argument_or_default(plot_kwargs, "m_name", "_masses")
    w_name = get_argument_or_default(plot_kwargs, "w_name", "_graph_edge_v")

    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']
    config = uni["cfg"]

    fig = hlpr.fig
    out_path = hlpr.out_path

    time_list = get_time_grps(grp)

    density_l = []
    for i in time_list:
        frame = grp[str(i)]

        G = get_graph_from_frame(frame, w_name, m_name)
        density_l.append(nx.density(G))

    hlpr.ax.plot(time_list, density_l)


# plot trait space for time (2d)
def plt_biomass_by_type(dm: DataManager, *,
                        uni: UniverseGroup,
                        out_path: str,
                        time: int = 0,
                        fmt: str = None,
                        xlim: tuple = None,
                        ylim: tuple = None,
                        save_kwargs: dict = None,
                        **plot_kwargs):
    """Plot biomasses of every species class

    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The selected universe data
        out_path (str): Where to store the plot to
        fmt (str, optional): the plt.plot format argument
        save_kwargs (dict, optional): kwargs to the plt.savefig function
        **plot_kwargs: Passed on to plt.plot
    """
    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']

    frame = grp[str(time)]

    types = frame["_types"]
    traits = frame["_traits"]
    masses = frame["_masses"]
    # cr = frame["_cr"]

    if not xlim:
        xlim = (-20, 20)
    elif not xlim == "auto":
        plt.xlim(xlim)

    if not ylim:
        ylim = (0.0, 1.0)
    elif not ylim == "auto":
        plt.ylim(ylim)

    for i in range(2):
        xy = sorted(zip(traits[types == i], masses[types == i]))
        # print(*zip(*xy))

        plt.scatter(*zip(*xy), **plot_kwargs)

    # Save and close figure
    save_and_close(out_path, save_kwargs=save_kwargs)


@is_plot_func(creator_type=UniversePlotCreator)
def plt_timestep_dev(dm: DataManager, *,
                     uni: UniverseGroup,
                     hlpr: PlotHelper,
                     fmt: str = None,
                     xlim: tuple = None,
                     ylim: tuple = None,
                     save_kwargs: dict = None,
                     cmap: str = "RdBu",
                     **plot_kwargs):
    """Plot biomasses of every species class

    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The selected universe data
        out_path (str): Where to store the plot to
        fmt (str, optional): the plt.plot format argument
        save_kwargs (dict, optional): kwargs to the plt.savefig function
        **plot_kwargs: Passed on to plt.plot
    """
    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']
    config = uni["cfg"]

    integrator = config["MuLAN_MA"]["integrator"]
    try:
        param = config["MuLAN_MA"]["error"]
    except:
        param = config["MuLAN_MA"]["dt"]

    time_list = []

    for t in grp:
        time_list.append(int(t))

    m = []
    dt_a = []
    ta = [0]
    cs = []
    c = 0
    for t in sorted(time_list):
        frame = grp[str(t)]
        dts = frame["dts"]
        # types = frame["_types"]
        traits = frame["_traits"]
        masses = frame["_masses"]
        m.append(masses[traits == 1])

        for dt in dts:
            dt_a.append(dt)
            ta.append(ta[-1] + dt)
            cs.append(c)

        # alternate c between 0 and 1 to mark dt2 timesteps
        c ^= 1

    # print(time_list)
    # plt.scatter(time_list, m)
    # plt.plot(dt_a)
    hlpr.ax.scatter(ta[:-1], dt_a, c=cs, marker=".", cmap=cmap)
    hlpr.ax.set_xlabel("t")
    hlpr.ax.set_ylabel("dt")
    hlpr.provide_defaults('set_title', title="Integrator: " + str(integrator) + " " + str(param))

    # Save and close figure
    # save_and_close(out_path, save_kwargs=save_kwargs)


@is_plot_func(creator_type=MultiversePlotCreator, use_dag=False)
def plt_timestep_multi(dm: DataManager, *,
                       hlpr: PlotHelper,
                       plt_dts: bool = False,
                       legend: bool = False,
                       save_kwargs: dict = None,
                       **plot_kwargs):
    mv = dm["multiverse"]

    xlim_r = 0
    for i in mv:
        uni = mv[i]
        config = uni["cfg"]
        data = uni["data"]["MuLAN_MA"]

        integrator = config["MuLAN_MA"]["integrator"]
        if integrator == "rkck":
            param = "Error: " + str(config["MuLAN_MA"]["error"])
        else:
            param = "dt: " + str(config["MuLAN_MA"]["dt"])

        time_list = get_time_grps(data)

        dt_a = []
        ta = [0]
        m2 = {}

        for t in time_list:
            frame = data[str(t)]
            dts = frame["dts"]
            types = frame["_types"]
            traits = frame["_traits"]
            error = frame["_error"]
            masses = frame["_masses"]
            # m.append(error[types == 1])

            for dt in dts:
                dt_a.append(dt)
                ta.append(ta[-1] + dt)

            for trait2 in traits[types == 1]:
                trait = trait2.values.item()

                if trait not in m2.keys():
                    m2[trait] = [[], []]
                m2[trait][0].append(error[types == 1][traits[types == 1] == trait])
                m2[trait][1].append(ta[-1])
                # ta2.append(ta[-1])

        # print(time_list)
        # plt.scatter(time_list, m)
        # plt.plot(dt_a)
        # print(plot_kwargs)
        if plt_dts:
            hlpr.ax.scatter(ta[:-1], dt_a, marker=".", label=integrator + "  " + param)
            hlpr.ax.set_xlabel("t")
            hlpr.ax.set_ylabel("dt")
        else:
            for k in sorted(m2.keys()):
                # hlpr.ax.scatter(ta2, m, marker=".", label=integrator + "  " + param)
                hlpr.ax.scatter(m2[k][1], m2[k][0], marker=".", s=2, label=integrator + "  " + param)
            hlpr.ax.set_xlabel("t")
            hlpr.ax.set_ylabel("error")

        if xlim_r == 0:
            xlim_r = ta[-1]
        else:
            xlim_r = min(xlim_r, ta[-1])
        # hlpr.provide_defaults('set_title', title="Integrator: " + str(integrator) + " " + str(param))
    hlpr.ax.set_xlim((0, xlim_r))
    if legend:
        hlpr.ax.legend()

# movie 2d both (pp and c) width by color
def mov_biomass_by_type(dm: DataManager, *,
                        uni: UniverseGroup,
                        out_path: str,
                        xlim: tuple = None,
                        ylim: tuple = None,
                        vlim: tuple = None,
                        plt_colorbar: bool = False,
                        save_kwargs: dict = None,
                        cmap: str = "autumn",
                        **plot_kwargs):
    """ Make a Movie of the biomass development

    Args:
        dm (DataManager): The data manager from which to retrieve the data
        uni (UniverseGroup): The selected universe data
        out_path (str): Where to store the plot to
        fmt (str, optional): the plt.plot format argument
        save_kwargs (dict, optional): kwargs to the plt.savefig function
        **plot_kwargs: Passed on to plt.plot
    """

    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "_masses")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_niche_w")

    x_label = get_argument_or_default(plot_kwargs, "x_label", "trait")
    y_label = get_argument_or_default(plot_kwargs, "y_label", "Biomass")
    z_label = get_argument_or_default(plot_kwargs, "z_label", "_niche_w")
    timestamp = get_argument_or_default(plot_kwargs, "timestamp", True)
    time_frame_l = get_argument_or_default(plot_kwargs, "timeframes", [])

    mov = get_argument_or_default(plot_kwargs, "mov", True)

    m_size = get_argument_or_default(plot_kwargs, "markersize", 15)
    m_size2 = get_argument_or_default(plot_kwargs, "markersize_pp", 5)

    font_size = get_argument_or_default(plot_kwargs, "font_size", 15)

    plt.rcParams.update({'font.size': font_size})

    # Get the group that all datasets are in
    grp = uni['data/MuLAN_MA']
    config = dm["cfg"]

    time_list = get_time_grps(grp)

    metadata = dict(title='MuLAN Mac Arthur Model', artist='Utopia',
                    comment='')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    fig = plt.figure(figsize=(16, 8))


    plt_l = [
        plt.plot([], [], 'b-o', markersize=m_size2)[0],
        # plt.plot([], [], 'ro', **plot_kwargs)[0]]

        # plt.scatter([], [], c=[], **plot_kwargs)[0],
        plt.scatter([], [], s=m_size, c=[], cmap=cmap)]

    if not xlim:
        xlim = (-20, 20)

    if xlim == "auto":
        xlim, ylim2, _ = get_maximums(grp, time_list)

    plt.xlim(xlim)

    if not ylim:
        ylim = (0.0, 1.0)

    if not vlim:
        vlim = (0.0, 20.0)

    plt_l[1].set_clim(vlim)

    if ylim == "auto":
        if xlim == "auto":
            ylim = ylim2
        else:
            xlim2, ylim, _ = get_maximums(grp, time_list)

    plt.ylim(ylim)

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if plt_colorbar:
        cbar = plt.colorbar(plt_l[1])
        cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
        cbar.ax.set_ylabel(z_label, rotation=270)

    if (timestamp):
        timestamp_t = fig.text(0.12, 0.9, "")

    def get_frame(time):

        frame = grp[str(time)]

        types = frame["_types"]
        traits = frame[x_name]
        masses = frame[y_name]
        nw = frame[z_name]

        xy = sorted(zip(traits[types == 0], masses[types == 0]))
        plt_l[0].set_data(*zip(*xy))

        plt_l[1].set_offsets(np.vstack((traits[types == 1], masses[types == 1])).T)
        plt_l[1].set_array(nw[types == 1])

        if (timestamp):
            dt2 = float(config["model"]["dt2"])

            timestamp_t.set_text("t=" + f'{time * dt2:.1f}' + " (" + str(time) + "dt)")


    if mov:
        with writer.saving(fig, out_path.split(".pdf")[0] + ".mp4", 100):
            for time in time_list:

                get_frame(time)

                writer.grab_frame()

    for time in time_frame_l:
        get_frame(time)
        plt.savefig(out_path.split(".pdf")[0] + "_" + str(time) + ".pdf", bbox_inches="tight")



# trait space over time, bifurcation diagram style
def plt_trait(dm: DataManager, *,
              uni: UniverseGroup,
              out_path: str,
              fmt: str = None,
              save_kwargs: dict = None,
              out_type: str = ".pdf",
              cmap: str = "YlOrBr",
              trait_name: str = "_traits",
              xlim: tuple = None,
              **plot_kwargs):

    fixed_size = get_argument_or_default(plot_kwargs, "fixed_size", None)

    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "Niche")
    zlabel = get_argument_or_default(plot_kwargs, "zlabel", "Biomass")

    font_size = get_argument_or_default(plot_kwargs, "font_size", 18)

    grp_main = uni['data/MuLAN_MA']

    dset_list = []
    for dset in grp_main:
        dset_list.append(int(dset))
    dset_list = sorted(dset_list)

    num_species_classes = grp_main.attrs["num_species_classes"][0]

    plots = {}
    mass = {}
    niche_width = {}
    min_m = {}
    max_m = {}

    for spec_class in range(num_species_classes):
        plots[spec_class] = []
        mass[spec_class] = []
        niche_width[spec_class] = []

    for dset in dset_list:
        grp = uni['data/MuLAN_MA/' + str(dset)]

        num_species = grp.attrs["num_species"]

        types = grp["_types"].values.reshape(num_species)
        traits = grp[trait_name].values.reshape(num_species)
        masses = grp["_masses"].values.reshape(num_species)
        # nw = grp["_niche_w"].reshape(num_species)

        for spec_class in range(num_species_classes):
            plots[spec_class].append(traits[types == spec_class])
            mass_v = masses[types == spec_class]
            mass[spec_class].append(mass_v)
            if spec_class in max_m.keys():
                max_m[spec_class] = np.maximum(max_m[spec_class], np.max(mass_v))
            else:
                max_m[spec_class] = np.max(mass_v)
            # niche_width[spec_class].append(nw[types == spec_class])

    out_path = out_path.split(".")

    for spec_class in range(num_species_classes):
        plt.clf()
        for xe, ye, ce in zip(dset_list, plots[spec_class], mass[spec_class]):
            if fixed_size:
                c_mappable = plt.scatter([xe] * len(ye), ye, s=fixed_size, c=ce, cmap=cmap)
            else:
                c_mappable = plt.scatter([xe] * len(ye), ye, ce, c=ce, cmap=cmap)

            plt.clim(0, max_m[spec_class])

        #plt.colorbar()
        plt.xlabel("time [dt]")
        plt.ylabel(ylabel)

        if xlim is not None:
            plt.xlim(xlim)

        cbar = plt.colorbar(c_mappable)
        cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
        cbar.ax.set_ylabel(zlabel, rotation=270)

        save_and_close(out_path[0] + "_" + str(spec_class) + out_type, save_kwargs=save_kwargs)


@is_plot_func(creator_type=UniversePlotCreator)
def plt_masses(dm: DataManager, *,
                   uni: UniverseGroup,
                   hlpr: PlotHelper,
                   fmt: str = None,
                   save_kwargs: dict = None,
                   out_type: str = ".png",
                   xlim: tuple = None,
                   ylim: tuple = None,
                   trait_name: str = "_masses",
                   spec_class = 1,
                   **plot_kwargs):

    grp_main = uni['data/MuLAN_MA']

    dset_list = get_time_grps(grp_main)
    num_species_classes = grp_main.attrs["num_species_classes"][0]

    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "time [dt]")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "Biomass")
    label = get_argument_or_default(plot_kwargs, "label", "_niche_w")
    llabel = get_argument_or_default(plot_kwargs, "llabel", "Width")

    mass = {}

    for dset in dset_list:
        grp = uni['data/MuLAN_MA/' + str(dset)]

        num_species = grp.attrs["num_species"]

        types = grp["_types"].values.reshape(num_species)
        masses = grp[trait_name].values.reshape(num_species)
        masses = masses[types == spec_class]
        labels = grp[label].values.reshape(num_species)
        labels = labels[types == spec_class]

        for i, j in zip(labels, masses):
            mass.setdefault(i, []).append(j)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    for l in mass:
        print(llabel + ": " + str(l))
        plt.plot(dset_list, mass[l], label=llabel + ": " + str(l))

    plt.legend()

# plot the total mass
def plt_total_mass(dm: DataManager, *,
              uni: UniverseGroup,
              out_path: str,
              fmt: str = None,
              save_kwargs: dict = None,
              out_type: str = "png",
              trait_name: str = "_masses",
              **plot_kwargs):

    grp_main = uni['data/MuLAN_MA']

    dset_list = []
    for dset in grp_main:
        dset_list.append(int(dset))
    dset_list = sorted(dset_list)

    num_species_classes = grp_main.attrs["num_species_classes"][0]

    plots = {}
    mass = {}
    niche_width = {}

    for spec_class in range(num_species_classes):
        plots[spec_class] = []
        mass[spec_class] = []
        niche_width[spec_class] = []

    for dset in dset_list:
        grp = uni['data/MuLAN_MA/' + str(dset)]

        num_species = grp.attrs["num_species"]

        types = grp["_types"].values.reshape(num_species)
        masses = grp[trait_name].values.reshape(num_species)

        for spec_class in range(num_species_classes):

            mass[spec_class].append(masses[types == spec_class])
            # niche_width[spec_class].append(nw[types == spec_class])

    out_path = out_path.split(".")



    for spec_class in range(num_species_classes):


        x = []
        y = []

        for xe, ce in zip(dset_list,  mass[spec_class]):
            x.append(xe)
            y.append(np.sum(ce))


        plt.clf()
        plt.xlabel("time")
        plt.ylabel("mass")

        plt.plot(x, y)

        save_and_close(out_path[0] + "_" + str(spec_class) + out_path[1], save_kwargs=save_kwargs)



def plt_biomass_dev(dm: DataManager, *,
              uni: UniverseGroup,
              out_path: str,
              fmt: str = None,
              save_kwargs: dict = None,
              out_type: str = ".pdf",
              plt_3d: bool = False,
              labels: dict = {},
              **plot_kwargs):


    grp_main = uni['data/MuLAN_MA']

    dset_list = []
    for dset in grp_main:
        dset_list.append(int(dset))
    dset_list = sorted(dset_list)

    num_species_classes = grp_main.attrs["num_species_classes"][0]

    mass = {}

    for spec_class in range(num_species_classes):
        mass[spec_class] = []

    for dset in dset_list:
        grp = uni['data/MuLAN_MA/' + str(dset)]

        num_species = grp.attrs["num_species"]

        types = grp["_types"].values.reshape(num_species)
        masses = grp["_masses"].values.reshape(num_species)
        # nw = grp["_niche_w"].reshape(num_species)

        for spec_class in range(num_species_classes):
            mass[spec_class].append(np.nansum(masses[types == spec_class]))

    out_path = out_path.split(".")

    for i in mass:
        plt.plot(dset_list, mass[i], label=labels[i])

    plt.xlabel("time")
    plt.ylabel("mass")


    save_and_close(out_path[0] + "_time_series." + out_type, save_kwargs=save_kwargs)

    if not get_argument_or_default(plot_kwargs, "phase_space", True):
        return

    if plt_3d and len(mass) >= 3:
        for comb in itertools.cobinations(mass, 2):

            plt.clf()
            ax = fig.add_subplot(111, projection="3d")
            for i in mass:
                ax.plot(mass[comb[0]], mass[comb[1]])

            if comb[0] in labels.keys():
                ax.set_xlabel(labels[comb[0]])
            else:
                ax.set_xlabel(comb[0])

            if comb[1] in labels.keys():
                ax.set_ylabel(labels[comb[1]])
            else:
                ax.set_ylabel(comb[1])

            if comb[2] in labels.keys():
                ax.set_zlabel(labels[comb[2]])
            else:
                ax.set_zlabel(comb[2])

            save_and_close(out_path[0] + "_phase_space" +
                           str(comb[0]) + "_" + str(comb[1]) +
                           "_" + str(comb[2]) +
                           out_type, save_kwargs=save_kwargs)
    else:
        for comb in itertools.combinations(mass, 2):

            plt.clf()
            for i in mass:
                plt.plot(mass[comb[0]], mass[comb[1]])

            if comb[0] in labels.keys():
                plt.xlabel(labels[comb[0]])
            else:
                plt.xlabel(comb[0])

            if comb[1] in labels.keys():
                plt.ylabel(labels[comb[1]])
            else:
                plt.ylabel(comb[1])

            save_and_close(out_path[0] + "_phase_space" +
                           str(comb[0]) + "_" + str(comb[1]) +
                           out_type, save_kwargs=save_kwargs)


@is_plot_func(creator_type=MultiversePlotCreator, use_dag=False)
def total_mass_fixed_y(dm: DataManager, *,
                              hlpr: PlotHelper,
                              save_kwargs: dict = None,
                              cmap: str = "Spectral",
                              xlim: tuple = None,
                              ylim: tuple = None,
                              zlim: tuple = None,
                              **plot_kwargs):
    fig = hlpr.fig
    out_path = hlpr.out_path
    mv = dm["multiverse"]
    config = dm["cfg"]


    # Get all the configuration options
    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "time [dt]")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "Niche width")
    zlabel = get_argument_or_default(plot_kwargs, "zlabel", "Total Biomass")

    spec_type = get_argument_or_default(plot_kwargs, "spec_type", 1)

    colorbar = get_argument_or_default(plot_kwargs, "colorbar", False)
    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "fixed_y")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_masses")

    azim = get_argument_or_default(plot_kwargs, "azim", -60)
    elev = get_argument_or_default(plot_kwargs, "elev", None)
    rot = get_argument_or_default(plot_kwargs, "rotation", 0)

    cmap = get_argument_or_default(plot_kwargs, "cmap", "Spectral")

    font_size = get_argument_or_default(plot_kwargs, "font_size", 15)

    plt.rcParams.update({'font.size': font_size})

    t = []
    m = []
    w = []

    for i in mv:

        grp_main = mv[i]["data"]["MuLAN_MA"]
        time_list = get_time_grps(grp_main)
        w_temp = mv[i]["cfg"]["MuLAN_MA"][y_name]
        for j in time_list:
            grp = grp_main[str(j)]

            num_species = grp.attrs["num_species"]

            types = grp["_types"].values.reshape(num_species)
            masses = grp[z_name].values.reshape(num_species)

            total_mass = np.sum(np.array(masses[types == spec_type]))

            t.append(j)
            m.append(total_mass)
            w.append(w_temp)

    #c_mappable = plt.scatter(t, m, c=w, marker=".", alpha=0.5)
    t = np.array(t)
    w = np.array(w)
    m = np.array(m)

    #max_t = np.max(t)
    #plt.plot(w[t == max_t], m[t == max_t])

    #plt.xlabel(ylabel)
    #plt.ylabel(zlabel)
    #return
    fig.clf()
    ax = fig.add_subplot(111, projection="3d")
    ax.view_init(azim=azim)

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel, rotation=270)

    #c_mappable = ax.scatter(t, m, w, marker=".")
    c_mappable = plt_3d_frame(ax, t, w, m, cmap, False, {})

    #if colorbar:
    #    cbar = fig.colorbar(c_mappable)
    #    cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
    #    cbar.ax.set_ylabel(ylabel, rotation=270)

@is_plot_func(creator_type=MultiversePlotCreator, use_dag=False)
def carrying_capacity_fixed_y(dm: DataManager, *,
                              hlpr: PlotHelper,
                              save_kwargs: dict = None,
                              cmap: str = "Spectral",
                              line_plot: bool = True,
                              mov: bool = True,
                              dpi: int = 300,
                              final_state: bool = True,
                              line_plot_args: dict = None,
                              plt_3d: bool = False,
                              # ngrid: tuple = (1000, 1000),
                              xlim: tuple = None,
                              ylim: tuple = None,
                              zlim: tuple = None,
                              **plot_kwargs):
    fig = hlpr.fig
    out_path = hlpr.out_path
    mv = dm["multiverse"]
    config = dm["cfg"]
    for key in config.keys():
        print(key, "  ", config[key].keys())
        print("---")
        try:
            for key2 in config[key].keys():
                print(key2, "  ", config[key][key2])
        except:
            pass
        print("-----------")

    # Get all the configuration options
    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "x")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "y")
    zlabel = get_argument_or_default(plot_kwargs, "zlabel", "mass")
    title = get_argument_or_default(plot_kwargs, "title", 'MuLAN Mac Arthur Model')
    artist = get_argument_or_default(plot_kwargs, "artist", 'Utopia')
    comment = get_argument_or_default(plot_kwargs, "comment", '')
    fps = get_argument_or_default(plot_kwargs, "fps", 15)
    start_frame = get_argument_or_default(plot_kwargs, "start_frame", 0)
    spec_type = get_argument_or_default(plot_kwargs, "spec_type", 1)
    plot_peaks = get_argument_or_default(plot_kwargs, "plot_peaks", False)
    contour_lines = get_argument_or_default(plot_kwargs, "contour_lines", None)
    colorbar = get_argument_or_default(plot_kwargs, "colorbar", False)
    log_z = get_argument_or_default(plot_kwargs, "log_z", False)
    azim = get_argument_or_default(plot_kwargs, "azim", -60)
    elev = get_argument_or_default(plot_kwargs, "elev", None)
    rot = get_argument_or_default(plot_kwargs, "rotation", 0)
    plot_axis = get_argument_or_default(plot_kwargs, "plot_axis", True)
    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "fixed_y")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_masses")

    timestamp = get_argument_or_default(plot_kwargs, "timestamp", True)
    font_size = get_argument_or_default(plot_kwargs, "font_size", 15)

    plt.rcParams.update({'font.size': font_size})

    # Cant plot both
    if plt_3d and plot_peaks:
        return

    # Init 3D axis
    if plt_3d:
        fig.clf()
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(azim=azim)
    else:
        ax = fig.gca()

    if (timestamp):
        timestamp_t = fig.text(0.12, 0.9, "")

    ### MOVIE ###
    if mov:

        metadata = dict(title=title, artist=artist,
                        comment=comment)
        writer = FFMpegWriter(fps=fps, metadata=metadata)

        with writer.saving(fig, out_path.split(".pdf")[0] + ".mp4", dpi):
            time = start_frame
            while True:

                cont, x, y, z, c = get_data_cc_fixed_y(mv,
                                                       x_name=x_name,
                                                       y_name=y_name,
                                                       z_name=z_name,
                                                       spec_type=spec_type,
                                                       time=time)

                # Calculate peaks
                if not plt_3d and plot_peaks:
                    peaks, _ = find_peaks(z)
                    x = x[peaks]
                    y = y[peaks]
                    z = z[peaks]

                # increment frame
                time += 1

                # Break if end of data reached
                if not cont:
                    break

                #
                if colorbar and not plot_peaks:
                    if plt_3d:
                        fig.clf()
                        ax = fig.add_subplot(111, projection="3d")
                        ax.view_init(azim=azim, elev=elev)
                        azim = azim + rot
                    else:
                        fig.clf()
                        ax = fig.add_subplot(111)
                else:
                    ax.clear()
                    if plt_3d:
                        ax.view_init(azim=azim, elev=elev)
                        azim = azim + rot

                # Cosmetics
                if not plot_axis:
                    ax.set_axis_off()

                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                if plt_3d and not colorbar:
                    ax.set_zlabel(zlabel, rotation=270)

                # Set limits
                if xlim is not None:
                    ax.set_xlim(xlim)
                if ylim is not None:
                    ax.set_ylim(ylim)

                try:
                    if plt_3d:
                        if zlim is not None:
                            ax.set_zlim(zlim)
                            z = np.clip(z, zlim[0], zlim[1])
                        c_mappable = plt_3d_frame(ax,
                                                  x,
                                                  y,
                                                  z,
                                                  cmap,
                                                  line_plot,
                                                  line_plot_args)

                        ax.relim()
                        ax.autoscale_view()
                    elif plot_peaks:
                        ax.scatter(x, y, z)
                    else:
                        c_mappable = plt_2d_frame(ax,
                                                  x,
                                                  y,
                                                  z,
                                                  cmap,
                                                  zlim,
                                                  contour_lines,
                                                  log_z)

                    if colorbar and not plot_peaks:
                        cbar = fig.colorbar(c_mappable)

                        cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
                        cbar.ax.set_ylabel(zlabel, rotation=270)
                        if log_z and not plt_3d:
                            tick_locator = ticker.LogLocator(base=10.0, numticks=contour_lines)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                    if (timestamp):
                        dt2 = float(config["model"]["dt2"])
                        # t_max = int(config["num_steps"])
                        write_every = int(config["update"]["parameter_space"]["write_every"])
                        time_s = time * write_every

                        timestamp_t.set_text("t=" + str(time_s * dt2) + " (" + str(time_s) + "dt)")


                except RuntimeError as e:
                    print("RuntimeError")
                    print(e)
                    pass
                except:
                    raise

                fig.tight_layout()
                writer.grab_frame()

    # Plot the final state
    if final_state:

        # Get azim in case of rotation
        azim = get_argument_or_default(plot_kwargs, "azim", -60)

        if colorbar and not plot_peaks:
            if plt_3d:
                fig.clf()
                ax = fig.add_subplot(111, projection="3d")
                ax.view_init(azim=azim, elev=elev)
            else:
                fig.clf()
                ax = fig.add_subplot(111)

            if (timestamp):
                timestamp_t = fig.text(0.12, 0.9, "")

        else:
            ax.clear()
            if plt_3d:
                ax.view_init(azim=azim, elev=elev)

        # Get the data
        cont, x, y, z, c = get_data_cc_fixed_y(mv,
                                               x_name=x_name,
                                               y_name=y_name,
                                               z_name=z_name,
                                               spec_type=spec_type,
                                               time=-1)

        if not plt_3d and plot_peaks:
            peaks, _ = find_peaks(z)
            x = x[peaks]
            y = y[peaks]
            z = z[peaks]

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if plt_3d and not colorbar:
            ax.set_zlabel(zlabel, rotation=270)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        if plt_3d:
            if zlim is not None:
                ax.set_zlim(zlim)
                z = np.clip(z, zlim[0], zlim[1])
            c_mappable = plt_3d_frame(ax, x, y, z, cmap, line_plot, line_plot_args)
            ax.relim()
            ax.autoscale_view()
        elif plot_peaks:
            ax.scatter(x, y, z)
        else:
            c_mappable = plt_2d_frame(ax, x, y, z, cmap, zlim, contour_lines, log_z)

        if colorbar and not plot_peaks:
            cbar = fig.colorbar(c_mappable)
            cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
            cbar.ax.set_ylabel(zlabel, rotation=270)
            if log_z and not plt_3d:
                tick_locator = ticker.LogLocator(base=10.0, numticks=contour_lines)
                cbar.locator = tick_locator
                cbar.update_ticks()

        if (timestamp):
            dt2 = float(config["model"]["dt2"])
            t_max = int(config["update"]["parameter_space"]["num_steps"])
            # write_every = int(config["write_every"])
            time_s = t_max

            timestamp_t.set_text("t=" + str(time_s * dt2) + " (" + str(time_s) + "dt)")


        # Save and close the figure
        # save_and_close(out_path, save_kwargs=save_kwargs)


@is_plot_func(creator_type=UniversePlotCreator)
def carrying_capacity_uni(dm: DataManager, *,
                          uni: UniverseGroup,
                          hlpr: PlotHelper,
                          cmap: str = "Spectral",
                          mov: bool = True,
                          dpi: int = 300,
                          final_state: bool = True,
                          plt_3d: bool = False,
                          # ngrid: tuple = (1000, 1000),
                          xlim: tuple = None,
                          ylim: tuple = None,
                          zlim: tuple = None,
                          contour_lines: int = None,
                          **plot_kwargs):
    fig = hlpr.fig
    out_path = hlpr.out_path
    config = uni["cfg"]

    xlabel = get_argument_or_default(plot_kwargs, "xlabel", "x")
    ylabel = get_argument_or_default(plot_kwargs, "ylabel", "y")
    zlabel = get_argument_or_default(plot_kwargs, "zlabel", "Biomass")
    title = get_argument_or_default(plot_kwargs, "title", 'MuLAN Mac Arthur Model')
    artist = get_argument_or_default(plot_kwargs, "artist", 'Utopia')
    comment = get_argument_or_default(plot_kwargs, "comment", '')
    fps = get_argument_or_default(plot_kwargs, "fps", 15)
    start_frame = get_argument_or_default(plot_kwargs, "start_frame", 0)
    spec_type = get_argument_or_default(plot_kwargs, "spec_type", 1)
    plot_peaks = get_argument_or_default(plot_kwargs, "plot_peaks", False)
    contour_lines = get_argument_or_default(plot_kwargs, "contour_lines", None)
    colorbar = get_argument_or_default(plot_kwargs, "colorbar", False)
    log_z = get_argument_or_default(plot_kwargs, "log_z", False)
    azim = get_argument_or_default(plot_kwargs, "azim", -60)
    elev = get_argument_or_default(plot_kwargs, "elev", None)
    rot = get_argument_or_default(plot_kwargs, "rotation", 0)
    plot_axis = get_argument_or_default(plot_kwargs, "plot_axis", True)

    x_name = get_argument_or_default(plot_kwargs, "x_name", "_traits")
    y_name = get_argument_or_default(plot_kwargs, "y_name", "_niche_w")
    z_name = get_argument_or_default(plot_kwargs, "z_name", "_masses")

    time_frame_l = get_argument_or_default(plot_kwargs, "timeframes", [])

    timestamp = get_argument_or_default(plot_kwargs, "timestamp", True)
    font_size = get_argument_or_default(plot_kwargs, "font_size", 15)

    plt.rcParams.update({'font.size': font_size})

    if plt_3d and plot_peaks:
        return

    if plt_3d:
        fig.clf()
        ax = fig.add_subplot(111, projection="3d")
        ax.view_init(azim=azim)

    else:
        ax = fig.gca()

    if (timestamp):
        timestamp_t = fig.text(0.12, 0.9, "")

    grp_main = uni['data/MuLAN_MA']

    if mov:

        metadata = dict(title='MuLAN Mac Arthur Model', artist='Utopia',
                        comment='')
        writer = FFMpegWriter(fps=15, metadata=metadata)

        with writer.saving(fig, out_path.split(".pdf")[0] + ".mp4", dpi):
            time = start_frame
            while True:

                cont, x, y, z = get_data_cc_uni(grp_main,
                                                spec_type=spec_type,
                                                x_name=x_name,
                                                y_name=y_name,
                                                z_name=z_name,
                                                time=time)

                if not plt_3d and plot_peaks:
                    peaks, _ = find_peaks(z)
                    x = x[peaks]
                    y = y[peaks]
                    z = z[peaks]

                time += 1

                if not cont:
                    break

                if x.size < 3:
                    continue

                if colorbar and not plot_peaks:
                    if plt_3d:
                        fig.clf()
                        ax = fig.add_subplot(111, projection="3d")
                        ax.view_init(azim=azim, elev=elev)
                        azim = azim + rot
                    else:
                        fig.clf()
                        ax = fig.add_subplot(111)
                else:
                    ax.clear()
                    if plt_3d:
                        ax.view_init(azim=azim, elev=elev)
                        azim = azim + rot

                if not plot_axis:
                    ax.set_axis_off()

                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                if plt_3d and not colorbar:
                    ax.set_zlabel(zlabel, rotation=270)

                if xlim is not None:
                    ax.set_xlim(xlim)

                if ylim is not None:
                    ax.set_ylim(ylim)

                try:
                    if plt_3d:
                        if zlim is not None:
                            ax.set_zlim(zlim)
                            z = np.clip(z, zlim[0], zlim[1])
                        c_mappable = plt_3d_frame(ax,
                                                  x,
                                                  y,
                                                  z,
                                                  cmap,
                                                  False,
                                                  {})

                        ax.relim()
                        ax.autoscale_view()
                    elif plot_peaks:
                        ax.scatter(x, y, z)
                    else:
                        c_mappable = plt_2d_frame(ax,
                                                  x,
                                                  y,
                                                  z,
                                                  cmap,
                                                  zlim,
                                                  contour_lines,
                                                  log_z)

                    if colorbar and not plot_peaks:
                        cbar = fig.colorbar(c_mappable)
                        cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
                        cbar.ax.set_ylabel(zlabel, rotation=270)
                        if log_z and not plt_3d:
                            tick_locator = ticker.LogLocator(base=10.0, numticks=contour_lines)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                    if (timestamp):
                        dt2 = float(config["MuLAN_MA"]["dt2"])
                        # t_max = int(config["num_steps"])
                        write_every = int(config["write_every"])
                        time_s = time * write_every

                        timestamp_t.set_text("t=" + str(time_s * dt2) + " (" + str(time_s) + "dt)")


                except RuntimeError as e:
                    print(e)
                    pass
                except:
                    raise

                writer.grab_frame()

    if final_state:
        azim = get_argument_or_default(plot_kwargs, "azim", -60)

        if colorbar and not plot_peaks:
            if plt_3d:
                fig.clf()
                ax = fig.add_subplot(111, projection="3d")
                ax.view_init(azim=azim, elev=elev)
            else:
                fig.clf()
                ax = fig.add_subplot(111)

            if (timestamp):
                timestamp_t = fig.text(0.12, 0.9, "")
        else:
            ax.clear()
            if plt_3d:
                ax.view_init(azim=azim, elev=elev)

        time_f = -1
        if len(time_frame_l) > 0:
            time_f = time_frame_l[0]

        cont, x, y, z = get_data_cc_uni(grp_main,
                                        spec_type=spec_type,
                                        x_name=x_name,
                                        y_name=y_name,
                                        z_name=z_name,
                                        time=time_f)

        if not plt_3d and plot_peaks:
            peaks, _ = find_peaks(z)
            x = x[peaks]
            y = y[peaks]
            z = z[peaks]

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if plt_3d and not colorbar:
            ax.set_zlabel(zlabel, rotation=270)

        if xlim is not None:
            ax.set_xlim(xlim)

        if ylim is not None:
            ax.set_ylim(ylim)

        if plt_3d:
            if zlim is not None:
                ax.set_zlim(zlim)
                z = np.clip(z, zlim[0], zlim[1])
            c_mappable = plt_3d_frame(ax,
                                      x,
                                      y,
                                      z,
                                      cmap,
                                      False,
                                      {})

            ax.relim()
            ax.autoscale_view()
        elif plot_peaks:
            ax.scatter(x, y, z)
        else:
            c_mappable = plt_2d_frame(ax,
                                      x,
                                      y,
                                      z,
                                      cmap,
                                      zlim,
                                      contour_lines,
                                      log_z)

        if colorbar and not plot_peaks:
            cbar = fig.colorbar(c_mappable)
            cbar.ax.get_yaxis().labelpad = 10 + font_size/2.0
            cbar.ax.set_ylabel(zlabel, rotation=270)
            if log_z and not plt_3d:
                tick_locator = ticker.LogLocator(base=10.0, numticks=contour_lines)
                cbar.locator = tick_locator
                cbar.update_ticks()

        if (timestamp):
            dt2 = float(config["MuLAN_MA"]["dt2"])
            t_max = int(config["num_steps"])
            # write_every = int(config["write_every"])
            time_s = t_max

            timestamp_t.set_text("t=" + str(time_s * dt2) + " (" + str(time_s) + "dt)")


        # Save and close the figure
        # save_and_close(out_path, save_kwargs=save_kwargs)
