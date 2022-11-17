import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation, LinearTriInterpolator
from matplotlib import ticker
from scipy.interpolate import griddata

from paramspace.paramdim import ParamDim

import numpy as np
import networkx as nx

from utopya import DataManager, UniverseGroup, MultiverseGroup


def save_and_close(out_path: str, *, save_kwargs: dict=None) -> None:
    """Save and close the figure, passing the kwargs

    Args:
        out_path (str): The output path to save the figure to
        save_kwargs (dict, optional): The additional save_kwargs
    """
    plt.savefig(out_path, **(save_kwargs if save_kwargs else {}))
    plt.close()


def print_cfg(cfg, depth=0):
    """Print a configuration group

    Args:
        cfg : the configuration group
        depth (int, optional): Depth in config tree. Defaults to 0.
    """
    if type(cfg) == type("str"):
        return
    for i in cfg:
        print(depth*"    " + i)
        try:
            print_cfg(cfg[i], depth=depth+1)
        except:
            pass

def save_np_data(*data, path=".", filename="data.npy"):
    """Save Data to npy binary file

    Args:
        data : arbitrary number of nump array
        path (str, optional): Save path. Defaults to ".".
        filename (str, optional): Filename. Defaults to "data.npy".
    """
    with open("/" + path.strip("/") + "/" + filename, "wb") as f:
        for d in data:
            np.save(f, d)


def get_argument_or_default(kwargs, argument, default):
    """Get Entries from Plot configuration if availabel or use default

    Args:
        kwargs (dict): The config dictionary
        argument (str): The name if the argument to assign
        default (any): Default value

    Returns:
        type(default) : The Argument
    """
    v = kwargs.get(argument)
    print(argument)
    print(v)
    print(type(v))
    if type(v) == ParamDim:
        v = v.default.value

    if type(v) == dict:
        v = v["default"]
    
    var = v if v is not None else default
    if default != None:
        return type(default)(var)
    else:
        return var

def get_time_grps(grp):
    """Get names of numbered data groups

    Args:
        grp : The parent group that contains the data groups for

    Returns:
        list: sorted list of all data group names (as ints)
    """
    time_list = []

    for t in grp:
        time_list.append(int(t))

    return sorted(time_list)

def get_maximums(grp, time_list, attx="_traits", atty="_masses", attz="_niche_w"):

    """Get the limits of a series of datasets

    Args:
        grp : The parent group that contains the data groups for
        time_list (list) : The names of the data groups 
            as list of ints or str 
        attx (str) : Name of 1st Dataset
        atty (str) : Name of 2nd Dataset

    Raises:
        ValueError: Raised if values are unreasonable
            Can happen if the data is empty or contains singularities

    Returns:
        (list, list): The Minimum and Maximum of attx and atty Datasets
    """

    mmx = [np.inf, -np.inf]
    mmy = [np.inf, -np.inf]
    mmz = [np.inf, -np.inf]


    for t in time_list:
        frame = grp[str(t)]

        arrx = frame[attx]
        arry = frame[atty]
        arrz = frame[atty]

        mmx[0] = np.minimum(mmx[0], np.amin(arrx))
        mmx[1] = np.maximum(mmx[1], np.amax(arrx))

        mmy[0] = np.minimum(mmy[0], np.amin(arry))
        mmy[1] = np.maximum(mmy[1], np.amax(arry))

        mmz[0] = np.minimum(mmz[0], np.amin(arrz))
        mmz[1] = np.maximum(mmz[1], np.amax(arrz))

    for i in range(2):
        if np.isinf(mmx[i]) or np.isinf(mmy[i]) or np.isinf(mmz[i]):
            raise ValueError("Could not calculate limits")

    return mmx, mmy, mmz


def get_data_cc_fixed_y(mv: MultiverseGroup,
                        x_name="_traits",
                        y_name="fixed_y",
                        z_name="_masses",
                        spec_type: int=1,
                        time: int=-1):
    
    y_val = []
    for i in mv:
        y_val.append(mv[i]["cfg"]["MuLAN_MA"][y_name])

    x = []
    y = []
    z = []
    c = []

    k = 0
    for i in mv:
        grp_main = mv[i]["data"]["MuLAN_MA"]

        iter = get_time_grps(grp_main)

        if len(iter) <= time:
            return False, np.array(x), np.array(y), np.array(z), np.array(c)

        grp = grp_main[str(iter[time])]

        num_species = grp.attrs["num_species"]

        types = grp["_types"].values.reshape(num_species)
        traits = grp[x_name].values.reshape(num_species)
        masses = grp[z_name].values.reshape(num_species)

        c_traits = np.array(traits[types == spec_type])
        c_masses = np.array(masses[types == spec_type])

        indi = c_traits.argsort()
        c_traits = np.array(c_traits[indi[::-1]])
        c_masses = np.array(c_masses[indi[::-1]])

        x.append(np.append(c_traits, np.nan))
        z.append(np.append(c_masses, np.nan))
        y.append(np.append(np.full_like(c_traits, y_val[k]), np.nan))
        c.append(np.append(np.full_like(c_traits, k), np.nan))

        k = k + 1

    return True, \
        np.concatenate(x), \
        np.concatenate(y), \
        np.concatenate(z), \
        np.concatenate(c)


def get_data_cc_uni(grp_main, x_name="_traits", y_name="_niche_w", z_name="_masses", spec_type: int=1, time: int = -1):

    """

    :param grp_main:
    :param x_name:
    :param y_name:
    :param z_name:
    :param spec_type:
    :param time:
    :return:
    """

    iter = get_time_grps(grp_main)

    if len(iter) <= time:
        return False, np.array([]), np.array([]), np.array([])

    grp = grp_main[str(iter[time])]

    num_species = grp.attrs["num_species"]

    types = grp["_types"].values.reshape(num_species)
    traits = grp[x_name].values.reshape(num_species)
    nw = grp[y_name].values.reshape(num_species)
    masses = grp[z_name].values.reshape(num_species)

    x = np.array(traits[types == spec_type])
    y = np.array(nw[types == spec_type])
    z = np.array(masses[types == spec_type])

    return True, x, y, z

def plt_3d_frame(ax, x, y, z, cmap, line_plot : bool, line_plot_args):

    """ Performs triangulation and plots Data in 3D ax

    Args:
        ax : Matplotlib axes(projection="3d")
        x : x Data (1D)
        y : y Data (1D)
        z : z Data (1D)
        cmap : colormap
        line_plot (bool) : Plot contour lines?
        line_plot_args : arguments for line plot (from config)

    """

    triang = Triangulation(x[~np.isnan(x)], y[~np.isnan(y)])

    if cmap == "None":
        c_mappable = ax.plot_trisurf(triang, z[~np.isnan(z)], linewidth=0.2, antialiased=True)
    else:
        c_mappable = ax.plot_trisurf(triang, z[~np.isnan(z)], cmap=cmap, linewidth=0.2, antialiased=True)

    if line_plot:
        ax.plot(x, y, z, **line_plot_args)
    
    return c_mappable


def plt_2d_frame(ax, x, y, z, cmap, zlim, contour_lines, log):

    '''
    xi = np.linspace(xlim[0], xlim[1], ngrid[0])
    yi = np.linspace(ylim[0], ylim[1], ngrid[1])

    # Perform linear interpolation of the data (x,y)
    # on a grid defined by (xi,yi)
    triang = Triangulation(x[~np.isnan(x)], y[~np.isnan(y)])
    interpolator = LinearTriInterpolator(triang, z[~np.isnan(z)])
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    ax.tricontour(Xi, Yi, zi, cmap=cmap)
    '''

    if log:
        mask = np.all([~np.isnan(x), z != 0], axis=0)
        x = x[mask]
        y = y[mask]
        z = z[mask]

    if contour_lines is not None and zlim is None:
        zlim = (np.nanmin(z), np.nanmax(z))

    if contour_lines is None and zlim is not None:
        contour_lines = 11

    if zlim is None:
        if log:

            return ax.tricontourf(x[~np.isnan(x)], y[~np.isnan(y)], z[~np.isnan(x)], 
                    cmap=cmap, locator=ticker.LogLocator())
        else:
            return ax.tricontourf(x[~np.isnan(x)], y[~np.isnan(y)], z[~np.isnan(x)], cmap=cmap)
    else:
        if log:

            # TODO: Find better solution
            if zlim[0] == 0:
                zlim[0] = 0.000001
            return ax.tricontourf(x[~np.isnan(x)], y[~np.isnan(y)], z[~np.isnan(x)],
                        levels=np.power(
                            10.0, 
                            np.linspace(np.log10(zlim[0]), np.log10(zlim[1]), 
                            contour_lines)),
                        cmap=cmap, locator=ticker.LogLocator())
        else:
            return ax.tricontourf(x[~np.isnan(x)], y[~np.isnan(y)], z[~np.isnan(x)],
                       levels=np.linspace(zlim[0], zlim[1], contour_lines), cmap=cmap)


def get_graph_from_frame(frame, w_name: str, m_name: str ):

    vertices = frame["_graph_vertex"]
    edges = frame["_graph_edges"]
    weights = frame[w_name]
    masses = frame[m_name]

    edge_l = []
    weight_l = []
    for i, j, k in zip(*edges, weights):
        edge_l.append((int(i), int(j), float(1)))

    node_l = []
    for i, j in zip(vertices, masses):
        v_i = int(i.data)
        node_l.append((v_i, {m_name: float(j)}))


    g = nx.DiGraph()
    g.add_weighted_edges_from(edge_l)
    g.add_nodes_from(node_l)

    return g