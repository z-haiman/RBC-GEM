"""
Contains functions to help visualize updates to the RBC-GEM content.
"""
import matplotlib as mpl
import numpy as np
import seaborn as sns


def visualize_comparison(df_comparison):
    if df_comparison.empty:
        raise ValueError("`df_comparison` should not be empty")

    # TODO update to have one location for values and colors
    values_to_colors = {
        "EMPTY": "xkcd:white",
        "NO CHANGE": "xkcd:light blue",
        "REMOVED": "xkcd:red",
        "CHANGED": "xkcd:yellow",
        "ADDED": "xkcd:light green",
    }
    values = list(values_to_colors.keys())
    ints = range(len(values_to_colors))
    colors = list(values_to_colors.values())
    cmap = mpl.colors.LinearSegmentedColormap.from_list("compare", colors, len(colors))
    norm = mpl.colors.BoundaryNorm(np.linspace(0, len(colors), len(colors) + 1), cmap.N)
    mapping = dict(zip(values, ints))
    df_comparison = df_comparison.map(lambda x: mapping[x])

    ax = sns.heatmap(
        df_comparison.T,
        xticklabels=[],
        cmap=cmap,
        norm=norm,
        cbar=True,
    )
    ax.set_yticks(
        ticks=np.linspace(0, len(df_comparison.T), len(df_comparison.T) + 1), minor=True
    )
    ax.grid(which="minor", axis="y", color="k", linewidth="2", linestyle="-", zorder=1)
    # Draw a frame
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(2)

    # Colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.dividers.set_color("black")
    colorbar.dividers.set_linewidth(2)
    colorbar.outline.set_color("black")
    colorbar.outline.set_linewidth(2)
    colorbar.set_ticks(list(np.array(ints) + 0.5), labels=values)
    return ax


def cmap_map(function, cmap):
    """Applies function (which should operate on vectors of shape 3: [r, g, b]), on colormap cmap.
    This routine will break any discontinuous points in a colormap.

    Taken from scipy cookbook: https://scipy-cookbook.readthedocs.io/items/Matplotlib_ColormapTransformations.html
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ("red", "green", "blue"):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step: np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(["red", "green", "blue"]):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j, i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1],), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector
    return mpl.colors.LinearSegmentedColormap("colormap", cdict, 1024)
