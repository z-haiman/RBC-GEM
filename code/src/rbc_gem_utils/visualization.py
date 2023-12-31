"""
Contains functions to help visualize updates to the RBC-GEM content.
"""
import matplotlib as mpl
import seaborn as sns
import numpy as np

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
    cmap = mpl.colors.LinearSegmentedColormap.from_list('compare', colors, len(colors))
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
    ax.set_yticks(ticks=np.linspace(0, len(df_comparison.T), len(df_comparison.T)+1), minor=True)
    ax.grid(which="minor", axis="y", color="k", linewidth="2",linestyle="-", zorder=1)
    # Draw a frame 
    for _, spine in ax.spines.items(): 
        spine.set_visible(True) 
        spine.set_linewidth(2) 

    # Colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.dividers.set_color('black')
    colorbar.dividers.set_linewidth(2)
    colorbar.outline.set_color("black")
    colorbar.outline.set_linewidth(2)
    colorbar.set_ticks(list(np.array(ints) + 0.5), labels=values)
    return ax