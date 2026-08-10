from matplotlib.colors import LinearSegmentedColormap


CTAO_COLORS = {
    # Primary colors
    "cherenkov_blue": "#00004A",
    "cherenkov_cyan": "#00E4D8",
    "light_gray": "#F5F5F5",
    # Secondary colors
    "white": "#FFFFFF",
    "black": "#000000",
    "cosmic_azure": "#007AFF",
    "interstellar_indigo": "#00009C",
}

CTAO_CMAP = LinearSegmentedColormap.from_list(
    "ctao_hist",
    [
        (0.00, CTAO_COLORS["black"]),
        (0.15, CTAO_COLORS["cherenkov_blue"]),
        (0.70, CTAO_COLORS["cherenkov_cyan"]),
        (1.00, CTAO_COLORS["white"]),
    ],
    N=256,
)

CTAO_CMAP_R = CTAO_CMAP.reversed(name="ctao_hist_r")
