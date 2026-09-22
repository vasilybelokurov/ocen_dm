"""Scientific comparison figures from saved dynamical diagnostics."""
from pathlib import Path

import numpy as np


def plot_evolutions(directories, output):
    from astropy.table import Table
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output = Path(output)
    if output.exists():
        raise FileExistsError(f"Figure already exists: {output}")
    figure, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    fields = ("dm_20pc_mass_msun", "dm_20pc_rho_msun_pc3", "dm_bound_mass_msun")
    labels = (r"$M_{\rm DM}(<20\,{\rm pc})\ [M_\odot]$",
              r"$\rho_{\rm DM}(20\,{\rm pc})\ [M_\odot\,{\rm pc}^{-3}]$",
              r"DM bound mass proxy $[M_\odot]$")
    try:
        for directory in map(Path, directories):
            table = Table.read(directory/"diagnostics.ecsv")
            label = directory.parent.name+" / "+directory.name.removeprefix("evolution_")
            for ax, field in zip(axes, fields):
                if field not in table.colnames:
                    raise ValueError("This figure requires 20 pc in the diagnostic radii")
                ax.plot(table["elapsed_myr"], table[field], label=label)
            low = np.asarray(table["dm_20pc_shell_neff"]) < 100
            if np.any(low):
                axes[1].scatter(table["elapsed_myr"][low], table[fields[1]][low],
                                marker="x", s=20, color=axes[1].lines[-1].get_color())
        for ax, label in zip(axes, labels):
            ax.set(xlabel="Elapsed time [Myr]", ylabel=label)
            ax.grid(alpha=.2)
        axes[0].legend(fontsize=7)
        axes[1].set_title(r"Shell estimate; × marks $N_{\mathrm{eff}}<100$", fontsize=9)
        figure.tight_layout()
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=160)
    finally:
        plt.close(figure)
    return output
