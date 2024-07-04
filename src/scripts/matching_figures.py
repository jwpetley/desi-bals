import numpy as np  
import paths
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import seaborn as sns

from tqdm import tqdm
from astropy.table import Table
from sklearn.neighbors import NearestNeighbors
from astropy.cosmology import Planck15
from scipy.interpolate import interp1d
from detection_fraction import AI_cut, BI_cut, flux_cut, redshift_cut, remove_bals, plot_detection_fraction
from matching import wise_detected, six_micron


if __name__ == "__main__":
    matched = Table.read(paths.data / "matched/matched_qsos.fits")
    bals = Table.read(paths.data / "matched/matched_bals.fits")
    matched_ai = Table.read(paths.data / "matched/matched_ai.fits")

    qsos = Table.read(paths.data / "matched_qsos.fits")
    ai_bals = AI_cut(bals, 1, 200000)

    qsos = wise_detected(qsos)

    qsos = remove_bals(qsos, bals)
    qsos = remove_bals(qsos, ai_bals)
    ai_bals = remove_bals(ai_bals, bals)
    bal = remove_bals(bals, ai_bals)

    qsos = redshift_cut(qsos, 1.5, 4)

    qsos["six_micron"] = six_micron(qsos)

    qsos = qsos[qsos["six_micron"] > 0]

    print(len(matched), len(bals), len(qsos))

    plt.figure()

    plt.scatter(bals["Z"], bals["six_micron"], label="BI BALs", s=1)
    plt.scatter(qsos["Z"], qsos["six_micron"], label="QSOs", s=1)
    plt.scatter(matched["Z"], matched["six_micron"], label="Matched QSOs", s=1)

    plt.xlabel(r"$z$")
    plt.ylabel(r"$\log_{10} L_{6 \mu m}$")
    plt.legend()
    plt.close()

    #Same as a kde plot
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
    axs = axs.flatten()
    
    axs[0].hist2d(bals["Z"], bals["six_micron"], bins=50, cmap="viridis")

    axs[1].hist2d(qsos["Z"], qsos["six_micron"], bins=50, cmap="viridis")

    axs[2].hist2d(matched["Z"], matched["six_micron"], bins=50, cmap="viridis")

    axs[3].hist2d(matched_ai["Z"], matched_ai["six_micron"], bins=50, cmap="viridis")

    axs[0].set_xlabel(r"$z$")
    axs[0].set_ylabel(r"$\log_{10} L_{6 \mu m}$")
    axs[1].set_xlabel(r"$z$")
    axs[2].set_xlabel(r"$z$")
    axs[1].set_ylabel(r"$\log_{10} L_{6 \mu m}$")
    axs[0].set_title("BI BALs")
    axs[1].set_title("QSOs")
    axs[2].set_title("Matched QSOs")

    plt.tight_layout()

    plt.savefig(paths.figures / "six_micron_z.pdf", bbox_inches = "tight",
                dpi = 300)

    plot_detection_fraction(bals, matched_ai, matched, "matched_detection_fraction.pdf")