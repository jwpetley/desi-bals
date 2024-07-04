import numpy as np   
import paths
import matplotlib.pyplot as plt

from astropy.table import Table


if __name__ == "__main__":
    bals = Table.read(paths.data / "matched_bals.fits")
    bals_ai = Table.read(paths.data / "matched_ai.fits")
    qsos = Table.read(paths.data / "matched_qsos.fits")


    plt.figure()

    plt.scatter(bals['Total_flux'], bals['BI_CIV'], label="BI BALs", s=1)

    plt.xlabel(r"$\log_{10} F_{\nu}$")
    plt.ylabel(r"$BI_{CIV}$")
    plt.legend()
    plt.savefig(paths.figures / "flux_bi.pdf", bbox_inches = "tight",
                dpi = 300)