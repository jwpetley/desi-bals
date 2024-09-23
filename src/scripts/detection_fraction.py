import matplotlib.pyplot as plt
import numpy as np
import paths

from astropy.table import Table



def remove_bals(data, bals):
    '''Remove any tagetid overlaps'''
    ids = bals["TARGETID"]
    removed = data[~np.isin(data["TARGETID"], ids)]
    print("Removed ", len(data) - len(removed), " targets")
    return removed


def radio_detected(data):
    return data[data["Total_flux"] >0.01]



def flux_cut(data, min, max):
    return data[(data["FLUX_G"] > min) & (data["FLUX_G"] < max)]

def redshift_cut(data, min, max):
    return data[(data["Z"] > min) & (data["Z"] < max)]

def BI_cut(data, min, max):
    return data[(data["BI_CIV"] > min) & (data["BI_CIV"] < max)]

def AI_cut(data, min, max):
    return data[(data["AI_CIV"] > min) & (data["AI_CIV"] < max)]


def detection_fraction(data, bins):
    '''Calculate the detection fraction as a function of redshift'''
    hist, _ = np.histogram(data["Z"], bins=bins)
    detected = radio_detected(data)
    detected_hist, _ = np.histogram(detected["Z"], bins=bins)
    return detected_hist/hist

def detection_fraction_err(data, bins):
    '''Calculate the detection fraction as a function of redshift'''
    hist, _ = np.histogram(data["Z"], bins=bins)
    detected = radio_detected(data)
    detected_hist, _ = np.histogram(detected["Z"], bins=bins)
    return detected_hist/hist, np.sqrt(detected_hist)/hist

def plot_detection_fraction(bals, ai_bals, qsos, fname):
    '''Create detection fraction figure'''

    bins = np.linspace(1.5, 4, 6)

    bals_frac, bals_frac_err = detection_fraction_err(bals, bins)
    ai_bals_frac, ai_bals_frac_err = detection_fraction_err(ai_bals, bins)
    qsos_frac, qsos_frac_err = detection_fraction_err(qsos, bins)

    # Two panel figure - Lower panel should have relative fractions
    fig , axs = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    axs[0].errorbar(bins[:-1], bals_frac, yerr=bals_frac_err, label="BI BALs",
                capsize=4)
    axs[0].errorbar(bins[:-1], ai_bals_frac, yerr=ai_bals_frac_err, label="AI BALs",
                capsize = 4)
    axs[0].errorbar(bins[:-1], qsos_frac, yerr=qsos_frac_err, label="QSOs",
                capsize = 4)
    axs[1].set_xlabel(r"$z$")
    axs[0].set_ylabel("LoTSS Detection fraction")
    axs[0].legend()

    axs[1].errorbar(bins[:-1], bals_frac/qsos_frac, yerr=bals_frac_err/qsos_frac,
                label="BI/QSOs", capsize = 4)
    axs[1].errorbar(bins[:-1], ai_bals_frac/qsos_frac, yerr=ai_bals_frac_err/qsos_frac,
                label="AI/QSOs", capsize = 4)
    
    axs[1].set_ylabel("Relative Fraction")
    axs[1].legend()
    plt.savefig(paths.figures / fname, bbox_inches = "tight",
                dpi = 300)
    

if __name__ == "__main__":
    # Read the data
    bals = paths.data / "desi_bals_radio.fits"
    qsos = paths.data / "desi_qso_radio.fits"

    bals = Table.read(bals)
    qsos = Table.read(qsos)

    #Define samples
    ai_bals = AI_cut(bals, 1, 20000)
    bals = BI_cut(bals, 1, 20000)
    qsos = remove_bals(qsos, bals)
    qsos = remove_bals(qsos, ai_bals)
    ai_bals = remove_bals(ai_bals, bals)
    bals = remove_bals(bals, ai_bals)

    #Redshift cut
    ai_bals = redshift_cut(ai_bals, 1.5, 5)
    bals = redshift_cut(bals, 1.5, 5)
    qsos = redshift_cut(qsos, 1.5, 5)



    plot_detection_fraction(bals, ai_bals, qsos, "detection_fraction_z.pdf")