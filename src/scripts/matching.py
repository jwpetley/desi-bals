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

def ergs(six_micron_lum):
    return np.log10(10**six_micron_lum*(3e8/6e-6)*1e7)


def six_micron(data):
    wise_wavs = np.array([4.6e-6,12e-6])
    zero_points = np.array([171.787, 31.674])
    ab_conversions = np.array([3.339, 5.174])
    
    luminosities = []
    
    for source in tqdm(data):
        mags = np.array([float(source['w2mpro']),
                         float(source['w3mpro'])], dtype = np.float64)
        
            
        ab_mags = mags + ab_conversions
        
        
        z = source['Z']
        
        lum_dist = Planck15.luminosity_distance(z).to(u.meter)
        
        
        fluxes = 10**((ab_mags - 8.926)/-2.5)*1e-26
        
        target_wav = 6e-6*(1+z)
        
        interp_points = np.log10(fluxes)
        
        
        #print(interp_points)
        
        line = interp1d( np.log10(wise_wavs), interp_points, fill_value="extrapolate", kind = 'linear' )
        
        
        interp_flux = 10**line(np.log10(target_wav))
        
        #print(interp_flux)
        
        # I have flux in SI units at this point
        
        # Consult https://ned.ipac.caltech.edu/level5/Sept02/Hogg/Hogg2.html (eq. 6 + 7)
        
        
        luminosity = (interp_flux*4*np.pi*(lum_dist**2))/(1+z) #What am I doing here??

        
        luminosities.append(np.log10(luminosity.value))

        # Change nan to zero
        if np.isnan(luminosities[-1]):
            luminosities[-1] = -99
       
        
    return ergs(np.array(luminosities))

def wise_detected(data):
    data = data[data["w3mpro"]!= "null"]
    data = data[data["w2mpro"]!= "null"]
    return data

def z_six_micron_match(data1, data2):
    """Use  kNN to match two datasets by redshift and 6 micron luminosity
    
    Return two tables with a new column "match" containing the TARGETID of the other dataset"""
    data1 = data1.copy()
    data2 = data2.copy()
    
    X1 = np.array([data1["Z"], data1["six_micron"]]).T
    X2 = np.array([data2["Z"], data2["six_micron"]]).T
    
    nbrs = NearestNeighbors(n_neighbors=1, algorithm="ball_tree", radius=0.2).fit(X2)
    distances, indices = nbrs.kneighbors(X1)
    
    data1["match"] = data2["TARGETID"][indices]
    
    matched = data2[indices.flatten()]

    return matched


if __name__ == "__main__":
    bals = Table.read(paths.data / "desi_bals_radio.fits")
    qsos = Table.read(paths.data / "desi_qso_radio.fits")

    bals = wise_detected(bals)
    qsos = wise_detected(qsos)

    
    
    ai_bals = AI_cut(bals, 1, 200000)
    bals = BI_cut(bals, 1, 200000)
    qsos = remove_bals(qsos, bals)
    qsos = remove_bals(qsos, ai_bals)
    ai_bals = remove_bals(ai_bals, bals)
    bal = remove_bals(bals, ai_bals)

    bals = redshift_cut(bals, 1.5, 4)
    ai_bals = redshift_cut(ai_bals, 1.5, 4)
    qsos = redshift_cut(qsos, 1.5, 4)

    bals["six_micron"] = six_micron(bals)
    qsos["six_micron"] = six_micron(qsos)
    ai_bals["six_micron"] = six_micron(ai_bals)

    bals = bals[bals["six_micron"] > 0]
    qsos = qsos[qsos["six_micron"] > 0]
    ai_bals = ai_bals[ai_bals["six_micron"] > 0]

    matched = z_six_micron_match(bals, qsos)

    matched_ai = z_six_micron_match(bals, ai_bals)

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
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    axs = axs.flatten()
    handles = []

    x = np.array(bals["Z"]).byteswap().newbyteorder()
    y = np.array(bals["six_micron"]).byteswap().newbyteorder()
    sns.kdeplot(x = x, y = y, label="BI BALs", ax=axs[0],
                color = "tab:orange")
    handles.append(mlines.Line2D([], [], color='tab:orange', label='BI BALs'))

    sns.kdeplot(x = np.array(qsos["Z"]).byteswap().newbyteorder(), y = np.array(qsos["six_micron"]).byteswap().newbyteorder(), label="QSOs", ax=axs[0],
                color = "tab:blue")
    handles.append(mlines.Line2D([], [], color='tab:blue', label='QSOs'))

    axs[0].set_xlabel(r"$z$")
    axs[0].set_ylabel(r"$\log_{10} L_{6 \mu m}$")

    axs[0].legend(handles=handles)

    handles = []

    sns.kdeplot(x = np.array(bals["Z"]).byteswap().newbyteorder(), y = np.array(bals["six_micron"]).byteswap().newbyteorder(), label="BI BALs", ax=axs[1],
                color = "tab:orange")
    handles.append(mlines.Line2D([], [], color='tab:orange', label='BI BALs'))
    
    sns.kdeplot(x = np.array(matched["Z"]).byteswap().newbyteorder(), y = np.array(matched["six_micron"]).byteswap().newbyteorder(), label="Matched QSOs",ax=axs[1],
                color = "tab:green")
    handles.append(mlines.Line2D([], [], color='tab:green', label='Matched QSOs'))
    
    axs[1].set_xlabel(r"$z$")
    axs[1].set_ylabel(r"$\log_{10} L_{6 \mu m}$")

    axs[1].legend(handles=handles)

    plt.savefig(paths.figures / "six_micron_z.pdf", bbox_inches = "tight",
                dpi = 300)

    plot_detection_fraction(bals, matched_ai, matched, "matched_detection_fraction.pdf")

    matched.write(paths.data / "matched/matched_qsos.fits", overwrite=True)
    matched_ai.write(paths.data / "matched/matched_ai.fits", overwrite=True)
    bals.write(paths.data / "matched/matched_bals.fits", overwrite=True)

    


