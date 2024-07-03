import numpy as np  
import paths
import astropy.units as u

from tqdm import tqdm
from astropy.table import Table
from sklearn.neighbors import NearestNeighbors
from astropy.cosmology import Planck15
from scipy.interpolate import interp1d

def six_micron(data):
    wise_wavs = np.array([4.6e-6,12e-6])
    zero_points = np.array([171.787, 31.674])
    ab_conversions = np.array([3.339, 5.174])
    
    luminosities = []
    
    for source in tqdm(data):
        if str(source['w2mpro']) == 'null':
            
            source['w2mpro'] = 99.0
        if str(source['w3mpro']) == 'null':
            
            source['w3mpro'] = 99.0
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
       
    
        
    return np.array(luminosities)

if __name__ == "__main__":
    bals = Table.read(paths.data / "desi_bals_radio.fits")
    qsos = Table.read(paths.data / "desi_qso_radio.fits")

    six_micron_bals = six_micron(bals)
