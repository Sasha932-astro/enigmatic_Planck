import numpy as np

class ReadData:
    def read_model_atmosphere(self, mh, teff, logg, set_type):
        file_name = f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_model_atmosphere.dat"
        return np.genfromtxt(file_name, skip_header=2, skip_footer=23)

    def read_clv_spectra(self, mh, teff, logg, set_type):
        file_name = f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_intensity_spectra.dat"
        data = np.loadtxt(file_name, skiprows=2)
        return data[:, 0], data[:, 1:]

    def read_disk_integrated_spectra(self, mh, teff, logg, set_type):
        file_name = f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_flux_spectra.dat"
        return np.loadtxt(file_name, skiprows=1, unpack=True)
    
    def read_mu_positions(self, mh, teff, logg, set_type):
        file_name = f"{set_type}/MH{mh}/teff{teff}/logg{logg}/mpsa_intensity_spectra.dat"
        f = open(file_name, "r")
        data = f.readlines()
        muval = data[1].split()[2:]
        return np.array(muval).astype(float)