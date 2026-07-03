
We present the MPS-ATLAS library of stellar model atmospheric structures, 
center-to-limb intensity spectra and disk-integrated spectra.

We provide a jupyter notebook which extract the required data from opur library for the given 
stellar parameters (effective temperature, surface gravity and metallicity) and
set of computations (set1 or set2 which stands for different abundances and mixing length parameter treatment).
The user can also choose the type of the output data ("model", "clv_spectra", "disk_integrated_spectra").

So, the set1 is computed with Grevesse and Sauval 1998 chemical abundances and a constant mixing length 
parameter of 1.25, and for the set2 computation we used chemical abundances from Asplund et al 2009 and 
different mixing length parameters  depending on stellar effective temperature, metallicity, and surface 
gravity of stars from Viani et al 2018. 

The jupyter notebook works with the unziped data, so set1.zip and set2.zip should be unzipped first. 
Note that the output data is selected from the grid for the closest stellar parameters.

ENVIRONMENT:
============
To use the script numpy library is necessary. 

INPUT:
======
To get data, the user should first select the input parameters from the range:
Teff = [3500 to 9000] K 
logg = [3.0 to 5.0]
M/H = [-5.0 to 1.5] 
set_type = ["set1", "set2"]
output_data = ["model", "clv_spectra", "disk_integrated_spectra"]

OUTPUT:
======
if output_data =  "model": 
   the output data is a 2D array of the shape (number_of_depth_points, 7), 
   where "7" describes different quantities:
   model[:, 0] stands for column mass in [g/cm^2]
   model[:, 1] stands for gas temperature in [K]
   model[:, 2] stands for gas pressure in [erg/cm^3]
   model[:, 3] stands for electron number density in [cm^-1]
   model[:, 4] stands for mean Rosseland opacity in [cm^2/g]
   model[:, 5] stands for radiation pressure in [erg/cm^3]
   model[:, 6] stands for turbulent velocity in [cm/s]

if output_data =  "clv_spectra": 
   for this case we provide 1D array of wavelengths, 1D array of position(mu), and 
   2D array of the specific intensity of the shape (number_of_wavelengths, number_of_positions)
   I_\nu in [erg * s^{-1} * cm^{-2} * Hz^{-1} * ster^{-1}]
   
if output_data =  "disk_integrated_spectra": 
   for this case we provide 1D array of wavelengths and 1D array of the disk_integrated flux at one AU from the Sun,
   F_\nu in [erg * s^{-1} * cm^{-2} * Hz^{-1}]

CITATION:
=========
Once the user makes a use of the providing data, we would ask to cite the following papers:

1. Witzke et al. 2021, MPS-ATLAS: A fast all-in-one code for synthesising stellar spectra, 2021A&A...653A..65W

@ARTICLE{2021A&A...653A..65W,
       author = {{Witzke}, V. and {Shapiro}, A.~I. and {Cernetic}, M. and {Tagirov}, R.~V. and {Kostogryz}, N.~M. and {Anusha}, L.~S. and {Unruh}, Y.~C. and {Solanki}, S.~K. and {Kurucz}, R.~L.},
        title = "{MPS-ATLAS: A fast all-in-one code for synthesising stellar spectra}",
      journal = {\aap},
     keywords = {stars: atmospheres, stars: late-type, radiative transfer, opacity, convection, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2021,
        month = sep,
       volume = {653},
          eid = {A65},
        pages = {A65},
          doi = {10.1051/0004-6361/202140275},
archivePrefix = {arXiv},
       eprint = {2105.13611},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021A&A...653A..65W},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

2. Kostogryz et al. 2022, Stellar limb darkening. A new MPS-ATLAS library for Kepler, TESS, CHEOPS, and PLATO passbands, 2022A&A...666A..60K
@ARTICLE{2022A&A...666A..60K,
       author = {{Kostogryz}, N.~M. and {Witzke}, V. and {Shapiro}, A.~I. and {Solanki}, S.~K. and {Maxted}, P.~F.~L. and {Kurucz}, R.~L. and {Gizon}, L.},
        title = "{Stellar limb darkening. A new MPS-ATLAS library for Kepler, TESS, CHEOPS, and PLATO passbands}",
      journal = {\aap},
     keywords = {radiative transfer, methods: numerical, Sun: atmosphere, stars: atmospheres, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Earth and Planetary Astrophysics},
         year = 2022,
        month = oct,
       volume = {666},
          eid = {A60},
        pages = {A60},
          doi = {10.1051/0004-6361/202243722},
archivePrefix = {arXiv},
       eprint = {2206.06641},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022A&A...666A..60K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}






