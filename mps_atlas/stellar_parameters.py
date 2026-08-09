
from typing import Tuple
import numpy as np 


class GridOfStellarParameters:
    def __init__(
        self,
        mh_input: float, 
        teff_input: float, 
        logg_input: float
    ) -> None:
        self.mh_input = mh_input
        self.teff_input = teff_input
        self.logg_input = logg_input
        
    def get_closest_mh_teff_logg(self) -> Tuple[float, int, float]:
        teff_grid = np.loadtxt("grid_teff.txt")
        mh_grid = np.loadtxt("grid_mh.txt")
        logg_grid = np.loadtxt("grid_logg.txt")
        index_teff = np.argmin(abs(teff_grid - self.teff_input))
        index_mh = np.argmin(abs(mh_grid - self.mh_input))
        index_logg = np.argmin(abs(logg_grid - self.logg_input))
        return float(mh_grid[index_mh]), int(teff_grid[index_teff]), float(logg_grid[index_logg])
        
    