import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
from ase.geometry.cell import cell_to_cellpar

from dataclasses import dataclass

@dataclass
class Surface_sp:
    get_H: np.ndarray
    get_O_wat: np.ndarray
    get_O_lattice: np.ndarray
    get_O_bridge: np.ndarray
    get_M_unsatur: np.ndarray

class surface_identify:
    def __init__(
            self, 
            atoms,
            cell: np.ndarray,
            metal_name: str,
            cutoff_O_M: float,
            cutoff_O_H: float,
            ):
        if type(atoms) == mda.core.universe.Universe:
            self.atoms = atoms.atoms
        elif type(atoms) == str:
            self.atoms = mda.Universe(atoms).atoms
        else:
            raise PermissionError("Unsupporting coordinate format. Only file path '*.xyz' or MDAnalysis.Universe is permitted.")
        
        if np.shape(cell) == (6, ):
            self.cell = cell
        elif np.shape(cell) == (3, 3):
            self.cell = cell_to_cellpar(cell)
        else:
            raise PermissionError("Unsupporting cell format. Only [x, y, z, alpha, beta, gamma] or 3*3 matrix is permitted.")
        self.M = metal_name
        self.cutoff_O_M = cutoff_O_M
        self.cutoff_O_H = cutoff_O_H
    def idx_elements(self):
        atoms = self.atoms
        M = self.M
        O_idx = atoms.select_atoms('name O').indices
        H_idx = atoms.select_atoms('name H').indices
        M_idx = atoms.select_atoms(f'name {M}').indices
        return O_idx, H_idx, M_idx
    
    def coordination_count(self):
        atoms = self.atoms
        atoms.dimensions = self.cell
        O_idx, H_idx, M_idx = self.idx_elements()
        pairs_OH, _ = capped_distance(reference=atoms[O_idx], configuration=atoms[H_idx], max_cutoff=self.cutoff_O_H, box=self.cell)
        pairs_OM, _ = capped_distance(reference=atoms[O_idx], configuration=atoms[M_idx], max_cutoff=self.cutoff_O_M, box=self.cell)

        cn_OH = np.bincount(pairs_OH[:, 0], minlength=O_idx.shape[0])
        cn_OM = np.bincount(pairs_OM[:, 0], minlength=O_idx.shape[0])
        return cn_OH, cn_OM
    
    def idx_oxy_wat_lattice(self):
        O_idx = self.idx_elements()[0]
        cn_OH, cn_OM = self.coordination_count()
        O_wat_idx = O_idx[(cn_OH >= 0) & (cn_OM <= 1)]
        O_lat_idx = O_idx[~np.isin(O_idx, O_wat_idx)]
        return O_wat_idx, O_lat_idx
    
    def surface_species(self, lattice_cn_MO) -> Surface_sp:
        O_idx, H_idx, M_idx = self.idx_elements()
        atoms = self.atoms
        cn_OH, cn_OM = self.coordination_count()
        O_wat_idx, O_lat_idx = self.idx_oxy_wat_lattice()

        O_bridge = np.array([O_idx[cn_OM == 2], atoms[O_idx[cn_OM == 2]].positions[:, 2]]).T
        O_bridge_z = O_bridge[:, 1]
        O_bridge_idx = np.array([O_bridge[O_bridge_z < np.mean(O_bridge_z)][:, 0].astype(int), O_bridge[O_bridge_z > np.mean(O_bridge_z)][:, 0].astype(int)])

        pairs_lat_M, _ = capped_distance(reference=atoms[M_idx], configuration=atoms[O_lat_idx], max_cutoff=self.cutoff_O_M, box=self.cell)
        cn_Olat_M = np.unique(pairs_lat_M[:, 0], return_counts=True)
        M_unsatur = M_idx[cn_Olat_M[0][cn_Olat_M[1] != lattice_cn_MO]]
        M_unsatur_z = atoms[M_unsatur].positions[:, 2]
        M_unsatur_idx = np.array([M_unsatur[M_unsatur_z < np.mean(M_unsatur_z)].astype(int), M_unsatur[M_unsatur_z > np.mean(M_unsatur_z)].astype(int)])
        
        return Surface_sp(
            get_H=self.idx_elements()[1],
            get_O_wat=O_wat_idx,
            get_O_lattice=O_lat_idx,
            get_O_bridge=O_bridge_idx,
            get_M_unsatur=M_unsatur_idx
        )
        

class surface_proton:
    def __init__(self, 
                 cell: np.ndarray,
                 metal_name: str,
                 atoms: mda.core.universe.Universe,
                 cutoff_O_M: float,
                 cutoff_O_H: float,
                 H_idx=None,
                 O_wat_idx=None,
                 O_bridge_idx=None,
                 M_unsatur_idx=None,     
                 ):
        self.atoms = atoms.atoms
        if np.shape(cell) == (6, ):
            self.cell = cell
        elif np.shape(cell) == (3, 3):
            self.cell = cell_to_cellpar(cell)
        else:
            raise PermissionError("Unsupporting cell format. Only [x, y, z, alpha, beta, gamma] or 3*3 matrix is permitted.")
        self.M = metal_name
        self.H_idx = H_idx
        self.O_wat_idx = O_wat_idx
        self.O_bridge_idx = O_bridge_idx
        self.M_unsatur_idx = M_unsatur_idx
        self.cutoff_O_M = cutoff_O_M
        self.cutoff_O_H = cutoff_O_H
    
    def bridge_protonic_charge(self): # Positive Charge
        pairs_Obr_H, _ = capped_distance(reference=self.atoms[self.O_bridge_idx], configuration=self.atoms[self.H_idx], max_cutoff=self.cutoff_O_H, box=self.cell)
        protonic_charge = pairs_Obr_H.shape[0]
        return protonic_charge
    
    def get_ads_water_O_idx(self):
        pairs_Owat_M, _ = capped_distance(reference=self.atoms[self.O_wat_idx], configuration=self.atoms[self.M_unsatur_idx], max_cutoff=self.cutoff_O_M, box=self.cell)
        cn_Owat_M = np.bincount(pairs_Owat_M[:, 0], minlength=self.O_wat_idx.shape[0])
        ads_idx = self.O_wat_idx[cn_Owat_M == 1]
        return ads_idx
    
    def water_hydroxide_charge(self): # Negative Charge
        pairs_Owat_H, _ = capped_distance(reference=self.atoms[self.get_ads_water_O_idx()], configuration=self.atoms[self.H_idx], max_cutoff=self.cutoff_O_H, box=self.cell)
        proton_neutral = 2 * self.get_ads_water_O_idx().shape[0] # 2 hydrogens per water molecule
        cn_Owat_H = np.bincount(pairs_Owat_H[:, 0], minlength=self.get_ads_water_O_idx().shape[0])
        OH_1 = np.count_nonzero(cn_Owat_H == 1)
        O_2 = np.count_nonzero(cn_Owat_H == 0)
        proton_non_neutral = pairs_Owat_H.shape[0]
        hydroxide_charge = proton_non_neutral - proton_neutral
        return hydroxide_charge, (OH_1, O_2)

        
