import ase
import ase.atom
import numpy as np
from ase import io
from ase.geometry.cell import cellpar_to_cell

class dipole_cell:
    def __init__(
            self,
            trajectory: list,
            cell: np.ndarray,
            charge_dict: dict, 
            ):
        self.trajectory = trajectory
        if np.shape(cell) == (6, ):
            cell = cellpar_to_cell(cell)
            self.cell = cell
        elif np.shape(cell) == (3, 3):
            self.cell = cell
        else:
            raise PermissionError("Unsupporting cell format. Only [x, y, z, alpha, beta, gamma] or 3*3 matrix is permitted.")
        self.charge_dict = charge_dict

    def get_dipole_x(self, wrap=False):
        trajectory = self.trajectory
        charge_dict = self.charge_dict
        tot_dipole = []
        if type(trajectory) == ase.atoms.Atoms:
            trajectory = [trajectory]   
        for frame in trajectory:
            if wrap == True:
                frame.set_cell(self.cell)
                frame.set_pbc(True)
                frame.wrap()
            elif wrap == False:
                pass 
            atoms = frame
            coord_x = atoms.get_positions()[:, 0]
            charge = np.array([charge_dict[s] for s in atoms.symbols])
            cell_dipole_x = coord_x*charge
            tot_dipole.append([int(trajectory.index(frame)), np.float64(np.sum(cell_dipole_x))])
            print(f'Frame {trajectory.index(frame)} Dipole: done.')
        return np.array(tot_dipole)
    
    def get_dipole_y(self, wrap=False):
        trajectory = self.trajectory
        charge_dict = self.charge_dict
        tot_dipole = []
        if type(trajectory) == ase.atoms.Atoms:
            trajectory = [trajectory]   
        for frame in trajectory:
            if wrap == True:
                frame.set_cell(self.cell)
                frame.set_pbc(True)
                frame.wrap()
            elif wrap == False:
                pass 
            atoms = frame
            coord_y = atoms.get_positions()[:, 1]
            charge = np.array([charge_dict[s] for s in atoms.symbols])
            cell_dipole_y = coord_y*charge
            tot_dipole.append([int(trajectory.index(frame)), np.float64(np.sum(cell_dipole_y))])
            print(f'Frame {trajectory.index(frame)} Dipole: done.')
        return np.array(tot_dipole)
    
    def get_dipole_z(self, wrap=False):
        trajectory = self.trajectory
        charge_dict = self.charge_dict
        tot_dipole = []
        if type(trajectory) == ase.atoms.Atoms:
            trajectory = [trajectory]   
        for frame in trajectory:
            if wrap == True:
                frame.set_cell(self.cell)
                frame.set_pbc(True)
                frame.wrap()
            elif wrap == False:
                pass 
            atoms = frame
            coord_z = atoms.get_positions()[:, 2]
            charge = np.array([charge_dict[s] for s in atoms.symbols])
            cell_dipole_z = coord_z*charge
            tot_dipole.append([int(trajectory.index(frame)), np.float64(np.sum(cell_dipole_z))])
            print(f'Frame {trajectory.index(frame)} Dipole: done.')
        return np.array(tot_dipole)
