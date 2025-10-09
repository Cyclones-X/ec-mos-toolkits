# Function 
import ase
import ase.atom
import numpy as np
from ase import io
from scipy import fft
from scipy import constants
from ase.geometry.cell import cellpar_to_cell

# AU
AU_TO_ANG = 5.29177208590000e-01
AU_TO_EV = 2.72113838565563e01
AU_TO_EV_EVERY_ANG = AU_TO_EV / AU_TO_ANG
DEBYE_TO_EA = 1.0 / 4.8

# epsilon [e/(V Angstrom)]
EPSILON = constants.epsilon_0 / constants.elementary_charge * constants.angstrom

class ElecPotential_fourier_solvation:
    def __init__(
            self,
            trajectory: list,
            cell: np.ndarray,
            z_coord: np.ndarray,
            spread_dict: dict, 
            charge_dict: dict, 
            n_k_vectors: int,
            ):
        self.trajectory = trajectory
        if np.shape(cell) == (6, ):
            cell = cellpar_to_cell(cell)
            self.cell = cell
        elif np.shape(cell) == (3, 3):
            self.cell = cell
        else:
            raise PermissionError("Unsupporting cell format. Only [x, y, z, alpha, beta, gamma] or 3*3 matrix is permitted.")
        self.z_coord = z_coord
        self.spread_dict = spread_dict
        self.charge_dict = charge_dict
        self.n_k_vectors = n_k_vectors
        self.XY_Plane_Area = np.linalg.norm(np.cross(cell[0], cell[1]))
        self.l_box = np.linalg.norm(cell[2])

    def Total_density(self):
        '''
        . Z-direction
        . Unit e/A^3
        '''
        l_box = self.l_box
        XY_Plane_Area = self.XY_Plane_Area
        trajectory = self.trajectory
        spread_dict = self.spread_dict
        charge_dict = self.charge_dict
        z_coord = self.z_coord
        rho = np.zeros(len(self.z_coord), dtype=np.float64) # Real rho
        if type(trajectory) == ase.atoms.Atoms:
            trajectory = [trajectory]            
        for frame in trajectory:
            atoms = frame
            mu_z = atoms.get_positions()[:, 2]
            spread = np.array([spread_dict[s] for s in atoms.symbols])
            charge = np.array([charge_dict[s] for s in atoms.symbols])
            # Calculate charge density with Fourier Series Expansion
            rho_frame = np.zeros(len(self.z_coord), dtype=np.complex128)
            for nk in np.arange(-self.n_k_vectors, self.n_k_vectors+1): # With k=0, considering the non-neutral conditions
                k_n = 2*np.pi/l_box*nk
                ft_trans_density = 1/XY_Plane_Area*charge*np.exp(-1j*k_n*mu_z)*np.exp(-k_n**2/4*(2*spread**2))
                ft_trans_density_k_n = np.sum(ft_trans_density)
                rho_kn = 1/l_box*np.exp(1j*k_n*z_coord)*ft_trans_density_k_n
                rho_frame += rho_kn
            rho_frame = np.real(rho_frame)
            rho += rho_frame
            print(f'Frame {trajectory.index(frame)} Density: done.')
        rho = rho/len(trajectory) # Average of trajectory. 
        return (z_coord, rho)
    
    def Electrostatic_potential(self):  # Solve Poisson Equation  -∇²φ = ρ/ε₀, and multiply an electronic charge (-e)
        '''
        . Z-direction
        . Hartree Potential Energy, which is the electrons-(electrons+ions) interactions, Unit eV 
        '''
        l_box = self.l_box
        XY_Plane_Area = self.XY_Plane_Area
        trajectory = self.trajectory
        spread_dict = self.spread_dict
        charge_dict = self.charge_dict
        z_coord = self.z_coord
        phi = np.zeros(len(self.z_coord), dtype=np.float64) # Real phi
        if type(trajectory) == ase.atoms.Atoms:
            trajectory = [trajectory]
        for frame in trajectory:
            atoms = frame
            mu_z = atoms.get_positions()[:, 2]
            spread = np.array([spread_dict[s] for s in atoms.symbols])
            charge = np.array([charge_dict[s] for s in atoms.symbols])
            # Calculate Electrostatic Potential with Fourier Series Expansion
            phi_frame = np.zeros(len(self.z_coord), dtype=np.complex128)
            # k should not equal to 0
            for nk in np.concatenate((np.arange(-self.n_k_vectors, 0), np.arange(1, self.n_k_vectors+1))):
                k_n = 2*np.pi/l_box*nk
                ft_trans_density = 1/XY_Plane_Area*charge*np.exp(-1j*k_n*mu_z)*np.exp(-k_n**2/4*(2*spread**2))
                ft_trans_density_k_n = 1/(EPSILON*k_n**2)*np.sum(ft_trans_density)
                phi_kn = 1/l_box*np.exp(1j*k_n*z_coord)*ft_trans_density_k_n
                phi_frame += phi_kn
            phi_frame = np.real(phi_frame)
            phi += phi_frame
            print(f'Frame {trajectory.index(frame)} Potential: done.')
        phi = phi/len(trajectory) # Average of trajectory, total electrostatic potential φ
        phi = -phi # Electrostatic Potential (φ) -> Hartree Potential (E), unit [V -> eV], etc.
        return (z_coord, phi)

class macro_avg:  # A macroscopic average method for microscopic electrostatic potential data. This method is writen by taking a part from CP2KDATA. Refernce: https://robinzyb.github.io/cp2kdata/_modules/cp2kdata/cube/cube.html#Cp2kCube.get_mav 
    def __init__(
            self,
            input_coord: np.ndarray,
            input_pot: np.ndarray,
            l_box: float,
            ):
        self.zcoord = input_coord
        self.potential = input_pot
        self.l_box = l_box
    
    def square_wave_filter(self, x, l: float, cell_z: float):
        half_l = l/2
        x_1st, x_2nd = np.array_split(x, 2)
        y_1st = np.heaviside(half_l - np.abs(x_1st), 0)/l
        y_2nd = np.heaviside(half_l - np.abs(x_2nd-cell_z), 0)/l
        y = np.concatenate([y_1st, y_2nd])
        return y
    
    def get_mav(self, l1, l2=0, ncov=1):
        length = self.l_box
        z_coord = self.zcoord
        theta_1_fft = fft.fft(self.square_wave_filter(z_coord, l1, length))
        pav_fft = fft.fft(self.potential)
        mav_fft = pav_fft*theta_1_fft*length/len(z_coord)
        if ncov == 2:
            theta_2_fft = fft.fft(self.square_wave_filter(z_coord, l2, length))
            mav_fft = mav_fft*theta_2_fft*length/len(z_coord)
        mav = fft.ifft(mav_fft)
        return z_coord, np.real(mav)
