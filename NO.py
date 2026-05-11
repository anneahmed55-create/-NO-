import numpy as np
from scipy.constants import h, c
import matplotlib.pyplot as plt



def param_Atmosphere(Atmosphere):
    
    data = np.loadtxt(Atmosphere, delimiter="," , skiprows=1)
    
    z = data[:, 0]
    T = data[:, 1]
    O = data[:, 2]
    NO = data[:, 3]
    
    return z, T, O, NO


class NO_flux:
   
    
    def __init__(self, z, T, O, NO):
       
        self.z = z
        self.z = z * 1e5 
        self.T = T
        self.O = O
        self.NO = NO
        
        
        self.h = h * 1e7  
        self.c = c * 100    
        
        
        self.A10 = 12.5    
        self.k_NO_O = 2.8e-11   
        self.S_E = 1.06e-4       
        
        self.l = 5.3e-4     # cm
        self.nu = self.c / self.l # Hz
        self.photon_energy = self.h * self.nu  # erg
        
    def compute_NO(self):
       
        collisional = self.k_NO_O * np.exp(-2700.0 / self.T) * self.O
        
        numerator = self.S_E + collisional
        denominator = self.A10 + self.k_NO_O * self.O
        
        
        denominator = np.maximum(denominator, 1e-20)
        
        
        NO = (numerator / denominator) * self.NO # NO*
        
        return NO
    
    def compute_L_NO(self):
        
        NO = self.compute_NO()
        L_NO = self.photon_energy * self.A10 * NO
        return L_NO, NO
    
    
    
    def compute_flux(self):
       
        L_NO = self.compute_L_NO()
        
        
        flux = np.trapz(L_NO, self.z)
        
        return flux, L_NO
    
    def get_flux(self, z_min, z_max):
    
        mask = (self.z >= z_min) & (self.z <= z_max)
        z_m = self.z[mask]
        L_NO_m = self.compute_L_NO()
        L_NO_m = L_NO_m[mask]
        
        flux = np.trapz(L_NO_m, z_m)
        return flux
    
    
    

if __name__ == "__main__":
    
    file= r"C:\Users\Amal\OneDrive\Desktop\mynumpy\Atmosphere.dat"
    
    z, T, O, NO = param_Atmosphere(file)
    
    
    calc = NO_flux(z, T, O, NO)
    total_flux, L_NO = calc.compute_flux()
    NO = calc.compute_NO()
    
    print("flux", total_flux)
    
    
   
   
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    
    axes[0].plot(NO, z, 'b-', label='NO')
    axes[0].plot(O, z, 'g-', label='O')
    axes[0].set_xlabel('Концентрация(cm$^{-3}$)')
    axes[0].set_ylabel('Z')
    axes[0].set_xscale('log')
    axes[0].legend()
    axes[0].set_title('Концентрация')
    axes[0].grid(True, alpha=0.3)

    
    axes[1].plot(NO, z, 'r-')
    axes[1].set_xlabel('[NO*] (cm$^{-3}$)')
    axes[1].set_ylabel('Z')
    axes[1].set_xscale('log')
    axes[1].set_title('NO*')
    axes[1].grid(True, alpha=0.3)
    
    
    
    axes[2].plot(L_NO[0], z, 'm-')
    axes[2].set_xlabel('L$_{NO}$ (erg cm$^{-3}$ s$^{-1}$)')
    axes[2].set_ylabel('Z')
    axes[2].set_xscale('log')
    axes[2].set_title('L_No')
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()