import numpy as np
import h5py

def lum_exclusion():
    
    lens_file = h5py.File("LRG_lum.h5")
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_rad = lens_file["radius"][:]*0.21 
    for i in range(3):

        source_file = h5py.File("source_zb_"+str(0.4+i*(0.2))+"_0.9.h5")  

        source_ra = source_file["ra"][:]
        source_dec = source_file["dec"][:]
	
        x =  source_ra
        y =  source_dec     
	rad_mask = np.where((x<np.any(lens_ra+lens_rad/3600.))&(x>np.all(lens_ra-lens_rad/3600.)))#&(y>np.all(lens_dec-lens_rad/3600.))&(y<np.all(lens_dec+lens_rad/3600.)))
  
        print rad_mask[0].shape

        return None
  
if __name__ == '__main__':

   lum_exclusion()
