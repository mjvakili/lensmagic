import numpy as np
import itertools
import h5py
import kmeans_radec

def lens_rand(zmin, zmax):
 
    lens_file = h5py.File('source_zb_'+str(zmax+0.1)+'_0.9.h5')
    print lens_file 

    z = lens_file["zb"][:]
    ra = lens_file["ra"][:]
    dec = lens_file["dec"][:]
    e1 = lens_file["e1"][:]
    e2 = lens_file["e2"][:]
    w = lens_file["w"][:]
    size = lens_file["rad"][:]
    snr = lens_file["snr"][:]
    
    mask_snr = snr > np.median(snr)    
    mask_size = size > np.median(size)    


    np.savetxt("esd/source_faint_zbin_"+str(zmin)+"_"+str(zmax) , \
    zip(ra[~mask_snr] , dec[~mask_snr] , e1[~mask_snr] , -1.*e2[~mask_snr], w[~mask_snr] , z[~mask_snr]), \
    fmt='%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f') 
    
    np.savetxt("esd/source_bright_zbin_"+str(zmin)+"_"+str(zmax) , \
    zip(ra[mask_snr] , dec[mask_snr] , e1[mask_snr] , -1.*e2[mask_snr], w[mask_snr] , z[mask_snr]), \
    fmt='%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f') 
    
    np.savetxt("esd/source_small_zbin_"+str(zmin)+"_"+str(zmax) , \
    zip(ra[~mask_size] , dec[~mask_size] , e1[~mask_size] , -1.*e2[~mask_size], w[~mask_size] , z[~mask_size]), \
    fmt='%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f') 
    
    np.savetxt("esd/source_large_zbin_"+str(zmin)+"_"+str(zmax) , \
    zip(ra[mask_size] , dec[mask_size] , e1[mask_size] , -1.*e2[mask_size], w[mask_size] , z[mask_size]), \
    fmt='%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f\t%1.8f') 
    
    return None


if __name__ == '__main__':

    zmin, zmax = 0.1, 0.3
    lens_rand(zmin, zmax) 
