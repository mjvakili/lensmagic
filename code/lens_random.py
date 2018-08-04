import numpy as np
import itertools
import h5py
import kmeans_radec

def lens_rand(zmin, zmax):
 
    lens_file = h5py.File('LRG_lum.h5')
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_mask = (lens_z>zmin)&(lens_z<zmax)
    
    lens_z = lens_z[lens_mask]
    lens_ra = lens_ra[lens_mask]     
    lens_dec = lens_dec[lens_mask]   
    
    nlens = len(lens_z)

    lens_input = np.vstack([np.arange(nlens), lens_ra , lens_dec, 
                            lens_z, np.ones((nlens))]).T

    
    #nrands = 80 * nlens
    
    random_file = h5py.File("randoms.h5")
    random_coord = random_file["random"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]
    random_file.close()
       
    nrands = len(random_ra)
    
    shuffle = np.random.choice(nrands, 10*nlens)
    random_ra = random_ra[shuffle]
    random_dec = random_dec[shuffle]
    nrands = len(random_ra)

    print "nrands/nlens" , nrands*1./nlens
    
    zdist = lens_z
    hist, bins = np.histogram(zdist, bins=100)
    bin_midpoints = bins[:-1] + np.diff(bins)/2
    cdf = np.cumsum(hist)
    cdf = cdf / float(cdf[-1])
    values = np.random.rand(nrands)
    value_bins = np.searchsorted(cdf, values)
    random_z = bin_midpoints[value_bins]
                                                   
    rands_input = np.stack([np.arange(nrands) , random_ra , random_dec , random_z, np.ones((nrands))]).T
    
    return lens_input , rands_input

def lens_rand_jk(zmin, zmax, jk):
 
    lens_file = h5py.File('LRG_lum_jk.h5')
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    coord_lens_jk = np.vstack([lens_ra, lens_dec]).T
    centers = np.loadtxt("centers.txt")
    labels_lens = kmeans_radec.find_nearest(coord_lens_jk, centers)

    label_mask = labels_lens == jk
    lens_mask = (lens_z>zmin)&(lens_z<zmax)&(~label_mask)
    
    lens_z = lens_z[lens_mask]
    lens_ra = lens_ra[lens_mask]     
    lens_dec = lens_dec[lens_mask]   
    
    nlens = len(lens_z)


    lens_input = np.vstack([np.arange(nlens), lens_ra , lens_dec, 
                            lens_z, np.ones((nlens))]).T
    
    
    random_file = h5py.File("randoms_jk.h5")
    random_coord = random_file["random_jk"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]



    centers = np.loadtxt("centers.txt")
    labels_random = kmeans_radec.find_nearest(random_coord, centers)

    label_mask = labels_random == jk
    random_ra = random_ra[~label_mask]
    random_dec = random_dec[~label_mask]
    random_file.close()
       
    nrands = len(random_ra)
    shuffle = np.random.choice(nrands, 10*nlens)
    random_ra = random_ra[shuffle]
    random_dec = random_dec[shuffle]
    nrands = len(random_ra)
    
    print "nrands/nlens" , nrands*1./nlens
    
    zdist = lens_z
    hist, bins = np.histogram(zdist, bins=100)
    bin_midpoints = bins[:-1] + np.diff(bins)/2
    cdf = np.cumsum(hist)
    cdf = cdf / float(cdf[-1])
    values = np.random.rand(nrands)
    value_bins = np.searchsorted(cdf, values)
    random_z = bin_midpoints[value_bins]
                                                   
    rands_input = np.stack([np.arange(nrands) , random_ra , random_dec , random_z, np.ones((nrands))]).T
    
    return lens_input , rands_input

if __name__ == '__main__':

    zmin, zmax = 0.1, 0.3
        
    lens_output, rand_output = lens_rand(zmin, zmax)
    np.savetxt("esd/lens_file_"+str(zmin)+"_"+str(zmax), lens_output , fmt='%i\t%1.8f\t%1.8f\t%1.8f\t%i')
    np.savetxt("esd/rand_file_"+str(zmin)+"_"+str(zmax), rand_output , fmt='%i\t%1.8f\t%1.8f\t%1.8f\t%i')

    for jk in range(100):

        lens_output, rand_output = lens_rand_jk(zmin, zmax, jk)
        np.savetxt("esd/lens_file_"+str(zmin)+"_"+str(zmax)+"_jk_"+str(jk), lens_output , fmt='%i\t%1.8f\t%1.8f\t%1.8f\t%i')
        np.savetxt("esd/rand_file_"+str(zmin)+"_"+str(zmax)+"_jk_"+str(jk), rand_output , fmt='%i\t%1.8f\t%1.8f\t%1.8f\t%i')
