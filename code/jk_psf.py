import h5py
import numpy as np
import kmeans_radec
import matplotlib.pyplot as plt
import treecorr
import seaborn as sns
from kmeans_radec import KMeans, kmeans_sample

random_file = h5py.File("randoms_jk.h5")
random_coord = random_file["random_jk"][:] 
random_ra = random_coord[:,0]
random_dec = random_coord[:,1]
random_file.close()

random_file_jk = h5py.File("randoms_jk.h5")
coord_rand_jk = random_file_jk["random_jk"][:] 
random_ra_jk = coord_rand_jk[:,0]
random_dec_jk = coord_rand_jk[:,1]
random_file_jk.close()
centers = np.loadtxt("centers.txt")
labels_random_jk = kmeans_radec.find_nearest(coord_rand_jk, centers)

lens_file_jk = h5py.File("LRG_dense_jk.h5")
dense_jk_mask = h5py.File("dense_jk_mask.h5")
dense_mask = dense_jk_mask["mask"][:]
dense_jk_mask.close()
lens_z_jk = lens_file_jk["redshift"][dense_mask]
lens_ra_jk = lens_file_jk["RA"][dense_mask]
lens_dec_jk = lens_file_jk["DEC"][dense_mask]
lens_file_jk.close()    
coord_lens_jk = np.vstack([lens_ra_jk, lens_dec_jk]).T 
centers = np.loadtxt("centers.txt")
labels_lens_jk = kmeans_radec.find_nearest(coord_lens_jk, centers)

random_ra = random_coord[:,0]
random_dec = random_coord[:,1]
random_file.close()
#shuffled = np.random.choice(len(random_ra) , 2000000)
#random_coord = random_coord[shuffled]
random_ra = random_coord[:,0]
random_dec = random_coord[:,1]

source_file = h5py.File("source_zb_0.4_0.9.h5")
source_jk_mask = h5py.File("source_mask_zb_0.4_0.9.h5")
source_mask = source_jk_mask["mask"][:]
source_jk_mask.close()
source_ra = source_file["ra"][source_mask]
source_dec = source_file["dec"][source_mask]
source_z = source_file["zb"][source_mask]
source_e1 = source_file["e1"][source_mask]
source_e2 = source_file["e2"][source_mask]
source_w = source_file["w"][source_mask]
source_size = source_file["snr"][source_mask]
source_file.close()

coord_source_jk = np.vstack([source_ra, source_dec]).T 
centers = np.loadtxt("centers.txt")
labels_source_jk = kmeans_radec.find_nearest(coord_source_jk, centers)

def xshear_lens_jk(zl1, zl2, jk_label):

    source_jk_mask = labels_source_jk == jk_label
    source_mask = (source_z>zl2+0.1)|(source_z==zl2+0.1)
    source_mask = (~source_jk_mask)&(source_mask)

    source_cat=treecorr.Catalog(x=source_ra[source_mask],y=source_dec[source_mask],g1=source_e1[source_mask],g2=-1.*source_e2[source_mask],w=source_w[source_mask],x_units='degree',y_units='degree')        
    
    mask_z_lens = (lens_z_jk>zl1)&(lens_z_jk<zl2)
    lens_mask = labels_lens_jk == jk_label
    lens_mask = (~lens_mask)&(mask_z_lens)
    
    lens_cat = treecorr.Catalog(x=lens_ra_jk[lens_mask], y=lens_dec_jk[lens_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 20, min_sep=0.5, max_sep=250, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat, num_threads = 20)
    r_lens_h, xt_lens_h, xx_lens_h , w_lens_h = ng.meanr , ng.xi , ng.xi_im, ng.npairs
    w_lens_h = w_lens_h / np.sum(w_lens_h)
   
    #print "lens", xt_lens_h

    random_mask = labels_random_jk == jk_label
    random_cat = treecorr.Catalog(x=random_ra_jk[~random_mask], y=random_dec_jk[~random_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 20, min_sep=0.5, max_sep=250, sep_units='arcmin', verbose=1)
    ng.process(random_cat, source_cat, num_threads = 20)
    r_rand_h, xt_rand_h, xx_rand_h , w_rand_h = ng.meanr , ng.xi, ng.xi_im, ng.npairs
    w_rand_h = w_rand_h / np.sum(w_rand_h)
  
    boost = 1. * len(random_ra_jk[~random_mask]) * w_lens_h / (len(lens_ra_jk[lens_mask]) * w_rand_h)

    #print "random", xt_rand_h 

    return r_lens_h, xt_lens_h, xx_lens_h , w_lens_h, r_rand_h, xt_rand_h, xx_rand_h , w_rand_h, boost 



def lum_all(nbins):

    xshear_file = h5py.File("psf_jk.h5" , "w")



    lum_jk_L1 = np.zeros((100,9,nbins))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, xx_lens_h , w_lens_h, r_rand_h, xt_rand_h, xx_rand_h , w_rand_h , boost = xshear_lens_jk(0.1, 0.3, i)
        
        lum_jk_L1[i,0,:] = r_lens_h
        lum_jk_L1[i,1,:] = xt_lens_h
        lum_jk_L1[i,2,:] = xx_lens_h
        lum_jk_L1[i,3,:] = w_lens_h
        lum_jk_L1[i,4,:] = r_rand_h
        lum_jk_L1[i,5,:] = xt_rand_h
        lum_jk_L1[i,6,:] = xx_rand_h
        lum_jk_L1[i,7,:] = w_rand_h
        lum_jk_L1[i,8,:] = boost
    
        print "done with", i
    xshear_file.create_dataset("lum_L1", (100,9,nbins) , data = lum_jk_L1)

    lum_jk_L2 = np.zeros((100,9,nbins))

    for i in xrange(100):
    
        
        r_lens_h, xt_lens_h, xx_lens_h , w_lens_h, r_rand_h, xt_rand_h, xx_rand_h , w_rand_h, boost = xshear_lens_jk(0.3, 0.5, i)
        
        lum_jk_L2[i,0,:] = r_lens_h
        lum_jk_L2[i,1,:] = xt_lens_h
        lum_jk_L2[i,2,:] = xx_lens_h
        lum_jk_L2[i,3,:] = w_lens_h
        lum_jk_L2[i,4,:] = r_rand_h
        lum_jk_L2[i,5,:] = xt_rand_h
        lum_jk_L2[i,6,:] = xx_rand_h
        lum_jk_L2[i,7,:] = w_rand_h
        lum_jk_L2[i,8,:] = boost
    
        print "done with", i
    xshear_file.create_dataset("lum_L2", (100,9,nbins) , data = lum_jk_L2)
    lum_jk_L3 = np.zeros((100,9,nbins))

    for i in xrange(100):
   
        r_lens_h, xt_lens_h, xx_lens_h , w_lens_h, r_rand_h, xt_rand_h, xx_rand_h , w_rand_h , boost= xshear_lens_jk(0.5, 0.7, i)
        
        lum_jk_L3[i,0,:] = r_lens_h
        lum_jk_L3[i,1,:] = xt_lens_h
        lum_jk_L3[i,2,:] = xx_lens_h
        lum_jk_L3[i,3,:] = w_lens_h
        lum_jk_L3[i,4,:] = r_rand_h
        lum_jk_L3[i,5,:] = xt_rand_h
        lum_jk_L3[i,6,:] = xx_rand_h
        lum_jk_L3[i,7,:] = w_rand_h
        lum_jk_L3[i,8,:] = boost
        print "done with", i

    xshear_file.create_dataset("lum_L3", (100,9,nbins) , data = lum_jk_L3)
 
    return None

if __name__ == '__main__':

   lum_all(20)
