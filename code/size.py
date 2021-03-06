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
lens_z_jk = lens_file_jk["redshift"][:]
lens_ra_jk = lens_file_jk["RA"][:]
lens_dec_jk = lens_file_jk["DEC"][:]
lens_file_jk.close()    
coord_lens_jk = np.vstack([lens_ra_jk, lens_dec_jk]).T 
centers = np.loadtxt("centers.txt")
labels_lens_jk = kmeans_radec.find_nearest(coord_lens_jk, centers)

random_ra = random_coord[:,0]
random_dec = random_coord[:,1]
random_file.close()
shuffled = np.random.choice(len(random_ra) , 2000000)
random_coord = random_coord[shuffled]
random_ra = random_coord[:,0]
random_dec = random_coord[:,1]

source_file = h5py.File("source_zb_0.4_0.9.h5")
source_ra = source_file["ra"][:]
source_dec = source_file["dec"][:]
source_z = source_file["zb"][:]
source_e1 = source_file["e1"][:]
source_e2 = source_file["e2"][:]
source_w = source_file["w"][:]
source_size = source_file["rad"][:]
source_file.close()

coord_source_jk = np.vstack([source_ra, source_dec]).T 
centers = np.loadtxt("centers.txt")
labels_source_jk = kmeans_radec.find_nearest(coord_source_jk, centers)


def xshear_lens_jk(zl1, zl2, jk_label):

    source_jk_mask = labels_source_jk == jk_label
    source_z_mask = (source_z>zl2+0.1)|(source_z==zl2+0.1)
    source_rad_mask = source_size > np.median(source_size)
    
    #### LARGE GALAXIES ####
    
    source_mask = (~source_jk_mask)&(source_z_mask)&(source_rad_mask)
    
    source_cat=treecorr.Catalog(x=source_ra[size_mask],y=source_dec[size_mask],g1=source_e1[size_mask],g2=-1.*source_e2[size_mask],w=source_w[size_mask],x_units='degree',y_units='degree')        
    
    
    mask_z_lens = (lens_z_jk>zl1)&(lens_z_jk<zl2)
    lens_mask = labels_lens_jk == jk_label
    lens_mask = (~lens_mask)&(mask_z_lens)
    
    lens_cat = treecorr.Catalog(x=lens_ra_jk[lens_mask], y=lens_dec_jk[lens_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    r_lens_h, xt_lens_h = ng.meanr , ng.xi 
    
    random_mask = labels_random_jk == jk_label
    random_cat = treecorr.Catalog(x=random_ra_jk[~random_mask], y=random_dec_jk[~random_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(random_cat, source_cat)
    r_rand_h, xt_rand_h = ng.meanr , ng.xi
    
    #### SMALL GALAXIES ####
    
    source_cat=treecorr.Catalog(x=source_ra[~size_mask],y=source_dec[~size_mask],g1=source_e1[~size_mask],g2=-1.*source_e2[~size_mask],w=source_w[~size_mask],x_units='degree',y_units='degree')        
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    r_lens_l, xt_lens_l = ng.meanr , ng.xi
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(random_cat, source_cat)
    r_rand_l, xt_rand_l = ng.meanr , ng.xi    
    
    
    return r_lens_h, xt_lens_h, r_rand_h, xt_rand_h, r_lens_l, xt_lens_l, r_rand_l, xt_rand_l 



def lum_all():

    xshear_file = h5py.File("size_jk.h5" , "w")



    lum_jk_L1 = np.zeros((100,8,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h, r_lens_l, xt_lens_l, r_rand_l, xt_rand_l  = xshear_lens_jk(0.1, 0.3, i)
        
        lum_jk_L1[i,0,:] = r_lens_h
        lum_jk_L1[i,1,:] = xt_lens_h
        lum_jk_L1[i,2,:] = r_rand_h
        lum_jk_L1[i,3,:] = xt_rand_h
        
        lum_jk_L1[i,4,:] = r_lens_l
        lum_jk_L1[i,5,:] = xt_lens_l
        lum_jk_L1[i,6,:] = r_rand_l
        lum_jk_L1[i,7,:] = xt_rand_l
        
    
        print "done with", i
    xshear_file.create_dataset("lum_L1", (100,8,25) , data = lum_jk_L1)

    lum_jk_L2 = np.zeros((100,8,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h, r_lens_l, xt_lens_l, r_rand_l, xt_rand_l  = xshear_lens_jk(0.3, 0.5, i)
        
        lum_jk_L2[i,0,:] = r_lens_h
        lum_jk_L2[i,1,:] = xt_lens_h
        lum_jk_L2[i,2,:] = r_rand_h
        lum_jk_L2[i,3,:] = xt_rand_h
        
        lum_jk_L2[i,4,:] = r_lens_l
        lum_jk_L2[i,5,:] = xt_lens_l
        lum_jk_L2[i,6,:] = r_rand_l
        lum_jk_L2[i,7,:] = xt_rand_l
        
    
        print "done with", i
    xshear_file.create_dataset("lum_L2", (100,8,25) , data = lum_jk_L2)
    lum_jk_L3 = np.zeros((100,8,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h, r_lens_l, xt_lens_l, r_rand_l, xt_rand_l  = xshear_lens_jk(0.5, 0.7, i)
        
        lum_jk_L3[i,0,:] = r_lens_h
        lum_jk_L3[i,1,:] = xt_lens_h
        lum_jk_L3[i,2,:] = r_rand_h
        lum_jk_L3[i,3,:] = xt_rand_h
        
        lum_jk_L3[i,4,:] = r_lens_l
        lum_jk_L3[i,5,:] = xt_lens_l
        lum_jk_L3[i,6,:] = r_rand_l
        lum_jk_L3[i,7,:] = xt_rand_l
        
    
        print "done with", i

    xshear_file.create_dataset("lum_L3", (100,8,25) , data = lum_jk_L3)
 
    return None

if __name__ == '__main__':

   lum_all()
