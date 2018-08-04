import h5py
import numpy as np
import kmeans_radec
import matplotlib.pyplot as plt
import treecorr
import seaborn as sns
from kmeans_radec import KMeans, kmeans_sample

random_file = h5py.File("randoms.h5")
random_coord = random_file["random"][:] 
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

lens_file_jk = h5py.File("LRG_lum_jk.h5")
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


def xshear_lens_jk(zl1, zl2, jk_label):

    source_file = h5py.File("source_zb_"+str(zl2+0.1)+"_0.9.h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["psf1"][:]
    source_e2 = source_file["psf2"][:]
    source_w = source_file["w"][:]
    source_size = source_file["snr"][:]
    source_file.close()
    
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,g1=source_e1,g2=-1.*source_e2,w=source_w,x_units='degree',y_units='degree')        
    
    
    mask_z_lens = (lens_z_jk>zl1)&(lens_z_jk<zl2)
    lens_mask = labels_lens_jk == jk_label
    lens_mask = (~lens_mask)&(mask_z_lens)
    
    lens_cat = treecorr.Catalog(x=lens_ra_jk[lens_mask], y=lens_dec_jk[lens_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(source_cat, lens_cat)
    r_lens_h, xt_lens_h = ng.meanr , ng.xi 
    
    random_mask = labels_random_jk == jk_label
    random_cat = treecorr.Catalog(x=random_ra_jk[~random_mask], y=random_dec_jk[~random_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(source_cat, random_cat)
    r_rand_h, xt_rand_h = ng.meanr , ng.xi
    
    return r_lens_h, xt_lens_h, r_rand_h, xt_rand_h



def lum_all():

    xshear_file = h5py.File("flipper_jk.h5" , "w")



    lum_jk_L1 = np.zeros((100,4,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h = xshear_lens_jk(0.1, 0.3, i)
        
        lum_jk_L1[i,0,:] = r_lens_h
        lum_jk_L1[i,1,:] = xt_lens_h
        lum_jk_L1[i,2,:] = r_rand_h
        lum_jk_L1[i,3,:] = xt_rand_h
    
        print "done with", i
    xshear_file.create_dataset("lum_L1", (100,4,25) , data = lum_jk_L1)

    lum_jk_L2 = np.zeros((100,4,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h = xshear_lens_jk(0.3, 0.5, i)
        
        lum_jk_L2[i,0,:] = r_lens_h
        lum_jk_L2[i,1,:] = xt_lens_h
        lum_jk_L2[i,2,:] = r_rand_h
        lum_jk_L2[i,3,:] = xt_rand_h
    
        print "done with", i
    xshear_file.create_dataset("lum_L2", (100,4,25) , data = lum_jk_L2)
    lum_jk_L3 = np.zeros((100,4,25))

    for i in xrange(100):
    
        r_lens_h, xt_lens_h, r_rand_h, xt_rand_h = xshear_lens_jk(0.5, 0.7, i)
        
        lum_jk_L3[i,0,:] = r_lens_h
        lum_jk_L3[i,1,:] = xt_lens_h
        lum_jk_L3[i,2,:] = r_rand_h
        lum_jk_L3[i,3,:] = xt_rand_h
    
        print "done with", i

    xshear_file.create_dataset("lum_L3", (100,4,25) , data = lum_jk_L3)
 
    return None

if __name__ == '__main__':

   lum_all()
