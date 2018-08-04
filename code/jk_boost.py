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


def xshear_lens_jk(zl1, zl2, jk_label):

    
    mask_z_lens = (lens_z_jk>zl1)&(lens_z_jk<zl2)
    lens_mask = labels_lens_jk == jk_label
    lens_mask = (~lens_mask)&(mask_z_lens)
    random_mask = labels_random_jk == jk_label
    
    
    return 1.* len(random_ra_jk[~random_mask])/len(lens_ra_jk[lens_mask])



def lum_all():

    xshear_file = h5py.File("boost_lum_jk.h5" , "w")



    lum_jk_L1 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.1, 0.3, i)
        lum_jk_L1[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L1", (100,) , data = lum_jk_L1)

    lum_jk_L2 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.3, 0.5, i)
        lum_jk_L2[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L2", (100,) , data = lum_jk_L2)

    lum_jk_L3 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.5, 0.7, i)
        lum_jk_L3[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L3", (100,) , data = lum_jk_L3)
 
    return None


def dense_all():

    xshear_file = h5py.File("boost_dense_jk.h5" , "w")

    lum_jk_L1 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.1, 0.3, i)
        lum_jk_L1[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L1", (100,) , data = lum_jk_L1)

    lum_jk_L2 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.3, 0.5, i)
        lum_jk_L2[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L2", (100,) , data = lum_jk_L2)

    lum_jk_L3 = np.zeros((100))

    for i in xrange(100):
    
        boost = xshear_lens_jk(0.5, 0.7, i)
        lum_jk_L3[i] = boost
        print "done with", i
    
    xshear_file.create_dataset("lum_L3", (100,) , data = lum_jk_L3)
 
    return None


if __name__ == '__main__':

   lum_all()
   dense_all()
