import pyfits as pf
import h5py
import numpy as np
import kmeans_radec
import matplotlib.pyplot as plt
import treecorr
import seaborn as sns
from kmeans_radec import KMeans, kmeans_sample
from astropy.stats import median_absolute_deviation
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.switch_backend("Agg")

#############

random_file = h5py.File("randoms_jk.h5")
random_coord = random_file["random_jk"][:] 
random_ra = random_coord[:,0]
random_dec = random_coord[:,1]
random_file.close()

#############

def xshear(zl1,zl2,zs1,zs2,lens_type):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_m = source_file["bias"][:]
    source_w = source_file["w"][:]
    source_file.close()
    
    #print "diagnosis" , source_z.min(), source_z.max()
    
    if lens_type == 'lum':
        lens_file = h5py.File("LRG_lum_jk.h5")
    
    elif lens_type == 'dense':
        lens_file = h5py.File("LRG_dense_jk.h5")
    
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_file.close()
    
    mask_z_lens = (lens_z>zl1)&(lens_z<zl2)
    
    lens_cat = treecorr.Catalog(x=lens_ra[mask_z_lens], y=lens_dec[mask_z_lens], x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,k =source_m, w=source_w,x_units='degree',y_units='degree')
    
    ng = treecorr.NKCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t,  weight = ng.meanr , ng.xi, ng.weight
    
    return r, xi_t, weight


#########################

def xshear_random(zs1,zs2):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_file.close()
    
    #print "diagnosis" , source_z.min(), source_z.max()
    
    random_file = h5py.File("randoms_jk.h5")
    random_coord = random_file["random_jk"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]
    
    #shuffled = np.random.choice(len(random_ra) , 2000000)
    
    lens_cat = treecorr.Catalog(x=random_ra, y=random_dec, x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,g1=source_e1,g2=-1.*source_e2,w=source_w,x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x,  weight = ng.meanr , ng.xi, ng.xi_im , ng.weight
    
    return r, xi_t, xi_x, weight
############################

def flipper():

    source_file = h5py.File("source_zb_0.1_0.3.h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_file.close()
    
    lens_file = h5py.File("LRG_dense_jk.h5")
    
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_file.close()
    
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    #print "diagnosis" , source_z.min(), source_z.max()
    
    source_mask = (source_z>0.1)&(source_z<0.3)

    source_cat=treecorr.Catalog(x=source_ra[source_mask],y=source_dec[source_mask],g1=source_e1[source_mask],g2=-1.*source_e2[source_mask],w=source_w[source_mask],x_units='degree',y_units='degree')        
    
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    lens_mask = mask_z_lens
    
    lens_cat = treecorr.Catalog(x=lens_ra[lens_mask], y=lens_dec[lens_mask], x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    r_lens_h, xt_lens_h, xx_lens_h , w_lens_h = ng.meanr , ng.xi , ng.xi_im, ng.weight
   
    print "lens", xt_lens_h

    random_cat = treecorr.Catalog(x=random_ra, y=random_dec, x_units='degree', y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(random_cat, source_cat)
    r_rand_h, xt_rand_h, xx_rand_h , w_rand_h = ng.meanr , ng.xi, ng.xi_im, ng.weight

    print "random", xt_rand_h 

    return r_lens_h, 10**5.*xt_lens_h, 10**5.*xx_lens_h , w_lens_h, r_rand_h, 10**5.*xt_rand_h, 10**5.*xx_rand_h , w_rand_h


############################

def xshear_flip():
    
    xshear_file = h5py.File("xshear_lens_henk.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw, rr, xtr, xxr, xwr = flipper()
    xshear_file.create_dataset("lum_L1", (8, nbins) , data = np.vstack([r, xt, xx, xw, rr, xtr, xxr, xwr]))
    xshear_file.close()
    
    return None

##############################



def xshear_all():
    
    xshear_file = h5py.File("xshear_lens_henk.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw = xshear(0.1,0.3,0.4,0.9,"lum")
    xshear_file.create_dataset("lum_L1", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
   
    r , xt, xx , xw  = xshear(0.3,0.5,0.6,0.9,"lum")
    xshear_file.create_dataset("lum_L2", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    r , xt, xx , xw = xshear(0.5,0.7,0.8,0.9,"lum")
    xshear_file.create_dataset("lum_L3", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
 
    ######### dense ###########
    
    r , xt, xx , xw = xshear(0.1,0.3,0.4,0.9,"dense")
    xshear_file.create_dataset("dense_L1", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
   
    r , xt, xx , xw  = xshear(0.3,0.5,0.6,0.9,"dense")
    xshear_file.create_dataset("dense_L2", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    r , xt, xx, xw  = xshear(0.5,0.7,0.8,0.9,"dense")
    xshear_file.create_dataset("dense_L3", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    xshear_file.close()
    
    return None


def xshear_random_all():
    
    nbins = 25
    xshear_file = h5py.File("xshear_random_henk.h5" , "w")

    ######### lum ###########
    
    r , xt, xx, xw  = xshear_random(0.4,0.9)
    xshear_file.create_dataset("rand_S1", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    r , xt, xx, xw  = xshear_random(0.6,0.9)
    xshear_file.create_dataset("rand_S2", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    r , xt, xx, xw  = xshear_random(0.8,0.9)
    xshear_file.create_dataset("rand_S3", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    xshear_file.close()


##################################

def plot_xshear_dense():
    
    dense_jk = h5py.File("dense_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]

    jk_L3 = dense_jk["lum_L3"][:]
   
    print jk_L1.shape
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["dense_L1"][:]

    L2 = xshear_file["dense_L2"][:]
  
    L3= xshear_file["dense_L3"][:]

    fig = plt.figure(figsize=(15,8))


    err = np.std((jk_L1[:,0,:]/100)**0.5 *(jk_L1[:,2,:]-jk_L1[:,6,:]) , axis = 0)
    #ax4.errorbar(L1[0,:] , (L1[0,:]/100)**0.5 * (L1[2,:] - rand_S1[2,:]),  yerr = err , elinewidth = 2, fmt='.',capsize = 2)
    
    for i in range(100):
        print jk_L1[i,2,:]-jk_L1[i,6,:]
        plt.plot(jk_L1[i,0,:] , 10**-5. * (jk_L1[i,2,:]-jk_L1[i,6,:]) , alpha = .1, lw = 1. , color = "r")
    print np.std(10**-5.*(jk_L1[:,2,:]-jk_L1[:,6,:]) , axis = 0)
    plt.xscale("log")

    plt.savefig("/home/vakili/public_html/dense.png")
    return None

def plot_xshear_lum_p():
   
    fct = (100.-1.)/(100. - 18 -1.)

    dense_jk = h5py.File("lum_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]
    jk_L3 = dense_jk["lum_L3"][:]
   
    print jk_L1.shape
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["lum_L1"][:]

    L2 = xshear_file["lum_L2"][:]
  
    L3= xshear_file["lum_L3"][:]

    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,8))
    gs = gridspec.GridSpec(2, 3, height_ratios=[3,1] , wspace = 0 , hspace = 0.0)

    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[0,2])
    ax4 = plt.subplot(gs[1,0])
    ax5 = plt.subplot(gs[1,1])
    ax6 = plt.subplot(gs[1,2])    
    #fig , ax = plt.subplots(nrows=2,ncols=3 , sharex='col', sharey='row')

    ############# L1 #################
    

    cov = fct* 99.**2. * np.cov((jk_L1[:,1,:]-jk_L1[:,5,:]).T)/100.
    err = np.diag(cov)**0.5
    print err.shape 
    print np.log(-20.0) 
    ax1.errorbar(L1[0,:] , L1[1,:] - rand_S1[1,:], err , fmt='o', capsize = 2 ,  label = r"$0.4<z_B<0.9$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    
    ax1.set_ylabel(r"$\gamma_{T}$" , fontsize = 30)
    ax1.set_ylim([5.*10**-6. , 0.07])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    cov =  fct* 99.**2. * np.cov((jk_L2[:,1,:]-jk_L2[:,5,:]).T)/100.
    err = np.diag(cov)**0.5
    
    ax2.errorbar(L2[0,:] , L2[1,:] - rand_S2[1,:], err,  fmt='o', capsize = 2 , label = r"$0.6<z_B<0.9$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")     
    ax2.legend(fontsize = 15)
   
    ax2.set_ylim([5.*10**-6. , 0.07])
    ax2.set_xlim([0.28 , 305])    
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    
     
    #for i in range(40):
        
    cov =  fct* 99.**2. * np.cov((jk_L3[:,1,:]-jk_L3[:,5,:]).T)/100.
    err = np.diag(cov)**0.5

    
    ax3.errorbar(L3[0,:] , L3[1,:] - rand_S3[1,:], err,  fmt='o', capsize = 2 , label = r"$0.8<z_B<0.9$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")    
    ax3.set_ylim([5.*10**-6. , 0.07])
    ax3.set_xlim([0.28 , 305])    
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.legend(fontsize = 15)

    
    #cov = 99.**2. * np.cov((jk_L1[:,2,:]-jk_L1[:,6,:]).T)/100.
    #cov = 99.**2. * np.cov((jk_L1[:,2,:]-jk_L1[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L1[:,2,:]-jk_L1[:,6,:])*((L1[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    print (L1[0,:]/100.)**0.5 *(L1[2,:]-rand_S1[2,:]) 
    ax4.errorbar(L1[0,:] , (L1[0,:]/100.)**0.5 * (L1[2,:]-rand_S1[2,:]) ,  err,  fmt='o', capsize = 2)
    
    ax4.set_xscale("log")
    ax4.tick_params(axis = "x", labelsize = 20)
    ax4.tick_params(axis = "y", labelsize = 20)

    ax4.set_ylabel(r"$ (\theta /100)^{0.5} \; \gamma_{\times}$" , fontsize = 20)
    ax4.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)

    ax4.set_ylim([-0.00035 , 0.00035])
    ax4.set_xlim([0.28 , 305])
    
    #err = np.std(jk_L2[:,2,:]-jk_L2[:,6,:] , axis = 0)
    #cov = 99.**2. * np.cov((jk_L2[:,2,:]-jk_L2[:,6,:]).T)/100.
    #cov = 99.**2. * np.cov((jk_L2[:,2,:]-jk_L2[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L2[:,2,:]-jk_L2[:,6,:])*((L2[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    ax5.errorbar(L2[0,:] , (L2[0,:]/100.)**0.5 * (L2[2,:]-rand_S2[2,:]) , err,  fmt='o', capsize = 2)
    ax5.set_xscale("log")

    chisq = np.dot((L2[2,2:]-rand_S2[2,2:]),np.linalg.solve(cov[2:,2:],(L2[2,2:]-rand_S2[2,2:])))
    print chisq/18.
    ax5.set_ylim([-0.00035 , 0.00035])
    ax5.set_xlim([0.28 , 305])
    ax5.set_yticklabels([])
    ax5.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax5.tick_params(axis = "x", labelsize = 20)

    #cov = 99.**2. * np.cov((jk_L3[:,2,:]-jk_L3[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L3[:,2,:]-jk_L3[:,6,:])*((L3[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    ax6.errorbar(L3[0,:] , (L3[0,:]/100.)**0.5 * (L3[2,:]-rand_S3[2,:]) , err,  fmt='o', capsize = 2)

    chisq = np.dot((L3[2,1:]-rand_S3[2,1:]),np.linalg.solve(cov[1:,1:],(L3[2,1:]-rand_S3[2,1:])))
    print chisq/18.

    ax6.set_xscale("log")
    ax6.set_ylim([-0.00035 , 0.00035])
    ax6.set_xlim([0.28 , 305])
    ax6.set_yticklabels([])
    ax6.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax6.tick_params(axis = "x", labelsize = 20)

    ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    ax2.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    ax3.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    ax4.axhline(0.0, color = "k" , linestyle="dashed")
    ax5.axhline(0.0, color = "k" , linestyle="dashed")
    ax6.axhline(0.0, color = "k" , linestyle="dashed")

    fig.tight_layout()
   
    plt.savefig("/home/vakili/public_html/lum.png")
    return None

def plot_xshear_dense_p():
   
    fct = (100.-1.)/(100. - 18 -1.)

    dense_jk = h5py.File("dense_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]*10**-5.
    jk_L2 = dense_jk["lum_L2"][:]*10**-5.
    jk_L3 = dense_jk["lum_L3"][:] *10**-5.
   
    print jk_L1.shape
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["dense_L1"][:]

    L2 = xshear_file["dense_L2"][:]
  
    L3= xshear_file["dense_L3"][:]

    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,8))
    gs = gridspec.GridSpec(2, 3, height_ratios=[3,1] , wspace = 0 , hspace = 0.0)

    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[0,2])
    ax4 = plt.subplot(gs[1,0])
    ax5 = plt.subplot(gs[1,1])
    ax6 = plt.subplot(gs[1,2])    
    #fig , ax = plt.subplots(nrows=2,ncols=3 , sharex='col', sharey='row')

    ############# L1 #################
    

    cov = fct* 99.**2. * np.cov((jk_L1[:,1,:]-jk_L1[:,5,:]).T)/100.
    err = np.diag(cov)**0.5
    print err.shape 
    print np.log(-20.0) 
    ax1.errorbar(L1[0,:] , L1[1,:] - rand_S1[1,:], err , fmt='o', capsize = 2 ,  label = r"$0.4<z_B<0.9$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    
    ax1.set_ylabel(r"$\gamma_{T}$" , fontsize = 30)
    ax1.set_ylim([5.*10**-6. , 0.07])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    cov =  fct* 99.**2. * np.cov((jk_L2[:,1,:]-jk_L2[:,5,:]).T)/100.
    err = np.diag(cov)**0.5
    
    ax2.errorbar(L2[0,:] , L2[1,:] - rand_S2[1,:], err,  fmt='o', capsize = 2 , label = r"$0.6<z_B<0.9$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")     
    ax2.legend(fontsize = 15)
   
    ax2.set_ylim([5.*10**-6. , 0.07])
    ax2.set_xlim([0.28 , 305])    
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    
     
    #for i in range(40):
        
    cov =  fct* 99.**2. * np.cov((jk_L3[:,1,:]-jk_L3[:,5,:]).T)/100.
    err = np.diag(cov)**0.5

    
    ax3.errorbar(L3[0,:] , L3[1,:] - rand_S3[1,:], err,  fmt='o', capsize = 2 , label = r"$0.8<z_B<0.9$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")    
    ax3.set_ylim([5.*10**-6. , 0.07])
    ax3.set_xlim([0.28 , 305])    
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.legend(fontsize = 15)

    
    #cov = 99.**2. * np.cov((jk_L1[:,2,:]-jk_L1[:,6,:]).T)/100.
    #cov = 99.**2. * np.cov((jk_L1[:,2,:]-jk_L1[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L1[:,2,:]-jk_L1[:,6,:])*((L1[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    print (L1[0,:]/100.)**0.5 *(L1[2,:]-rand_S1[2,:]) 
    ax4.errorbar(L1[0,:] , (L1[0,:]/100.)**0.5 * (L1[2,:]-rand_S1[2,:]) ,  err,  fmt='o', capsize = 2)
    
    ax4.set_xscale("log")
    ax4.tick_params(axis = "x", labelsize = 20)
    ax4.tick_params(axis = "y", labelsize = 20)

    ax4.set_ylabel(r"$ (\theta /100)^{0.5} \; \gamma_{\times}$" , fontsize = 20)
    ax4.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)

    ax4.set_ylim([-0.00035 , 0.00035])
    ax4.set_xlim([0.28 , 305])
    
    #err = np.std(jk_L2[:,2,:]-jk_L2[:,6,:] , axis = 0)
    #cov = 99.**2. * np.cov((jk_L2[:,2,:]-jk_L2[:,6,:]).T)/100.
    #cov = 99.**2. * np.cov((jk_L2[:,2,:]-jk_L2[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L2[:,2,:]-jk_L2[:,6,:])*((L2[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    ax5.errorbar(L2[0,:] , (L2[0,:]/100.)**0.5 * (L2[2,:]-rand_S2[2,:]) , err,  fmt='o', capsize = 2)
    ax5.set_xscale("log")

    chisq = np.dot((L2[2,2:]-rand_S2[2,2:]),np.linalg.solve(cov[2:,2:],(L2[2,2:]-rand_S2[2,2:])))
    print chisq/18.
    ax5.set_ylim([-0.00035 , 0.00035])
    ax5.set_xlim([0.28 , 305])
    ax5.set_yticklabels([])
    ax5.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax5.tick_params(axis = "x", labelsize = 20)

    #cov = 99.**2. * np.cov((jk_L3[:,2,:]-jk_L3[:,6,:]).T)/100.
    cov =  fct* 99.**2. * np.cov(((jk_L3[:,2,:]-jk_L3[:,6,:])*((L3[0,:][None,:]/100.)**0.5)).T)/100.
    err = np.diag(cov)**0.5
    ax6.errorbar(L3[0,:] , (L3[0,:]/100.)**0.5 * (L3[2,:]-rand_S3[2,:]) , err,  fmt='o', capsize = 2)

    chisq = np.dot((L3[2,1:]-rand_S3[2,1:]),np.linalg.solve(cov[1:,1:],(L3[2,1:]-rand_S3[2,1:])))
    print chisq/18.

    ax6.set_xscale("log")
    ax6.set_ylim([-0.00035 , 0.00035])
    ax6.set_xlim([0.28 , 305])
    ax6.set_yticklabels([])
    ax6.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax6.tick_params(axis = "x", labelsize = 20)

    ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    ax2.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    ax3.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    ax4.axhline(0.0, color = "k" , linestyle="dashed")
    ax5.axhline(0.0, color = "k" , linestyle="dashed")
    ax6.axhline(0.0, color = "k" , linestyle="dashed")

    fig.tight_layout()
   
    plt.savefig("/home/vakili/public_html/dense.png")
    return None

##########################################################################


def plot_xshear_lum():
    
    dense_jk = h5py.File("lum_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]
    jk_L3 = dense_jk["lum_L3"][:]    
    
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["lum_L1"][:]

    L2 = xshear_file["lum_L2"][:]
  
    L3= xshear_file["lum_L3"][:]

    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,8))
    gs = gridspec.GridSpec(2, 3, height_ratios=[3,1] , wspace = 0 , hspace = 0.0)

    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[0,1])
    ax3 = plt.subplot(gs[0,2])
    ax4 = plt.subplot(gs[1,0])
    ax5 = plt.subplot(gs[1,1])
    ax6 = plt.subplot(gs[1,2])    
    #fig , ax = plt.subplots(nrows=2,ncols=3 , sharex='col', sharey='row')

    ############# L1 #################
    
    

    err = np.std(jk_L1[:,1,:]-jk_L1[:,5,:] , axis = 0)

    
    ax1.errorbar(L1[0,:] , L1[1,:] - rand_S1[1,:], err , fmt='o', capsize = 2 , label = r"$0.4<z_B<0.9$")

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    
    ax1.set_ylabel(r"$\gamma_{T}$" , fontsize = 30)
    ax1.set_ylim([5.*10**-6. , 0.07])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    err = np.std(jk_L2[:,1,:]-jk_L2[:,5,:] , axis = 0)

  
    
    ax2.errorbar(L2[0,:] , L2[1,:] - rand_S2[1,:], err,  fmt='o', capsize = 2 , label = r"$0.6<z_B<0.9$")
    ax2.set_xscale("log")
    ax2.set_yscale("log")     
    ax2.legend(fontsize = 15)
   
    ax2.set_ylim([5.*10**-6. , 0.07])
    ax2.set_xlim([0.28 , 305])    
    ax2.set_yticklabels([])
    ax2.set_xticklabels([])
    
     
    #for i in range(40):
        
    err = np.std(jk_L3[:,1,:]-jk_L3[:,5,:] , axis = 0)

    
    ax3.errorbar(L3[0,:] , L3[1,:] - rand_S3[1,:], err,  fmt='o', capsize = 2 , label = r"$0.8<z_B<0.9$")
    ax3.set_xscale("log")
    ax3.set_yscale("log")    
    ax3.set_ylim([5.*10**-6. , 0.07])
    ax3.set_xlim([0.28 , 305])    
    ax3.set_yticklabels([])
    ax3.set_xticklabels([])
    ax3.legend(fontsize = 15)

    err = np.std(jk_L1[:,2,:]-jk_L1[:,6,:] , axis = 0)
    ax4.errorbar(L1[0,:] , (L1[0,:]/100)**0.5 * (L1[2,:] - rand_S1[2,:]), (L1[0,:]/100)**0.5 * err , err, fmt='o',capsize = 2)
    ax4.set_xscale("log")
    ax4.tick_params(axis = "x", labelsize = 20)
    ax4.tick_params(axis = "y", labelsize = 20)

    ax4.set_ylabel(r"$ (\theta /100)^{0.5} \; \gamma_{\times}$" , fontsize = 20)
    ax4.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)

    ax4.set_ylim([-0.00035 , 0.00035])
    ax4.set_xlim([0.28 , 305])
    
    err = np.std(jk_L2[:,2,:]-jk_L2[:,6,:] , axis = 0)
    ax5.errorbar(L2[0,:] , (L2[0,:]/100)**0.5 * (L2[2,:] - rand_S2[2,:]), (L2[0,:]/100)**0.5 * err, fmt='o',capsize = 2 )
    ax5.set_xscale("log")

    ax5.set_ylim([-0.00035 , 0.00035])
    ax5.set_xlim([0.28 , 305])
    ax5.set_yticklabels([])
    ax5.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax5.tick_params(axis = "x", labelsize = 20)

    err = np.std(jk_L3[:,2,:]-jk_L3[:,6,:] , axis = 0)
    ax6.errorbar(L3[0,:] , (L3[0,:]/100)**0.5 * (L3[2,:] - rand_S3[2,:]), (L3[0,:]/100)**0.5 * err,fmt='o', capsize = 2)
    ax6.set_xscale("log")

    
    ax6.set_ylim([-0.00035 , 0.00035])
    ax6.set_xlim([0.28 , 305])
    ax6.set_yticklabels([])
    ax6.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax6.tick_params(axis = "x", labelsize = 20)

    ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    ax2.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    ax3.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None

############################################################

def plot_boost():
    
    
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,8))
    gs = gridspec.GridSpec(1, 3 , wspace = 0 , hspace = 0.0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
 
    
    ############## dense #####################
    
    lens_file = h5py.File("LRG_dense.h5")
    lens_z = lens_file["redshift"][:]

    dense_jk = h5py.File("dense_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]
    jk_L3 = dense_jk["lum_L3"][:]    
    
    boost_jk = h5py.File("boost_dense_jk.h5")
    boost_jk_L1 = boost_jk["lum_L1"][:]
    boost_jk_L2 = boost_jk["lum_L2"][:]
    boost_jk_L3 = boost_jk["lum_L3"][:]
    
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    L1 = xshear_file["dense_L1"][:]
    L2 = xshear_file["dense_L2"][:]
    L3= xshear_file["dense_L3"][:]

    ############# L1 ###############
    mask_z_lens = (lens_z>0.1)&(lens_z<0.3)
    
    
    errs = boost_jk_L1[:,None]*jk_L1[:,3,:]/jk_L1[:,7,:]
    err = np.std(errs, axis = 0)
    ax1.errorbar(L1[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L1[3,:]/rand_S1[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{dense}\;\;0.4<z_B<0.9$")

    ax1.set_xscale("log")
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"Boost factor" , fontsize = 30)
    ax1.set_ylim([0.98,1.2])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    ############ L2 ##########################
    mask_z_lens = (lens_z>0.3)&(lens_z<0.5)
    errs = boost_jk_L2[:,None]*jk_L2[:,3,:]/jk_L2[:,7,:]
    err = np.std(errs, axis = 0)
    ax2.errorbar(L2[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L2[3,:]/rand_S2[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{dense}\;\;0.6<z_B<0.9$")
    ax2.set_xscale("log")
    ax2.legend(fontsize = 15)
 
    ax2.set_ylim([0.98,1.2])
    ax2.set_xlim([0.28 , 305])    
    ax2.set_yticklabels([])
    
    #############L3###########################    
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    errs = boost_jk_L3[:,None]*jk_L3[:,3,:]/jk_L3[:,7,:]
    err = np.std(errs, axis = 0)
    ax3.errorbar(L3[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L3[3,:]/rand_S3[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{dense}\;\;0.8<z_B<0.9$")
    ax3.set_xscale("log")
    ax3.legend(fontsize = 15)
 
    ax3.set_ylim([0.97,1.2])
    ax3.set_xscale("log")
    ax3.set_xlim([0.28 , 305])    
    ax3.set_yticklabels([])
    ax3.legend(fontsize = 15)


    ############## lum #####################
    
    lens_file = h5py.File("LRG_lum.h5")
    lens_z = lens_file["redshift"][:]

    lum_jk = h5py.File("lum_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]
    jk_L3 = dense_jk["lum_L3"][:]    
    
    boost_jk = h5py.File("boost_lum_jk.h5")
    boost_jk_L1 = boost_jk["lum_L1"][:]
    boost_jk_L2 = boost_jk["lum_L2"][:]
    boost_jk_L3 = boost_jk["lum_L3"][:]
    
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    L1 = xshear_file["lum_L1"][:]
    L2 = xshear_file["lum_L2"][:]
    L3= xshear_file["lum_L3"][:]    

    
    ############# L1 ###############
    mask_z_lens = (lens_z>0.1)&(lens_z<0.3)
    errs = boost_jk_L1[:,None]*jk_L1[:,3,:]/jk_L1[:,7,:]
    err = np.std(errs, axis = 0)
    ax1.errorbar(L1[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L1[3,:]/rand_S1[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{luminous}\;\;0.4<z_B<0.9$")

    ax1.set_xscale("log")
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"Boost factor" , fontsize = 30)
    ax1.set_ylim([0.97,1.2])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    ############ L2 ##########################
    mask_z_lens = (lens_z>0.3)&(lens_z<0.5)
    errs = boost_jk_L2[:,None]*jk_L2[:,3,:]/jk_L2[:,7,:]
    err = np.std(errs, axis = 0)
    ax2.errorbar(L2[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L2[3,:]/rand_S2[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{luminous}\;\;0.6<z_B<0.9$")
    ax2.set_xscale("log")
    ax2.legend(fontsize = 15)
 
    ax2.set_ylim([0.97,1.2])
    ax2.set_xlim([0.28 , 305])    
    ax2.set_yticklabels([])
    
    #############L3###########################    
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    errs = boost_jk_L3[:,None]*jk_L3[:,3,:]/jk_L3[:,7,:]
    err = np.std(errs, axis = 0)
    ax3.errorbar(L3[0,:],(len(random_ra)*1.0/len(lens_z[mask_z_lens]))*L3[3,:]/rand_S3[3,:],err,fmt='o',capsize =2, label = r"$\mathrm{luminous}\;\; 0.8<z_B<0.9$")
    ax3.set_xscale("log")
    ax3.legend(fontsize = 15)
 
    ax1.axhline(1, color="k", linestyle = "dashed")
    ax2.axhline(1, color="k", linestyle = "dashed")
    ax3.axhline(1, color="k", linestyle = "dashed")

    ax1.tick_params(axis = "x", labelsize = 20)
    ax2.tick_params(axis = "x", labelsize = 20)
    ax3.tick_params(axis = "x", labelsize = 20)
    
    ax1.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax2.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax3.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)

    ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    ax2.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    ax3.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None


def xshear_lens_snr(zl1,zl2,zs1,zs2,lens_type):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_snr = source_file["snr"][:]
    snr_mask = source_snr > np.median(source_snr)
    source_file.close()
   
    if lens_type == 'lum':
        lens_file = h5py.File("LRG_lum.h5")
    
    elif lens_type == 'dense':
        lens_file = h5py.File("LRG_dense.h5")
    
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_file.close()
    
    mask_z_lens = (lens_z>zl1)&(lens_z<zl2)
    
    lens_cat = treecorr.Catalog(x=lens_ra[mask_z_lens], y=lens_dec[mask_z_lens], x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra[snr_mask],y=source_dec[snr_mask],g1=source_e1[snr_mask],g2=-1.*source_e2[snr_mask],w=source_w[snr_mask],x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x , weight = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    source_cat=treecorr.Catalog(x=source_ra[~snr_mask],y=source_dec[~snr_mask],g1=source_e1[~snr_mask],g2=-1.*source_e2[~snr_mask],w=source_w[~snr_mask],x_units='degree',y_units='degree')
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    rl, xi_tl, xi_xl , weightl = ng.meanr , ng.xi, ng.xi_im, ng.weight

    
    return r, xi_t, xi_x, weight, rl, xi_tl, xi_xl, weightl

def xshear_random_snr(zs1,zs2):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_snr = source_file["snr"][:]
    source_file.close()
    
    snr_mask = source_snr > np.median(source_snr)
    source_file.close()
        
    random_file = h5py.File("randoms.h5")
    random_coord = random_file["random"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]
        
    lens_cat = treecorr.Catalog(x=random_ra, y=random_dec, x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra[snr_mask],y=source_dec[snr_mask],g1=source_e1[snr_mask],g2=-1.*source_e2[snr_mask],w=source_w[snr_mask],x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x , weight = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    source_cat=treecorr.Catalog(x=source_ra[~snr_mask],y=source_dec[~snr_mask],g1=source_e1[~snr_mask],g2=-1.*source_e2[~snr_mask],w=source_w[~snr_mask],x_units='degree',y_units='degree')
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)

    ng.process(lens_cat, source_cat)
    rl, xi_tl, xi_xl , weightl = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    return r, xi_t, xi_x, weight, rl, xi_tl, xi_xl, weightl


def xshear_snr_all():
    
    xshear_file = h5py.File("xshear_lens_henk_snr.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_snr(0.1,0.3,0.4,0.9,"lum")
    xshear_file.create_dataset("lum_L1", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
   
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_snr(0.3,0.5,0.6,0.9,"lum")
    xshear_file.create_dataset("lum_L2", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_snr(0.5,0.7,0.8,0.9,"lum")
    xshear_file.create_dataset("lum_L3", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
 
    xshear_file.close()
    
    return None

def xshear_snr_random_all():
    
    xshear_file = h5py.File("xshear_random_henk_snr.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_snr(0.4,0.9)
    xshear_file.create_dataset("lum_L1", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
   
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_snr(0.6,0.9)
    xshear_file.create_dataset("lum_L2", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_snr(0.8,0.9)
    xshear_file.create_dataset("lum_L3", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
 
    xshear_file.close()
    
    return None

def xshear_lens_size(zl1,zl2,zs1,zs2,lens_type):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_snr = source_file["rad"][:]
    snr_mask = source_snr > np.median(source_snr)
    source_file.close()
    
    #print "diagnosis" , source_z.min(), source_z.max()
    
    if lens_type == 'lum':
        lens_file = h5py.File("LRG_lum.h5")
    
    elif lens_type == 'dense':
        lens_file = h5py.File("LRG_dense.h5")
    
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_file.close()
    
    mask_z_lens = (lens_z>zl1)&(lens_z<zl2)
    
    lens_cat = treecorr.Catalog(x=lens_ra[mask_z_lens], y=lens_dec[mask_z_lens], x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra[snr_mask],y=source_dec[snr_mask],g1=source_e1[snr_mask],g2=-1.*source_e2[snr_mask],w=source_w[snr_mask],x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x , weight = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    source_cat=treecorr.Catalog(x=source_ra[~snr_mask],y=source_dec[~snr_mask],g1=source_e1[~snr_mask],g2=-1.*source_e2[~snr_mask],w=source_w[~snr_mask],x_units='degree',y_units='degree')
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    rl, xi_tl, xi_xl , weightl = ng.meanr , ng.xi, ng.xi_im, ng.weight

    
    return r, xi_t, xi_x, weight, rl, xi_tl, xi_xl, weightl

def xshear_random_size(zs1,zs2):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["e1"][:]
    source_e2 = source_file["e2"][:]
    source_w = source_file["w"][:]
    source_snr = source_file["rad"][:]
    source_file.close()
    
    snr_mask = source_snr > np.median(source_snr)
    source_file.close()
        
    random_file = h5py.File("randoms.h5")
    random_coord = random_file["random"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]
        
    lens_cat = treecorr.Catalog(x=random_ra, y=random_dec, x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra[snr_mask],y=source_dec[snr_mask],g1=source_e1[snr_mask],g2=-1.*source_e2[snr_mask],w=source_w[snr_mask],x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x , weight = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    source_cat=treecorr.Catalog(x=source_ra[~snr_mask],y=source_dec[~snr_mask],g1=source_e1[~snr_mask],g2=-1.*source_e2[~snr_mask],w=source_w[~snr_mask],x_units='degree',y_units='degree')
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    rl, xi_tl, xi_xl , weightl = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    return r, xi_t, xi_x, weight, rl, xi_tl, xi_xl, weightl

def xshear_size_all():
    
    xshear_file = h5py.File("xshear_lens_henk_size.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_size(0.1,0.3,0.4,0.9,"lum")
    xshear_file.create_dataset("lum_L1", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
   
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_size(0.3,0.5,0.6,0.9,"lum")
    xshear_file.create_dataset("lum_L2", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_lens_size(0.5,0.7,0.8,0.9,"lum")
    xshear_file.create_dataset("lum_L3", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
 
    xshear_file.close()
    
    return None


def xshear_size_random_all():
    
    xshear_file = h5py.File("xshear_random_henk_size.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_size(0.4,0.9)
    xshear_file.create_dataset("lum_L1", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
   
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_size(0.6,0.9)
    xshear_file.create_dataset("lum_L2", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
    
    r , xt, xx , xw, rl , xtl, xxl , xwl = xshear_random_size(0.8,0.9)
    xshear_file.create_dataset("lum_L3", (8, nbins) , data = np.vstack([ r , xt, xx , xw, rl , xtl, xxl , xwl]))
 
    xshear_file.close()
    
    return None

###############################################################################################

def plot_snr():
    
    
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(15,5))
    gs = gridspec.GridSpec(1, 3 , wspace = 0 , hspace = 0.0)

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])

    ############## lum #####################
    
    lens_file = h5py.File("LRG_lum.h5")
    lens_z = lens_file["redshift"][:]

    lum_jk = h5py.File("lum_jk.h5")
    jk_L1 = dense_jk["lum_L1"][:]
    jk_L2 = dense_jk["lum_L2"][:]
    jk_L3 = dense_jk["lum_L3"][:]   
    
    
    
    #snr_jk = h5py.File("snr_jk.h5")
    #jk_snr_L1 = dense_jk["lum_L1"][:]
    #jk_snr_L2 = dense_jk["lum_L2"][:]
    #jk_snr_L3 = dense_jk["lum_L3"][:]   
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    L1 = xshear_file["lum_L1"][:]
    L2 = xshear_file["lum_L2"][:]
    L3= xshear_file["lum_L3"][:] 
    
    xshear_snr_file = h5py.File("xshear_lens_henk_snr.h5")
    snr_L1 = xshear_snr_file["lum_L1"][:]
    snr_L2 = xshear_snr_file["lum_L2"][:]
    snr_L3 = xshear_snr_file["lum_L3"][:]
 
    random_snr_file = h5py.File("xshear_random_henk_snr.h5")
    snr_rand_S1 = random_snr_file["lum_L1"][:]
    snr_rand_S2 = random_snr_file["lum_L2"][:]
    snr_rand_S3 = random_snr_file["lum_L3"][:]
    
    ############# L1 ###############
    mask_z_lens = (lens_z>0.1)&(lens_z<0.3)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    
    x = L1[0,:]
    yall = L1[1,:] - rand_S1[1,:]
    wall = rand_S1[3,:]/L1[3,:]
    yall = yall*wall
    
    yh = snr_L1[1,:] - snr_rand_S1[1,:]
    wh = snr_rand_S1[3,:]/snr_L1[3,:]
    mask_z_lens = (lens_z>0.1)&(lens_z<0.3)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yh = yh * ratio / wh
    
    yl = snr_L1[5,:] - snr_rand_S1[5,:]
    wl = snr_rand_S1[7,:]/snr_L1[7,:]
    mask_z_lens = (lens_z>0.1)&(lens_z<0.3)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yl = yl * ratio / wl

    err_all = np.std(jk_L1[:,1,:] - jk_L1[:,5,:], axis =0)
    

    var_ha = (yh/yall)**2. * (err_all**2/yh**2 + err_all**2/yall**2)
    var_la = (yl/yall)**2. * (err_all**2/yl**2 + err_all**2/yall**2)
    var_tot = var_ha + var_la
    err = var_tot**0.5
    
    #ax1.errorbar(x[yh<0],np.abs(yh[yh<0]),err_all[yh<0],fmt='o',capsize =2)
    #ax1.errorbar(x[yl<0],np.abs(yl[yl<0]),err_all[yl<0],fmt='o',capsize =2)
    
    ax1.errorbar(x,yh,err_all,fmt='o',capsize =2)
    ax1.errorbar(x,yl,err_all,fmt='o',capsize =2)
    #ax1.errorbar(x,yall,err,fmt='o',capsize =2)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"Boost factor" , fontsize = 30)
    #ax1.set_ylim([-2,2])
    ax1.set_xlim([0.28 , 305])
    ax1.tick_params(axis = "y", labelsize = 20)
    
    ############ L2 ##########################
    x = L2[0,:]
    
    yall = L2[1,:] - rand_S2[1,:]
    wall = rand_S2[3,:]/L2[3,:]
    yall = yall*wall
    
    yh = snr_L2[1,:] - snr_rand_S2[1,:]
    wh = snr_rand_S2[3,:]/snr_L2[3,:]    
    mask_z_lens = (lens_z>0.3)&(lens_z<0.5)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yh = yh * ratio / wh
    
    
    yl = snr_L2[5,:] - snr_rand_S2[5,:]
    wl = snr_rand_S2[7,:]/snr_L2[7,:]    
    mask_z_lens = (lens_z>0.3)&(lens_z<0.5)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yl = yl * ratio / wl
    
    err_all = np.std(jk_L2[:,1,:] - jk_L2[:,5,:], axis =0)
    var_ha = (yh/yall)**2. * (err_all**2/yh**2 + err_all**2/yall**2)
    var_la = (yl/yall)**2. * (err_all**2/yl**2 + err_all**2/yall**2)
    var_tot = var_ha + var_la
    err = var_tot**0.5
    
    print yh
    ax2.errorbar(x,yh,err_all,fmt='o',capsize =2)
    ax2.errorbar(x,yl,err_all,fmt='o',capsize =2)
    #ax1.errorbar(x,yall,err,fmt='o',capsize =2)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.legend(fontsize = 15)
    ax2.set_yticks([])
    #ax2.set_ylim([-2,2])
    ax2.set_xlim([0.28 , 305])
    #ax2.tick_params(axis = "y", labelsize = 20)
    
    ############ L3 ##########################
    x = L3[0,:]
    
    yall = L3[1,:] - rand_S3[1,:]
    wall = rand_S3[3,:]/L3[3,:]
    yall = yall*wall
    
    yh = snr_L3[1,:] - snr_rand_S3[1,:]
    wh = snr_rand_S3[3,:]/snr_L3[3,:]
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yh = yh * ratio / wh
        
    yl = snr_L3[5,:] - snr_rand_S3[5,:]
    wl = snr_rand_S3[7,:]/snr_L3[7,:]
    mask_z_lens = (lens_z>0.5)&(lens_z<0.7)
    ratio = 1. * len(random_ra)/len(lens_z[mask_z_lens])
    yl = yl * ratio / wl
        
    err_all = np.std(jk_L3[:,1,:] - jk_L3[:,5,:], axis =0)
    var_ha = (yh/yall)**2. * (err_all**2/yh**2 + err_all**2/yall**2)
    var_la = (yl/yall)**2. * (err_all**2/yl**2 + err_all**2/yall**2)
    var_tot = var_ha + var_la
    err = var_tot**0.5
   
    ax3.errorbar(x,yh,err_all,fmt='o',capsize =2)
    ax3.errorbar(x,yl,err_all,fmt='o',capsize =2)
    #ax1.errorbar(x,yall,err,fmt='o',capsize =2)
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.legend(fontsize = 15)
    #ax3.set_ylim([-2,2])
    ax3.set_xlim([0.28 , 305])
    ax3.set_yticks([])
    #ax3.tick_params(axis = "y", labelsize = 20)
    #ax1.axhline(1, color="k", linestyle = "dashed")
    #ax2.axhline(1, color="k", linestyle = "dashed")
    #ax3.axhline(1, color="k", linestyle = "dashed")

    ax1.tick_params(axis = "x", labelsize = 20)
    ax2.tick_params(axis = "x", labelsize = 20)
    ax3.tick_params(axis = "x", labelsize = 20)
    
    ax1.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax2.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)
    ax3.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 20)

    ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    ax2.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    ax3.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None

def xshear_psf(zl1,zl2,zs1,zs2,lens_type):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["psf1"][:]
    source_e2 = source_file["psf2"][:]
    source_w = source_file["w"][:]
    source_file.close()
    
    #print "diagnosis" , source_z.min(), source_z.max()
    
    if lens_type == 'lum':
        lens_file = h5py.File("LRG_lum.h5")
    
    elif lens_type == 'dense':
        lens_file = h5py.File("LRG_dense.h5")
    
    lens_z = lens_file["redshift"][:]
    lens_ra = lens_file["RA"][:]
    lens_dec = lens_file["DEC"][:]
    lens_file.close()
    
    mask_z_lens = (lens_z>zl1)&(lens_z<zl2)
    
    lens_cat = treecorr.Catalog(x=lens_ra[mask_z_lens], y=lens_dec[mask_z_lens], x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,g1=source_e1,g2=-1.*source_e2,w=source_w,x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x , weight = ng.meanr , ng.xi, ng.xi_im, ng.weight
    
    return r, xi_t, xi_x, weight

def xshear_psf_random(zs1,zs2):
    
    source_file = h5py.File("source_zb_"+str(zs1)+"_"+str(zs2)+".h5")
    source_ra = source_file["ra"][:]
    source_dec = source_file["dec"][:]
    source_z = source_file["zb"][:]
    source_e1 = source_file["psf1"][:]
    source_e2 = source_file["psf2"][:]
    source_w = source_file["w"][:]
    source_file.close()
    
    #print "diagnosis" , source_z.min(), source_z.max()
    
    random_file = h5py.File("randoms.h5")
    random_coord = random_file["random"][:] 
    random_ra = random_coord[:,0]
    random_dec = random_coord[:,1]
    
    shuffled = np.random.choice(len(random_ra) , 2000000)
    
    lens_cat = treecorr.Catalog(x=random_ra[shuffled], y=random_dec[shuffled], x_units='degree', y_units='degree')
    source_cat=treecorr.Catalog(x=source_ra,y=source_dec,g1=source_e1,g2=-1.*source_e2,w=source_w,x_units='degree',y_units='degree')
    
    ng = treecorr.NGCorrelation(nbins = 25, min_sep=0.3, max_sep=300, sep_units='arcmin', verbose=1)
    ng.process(lens_cat, source_cat)
    
    r, xi_t, xi_x,  weight = ng.meanr , ng.xi, ng.xi_im , ng.weight
    
    return r, xi_t, xi_x, weight

def xshear_psf_all():
    
    xshear_file = h5py.File("xpsf_lens_henk.h5" , "w")
    nbins = 25
    
    ######### lum ###########
    
    r , xt, xx , xw = xshear_psf(0.1,0.3,0.4,0.9,"dense")
    xshear_file.create_dataset("lum_L1", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
   
    r , xt, xx , xw  = xshear_psf(0.3,0.5,0.6,0.9,"dense")
    xshear_file.create_dataset("lum_L2", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    r , xt, xx , xw = xshear_psf(0.5,0.7,0.8,0.9,"dense")
    xshear_file.create_dataset("lum_L3", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    xshear_file.close()
    
    return None


def xshear_psf_random_all():
    
    nbins = 25
    xshear_file = h5py.File("xpsf_random_henk.h5" , "w")

    ######### lum ###########
    
    r , xt, xx, xw  = xshear_psf_random(0.4,0.9)
    xshear_file.create_dataset("rand_S1", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    r , xt, xx, xw  = xshear_psf_random(0.6,0.9)
    xshear_file.create_dataset("rand_S2", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    r , xt, xx, xw  = xshear_psf_random(0.8,0.9)
    xshear_file.create_dataset("rand_S3", (4, nbins) , data = np.vstack([r,xt,xx, xw ]))
    
    xshear_file.close()

def plot_xshear_psf_lum():
    
    rand_shear = h5py.File("xpsf_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xpsf_lens_henk.h5")
    
    L1 = xshear_file["lum_L1"][:]

    L2 = xshear_file["lum_L2"][:]
  
    L3= xshear_file["lum_L3"][:]
    
    
    jk_file = h5py.File("psf_jk.h5")
    
    psf1 = jk_file["lum_L1"][:]
    psf2 = jk_file["lum_L2"][:]
    psf3 = jk_file["lum_L3"][:]

    psf1 = psf1[:,1,:] - psf1[:,3,:]
    psf2 = psf2[:,1,:] - psf2[:,3,:]
    psf3 = psf3[:,1,:] - psf3[:,3,:]

    
    
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(1, 1)

    ax1 = plt.subplot(gs[0,0])
   
    ############# L1 #################
    
    
    err = np.std(psf1, axis = 0)
    
    factor = (100-25-1)/(100.-1.)
    C1 = np.cov(psf1.T)
    print C1.shape
    C1inv = factor * np.linalg.inv(C1)
    chisq = factor * np.dot(L1[1,:]  - rand_S1[1,:], np.dot(C1inv, L1[1,:]  - rand_S1[1,:]))
    #chisq = np.sum((L1[1,:]  - rand_S1[1,:])**2/err**2)
    print "chisq", chisq
    ax1.errorbar(L1[0,:] , L1[1,:], err, fmt='o', capsize = 2 , label = r"$0.1<z_l<0.3\;, \; 0.4<z_B<0.9$")
    ax1.set_xscale("log")
    #ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"$\gamma_{T} \; \; (\mathrm{PSF \; \; leakage})$" , fontsize = 30)
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])
    #ax1.set_xticks([0.1,1,10,100,300])
    
    err = np.std(psf2, axis = 0)
    chisq = np.sum((L2[1,:]  - rand_S2[1,:])**2/err**2)
    print "chisq", chisq

    ax1.errorbar(L2[0,:] , L2[1,:], err, fmt='o', capsize = 2, label = r"$0.3<z_l<0.5\;, \; 0.6<z_B<0.9$")
    ax1.set_xscale("log")
    #ax2.set_yscale("log")     
    ax1.legend(fontsize = 15)

    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    
    err = np.std(psf3, axis = 0)
    chisq = np.sum((L3[1,10:19]  - rand_S3[1,10:19])**2/err[10:19]**2)
    print "chisq", chisq
    ax1.errorbar(L3[0,:] , L3[1,:] , err, fmt='o', capsize = 2 ,label = r"$0.5<z_l<0.7\; , \; 0.8<z_B<0.9$")
    ax1.set_xscale("log")
    plt.axhline(y = 0.0 , color = "k", linestyle = 'dashed')
    #ax3.set_yscale("log")    
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    ax1.legend(fontsize = 15)
    
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)

    ax1.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 30)


    #ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    #ax1.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    #ax1.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None




###############################################################

def plot_cross_lum():
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["lum_L1"][:]

    L2 = xshear_file["lum_L2"][:]
  
    L3= xshear_file["lum_L3"][:]
    
    
    jk_file = h5py.File("lum_jk.h5")
    
    psf1 = jk_file["lum_L1"][:]
    psf2 = jk_file["lum_L2"][:]
    psf3 = jk_file["lum_L3"][:]

    psf1 = psf1[:,2,:] - psf1[:,6,:]
    psf2 = psf2[:,2,:] - psf2[:,6,:]
    psf3 = psf3[:,2,:] - psf3[:,6,:]

    
    
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(1, 1)

    ax1 = plt.subplot(gs[0,0])
   
    ############# L1 #################
    
    
    err = np.std(psf1, axis = 0)
    
    factor = (100-25-1)/(100.-1.)
    C1 = np.cov(psf1.T)
    print C1.shape
    C1inv = factor * np.linalg.inv(C1)
    chisq = factor * np.dot(L1[2,:]  - rand_S1[2,:], np.dot(C1inv, L1[2,:]  - rand_S1[2,:]))
    #chisq = np.sum((L1[1,:]  - rand_S1[1,:])**2/err**2)
    print "chisq", chisq
    ax1.errorbar(L1[0,:] , L1[2,:]-rand_S1[2,:] , err, fmt='o', capsize = 2 , label = r"$0.1<z_l<0.3\;, \; 0.4<z_B<0.9$")
    ax1.set_xscale("log")
    #ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"$\gamma_{\times}$" , fontsize = 30)
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])
    #ax1.set_xticks([0.1,1,10,100,300])
    
    err = np.std(psf2, axis = 0)
    chisq = np.sum((L2[1,:]  - rand_S2[1,:])**2/err**2)
    print "chisq", chisq

    ax1.errorbar(L2[0,:] , L2[2,:]-rand_S2[1,:] , err, fmt='o', capsize = 2, label = r"$0.3<z_l<0.5\;, \; 0.6<z_B<0.9$")
    ax1.set_xscale("log")
    #ax2.set_yscale("log")     
    ax1.legend(fontsize = 15)

    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    
    err = np.std(psf3, axis = 0)
    chisq = np.sum((L3[2,10:]  - rand_S3[2,10:])**2/err[10:]**2)
    print "chisq", chisq
    ax1.errorbar(L3[0,:] , L3[2,:]-rand_S3[2,:] , err, fmt='o', capsize = 2 ,label = r"$0.5<z_l<0.7\; , \; 0.8<z_B<0.9$")
    ax1.set_xscale("log")
    plt.axhline(y = 0.0 , color = "k", linestyle = 'dashed')
    #ax3.set_yscale("log")    
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    ax1.legend(fontsize = 15)
    
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)

    ax1.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 30)


    #ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    #ax1.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    #ax1.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None

########################################################################################################


from astropy.stats import median_absolute_deviation

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def plot_random_lum():
    
    rand_shear = h5py.File("xshear_random_henk.h5")
    
    rand_S1 = rand_shear["rand_S1"][:]
    rand_S2 = rand_shear["rand_S2"][:]
    rand_S3 = rand_shear["rand_S3"][:]
    
    xshear_file = h5py.File("xshear_lens_henk.h5")
    
    L1 = xshear_file["lum_L1"][:]

    L2 = xshear_file["lum_L2"][:]
  
    L3= xshear_file["lum_L3"][:]
    
    
    jk_file = h5py.File("lum_jk.h5")
    
    psf1 = jk_file["lum_L1"][:]
    psf2 = jk_file["lum_L2"][:]
    psf3 = jk_file["lum_L3"][:]

    psf1 = psf1[:,5,:]
    psf2 = psf2[:,5,:] 
    psf3 = psf3[:,5,:] 

    
    
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(1, 1)

    ax1 = plt.subplot(gs[0,0])
   
    ############# L1 #################
    
    
    err = np.std(psf1, axis = 0)
    
    factor = (100-25-1)/(100.-1.)
    C1 = np.cov(psf1.T)
    print C1.shape
    C1inv = factor * np.linalg.inv(C1)
    chisq = factor * np.dot(rand_S1[1,:], np.dot(C1inv, rand_S1[1,:]))
    #chisq = np.sum((L1[1,:]  - rand_S1[1,:])**2/err**2)
    print "chisq", chisq
    ax1.errorbar(L1[0,:] , rand_S1[1,:] , err, fmt='o', capsize = 2 , label = r"$0.1<z_l<0.3\;, \; 0.4<z_B<0.9$")
    ax1.set_xscale("log")
    #ax1.set_yscale("log")
    
    ax1.legend(fontsize = 15)
    ax1.set_ylabel(r"$\gamma_{T} \; \; (\mathrm{random})$" , fontsize = 30)
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])
    #ax1.set_xticks([0.1,1,10,100,300])
    
    err = np.std(psf2, axis = 0)
    chisq = np.sum((rand_S2[1,:])**2/err**2)
    print "chisq", chisq

    ax1.errorbar(L2[0,:] , rand_S2[1,:] , err, fmt='o', capsize = 2, label = r"$0.3<z_l<0.5\;, \; 0.6<z_B<0.9$")
    ax1.set_xscale("log")
    #ax2.set_yscale("log")     
    ax1.legend(fontsize = 15)

    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    
    err = np.std(psf3, axis = 0)
    chisq = np.sum((rand_S3[1,10:20])**2/err[10:20]**2)
    print "chisq", chisq
    ax1.errorbar(L3[0,:] ,rand_S3[1,:] , err, fmt='o', capsize = 2 ,label = r"$0.5<z_l<0.7\; , \; 0.8<z_B<0.9$")
    ax1.set_xscale("log")
    plt.axhline(y = 0.0 , color = "k", linestyle = 'dashed')
    #ax3.set_yscale("log")    
    ax1.set_ylim([-5*10**-4. , 5*10**-4.])
    ax1.set_xlim([0.28 , 305])    
    #ax1.set_yticklabels([])
    #ax1.set_xticklabels([])
    ax1.legend(fontsize = 15)
    
    plt.xticks(fontsize = 20)
    plt.yticks(fontsize = 20)

    ax1.set_xlabel(r"$\mathrm{separation \; [arcmin]}$" , fontsize = 30)


    #ax1.set_title(r"$0.1<z_l<0.3$" , fontsize = 20)
    #ax1.set_title(r"$0.3<z_l<0.5$" , fontsize = 20)
    #ax1.set_title(r"$0.5<z_l<0.7$" , fontsize = 20)    
    
    
    fig.tight_layout()
    
    return None




if __name__ == '__main__':

   xshear_flip()
   #xshear_all()
   #xshear_random_all()
   #plot_xshear_dense_p()
   #plot_xshear_lum_p()
   #plot_xshear_lum()
   #plot_boost()
   #xshear_snr_all()
   #xshear_snr_random_all()
   #xshear_size_all()
   #xshear_size_random_all()
   #plot_snr()
   #xshear_psf_all()
   #xshear_psf_random_all()
   #plot_xshear_psf_lum()
   #plot_cross_lum()
   #plot_random_lum()
