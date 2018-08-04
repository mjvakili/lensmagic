import pyfits as pf
import h5py
import numpy as np
import kmeans_radec
import matplotlib.pyplot as plt
plt.switch_backend("Agg")
import treecorr
import seaborn as sns

shear_data = ["../../redsequence/data/KiDS_DR3.1_G9_ugri_shear.fits", 
 "../../redsequence/data/KiDS_DR3.1_G12_ugri_shear.fits",
 "../../redsequence/data/KiDS_DR3.1_G15_ugri_shear.fits",
 "../../redsequence/data/KiDS_DR3.1_G23_ugri_shear.fits",
 "../../redsequence/data/KiDS_DR3.1_GS_ugri_shear.fits"]


kids_cat = pf.open("../../redsequence/data/KiDS_DR3.1_G9_ugri_shear.fits")
e1 = kids_cat[1].data['e1']
e2 = kids_cat[1].data['e2']
w = kids_cat[1].data['weight']
z_b = kids_cat[1].data['Z_B']
sg = kids_cat[1].data['SG_FLAG']
mask = kids_cat[1].data['MASK']
fitclass = kids_cat[1].data['fitclass']
ra = kids_cat[1].data['RAJ2000']
dec = kids_cat[1].data['DECJ2000']
snr = kids_cat[1].data['model_SNratio']
rm = kids_cat[1].data['bias_corrected_scalelength']
psf1 = kids_cat[1].data['PSF_e1']
psf2 = kids_cat[1].data['PSF_e2']
magr = kids_cat[1].data['MAG_r']
bias = kids_cat[1].data['m']
mask = kids_cat[1].data['MASK']
flag = kids_cat[1].data['Flag']

#sns.distplot(kids_cat[1].data['MASK'], kde= False)
#plt.xscale("log")
#plt.savefig("/home/vakili/public_html/mask.png")
#plt.close()
#sns.distplot(kids_cat[1].data['Flag'], kde= False)
#plt.savefig("/home/vakili/public_html/flag.png")
#plt.close()

star_mask = np.where((sg==1)&(magr>20)&(fitclass==0)&(rm>0.5)&(mask==0)&(flag==0))[0]

e1, e2, w, z_b, ra, dec , snr, rm , psf1, psf2, magr, bias = e1[star_mask], e2[star_mask], w[star_mask], z_b[star_mask], ra[star_mask], dec[star_mask], snr[star_mask],rm[star_mask], psf1[star_mask], psf2[star_mask], magr[star_mask], bias[star_mask]

for i in range(1,5): 
    kids_cat = pf.open(shear_data[i])

    e1c = kids_cat[1].data['e1']
    e2c = kids_cat[1].data['e2']
    wc= kids_cat[1].data['weight']
    z_bc = kids_cat[1].data['Z_B']
    sgc = kids_cat[1].data['SG_FLAG']
    maskc = kids_cat[1].data['MASK']
    fitclassc = kids_cat[1].data['fitclass']
    rac = kids_cat[1].data['RAJ2000']
    decc = kids_cat[1].data['DECJ2000']
    snrc = kids_cat[1].data['model_SNratio']
    rmc = kids_cat[1].data['bias_corrected_scalelength']
    psf1c = kids_cat[1].data['PSF_e1']
    psf2c = kids_cat[1].data['PSF_e2']
    magrc = kids_cat[1].data['MAG_r']
    biasc = kids_cat[1].data['m']
    maskc = kids_cat[1].data['MASK']
    flagc = kids_cat[1].data['Flag']

    star_maskc = np.where((sgc==1)&(magrc>20)&(fitclassc==0)&(maskc==0)&(flagc==0))[0]
    e1c, e2c, wc, z_bc, rac, decc, snrc, rmc, psf1c, psf2c , magrc, biasc = e1c[star_maskc], e2c[star_maskc], wc[star_maskc], z_bc[star_maskc], rac[star_maskc], decc[star_maskc], snrc[star_maskc], rmc[star_maskc], psf1c[star_maskc], psf2c[star_maskc], magrc[star_maskc], biasc[star_maskc]

    e1 = np.hstack([e1,e1c])  
    e2 = np.hstack([e2,e2c])        
    w = np.hstack([w,wc])        
    z_b = np.hstack([z_b,z_bc])        
    ra = np.hstack([ra,rac])        
    dec = np.hstack([dec,decc]) 
    snr = np.hstack([snr,snrc])
    rm = np.hstack([rm,rmc])
    psf1 = np.hstack([psf1,psf1c])
    psf2 = np.hstack([psf2,psf2c])
    bias = np.hstack([bias,biasc])
    magr = np.hstack([magr,magrc])

kids_cat[1].header

for i in range(4):
    print i
    result_file = h5py.File("source_zb_"+str(0.1+i*(0.2))+"_"+str(0.1+(i+1)*(0.2))+".h5" , 'w')
    zmask = (z_b>0.1+i*(0.2))&((z_b<0.1+(i+1)*(0.2))|(z_b==0.1+(i+1)*(0.2)))
    ns = len(ra[zmask])
    result_file.create_dataset("ra", (ns, ) , data = ra[zmask])
    result_file.create_dataset("dec",(ns, ) , data = dec[zmask])
    result_file.create_dataset("e1", (ns, ) , data = e1[zmask])
    result_file.create_dataset("e2", (ns, ) , data = e2[zmask])
    result_file.create_dataset("zb", (ns, ) , data = z_b[zmask])
    result_file.create_dataset("w", (ns, ) , data = w[zmask])
    result_file.create_dataset("rad", (ns,) , data= rm[zmask])
    result_file.create_dataset("snr", (ns ,) , data = snr[zmask])
    result_file.create_dataset("bias", (ns ,) , data = bias[zmask])

    result_file.close()


for i in range(3):
    print i
    result_file = h5py.File("source_zb_"+str(0.4+i*(0.2))+"_0.9.h5" , 'w')
    
    bad_north = ((dec>-10)&(((ra>145)&(ra<171))|((ra>195)&(ra<210))|(ra>227)))
    bad_south = (dec<-28)&(ra<30)
    bad = (bad_south)|(bad_north)
    
    
    zmask = (z_b>0.4+i*(0.2))&((z_b<0.9)|(z_b==0.9))
    zmask = zmask&(~bad)
    ns = len(ra[zmask])
    result_file.create_dataset("ra", (ns, ) , data = ra[zmask])
    result_file.create_dataset("dec",(ns, ) , data = dec[zmask])
    result_file.create_dataset("e1", (ns, ) , data = e1[zmask])
    result_file.create_dataset("e2", (ns, ) , data = e2[zmask])
    result_file.create_dataset("zb", (ns, ) , data = z_b[zmask])
    result_file.create_dataset("w", (ns, ) , data = w[zmask])
    result_file.create_dataset("rad", (ns,) , data= rm[zmask])
    result_file.create_dataset("snr", (ns ,) , data = snr[zmask])
    result_file.create_dataset("psf1", (ns,) , data= psf1[zmask])
    result_file.create_dataset("psf2", (ns ,) , data = psf2[zmask])
    result_file.create_dataset("bias", (ns ,) , data = bias[zmask])
    result_file.close()


