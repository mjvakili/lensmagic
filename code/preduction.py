import pyfits as pf
import numpy as np
import astropy.io.fits as fits

filename = '../data/source64.fits'
cat = fits.open(filename,memmap=True)

data = cat[1].data
header = cat[1].header
	
zbins = np.linspace(0.9,2.0,11)

def source(iz):
    zmin , zmax = zbins[iz] , zbins[iz+1]
    print zmin , zmax
    mask = data['true_redshift_gal']>zmin
    return data[mask]


def write(iz):

    source_gals = source(iz)
    shape = source_gals.shape[0]
    print "shape" , shape
    np.savetxt("/net/vuntus/data2/vakili/flagship/code/source64_file_zbin"+str(iz) , \
    zip(source_gals['ra_gal'] , source_gals['dec_gal'] , source_gals['gamma1'] , -1.*source_gals['gamma2'], np.ones((shape)) , source_gals['true_redshift_gal']), \
    fmt='%1.8f\t%1.8f\t%1.8f\t%1.8f\t%i\t%1.8f')
    
    print "finished writing into the source file" , iz

    return None

if __name__ == '__main__':

    import multiprocessing
    from multiprocessing import Pool
    Nthreads = 10
    pool = Pool(Nthreads)
    mapfn = pool.map
    arglist = [None] * Nthreads
    for i in range(Nthreads):
       arglist[i] = i
    result = list(mapfn(write, [ars for ars in arglist]))   
    pool.close()

