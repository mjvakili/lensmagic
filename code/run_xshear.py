import os
import sys
import stat

def bash_writer(i):

    bash_file = 'sh.xshear_runner'+str(i)
    bash_f = open(bash_file, 'w')
    bash_file_content = '\n'.join([
                        "#! /bin/csh -f",
  			"set echo",
			"set pp = /net/vuntus/data2/vakili/xshear/src",  #the directory where the executable lives
  			"$cat ../../flagship/code/source_file_zbin"+str(i)+" | xshear pointz.cfg ../../flagship/code/lens_file_zbin"+str(i)+" > output_file_zbin"+str(i),
	                 "exit",
                        ])
    bash_f.write(bash_file_content)
    bash_f.close()
    st = os.stat(bash_file)
    os.chmod(bash_file, st.st_mode | stat.S_IEXEC) 
    
    return None

def xshear_executer(i):
    '''
    Execute xshear bash script
    '''
    bash_writer(i)
    bash_file = 'sh.xshear_runner'+str(i)
    bash_cmd = './'+bash_file
    os.system(bash_cmd)

    return bash_cmd    

if __name__ == '__main__':

    import multiprocessing
    from multiprocessing import Pool
    Nthreads = 2
    pool = Pool(Nthreads)
    mapfn = pool.map
    arglist = [None] * Nthreads
    for i in range(Nthreads):
       arglist[i] = i
    result = list(mapfn(xshear_executer, [ars for ars in arglist]))   
    pool.close()
