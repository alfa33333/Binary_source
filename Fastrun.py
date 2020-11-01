import os
import sys
#import mkl
#mkl.set_num_threads(1)
os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1" 
########################################
# comment this block to read your data in your own way
########################################
import Read_Local_re as read
import Read_kmt as readkmt
import Read_Local_OGLE as readOGLE
#######################################

from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt
import SingleLensFitter as slf
from tlib import *
from time import time
#


if __name__ == '__main__':
    #####################
    #This block only reads data. You can include your own as long as
    #it is a dictionary of each source containing a tupple of (date,flux,err_flux)
    #####################
    trange = (8200,8260)
    
    # renormdic = readdict("./renorm/renormalise_values_test.dat")
    # for i in renormdic.keys():
    #     print('source {} : [ {} {} ]'.format(i,renormdic[i][0],renormdic[i][1]))
    # data = read.read_data('./Data/',trange,renorm=renormdic,minmag=False)
    data = read.read_data('./Data/',trange) 
    #data = readkmt.read_data('./Data/',trange)
    #rmdpoint(data,'KMTS37-I',8198.50797)
    # removedict = removesignal('./Data/filter', trange)
    # print('Filtering data ...')
    # for y in data.keys():
    #     rmdpoint(data,y, removedict[y])
    # dataOGLE = readOGLE.read_data('./Data/blg-0677/',trange,renorm=renormdic,minmag=False)
    # data[list(dataOGLE.keys())[0]] = dataOGLE[list(dataOGLE.keys())[0]]
    # dataOGLE = readOGLE.read_data('./Data/blg-0680/',trange,renorm=renormdic,minmag=False)
    # data[list(dataOGLE.keys())[0]] = dataOGLE[list(dataOGLE.keys())[0]]
    

    #######################################
    fitter = slf.SingleLensFitter(data,[1.07056800e-01, 8.22955366e+03, 5.00686788e+00],reference_source='CTIOI01')  


    folder = './output/Finalbinary0677D-renom2/'
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_prefix = 'binary-source-dynesty'
    fitter.add_binary_source(system='both')

    #Nested parameters for the  dnesty option
    fitter.nlive = 500
    fitter.tol = 0.1
    fitter.dynestyparallel = True
    fitter.sampletype = 'rslice'
    fitter.cpu = 10
    
    #emcee parameters for the Metropolist option
    fitter.nwalkers = 100
    fitter.nsteps = 50
    fitter.nsteps_production = 200 
    fitter.gen_optimize_first = True
    fitter.emcee_optimize_first = True 

    #limits 
    fitter.t0t1extlim = True
    fitter.t0_limits = (8220, 8240)
    fitter.t02_limits = (8228, 8232)
    
    #plotting
    fitter.plotrange = (8228, 8231)
    
    #fitter.plotprefix = folder+file_prefix
    fitter.p = [0.10304,8229.5414,4.93678,0.00019,8229.86480,0.00021] # irrelevant if using dnesty
    #Comments the following 4 lines if you want to use emcee only.
    t0 = time()
    fitter.dnesty()
    t1 = time()
    print('Sampler time: {:f}'.format(t1-t0))

    #It run emcee for a short convergence
    fitter.plotprefix = folder+'binext2-lens-emcee'
    t0 = time()
    fitter.fitparallel(cpu_cores=10, max_opt_iterations=1)
    t1 = time()
    
    print('Sampler time: {:f}'.format(t1-t0))
    # fitter.pmin2 = fitter.p
    chi2, _ = fitter.chi2_calc(p=fitter.pmin2)
    np.save(fitter.plotprefix+'p_minchi2',fitter.pmin2)
    with open(folder+file_prefix+'minchi2.dat','w') as fid:
        fid.write('The minchi2 is {:f} \n'.format(chi2))
    fid.close()
    print('The minchi2 is {:f} \n'.format(chi2))
