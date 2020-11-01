import sys
import os
import numpy as np


def rmdpoint(data,site,timepoint):
	''' Help method to remove manually outliers in the loaded data.

	Parameters :
	------------
	data = dictionary 
	 Full dictionary of data.
	site = str
	 desire site to for removal.
	timepoint = int, float or list of int or floats
	 exact point(s) in time to remove. 


	Returns:
	-----------
	out:
	In place removal of the data library site selected.	

	'''
	t = data[site][0]
	y = data[site][1]
	dy = data[site][2]
	npoints = 1
	ltimepoint = [timepoint]  
	# Checking if is an array or point.
	if isinstance(timepoint, list):
		npoints = len(timepoint)
		ltimepoint = timepoint
		print('list points ...')
	for i in range(npoints):
		outl = np.where(t == ltimepoint[i])[0]
		t = np.delete(t, outl)
		y = np.delete(y, outl)
		dy = np.delete(dy, outl)

	data[site] = ( t,y,dy )

###########
def removeq(data_directory, t_range,threshold = 1.5):
	''' Help function to remove datapoint above a threshold of q value.

	Parameters:
	----------
	data_directory: str
		Address to the directory where the data is located.
	t_range : tupla or list 
		Initial and final values in which the data will be contained.
	threshold: float
		Discriminator for q value.

	Returns:
	--------
	remove: Dictionary
		Dictionary for the datafiles with the values above the the 
		threshold. 
	
	'''	

	remove = {}
	files = [f for f in os.listdir(data_directory)]

	for f in files:
		source = f.split('.')[0]

		datafile = data_directory+'/'+f

		if os.path.exists(datafile):
			d = np.loadtxt(datafile,comments=['<','#'],usecols=(0,5))
			t = d[:,0]
			if t[0] > 2440000:
				t -= 2450000
			points = np.where((t_range[0] < t) & (t< t_range[1]) )[0]
			q = np.where(d[points,1] > threshold)
			tq = t[points[q]]
			tsave = tq.tolist()
			remove[source] = tsave
	return remove

##########
def removesignal(data_directory, t_range,threshold = 200.0):
	''' Help function to remove datapoint below a threshold of signal value.

	Parameters:
	----------
	data_directory: str
		Address to the directory where the data is located.
	t_range : tupla or list 
		Initial and final values in which the data will be contained.
	threshold: float
		Discriminator for q value.

	Returns:
	--------
	remove: Dictionary
		Dictionary for the datafiles with the values below the the 
		threshold. 
	
	'''	

	remove = {}
	files = [f for f in os.listdir(data_directory)]

	for f in files:
		source = f.split('.')[0]

		datafile = data_directory+'/'+f

		if os.path.exists(datafile):
			d = np.loadtxt(datafile,comments=['<','#'],usecols=(0,9))
			t = d[:,0]
			if t[0] > 2440000:
				t -= 2450000
			points = np.where((t_range[0] < t) & (t< t_range[1]) )[0]
			q = np.where(d[points,1] < threshold)
			tq = t[points[q]]
			tsave = tq.tolist()
			remove[source] = tsave
	return remove

###########
def writelist(name,listval):
	''' Custom method to write values to a file.

	Parameters:
	----------
	name: str
		Address to the directory and name where the file is going to be writen.
	listval : list 
		list of values.

	Returns:
	--------
	
	'''
	f = open(name,"w+")
	for i in range(len(listval)):
		f.write(listval[i]+'\n')
	f.close()

def writedict(fname,dict):
	''' Custom method to write values of a dictionary to a file.

	Parameters:
	----------
	fname: str
		file name handle. It can include the full path.
	dict : dictionary 
		dictionary to be written, it only supports single elements (number or strings), lists of numbers and numpy arrays.

	Returns:
	--------
	
	'''
	f = open(fname,"w+")
	for i in dict.keys():
		f.write(str(i) + '\t')
		if isinstance(i,str):
			f.write('str'+ '\t')
		else:
			f.write('num'+ '\t')
		if isinstance(dict[i],list):
			for j in dict[i]:
				f.write('{}'.format(j)+'\t')
		elif isinstance(dict[i],np.ndarray):
			for j in dict[i]:
				f.write('{}'.format(j)+'\t')
		else:
			f.write('{}'.format(dict[i])+'\t')
		f.write('\n')
	f.close()

def readdict(fname):
	''' Custom method to read values saved by writedict(). If the original dict had numpy objects these become lists.

	Parameters:
	----------
	fname: str
		file name handle. It can include the full path.

	Returns:
	--------
	dict: dictionary
		dictionary from file.
	
	'''
	retdict = {}
	with open(fname,"r") as f:
		for line in f:
			data = line.split('\t')
			if data[1] == 'str':
				key = data[0]
			elif data[1] == 'num':
				try:
					key = int(data[0])
				except:
					key = float(data[0])
			else:
				raise Exception('Keys have unknown format.')
			if len(data) > 4 :
				temp = data[2:len(data)-1]
				try:
					retdict[key] = [int(i) for i in temp]
				except:
					try:
						retdict[key] = [float(i) for i in temp]
					except:
						raise Exception('List values are not numbers.')
			elif data[2] == '\n':
				retdict[key] = []
			else:
				try:
					retdict[key] = int(data[2]) 
				except:
					try:
						retdict[key] = float(data[2])
					except:
						if isinstance(data[2],str):
							retdict[key] = data[2]
						else:
							raise Exception('Single value is a recognized type.')
				
										

	return retdict

############3
def readlist(name):
	''' Custom function to read a list created with writelist method.

	Parameters:
	----------
	name: str
		Address to the directory and name where the file is located.
	Returns:
	--------
	val: list
		List of strings of the list of values on the file.
	
	'''
	with open(name,'r') as f:
		val = f.read()
	val = val.split()
	f.close()
	return val

#####################
def renormalise_data_uncertainties_try1(self,p=None,source=None):
		import scipy
		from scipy.optimize import least_squares

		"""Adjust data uncertainties to force reduced chi^2 to 1 for each data source."""

		print('Renormalising data uncertainties')
		savedict = {}

		if p is None:
			p = self.p

		for site in self.data:

			if source is None or site == source:
				magnitude = self.zp - 2.5*np.log10(self.data[site][1])
				err_mag = self.zp - 2.5*np.log10(self.data[site][1]+self.data[site][2]) - magnitude

				flux = 10**(0.4*(self.zp-magnitude))
				err_flux = 10**(0.4*(self.zp - magnitude - err_mag))- flux

				chi2, _, _, _, _, chi2_elements = self.chi2_calc(p,source=site)

				chi2_elements_site = chi2_elements[site]

				mag = self.magnification(self.data[site][0])
				points = np.argsort(mag)
				tpoints = self.data[site][0][points]
				print('tpoints=', len(tpoints))


				ls = p[5]-0.5
				rs = p[5]+0.5

				def func1(x,dmsq_high, dmsq_low, err_mag_high, err_mag_low, magnitude, flux_high_mag_points, flux_low_mag_points):
					scale0 = x[0]
					eps = x[1]

					err_mag_high_new = scale0*np.sqrt(err_mag_high**2 + eps**2)
					err_mag_low_new = scale0*np.sqrt(err_mag_low**2 + eps**2)

					err_flux_high_new = 10.0**(0.4*(self.zp - magnitude[high_mag_points] - err_mag_high_new)) - flux_high_mag_points
					err_flux_low_new = 10.0**(0.4*(self.zp - magnitude[low_mag_points] - err_mag_low_new)) - flux_low_mag_points	
				
					diff1 = (np.sum(dmsq_high/(err_flux_high_new**2)) / (len(high_mag_points)))-1.0

					diff2 = (np.sum(dmsq_low/(err_flux_low_new**2)) / (len(low_mag_points)))-1.0

					return diff1, diff2


				#high_mag_points = np.where((ls<=tpoints) & (tpoints>=rs))[0]

				high_range = np.linspace(0.1,0.5,10)

				for i in range(len(high_range)):
					high_range_point = high_range[i]
					high_mag_points = np.where(mag>=np.max(mag)*high_range_point)[0]
					if len(high_mag_points)>=self.dims:
						high_mag_points = high_mag_points

				low_mag_points = np.setdiff1d(points,high_mag_points)

				chi2_elements_high = chi2_elements_site[high_mag_points]
				flux_high_mag_points = flux[high_mag_points]
				err_mag_high = err_mag[high_mag_points]
				err_flux_high = 10**(0.4*(self.zp - magnitude[high_mag_points] - err_mag_high))- flux_high_mag_points
				dmsq_high = chi2_elements_high*(err_flux_high**2)


				chi2_elements_low = chi2_elements_site[low_mag_points]
				err_mag_low = err_mag[low_mag_points]
				flux_low_mag_points = flux[low_mag_points]
				err_flux_low = 10**(0.4*(self.zp - magnitude[low_mag_points] - err_mag_low))- flux_low_mag_points
				dmsq_low = chi2_elements_low*(err_flux_low**2)
        
				if self.dims>len(points):
					scale = np.sqrt(chi2/(len(points)-self.dims-2))
				else:
					scale = np.sqrt(chi2/(len(points)))

				eps = 0.0001
				x0 = (scale,eps)
				
				value = scipy.optimize.least_squares(func1,x0,bounds=([0.0,0.0], [8.0,1.0]),args=(dmsq_high, dmsq_low, err_mag_high, err_mag_low, magnitude,flux_high_mag_points, flux_low_mag_points),max_nfev=10000000)
				print('value=', value.x, site)
				scalev = (value.x[0])
				epsv = (value.x[1])
				savedict[site] = [scalev,epsv]
				writedict('renormalise_values.dat',savedict)

				sigma = scalev*np.sqrt(err_mag**2+epsv**2)
				y = 10.0**(0.4*(self.zp-magnitude))
				dy = 10.0**(0.4*(self.zp-magnitude+sigma))-y
				self.data[site] = (self.data[site][0],y,dy)
				self.reformat_data()
				chi2_new,_,_,_,_,_ = self.chi2_calc(p,source=site)
				print('chi2_new=', chi2_new)


				
def renormalise_data_uncertainties_single(self,p=None,source=None):
		"""Adjust data uncertainties to force reduced chi^2 to 1 for each data source. version for single lens"""
		import scipy
		from scipy.optimize import least_squares

		

		print('Renormalising data uncertainties')
		savedict ={}

		if p is None:
			p = self.p

		for site in self.data:

			if source is None or site == source:
				magnitude = self.zp - 2.5*np.log10(self.data[site][1])
				err_mag = self.zp - 2.5*np.log10(self.data[site][1]+self.data[site][2]) - magnitude

				flux = 10**(0.4*(self.zp-magnitude))
				err_flux = 10**(0.4*(self.zp - magnitude - err_mag))- flux

				chi2, chi2_elements = self.chi2_calc(p,source=site)

				chi2_elements_site = chi2_elements[site]

				mag = self.magnification(self.data[site][0])
				points = np.argsort(mag)
				tpoints = self.data[site][0][points]
				print('tpoints=', len(tpoints))


				def func1(x,dmsq_high, dmsq_low, err_mag_high, err_mag_low, magnitude, flux_high_mag_points, flux_low_mag_points):
					scale0 = x[0]
					eps = x[1]

					err_mag_high_new = scale0*np.sqrt(err_mag_high**2 + eps**2)
					err_mag_low_new = scale0*np.sqrt(err_mag_low**2 + eps**2)

					err_flux_high_new = 10.0**(0.4*(self.zp - magnitude[high_mag_points] - err_mag_high_new)) - flux_high_mag_points
					err_flux_low_new = 10.0**(0.4*(self.zp - magnitude[low_mag_points] - err_mag_low_new)) - flux_low_mag_points	
				
					diff1 = (np.sum(dmsq_high/(err_flux_high_new**2)) / (len(high_mag_points)))-1.0

					diff2 = (np.sum(dmsq_low/(err_flux_low_new**2)) / (len(low_mag_points)))-1.0

					return diff1, diff2


				#high_mag_points = np.where((ls<=tpoints) & (tpoints>=rs))[0]

				high_range = np.linspace(0.1,0.5,10)

				for i in range(len(high_range)):
					high_range_point = high_range[i]
					high_mag_points = np.where(mag>=np.max(mag)*high_range_point)[0]
					if len(high_mag_points)>=self.ndim:
						high_mag_points = high_mag_points

				low_mag_points = np.setdiff1d(points,high_mag_points)

				chi2_elements_high = chi2_elements_site[high_mag_points]
				flux_high_mag_points = flux[high_mag_points]
				err_mag_high = err_mag[high_mag_points]
				err_flux_high = 10**(0.4*(self.zp - magnitude[high_mag_points] - err_mag_high))- flux_high_mag_points
				dmsq_high = chi2_elements_high*(err_flux_high**2)


				chi2_elements_low = chi2_elements_site[low_mag_points]
				err_mag_low = err_mag[low_mag_points]
				flux_low_mag_points = flux[low_mag_points]
				err_flux_low = 10**(0.4*(self.zp - magnitude[low_mag_points] - err_mag_low))- flux_low_mag_points
				dmsq_low = chi2_elements_low*(err_flux_low**2)

				if self.ndim>len(points):
					scale = np.sqrt(chi2/(len(points)-self.ndim-2))
				else:
					scale = np.sqrt(chi2/(len(points)))

				eps = 0.0001
				x0 = (scale,eps)
				
				value = scipy.optimize.least_squares(func1,x0,bounds=([0.0,0.0], [8.0,1.0]),args=(dmsq_high, dmsq_low, err_mag_high, err_mag_low, magnitude,flux_high_mag_points, flux_low_mag_points),max_nfev=10000000)
				print('value=', value.x, site)
				scalev = (value.x[0])
				epsv = (value.x[1])
				savedict[site] = [scalev,epsv]
				writedict('renormalise_values.dat',savedict)

				sigma = scalev*np.sqrt(err_mag**2+epsv**2)
				y = 10.0**(0.4*(self.zp-magnitude))
				dy = 10.0**(0.4*(self.zp-magnitude+sigma))-y
				self.data[site] = (self.data[site][0],y,dy)
				chi2_new,_ = self.chi2_calc(p,source=site)
				print('chi2_new=', chi2_new)

