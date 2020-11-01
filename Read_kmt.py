
import sys
import os
import numpy as np 

def read_data(data_directory,t_range=None,minimum_points=5,max_uncertainty=0.05, renorm=None):

	"""Load and assemble data.

		This is a personal file to load KMT data only uses (date,mag,err_mag).

		t_range (optional) is a tupple of (start_date,end_date).

		0.01 mag is added in quadrature to all uncertainties.

		Returns a dictionary of (source,source_data) pairs, where each source_data is a tupple of (date,flux,err_flux).
	"""
	print('Reading data:')

	if t_range is None:
		t_range = (-1.e6,1.e6)

	data = {}

	files = [f for f in os.listdir(data_directory) if f.split('.')[-1]=='pysis' and f.split('.')[0][-1] != 'V' ]
	print(files)

	for f in files:

		source = f[:8]
		source = source.replace('_','-')
		datafile = data_directory+'/'+f

		if os.path.exists(datafile):
			d = np.loadtxt(datafile,comments=['<','#'],usecols=(0,3,4))
			t = d[:,0]
			d2 = np.copy(d[:,2])
			if t[0] > 2440000:
				t -= 2450000
			points = np.where((t_range[0] < t) & (t< t_range[1]))[0]

			t = t[points]
			if renorm == None:
				y = 10.0**(0.4*(25 - d[points,1]))
				d[points,2] = np.sqrt(d[points,2]**2 + (0.01)**2)
				dy = 10.0**(0.4*(25 - d[points,1] + d[points,2])) - y
				points = np.where(dy/y < max_uncertainty)[0]
			else:
				y = 10.0**(0.4*(25 - d[points,1]))
				d[points,2] = renorm[source][0]*np.sqrt((d[points,2]**2 + renorm[source][1]**2))
				dy = 10.0**(0.4*(25 - d[points,1] + d[points,2])) - y
				d2[points] = np.sqrt(d2[points]**2 + (0.01)**2)
				d2y = 10.0**(0.4*(25 - d[points,1] + d2[points])) - y
				points = np.where(d2y/y < max_uncertainty)[0]

			t = t[points]
			y = y[points]
			dy = dy[points]
			if len(t) >= minimum_points:
				data[source] = (t,y,dy)

	return data




