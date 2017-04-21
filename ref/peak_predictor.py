#!/usr/bin/env python
'''
Generate predicted MS peak distribution from atom list
through iterative matrix multiplication of relevant isotope data
'''

from numpy import *
import matplotlib as mpl
mpl.use('TkAgg')
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt
#get_current_fig_manager().window.raise_()



class spectrum(object):
	def __init__(self):
		self.masses = []
		self.intensities = []
	

def check_vals(s):
	for i in  sorted([(s.masses[x],s.intensities[x]) for x in range(len(s.masses))]):
		print i


def plot_data(s):
	fig, ax = plt.subplots()
	ax.stem(s.masses, s.intensities, markerfmt=' ')
	plt.show()


def ascii_print(s):
	pass


def predict(atoms,elements):

	s = spectrum()
	s.masses = array(elements[atoms[0]].masses)
	s.intensities = array(elements[atoms[0]].abundances)

	for atom in range(1,len(atoms)):
		'''	Iteratively compute predicted spectrum for the addition of each atom. 
			Adding "ones" column/row placeholders to mass vectors allows matrix multiplication 
			to generate complete set of arithmetic mass permutations, based on each new element 
			added to the previously-condensed spectrum.
			The outer product of intensity matrices is taken to generate an intensity matrix
			for new peaks. 
			Mass and intensity matrices are flattened and condensed (the intensities of peaks 
			with equal masses are added together, and any peaks below 0.5% threshold are removed). 
		
		'''
		
		
		m1 = matrix( vstack(( s.masses, ones(shape(s.masses)) )) ).T	
		m2 = elements[atoms[atom]].masses
		m2 = matrix( vstack(( ones(shape(m2)) , m2 )) )
			
		masses = [round(x,6) for x in list((m1*m2).flat)]
		
		i2 = elements[atoms[atom]].abundances
		intens = list(outer(s.intensities, i2).flat)
	
		s = spectrum()
		for data_el in set(masses):
			s.masses.append(data_el)
			s.intensities.append( sum( [intens[x] for x in range(len(intens)) if data_el == masses[x]] ) )
			if s.intensities[-1] < .005:
				s.masses.pop()
				s.intensities.pop()
		
	
	check_vals(s)	
		
	plot_data(s)	
	
	ascii_print(s)
	
	
def main():
	pass
if __name__ == '__main__':
	main()