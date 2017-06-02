#!/usr/bin/env python
'''
Generate predicted mass spec peak distribution from atom_list
through iterative matrix multiplication of applicable isotope data

Returns a spectrum object s, whose attributes can be accessed as: 

	>>> s.name
	OAc
	
	>>> s.formula
	O2CCH3
	
	>>> s.masses 
	[59.013304, 60.016659, 60.019581, 61.01755]
	
	>>> s.intensities
	[97.413, 2.167, 0.029, 0.391]
'''

from numpy import *
import sys


class spectrum(object):
	def __init__(self):
		self.masses 		= []
		self.intensities 	= []


def predict( atom_list, elements, formula , name, spectra ):
	'''Initialize spectrum and populate with first atom/frag in list'''
	
	if len(atom_list) > 1:
		s = spectrum()
		
		if atom_list[0] in spectra:
			data_source = spectra[ atom_list[0] ]
		elif atom_list[0] in elements:
			data_source = elements[ atom_list[0] ]
		else:
			print 'Unknown fragment/atom: %s' % atom_list[0]
			sys.exit()
		
		s.masses = array( data_source.masses )
		s.intensities = array( data_source.intensities )
		
		for atom in range(1,len(atom_list)):
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
			m2 = elements[atom_list[atom]].masses
			m2 = matrix( vstack(( ones(shape(m2)) , m2 )) )
			
			masses = [round(x,6) for x in list((m1*m2).flat)]
		
			i2 = elements[atom_list[atom]].intensities
			intens = list(outer(s.intensities, i2).flat)
	
			s = spectrum()
			for data_el in set(masses):
				s.masses.append(data_el)
				s.intensities.append( sum( [intens[x] for x in range(len(intens)) if data_el == masses[x]] ) )
				if s.intensities[-1] < .005:
					s.masses.pop()
					s.intensities.pop()
		
		s.intensities = [round(x,6) for x in s.intensities]
		
	else:
		s = spectrum()
		s.masses = elements[atom_list[0]].masses
		s.intensities = elements[atom_list[0]].intensities
	
	''' Sort spectral peaks (both of masses/intensities lists), by mass; normalize intensities to 100. '''
	z = sorted( zip( s.masses, s.intensities ) )
	s.masses = [x[0] for x in z]
	s.intensities = [ round( x[1]/sum(s.intensities)*100. , 3) for x in z]
	
	s.formula = formula
	s.name = name
	
	return s
	
	
	
def main():
	pass
if __name__ == '__main__':
	main()