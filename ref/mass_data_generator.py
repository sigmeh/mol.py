#!/usr/bin/env python
'''
Parse mass data from mass_data.txt into element objects 
 *changed abundance/abundances to intensity/intensities for consistency across modules (when building fragments

>>> elements['Os'].atomic_num
76

>>> elements['Ru'].masses
[95.907597681, 97.905287111, 98.905939307, 99.904219664, 100.905582219, 101.904349503, 103.905430145]

>>> elements['Si'].intensities
[0.9223, 0.0467, 0.031]

'''
import numpy as np
import sys

elements = {}

class element(object):
	def __init__(self, atomic_num, mass, intensity):
		self.atomic_num = atomic_num
		self.masses = [mass]
		self.intensities = [intensity]
	

def generate_mass_data():

	with open('ref/mass_data.txt','r') as f:
		mass_data = [x.strip() for x in f.read().split('\r') if x]
	
	els = {}
	
	for line in mass_data:

		i = line.split(' ')
		atomic_num = int(i[0])
		el = i[2]
		mass = float(i[3].strip('#'))
		intensity = float(i[5])/100.	
		
		if el not in elements:
			elements[el] = element(atomic_num=atomic_num, mass=mass, intensity=intensity)
		else:
			elements[el].masses.append(mass)
			elements[el].intensities.append(intensity)

	for el in elements:
		if len(elements[el].masses) == 1:
			elements[el].molar_mass = elements[el].masses[0]
		elif all([x == 0. for x in elements[el].intensities]):
			elements[el].molar_mass = sum(elements[el].masses)/len(elements[el].masses)
		else:
			elements[el].molar_mass = sum(np.array(elements[el].masses) * np.array(elements[el].intensities))
	
	return elements


def main():
	return generate_mass_data()


if __name__ == '__main__':
	main()
	
