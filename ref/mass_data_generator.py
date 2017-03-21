#!/usr/bin/env python
'''
Parse mass data from mass_data.txt into element objects 
'''
import numpy as np
elements = {}

class element(object):
	def __init__(self,atomic_num,mass,abundance):
		self.atomic_num = atomic_num
		self.masses = [mass]
		self.abundances = [abundance]
	

def generate_mass_data():

	with open('ref/mass_data.txt','r') as f:
		mass_data = [x.strip() for x in f.read().split('\r') if x]
	
	els = {}
	
	for line in mass_data:
		i = line.split(' ')
		atomic_num = int(i[0])
		el = i[2]
		mass = float(i[3].strip('#'))
		abundance = float(i[5])/100.	
		
		if el not in elements:
			elements[el] = element(atomic_num=atomic_num,mass=mass,abundance=abundance)
		else:
			elements[el].masses.append(mass)
			elements[el].abundances.append(abundance)

	for el in elements:
		if len(elements[el].masses) == 1:
			elements[el].molar_mass = elements[el].masses[0]
		elif all([x == 0. for x in elements[el].abundances]):
			elements[el].molar_mass = sum(elements[el].masses)/len(elements[el].masses)
		else:
			elements[el].molar_mass = sum(np.array(elements[el].masses) * np.array(elements[el].abundances))
	return elements


def main():
	return generate_mass_data()


if __name__ == '__main__':
	main()
	
