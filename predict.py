#!/usr/bin/env python
'''
Generate spectra objects for all permutations of elements/fragments, based on min/max values defined in frag_chart below
frag_chart can be modified directly; the result is a dictionary object called spectra, indexed by fragment name 
	(frag_name is given in frag_chart, ref/known_fragments.txt, or an alphabetic string concatenation of the formula)

Usage:
	>>> import predict.py
	>>> predict.main()

'''
from numpy import *
from ref import mass_data_generator
from ref import plot_mass_data
from ref import atom_list_from_formula
from ref import peak_predictor
from ref import plot_multiple_spectra
import mol
import sys
import json
import subprocess as sp
import time

start=time.time()
ca=0


frag_chart = '''
		#Input desired elements, with min/max numbers of atoms/fragments for this test (comma-separated)
		#Commented lines are ignored
		#element, 	min, 	max
	
		Ru, 		2, 		2,
		OPiv, 		4, 		4,	
		#H, 		0, 		1,	
		#HOPiv, 	0, 		3, 
		#OAc, 		0, 		6,
		#HOAc, 		0, 		6
	'''	

spectra = {}
class spectrum(object):
	def __init__(self,frag_name):
		self.name = frag_name
		self.masses = []
		self.intensities = []


def handle_error(error):
	if error:
		print 'Error:',error
		sys.exit()
		

def get_known_fragments():
	'''known_fragments.txt contains references to atomic constitution of chemical fragments' shorthand notation; e.g., 'OAc' == O2CCH3  '''
	with open ('ref/known_fragments.txt','r') as f:
		known_fragments = {x.split(',')[0].strip():x.split(',')[1].strip() for x in f.read().split('\n') if x and not x.startswith('#')}
	return known_fragments


def main():
	
	'''----------elements-----------'''
	elements = mass_data_generator.generate_mass_data()
	'''
	try:
		with open('ref/precomputed_fragments.txt','r') as f:
			data = json.loads(f.read())
	except: 
		data = {}
	'''
	'''----------elements-----------'''
	
	
	'''----------known_fragments-----------'''	
	'''get atom list from formula for all known_fragments'''
	known_fragments_dict = get_known_fragments()	
	for known_frag in known_fragments_dict:
		formula = known_fragments_dict[known_frag]
		atom_list, formula_dict, formatted_formula, error = atom_list_from_formula.get(formula, elements)
	
		s = peak_predictor.predict(
			atom_list		= atom_list,
			elements		= elements, 
			formula 		= formula, 
			name 			= known_frag, 
			spectra			= spectra
		)
		
		spectra[s.name] = s
		
		''' To access class attributes: 
			>>> spectra['OAc'].masses 
			[59.013304, 60.016659, 60.019581, 61.01755] '''
			
	'''----------known_fragments-----------'''	
	

	
	frags = {}
	
	
	class frag_obj(object):
		'''	Define fragment objects for predictor book-keeping; includes min and max values (fragment quantity), if any'''
		'''	frags dictionary allows direct method calls; e.g., frags['Ru'].min'''
		
		def __init__(self, frag_data):
			self.name = frag_data[0]
			self.min = int( frag_data[1] )
			self.max = int( frag_data[2] )
			
			frags[frag_data[0]] = self
	
	frag_split = [ frag_obj( [ y.strip() for y in x.split(',')]) for x in frag_chart.replace('\t','').split('\n') if x and not x.startswith('#')]
	

	
	frags_to_test = [x for x in sorted(frags)]	
	

	
	all_frags = {}	#keep running count of unique fragments for testing purposes
	def conf_frag(fraglist):
		'''Count frag instances'''
		new_entry = ''.join([x+str(fraglist.count(x)) for x in set(fraglist)])
		if new_entry in all_frags:
			all_frags[new_entry] +=1
		else:
			all_frags[new_entry] = 1
		return new_entry
	
	
	def frags_parameters_ok(frag_list):
		''' 1. Ensure frag count is within threshold set by min/max parameters for each frag'''
		''' 2. Ensure all fragments are present which are required to be, as per min parameter > 0'''
		unique_frags = set(frag_list)
		for i in unique_frags:
			if frag_list.count(i) < frags[i].min or frag_list.count(i) > frags[i].max:
				return False
		for i in frags:
			if i not in unique_frags and frags[i].min > 0:
				return False
		return True

	
	frag_list = []
	count = 0

	def gen_spec(frag_list,count):
		print 'now'
		if frag_list:	
			lastfrag = frag_list[-1]
			endcount = frag_list.count(lastfrag)
			if endcount > frags[lastfrag].max: 	
				frag_list.pop()
				return frag_list,count
		for i in range(count,len(frags_to_test)):	
			
			frag_list.append(frags_to_test[i])

			frag_name = conf_frag(frag_list)
			
			global ca 
			ca+=1
			if frags_parameters_ok(frag_list):
				formula = ''.join([ x + str(frag_list.count(x)) if x in elements else '(%s)%d' % (x,frag_list.count(x) ) for x in set(frag_list) ])
				name = formula
				atom_list, formula_dict, formatted_formula, error = atom_list_from_formula.get(formula, elements)
				handle_error(error)
				
				''' get spectrum from atom_list '''
				#'''
				s = peak_predictor.predict(
					atom_list		= atom_list,
					elements		= elements, 
					formula 		= formula, 
					name 			= name, 
					spectra			= spectra
				)
		
				spectra[s.name] = s
				#'''
			gen_spec(frag_list,count)
			count+=1
		else:
			if frag_list:
				frag_list.pop()
			return frag_list,count
	
	
	gen_spec(frag_list,count)


	#s = sorted(spectra)[0]
	#plot_mass_data.plot_data(spectra[s],spectra[s].name)
	
	#plot_multiple_spectra.plot(spectra)
	
	#for i in sorted(spectra):
	#	print i
	
	#print ca, len(spectra)
	return spectra
	
	
	pass
if __name__ == '__main__':
	main()

print 'Elapsed time: ',round(time.time() - start,3),'s'