#!/usr/bin/env python
import re
import string
import sys
from ref import mass_data_generator
from ref import peak_predictor_v01

'''
Determine molecular weight from input formula based on iterative term expansion using regular expressions.
Formula may contain any of "( { [ ] } )" separators and spaces. 

Usage:
	$ python mol.py [formula]
	
	$ python
	>>> import mol
	>>> mol.main( formula_here )
	
	$ python mol.py 	#display prompt when invoked without arguments
	Enter formula: 	


Notes:
	During term expansion, ref/known_fragments.txt is consulted, and any instances of fragments in the formula,
	whose definitions are contained in the known_fragments.txt file, are replaced first. 
	
'''

def get_known_fragments():
	with open ('ref/known_fragments.txt','r') as f:
		known_fragments = {x.split(',')[0].strip():x.split(',')[1].strip() for x in f.read().split('\n') if x and not x.startswith('#')}
	return known_fragments

def calc(f):
	'''	1. Convert separators to parentheses and remove spaces
		2. Ensure formula contains only valid characters
		3. Employ term expansion 
			a. Expand direct terms (e.g., O2 --> OO)
			b. Expand parenthesized terms (e.g.,  (CHHH)2 --> CCHHHHHH  ) iteratively, until no nested terms remain. 
		4. Generate molecular formula dictionary
		5. Test against valid elements and calculate molecular weight
	'''
	
	f = re.sub( '\{|\[' , '(', re.sub( '\}|\]' , ')' , f ) ).replace(' ','').replace('()','')	
	valid_chars = string.ascii_letters + '0123456789()'
	
	if f.count('(') != f.count(')') or any([x not in valid_chars for x in f]): 
		print 'invalid formula'
		sys.exit()
	
	known_fragments = get_known_fragments()
	for frag in known_fragments:
		f = re.sub( frag, known_fragments[frag] , f )
	
	regex = [	'[A-Z][a-z]?\d+',
				'\(\w*\)\d*'	]
			
	r = 0
	while r < len(regex): 
		frags = re.findall( regex[r], f )
		if not frags:
			r+=1
			continue
				
		for x in frags:
			n = re.search( '\d+', x )
			if n: 
				n = n.group()
			else:
				n = 1	
				
			subfrags_expansion = ''.join([ int(n)*y for y in re.findall( '[A-Z][a-z]?', x ) ])
			escaped_x = ''.join([ '\\'+y if y in '()' else y for y in x ])	# Need to escape any parentheses for regex substitution (on next line)
			f = re.sub( escaped_x, subfrags_expansion, f, count=1 )
	
	all_els = re.findall('[A-Z][a-z]?',f)
	if not all_els:
		print 'formula has no valid elements'
		sys.exit()
	els = list(set(all_els))
	mol_f = {i:all_els.count(i) for i in els}
	mol_f_formatted = ' '.join([x+str(mol_f[x]) for x in mol_f])
	#print mol_f
	return mol_f, all_els


def main(*args):
	'''
	1. Obtain formula via arguments or prompt
	2. Generate element:element-count dictionary
	3. Generate element mass standard data from mass_data_generator.py
	4. Ensure all elements in formula are valid
	5. Determine molecular weight (average, based on isotope abundance)
	6. Predict exact isotopic distribution from known element masses using peak_predictor.py
	'''
	if args:
		f = args
	elif len(sys.argv) == 2:
		f = sys.argv[1]
	else:
		
		# testing area
		f = 'OPiv'
		#f = 'OPiv'
		##
		##
		#f = str(raw_input('Enter formula: '))
		##
		##
		
	print 'Formula input:',f
	mol_f,all_els = calc(f)
	'''mol_f is molecular formula dictionary (key = atom, value = # atoms); all_els is expanded list of atoms (containing repeats)'''
	
	elements = mass_data_generator.generate_mass_data()
	#print elements
	''' elements['Ru'].masses is list of exact masses'''
	
	for x in mol_f:
		if x not in elements:
			print 'Formula contains invalid elements:','\"%s\"'%x
			sys.exit()

	#if any([x not in elements for x in mol_f]):
	#	print 'formula contains invalid elements:',x
	#	sys.exit()
	
	if mol_f:
		print 'Molecular formula:',' '.join([x+str(mol_f[x]) for x in mol_f])
	
	mw = sum([mol_f[el]*elements[el].molar_mass for el in mol_f])
	print 'Molecular weight:',mw
	
	#print all_els
	#print mol_f
	formatted_mol_f = ''.join(sorted([x+str(mol_f[x]) for x in mol_f]))
	
	
	peak_profile = peak_predictor_v01.predict(atoms=all_els,elements=elements,input_formula=f,formatted_mol_f=formatted_mol_f)
	print formatted_mol_f
	sys.exit()
	print 'peak_profile:',peak_profile
	
	

if __name__ == '__main__':
	main()