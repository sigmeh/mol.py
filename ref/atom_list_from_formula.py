#!/usr/bin/env python
'''
Takes as input a potential chemical formula, cleans up formula if necessary/possible and evaluates formula's validity
Returns atom list (total list of all individual atom identities)

'''
import re
import sys
import string
from ref import mass_data_generator

elements = mass_data_generator.generate_mass_data()

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
	
	
	''' 1. Convert separators to parentheses and remove spaces'''
	f = re.sub( '\{|\[' , '(', re.sub( '\}|\]' , ')' , f ) ).replace(' ','').replace('()','')	
	valid_chars = string.ascii_letters + '0123456789()'
	
	
	'''2. Ensure all parentheses are closed and formula contains only valid characters'''
	if f.count('(') != f.count(')') or any([x not in valid_chars for x in f]): 
		error = 'invalid formula'
		return None, None, error
	
	
	'''3.0 Substitute any known fragments in formula with their known formula from ./known_fragments.txt '''
	known_fragments = get_known_fragments()
	for frag in known_fragments:
		f = re.sub( frag, known_fragments[frag] , f )
	
	
	'''Regex expressions are 3.a and 3.b selectors, respectively'''
	regex = [	'[A-Z][a-z]?\d+',
				'\(\w*\)\d*'	]
	
	
	''' For each regex expression, iteratively find all matching terms; this will take care of nested fragments by iterative term expansion'''		
	r = 0
	while r < len(regex): 
		print f
		frags = re.findall( regex[r], f )
		if not frags:
			r+=1
			continue
				
		for x in frags:
			''' Obtain number of preceding atoms; n=1 if no number is found'''
			n = re.search( '\d+', x )
			if n: 
				n = n.group()
			else:
				n = 1	
			
			''' 3.a. Expand number of preceding atoms based on n (found above); return as string'''	
			subfrags_expansion = ''.join([ int(n)*y for y in re.findall( '[A-Z][a-z]?', x ) ])
			
			''' Must escape any parentheses for regex substitution; replace first instance (x) with '''
			escaped_x = ''.join([ '\\'+y if y in '()' else y for y in x ])	
			
			''' Replace each x (escaped) in frags with corresponding expanded term; e.g., O2 --> OO'''
			f = re.sub( escaped_x, subfrags_expansion, f, count=1 )
	
	atom_list = re.findall('[A-Z][a-z]?',f)
	if not atom_list:
		error = 'Formula has no valid elements'
		return None, None, error
	
	''' Create dictionary mol_f containing entries for element: element-count'''
	els = list(set(atom_list))
	formula_dict = { i : atom_list.count(i) for i in els}
	
	formatted_formula_dict = ' '.join( [ x + str(formula_dict[x]) for x in formula_dict ] )
	
	print formatted_formula_dict
	sys.exit()
	
	for x in formula_dict:
		if x not in elements:
			error = 'Formula contains invalid elements:','\"%s\"'%x
			return None, None, error
	
	return atom_list, formula_dict, formatted_formula_dict, None


def get(formula):
	
	
	
	atom_list, formula_dict, formatted_formula_dict, error = calc(formula)
	
	
	return atom_list
	
	
	
	

def main():
	pass
	
if __name__ == '__main__':
	main()