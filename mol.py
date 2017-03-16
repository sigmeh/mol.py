#!/usr/bin/env python
import re
import string
import sys
'''
Determine molecular weight from input formula based on iterative term expansion using regular expressions.
Formula may contain any of "( { [ ] } )" separators and spaces. 

'''

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
	valid_chars = string.ascii_letters + string.octdigits + '()'
	if f.count('(') != f.count(')') or any([x not in valid_chars for x in f]): 
		print 'invalid formula'
		return None

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
				
			subfrags_exp = ''.join([ int(n)*y for y in re.findall( '[A-Z][a-z]?', x ) ])
			esc_x = ''.join([ '\\'+y if y in '()' else y for y in x ])	# Need to escape any parentheses for regex substitution (on next line)
			f = re.sub( esc_x, subfrags_exp, f, count=1 )
	
	mol_f = { i:f.count(i) for i in list( set( re.findall('[A-Z][a-z]?',f) ) ) }
	return mol_f

def main():
	
	
	#f = 'Fe4{O2 CCH3 }3 Cl[ HO2C(CH2)C(O)CH3 ]2 F2I2'
	
	#'''
	if len(sys.argv) == 2:	
		f = str(sys.argv[1])
	else:
		f = str(raw_input('Enter formula: '))
	#'''	
		
	mol_f = calc(f)
	if emp_f:
		print 'Molecular formula:',' '.join([x+str(emp_f[x]) for x in emp_f])
	

if __name__ == '__main__':
	main()