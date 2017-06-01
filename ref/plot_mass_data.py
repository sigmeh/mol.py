#!/usr/bin/env python
from numpy import *
import matplotlib as mpl
mpl.use('TkAgg')
#mpl.use('Qt4Agg')
import matplotlib.pyplot as plt


def plot_data(s,formula):
	fig, ax = plt.subplots()
	ax.stem(s.masses, s.intensities, markerfmt=' ')
	ax.set_title('Predicted isotopic distribution for %s' %formula)
	ax.set_xlabel('m/z')
	ax.set_ylabel('intensities')
	plt.show()


def main():
	plot_data(s,formula)
	
	
	
	
	
	
	
	
	pass
if __name__ == '__main__':
	main()