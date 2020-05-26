from matplotlib import pyplot as plt 
import scipy.stats as stats
import numpy as np  
import math, glob, os

def plot_histograms_within_a_multisim(startpath, box, value):
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/*/%s_BOX_0.dat" % (startpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/%s_BOX_1.dat" % (startpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_histogram.png" % (startpath, value)
	plt.title(startpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_individual_replicas(startpath, box, value):
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/%s_BOX_0.dat" % (startpath, value)
	elif box == 1:
		checkIfEmpty = "%s/%s_BOX_1.dat" % (startpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	
	for name in glob.glob(checkIfEmpty):
		tokens = name.split(os.sep)
		tokens.remove(tokens[len(tokens)-1])
		newName = os.sep.join(tokens)
		filename = "%s/%s_BOX_%s_histogram_exchange_vs_nonexchange.png" % (newName, value, box)
		#print filename	
		plt.title(startpath) 
		plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_all_replicas(startpath, box, value):
	#	print startpath, box
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/*/*/%s_BOX_0.dat" % (startpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/*/%s_BOX_1.dat" % (startpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
#	print "Right before for loop"
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_BOX_%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (startpath, value, box)
	plt.title(startpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_within_a_multisim_num_molecules(startpath, box, value):
	global checkIfEmpty
	floats = []
	checkIfEmpty = "%s/*/%s" % (startpath, value)
	
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_histogram.png" % (startpath, value)
	plt.title(startpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_individual_replicas_num_molecules(startpath, box, value):
	global checkIfEmpty
	floats = []
	checkIfEmpty = "%s/%s" % (startpath, value)

	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	
	for name in glob.glob(checkIfEmpty):
		tokens = name.split(os.sep)
		tokens.remove(tokens[len(tokens)-1])
		newName = os.sep.join(tokens)
		filename = "%s/%s_histogram_exchange_vs_nonexchange.png" % (newName, value)
		#print filename	
		plt.title(startpath) 
		plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_all_replicas_num_molecules(startpath, box, value):
	#	print startpath, box
	global checkIfEmpty
	floats = []

	checkIfEmpty = "%s/*/*/%s" % (startpath, value)

	index = 0
#	print "Right before for loop"
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		ag,bg,cg = stats.gamma.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.gamma.pdf(lnspc, ag, bg,cg)
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (startpath, value)
	plt.title(startpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()
