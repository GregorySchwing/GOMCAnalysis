from matplotlib import pyplot as plt 
import scipy.stats as stats
import scipy.optimize as opt
import numpy as np  
import math, glob, os
import re
from itertools import islice
import sympy as sy
from scipy.integrate import simps
import itertools    
import itertools as IT
import numdifftools as nd
from lmfit import Minimizer, Parameters
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D


def visual_inspection_plots(dottedtpath, box, value):
	# empty dictionary
	my_dict = {}
	global checkIfEmpty
	floats = []
	mins = []
	temperatures = []
	betas = []
	if box == 0:
		checkIfEmpty = "%s/*/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/%s_BOX_1.dat" % (dottedtpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		temperatures.append(re.findall("_(\d*[.,]?\d*)/", name))
		floats.append(np.fromfile(name, dtype='float64', count=-1, sep='\n'))
		index = index+1

	for arrays in floats:
		mins.append(np.min(arrays))
	sort_index = np.argsort(mins)	
	print sort_index
	#print floats[1][floats[1] > -6239910.23]
	index = 0
	#sorted_dict = list(sorted(my_dict.items()))
	temp = iter(sort_index) 
	#print sorted_dict	
	for i in temp:	
		print i
		overlaps = []
		try:
			next_i = next(temp)
			a_greater_than_right_replicas_min = floats[i][floats[i]>=mins[next_i]]
			a_less_than_left_replicas_max = a_greater_than_right_replicas_min[a_greater_than_right_replicas_min<=np.max(floats[i])]
			b_greater_than_right_replicas_min = floats[next_i][floats[next_i]>=mins[next_i]]
			b_less_than_left_replicas_max = b_greater_than_right_replicas_min[b_greater_than_right_replicas_min<=np.max(floats[i])]

			overlaps.append(a_less_than_left_replicas_max)
			overlaps.append(b_less_than_left_replicas_max)

			#print temperatures
			
			#betas.append(1./(float(temperatures[i][0]) * 1.380649E-23))
			#betas.append(1./(float(temperatures[next_i][0]) * 1.380649E-23))
			betas.append(1./(float(temperatures[i][0])))
			betas.append(1./(float(temperatures[next_i][0])))
			
			index = index+1

			##### Calculate the number of bins to use, by taking the average of two auto bins ####
### Bring this back once you figure our how to set up the for loops
			sample_edges1 = np.histogram_bin_edges(overlaps[index-1], bins='auto', range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))
			sample_edges2 = np.histogram_bin_edges(overlaps[index],  bins='auto', range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))

			numBins = sample_edges1.size + sample_edges2.size / 2


			########################################################################################

			##### Histogram the data ###############################################

			unnormedHist1, unnormedBin_edges1 = np.histogram(overlaps[index-1], bins=numBins, range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))

			unnormedHist2, unnormedBin_edges2 = np.histogram(overlaps[index],  bins=numBins, range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))

			########################################################################

			##### Let's find the bin with the smallest difference between the two histograms
			##### This will be our 'middle' bin ############################

			middleBin = np.argmin(np.absolute(np.subtract(unnormedHist1, unnormedHist2)))

			################################################################

			##### Let's digitize our inputs to make for easier trimming of data

			binIndicesOfLeftReplicaPoints = np.digitize(overlaps[index-1], bins=unnormedBin_edges1)
			binIndicesOfRightReplicaPoints = np.digitize(overlaps[index], bins=unnormedBin_edges2)

			##################################################################

			##### Lets find the center of the overlaps ###

			kSetsHist1 = np.array(list(getAllWindows(unnormedHist1)))
			kSetsHist2 = np.array(list(getAllWindows(unnormedHist2)))

			z_val = []
			y_val = []
			x_val = []
			minimum_index = 0
			minimum_diff = float("inf")
			for i in range(kSetsHist1.size):
    				y_val.append(np.true_divide(np.sum(np.absolute(np.subtract(kSetsHist1[i], kSetsHist2[i]))), kSetsHist1[i].size))
				x_val.append(abs(numBins - kSetsHist1[i].size))
				z_val.append(i)
				final_val = y_val[i]
				if (final_val < minimum_diff and kSetsHist1[i].size > numBins/4):
					minimum_diff = final_val
					minimum_index = i


			#area_linear = simps(y_val, x=x_val)

			indices1 = find_sub_list(kSetsHist1[minimum_index].tolist(), unnormedHist1.tolist())
			indices2 = find_sub_list(kSetsHist2[minimum_index].tolist(), unnormedHist2.tolist())

			intersectionOfResultArrays = np.intersect1d(indices1, indices2)	
			print intersectionOfResultArrays

			lowerBound = intersectionOfResultArrays[0]
			
			upperBound = intersectionOfResultArrays[1]

			##### Lets remove points that fall outside our overlap bins ###

			trimmedLeftReplica = overlaps[index-1][(binIndicesOfLeftReplicaPoints>lowerBound)&(binIndicesOfLeftReplicaPoints<upperBound)]	
			trimmedRightReplica = overlaps[index][(binIndicesOfRightReplicaPoints>lowerBound)&
(binIndicesOfRightReplicaPoints<upperBound)]	

			################################################################

			##### Let's histogram the trimmed data in the original histogram ranges to make sure we trimmed the edge blocks

			unnormedhist1, unnormedbin_edges1 = np.histogram(trimmedLeftReplica, bins=unnormedBin_edges1, range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))
			unnormedhist2, unnormedbin_edges2 = np.histogram(trimmedRightReplica, bins=unnormedBin_edges2, range=(mins[next_i], np.amax(a_less_than_left_replicas_max)))

			print "Notice the 0's on the edges.  Trimming worked."
			print unnormedhist1
			print unnormedhist2

			################################################################

			##### Let's remove the edge blocks from our range we use to histogram our data so we don't divide by 0/0 ###
	
			trimmedEdges1 = unnormedBin_edges1[lowerBound:upperBound]
			trimmedEdges2 = unnormedBin_edges2[lowerBound:upperBound]

			unnormedTrimmedhist1, unnormedTrimmedbin_edges1 = np.histogram(trimmedLeftReplica, bins=30)
			unnormedTrimmedhist2, unnormedTrimmedbin_edges2 = np.histogram(trimmedRightReplica,  bins=30)

			print "Notice how the histograms are identical to above, except the 0 bins are gone."
			print unnormedTrimmedhist1
			print unnormedTrimmedhist2


			################################################################

			bincenters1 = np.array(0.5*(unnormedTrimmedbin_edges1[1:]+unnormedTrimmedbin_edges1[:-1]))
			bincenters2 = np.array(0.5*(unnormedTrimmedbin_edges2[1:]+unnormedTrimmedbin_edges2[:-1]))	

			###########################################################

			log_final_hist = np.log(np.true_divide(unnormedTrimmedhist2, unnormedTrimmedhist1))

			################	Figure 1a	###########################################################

			################### linear -(B2-B1)E ##########################

			betaTerm = -(betas[index]-betas[index-1])
				
			xcellMethodRun = bincenters1[-1] - bincenters1[0] 
			xcellMethodRunTimesBeta = xcellMethodRun * betaTerm
			xcellMethodEndVal = log_final_hist[0] + xcellMethodRunTimesBeta

			lnspcy = np.linspace(log_final_hist[0], xcellMethodEndVal, log_final_hist.size)

			#########################################################

			########## Linear LSqu Weights: variance in the logarithm of the ratio ln(rk) = ln(pk,2 / pk,1) Eq 33

			varianceOfRatioInLogSpace = np.array(np.true_divide(1, unnormedTrimmedhist1) - np.true_divide(1, unnormedTrimmedhist1.sum()) + np.true_divide(1, unnormedTrimmedhist2) - np.true_divide(1, unnormedTrimmedhist2.sum()))

			print "variance in the logarithm of the ratio ln(rk) = ln(pk,2 / pk,1) Eq 33"
			print varianceOfRatioInLogSpace
			#########################################################################

			################ Linear Least Sq ###############################

			guess = [betaTerm, 1]
			popt, pcov = opt.curve_fit(f = linearFunc, xdata = bincenters1-bincenters1[0], ydata = log_final_hist, p0 = guess, sigma = np.sqrt(varianceOfRatioInLogSpace))

			#######################################################

			############ Plot histogram ratio in logspace, -(B2-B1)E, & Linear Least Squares Fit

			# histogram ratio in logspace
			plt.errorbar(
    			x = bincenters1,
    			y = log_final_hist,
    			yerr = 2*np.sqrt(varianceOfRatioInLogSpace),
			color='red'
			)
			# -(B2-B1)E
			plt.errorbar(
    			x = bincenters1,
    			y = lnspcy,
			color='green'
			)

			# Linear Least Squares Fit
			plt.errorbar(
    			x = bincenters1,
    			y = linearFunc(bincenters1-bincenters1[0], *popt),
			color='blue'
			)

			colors = ['red', 'green', 'blue']
			lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
			labels = [r'ln$\dfrac{P_%d(E)}{P_%d(E)}$' % (index, index-1), r'$-(\beta_%d-\beta_%d)\mathrm{E}$' % (index, index-1), r'Fit to $y=b+aE$']
			plt.legend(lines, labels)


			plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			ylabe = r'$\mathrm{ln}[\mathrm{P}_%d(E)/\mathrm{P}_%d(E)$]' % (index, index-1)
			special = "OverlapOfReplica_{}_AndReplica_{}_linear_bin_{}_to_bin_{}".format(index, index-1, lowerBound, upperBound)
			filename = "%s/%s_BOX_%s_%s.png" % (dottedtpath, value, box, special)
			plt.tick_params(axis='both', which='major', labelsize=16)
			plt.xlabel(r'E($k_B$T)', fontsize=18)
			plt.ylabel(ylabe, fontsize=18)
			plt.title("Linear Visual Inspection : Replica %d vs Replica %d" % (index, index-1), fontsize=20, pad=15)
			plt.savefig(filename, bbox_inches='tight', pad_inches=0.5)
			plt.clf()
			plt.cla()

			################	Figure 1b	###########################################################

			ratio_final_hist = np.true_divide(unnormedTrimmedhist2, unnormedTrimmedhist1)

			################### exp( -(B2-B1)E ) ##########################

			betaTerm = -(betas[index]-betas[index-1])
				
			energyRun = bincenters1[-1] - bincenters1[0] 
			energyRunTimesBeta = energyRun * betaTerm


			energyRunTimesBetaEndVal = ratio_final_hist[0] + energyRunTimesBeta
			lnspcy = np.linspace(np.log(ratio_final_hist[0]), np.log(energyRunTimesBetaEndVal), ratio_final_hist.size, endpoint = 'true')

			#########################################################

			########### Nonlinear LSqu Weights - variance in the ratio of the histograms themselves Eq 34 ####


			varianceOfRatios = np.multiply(np.square(np.true_divide(unnormedTrimmedhist2 * unnormedTrimmedhist1.sum(), unnormedTrimmedhist1 * unnormedTrimmedhist2.sum())), varianceOfRatioInLogSpace)

			print "variance in the ratio of the histograms themselves Eq 34"
			print varianceOfRatios


			##########################################################################

			############### Nonlinear Least Sq ############################


			guess = [betaTerm, 1]
			popt, pcov = opt.curve_fit(f = exponentialFunc, xdata = bincenters1-bincenters1[0], ydata = ratio_final_hist, p0 = guess, sigma = np.sqrt(varianceOfRatios))
		

			#######################################################

			############ Plot histogram ratio in realspace, exp( -(B2-B1)E ), & Non-Linear Least Squares Fit

			# histogram ratio in realspace
			plt.errorbar(
    			x = bincenters1,
    			y = ratio_final_hist,
    			yerr = 2*np.sqrt(varianceOfRatios),
			color='red'
			)

			# exp( -(B2-B1)E )
			plt.errorbar(
    			x = bincenters1,
    			y = np.exp(lnspcy),
			color='green'
			)

			#Non-Linear Least Squares Fit
			plt.errorbar(
    			x = bincenters1,
    			y = exponentialFunc(bincenters1-bincenters1[0], *popt),
			color='blue'
			)

			filename = "%s/%s_BOX_%s_%s_data.txt" % (dottedtpath, value, box, special)
			np.savetxt(filename, list(zip(bincenters1, ratio_final_hist)), delimiter=',')

			colors = ['red', 'green', 'blue']
			lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='-') for c in colors]
			labels = [r'$\dfrac{P_%d(E)}{P_%d(E)}$' % (index, index-1), r'exp$(-(\beta_%d-\beta_%d)\mathrm{E})$' % (index, index-1), r'Fit to $y=\mathrm{exp}(b+aE)$']
			plt.legend(lines, labels)


			ylabe = r'$\mathrm{P}_%d(E)/\mathrm{P}_%d(E)$' % (index, index-1)
			special = "OverlapOfReplica_{}_AndReplica_{}_exponential_bin_{}_to_{}_bin".format(index, index-1, lowerBound, upperBound)
			filename = "%s/%s_BOX_%s_%s.png" % (dottedtpath, value, box, special)
			plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			plt.tick_params(axis='both', which='major', labelsize=16)
			plt.xlabel(r'E($k_B$T)', fontsize=18)
			plt.ylabel(ylabe, fontsize=20)
			plt.title("Nonlinear Visual Inspection : Replica %d vs Replica %d" % (index, index-1) , fontsize=20, pad = 15)
			plt.savefig(filename, bbox_inches='tight', pad_inches=0.5)
			plt.clf()
			plt.cla()

		except StopIteration:
			print "Reached last sim"

def single_replica_MLE(dottedtpath, box, value):
	# empty dictionary
	my_dict = {}
	global checkIfEmpty
	floats = []
	sortedtemps = []
	temperatures = []
	betas = []
	if box == 0:
		checkIfEmpty = "%s/*/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/%s_BOX_1.dat" % (dottedtpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		temperatures.append(re.findall("_(\d*[.,]?\d*)/", name))
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		
	for i in temperatures:
		print float(i[0])

	tempFloats = np.array(temperatures, dtype = float)
	for i in tempFloats:
		print i

	tempFloats = np.reshape(tempFloats, tempFloats.size)

	for i in tempFloats:
		print i

	sort_index = np.argsort(tempFloats)	
	print sort_index
	index = 0
	temp = iter(sort_index) 
	for i in temp:
		mean = np.mean(floats[i])
		var = np.var(floats[i])
		likelihood = 0
		for z in floats[i]:
			print "mean : %f ; var : %f, likelihood of val %f : %f )" % (mean, var, z, np.log(g_theta(z, mean, var)))
			likelihood+=np.log(g_theta(z, mean, var))
		print "likelihood of %d = %f" % (i, likelihood)
	
	

def combined_max_likelihood(dottedtpath, box, value):
	# empty dictionary
	my_dict = {}
	global checkIfEmpty
	floats = []
	sortedtemps = []
	temperatures = []
	betas = []
	if box == 0:
		checkIfEmpty = "%s/*/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/%s_BOX_1.dat" % (dottedtpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		temperatures.append(re.findall("_(\d*[.,]?\d*)/", name))
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		
	for i in temperatures:
		print float(i[0])

	tempFloats = np.array(temperatures, dtype = float)
	for i in tempFloats:
		print i

	tempFloats = np.reshape(tempFloats, tempFloats.size)

	for i in tempFloats:
		print i

	sort_index = np.argsort(tempFloats)	
	print sort_index
	index = 0
	temp = iter(sort_index) 
	for i in temp:
		print i
		overlaps = []
		try:
			next_i = next(temp)
			betas.append(1./(float(temperatures[i][0])* 1.380649E-23))
			betas.append(1./(float(temperatures[next_i][0])* 1.380649E-23))

			min_i = np.min(floats[i])
			min_next_i = np.min(floats[next_i])

			smallest = min(min_i, min_next_i)
			print smallest

			betaAverage = 0.5 * (1./(float(temperatures[i][0])) + 1./(float(temperatures[next_i][0])))
			deltaBeta = -(1./(float(temperatures[next_i][0])) - 1./(float(temperatures[i][0])))
			#deltaBeta = 0
			print np.min(np.subtract(floats[i], smallest))
			print np.min(np.subtract(floats[next_i], smallest))

			x_in = np.array([floats[i].size])			
			x_in_int = np.concatenate((x_in, np.subtract(floats[i], smallest)), axis = None)
			x_in_final = np.concatenate((x_in_int, np.subtract(floats[next_i], smallest)), axis = None)
			x_in_final_w_betaAv = np.concatenate((x_in_final, betaAverage), axis = None)


			aGuess = 1
			guessStart = [deltaBeta, aGuess]
			guessInt = np.array(guessStart)

			lowerBounds = guessInt.tolist()
			upperBounds = guessInt.tolist()

			lowerBounds[1] = .1
			upperBounds[1] = 'inf'	
			bnds = opt.Bounds(lowerBounds, upperBounds) 
			print bnds
			#print bnds.shape
			#print guess.size
			lik_model = opt.minimize(lik, x0 = guessStart, args = x_in_final_w_betaAv,method='Nelder-Mead')
			#plt.scatter(x,y)
			print lik_model

			#plt.show()
			#guess = [deltaBeta, 1]
			#popt, pcov = opt.curve_fit(f = linearFunc, xdata = floats[i], ydata = floats[i]*deltaBeta, p0 = guess)

			# Linear Least Squares Fit
			plt.errorbar(
    			x = floats[i],
    			y = floats[i]*deltaBeta,
			yerr = 2*np.sqrt(np.var(floats[i]*deltaBeta)),
			color='blue'
			)
			# -(B2-B1)E
			plt.errorbar(
    			x = floats[i],
    			y = lik_model.x[0]*floats[i],
			color='green'
			)

			print "delBeta true: %f" % deltaBeta
			print "delBeta max like : %f" % lik_model.x[0]

			plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			ylabe = r'$\mathrm{ln}[\mathrm{P}_%d(E)/\mathrm{P}_%d(E)$]' % (index, index-1)
			special = "LinearMaxLikeReplica_{}_AndReplica_{}_linear".format(index, index-1)
			filename = "%s/%s_BOX_%s_%s.png" % (dottedtpath, value, box, special)
			plt.tick_params(axis='both', which='major', labelsize=16)
			plt.xlabel(r'E($k_B$T)', fontsize=18)
			plt.ylabel(ylabe, fontsize=18)
			plt.title("Linear Visual Inspection : Replica %d vs Replica %d" % (index, index-1), fontsize=20, pad=15)
			plt.savefig(filename, bbox_inches='tight', pad_inches=0.5)
			plt.clf()
			plt.cla()

			index = index+1
		except StopIteration:
			print "Reached last sim"


def plot_histograms_across_two_multisims_individual_replicas(dottedtpath, box, value):
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/%s_BOX_1.dat" % (dottedtpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), (np.amax(floats[index]) - np.amin(floats[index]))/2)
		pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	
	for name in glob.glob(checkIfEmpty):
		tokens = name.split(os.sep)
		tokens.remove(tokens[len(tokens)-1])
		newName = os.sep.join(tokens)
		filename = "%s/%s_BOX_%s_histogram_exchange_vs_nonexchange.png" % (newName, value, box)
		#print filename	
		plt.title(dottedtpath) 
		plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_all_replicas(dottedtpath, box, value):
	#	print dottedtpath, box
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/*/*/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/*/%s_BOX_1.dat" % (dottedtpath, value)
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
		parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), (np.amax(floats[index]) - np.amin(floats[index]))/2)
		pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		if "dev_multisim" in name:
			plt.plot(lnspc, pdf_gamma, linestyle="dashed", linewidth=2)
		elif "nonMPInonPT" in name:
			plt.plot(lnspc, pdf_gamma, linestyle="dotted", linewidth=2)
		else:
			plt.plot(lnspc, pdf_gamma, linewidth=1)		
		index = index+1
	filename = "%s/%s_BOX_%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (dottedtpath, value, box)
	#plt.title("NVT Butane") 
	#plt.annotate('solid: development branch after merge compiled with mpi run as multisim', xy=(0.625, 0.96), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashed: development branch before merge', xy=(0.625, 0.93), xycoords='axes fraction', size=7.3)
	#plt.annotate('dotted: development branch after merge compiled without mpi', xy=(0.625, 0.90), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashdot: No Parallel Tempering', xy=(0.625, 0.86), xycoords='axes fraction', size=7.3)
	plt.ylabel('Probability')
	plt.xlabel(value)
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_within_a_multisim_num_molecules(dottedtpath, box, value):
	global checkIfEmpty
	floats = []
	checkIfEmpty = "%s/*/%s" % (dottedtpath, value)
	
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		parameters = stats.norm.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_histogram.png" % (dottedtpath, value)
	plt.title(dottedtpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_individual_replicas_num_molecules(dottedtpath, box, value):
	global checkIfEmpty
	floats = []
	checkIfEmpty = "%s/%s" % (dottedtpath, value)

	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		parameters = stats.norm.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	
	for name in glob.glob(checkIfEmpty):
		tokens = name.split(os.sep)
		tokens.remove(tokens[len(tokens)-1])
		newName = os.sep.join(tokens)
		filename = "%s/%s_histogram_exchange_vs_nonexchange.png" % (newName, value)
		#print filename	
		plt.title(dottedtpath) 
		plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_all_replicas_num_molecules(dottedtpath, box, value):
	#	print dottedtpath, box
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/*/*/n1dis5a.dat" % (dottedtpath)
	elif box == 1:
		checkIfEmpty = "%s/*/*/n1dis5a.dat" % (dottedtpath)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
#	print "Right before for loop"
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		data = np.loadtxt(name)
		#print floats[index]
		#parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
		#pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		if "nonexchange" in name:
			plt.plot(data[:,0], data[:,1], linestyle="dashed", linewidth=2)
		elif "reference" in name:
			plt.plot(data[:,0], data[:,1], linestyle="dotted", linewidth=2)
		else:
			plt.plot(data[:,0], data[:,1], linewidth=1)		
		#index = index+1
	filename = "%s/%s_BOX_%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (dottedtpath, value, box)
	plt.title(value) 
	plt.annotate('solid: MPI-multisim MPIPT branch; dashed: non-MPI MPIPT branch; dotted: development branch', xy=(0.015, 0.97), xycoords='axes fraction', size=7.3)
	plt.ylabel("Number of Samples With X Number of Molecules")
	plt.xlabel("Number of Molecules")
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def testEq6WithinAMultiSim(dottedtpath, box, value):
	global checkIfEmpty
	floats = []
	if box == 0:
		checkIfEmpty = "%s/*/%s_BOX_0.dat" % (dottedtpath, value)
	elif box == 1:
		checkIfEmpty = "%s/*/%s_BOX_1.dat" % (dottedtpath, value)
	else:
		print "Box must be either 0 or 1."
		return
	index = 0
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		#print floats[index]
		parameters = stats.norm.fit(floats[index])  
		lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), np.amax(floats[index]) - np.amin(floats[index]))
		pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
		lab =   "Gamma_"+str(index)
		plt.plot(lnspc, pdf_gamma, label=lab)
		index = index+1
	filename = "%s/%s_BOX_%s_histogram.png" % (dottedtpath, value, box)
	plt.title(dottedtpath) 
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def findOverlap(a, b, rtol = 1e-05, atol = 1e-08):
    ovr_a = []
    ovr_b = []
    start_b = 0
    for i, ai in enumerate(a):
        for j, bj in islice(enumerate(b), start_b, None):
            if np.isclose(ai, bj, rtol=rtol, atol=atol, equal_nan=False):
                ovr_a.append(i)
                ovr_b.append(j)
            elif bj > ai: # (more than tolerance)
                break # all the rest will be farther away
            else: # bj < ai (more than tolerance)
                start_b += 1 # ignore further tests of this item
    return (ovr_a, ovr_b)

def getAllWindows(L):
    for w in range(1, len(L)+1):
        for i in range(len(L)-w+1):
            yield L[i:i+w]

def find_sub_list(sl,l):
    results=[]
    sll=len(sl)
    for ind in (i for i,e in enumerate(l) if e==sl[0]):
        if l[ind:ind+sll]==sl:
            results.append((ind,ind+sll-1))

    return results

def linearFunc(x, a, b):
    return b + a*x

def exponentialFunc(x, a, b):
    return np.exp(b + a*x)

def fermiFunc(x, m, betaAvg, delBeta, a):
    #print "m = %d" % m
    #print "betaAvg*a %.10f*%f = %f" % (betaAvg, a, betaAvg*a)
    #print "delBeta*x %.10f*%f = %f" % (delBeta, x, delBeta*x)
    #print "sum = %f" % (betaAvg*a + delBeta*x)
    #print "np.exp(m + sum) = %f" % (np.exp(m+ betaAvg*a + delBeta*x))
    #print "1./1+np.exp(%d + %f*%f + %.10f*%f) = %f" % (m, betaAvg, delBeta, a, x, 1./(1.0+np.exp(m + betaAvg*a + delBeta*x)))
    return 1./float(1.0+np.exp(-(m + betaAvg*a + delBeta*x)))
    #return 1./float(1.0+np.exp((m + betaAvg*a + delBeta*x)))

def lik(x, parameters):

    delBeta = x[0]
    a = x[1]
    energyAndBetaAv = np.array(parameters)

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = energyAndBetaAv[1:int(energyAndBetaAv[0])+1]
    N2 = energyAndBetaAv[int(energyAndBetaAv[0]+1):-1]
    betaAvg = energyAndBetaAv[-1]

    sum1 = 0
    sum2 = 0
    M = np.log(float(N1.size)/float(N2.size))

    #print "M :"
    #print M

    #print "betaAvg: "
    #print betaAvg

    #print "delBeta: "
    #print delBeta

    #print "a guess: "
    #print a

    #print i
    index = 0
    for i in N1:
	# delBeta and a are our parameters to estimate, so keep them in the same sign in both calls
	#print np.log(fermiFunc(i, -M, -betaAvg, -delBeta, a))
        sum1 += np.log(fermiFunc(-i, -M, -betaAvg, delBeta, a))

    index = 0
    for i in N2:
	#print np.log(fermiFunc(i, M, betaAvg, delBeta, a))
	sum2 += np.log(fermiFunc(i, M, betaAvg, delBeta, a))

    print "a : %f" % a
    print "delBeta : %f" % delBeta
    print "E: "
    print sum1 + sum2
    return -(sum1 + sum2)

def g_theta(z, mean, var):
    return (1./np.sqrt(2*np.pi*np.sqrt(var)))*np.exp(-0.5*np.square(z-mean)/var)

