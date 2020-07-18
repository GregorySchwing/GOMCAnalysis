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
import numdifftools as ndt
from lmfit import Minimizer, Parameters
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

def mleExpr():
	sy.init_printing()

	#a = sy.Symbol("a")
	#delBeta = sy.Symbol("dB")
	#m = sy.Symbol("m")
	#betaAvg = sy.Symbol("betaAvg")
	x1 = sy.Symbol("x1")
	#i = sy.Symbol("i")
	#expr = sy.Sum(2.0*x1)
	#x2 = sy.Symbol("x2")
	#expr = sy.Sum(sy.ln(1./(1.0+sy.exp(x1)))# + sy.ln(1./(1.0+sy.exp(-(delBeta*x2)))))
	expr = 2.0*x1# + sy.ln(1./(1.0+sy.exp(-(delBeta*x2))))
	#expr = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x1)))) + sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))
	print "the expression"
	print expr
	#print "the jacobian"
	#print sy.diff(expr, a, delBeta)
	#print "Jac"	
	#print list(sy.derive_by_array(expr, (a, delBeta)))
	#print "Hess"	
	#print (sy.derive_by_array(sy.derive_by_array(expr, (a, delBeta)), (a, delBeta))).tomatrix()
	#s=Sum(Indexed('x1',i),(i,1,3))
	#f = lambda x: Subs(s.doit(), [s.function.subs(s.variables[0], j)
	#for j in range(s.limits[0][1], s.limits[0][2] + 1)], x).doit()
	#print f((30,10,2))
	#x = sy.Function('x')
	#sy.Sum(x1(i),(i,1,3)).doit()
	#h = sy.lambdify(x1, expr)
	#print h(0)
	#print h(1, 1)



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

			print "Untrimmed Histograms of overlaps"
			print unnormedHist1
			print unnormedHist2


			########################################################################

			##### Let's find the bin with the smallest difference between the two histograms
			##### This will be our 'middle' bin ############################

			#middleBin = np.argmin(np.absolute(np.subtract(unnormedHist1, unnormedHist2)))
			metrics = np.sqrt(np.add(np.square(unnormedHist1), np.square(unnormedHist2)))
			print "metrics"
			print metrics

			angles = np.rad2deg(np.arctan(np.true_divide(unnormedHist1, unnormedHist2)))
			print "angles"
			print angles

			metricScale = np.absolute(np.subtract(45.0, angles))/45.0
			print "metricScale"
			print metricScale

			middleBin = np.argmax(np.multiply(metrics, metricScale))
			# Lets make this skew proof
			print "middle bin: "			
			print middleBin
			################################################################

			##### Let's digitize our inputs to make for easier trimming of data

			binIndicesOfLeftReplicaPoints = np.digitize(overlaps[index-1], bins=unnormedBin_edges1)
			binIndicesOfRightReplicaPoints = np.digitize(overlaps[index], bins=unnormedBin_edges2)

			##################################################################

			##### Lets find the center of the overlaps ###

			#kSetsHist1 = np.array(list(getAllWindows(unnormedHist1)))
			#kSetsHist2 = np.array(list(getAllWindows(unnormedHist2)))

			#z_val = []
			#y_val = []
			#x_val = []
			#minimum_index = 0
			#minimum_diff = float("inf")
			#for i in range(kSetsHist1.size):
    			#	y_val.append(np.true_divide(np.sum(np.absolute(np.subtract(kSetsHist1[i], kSetsHist2[i]))), kSetsHist1[i].size))
			#	x_val.append(abs(numBins - kSetsHist1[i].size))
			#	z_val.append(i)
			#	final_val = y_val[i]
			#	if (final_val < minimum_diff and kSetsHist1[i].size > numBins/4):
			#		minimum_diff = final_val
			#		minimum_index = i


			#area_linear = simps(y_val, x=x_val)

			#indices1 = find_sub_list(kSetsHist1[minimum_index].tolist(), unnormedHist1.tolist())
			#indices2 = find_sub_list(kSetsHist2[minimum_index].tolist(), unnormedHist2.tolist())

			#intersectionOfResultArrays = np.intersect1d(indices1, indices2)	
			#print intersectionOfResultArrays

			#lowerBound = intersectionOfResultArrays[0]
			
			#upperBound = intersectionOfResultArrays[1]

			lowerBound = middleBin - 8
			upperBound = middleBin + 8

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

			unnormedTrimmedhist1, unnormedTrimmedbin_edges1 = np.histogram(trimmedLeftReplica, bins=20)
			unnormedTrimmedhist2, unnormedTrimmedbin_edges2 = np.histogram(trimmedRightReplica,  bins=20)

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

def left_skewed_visual_inspection_plots(dottedtpath, box, value):
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
		b4 = np.fromfile(name, dtype='float64', count=-1, sep='\n')

		smallest = np.min(b4)

		scaledb4 = np.add(1-smallest, b4)

		largest = np.max(scaledb4)
		transformed = np.sqrt(np.subtract(largest, b4))
		floats.append(transformed)
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

			print "Untrimmed Histograms of overlaps"
			print unnormedHist1
			print unnormedHist2


			########################################################################

			##### Let's find the bin with the smallest difference between the two histograms
			##### This will be our 'middle' bin ############################

			middleBin = np.argmax(np.multiply(np.sqrtnp.add(np.square(unnormedHist1), np.square(unnormedHist2)),(1-np.arctan(unnormedHist1, unnormedHist2))))
			# Lets make this skew proof
			print "middle bin: "			
			print middleBin

			################################################################

			##### Let's digitize our inputs to make for easier trimming of data

			binIndicesOfLeftReplicaPoints = np.digitize(overlaps[index-1], bins=unnormedBin_edges1)
			binIndicesOfRightReplicaPoints = np.digitize(overlaps[index], bins=unnormedBin_edges2)

			##################################################################

			##### Lets find the center of the overlaps ###

			#kSetsHist1 = np.array(list(getAllWindows(unnormedHist1)))
			#kSetsHist2 = np.array(list(getAllWindows(unnormedHist2)))

			#z_val = []
			#y_val = []
			#x_val = []
			#minimum_index = 0
			#minimum_diff = float("inf")
			#for i in range(kSetsHist1.size):
    			#	y_val.append(np.true_divide(np.sum(np.absolute(np.subtract(kSetsHist1[i], kSetsHist2[i]))), kSetsHist1[i].size))
			#	x_val.append(abs(numBins - kSetsHist1[i].size))
			#	z_val.append(i)
			#	final_val = y_val[i]
			#	if (final_val < minimum_diff and kSetsHist1[i].size > numBins/4):
			#		minimum_diff = final_val
			#		minimum_index = i


			#area_linear = simps(y_val, x=x_val)

			#indices1 = find_sub_list(kSetsHist1[minimum_index].tolist(), unnormedHist1.tolist())
			#indices2 = find_sub_list(kSetsHist2[minimum_index].tolist(), unnormedHist2.tolist())

			#intersectionOfResultArrays = np.intersect1d(indices1, indices2)	
			#print intersectionOfResultArrays

			#lowerBound = intersectionOfResultArrays[0]
			
			#upperBound = intersectionOfResultArrays[1]

			lowerBound = middleBin - 5
			upperBound = middleBin + 5

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
		floats.append(np.fromfile(name, dtype=np.float64, count=-1, sep='\n'))
		
	floats = np.float64(floats)
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
			betas.append(1./(float(temperatures[i][0])))
			betas.append(1./(float(temperatures[next_i][0])))

			min_i = np.min(floats[i])
			min_next_i = np.min(floats[next_i])

			smallest = min(min_i, min_next_i)
			print smallest

			betaAverage = 0.5 * (1./(float(temperatures[i][0])) + 1./(float(temperatures[next_i][0])))
			deltaBeta = -(1./(float(temperatures[next_i][0])) - 1./(float(temperatures[i][0])))
			#deltaBeta = 0
			#print np.min(np.subtract(floats[i], smallest))
			#print np.min(np.subtract(floats[next_i], smallest))

			x_in = np.array([floats[i].size])			
			x_in_int = np.concatenate((x_in, np.subtract(floats[i], smallest)), axis = None)
			x_in_final = np.concatenate((x_in_int, np.subtract(floats[next_i], smallest)), axis = None)
			x_in_final_w_betaAv = np.concatenate((x_in_final, betaAverage), axis = None)
			x_in_final_w_betaAv = np.array(x_in_final_w_betaAv)

			aGuess = -1000
			guessStart = np.array([deltaBeta-0.005, aGuess])
			guessInt = np.array(guessStart)

			lowerBounds = guessInt.tolist()
			upperBounds = guessInt.tolist()

			lowerBounds[1] = .1
			upperBounds[1] = 'inf'	
			bnds = opt.Bounds(lowerBounds, upperBounds) 
			print bnds
			#print bnds.shape
			#print guess.size
			#Hessian(guessStart, x_in_final_w_betaAv)
			#lik_model = opt.minimize(lik, x0 = guessStart, args = x_in_final_w_betaAv, jac=Jacobian, hess=Hessian, method='COBYLA', options={'xtol': 1e-8, 'disp': True})
			#print opt.fmin_cg(f = lik, x0 = guessStart, args = (x_in_final_w_betaAv,), fprime=Jacobian, full_output = True)
			#print "likelihood out"			
			#print lik_model
			#print opt.check_grad(symlik, Jacobian, guessStart, x_in_final_w_betaAv)
			print Jacobian(guessStart, x_in_final_w_betaAv)

			#hessian_ndt, info = ndt.Hessian(lik, full_output=True)(guessStart, x_in_final_w_betaAv)
			#print "Hessian"
			#print hessian_ndt
			#print "info"
			#print info
			#hessian_ndt, info = Hfun(lik_model)

			#print "inverse"
			#print np.linalg.inv(hessian_ndt)

			#se = np.sqrt(np.diag(np.linalg.inv(hessian_ndt)))

			#print "se"
			#print se
			#results = pd.DataFrame({'parameters':res['x'],'std err':se})
			#results.index=['Intercept','Slope','Sigma']   
			#results.head() 

			#plt.scatter(x,y)
			#print lik_model

			#plt.show()
			#guess = [deltaBeta, 1]
			#popt, pcov = opt.curve_fit(f = linearFunc, xdata = floats[i], ydata = floats[i]*deltaBeta, p0 = guess)

			# Linear Least Squares Fit
			#plt.errorbar(
    			#x = floats[i],
    			#y = floats[i]*deltaBeta,
			#yerr = 2*np.sqrt(np.var(floats[i]*deltaBeta)),
			#color='blue'
			#)
			# -(B2-B1)E
			#plt.errorbar(
    			#x = floats[i],
    			#y = lik_model.x[0]*floats[i],
			#yerr = 2*se,
			#color='green'
			#)

			#print "delBeta true: %f" % deltaBeta
			#print "delBeta max like : %f" % lik_model.x[0]

			#plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
			#ylabe = r'$\mathrm{ln}[\mathrm{P}_%d(E)/\mathrm{P}_%d(E)$]' % (index, index-1)
			#special = "LinearMaxLikeReplica_{}_AndReplica_{}_linear".format(index, index-1)
			#filename = "%s/%s_BOX_%s_%s.png" % (dottedtpath, value, box, special)
			#plt.tick_params(axis='both', which='major', labelsize=16)
			#plt.xlabel(r'E($k_B$T)', fontsize=18)
			#plt.ylabel(ylabe, fontsize=18)
			#plt.title("Linear Visual Inspection : Replica %d vs Replica %d" % (index, index-1), fontsize=20, pad=15)
			#plt.savefig(filename, bbox_inches='tight', pad_inches=0.5)
			#plt.clf()
			#plt.cla()

			index = index+1
		except StopIteration:
			print "Reached last sim"

def left_skewed_combined_max_likelihood(dottedtpath, box, value):
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
	preTransformed = []
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		temperatures.append(re.findall("_(\d*[.,]?\d*)/", name))
		b4 = np.fromfile(name, dtype='float64', count=-1, sep='\n')

	maxOfAll = 0.0
	for array in b4:
    		max_N1 = np.amax(array)
		if (max_N1 > maxOfAll):
			maxOfAll = max_N1
		
	for array in b4:
    		transformed = np.sqrt(np.subtract(maxOfAll, array))

		mean = np.mean(transformed)
		standard_deviation = np.std(transformed)
		distance_from_mean = abs(transformed - mean)
		max_deviations = 2
		not_outlier = distance_from_mean < max_deviations * standard_deviation
		no_outliers = transformed[not_outlier]
		floats.append(transformed)

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


			aGuess = -1000
			#guessStart = [-9000, aGuess]
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
			lik_model = opt.minimize(symlik, x0 = guessStart, args = x_in_final_w_betaAv,method='l-bfgs-b')
			#plt.scatter(x,y)
			print lik_model

			hessian_ndt, info = ndt.Hessian(lik, full_output=True)(guessStart, x_in_final_w_betaAv)
			print "Hessian"
			print hessian_ndt
			print "info"
			print info
			#hessian_ndt, info = Hfun(lik_model)

			print "inverse"
			print np.linalg.inv(hessian_ndt)

			se = np.sqrt(np.diag(np.linalg.inv(hessian_ndt)))

			print "se"
			print se

			#plt.show()
			#guess = [deltaBeta, 1]
			#popt, pcov = opt.curve_fit(f = linearFunc, xdata = floats[i], ydata = floats[i]*deltaBeta, p0 = guess)

			# Linear Least Squares Fit
			plt.errorbar(
    			x = floats[i],
    			y = floats[i]*deltaBeta,
			#yerr = 2*np.sqrt(np.var(floats[i]*deltaBeta)),
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


def log_transformed_combined_max_likelihood(dottedtpath, box, value):
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
	preTransformed = []
	for name in glob.glob(checkIfEmpty):
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		temperatures.append(re.findall("_(\d*[.,]?\d*)/", name))
		b4 = np.fromfile(name, dtype='float64', count=-1, sep='\n')

	minOfAll = 0.0
	for array in b4:
    		min_N1 = np.amin(array)
		if (min_N1 < minOfAll):
			minOfAll = min_N1
		
	for array in b4:
    		transformed = np.log(np.add(array, np.absolute(1.0-minOfAll)))
		floats.append(transformed)

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


			aGuess = -1000
			#guessStart = [-9000, aGuess]
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
			lik_model = opt.minimize(lik, x0 = guessStart, args = x_in_final_w_betaAv,method='l-bfgs-b')
			#plt.scatter(x,y)
			print lik_model

			hessian_ndt, info = ndt.Hessian(lik, full_output=True)(guessStart, x_in_final_w_betaAv)
			print "Hessian"
			print hessian_ndt
			print "info"
			print info
			#hessian_ndt, info = Hfun(lik_model)

			print "inverse"
			print np.linalg.inv(hessian_ndt)

			se = np.sqrt(np.diag(np.linalg.inv(hessian_ndt)))

			print "se"
			print se

			#plt.show()
			#guess = [deltaBeta, 1]
			#popt, pcov = opt.curve_fit(f = linearFunc, xdata = floats[i], ydata = floats[i]*deltaBeta, p0 = guess)

			# Linear Least Squares Fit
			plt.errorbar(
    			x = floats[i],
    			y = floats[i]*deltaBeta,
			#yerr = 2*np.sqrt(np.var(floats[i]*deltaBeta)),
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
	fig,(ax1)=plt.subplots(1,1)
	for name in glob.glob(checkIfEmpty):
		tokens = name.split(os.sep)
		print "I think this is temp"
		lab = tokens[len(tokens)-3] + " : " + tokens[len(tokens)-2]
		print lab
		print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		floats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		unnormedHist1, unnormedBin_edges1 = np.histogram(floats[index], bins='auto')
		bincenters1 = np.array(0.5*(unnormedBin_edges1[1:]+unnormedBin_edges1[:-1]))
		#print floats[index]



		# histogram ratio in realspace
		ax1.errorbar(
    		x = bincenters1,
    		y = unnormedHist1,
		label = lab
		)
		index = index+1
	
	filename = "%s/%s_BOX_%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (dottedtpath, value, box)
	# get handles
	handles, labels = ax1.get_legend_handles_labels()
	# remove the errorbars
	handles = [h[0] for h in handles]
	# use them in the legend
	ax1.legend(handles, labels, loc='upper left',numpoints=1)
	#plt.title("NVT Butane") 
	#plt.annotate('solid: development branch after merge compiled with mpi run as multisim', xy=(0.625, 0.96), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashed: development branch before merge', xy=(0.625, 0.93), xycoords='axes fraction', size=7.3)
	#plt.annotate('dotted: development branch after merge compiled without mpi', xy=(0.625, 0.90), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashdot: No Parallel Tempering', xy=(0.625, 0.86), xycoords='axes fraction', size=7.3)

	#swapfloats = []
	#index = 0
	#checkIfEmpty = "%s/*/*/swaps.dat" % (dottedtpath)
	#for name in glob.glob(checkIfEmpty):
	#	print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
	#	swapfloats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
	#	plt.plot(swapfloats[index], np.zeros(swapfloats[index].size))
		#print floats[index]
	#	parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
	#	lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), (np.amax(floats[index]) - np.amin(floats[index]))/2)
	#	pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
	#	if "dev_multisim" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dashed", linewidth=2)
	#	elif "nonMPInonPT" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dotted", linewidth=2)
	#	else:
	#		plt.plot(lnspc, pdf_gamma, linewidth=1)		
	#	index = index+1


	plt.ylabel('Probability')
	plt.xlabel(value)
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def plot_histograms_across_two_multisims_all_replicas_left_skewed(dottedtpath, box, value):
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

	maxOfAll = 0.0
	for array in floats:
    		max_N1 = np.amax(array)
		if (max_N1 > maxOfAll):
			maxOfAll = max_N1
		
	for array in floats:
    		transformed = np.sqrt(np.subtract(maxOfAll, array))

		mean = np.mean(transformed)
		standard_deviation = np.std(transformed)
		distance_from_mean = abs(transformed - mean)
		max_deviations = 2
		not_outlier = distance_from_mean < max_deviations * standard_deviation
		no_outliers = transformed[not_outlier]


		unnormedHist1, unnormedBin_edges1 = np.histogram(no_outliers, bins='auto')
		bincenters1 = np.array(0.5*(unnormedBin_edges1[1:]+unnormedBin_edges1[:-1]))
		#print floats[index]



		# histogram ratio in realspace
		plt.errorbar(
    		x = bincenters1,
    		y = unnormedHist1
		)
		index = index+1
	
	filename = "%s/%s_BOX_%s_histograms_of_all_replicas_exchange_vs_nonexchange.png" % (dottedtpath, value, box)
	#plt.title("NVT Butane") 
	#plt.annotate('solid: development branch after merge compiled with mpi run as multisim', xy=(0.625, 0.96), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashed: development branch before merge', xy=(0.625, 0.93), xycoords='axes fraction', size=7.3)
	#plt.annotate('dotted: development branch after merge compiled without mpi', xy=(0.625, 0.90), xycoords='axes fraction', size=7.3)
	#plt.annotate('dashdot: No Parallel Tempering', xy=(0.625, 0.86), xycoords='axes fraction', size=7.3)

	#swapfloats = []
	#index = 0
	#checkIfEmpty = "%s/*/*/swaps.dat" % (dottedtpath)
	#for name in glob.glob(checkIfEmpty):
	#	print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
	#	swapfloats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
	#	plt.plot(swapfloats[index], np.zeros(swapfloats[index].size))
		#print floats[index]
	#	parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
	#	lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), (np.amax(floats[index]) - np.amin(floats[index]))/2)
	#	pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
	#	if "dev_multisim" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dashed", linewidth=2)
	#	elif "nonMPInonPT" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dotted", linewidth=2)
	#	else:
	#		plt.plot(lnspc, pdf_gamma, linewidth=1)		
	#	index = index+1


	plt.ylabel('Probability')
	plt.xlabel(value)
	plt.savefig(filename, bbox_inches='tight')
	plt.clf()
	plt.cla()

def fit_data_to_norm_dist_across_two_multisims_all_replicas(dottedtpath, box, value):
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

	swapfloats = []
	index = 0
	checkIfEmpty = "%s/*/*/swaps.dat" % (dottedtpath)
	for name in glob.glob(checkIfEmpty):
	#	print "For BOX_%s %s histograms, plotting %s" % (box, value, name)
		#print temp
		swapfloats.append(np.fromfile(name, dtype=float, count=-1, sep='\n'))
		plt.plot(swapfloats[index], np.zeros(swapfloats[index].size))
		#print floats[index]
	#	parameters = stats.norm.fit(floats[index])  
		#lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), floats[index].size)
	#	lnspc = np.linspace(np.amin(floats[index]), np.amax(floats[index]), (np.amax(floats[index]) - np.amin(floats[index]))/2)
	#	pdf_gamma = stats.norm.pdf(lnspc,loc = parameters[0],scale = parameters[1])
	#	if "dev_multisim" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dashed", linewidth=2)
	#	elif "nonMPInonPT" in name:
	#		plt.plot(lnspc, pdf_gamma, linestyle="dotted", linewidth=2)
	#	else:
	#		plt.plot(lnspc, pdf_gamma, linewidth=1)		
		index = index+1


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
    return 1./float(1.0+np.exp(-(m + betaAvg*a + delBeta*x)))

def symlik(x, parameters):

    energyAndBetaAv = np.array(parameters)

    a = sy.Symbol("a")
    delBeta = sy.Symbol("dB")
    m = sy.Symbol("m")
    betaAvg = sy.Symbol("betaAvg")
    x1 = sy.Symbol("x1")
    x2 = sy.Symbol("x2")
    expr1 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x1)))) 
    expr2 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))
    expr3 = sy.ln(1./(1.0+sy.exp(-(-m - betaAvg*a - delBeta*x1)))) + sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))

    N1_mle = sy.lambdify((x1, m, betaAvg, delBeta, a), expr1, modules=['numpy', 'sympy'])
    N2_mle = sy.lambdify((x2, m, betaAvg, delBeta, a), expr2, modules=['numpy', 'sympy'])

    Ncomb = sy.lambdify((x1, x2, m, betaAvg, delBeta, a), expr3, modules=['numpy', 'sympy'])

    #print NSimple

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = np.float64(energyAndBetaAv[1:int(energyAndBetaAv[0])+1])
    N2 = np.float64(energyAndBetaAv[int(energyAndBetaAv[0]+1):-1])
    delBetaArg = np.float64(x[0])
    aArg = np.float64(x[1])  
    betaAvgArg = np.float64(energyAndBetaAv[-1])

    sum1 = np.float64(0.0)
    sum2 = np.float64(0.0)

    #max_N1 = np.amax(N1)
    #max_N2 = np.amax(N2)

    #largest = max(max_N1, max_N2)

   # np.sqrt(np.subtract(largest, N1))
    #np.sqrt(np.subtract(largest, N2))

    M = np.float64(np.log(float(N1.size)/float(N2.size)))


    for i in range(N1.size):
	# delBeta and a are our parameters to estimate, so keep them in the same sign in both calls
	#print np.log(fermiFunc(i, -M, -betaAvg, -delBeta, a))
        #sum1 += np.log(fermiFunc(-i, -M, -betaAvg, delBeta, a))
	#sum1 += Ncomb(np.float64(-N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)
	sum1 += Ncomb(np.float64(-N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)


    print "a : %f" % aArg
    print "delBeta : %f" % delBetaArg
    print "E: "
    print sum1 + sum2
    return sum1

def lik(x, parameters):

    delBeta = x[0]
    a = x[1]
    energyAndBetaAv = np.array(parameters)

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = energyAndBetaAv[1:int(energyAndBetaAv[0])+1]
    N2 = energyAndBetaAv[int(energyAndBetaAv[0]+1):-1]
    betaAvg = energyAndBetaAv[-1]

    sum1 = 0.0
    sum2 = 0.0

    M = np.log(float(N1.size)/float(N2.size))

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

def Hessian(x, parameters):

    energyAndBetaAv = np.array(parameters)

    a = sy.Symbol("a")
    delBeta = sy.Symbol("dB")
    m = sy.Symbol("m")
    betaAvg = sy.Symbol("betaAvg")
    x1 = sy.Symbol("x1")
    x2 = sy.Symbol("x2")
    expr1 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x1)))) 
    expr2 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))
    expr3 = sy.ln(1./(1.0+sy.exp(-(-m - betaAvg*a - delBeta*x1)))) + sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))

    N1_mle = sy.lambdify((x1, m, betaAvg, delBeta, a), expr1, modules=['numpy', 'sympy'])
    N2_mle = sy.lambdify((x2, m, betaAvg, delBeta, a), expr2, modules=['numpy', 'sympy'])

    h = sy.derive_by_array(sy.derive_by_array(expr3, (a, delBeta)), (a, delBeta))

    Ncomb = sy.lambdify((x1, x2, m, betaAvg, delBeta, a), h, modules=['numpy', 'sympy'])

    print Ncomb

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = np.float64(energyAndBetaAv[1:int(energyAndBetaAv[0])+1])
    N2 = np.float64(energyAndBetaAv[int(energyAndBetaAv[0]+1):-1])
    delBetaArg = np.float64(x[0])
    aArg = np.float64(x[1])  
    betaAvgArg = np.float64(energyAndBetaAv[-1])


    sumHes = np.float64(np.zeros((2, 2)))

    #max_N1 = np.amax(N1)
    #max_N2 = np.amax(N2)

    #largest = max(max_N1, max_N2)

   # np.sqrt(np.subtract(largest, N1))
    #np.sqrt(np.subtract(largest, N2))

    M = np.float64(np.log(float(N1.size)/float(N2.size)))


    for i in range(N1.size):
	# delBeta and a are our parameters to estimate, so keep them in the same sign in both calls
	#print np.log(fermiFunc(i, -M, -betaAvg, -delBeta, a))
        #sum1 += np.log(fermiFunc(-i, -M, -betaAvg, delBeta, a))
	#sum1 += Ncomb(np.float64(-N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)
	sumHes+= Ncomb(np.float64(N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)

    print "a : %f" % aArg
    print "delBeta : %f" % delBetaArg
    print "E: "
    print sumHes
    return sumHes

def Jacobian(x, parameters):

    energyAndBetaAv = np.array(parameters)

    a = sy.Symbol("a")
    delBeta = sy.Symbol("dB")
    m = sy.Symbol("m")
    betaAvg = sy.Symbol("betaAvg")
    x1 = sy.Symbol("x1")
    x2 = sy.Symbol("x2")
    expr1 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x1)))) 
    expr2 = sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))
    expr3 = sy.ln(1./(1.0+sy.exp(-(-m - betaAvg*a - delBeta*x1)))) + sy.ln(1./(1.0+sy.exp(-(m + betaAvg*a + delBeta*x2))))

    N1_mle = sy.lambdify((x1, m, betaAvg, delBeta, a), expr1, modules=['numpy', 'sympy'])
    N2_mle = sy.lambdify((x2, m, betaAvg, delBeta, a), expr2, modules=['numpy', 'sympy'])

    j = sy.derive_by_array(expr3, (a, delBeta))
    print "Jacobian"
    print j

    Ncomb = sy.lambdify((x1, x2, m, betaAvg, delBeta, a), j, modules=['numpy', 'sympy'])

    print Ncomb

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = np.float64(energyAndBetaAv[1:int(energyAndBetaAv[0])+1])
    N2 = np.float64(energyAndBetaAv[int(energyAndBetaAv[0]+1):-1])
    delBetaArg = np.float64(x[0])
    aArg = np.float64(x[1])  
    betaAvgArg = np.float64(energyAndBetaAv[-1])


    sumJac = np.zeros(2)

    #max_N1 = np.amax(N1)
    #max_N2 = np.amax(N2)

    #largest = max(max_N1, max_N2)

   # np.sqrt(np.subtract(largest, N1))
    #np.sqrt(np.subtract(largest, N2))

    M = np.float64(np.log(float(N1.size)/float(N2.size)))


    for i in range(N1.size):
	# delBeta and a are our parameters to estimate, so keep them in the same sign in both calls
	#print np.log(fermiFunc(i, -M, -betaAvg, -delBeta, a))
        #sum1 += np.log(fermiFunc(-i, -M, -betaAvg, delBeta, a))
	#sum1 += Ncomb(np.float64(-N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)
	sumJac+= Ncomb(np.float64(N1[i]), np.float64(N2[i]), M, betaAvgArg, delBetaArg, aArg)

    print "a : %f" % aArg
    print "delBeta : %f" % delBetaArg
    print "E: "
    print sumJac
    return sumJac

def prodlik(x, parameters):

    delBeta = x[0]
    a = x[1]
    energyAndBetaAv = np.array(parameters)

    # x array is setup like this:
    # N1 size, N1, N2
    N1 = energyAndBetaAv[1:int(energyAndBetaAv[0])+1]
    N2 = energyAndBetaAv[int(energyAndBetaAv[0]+1):-1]
    betaAvg = energyAndBetaAv[-1]

    sum1 = 1.0
    sum2 = 1.0
    M = np.log(float(N1.size)/float(N2.size))

    index = 0
    for i in N1:
	# delBeta and a are our parameters to estimate, so keep them in the same sign in both calls
	#print np.log(fermiFunc(i, -M, -betaAvg, -delBeta, a))
        sum1 *= fermiFunc(-i, -M, -betaAvg, delBeta, a)

    index = 0
    for i in N2:
	#print np.log(fermiFunc(i, M, betaAvg, delBeta, a))
	sum2 *= fermiFunc(i, M, betaAvg, delBeta, a)

    print "a : %f" % a
    print "delBeta : %f" % delBeta
    print "E: "
    print sum1 * sum2
    return -(sum1 * sum2)


def g_theta(z, mean, var):
    return (1./np.sqrt(2*np.pi*np.sqrt(var)))*np.exp(-0.5*np.square(z-mean)/var)

