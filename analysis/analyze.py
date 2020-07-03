import sys, os, glob
import numpy as np  
from simpleHist import visual_inspection_plots
from simpleHist import combined_max_likelihood
from simpleHist import single_replica_MLE
from simpleHist import plot_histograms_across_two_multisims_individual_replicas
from simpleHist import plot_histograms_across_two_multisims_all_replicas
from simpleHist import fit_data_to_norm_dist_across_two_multisims_all_replicas
from subprocess import check_output
from simpleHist import plot_histograms_within_a_multisim_num_molecules
from simpleHist import plot_histograms_across_two_multisims_individual_replicas_num_molecules
from simpleHist import plot_histograms_across_two_multisims_all_replicas_num_molecules
from simpleHist import testEq6WithinAMultiSim

def generate_energy_files(startpath):
    for root, dirs, files in os.walk(startpath):
	checkIfEmpty_Ensemble = "%s/*/*/Blk_*BOX_0.dat" % root
	if len(glob.glob(checkIfEmpty_Ensemble)) > 0:
		tokens = root.split(os.sep)
       		ensemble = tokens[len(tokens)-1]
########## ONLY IF A Blk*BOX_0.dat EXISTS ################
        checkIfEmpty_0 = "%s/Blk_*BOX_0.dat" % root
        for name in glob.glob(checkIfEmpty_0):
            print ensemble
            for value in my_dict[ensemble]:
            	# To obtain a file with only the energy values, assumed to be column 2 for now
            	systemCall1 =  "cat %s | awk \'{print $2}\' > %s/%s_BOX_0.dat" % (name, root, value)
            	os.popen(systemCall1)
            	# To remove equilibration outputs, hardcoded as 5000 lines for now
            	systemCall3 = "sed -i '1,1d' %s/%s_BOX_0.dat" % (root, value)
            	os.popen(systemCall3)
########## ONLY IF A Blk*BOX_0.dat EXISTS ################

########## ONLY IF A Blk*BOX_1.dat EXISTS ################
        checkIfEmpty_1 = "%s/Blk_*BOX_1.dat" % root
        for name in glob.glob(checkIfEmpty_1):
            #print name
            tokens = name.split(os.sep)
            for value in my_dict[ensemble]:
#	    print "creating %s" % name
            # To obtain a file with only the energy values, assumed to be column 2 for now
            	systemCall2 =  "cat %s | awk \'{print $2}\' > %s/%s_BOX_1.dat" % (name, root, value)
            	os.popen(systemCall2)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            	systemCall4 = "sed -i '1,5000d' %s/%s_BOX_1.dat" % (root, value)
            	os.popen(systemCall4)
########## ONLY IF A Blk*BOX_1.dat EXISTS ################

def generate_energy_swaps(startpath):
    for root, dirs, files in os.walk(startpath):
	checkIfEmpty_Ensemble = "%s/*/*/ConsoleOut.dat" % root
	if len(glob.glob(checkIfEmpty_Ensemble)) > 0:
	    tokens = root.split(os.sep)
       	    ensemble = tokens[len(tokens)-1]
        checkIfEmpty_0 = "%s/ConsoleOut.dat" % root
        for name in glob.glob(checkIfEmpty_0):
      	    systemCall1 =  "cat %s | grep swap | awk \'{print $2}\' > %s/swaps.dat" % (name, root)
            os.popen(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,1d' %s/swaps.dat" % (root)
            os.popen(systemCall3)



def generate_density_NPT_files(startpath):
    for root, dirs, files in os.walk(startpath):

########## ONLY IF A Blk*BOX_0.dat EXISTS ################
        checkIfEmpty_0 = "%s/Blk_*BOX_0.dat" % root
        for name in glob.glob(checkIfEmpty_0):
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall1 =  "cat %s | awk \'{print $11}\' > %s/density_BOX_0.dat" % (name, root)
            os.system(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,5000d' %s/density_BOX_0.dat" % root
            os.system(systemCall3)
########## ONLY IF A Blk*BOX_0.dat EXISTS ################

########## ONLY IF A Blk*BOX_1.dat EXISTS ################
        checkIfEmpty_1 = "%s/Blk_*BOX_1.dat" % root
        for name in glob.glob(checkIfEmpty_1):
#	    print "creating %s" % name
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall2 =  "cat %s | awk \'{print $11}\' > %s/density_BOX_1.dat" % (name, root)
            os.system(systemCall2)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall4 = "sed -i '1,5000d' %s/density_BOX_1.dat" % root
            os.system(systemCall4)
########## ONLY IF A Blk*BOX_1.dat EXISTS ################

def generate_pressure_NVT_files(startpath):
    for root, dirs, files in os.walk(startpath):

########## ONLY IF A Blk*BOX_0.dat EXISTS ################
        checkIfEmpty_0 = "%s/Blk_*BOX_0.dat" % root
        for name in glob.glob(checkIfEmpty_0):
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall1 =  "cat %s | awk \'{print $11}\' > %s/pressure_BOX_0.dat" % (name, root)
            os.system(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,5000d' %s/pressure_BOX_0.dat" % root
            os.system(systemCall3)
########## ONLY IF A Blk*BOX_0.dat EXISTS ################

########## ONLY IF A Blk*BOX_1.dat EXISTS ################
        checkIfEmpty_1 = "%s/Blk_*BOX_1.dat" % root
        for name in glob.glob(checkIfEmpty_1):
#	    print "creating %s" % name
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall2 =  "cat %s | awk \'{print $11}\' > %s/pressure_BOX_1.dat" % (name, root)
            os.system(systemCall2)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall4 = "sed -i '1,5000d' %s/pressure_BOX_1.dat" % root
            os.system(systemCall4)
########## ONLY IF A Blk*BOX_1.dat EXISTS ################

def generate_density_GEMC_files(startpath):
    for root, dirs, files in os.walk(startpath):

########## ONLY IF A Blk*BOX_0.dat EXISTS ################
        checkIfEmpty_0 = "%s/Blk_*BOX_0.dat" % root
        for name in glob.glob(checkIfEmpty_0):
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall1 =  "cat %s | awk \'{print $13}\' > %s/density_BOX_0.dat" % (name, root)
            os.system(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,5000d' %s/density_BOX_0.dat" % root
            os.system(systemCall3)
########## ONLY IF A Blk*BOX_0.dat EXISTS ################

########## ONLY IF A Blk*BOX_1.dat EXISTS ################
        checkIfEmpty_1 = "%s/Blk_*BOX_1.dat" % root
        for name in glob.glob(checkIfEmpty_1):
#	    print "creating %s" % name
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall2 =  "cat %s | awk \'{print $13}\' > %s/density_BOX_1.dat" % (name, root)
            os.system(systemCall2)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall4 = "sed -i '1,5000d' %s/density_BOX_1.dat" % root
            os.system(systemCall4)
########## ONLY IF A Blk*BOX_1.dat EXISTS ################

def generate_density_GEMC_NVT_binary_molfrac_files(startpath):
    for root, dirs, files in os.walk(startpath):

########## ONLY IF A Blk*BOX_0.dat EXISTS ################
        checkIfEmpty_0 = "%s/Blk_*BOX_0.dat" % root
        for name in glob.glob(checkIfEmpty_0):
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall1 =  "cat %s | awk \'{print $16}\' > %s/MOLFRACT_AR_BOX_0.dat" % (name, root)
            os.system(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,5000d' %s/MOLFRACT_AR_BOX_0.dat" % root
            os.system(systemCall3)

            systemCall1 =  "cat %s | awk \'{print $17}\' > %s/MOLFRACT_KR_BOX_0.dat" % (name, root)
            os.system(systemCall1)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall3 = "sed -i '1,5000d' %s/MOLFRACT_KR_BOX_0.dat" % root
            os.system(systemCall3)
########## ONLY IF A Blk*BOX_0.dat EXISTS ################

########## ONLY IF A Blk*BOX_1.dat EXISTS ################
        checkIfEmpty_1 = "%s/Blk_*BOX_1.dat" % root
        for name in glob.glob(checkIfEmpty_1):
#	    print "creating %s" % name
            # To obtain a file with only the energy values, assumed to be column 2 for now
            systemCall2 =  "cat %s | awk \'{print $17}\' > %s/MOLFRACT_KR_BOX_1.dat" % (name, root)
            os.system(systemCall2)
            # To remove equilibration outputs, hardcoded as 5000 lines for now
            systemCall4 = "sed -i '1,5000d' %s/MOLFRACT_KR_BOX_1.dat" % root
            os.system(systemCall4)
########## ONLY IF A Blk*BOX_1.dat EXISTS ################


def generate_energy_histograms_within_a_multisim(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		checkIfEmpty_2 = "%s/*/Blk_*BOX_0.dat" % dirname
		if len(glob.glob(checkIfEmpty_2)) > 0:
			tokens = dirname.split(os.sep)
			for value in my_dict[tokens[len(tokens)-2]]:
				visual_inspection_plots(dirname, 0, value)
				#combined_max_likelihood(dirname, 0, value)
				#single_replica_MLE(dirname, 0, value)		
				print dirname+" "+value+" Box 0 histogram generated"
				

		checkIfEmpty_3 = "%s/*/Blk_*BOX_1.dat" % dirname
		if len(glob.glob(checkIfEmpty_3)) > 0:
			tokens = dirname.split(os.sep)
			for value in my_dict[tokens[len(tokens)-2]]:
				visual_inspection_plots(dirname, 1, value)
				print dirname+" "+value+" Box 1 histogram generated"


def generate_energy_histograms_across_two_multisims_individual_replicas(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		########## ONLY IF A Blk*BOX_0.dat EXISTS ################
		checkIfEmpty_2 = "%s/*/*/Blk_*BOX_0.dat" % dirname
		ensembleDirectory = dirname
		filePath = glob.glob(checkIfEmpty_2)
		if len(filePath) > 0:
			for files in glob.glob(checkIfEmpty_2):
				tokens = files.split(os.sep)
				tokens[len(tokens)-3] = '*'
				tokens.remove(tokens[len(tokens)-1])
				newName = os.sep.join(tokens)
				tokensForValues = dirname.split(os.sep)
				for value in my_dict[tokensForValues[len(tokensForValues)-1]]:
					plot_histograms_across_two_multisims_individual_replicas(newName, 0, value)	

		########## ONLY IF A Blk*BOX_0.dat EXISTS ################	

		########## ONLY IF A Blk*BOX_1.dat EXISTS ################
		checkIfEmpty_2 = "%s/*/*/Blk_*BOX_1.dat" % dirname
		ensembleDirectory = dirname
		filePath = glob.glob(checkIfEmpty_2)
		if len(filePath) > 0:
			for files in glob.glob(checkIfEmpty_2):
				tokens = files.split(os.sep)
				tokens[len(tokens)-3] = '*'
				tokens.remove(tokens[len(tokens)-1])
				newName = os.sep.join(tokens)
				print newName
				tokensForValues = dirname.split(os.sep)
				for value in my_dict[tokensForValues[len(tokensForValues)-1]]:
					plot_histograms_across_two_multisims_individual_replicas(newName, 1, value)	

		########## ONLY IF A Blk*BOX_1.dat EXISTS ################


def generate_energy_histograms_across_two_multisims_all_replicas(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		########## ONLY IF A Blk*BOX_0.dat EXISTS ################
		checkIfEmpty_2 = "%s/*/*/Blk_*BOX_0.dat" % dirname
		if len(glob.glob(checkIfEmpty_2)) > 0:
			tokensForValues = dirname.split(os.sep)
			for value in my_dict[tokensForValues[len(tokensForValues)-1]]:
				plot_histograms_across_two_multisims_all_replicas(dirname, 0, value)
				print dirname+" "+value+" Box 0 histogram generated"

		########## ONLY IF A Blk*BOX_0.dat EXISTS ################	

		########## ONLY IF A Blk*BOX_1.dat EXISTS ################
		checkIfEmpty_3 = "%s/*/*/Blk_*BOX_1.dat" % dirname
		if len(glob.glob(checkIfEmpty_3)) > 0:
			tokensForValues = dirname.split(os.sep)
			for value in my_dict[tokensForValues[len(tokensForValues)-1]]:
				plot_histograms_across_two_multisims_all_replicas(dirname, 1, value)
				print dirname+" "+value+" Box 1 histogram generated"

		########## ONLY IF A Blk*BOX_1.dat EXISTS ################

def generate_num_molecules_histograms_within_a_multisim(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		checkIfEmpty_2 = "%s/*/n1dis5a.dat" % dirname
		if len(glob.glob(checkIfEmpty_2)) > 0:
			plot_histograms_within_a_multisim_num_molecules(dirname, 0, 'n1dis5a.dat')
			print dirname+" num molecules histogram generated"
				


def generate_num_molecules_histograms_across_two_multisims_individual_replicas(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		########## ONLY IF A Blk*BOX_0.dat EXISTS ################
		checkIfEmpty_2 = "%s/*/*/n1dis5a.dat" % dirname
		ensembleDirectory = dirname
		filePath = glob.glob(checkIfEmpty_2)
		if len(filePath) > 0:
			for files in glob.glob(checkIfEmpty_2):
				plot_histograms_across_two_multisims_individual_replicas_num_molecules(newName, 0, 'n1dis5a.dat')
				print dirname+" num molecules histogram generated"	

		########## ONLY IF A Blk*BOX_0.dat EXISTS ################	


def generate_num_molecules_histograms_across_two_multisims_all_replicas(startpath, my_dict):
	for dirname, dirnames, filenames in os.walk('.'):
    	# print path to all subdirectories first.
		########## ONLY IF A Blk*BOX_0.dat EXISTS ################
		checkIfEmpty_2 = "%s/*/*/n1dis5a.dat" % dirname
		if len(glob.glob(checkIfEmpty_2)) > 0:
			plot_histograms_across_two_multisims_all_replicas_num_molecules(dirname, 0, "Num Molecules")
			print dirname+" num molecules histogram generated"	

		########## ONLY IF A Blk*BOX_0.dat EXISTS ################	


print "Does ./validation exist? : \t\t", os.path.isdir("./validation")
print "Does ./validation/GCMC exist? :\t\t", os.path.isdir("./validation/GCMC")
print "Does ./validation/NPT exist? : \t\t", os.path.isdir("./validation/NPT")
print "Does ./validation/NVT exist? : \t\t", os.path.isdir("./validation/NVT")
print "Does ./validation/GEMC_NVT exist? : \t", os.path.isdir("./validation/GEMC_NVT")
print "Does ./validation/GEMC_NPT exist? : \t", os.path.isdir("./validation/GEMC_NPT")

listOfNVTObservables = ['energy']#, 'pressure']
listOfNPTObservables = ['energy', 'density']
listOfGCMCObservables = ['energy']
listOfGEMC_NVT_BinaryObservables = ['energy', 'density', 'MOLFRACT_AR', 'MOLFRACT_KR']
######### Note: need to write a generate vapor pressure method ###########
listOfGEMC_NVT_PureObservables = ['energy', 'pressure']#, 'vapor_pressure']
listOfGEMC_NPTObservables = ['energy', 'density']

my_dict = {'mu':listOfGCMCObservables, 'mu_temp':listOfGCMCObservables, 'temp':listOfGCMCObservables, 'NVT':listOfNVTObservables, 'NPT':listOfNPTObservables, 'GEMC_NPT':listOfGEMC_NPTObservables, 'binary_mixture':listOfGEMC_NVT_BinaryObservables, 'pure_fluid':listOfGEMC_NVT_PureObservables}
#print my_dict['GCMC']
#for files in my_dict['GCMC']:
#	print files



generate_energy_files(os.getcwd())
generate_energy_swaps(os.getcwd())
#os.chdir(homedir)
#os.chdir("./validation_9_17_19/NPT")
#generate_density_NPT_files(os.getcwd())
#os.chdir(homedir)
#os.chdir("./validation_9_17_19/NVT")
#generate_pressure_NVT_files(os.getcwd())
#os.chdir(homedir)
#os.chdir("./validation_9_17_19/GEMC_NVT/pure_fluid")
#generate_pressure_NVT_files(os.getcwd())
#os.chdir(homedir)
#os.chdir("./validation_9_17_19/GEMC_NVT/binary_mixture")
#generate_density_GEMC_NVT_binary_molfrac_files(os.getcwd())
#generate_density_GEMC_files(os.getcwd())
#os.chdir(homedir)
#os.chdir("./validation_9_17_19/GEMC_NPT")
#generate_density_GEMC_files(os.getcwd())
#os.chdir(homedir)

#generate_num_molecules_histograms_within_a_multisim(os.getcwd(), my_dict)
#generate_num_molecules_histograms_across_two_multisims_individual_replicas(os.getcwd(), my_dict)
#generate_num_molecules_histograms_across_two_multisims_all_replicas(os.getcwd(), my_dict)
#generate_energy_histograms_within_a_multisim(os.getcwd(), my_dict)

#generate_energy_histograms_across_two_multisims_individual_replicas(os.getcwd(), my_dict)
generate_energy_histograms_across_two_multisims_all_replicas(os.getcwd(), my_dict)
#generate_num_molecules_histograms_across_two_multisims_all_replicas(os.getcwd(), my_dict)

