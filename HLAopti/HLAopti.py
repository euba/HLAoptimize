import math
import copy
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

#################################################################################
############# Functions relevant to Minimal analysis
#################################################################################

def convertVeldata(csvfile, virus, immunocut):
	### cast inputs to respective data types
	csvfile = str(csvfile)
	virus = str(virus)
	immunocut = float(immunocut)
    ### subselect the data of interest
	hladf = pd.read_csv(csvfile) # import the csv file
	virusdf = hladf[hladf["Virus Strain"] == virus] # subselect the virus of interest
	virusdf_cov = virusdf[virusdf["Immunogenicity Score"] > immunocut] # subselect the HLA that are covered
	### calculate important metrics per vaccine element
	vels = virusdf_cov['Vaccine Element'].unique() # get vaccine elements to iterate through them
	veldata = [] # create empty data structure to fill in
	for vel in vels:
		veldf = virusdf_cov[virusdf_cov['Vaccine Element'] == vel] # select the right element
		row = [vel,
			round(np.median(veldf["Immunogenicity Score"]), 2), # calculate median of immunogenicity
			round(np.median(veldf["Error Rate"]), 2), # calculate median of error rate
			len(veldf.index)] # calcute number of covered HLAs
		veldata.append(row)
	veldata = pd.DataFrame(veldata,
			columns=['Element','Immunogenicity','Error','HLA coverage']) # convert to data frame
	veldata = veldata.sort_values(by=['HLA coverage','Immunogenicity'], 
		ignore_index=True, ascending=False) # sort values by HLA coverage and Immunogenicity
	return(veldata)

def scoreBars(veldata, outplot, title):
	veldata_cnts = veldata["HLA coverage"] # taking out the HLA coverage for downstream calculation
	### canvas parameters
	fig, ax = plt.subplots(figsize=(18, len(veldata.index)/4)) # scaling the height of plot by number of vaccine elements
	cols = ['grey' if (x < max(veldata_cnts)) else 'red' for x in veldata_cnts] # color the bars based on the top picks
	### setting up the barplot
	snsplot = sns.barplot(x="HLA coverage", y="Element", 
			palette=cols, data=veldata, ax=ax) # main plot function
	snsplot.set_title(title, size=50) # setting the plot title based on selected virus
	snsplot.set_xlabel('HLA coverage', size=20) # setting size of x-axis label
	snsplot.set_ylabel('Element', size=20) # setting size of y-axis label
	snsplot.figure.savefig(outplot) # saving the plot

def score2D(veldata, outplot, title):
	veldata_cnts = veldata["HLA coverage"] # taking out the HLA coverage for downstream calculation
	### canvas parameters
	fig, ax = plt.subplots(figsize=(10, 10)) # scaling the height/width of plot
	veldata["Pick"] = ['other' if (x < max(veldata_cnts)) else 'top' for x in veldata_cnts] # color the points based on the top picks
	### setting up the scatterplot
	scat2D = sns.scatterplot(x="HLA coverage", y="Immunogenicity", size="Error",
			hue="Pick", palette=['red','grey'], 
			sizes=(20, 200), data=veldata, ax=ax) # main plotting function
	scat2D.set_title(title, size=18) # setting the plot title based on selected virus
	scat2D.set_xlabel('HLA coverage', size=15) # setting size of x-axis label
	scat2D.set_ylabel('Immunogenicity (Median)', size=15) # setting size of y-axis label
	sns.move_legend(scat2D, "upper left", bbox_to_anchor=(1, 1)) # shifting the legend to the left
	scat2D.figure.savefig(outplot) # saving the plot

#################################################################################
############# Functions relevant to Complementary analysis
#################################################################################

def compareHLAdata(csvfile, immunocut):
	### cast inputs to respective data types
	csvfile = str(csvfile)
	immunocut = float(immunocut)
	### subselect the data of interest based on the vaccine elements
	hladf = pd.read_csv(csvfile) # import the csv file
	hladf_cov = hladf[hladf["Immunogenicity Score"] > immunocut] # subselect the HLA that are covered
	### calculate important metrics per vaccine element
	vels = hladf_cov['Vaccine Element'].unique() # get vaccine elements to iterate through them
	veldata = [] # create empty data structure to fill in
	for vel in vels:
		veldf = hladf_cov[hladf_cov['Vaccine Element'] == vel] # select the right element
		row = [vel,
			# calculate median of immunogenicity over all viruses:
			np.sum(veldf["Immunogenicity Score"]),
			# calculate median of binding over all viruses:
			np.sum(veldf["Binding"]),
			# calculate median of error rate over all viruses:
			np.sum(veldf["Error Rate"]), 
			# calcute number of covered HLAs:
			len(veldf.index)]
		veldata.append(row)
	veldata = pd.DataFrame(veldata,
			columns=['Element','Immunogenicity','Binding','Error','HLA coverage']) # convert to data frame
	### scale the data by min and max to make values of different vaccine elements comparable
	scaler = MinMaxScaler() # use min max scaler from scikit package
	veldata_scaled = pd.DataFrame(scaler.fit_transform(veldata.drop('Element', axis=1)),
		columns=['Immunogenicity','Binding','Error','HLA coverage']) # scale only numeric elements
	veldata_scaled['Error'] = 1-veldata_scaled['Error'] # inverse error rate to allow maximization
	veldata_scaled.index = veldata['Element'] # rename index by vaccine elements
	return([veldata, veldata_scaled])

def calcFitness(velscale, velpheno, velweight):
	fitness_mat = velscale.loc[velpheno] # select a subset of the scaled vaccine element features based on minimal set
	fitness_all = np.array(fitness_mat.mean(axis=0)) # calculate array of individual fitness functions
	len_fitness = 1-len(velpheno)/len(velscale.index) # calculate inverse of vaccine element list, which can be maximized (e.g., minimized)
	fitness_all = np.append(len_fitness, fitness_all) # add length fitness to array
	fitness = np.sum(fitness_all * velweight)/np.sum(velweight) # calculated overall fitness as weighted average of individual fitness functions
	return(fitness)

def geno2Pheno(velindex, velgenome): # helper function to convert array of presence/absence (0/1) to list of vaccine elements 
	velpheno = list(np.take(velindex, np.where(velgenome == 1)[0]))
	return(velpheno)

def calcPopFitness(population, velscale, velweight):
	popfit = []
	for i in range(len(population)): # calculate fitness value per individual solution in population
		popfit.append(calcFitness(velscale, 
			geno2Pheno(velscale.index, population[i]),
			velweight))
	return(popfit)

def mutation(velgenome, mutrate): # helper function to mutate random positions in genome
	for i in range(1,mutrate): # mutation rate is implemented as number of positions mutated
		randind = random.randint(0, len(velgenome)-1) # select a random position
		if velgenome[randind] == 1: # invert the boolean integer of the position
			velgenome[randind] = 0
		else:
			velgenome[randind] = 1
	if np.sum(velgenome) == 0: # in case mutation creates an empty solution set
		velgenome[random.randint(0, len(velgenome)-1)] = 1
	return(velgenome)

def crossOver(velgenome1, velgenome2): # helper function to cross over two genomes with each other
	randind = random.randint(0, len(velgenome1)-1) # select random position to cross over in
	velgenome_new = np.concatenate((velgenome1[:randind], velgenome2[randind:])) # cross over both genomes
	if np.sum(velgenome_new) == 0: # in case cross over creates an empty solution set
		velgenome_new[random.randint(0, len(velgenome_new)-1)] = 1
	return(velgenome_new)

def initPop(positions, popsize): # helper function to initialize a set of random solutions
	popinit = np.ones((popsize, positions)).astype(int) # create a initial matrix of integers
	for i in range(0,popsize-1): # fill each individual in solution population with 0
		randpos = random.sample(range(0, positions-1),
			positions-math.ceil(positions/2)) # leave at least half vaccine elements in the solution
		popinit[i][randpos] = 0
	return(popinit)

def selectPop(population, velscale, mutrate, velweight): # helper function to select high fitness individuals
	popfit = calcPopFitness(population, velscale, velweight) # calculate population wide fitness per individual solution
	fitind = copy.deepcopy(population[np.argmax(popfit)]) # find solution with highest fitness
	for i in range(math.ceil(len(population)/3)): # replace low fitness solutions (one third of population) with highest fitness
		deadind = np.argmin(popfit) # find the lowest fitness individual
		population[deadind] = fitind # replace low with highest fitness individual
		popfit[deadind] = max(popfit) # replace fitness value
	for i in range(len(population)): # go trough each individual solution and perform mutation
		population[i] = mutation(population[i], mutrate)
	population[np.argmax(popfit)] = crossOver(population[np.argmax(popfit)], 
		population[random.randint(0, len(population)-1)]) # perform cross over with highest fitness and one random individual
	population[0] = fitind # keep at least one individual with original (unmutated) solution
	return(population)

def geneticAlgo(velscale, popsize, mutrate, generations, randfitness=False): # main loop of genetic algorithm
	pop = initPop(len(velscale.index), popsize) # initialize population with random solutions
	fitlist = [] # create a list with max fitness value per generation
	velweight = [0.5,1,1,1,1] # use this as default weights for objective function
	for x in range(generations): # interate through the specified number of generations
		if randfitness: # flag determines if weighting of objective function should be randomized
			velweight = [random.uniform(0.1, 2),random.uniform(0.1, 2),random.uniform(0.1, 2),random.uniform(0.1, 2),random.uniform(0.1, 2)]
		pop = selectPop(pop, velscale, mutrate, velweight) # perform selection process
		fitlist.append(np.max(calcPopFitness(pop, velscale, velweight))) # fill in the fitness list
		print(fitlist[x])
	best_sol = geno2Pheno(velscale.index,pop[np.argmax(calcPopFitness(pop, velscale, velweight))]) # select the best solution at the end of the simulation
	return([best_sol, fitlist]) # export best solution and list of max fitness values

def plotFitness(fitdat, outplot):
	fitdf = pd.DataFrame(fitdat) # convert to list of lists to data frame
	fitdf = fitdf.transpose() # transpose the data
	fitdf['Generations'] = range(1, len(fitdf.index)+1) # set the number of generations
	fitdfl = pd.melt(fitdf, ["Generations"]) # restructure data for plotting
	fitdfl = fitdfl.rename(columns={'variable': 'replicate'}) # rename variable column
	fig, ax = plt.subplots(figsize=(12, 10)) # set up plot
	lp = sns.lineplot(data=fitdfl, x='Generations', y='value', hue='replicate', ax=ax) # plot line plot
	lp.set_title("GA of 10 replicates", size=18) # setting the plot title
	lp.set_xlabel('Generations', size=15) # setting size of x-axis label
	lp.set_ylabel('Max Fitness', size=15) # setting size of y-axis label
	sns.move_legend(lp, "upper left", bbox_to_anchor=(1, 1)) # shifting the legend to the left
	lp.figure.savefig(outplot) # save the plot

