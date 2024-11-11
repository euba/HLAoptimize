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
	return(veldata_scaled)

def calcFitness(velscale, velgenome, velweight=[0.5,2,2,1,2]):
	fitness_mat = velscale.loc[velgenome] # select a subset of the scaled vaccine element features based on minimal set
	fitness_all = np.array(fitness_mat.mean(axis=0)) # calculate array of individual fitness functions
	len_fitness = 1-len(velgenome)/len(velscale.index) # calculate inverse of vaccine element list, which can be maximized (e.g., minimized)
	fitness_all = np.append(len_fitness, fitness_all) # add length fitness to array
	fitness = np.sum(fitness_all * velweight)/np.sum(velweight) # calculated overall fitness as weighted average of individual fitness functions
	return(fitness)





