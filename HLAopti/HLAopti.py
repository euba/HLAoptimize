import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
			np.median(veldf["Immunogenicity Score"]), # calculate median of immunogenicity
			np.median(veldf["Error Rate"]), # calculate median of error rate
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
	snsplot = sns.barplot(x="HLA coverage", y="Element", palette=cols, data=veldata, ax=ax) # main plot function
	snsplot.set_title(title, size=50) # setting the plot title based on selected virus
	snsplot.set_xlabel('HLA coverage', size=20) # setting size of x-axis label
	snsplot.set_ylabel('Element', size=20) # setting size of y-axis label
	snsplot.figure.savefig(outplot) # saving the plot

