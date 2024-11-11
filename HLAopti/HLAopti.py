import numpy as np
import pandas as pd

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

def some_func():
    return 42

def some_func2(input, output):
    return ("Input: "+input+", Output: "+output)