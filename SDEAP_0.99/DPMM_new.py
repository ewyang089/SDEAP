import csv
import sys
import itertools
import numpy as np

from scipy import linalg
from sklearn import mixture

from numpy import genfromtxt

X = genfromtxt(sys.argv[1], delimiter='\t', skip_header=1)
X=X[:,1:]
X=X.T
#print X
X=X+0.001
#print X

with open( sys.argv[1]+".cluster", "w") as myfile:

	numrows = len(X)
	numcols = len(X[0])
	
		
	for i in range(numcols) :
		#m = np.amax(X[:,i])
       		#Q = 10*X[:,i]/m
       		
                m =  np.percentile(X[:,i], 50)
                Q = 10*X[:,i]/m
                
		dpmm = mixture.BayesianGaussianMixture(n_components=numrows/6,weight_concentration_prior_type='dirichlet_process')
		dpmm.fit(Q.reshape(-1, 1))
		X_label = dpmm.predict(Q.reshape(-1, 1))
		for c in X_label[:-1] :
			myfile.write("%s\t" % c)
		myfile.write("%s" % X_label[-1])        	
                myfile.write("\n")
		print X_label