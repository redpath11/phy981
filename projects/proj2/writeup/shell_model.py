import numpy as np 
from decimal import Decimal
import os
from scipy.misc import comb
from numpy import linalg as LA
import itertools

# bitCount() counts the number of bits set (not an optimal function)

def bitCount(int_type):
    """ Count bits set in integer """
    count = 0
    while(int_type):
        int_type &= int_type - 1
        count += 1
    return(count)

# allperms() generates a .jpg file displaying a composite image of all
# the WORST hairstyles from the '80s
#
# actually it generates all possible <Nbits> words as a list of strings

def allperms(Nbits):
	result = ["".join(seq) for seq in itertools.product("01", repeat=Nbits)]
	return result



def slaterdet(numLevels, numPairs):
#	file = open("slater_dets.list", "w")
#	file.write("This is my output! You have " + str(numParticles) + " particles! \n")
	# This is where you can do your calculation. You can nest the file.write line inside a loop, if you need to
#	nLvls=numStates/2
#	nPairs=numParticles/2

	numSDs=int(round(comb(numLevels, numPairs, exact=False)))
	SD = np.zeros([numSDs,numLevels])

	ap=allperms(numLevels)
	apmax=len(ap)
	# loop over generated list of <nLvls>-bit words and select ones with <nPairs> bits set
	# i.e. eliminate states that don't have the specified number of pairs
	idxSD=0
	idxLVL=0
	for i in range(apmax):
		state=int(ap[i],2)
		# if word has right number of pairs
		# we add a new SD
		if (bitCount(state) == numPairs):
			idxLVL=0
			for j in range(numLevels):
				SD[idxSD][idxLVL]=ap[i][j]
				idxLVL+=1
			idxSD+=1



#	file.close()
	return SD

def hamiltonian(numLevels, numPairs, dim_H, G,dE):
	# The following lines should give you a NumPy array, where each row corresponds to a Slater determinant
#	SDarray = np.genfromtxt('slater_dets.list-sample1')

	SDarray = slaterdet(numLevels, numPairs)
	b = slaterdet(numLevels, numPairs)
	print 'My array of Slater determinants is \n', SDarray

	# This next line should let you check that the data was imported correctly
	print 'It contains', np.size(SDarray,0), 'Slater determinants constructed from', np.size(SDarray,1), 'single particle states.'

	H = np.zeros([dim_H, dim_H])
	ne=0
	for i in range(dim_H):
            for q in range(numLevels):
                if SDarray[i][q]==1:
                    for p in range(numLevels):
                        if SDarray[i][p]==0:
                            b[i][q]=0
                            b[i][p]=1
                            for j in range(dim_H):
                                if np.array_equal(SDarray[j],b[i]):
#				    print "SDarray,b: ", SDarray[j],b[i]
                                    H[i][j]=H[i][j]-G
                            b[i][p]=0
                            b[i][q]=1
        print(H*G)

        for i in range(dim_H):
            for j in range(numLevels):
                if SDarray[i][j]==1:
                    ne=ne+(numLevels-j-1)
            H[i][i]=H[i][i]+2*ne*dE-G*numPairs
            ne=0
        print H
	
	return H

if __name__ == '__main__':
	Nlevels = int(raw_input('Number of Levels = '))
#	Nstates = 2*Nlevels
	Npairs = int(raw_input('Number of Pairs = '))
	size_H = comb(Nlevels, Npairs, exact=False)
#	g=1
	dE=1 #0
#	h = hamiltonian(Nlevels, Npairs, size_H, g,dE)
#	EigValues, EigVectors = LA.eigh(h)
#	permute = EigValues.argsort()
#	print '\n For g =', g, 'the eigenvalues of the Hamiltonian are \n', EigValues[permute]
#	print 'and the corresponding eigenvectors are \n', EigVectors[:,permute], '\n'
	for g in np.linspace(-1, 1, num=5):
#		slaterdet(Nstates, Nparticles)
#		h = hamiltonian(Nstates, Nparticles, size_H, g)
#		EigValues, EigVectors = LA.eigh(h)

		h = hamiltonian(Nlevels, Npairs, size_H, g,dE)
		EigValues, EigVectors = LA.eigh(h)

		permute = EigValues.argsort()
		print '\n For g =', g, 'the eigenvalues of the Hamiltonian are \n', EigValues[permute]
		print 'and the corresponding eigenvectors are \n', EigVectors[:,permute], '\n'
