import numpy as np 
import time

# set numpy print options
np.set_printoptions(suppress=True,precision=4)

# expectation value for the one body part, Harmonic oscillator in three dimensions
def onebody(i, nn, il):
    hbar = 1
    omega = 1
    return hbar*omega*(2*nn[i] + il[i] + 1.5)

# get energy from eigenvalues
def get_energy(epsilon):
    return np.sum(epsilon)

# map nucleon number to state 
def state(alpha, iidx):
#    return int(np.floor(alpha/2.))
    return iidx[alpha]-1



#if __name__ == '__main__':

Norbitals = 40 # number of orbits to include
Nocc = 8 # number of occupied orbitals

print "Norbitals = ",Norbitals," ",Nocc," occupied."
    
#""" Read quantum numbers from file """
idx = []
n = []
l = []
j = []
mj = []
tz = []
states = 0
with open("nucleispnumbers.dat", "r") as qnumfile:
    for line in qnumfile:
        nums = line.split()
        # now only information for neutrons is being read in
        if len(nums) != 0 and int(nums[5]) == 1:
#        if len(nums) != 0:
            idx.append(int(nums[0]))
            n.append(int(nums[1]))
            l.append(int(nums[2]))
            j.append(int(nums[3]))
            mj.append(int(nums[4]))
            tz.append(int(nums[5]))
            states += 1
            
print "states: ", states
print 'length of idx: ', len(idx)
print ""


#""" Calculate one-electron integral """
one_nucleon_I = np.zeros(states)
for i in range(states):
    one_nucleon_I[i] = onebody(i, n, l)
    print idx[i],n[i], l[i], j[i], mj[i], tz[i]

#""" Read ALL two-nucleon integrals from file """
itest=0
states*=2

print "Now states is ",states
print ""

two_nucleon_I = np.zeros([states, states, states, states])
with open("twobodyIN.dat", "r") as infile:
    for line in infile:
        l = line.split()
        a = int(l[0]) - 1
        b = int(l[1]) - 1
        c = int(l[2]) - 1
        d = int(l[3]) - 1
        itest+=1
        if itest < 10:
            print a, b, c, d, float(l[4])
            itest+=1
        two_nucleon_I[a][b][c][d] = float(l[4])
                

#        """ Run HF-iterations """
hf_count = 0
E_old = 0
E_new = 1
C = np.ones([Norbitals, Norbitals]) # HF coefficients, all initialized to 1

#while hf_count < 24:
while abs(E_new - E_old) >= 1:
    print "############### Iteration %i ###############" % hf_count

    h = np.zeros([Norbitals, Norbitals])
    F = np.zeros([Norbitals, Norbitals])

    for i in range(Norbitals):
#        h[i][i] = C[i][i]*one_nucleon_I[state(i)]
        h[i][i] = C[i][i]*one_nucleon_I[i] 
        for j in range(Norbitals):

                    # F[i][j] += C[i][j]*(two_electron_I[state(i)][state(j)][state(i)][state(j)] - \
                    # 							two_electron_I[state(i)][state(j)][state(j)][state(i)])

            alpha = i
            gamma = j

            for p in range(Nocc):
                for beta in range(Norbitals):
                    for delta in range(Norbitals):
                        a, b, d, g = state(alpha,idx), state(beta,idx), state(delta,idx), state(gamma,idx)
#                        a, b, d, g = alpha, beta, delta, gamma
                        F[alpha][gamma] += C[p][beta]*C[p][delta]* \
                                        (two_nucleon_I[a][b][g][d] - two_nucleon_I[a][b][d][g])

    print "F:"
    print F
    print "h:"
    print h

    epsilon, C = np.linalg.eigh(h + F)
    E_old = E_new
    E_new = get_energy(epsilon)

    print "new eps"
    print epsilon
    print "new C"
    print C
    print 
    print "dE = ", abs(E_new - E_old)
    print "Energy = ", E_new
    print
    print

    hf_count += 1

