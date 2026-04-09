
import numpy as np
from random import sample, seed, choices, randint
import math
import sys
from time import time, sleep
import networkx as nx
import scipy 
import os
import matplotlib.pyplot as plt
import pickle
from sklearn.neighbors import KDTree
from numba import cuda,njit

prcNr = 4000
aggNr = 300
RandomSeed = 103317
vol=0.1
L0 = 750
L0_z_factor= 0.1
Mu = 5
Sigma = 0
timeScale=0.01
totaltime=100
timeIni=0
density=Mu**3/2700
damp=0.1
VGradient=20
MixingDire=0
FixedBound=2
ovlp=0.05
BalancePos=0.5
jumpSlice=50
cutoff=20*Mu

seed(RandomSeed)
np.random.seed(RandomSeed)

perB=np.array([1,1,0])

curloc=os.path.abspath(os.curdir)


def sphereVolume(r):
    return 4/3 * math.pi * r ** 3  

def initFunc(vol):  # only under all PBC be used
    global prcNr,aggNr,PL,PS,PV,PF
    
    prcNr=L0**3*L0_z_factor*vol/sphereVolume(Mu)
    prcNr = int(prcNr)
    aggNr=prcNr
    aggNr = int(aggNr)
    print('after calculation there will be %i particles and %i agglomerates to be generated' %(prcNr,aggNr))
    rng = np.random.RandomState(RandomSeed)
    PL1=rng.uniform(0,L0,(prcNr,2))
    PL2=rng.uniform(0,L0*L0_z_factor,(prcNr,1))
    PL=np.hstack((PL1,PL2))
    PS1=abs(rng.normal(Mu,Sigma,(prcNr//10,1)))
    PS2=abs(rng.normal(Mu,Sigma,(prcNr-prcNr//10,1)))
    PS=np.vstack((PS1,PS2))
    
    #PV=rng.uniform(-L0/10/2,L0/10/2,(prcNr, 3))
    PV=rng.uniform(-L0/10/2,L0/10/2,(prcNr, 3))
    
    for i in range(prcNr):  # initial speed
        PV[i][0]+=PL[i][2]*0.5
    PF=np.zeros((prcNr,3))
    return PL,PS,PV,PF

class ProgressBar():

    def __init__(self, max_steps):
        self.max_steps = max_steps
        self.current_step = 0
        self.progress_width = 50

    def update(self, step=None):
        self.current_step = step

        num_pass = int(self.current_step * self.progress_width / self.max_steps) + 1
        num_rest = self.progress_width - num_pass 
        percent = (self.current_step+1) * 100.0 / self.max_steps 
        progress_bar = '[' + '\u25a0' * (num_pass-1) + '|' + '-' * num_rest + ']'
        progress_bar += '%.2f' % percent + '%' 
        if self.current_step < self.max_steps - 1:
            progress_bar += '\r' 
        else:
            progress_bar += '\n' 
        sys.stdout.write(progress_bar) 
        sys.stdout.flush()
        if self.current_step >= self.max_steps:
            self.current_step = 0
            print
    
def reactionForce(r,a,b):
    dist = r+a+b
    if r<0:
        return -((a+b)/(dist-0.1*(1-ovlp)*(a+b))-1)*500*(min(a,b)/Mu*2)**3-10
    elif r<cutoff:
        return -1000*(((a+b)/(dist))**6-2*((a+b)/(dist))**3)
    else:
        return 0.0
        
def xyzout(PI,PR,label):
    packet_file = 'data_set_shear/distribution_%i.xyz' % label
    os.path.join(curloc, packet_file)
    file_xyz=open(packet_file,'w')
    file_xyz.write('%d\n'% len(PI))
    file_xyz.write('Lattice="%.2f 0.0 0.0 0.0 %.2f 0.0 0.0 0.0 %.2f" Properties=species:S:1:pos:R:3:Radius:R:1\n'%(L0,L0,L0*L0_z_factor))
    for i in range(len(PI)//10):
        file_xyz.write('huge_Particle ')
        file_xyz.write(' %.4f  %.4f  %.4f  %.4f\n'%(PI[i][0],PI[i][1],PI[i][2],PR[i][0]))
        #file_xyz.write(' %.4f  %.4f  %.4f  %.4f\n'%(PI[i][0],PI[i][1],PI[i][2],Mu))
    for i in range(len(PI)//10,len(PI),1):
        file_xyz.write('small_Particle ')
        #file_xyz.write(' %.4f  %.4f  %.4f  %.4f\n'%(PI[i][0],PI[i][1],PI[i][2],Mu))
        file_xyz.write(' %.4f  %.4f  %.4f  %.4f\n'%(PI[i][0],PI[i][1],PI[i][2],PR[i][0]))
    file_xyz.close()    

def main():
    global prcNr,aggNr,timeIni
    sliceNr=0
    PL,PS,PV,PF=initFunc(vol)
    progress_bar=ProgressBar(totaltime/timeScale)
    while timeIni<totaltime:
        progress_bar.update(sliceNr)
        flag=[]
        PL__=PL.copy()
        PS__=PS.copy()  
        tree = KDTree(PL__, leaf_size=2) 
        s = pickle.dumps(tree) 
        tree_copy = pickle.loads(s)  
        PF__=np.array([0,0,0])
        PV__=np.array([0,0,0])
        #print('TimeStep',timeIni)
        for i in range(prcNr):
            dist, ind = tree_copy.query([PL__[i]], k=8)
            PF__=np.zeros(3)
            for j,r in zip(ind[0],dist[0]):
                if j>=prcNr:
                    j=flag[j-prcNr]
                if r!=0 and r<cutoff: 
                    PF__+=10*reactionForce(r-PS[j][0]-PS[i][0],PS[i][0],PS[j][0])*(PL[j]-PL[i])/r
            PF__[MixingDire]+=PL[i][FixedBound]*VGradient-BalancePos*L0*L0_z_factor*VGradient
            PV__=PV[i]+0.5*timeScale*(PF__+PF[i])/(PS[i]**3/density)
            PL[i]=PL[i]+0.5*timeScale*(PV[i]+PV__)*(1-damp)
            PV[i]=PV__*(1-damp)
            #PV[i][MixingDire]+=PL[i][FixedBound]*VGradient-BalancePos*L0*VGradient
            if PL[i][FixedBound]>L0-Mu:# velosity reflection when the particle reach the upper or bottom boundaries
                PV[i][FixedBound]=-abs(PV[i][FixedBound])
            if PL[i][FixedBound]<0+Mu:
                PV[i][FixedBound]=abs(PV[i][FixedBound])
            PL[i]=PL[i]*perB%L0-PL[i]*(perB-1)
            PF[i]=PF__

        for i in range(prcNr):
            pb=[0,0,0]
            __PL__=PL[i].copy()
            if PL[i][0]<cutoff:
                __PL__[0]=PL[i][0]+L0
                pb[0]=1*perB[0]
            elif PL[i][0]>L0-cutoff:
                __PL__[0]=PL[i][0]-L0
                pb[0]=1*perB[0]
            if PL[i][1]<cutoff:
                __PL__[1]=PL[i][1]+L0
                pb[1]=1*perB[1]
            elif PL[i][1]>L0-cutoff:
                __PL__[1]=PL[i][1]-L0
                pb[1]=1*perB[1]
            if PL[i][2]<cutoff:#
            #if PL[i][2]<0:
                __PL__[2]=PL[i][2]+L0#
                #PL[i][2]=-PL[i][2]
                pb[2]=1*perB[2]
            elif PL[i][2]>L0-cutoff:#
            #elif PL[i][2]>L0:
                __PL__[2]=PL[i][2]-L0#
                #PL[i][2]=2*L0-PL[i][2]
                pb[2]=1*perB[2]
            if pb[0]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([__PL__[0],PL[i][1],PL[i][2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[1]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([PL[i][0],__PL__[1],PL[i][2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[2]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([PL[i][0],PL[i][1],__PL__[2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[0]==1 and pb[1]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([__PL__[0],__PL__[1],PL[i][2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[0]==1 and pb[2]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([__PL__[0],PL[i][1],__PL__[2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[1]==1 and pb[2]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([PL[i][0],__PL__[1],__PL__[2]])))
                PS__=np.vstack((PS__,PS[i][0]))
            if pb[0]==1 and pb[1]==1 and pb[2]==1:
                flag.append(i)
                PL__=np.vstack((PL__,np.array([__PL__[0],__PL__[1],__PL__[2]])))
                PS__=np.vstack((PS__,PS[i][0]))
        if sliceNr%jumpSlice==0:
            xyzout(PL,PS,sliceNr)
        sliceNr+=1
        timeIni+=timeScale
    
if __name__=="__main__":
    start=time()
    main()
    print('the dynamic simulation takes %.1fs'%(time()-start))