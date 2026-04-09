
import numpy as np
from random import sample, seed, choices, randint
import math
import sys
from time import time, sleep
import os
import matplotlib.pyplot as plt


L0=2000
aspRatio=150
vol=0.005
num=10
radius=5
segNr=100

RandomSeed = 5838317

curloc=os.path.abspath(os.curdir)

def volConverter():
    global num
    num=L0**3*vol/(3.14*radius**2*aspRatio*2)
    num=int(num)
    print('the volume ratio of the filling is %.2f%s with %i fibers.' %((3.14*radius**2*aspRatio*2*num)/L0**3*100,'%',num))
    print('the bond coeff should be set as %.2f.'%(aspRatio*radius*2/segNr))

def dataout(PLAtoms,PLBonds,PLAngles):
    packet_file = 'data_set/test.data'
    os.path.join(curloc, packet_file)
    file_data=open(packet_file,'w')
    file_data.write('first line is always a comment\n')
    file_data.write('\n')
    file_data.write('%i atoms\n'% len(PLAtoms))
    file_data.write('%i bonds\n'% len(PLBonds))
    file_data.write('%i angles\n'% len(PLAngles))
    file_data.write('\n')
    file_data.write('1 atom types\n')
    file_data.write('1 bond types\n')
    file_data.write('1 angle types\n')
    file_data.write('\n')
    file_data.write('0.0 %.2f xlo xhi\n'%(L0))
    file_data.write('0.0 %.2f ylo yhi\n'%(L0))
    file_data.write('0.0 %.2f zlo zhi\n'%(L0))
    file_data.write('\n')
    file_data.write('Masses\n')
    file_data.write('\n')
    file_data.write('1 16\n')
    file_data.write('\n')
    file_data.write('Atoms\n')
    file_data.write('\n')

    for i in range(len(PLAtoms)):
        file_data.write('%i %i 1 '%(i+1, (i//(segNr+1)+1)))        
        file_data.write(' %.4f  %.4f  %.4f\n'%(PLAtoms[i][0][0],PLAtoms[i][0][1],PLAtoms[i][0][2]))
    
    file_data.write('\n')
    file_data.write('Bonds\n')    
    file_data.write('\n')
    
    for i in range(len(PLBonds)):
        file_data.write('%i 1   '%(i+1))
        file_data.write(' %i  %i  \n'%(PLBonds[i][0]+1,PLBonds[i][1]+1))
        
    file_data.write('\n')
    file_data.write('Angles\n')
    file_data.write('\n')
    
    for i in range(len(PLAngles)):
        file_data.write('%i %i '%(i+1,1))
        file_data.write(' %i  %i  %i  \n'%(PLAngles[i][0]+1,PLAngles[i][1]+1,PLAngles[i][2]+1))
    
    file_data.close()

def main():
    global num
    
    volConverter()
    
    rng = np.random.RandomState(RandomSeed)
    
    PLIni=rng.uniform(0,L0,(num,3))
    
    PLAtoms=[]
    PLBonds=[]
    PLAngles=[]
    orien=rng.uniform(-1,1,(1,3))
    orien=orien/np.linalg.norm(orien)
    for k in range(num):
        
        
        PLAtoms.append([PLIni[k]])
        num_before=k*(segNr+1)
        sys.stdout.write("\r{0}".format((k+1)/num))
        sys.stdout.flush()
        for i in range(1,segNr+1):
            PLAtoms.append(PLIni[k]+(i)*orien*aspRatio*radius*2/segNr)
            PLBonds.append([num_before+i-1,num_before+i,k])
            if i!=segNr:
                PLAngles.append([num_before+i-1,num_before+i,num_before+i+1,k])
                #print(i)
                                   
    print('\n')
    dataout(PLAtoms,PLBonds,PLAngles)
if __name__=="__main__":
    
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    