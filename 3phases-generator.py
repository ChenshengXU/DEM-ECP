
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

safeVol_CNT=0.0001

#CB info
sphVol=0.001
segNr=150
sphNr=0
sphRa=20

#gnp info
gnp_num=50
width=400
length=480
num_w=50
num_l=60
gnp_nodes=num_w*num_l*gnp_num
orien_m=np.array([1,3,2])
orien_s=np.array([2,4,5])

gnp_vol=(width/num_w/2)*width*length*gnp_num/L0**3

safeVol_GNP=safeVol_CNT*segNr/num_w/num_l

RandomSeed = 5838317

curloc=os.path.abspath(os.curdir)

def orienConverter(ori1,ori2):
    global orien_m,orien_s
    ori2=np.cross(ori1,ori2)
    ori1=ori1/np.linalg.norm(ori1)*width/num_w
    ori2=ori2/np.linalg.norm(ori2)*length/num_l
    return ori1,ori2

def volConverter():
    global num,sphNr,L0
    num=L0**3*vol/(3.14*radius**2*aspRatio*2)
    num=int(num)
    sphNr=int(sphVol*L0**3/(4/3*3.14*sphRa**3))
    
    
    
    print('the volume ratio of the filling is %.2f%s with %i fibers.' %((3.14*radius**2*aspRatio*2*num)/L0**3*100,'%',num))
    print('the bond coeff should be set as %.2f'%(aspRatio*radius*2/segNr))
    print('the volume ratio of the sphere is %.2f%s with %i particles' %((3.14*sphRa**3*4/3*100*sphNr/L0**3),'%',sphNr))
    print('the pair coeff between 2 and 3 should be %.2f'%(radius+sphRa))
    print('the GNP vol ratio is %.2f%s, the pair coeff between 1 and 2 should be %.2f, the pair coeff between 1 and 3 should be %.2f' %(gnp_vol*100,'%', radius+width/num_w/2, width/num_w/2+sphRa))
    print('the radius of the matrix particle is %.2f' %(L0/200))
    
    L0=L0*(max(vol/safeVol_CNT,gnp_vol/safeVol_GNP)**(1/3))
    print('the box should be compressed from %.2f to %.2f'%(L0,L0/(max(vol/safeVol_CNT,gnp_vol/safeVol_GNP)**(1/3))))
    

def dataout(PLAtoms,PLBonds,PLAngles,CNTNode,PLDihedrals):
    packet_file = 'test.data'
    os.path.join(curloc, packet_file)
    file_data=open(packet_file,'w')
    file_data.write('first line is always a comment\n')
    file_data.write('\n')
    file_data.write('%i atoms\n'% (len(PLAtoms)))
    print('there are totally %i of particles need to be generated' %len(PLAtoms))
    file_data.write('%i bonds\n'% len(PLBonds))
    file_data.write('%i angles\n'% len(PLAngles))
    file_data.write('%i dihedrals\n'% len(PLDihedrals))
    file_data.write('\n')
    if sphVol*vol!=0:
    	file_data.write('3 atom types\n')
    else:
    	file_data.write('2 atom types\n')
    file_data.write('\n')
    if vol!=0:
        file_data.write('2 bond types\n')
        file_data.write('2 angle types\n')
        file_data.write('1 dihedral types\n')
    file_data.write('\n')
    file_data.write('0.0 %.2f xlo xhi\n'%(L0))
    file_data.write('0.0 %.2f ylo yhi\n'%(L0))
    file_data.write('0.0 %.2f zlo zhi\n'%(L0))
    file_data.write('0 0 0 xy xz yz\n')
    file_data.write('\n')
    file_data.write('Masses\n')
    file_data.write('\n')
    
    if vol!=0:
        file_data.write('2 16\n')
        if sphVol!=0:
            file_data.write('3 2304.0\n')
            #file_data.write('4 150\n')
            file_data.write('1 1.024\n')
    else:
        file_data.write('2 2304.0\n')
    file_data.write('\n')
    file_data.write('Atoms\n')
    file_data.write('\n')

    for i in range(gnp_nodes):
        file_data.write('%i %i 1 '%(i+1, (i//(num_l*num_w)+1)))
        file_data.write(' %.4f  %.4f  %.4f\n'%(PLAtoms[i][0],PLAtoms[i][1],PLAtoms[i][2]))

    for i in range(gnp_nodes,gnp_nodes+CNTNode):
        file_data.write('%i %i 2 '%(i+1, gnp_num+1+(i-gnp_nodes)//(segNr+1)))        
        file_data.write(' %.4f  %.4f  %.4f\n'%(PLAtoms[i][0][0],PLAtoms[i][0][1],PLAtoms[i][0][2]))

    if vol!=0:
    	mark=3
    else:
    	mark=2
    for i in range(gnp_nodes+CNTNode, len(PLAtoms)):
    	file_data.write('%i %i %i '%(i+1, gnp_num+num+1+i-gnp_nodes-CNTNode, mark))        
    	file_data.write(' %.4f  %.4f  %.4f\n'%(PLAtoms[i][0][0],PLAtoms[i][0][1],PLAtoms[i][0][2]))
    
    length=len(PLAtoms)
    aktNr=CNTNode//(segNr+1)+length-CNTNode
    '''
    for i in range(100):
        for j in range(100):
            for k in range(100):
                aktNr+=1
                file_data.write('%i %i 4 '%(length+k*10000+j*100+i+1,aktNr++k*10000+j*100+i))        
                file_data.write(' %.4f  %.4f  %.4f\n'%(L0/200*(2*i+1),L0/200*(2*j+1),L0/200*(2*k+1)))
    '''
    if vol!=0:
        file_data.write('\n')
        file_data.write('Bonds\n')    
        file_data.write('\n')
        
        for i in range(gnp_bonds):
            file_data.write('%i 1   '%(i+1))
            file_data.write(' %i  %i  \n'%(PLBonds[i][0]+1,PLBonds[i][1]+1))

        for i in range(gnp_bonds,len(PLBonds)):
            file_data.write('%i 2   '%(i+1))
            file_data.write(' %i  %i  \n'%(PLBonds[i][0]+1,PLBonds[i][1]+1))
        
        file_data.write('\n')
        file_data.write('Angles\n')
        file_data.write('\n')
        
        for i in range(gnp_angles):
            file_data.write('%i %i '%(i+1,1))
            file_data.write(' %i  %i  %i  \n'%(PLAngles[i][0]+1,PLAngles[i][1]+1,PLAngles[i][2]+1))
        
        for i in range(gnp_angles,len(PLAngles)):
            file_data.write('%i %i '%(i+1,2))
            file_data.write(' %i  %i  %i  \n'%(PLAngles[i][0]+1,PLAngles[i][1]+1,PLAngles[i][2]+1))
        
        file_data.write('\n')
        file_data.write('Dihedrals\n')    
        file_data.write('\n')    
    
        for i in range(len(PLDihedrals)):
            file_data.write('%i  1'%(i+1))
            file_data.write(' %i  %i  %i  %i\n'%(PLDihedrals[i][0]+1,PLDihedrals[i][1]+1,PLDihedrals[i][2]+1,PLDihedrals[i][3]+1))
    
        file_data.close()

def main():
    global num, orien_m,orien_s,gnp_bonds,gnp_angles
    
    orien_m,orien_s=orienConverter(orien_m,orien_s)
    rng = np.random.RandomState(RandomSeed)
    
    volConverter()
    
    rng = np.random.RandomState(RandomSeed)
    
    PLIni=rng.uniform(0,L0,(gnp_num+num+sphNr,3))
    
    PLAtoms=[]
    PLBonds=[]
    PLAngles=[]
    PLDihedrals=[]
    
    for k in range(gnp_num):
        orien_m=rng.uniform(0,L0,(1,3))
        orien_s=rng.uniform(0,L0,(1,3))
        orien_m,orien_s=orienConverter(orien_m,orien_s)
        #print(orien_m,orien_s)
        for j in range(num_l):
            for i in range(num_w):
                PLAtoms.append(PLIni[k]+i*orien_m[0]+j*orien_s[0])
                num_total=num_w*num_l*k
                
                if i!=num_w-1:
                    PLBonds.append([i+j*num_w+num_total,i+j*num_w+1+num_total,k])
                    #print([i+j*num_w+num_total,i+j*num_w+1+num_total], i, j)
                    if j!=num_l-1:
                        PLAngles.append([i+(j+1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w+1+num_total,k])
                        PLDihedrals.append([i+(j+1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w+1+num_total,i+(j+1)*num_w+num_total+1,k])
                    if j!=0:
                        PLAngles.append([i+(j-1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w+1+num_total,k])
                        PLDihedrals.append([i+(j-1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w+1+num_total,i+(j-1)*num_w+num_total+1,k])
                if j!=num_l-1:
                    PLBonds.append([i+j*num_w+num_total,i+(j+1)*num_w+num_total,k])
                    if i!=0:
                        PLAngles.append([i+(j+1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w-1+num_total,k])
                        PLDihedrals.append([i+(j+1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w-1+num_total,i+(j+1)*num_w+num_total-1,k])
                if i!=0 and j!=0:
                    PLAngles.append([i+(j-1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w-1+num_total,k])   
                    PLDihedrals.append([i+(j-1)*num_w+num_total,i+j*num_w+num_total,i+j*num_w-1+num_total,i+(j-1)*num_w+num_total-1,k]) 
    
    gnp_bonds=len(PLBonds)
    gnp_angles=len(PLAngles)
    
    for k in range(gnp_num,gnp_num+num):
        orien=rng.uniform(-1,1,(1,3))
        orien=orien/np.linalg.norm(orien)
        
        PLAtoms.append([PLIni[k]])
        num_before=(k-gnp_num)*(segNr+1)+gnp_nodes
        
        for i in range(1,segNr+1):
            PLAtoms.append(PLIni[k]+(i)*orien*aspRatio*radius*2/segNr)
            PLBonds.append([num_before+i-1,num_before+i,k])
            if i!=segNr:
                PLAngles.append([num_before+i-1,num_before+i,num_before+i+1,k])
                #print(i)
    CNTNode=len(PLAtoms)-gnp_nodes
    
    
    for j in range(gnp_num+num,gnp_num+num+sphNr):
    	PLAtoms.append([PLIni[j]])

    dataout(PLAtoms,PLBonds,PLAngles,CNTNode,PLDihedrals)
if __name__=="__main__":
    
    main()
    
    