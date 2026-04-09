#!/usr/bin/env python

import numpy as np
import os
import networkx as nx
from scipy import spatial
import math
from time import time
import multiprocessing as mp
from pathlib import Path
import sys
import socket

rate=1
num_cores = int(mp.cpu_count()*rate) # number of cores in the computer
num_cores=min(10,num_cores)
used_cores=200
num_cores=used_cores
#used_cores=200 #num_cores # default use all the cores in the computer



curloc=os.path.abspath(os.curdir)
#curloc='/home/st/st_st/st_ac136602/ecp'
k=0
tstep=0
L0=5000

sphRa=20
vols=0.0126
radius=5
volf=0.01

tou_1 = 22.9
tou_2 = 17.7
ConduQuan = 7.748091729e-5
M = 400
rtunnel=14

polyCondu=10**-17
fiCondu=4000e-9
CNTCondu=0.001

percoDirection=0

G =[]
CN=[]

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

def init_list_of_objects(size):
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() ) #different object reference each time
    return list_of_objects

def dumpsplit(fileName):

    batch_file = '%s.dump'%fileName
    os.path.join(curloc, batch_file)
    file=open(batch_file,'r')
    steplist=[]
    for lines in file.readlines():
        #print(lines)
        if "ITEM: TIMESTEP" in lines:
            print(lines)
            k=1
            continue
        if k==1:
            tstep=int(lines.strip('\n'))
            output='data/timestep_%i.xyz'%tstep
            os.path.join(curloc, output)
            output=Path(output)
            output.touch(exist_ok=True)
            file_output=open(output,'w+')
            file_output.write(lines)
            steplist.append(tstep)
            k=2
        if k==2:
            if "ITEM: TIMESTEP" not in lines:
                file_output.write(lines)
            else:
                file_output.close()
                k=0
    file.close()
    return tstep,steplist

def lastTimestepReader(timeStep):
    final_file_name='data/timestep_%i.xyz'%timeStep
    final_file_name=os.path.join(curloc, final_file_name)
    #print(curloc,final_file_name)
    final_file=open(final_file_name,'r')
    k=0
    #PI=init_list_of_objects(AgggNr)
    PI=[]
    
    for lines in final_file.readlines():
        if "ITEM: ATOMS id" in lines:
            k=1
            continue
        if k==1:
            nlist=lines.strip('\n').split()
            nlist=list(map(float,nlist))
            PI.append(nlist)
        #print(lines)    
    return PI


def resis(distance,r1,r2):
    #R_in=fiCondu*Mu*2/(Mu**2*0.01*3.14)
    distance=distance-r1-r2
    r=min(r1,r2)
    pix=1
    #R_in=0.195/fiCondu/r
    R_in=0.01/CNTCondu/r
    
    if -r1-r2<distance <= 0:
        R = R_in*(r1+r2+distance)/(r1+r2)-3.14*(0.66*math.log(100*r+1))**-2*(-distance)*CNTCondu # Soft-core nur angewendet für den Startpunkt in Boundary    
    elif distance <= 0.34:
        conduc=M * ConduQuan * math.exp(-tou_1 * 0.34)
        R = 1/conduc+R_in
    elif distance < 0.6:
        conduc=M * ConduQuan * math.exp(-tou_1 * distance)
        R = 1/conduc+R_in  
    elif distance <= 1.4:
        conduc = M * ConduQuan * math.exp(-tou_2 * distance)
        R = 1/conduc+R_in 
    elif 0<=distance <= pix:
        R_inter = 0.1*(1.5*distance)**2/(CNTCondu*pix/distance*0.08*4/3*3.14*(r1**3+r2**3))
        R = R_inter+R_in  
    elif pix<=distance <= (r1+r2)/20:
        R=1/(-np.log(distance/r)*polyCondu*3.14*r)+R_in
    else:
        R = 'Inf'
    return R

def moleDivi(PI):
    prcNr=len(PI)
    PI = np.asarray(PI)
    PI=np.asarray(PI[:,:6])
    PI = sorted(PI, key=lambda x:x[0])
    PI_m=np.asarray(PI[:])
    PI=PI_m[:,2:6]

    return PI, prcNr


def get_key (dict):
    return [k for k, v in dict.items() if v != 0]

def distMatrix(PI):
    global CN,G
    nested_dict = {}
    param_dict={}
    kd_tree = spatial.cKDTree(PI)
    sdm = kd_tree.sparse_distance_matrix(kd_tree, rtunnel,output_type='dict')
    k=get_key(sdm) #(nodeID1+1,nodeID2)
    pool=mp.Pool(used_cores)
    num_task=len(k) ###-1
    # ppCal(PI,sdm,k)
    for i in range(used_cores):
        nested_dict[i]=[]
        try:
            nested_dict[i]=k[(num_task//used_cores+1)*i:(num_task//used_cores+1)*(i+1)][:]
            param_dict.update({i: list(range((num_task//used_cores+1)*i,(num_task//used_cores+1)*(i+1)))})
        except IndexError:
            nested_dict[i]=k[(num_task//used_cores+1)*i:num_task][:]
            param_dict.update({i: list(range((num_task//used_cores+1)*i,num_task))})

    # manager = mp.Manager()
    # managed_locker = manager.Lock()
    # managed_dict = manager.dict()
    for name, param in param_dict.items():
        pool.apply_async(ppCal, args=(PI,sdm,nested_dict[name]),callback=result_app)
    pool.close()
    pool.join()
    

def ppCal(PI,sdm,k):
    global G,CN
    R1=[]
    R2=[]
    R3=[]
    for num,paar in enumerate(k): # nodeID1 is paar[0]+1; node ID2 is paar[1] 
        if PI[paar[0]+1][-1]!=PI[paar[1]][-1]:
            dist=sdm[paar]
            Re=resis(dist,radius,radius)
            if Re!='Inf': #check if the connected one is from the other side of the fibre
                print('/')
                m_id1=int(PI[paar[0]+1][-1])
                m_id2=int(PI[paar[1]][-1])
                if Re<0:
                    Re=0.01
                G.append(['r0%i%i' % (paar[0]+2,paar[1]+1), paar[0]+2,paar[1]+1, Re])
                CN.append([m_id1,paar[0]+1])
                CN.append([m_id2,paar[1]])
    return [R1,R2,R3]

def result_app(result):
    global G,CN
    G+=result[0]
    CN+=result[1]
    CN+=result[2]

def inPartRe():
    global G,CN
    CN=sorted(CN, key=lambda x:x[1])
    for i in range(1,len(CN)):
        if CN[i][0]==CN[i-1][0] and abs(CN[i][1]-CN[i-1][1])>1:
            Re=abs(CN[i][1]-CN[i-1][1])*2/(math.pi*radius*CNTCondu)
            G.append(['r00%i%i' % (CN[i-1][1]+1,CN[i][1]+1), CN[i-1][1]+1, CN[i][1]+1, Re])

def BC(PI,prcNr):
    global G,CN
    for i in range(1,len(PI)-1): 
        if PI[i][percoDirection]>L0 and PI[i+1][percoDirection]<L0 and PI[i][-1]==PI[i+1][-1]:
            G.append(['r%i' % len(G), i+3, prcNr*2, '1']) 
            G.append(['r%i' % len(G), i+2,1,'1'])
            CN.append([PI[i][-1],i])
            CN.append([PI[i][-1],i+1])
        elif PI[i][percoDirection]<L0 and PI[i+1][percoDirection]>L0 and PI[i][-1]==PI[i+1][-1]:
            G.append(['r%i' % len(G), i+2, prcNr*2, '1']) 
            G.append(['r%i' % len(G), i+3,1,'1'])
            CN.append([PI[i][-1],i])
            CN.append([PI[i][-1],i+1])
        if PI[i][percoDirection]>0 and PI[i+1][percoDirection]<0 and PI[i][-1]==PI[i+1][-1]:
            G.append(['r%i' % len(G), i+2, prcNr*2, '1']) 
            G.append(['r%i' % len(G), i+3,1,'1'])
            CN.append([PI[i][-1],i])
            CN.append([PI[i][-1],i+1])
        elif PI[i][percoDirection]<0 and PI[i+1][percoDirection]>0 and PI[i][-1]==PI[i+1][-1]:
            G.append(['r%i' % len(G), i+3, prcNr*2, '1']) 
            G.append(['r%i' % len(G), i+2,1,'1'])
            CN.append([PI[i][-1],i])
            CN.append([PI[i][-1],i+1])
        
    PI[:,percoDirection]=PI[:,percoDirection]%L0
    return PI       

def output(G,prcNr,steppoint):
    output=Path('data_set/data_step_%i.cir'%steppoint)
    output.touch(exist_ok=True)
    file=open('data_set/data_step_%i.cir'%steppoint,'w+')
    file.write('Simulation using vdd and resistors\n')
    file.write('vdd %s %s dc 2000v\n' %(prcNr*2, 1))
    for i in range(len(G)):
        file.write(" ".join(str(i) for i in G[i]))
        file.write('\n')
    file.write('.op i(1)\n')
    file.write('.end')
    file.close()

def main():
    tstep,steplist=dumpsplit('Regular')
    #steplist=range(0,2550000,50000)
    prcNrlist=[]
    for steppoint in steplist:
        print('prepare circuit file for step %i'%steppoint)
        PI=lastTimestepReader(steppoint)
        PI, prcNr=moleDivi(PI)
        PI=BC(PI,prcNr)
        distMatrix(PI)
        inPartRe()
        output(G,prcNr,steppoint)
        prcNrlist.append(prcNr)
    return steplist
    
if __name__=="__main__":
    start=time()
    #print(socket.gethostname())
    main()
    print("takes %.2fs" %(time()-start))
