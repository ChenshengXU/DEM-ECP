
import numpy as np
import os
import networkx as nx
from scipy import spatial
import math
from time import time

curloc=os.path.abspath(os.curdir)

k=0
tstep=0
L0=2000

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
CNTCondu=0.1

percoDirection=0

G =[]

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
    
    output='data/Regular.dump'
    os.path.join(curloc, output)
    file_output=open(output,'w')
    m=0
    for lines in file.readlines():
        if "ITEM: TIMESTEP" in lines:
            print(lines)
            file_output.write(lines)
            k=1
            continue
        if k==1:
            tstep=int(lines.strip('\n'))
            #file_output.write(lines)
            k=2
        if k==2:
            if "ITEM: TIMESTEP" not in lines:
                
                if "ITEM: BOX BOUNDS xy xz yz pp pp pp" in lines:
                    print(lines)
                    file_output.write(lines)
                    m=1
                    continue
                if m==1:
                    file_output.write('%i %i %i\n' %(-0/200,2000,-0/200))
                    m=2
                    continue
                if m==2:
                    if "ITEM: BOX BOUNDS xy xz yz pp pp pp" not in lines:
                        file_output.write(lines)
                        continue
                file_output.write(lines)
            else:
                file_output.close()
                k=0

    file.close()
    return tstep

def lastTimestepReader(timeStep):
    final_file='data/timestep_%i.xyz'%timeStep
    os.path.join(curloc, final_file)
    final_file=open(final_file,'r')
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

'''   
def DistMatrix(PI, CP):
    global G
    PI = sorted(PI, key=lambda x:x[0])
    
    partDist(PI)
    CI=[] #contact index
    for i in range(len(PI)):
        
        if len(PI[i]) == 1:#sphere
            PI[i].append(sphRa)
        else: #CNT
            PI[i].append(radius)
    
    prcNr = len(PI)
    PI=np.asarray(PI)
    #G.append(['r0000', 0,prcNr*2,(1-vol)*polyCondu**-1/L0])
    distMat=cal(PI[:,6:]/1000,PI[:,6:]/1000)
    #print(distMat[-1])
    indexID=np.argwhere(distMat!=10000)
    progress_bar = ProgressBar(len(indexID))
 
    for i in range(len(indexID)):
        progress_bar.update(i)
        if indexID[i][0]!=indexID[i][1]:
            re=2*resistance(distMat[indexID[i][0]][indexID[i][1]],PI[indexID[i][0]][-1],PI[indexID[i][1]][-1])
            if re!='InfInf' and PI[indexID[i][0]][5]!=PI[indexID[i][1]][5]:
                G.append(['r%s000%s%i'%(str(PI[indexID[i][0]][0]),str(PI[indexID[i][1]][0]),len(G)),PI[indexID[i][0]][0],PI[indexID[i][1]][0],re])
                CI.append([PI[indexID[i][0]][0],PI[indexID[i][0]][5]])
                CI.append([PI[indexID[i][1]][0],PI[indexID[i][1]][5]])
    
    CI= sorted(CI, key=lambda x:x[0])
    CI= sorted(CI, key=lambda x:x[1])
    CP= sorted(CP, key=lambda x:x[0])
    CP= sorted(CP, key=lambda x:x[1])
    Mid=CI[0][1]
    
    for i in range(1,len(CI)):
        if CI[i][0]!=Mid:
            Mid=CI[i][0]
            continue
        else:
            for j in CP:
                if j[1]==CI[i][1] and CI[i-1][1]<j[0]<CI[i][1]:
                    continue
                else:
                    re=(CI[i][1]-CI[i-1][1])*2/3.1416/radius/CNTCondu
                    G.append(['r%s000%s%i'%(str(CI[i-1][1]),str(CI[i][1]),CI[i][0]),CI[i-1][1],CI[i][1],re])
'''

def moleDivi(PI):
    prcNr=len(PI)
    
    PI = sorted(PI, key=lambda x:x[-4])
    PI_m=[]
    mid=0
    for atom in PI:
        if atom[-4]!=mid:
            mid=atom[-4]
            PI_m.append([atom])
        else:
            PI_m[-1].append(atom)
    for i in range(len(PI_m)):
        PI_m[i] = sorted(PI_m[i], key=lambda x:x[0])
    segNr=len(PI_m)
    return PI_m, prcNr, segNr

def distOf2Group(A1,A2):
    tree=spatial.cKDTree(A1[:,2:5])
    mat_mindist, mat_minid=tree.query(A2[:,2:5])
    mindist=min(mat_mindist)
    a1=mat_mindist.tolist()
    id1=a1.index(mindist)
    id2=mat_minid[id1]
    id1=int(A1[id1][0])
    id2=int(A2[id2][0])
    return mindist, (id1,id2)

def distMatrix(PI, CP, segNr):
    global G
    CN=CP.copy() #connected points (molecular id, particle id)
    mindist=0
    for i in range(len(PI)-1):
        print(i/len(PI))
        for j in range(i+1,len(PI)):
            A1=np.asarray(PI[i])
            A2=np.asarray(PI[j])
            mindist, minid=distOf2Group(A1,A2)
            if mindist>rtunnel:
                continue
            
            Re=resis(mindist,radius,radius)
            if Re!='Inf':
                m_id1=i+1
                if i+1 in CP[:][0]:
                    checkID1=CP[:][0].index(i+1)
                    if CP[checkID1][1]*CP[checkID1][2]<minid[0]*CP[checkID1][2]:  # need to be checked
                        m_id1=i+1+segNr
                m_id2=j+1        
                if j+1 in CP[:][0]:    
                    checkID2=CP[:][0].index(j+1)
                    if CP[checkID2][1]*CP[checkID2][2]<minid[1]*CP[checkID2][2]:# need to be checked
                        m_id2=j+1+segNr
                G.append(['r0%i%i' % (minid[0],minid[1]), minid[0], minid[1], Re])
                CN.append([m_id1,minid[0]])
                CN.append([m_id2,minid[1]])
    CN=sorted(CN, key=lambda x:x[1])
    return CN

def inPartRe(CN):
    global G
    for i in range(1,len(CN)):
        if CN[i][0]==CN[i-1][0] and CN[i][1]!=CN[i-1][1]:
            Re=abs(CN[i][1]-CN[i-1][1])*2/(math.pi*radius*CNTCondu)
            G.append(['r00%i%i' % (CN[i-1][1],CN[i][1]), CN[i-1][1], CN[i][1], Re])

def BC(PI,prcNr,segNr):
    global G 
    CP=[] #cut point
    for i in range(len(PI)):
        for j in range(len(PI[i])):

            PI[i][j][2]=PI[i][j][2]%L0
            PI[i][j][3]=PI[i][j][3]%L0
            PI[i][j][4]=PI[i][j][4]%L0
    
    for i in range(0,len(PI)):
        for j in range(0, len(PI[i])):
            
            if j>0 and abs(PI[i][j][percoDirection+2]-PI[i][j-1][percoDirection+2])>20:
                if PI[i][j][percoDirection+2]<PI[i][j-1][percoDirection+2]:
                    CP.append([PI[i][j][0],PI[i][j][5],1]) # the fiber goes through the boundary along the axis, the small number of particle are on the right side
                else:
                    CP.append([PI[i][j][0],PI[i][j][5],-1]) # the fiber goes through the boundary against the axis, the small number of particle are on the left side
                G.append(['r%i' % len(G), int(PI[i][j][0]+segNr), prcNr*2, '0.01k']) # right side will get connected with the fiber part with the id (m_id + total number of fiber)
                G.append(['r%i' % len(G), int(PI[i][j-1][0]),0,'0.01k'])
    return CP,PI       

def output(G,prcNr):
    file=open('data_set/data.cir','w')
    file.write('Simulation using vdd and resistors\n')
    file.write('vdd %s %s dc 2000v\n' %(prcNr*2, 0))
    for i in range(len(G)):
        file.write(" ".join(str(i) for i in G[i]))
        file.write('\n')
    file.write('.op i(1)\n')
    file.write('.end')
    file.close()

def main():
    tstep=dumpsplit('poly_r_2CB')
    
if __name__=="__main__":
    start=time()
    main()
    print("takes %.2fs" %(time()-start))