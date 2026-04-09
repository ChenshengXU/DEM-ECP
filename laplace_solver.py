
#import reader_pp as rp
import numpy as np
import networkx as nx
import os
from time import time
import multiprocessing as mp
import sys


class ProgressBar():

    def __init__(self, max_steps):
        self.max_steps = max_steps
        self.current_step = 0
        self.progress_width = 50

    def update(self, step=None):
        self.current_step = step

        num_pass = int(self.current_step * self.progress_width / self.max_steps) + 1
        num_rest = self.progress_width - num_pass
        percent = (self.current_step + 1) * 100.0 / self.max_steps
        progress_bar = '[' + '\u25a0' * (num_pass - 1) + '|' + '-' * num_rest + ']'
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


#num_cores = int(mp.cpu_count() * rp.rate)
#num_cores = int(rp.num_cores)
num_cores=4
used_cores = num_cores
#print(num_cores,'cores will be used')
epsilon = 0.001
x_k1 = np.array([])
b_k1 = np.array([])
Tran=np.array([])

curloc = os.path.abspath(os.curdir)


def readCIR(fileName):
    G = []
    batch_file = 'data_set/%s.cir' % fileName
    os.path.join(curloc, batch_file)
    file = open(batch_file, 'r')
    k = 0
    print('reading start')
    for lines in file.readlines():
        
        if "dc" in lines:
            k = 1
            slist = lines.strip('\n').split()
            prcNr = int(slist[1])
            continue
        if 'op' in lines:
            k = 0
            continue
        if k == 1:
            nlist = lines.strip('\n').split()
            G.append([int(nlist[1]), int(nlist[2]), float(nlist[3])])
    print('reading finished')
    return G, prcNr


def G2GP(G):
    print('converting start')
    prcNr = max(max(G))
    GP = nx.Graph()
    GP.add_node(0)
    GP.add_node(prcNr)
    for i in range(len(G)):
        
        nx.add_path(GP, [G[i][0],G[i][1]])
    print('converting finished')
    return GP, prcNr


def nodesMapping(G, nodeNr,node_degree):
    mapping = np.zeros((nodeNr-1, 2), dtype=int)
    pointer = 0
    #print(node_degree)
    print('mapping start')
    '''
    for line in G:
        print(pointer/nodeNr)
        if line[0] not in mapping[:, 0]:
            mapping[pointer][0] = line[0]
            mapping[pointer][1] += 1
            pointer += 1
        else:
            idNr = np.where(mapping[:, 0] == line[0])
            idNr = idNr[0][0]
            mapping[idNr][1] += 1
        if line[1] not in mapping[:, 0]:
            mapping[pointer][0] = line[1]
            mapping[pointer][1] += 1
            pointer += 1
        else:
            idNr = np.where(mapping[:, 0] == line[1])
            idNr = idNr[0][0]
            mapping[idNr][1] += 1
    '''

    for info in node_degree:
        
        mapping[pointer-1][0]=info[0]
        mapping[pointer-1][1]=info[1]
        pointer+=1
        print(1,pointer/nodeNr, end = "\r")
    # diagonal matrix [Node ID, degree]
    D = mapping.copy()

    G = np.asarray(G)
    '''
    delID=np.where(D[:,0]==0)[0]
    print(delID)
    blindBranch=np.where(D[:,0]==1)
    print(blindBranch)
    for idNr in blindBranch:
    	print(idNr/nodeNr)
    	for blindBranchidNr in np.where(G == D[idNr][0]):
    		print(np.where(G == D[idNr][0]))
    		pairNr = 1 - blindBranchidNr[1]
    		mapNr = np.where(mapping == G[idNr[pairNr]])[0][0]
    		D[mapNr][1]-=1
    '''
    '''
    for i in range(len(D)):
        print(i/nodeNr)
        if D[i][1] == 1:
            for blindBranchidNr in np.where(G == D[i][0])[0]:
                pairNr = 1 - np.where(G == D[i][0])[1][0]
                mapNr = np.where(mapping == G[blindBranchidNr][pairNr])[0][0]
                mapping[mapNr][1] -= 1
            delID.append(i)
        if D[i][1] == 0:
            delID.append(i)
    
    for ID in range(len(delID) - 1, -1, -1):
        D = np.delete(D, (delID[ID]), axis=0)
    nodeNr = len(D)
    '''
    A = {}
    for pair in range(nodeNr-1):
        print(2,pair/nodeNr, end = "\r")
        A[int(D[pair][0])] = []
        '''
        A[int(D[pair][0])] = np.zeros(nodeNr)
        A[int(D[pair][0])][pair] = int(D[pair][1])
        '''
    print('mapping finished')
    return mapping, D, nodeNr, A


def condM(G, D):
    global Tran
    print('connectivity building')
    #used_cores=1
    num_task=len(G)
    Tran=np.zeros((num_task,5))
    nested_dict = {}
    param_dict={}
    pool1 = mp.Pool(used_cores)
    print('mp scanning')

    for i in range(used_cores):
        
        nested_dict[i]=[]
        try:
            nested_dict[i]=G[(num_task//used_cores+1)*i:(num_task//used_cores+1)*(i+1)][:]
            param_dict.update({i: list(range((num_task//used_cores+1)*i,(num_task//used_cores+1)*(i+1)))})
        except IndexError:
            nested_dict[i]=G[(num_task//used_cores+1)*i:num_task][:]
            param_dict.update({i: list(range((num_task//used_cores+1)*i,num_task))})
    manager = mp.Manager()
    
    for name, param in param_dict.items():
        pool1.apply_async(condPP, args=(name,used_cores,nested_dict[name],D, num_task,param,), callback=result_app_cond)
        
    pool1.close()
    pool1.join()

    #result_app_cond(condPP(0,used_cores,G,D, num_task,range(0,num_task)))
    print('scanning finished')

def condPP(name,used_cores,G,D,num_task, param):
    resultxx=np.zeros((num_task//used_cores+1,6))
    for i in param:
        
        id=i-name*(num_task//used_cores+1)
        #id=i
        print('2.5 PP',4*id/num_task, end = "\r")
        id1 = np.where(D[:, 0] == G[id][0])
        id1 = id1[0][0]
        id2 = np.where(D[:, 0] == G[id][1])
        id2 = id2[0][0]
        resi = G[id][-1]
        resultxx[id]=np.array([i,id1,id2,G[id][0],G[id][1],1 / resi])
    return resultxx

def result_app_cond(result):
    global Tran
    print('try to output')
    for line in result:
        if line[1]!=0:
            Tran[int(line[0])]=line[1:]

def laterAc(A):    
    n=0
    lend=len(Tran)
    
    for i in range(lend):
        if Tran[i][2]*Tran[i][3]!=0:         
            A[int(Tran[i][2])].append([int(Tran[i][3]),Tran[i][4]])
            A[int(Tran[i][3])].append([int(Tran[i][2]),Tran[i][4]])
    for i in A.copy():
        if not A[i]:
            A.pop(i)
            #print(3,(i+1)/lend)
            #print(A[int(Tran[i][2])])
    '''
    for line in G:
        n+=1
        lend=len(G)
        print(3,n/lend)
        if line[0] in D[:, 0] and line[1] in D[:, 0]:
            id1 = np.where(D[:, 0] == line[0])
            id1 = id1[0][0]
            id2 = np.where(D[:, 0] == line[1])
            id2 = id2[0][0]
            resi = line[-1]
    
            #A[line[0]][id2] = 1 / resi
            #A[line[1]][id1] = 1 / resi
    
            A[line[0]].append([id2,1 / resi])
            A[line[1]].append([id1,1 / resi])
    '''
    
    print('connectivity finished')
    return A



def jacNewIt(b, D, A, nodeNr, id_a, id_c):
    global x_k1, b_k1
    #x_k1 = np.random.rand(nodeNr, 1)  # Nr.0 Itration
    x_k1 =np.zeros((nodeNr,1))
    x_k =np.zeros((nodeNr,1))
    #x_k = np.random.rand(nodeNr, 1)
    x_k1[id_a] -= 4
    x_k1[id_c] += 4
    err = 1
    b_tran = 0
    ItNr=0
    while err > epsilon :
        b_k1 = np.zeros(nodeNr)
        x_k = x_k1.copy()
        
        jacIter(x_k, D, b, A, nodeNr)  # Resterror estimation

        err=0
        for i in range(len(b)):
            err+=abs(b[i]-b_k1[i])
        #err = (np.linalg.norm(b) - np.linalg.norm(b_k1))
        EquRe = (x_k1[id_a] - x_k1[id_c]) 
        save_file = 'data_set/his_err_pp.csv'
        os.path.join(curloc, save_file)
        file = open(save_file, 'a+')
        file.writelines('%.5f   %.5f Omu\n' % (err, EquRe))
        print('err:', err)
        ItNr+=1
        if ItNr>10000:
            EquRe=0
            break
    return EquRe


def jacIter(x_k, D, b, A, nodeNr):
    global x_k1, b_k1
    param_dict = {}

    pool = mp.Pool(used_cores)
    num_task = nodeNr-1
    nested_dict = {}
    
    for i in range(used_cores):
        param_dict.update({i: list(range(i, num_task, used_cores))})
        nested_dict[i] = {}
        for id_dict in range(i, num_task, used_cores):
            nested_dict[i].update({D[id_dict][0]: A.get(D[id_dict][0])})
    manager = mp.Manager()
    for name, param in param_dict.items():
        pool.apply_async(ppCal, args=(x_k, D, b, nested_dict[name], nodeNr, param,), callback=result_app)
    pool.close()
    pool.join()
    
    result_app(ppCal(x_k, D, b, A, nodeNr, range(0,len(A))))

def ppCal(x_k, D, b, A, nodeNr, param):
    x_k1__ = np.zeros(nodeNr)
    b__ = np.zeros(nodeNr)
    outcome=[]
    for name,param in A.items():
        
        summe = 0
        D_tran = 0
        id_m = np.where(D[:, 0] == name)
        id_m = id_m[0][0]
        for sp in param:
            id_n = np.where(D[:, 0] == sp[0])
            id_n = id_n[0][0]
            
            summe-=x_k[id_n][0]/sp[1]
            D_tran+=1/sp[1]
        
        #b__[num] = (summe + A[D[num][0]][num] * x_k[num][0])
        b__[id_m]=summe+D_tran*x_k[id_m][0]
        
        b_tran = b[id_m] - summe
        x_k1__[id_m] = ((1 / D_tran) * b_tran-x_k[id_m])*2+x_k[id_m]
        outcome.append(id_m)
        # print('result:', [num, x_k1__[-1], b__[-1]])
    return [outcome, x_k1__, b__]


def result_app(result):
    global x_k1, b_k1
    for num in result[0]:
        x_k1[num] = result[1][num]
        b_k1[num] = result[2][num]


def main():  # mainly to solve the laplance matrix with jacobi-iteration and newton-raphon-iteration
    #G, prcNr = readCIR('data')
    #steplist=rp.main()
    steplist=[2000000]
    for steppoint in steplist:
        G, prcNr = readCIR('data_step_%i'%steppoint)
        GP, prcNr = G2GP(G)
        
        EquRe=1e12
        if nx.has_path(GP, 1, prcNr):  # only when it's connected, it can progress further
            nodeNr = GP.number_of_nodes()
            node_list=GP.nodes()
            node_degree=GP.degree(node_list)
            NP, D, nodeNr, A = nodesMapping(G,
                                            nodeNr, node_degree)  # obtain the node ID mapping and corresponding degree of the nodes, which are more than 1
            
            condM(G, NP)  # A is the 2D conductance matrix(present as dictionary, with the row number is the key); D is the diagonal matrix of conductance matrix, notice the element should be inversed; R is the rest part
            A=laterAc(A)
            # wehere the D[_] is -1, the corresponding cow and culomn will be removed and the other D should be recalculated; here the node number will be that one exclude the blind connection
            b = np.zeros(nodeNr)  # current vector
            id_a = np.where(D[:, 0] == 1)
            id_a = id_a[0][0]
            id_c = np.where(D[:, 0] == prcNr)
            id_c = id_c[0][0]
            b[id_a] = -1
            b[id_c] = 1
            EquRe = jacNewIt(b, D, A, nodeNr, id_a, id_c)
            
        save_file = 'data_set/step_cond_coo.csv'
        save_file=os.path.join(curloc, save_file)
        with open(save_file,'a') as file:
            file.writelines('step Nr. %i   %.5f Omu\n' % (steppoint, EquRe))
            
            #print(EquRe)


if __name__ == "__main__":
    start = time()
    main()
    print('the simulation takes %.1fs' % (time() - start))

