import sys
import numpy as np

def dist_3D_Segment_to_Segment(Segment_S1, Segment_S2):

    S1_P0 = Segment_S1[0, :]
    S1_P1 = Segment_S1[1, :]
    S2_P0 = Segment_S2[0, :]
    S2_P1 = Segment_S2[1, :]
    u = S1_P1 - S1_P0
    v = S2_P1 - S2_P0
    w = S1_P0 - S2_P0
    a = np.dot(u, u)
    b = np.dot(u, v)
    c = np.dot(v, v)
    d = np.dot(u, w)
    e = np.dot(v, w)
    D = a * c - b * b
    sc = D
    sN = D
    sD = D
    tc = D
    tN = D
    tD = D
    
    # compute the line parameters of the two closest points
    if D < sys.float_info.epsilon: # epsilon is a very small number, here the lines are almost parallel
        sN = 0.0         # force using point P0 on segment S1
        sD = 1.0         # to prevent possible division by 0.0 later
        tN = e
        tD = c
    
    else:  # get the closest points on the infinite lines
        sN = (b * e - c * d)
        tN = (a * e - b * d)
        if sN < 0.0:       # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif sN > sD:  # sc > 1  => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c
    
    if tN < 0.0:            # tc < 0 => the t=0 edge is visible
        tN = 0.0
        
        # recompute sc for this edge
        if -d < 0.0:
            sN = 0.0
        elif -d > a:
            sN = sD
        else:
            sN = -d
            sD = a
    
    elif tN > tD:     # tc > 1  => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if (-d + b) < 0.0:
            sN = 0
        elif (-d + b) > a:
            sN = sD
        else:
            sN = (-d +  b)
            sD = a
    
    # finally do the division to get sc and tc
    if abs(sN) < sys.float_info.epsilon:
        sc = 0.0
    else:
        sc = sN / sD
    if abs(tN) < sys.float_info.epsilon:
        tc = 0.0
    else:
        tc = tN / tD
    
    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v) # =  S1(sc) - S2(tc)
    
    return np.linalg.norm(dP)

if __name__ == '__main__':
    Segment_S1 = np.array([[0.0,1.0,0.0],[1.0,0.0,0.0]])
    Segment_S2 = np.array([[0.0,0.0,0.5],[0.0,0.0,1.0]])
    dist = dist_3D_Segment_to_Segment(Segment_S1,Segment_S2)
    print(dist)