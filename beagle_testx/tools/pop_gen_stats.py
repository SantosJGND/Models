from sklearn.metrics import pairwise_distances
import numpy as np

## GGVE p. 62 eq. 2.38
def PI_Taj(dataM):
    
    n= dataM.shape[0]
    pairDiff= pairwise_distances(dataM, metric='manhattan')
    mu= np.tril_indices(n)
    pair_Diff= sum(pairDiff[mu])
    
    pi= (2 / n / (n-1)) * pair_Diff
    
    return pi

## GGVE p. 62 eq. 2.40
def Watt_est(dataM, L= 0):
    
    Sn= np.sum(dataM,axis= 0)
    Sn= [x for x in Sn if x > 0]
    Sn= len(Sn)
    
    An= [1 / j for j in range(1,dataM.shape[0])]
    An= sum(An)
    
    Tw= Sn / An
    
    return Tw

## GGVE p. 62 eq. 2.42
def TajD(dataM):
    
    n= dataM.shape[0]
    
    Pit= PI_Taj(dataM)
    Watt= Watt_est(dataM)
    
    ##
    Sn= np.sum(dataM,axis= 0)
    Sn= [x for x in Sn if x >= 0]
    Sn= len(Sn)
    
    An= [1 / j for j in range(1,n)]
    An= sum(An)
    
    Bn= [1 / j**2 for j in range(1,n)]
    Bn= sum(Bn)
    
    E1= (n + 1) / (3 * An * (n-1)) - 1 / An**2
    E2= 2 * (n**2 + n + 3) / 9*n / (n-1)
    E2= E2 - (n+2)/ n / An + Bn / An**2
    
    E2= E2 / (An**2 + Bn)
    
    ##
    denom= E1 * Sn + E2 * Sn * (Sn - 1)
    denom= np.sqrt(denom)
    
    D= Pit - Watt
    D= D / denom
    
    return D