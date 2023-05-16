# -*- coding: utf-8 -*-
## Greedy Algorithm for NIRB
## Elise Grosjean
## 01/2021


from BasicTools.FE import FETools as FT
import numpy as np
from scipy import linalg

def orthogonality_check(Matrix,CorrelationMatrix):
    """
    This function check for the pairwise orthogonality of the new basis
    """
    list_ = list(Matrix)
    dot_matrix = np.array([[np.dot(CorrelationMatrix.dot(item1), item2.T) for item1 in list_] for item2 in list_])
    if (dot_matrix - np.eye(dot_matrix.shape[0]) < 1e-10).all():
        return True
    else:
        error = dot_matrix - np.eye(dot_matrix.shape[0])
        print("max error with identity: ",np.max(np.abs(error)))
        return False


def greedy_algorithm(snapshots,NumerOfTimesSteps,snapshotCorrelationOperator,h1ScalarProducMatrix, tolerance, NumberOfModes):
    """
    Compute the greedy reduced basis using a greedy algorithm.
    """
   
    TolGreedy=1e-6
    DegreesOfFreedom=np.shape(snapshotCorrelationOperator)[0]
    NbParam=len(snapshots)
    
    # initialization: First parameter
    ind=0 #first parameter (orthonormalized, so norm(snapshot) > 1e-10)

    ListeIndex=[ind] 
    basis_vectors=[] # basis functions list
    snapshotsList=snapshots.copy()
    snapshotsArray=snapshotsList.pop(ind) 
    
    residuals=np.array(snapshotsArray) #all time steps for first parameter 
    n,m = residuals.shape
            
    tol=1
    cptmodes=1
    cpt=0
    matVecProduct0=np.zeros((m,n))
    NormTestMax0=np.zeros(m)
    #for k in range(1):
    tolold=tol
    while(tol >= TolGreedy and NumberOfModes+1> cptmodes):
        matVecProduct= snapshotCorrelationOperator.dot(residuals.T)
        NormTestMax=np.sqrt(np.diag(np.dot(residuals,matVecProduct)))
        
        if cpt==0:
            matVecProduct0=matVecProduct 
            NormTestMax0=NormTestMax
            ind=np.argmax(NormTestMax0)
            tol=NormTestMax[ind]
            basis_vectors.append((np.reshape(residuals[ind,:],((m,1)))/NormTestMax[ind])) #first time index=random
            print(cptmodes," / ",NumberOfModes)
            print(NormTestMax[ind],"/",TolGreedy ,"(treshold)")
        else:
            TimeIndex= np.argmax(NormTestMax)
            tol=np.max(NormTestMax)
            print(cptmodes," / ",NumberOfModes, "(same parameter, new time step)")
            print(tol,"/",TolGreedy ,"(treshold)")

            basis_vectors.append(np.reshape(residuals[TimeIndex,:]/tol,((m,1))))
            if(tol > tolold):
                break;
        residuals -= np.outer(basis_vectors[-1],np.dot(basis_vectors[-1].T,matVecProduct0)).T
        cptmodes+=1
        cpt+=1
        
    # greedy on the parameters

    matrix=np.array(snapshotsList)
    cptParam=1
    
    while(NumberOfModes+1>cptmodes):
        residuals = np.reshape(matrix,(((NbParam-cptParam)*n,m)))
        
        matVecProduct0 = snapshotCorrelationOperator.dot(residuals.T)
        for i in range(len(basis_vectors)):
            residuals -= np.outer(basis_vectors[i],np.dot(basis_vectors[i].T,matVecProduct0)).T
           
        matVecProduct = snapshotCorrelationOperator.dot(residuals.T)
        NormTestMax = np.sqrt(np.diag(np.dot(residuals,matVecProduct)))
        Index = np.argmax(NormTestMax)
        IndexParam = int(Index/n)  
        tol =np.max(NormTestMax)
    
        basis_vectors.append(np.reshape(residuals[Index,:]/tol,((m,1))))
        print(cptmodes," / ",NumberOfModes, "(new parameter)")
        print(tol,"/",TolGreedy,"(treshold)")
        cptmodes+=1
        snapshotsArray=snapshotsList.pop(IndexParam)
        matrix=np.array(snapshotsList)
        residuals=np.array(snapshotsArray)
        cpt=0
        tol=1
        oldtol=tol
        while(tol >= TolGreedy and NumberOfModes >= cptmodes): #greedy on time steps
         
            if cpt==0:
                matVecProduct0=snapshotCorrelationOperator.dot(residuals.T)
                for i in range(len(basis_vectors)):
                    residuals -= np.outer(basis_vectors[i],np.dot(basis_vectors[i].T,matVecProduct0)).T
                   
            matVecProduct = snapshotCorrelationOperator.dot(residuals.T)
            NormTestMax = np.sqrt(np.diag(np.dot(residuals,matVecProduct)))
            TimeIndex = np.argmax(NormTestMax)
            tol = np.max(NormTestMax)
            if(tol > oldtol):
                break;
            basis_vectors.append(np.reshape(residuals[TimeIndex,:]/tol,((m,1))))
            residuals -= np.outer(basis_vectors[-1],np.dot(basis_vectors[-1].T,matVecProduct0)).T
            print(cptmodes," / ",NumberOfModes, "(same parameter, new time step)")
            print(tol,"/",TolGreedy,"(treshold)")
            cptmodes+=1
            cpt+=1
        
        cptParam+=1
    basis_vectors=np.reshape(basis_vectors,((NumberOfModes,m)))
    reducedOrderBasisU=np.array(np.reshape(basis_vectors,((NumberOfModes,m))))
    
    orthogonality=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    cpt=0
    while (orthogonality==False and 10>=cpt):       #Gram-schmidt 
        for i in range(1,NumberOfModes):
              
            reducedOrderBasisU[i,:]=basis_vectors[i]-sum((reducedOrderBasisU[k,:]*np.dot(snapshotCorrelationOperator.dot(reducedOrderBasisU[k]),basis_vectors[i]) for k in range(i)))            
            reducedOrderBasisU[i]/=np.sqrt(np.dot((snapshotCorrelationOperator.dot(reducedOrderBasisU[i])),reducedOrderBasisU[i]))
        for j in range(1,NumberOfModes):
            basis_vectors[j]=reducedOrderBasisU[j,:]
        cpt+=1
        orthogonality=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
        print(orthogonality)

       ### H1 Orthogonalization
    
    K=np.zeros((NumberOfModes,NumberOfModes)) #rigidity matrix
    M=np.zeros((NumberOfModes,NumberOfModes)) #mass matrix
    for i in range(NumberOfModes):
        matVecH1=h1ScalarProducMatrix.dot(reducedOrderBasisU[i,:])
        matVecL2=snapshotCorrelationOperator.dot(reducedOrderBasisU[i,:])
        for j in range(NumberOfModes):
            if i>=j:
            
                K[i,j]=np.dot(matVecH1,reducedOrderBasisU[j,:])
                M[i,j]=np.dot(matVecL2,reducedOrderBasisU[j,:])
                K[j,i]=K[i,j]
                M[j,i]=M[i,j]
    
    
    # on resoud Kv=lambd Mv
    eigenValues,vr=linalg.eig(K, b=M) #eigenvalues + right eigenvectors
    idx = eigenValues.argsort()[::-1]
    eigenValues = eigenValues[idx]
    print("EV : ",eigenValues)
    print("sqrt(EV) : ",np.sqrt(eigenValues))
    eigenVectors = vr[:, idx]
    reducedOrderBasisU=np.dot(eigenVectors.transpose(),reducedOrderBasisU)

    for i in range(NumberOfModes):
        reducedOrderBasisNorm=np.sqrt(reducedOrderBasisU[i,:]@(snapshotCorrelationOperator@reducedOrderBasisU[i,:]))
        reducedOrderBasisU[i,:]/=reducedOrderBasisNorm#np.sqrt(M[i,i]) #L2 orthonormalization
    
    orthogonality=orthogonality_check(reducedOrderBasisU,snapshotCorrelationOperator)
    print(orthogonality)
    return reducedOrderBasisU



