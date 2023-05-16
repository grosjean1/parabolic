# -*- coding: utf-8 -*-
## NIRB parabolic test with OFFLINE/ONLINE DECOMPOSITION

## Elise Grosjean
## 01/2022



import numpy as np
import sys
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import os
import os.path as osp

import Readers as MR
import SolutionVTKWriter as SVTKW

from BasicTools.FE import FETools as FT
#import pickle

import Greedy as GD
#import SVD

from scipy import linalg
from scipy import interpolate
from scipy.interpolate import griddata 
from scipy.interpolate import interp1d #time interpolation

from scipy.spatial import cKDTree #space interpolation
from scipy.sparse import coo_matrix
        
############################################################
"""          Initialization                              """

############################################################

onlineParameter='9'

## Directories
currentFolder=os.getcwd()
dataFolder=currentFolder
FinedataFolder=osp.join(dataFolder,'FineSnapshots/'+sys.argv[2]) #for fine snapshotsS
CoarsedataFolder=osp.join(dataFolder,'CoarseSnapshots/'+sys.argv[3]+'/'+sys.argv[4]) #for coarse snapshots
print("fine folder: ", FinedataFolder)
print("coarse folder: ", CoarsedataFolder)


ns=0 #number of snapshots
count1=0
for _, folder, _ in os.walk(FinedataFolder): #number of folders in FineData
    count1 += len(folder)

ns=count1-1 #-1 because of the online parameter not included in offline snapshots
print("Number of snapshots: ",ns)
#ns=18


nev=int(sys.argv[1])   #nombre de modes
print("nev: ",nev)

time=0.0 #init
dimension=2 #2D
           
TF=int(len([name for name in os.listdir(FinedataFolder+"/1/")]))
print("Number of fine time steps: ",TF)

TG=int(len([name for name in os.listdir(CoarsedataFolder+"/1/")]))
print("Number of coarse time steps: ",TG)

dtF=float(sys.argv[2]) #fine time steps
dtG=float(sys.argv[4]) #coarse time steps

for time in np.arange(0, 1.0001, dtF):
    if time>=0: #.9999: #NIRB sur t0=]1, 2]
        t0f=time
        break
for time in np.arange(0, 1.0001, dtG):
    if time>=0: #.9999:
        t0g=time
        break
    
#for time interpolation
oldtime=np.arange(t0g, 1.0001, dtG)
newtime=np.arange(t0f, 1.0001, dtF)

print(np.shape(oldtime),len(newtime)," ", TF," ",TG)
print(oldtime," ",newtime)
"""
-------------------
###  Read fine mesh
------------------- 
"""
meshFileName = FinedataFolder + "/1/Snapshoth_1.vtu";
mesh=MR.Readmesh(meshFileName)
mesh.nodes= mesh.nodes[:,:2] #2D

print("Fine mesh defined in " + meshFileName + " has been read")
nbeOfComponentsPrimal = 1 # 1 field 
numberOfNodes = mesh.GetNumberOfNodes()
print("DoF fine mesh ", numberOfNodes)
"""
-------------------
###  Read coarse mesh
------------------- 
"""

meshFileName2 = CoarsedataFolder + "/1/Snapshoth_1.vtu";#//FineMesh/mesh1.msh"
mesh2=MR.Readmesh(meshFileName2)
mesh2.nodes = mesh2.nodes[:,:2] #CAS 2D

print("Coarse mesh defined in " + meshFileName2 + " has been read")
numberOfNodes2 = mesh2.GetNumberOfNodes()
print("DoF coarse mesh ", numberOfNodes2)

"""
-------------------
###  mesh space interpolation
------------------- 
"""

inputnodes=mesh2.nodes
outputnodes=mesh.nodes
"""
kdt = cKDTree(inputnodes)
nbtp = outputnodes.shape[0]
_, ids = kdt.query(outputnodes)
cols=ids
row = np.arange(nbtp)
data = np.ones(nbtp)
Operator=coo_matrix((data, (row, cols)), shape=(nbtp , inputnodes.shape[0]))
"""
"""
-------------------
###  read all snapshots ...
------------------- 
"""

parameters=[]
for i in range(1,ns+2):
    if(float(0.5*i)%1>1e-3):
        parameters.append(str(0.5*i))
    else:
        parameters.append(str(int(0.5*i)))

print("param:",parameters)
parameters.remove(onlineParameter) #online parameter mu=1
print("parameters :",parameters)


snapshots=[]

for e,i in enumerate(parameters):
    snapshotsTime=[]
    for time in range(0,TF):   
        snapshot =MR.VTKReadToNp("Velocity",FinedataFolder+"/"+i+"/Snapshoth_",time)
        snapshotsTime.append(snapshot)
    snapshots.append(snapshotsTime)

    
snapshotsH=[]

for e,i in enumerate(parameters):
    snapshotsHTime=[]
    for time in range(0,TG):

        snapshotH =MR.VTKReadToNp("Velocity",CoarsedataFolder+"/"+i+"/Snapshoth_",time)

        #Compute the projected data using the projection operator
        snapshotHSpaceinterpolated = griddata(inputnodes,snapshotH,outputnodes,method='linear')#//Operator.dot(snapshotH)
        snapshotsHTime.append(snapshotHSpaceinterpolated)
    print(i," ",len(oldtime),len(snapshotsHTime))
    interp  = interp1d(oldtime,snapshotsHTime,kind='quadratic',axis=0,fill_value="extrapolate")

    solutionUHI=interp(newtime) #time and space interpolation
    snapshotsH.append(solutionUHI)
  
    
snapshotsHarraytl=np.array(snapshotsH)
#print("Temps,param ",np.shape(snapshotsHarraytl))


############################################################
"""          Greedy                                      """
############################################################

print("ComputeL2ScalarProducMatrix ...")

l2ScalarProducMatrix = FT.ComputeL2ScalarProducMatrix( mesh, nbeOfComponentsPrimal)
h1ScalarProducMatrix = FT.ComputeH10ScalarProductMatrix(mesh, nbeOfComponentsPrimal)

##### ALGO (full) GREEDY
#reducedOrderBasisU,_=GD.Greedy(snapshots,TF,l2ScalarProducMatrix,h1ScalarProducMatrix=None,NumberOfModes=nev,Tol=1e-6)

reducedOrderBasisU=GD.greedy_algorithm(snapshots,TF,l2ScalarProducMatrix,h1ScalarProducMatrix,1e-6,nev)

print("Number of modes after greedy",nev)

############################################################
"""          Rectification                               """
############################################################

RI=np.zeros((TF,nev,nev)) #global rectification over the time steps
for time in range(TF):

    alpha=np.zeros((nev,ns)) #fine coefficients
    beta=np.zeros((ns,nev)) #coarse coefficients
    betainv=np.zeros((nev,ns)) 
    R=np.zeros((nev,nev)) #rectification matrix

    for j,elt in enumerate(parameters):

        u1PT = snapshots[j][time]
        u1T = snapshotsH[j][time]
        
        for i in range(nev):
            alpha[i,j]=u1PT@(l2ScalarProducMatrix@reducedOrderBasisU[i,:])
            beta[j,i]=u1T@(l2ScalarProducMatrix@reducedOrderBasisU[i,:])
            betainv[i,j]=beta[j,i]
        
    lambd=1e-13
    Rns=np.zeros((nev,ns))
    
    for j in range(ns):
        Rns[:,j]=np.linalg.inv(beta.transpose()@beta+lambd*np.eye(nev))@betainv[:,j]
        
    for i in range(nev):
        R[i,:]=Rns@alpha[i,:]

    RI[time,:,:]=R

############################################################
"""          Online part                                 """
############################################################

for i in [onlineParameter]:
    
    snapshotsHTime=[]
    
    for time in range(0,TG):
        snapshotH =MR.VTKReadToNp("Velocity",CoarsedataFolder+"/"+i+"/Snapshoth_",time)
   
        #Compute the projected data using the projection operator
        snapshotHSpaceinterpolated = snapshotHSpaceinterpolated = griddata(inputnodes,snapshotH,outputnodes,method='linear')#//Operator.dot(snapshotH)
        snapshotsHTime.append(snapshotHSpaceinterpolated)
    
    interp  = interp1d(oldtime,snapshotsHTime,kind='quadratic',axis=0,fill_value="extrapolate")
    solutionUHI=interp(newtime) #time and space interpolation

    for time in range(0,TF):
        
        R=RI[time,:,:]
        #u1PT=solutionUHI[time]
        u1PT=MR.VTKReadToNp("Velocity",FinedataFolder+"/"+str(i)+"/Snapshoth_",time)
        coef=np.zeros(nev)
        CompressedSolutionUj=np.zeros(nev)
        for k in range(nev):
            coef[k]=0
            CompressedSolutionUj[k]=u1PT@(l2ScalarProducMatrix@reducedOrderBasisU[k,:])
            #print("k: ", k, CompressedSolutionUj[k])
            for j in range(nev):
                coef[k]+=R[k,j]*u1PT@(l2ScalarProducMatrix@reducedOrderBasisU[j,:])#CompressedSolutionUj[j]
            #print(coef[k])
        
        #reconstructedCompressedSolution =MR.VTKReadToNp("Velocity",CoarsedataFolder+"/"+str(i)+"/Snapshoth_",time)
        reconstructedCompressedSolution = np.dot(coef, reducedOrderBasisU) #rectified nirb
        #reconstructedCompressedSolution = np.dot(CompressedSolutionUj, reducedOrderBasisU) #rectified nirb
        
        ##################################################
        #######   saving solution in VTK ############
        VTKBase = MR.VTKReadmesh(meshFileName)
        SVTKW.numpyToVTKWrite(VTKBase,reconstructedCompressedSolution,"NIRB_approximation_"+str(time)+"_"+str(nev)+".vtu")
        ##################################################

        
