# -*- coding: utf-8 -*-
"""
---------------- OPTIMAL SENSOR AND ACTUATOR PLACEMENT -----------------------
-----------------------------ISOTROPIC PLATE----------------------------------
Created on Wed May 13 19:46:00 2020

@author: Sergio Cantero Chinchilla
"""
from IPython import get_ipython
get_ipython().magic('reset -sf') 
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Palatino']})
rc('text', usetex=True)

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator
from matplotlib import path
from pathlib import Path

# Function: Generate prior samples for isotropic materials
def priorSmpl(N_smpl,X_d_v,Y_d_v,V_mean,V_std):
###############################################################################
#   Input paramters:
#       N_smpl          -   Number of model parameter samples to evaluate
#                           the objective function
#       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
#                           area of possible damage occurrence
#       V_mean, V_std   -   Mean and standard deviation of the prior PDF of
#                           the wave propagation velocity
#
#   Return (save):
#       prior_smpl      -   Set of N_smpl samples obtained from the prior
#                           distribution of the model parameters
###############################################################################
    xd=np.zeros(N_smpl)
    yd=np.zeros(N_smpl)
    aux_vect=np.array([X_d_v.transpose(),Y_d_v.transpose()]).transpose()
    p = path.Path(aux_vect)
    for i in np.arange(0,N_smpl):
        flagIsIn=0
        while np.logical_not(flagIsIn):
            xd[i]=np.min(X_d_v) + np.multiply((np.max(X_d_v) - \
              np.min(X_d_v)),np.random.rand(1))
            yd[i]=np.min(Y_d_v) + np.multiply((np.max(Y_d_v) - \
              np.min(Y_d_v)),np.random.rand(1))
            flagIsIn=p.contains_points(np.array([xd[i],yd[i]]).reshape(1,2))
            
    V=np.random.normal(V_mean, V_std, N_smpl)
    
    prior_smpl=np.squeeze(np.array([[xd],[yd],[V]]))
    
    np.save('./dat/prior_smpl', prior_smpl)
    
# Function: Generate grid of possible sensor/actuator locations
def gridPoints(nSubGrid,minDist,X_d_v,Y_d_v,X_g_v,Y_g_v):
###############################################################################
#   Input paramters:
#       nSubGrid        -   Number of rings around the area of possible 
#                           damage occurrence for the grid of points
#       minDist         -   Approximate distance between consecutive grid
#                           points (possible sensor/actuator locations)
#       X_d_v, Y_d_v    -   X and Y coordinates of the vertices of the
#                           area of possible damage occurrence
#       X_g_v, Y_g_v    -   X and Y coordinates of the vertices of the
#                           outer geometry of the plate
#
#   Return (save):
#       GridPtX,GridPtY -   X and Y coordinates of the grid points of
#                           possible sensor/actuator locations
###############################################################################
    # Create the corner points
    auxGrdPtX=np.zeros((nSubGrid,len(X_d_v)+1))
    auxGrdPtY=np.zeros((nSubGrid,len(X_d_v)+1))
    for i in np.arange(0,nSubGrid):
        for j in np.arange(0,len(X_d_v)):
            auxGrdPtX[i,j]=X_g_v[j]-(i+1)*(X_g_v[j]-X_d_v[j])/(nSubGrid+1)
            auxGrdPtY[i,j]=Y_g_v[j]-(i+1)*(Y_g_v[j]-Y_d_v[j])/(nSubGrid+1)
        auxGrdPtX[i,j+1]=auxGrdPtX[i,0]
        auxGrdPtY[i,j+1]=auxGrdPtY[i,0]
        
    # Generate the grid points with approx minDist between them:
    GridPtX=np.array([])
    GridPtY=np.array([])
    for j in np.arange(0,nSubGrid):
        for i in np.arange(0,len(X_d_v)-1):
            
            auxDist=np.sqrt((auxGrdPtX[j,i+1]-auxGrdPtX[j,i])**2+\
                            (auxGrdPtY[j,i+1]-auxGrdPtY[j,i])**2)
            auxNpt=np.around(auxDist/minDist)
            
            if auxGrdPtX[j,i]==auxGrdPtX[j,i+1]:
                GrdPtXtemp=np.linspace(auxGrdPtX[j,i],auxGrdPtX[j,i+1],auxNpt.astype(int))
            else:
                GrdPtXtemp=np.arange(auxGrdPtX[j,i],auxGrdPtX[j,i+1],\
                                 (auxGrdPtX[j,i+1]-auxGrdPtX[j,i])/auxNpt)
            
            if auxGrdPtY[j,i]==auxGrdPtY[j,i+1]:
                GrdPtYtemp=np.linspace(auxGrdPtY[j,i],auxGrdPtY[j,i+1],auxNpt.astype(int))
            else:
                GrdPtYtemp=np.arange(auxGrdPtY[j,i],auxGrdPtY[j,i+1],\
                                 (auxGrdPtY[j,i+1]-auxGrdPtY[j,i])/auxNpt)
            
            GridPtX=np.hstack((GridPtX,GrdPtXtemp))
            GridPtY=np.hstack((GridPtY,GrdPtYtemp))
            
    np.savez('./dat/GridPoints', GridPtX=GridPtX, GridPtY=GridPtY)
    
# Function: P matrices for isotropic materials
def P_fun_py(xs, ys, xa, ya, xd, yd, V):
###############################################################################
#   Input paramters:
#       xs, ys          -   X and Y coordinates of the sensor location
#       xa, ya          -   X and Y coordinates of the actuator location
#       xd, yd          -   X and Y coordinates of the damage location
#       V               -   Wave propagation velocity
#
#   Return:
#       Multiplication of gradient of the time-of-flight model in a matrix form
###############################################################################
    return np.array([[(- ((xa - xd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) -(xs - xd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2)))  ** 2,\
        (- ((xa - xd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) -(xs - xd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- ((ya - yd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (ys - yd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))),\
        (- ((xa - xd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) -(xs - xd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- (np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2) / V ** 2) - np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2) / V ** 2)],\
        [(- ((xa - xd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (xs - xd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- ((ya - yd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (ys - yd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))),\
        (- ((ya - yd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (ys - yd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) ** 2,\
         (- ((ya - yd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (ys - yd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- (np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2) / V ** 2) - np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2) / V ** 2)],\
         [(- ((xa - xd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) -(xs - xd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- (np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2) / V ** 2) - np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2) / V ** 2),\
         (- ((ya - yd) / (V*np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2))) - (ys - yd) / (V*np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2))) * (- (np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2) / V ** 2) - np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2) / V ** 2),\
        (- (np.sqrt((xa - xd) ** 2 + (ya - yd) ** 2) / V ** 2) - np.sqrt((xs - xd) ** 2 + (ys - yd) ** 2) / V ** 2) ** 2]])

# Function: Evaluation of P matrices
def P_eval(GridPtX,GridPtY,prior_smpl,N_smpl):
###############################################################################
#   Input paramters:
#       GridPtX,GridPtY -   X and Y coordinates of the grid points of
#                           possible sensor/actuator locations
#       prior_smpl      -   Set of N_smpl samples obtained from the prior
#                           distribution of the model parameters
#       N_smpl          -   Number of model parameter samples to evaluate
#                           the objective function
#
#   Return (save):
#       P_sns_eval      -   Model-gradient multiplication matrices
###############################################################################
    NmaxSns=len(GridPtX)
    NmaxAct=len(GridPtX)
    P_sns_eval=P_fun_py(np.array([np.tile(GridPtX[j],(NmaxAct,N_smpl)) for j in np.arange(NmaxSns)]), np.array([np.tile(GridPtY[j],(NmaxAct,N_smpl)) for j in np.arange(NmaxSns)]), np.tile(np.array([np.repeat(GridPtX[l],N_smpl) for l in np.arange(NmaxAct)]),(NmaxSns,1,1)), np.tile(np.array([np.repeat(GridPtY[l],N_smpl) for l in np.arange(NmaxAct)]),(NmaxSns,1,1)), np.tile(prior_smpl[0,:],(NmaxSns,NmaxAct,1)), np.tile(prior_smpl[1,:],(NmaxSns,NmaxAct,1)), np.tile(prior_smpl[2,:],(NmaxSns,NmaxAct,1)))
    np.save('./dat/P_sns_eval', P_sns_eval)

# Function: Entropy-based convex objective function
def objFun(candSol,NmaxSns,NmaxAct,N_smpl,Nvar,x_cost,y_cost,P_sns_eval):
###############################################################################
#   Input paramters:
#       candSol         -   candidate solution
#       NmaxAct         -   Maximum number of actuators
#       NmaxSns         -   Maximum number of sensors
#       N_smpl          -   Number of model parameter samples to evaluate
#                           the objective function
#       N_var           -   Number of model parameters
#       x_cost, y_cost  -   X and Y coordinates of the interpolating points
#                           used to build the cost function (by pchip())
#       P_sns_eval      -   Model-gradient multiplication matrices
#
#   Output parameters:
#       f               -   Objective function evaluation
#       g               -   Gradient eval. of the the objective function
###############################################################################
    z=np.array(candSol[0:NmaxAct]) # Vector for the sensors - z
    w=np.array(candSol[NmaxAct:2*NmaxAct]) # Vector for the actuators - w
    n=np.around(candSol[-1],decimals=4) # Number of sensors + actuators
    
    f_aux=np.zeros(N_smpl)
    g_aux1=np.zeros((NmaxAct,N_smpl))
    g_aux2=np.zeros((NmaxSns,N_smpl))
    for i in np.arange(0,N_smpl):
        P_aux=np.array(P_sns_eval[:,:,:,:,i])
        Q=np.sum(np.reshape((z*np.reshape(np.transpose(np.sum(np.reshape((w*np.reshape(np.transpose(P_aux,(2,0,1,3)),(NmaxAct,Nvar*Nvar*NmaxSns)).T).T,(NmaxAct,Nvar,Nvar,NmaxSns)),0),(2,0,1)),(NmaxSns,Nvar*Nvar)).T).T,(NmaxAct,Nvar,Nvar)),0)
        f_aux[i]=np.log(np.linalg.det(Q))
        # For loops to address the gradient (Jacobian)
        for l in np.arange(0,NmaxAct):
            g_aux1[l,i]=np.trace(np.linalg.solve(Q,np.sum(np.transpose(np.reshape((w*np.reshape(np.transpose(P_aux[:,:,:,l],[2,0,1]),(NmaxSns,Nvar*Nvar)).T).T,(NmaxSns,Nvar,Nvar)),[1,2,0]),2)))
        for j in np.arange(0,NmaxSns):
            g_aux2[j,i]=np.trace(np.linalg.solve(Q,np.sum(np.transpose(np.reshape((z* np.reshape(np.transpose(P_aux[:,:,j,:],[2,0,1]),(NmaxSns,Nvar*Nvar)).T).T,(NmaxSns,Nvar,Nvar)),[1,2,0]),2)))
    
    # Cost function evaluation
    costPchip=PchipInterpolator(x_cost,y_cost)
    g_eval=costPchip(np.array([n]))
    
    # Monte Carlo approximation of the objective function:
    f_MC = -1/N_smpl * np.sum(f_aux)
    f = f_MC + np.abs(f_MC) * g_eval
    
    # Derivative of the cost function
    grdPchip=np.arange(n-1,n+1,0.01)
    dgAux=costPchip(grdPchip)
    slp= np.diff(dgAux) / np.diff(grdPchip)
    slpGrd=grdPchip[:-1]
    dgEval=np.interp(n,slpGrd,slp)
    
    # Jacobian (gradient) 
    g = np.array([-1/N_smpl * np.sum(g_aux1,1), -1/N_smpl * np.sum(g_aux2,1)]) * \
        (1+g_eval*f_MC/np.abs(f_MC))
    g = np.append(g, np.abs(f_MC)*dgEval)
    return f, g

# Function: Table formatting
def print_table(table, text_file):
###############################################################################
#   Input paramters:
#       table           -   List with the headings and optimal solution
#       text_file       -   Text file where the solution is written
#
#   Return (print):
#       Write formatted table with the optimal solution
###############################################################################
    longest_cols = [
        (max([len(str(row[i])) for row in table]) + 3)
        for i in range(len(table[0]))
    ]
    row_format = "".join(["|{:^" + str(longest_col) + "}|" for longest_col in longest_cols])
    for row in table:
        print(row_format.format(*row), file=text_file)


"""
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        INPUT PARAMETERS FOR THE SOFTWARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""  
# If you are going to change any of the input data, remove the files in /dat
###############################################################################
# 00 - Geometry of the plate
X_g_v=[- 1  ,   1  ,   1   , - 0.9, - 1   ]
Y_g_v=[- 0.5, - 0.5, - 0.15,   0.5, - 0.5 ]
###############################################################################
# 01 - Prior samples of the ToF model parameters (X_d,Y_d, and V)
# X_g_v and X_d_v need to have the same number of vertices
N_smpl = 10 # Number of samples for the MC integrals
X_d_v=np.array([- 0.7 ,  0.7 ,  0.7  ,- 0.65,- 0.7 ])
Y_d_v=np.array([- 0.35,- 0.35,- 0.23 ,  0.22,- 0.35])
V_mean, V_std= 2800 , 60

# Calculate or load the samples from the prior PDF
my_file = Path("./dat/prior_smpl.npy")
if my_file.is_file()==0:
    priorSmpl(N_smpl,X_d_v,Y_d_v,V_mean,V_std)
    prior_smpl=np.load(my_file)
else:
    prior_smpl=np.load(my_file)

###############################################################################
# 02 - Grid of possible sensor locations
nSubGrid=2 # Number of concentric grids

minDist=0.2 # Distance between grid points

# Calculate or load the grid of possible sensor locations
my_file = Path("./dat/GridPoints.npz")
if my_file.is_file()==0:
    gridPoints(nSubGrid,minDist,X_d_v,Y_d_v,X_g_v,Y_g_v)
    auxFilez=np.load(my_file)
    GridPtX=auxFilez['GridPtX']
    GridPtY=auxFilez['GridPtY']
else:
    auxFilez=np.load(my_file)
    GridPtX=auxFilez['GridPtX']
    GridPtY=auxFilez['GridPtY']

###############################################################################
# 03 - Evaluation of the P matrices - Isotropic material

# Calculate or load the grid of possible sensor locations
my_file = Path("./dat/P_sns_eval.npy")
if my_file.is_file()==0:
    P_eval(GridPtX,GridPtY,prior_smpl,N_smpl)
    P_sns_eval=np.load(my_file)
else:
    P_sns_eval=np.load(my_file)
###############################################################################
# 04 - Definition of the objective function

NmaxSns=len(GridPtX)     # No sensors
# N_smpl                 # No samples
NmaxAct=len(GridPtX)     # No actuators
Nvar=P_sns_eval.shape[1] # No variables

# Cost function definition by interpolating points (for Pchip) in [0,1]
x_cost=np.array([0  , 30 , 55, 60])
y_cost=np.array([0.0, 0.3, 0.95,1 ])

# Generate the inequality linear constraints
A=np.zeros((2*(NmaxAct+NmaxSns+1),NmaxAct+NmaxSns+1))
for i in np.arange(0,NmaxAct+NmaxSns+1):
    for j in np.arange(0,NmaxAct+NmaxSns+1):
        if i==j:
            A[i,j]=1
for i in np.arange(NmaxAct+NmaxSns+1,2*(NmaxAct+NmaxSns+1)):
    auxA=i-(NmaxAct+NmaxSns+1)
    for j in np.arange(0,NmaxAct+NmaxSns+1):
        if auxA==j:
            A[i,j]=-1
A[:,-1]=0
A[-1,:]=0
b=np.concatenate((np.zeros(NmaxAct+NmaxSns+1), -1* np.ones(NmaxAct+NmaxSns+1)), axis=None)
b[-1]=0


# Inequality contraints for the number of sensors
Aeq = np.concatenate((np.ones(NmaxAct+NmaxSns), -1), axis=None);
beq = np.array([1e-5])

# Both constraints
cons=[{'type': 'ineq', 'fun': lambda x: A @ x - b},
      {'type': 'ineq', 'fun': lambda x: -1 * Aeq @ x + beq}]

# Initial point by number of sensors and actuators
Nini=12
candSol = np.concatenate((np.ones(NmaxAct)/NmaxAct*Nini/2,\
                          np.ones(NmaxSns)/NmaxSns*Nini/2,\
                          Nini), axis=None)

# Minimization of the objective function considering the gradient (Jacobian)
print("COST-BENEFIT ANALYSIS")
res = minimize(objFun,candSol,\
               args=(NmaxSns,NmaxAct,N_smpl,Nvar,x_cost,y_cost,P_sns_eval),\
               jac=True,constraints=cons, tol=1e-7,\
               options={'disp': True, 'iprint':3, 'maxiter':100})

np.save('./res/Opt_Sol', res)
###############################################################################
# 05 - Plot the prior samples of damage locations
fig=plt.figure(1)
for i in np.arange(0,len(X_d_v)-1):
    plt.plot(np.array([X_d_v[i],X_d_v[i+1]]),np.array([Y_d_v[i],Y_d_v[i+1]]),\
         color = (0.3,0.3,0.3), label='_nolegend_')
    plt.plot(np.array([X_g_v[i],X_g_v[i+1]]),np.array([Y_g_v[i],Y_g_v[i+1]]),\
         color = (0.,0.,0.), label='_nolegend_')
plt.plot(prior_smpl[0,],prior_smpl[1,],'.',color = (0.7,0.7,0.7), label='_nolegend_')
plt.plot(GridPtX,GridPtY,'.',alpha=.3, label='_nolegend_')
plt.scatter(GridPtX,GridPtY,s=res.x[:(NmaxSns)]*100,alpha=1,marker=8,label='Sensors')
plt.scatter(GridPtX,GridPtY,s=res.x[(NmaxSns):-1]*100,alpha=1,marker=9,label='Actuators')
plt.ylabel('Y coordinate [m]')
plt.xlabel('X coordinate [m]')
plt.legend()
plt.show()
fig.savefig("./res/OptSol.pdf", bbox_inches='tight')
###############################################################################
# 06 - Table for the output
tabOutp=[["PZT #", "X", "Y", "z", "w"],\
         ["------", "------", "------", "------", "------", "------"]]

n_opt=np.around(res.x[-1],0).astype(int)
z=np.array(res.x[0:NmaxAct]) # Vector for the sensors
w=np.array(res.x[NmaxAct:2*NmaxAct]) # Vector for the actuators
maxRow=np.around(np.max([n_opt,len(z[z>.2]),len(w[w>.2])])/2+1).astype(int)

iDx=np.argsort(-z)
nDecil=3
for i in np.arange(0,maxRow):
    tabOutp.append([i,np.around(GridPtX[iDx[i]],nDecil),np.around(GridPtY[iDx[i]],nDecil),\
                    np.around(z[iDx[i]],nDecil), np.around(w[iDx[i]],nDecil)])
with open("./res/Output.txt", 'w') as text_file :
    print('ENTROPY-BASED CONVEX OPTIMIZATION - ISOTROPIC PLATE\nOptimal sensor and actuator solution\n',file=text_file)
    print_table(tabOutp,text_file)    
    print('\nNOMENCLATURE',file=text_file)
    print('X   - X coordinate of the sensor/actuator [m]',file=text_file)
    print('Y   - Y coordinate of the sensor/actuator [m]',file=text_file)
    print('z   - Continuous decision variable for sensors',file=text_file)
    print('w   - Continuous decision variable for actuators\n',file=text_file)
    print('--------------------------------------------------\n',file=text_file)
    print('\nThe optimal number of sensors and actuators is: %.2f' % res.x[-1],\
          file=text_file)
###############################################################################
# 07 - Flag when finished
print("Finished")
###############################################################################
get_ipython().magic('reset -sf') 
