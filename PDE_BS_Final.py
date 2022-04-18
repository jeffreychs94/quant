import numpy as np
from numpy.linalg import inv
import scipy.sparse
import xlwings as xw
from math import *

####PDE
####https://towardsdatascience.com/option-pricing-using-the-black-scholes-model-without-the-formula-e5c002771e2f


#https://stackoverflow.com/questions/809362/how-to-calculate-cumulative-normal-distribution

def Norm_CDF(x):
    y = x / sqrt(2.0)
    return (1.0 + erf(y)) / 2.0

  
class BS_Inputs:
    
    def __init__(self,call_put,Eur_Amer,st0,K,vol,T_mat,r,q,st_min,st_max,stock_grid,time_grid):
                
            self.C_P = call_put
            self.Amer = Eur_Amer
            self.st0 = st0                      #Stock Spot Price
            self.K = K            
            self.sigma = vol              #Strike
            self.T_mat = T_mat
            self.r = r
            self.q = q
            self.st_min = st_min
            self.st_max = st_max
            self.time_grid = int(time_grid)
            self.stock_grid = int(stock_grid)
            
    def init_grid(self,vec_T,vec_S): ####Create Grid with Boundaries set
        
        #print("(PDE) Initiating Grid With Boundaries")
        
        Arr_OptVal = np.zeros((self.time_grid+1,self.stock_grid+1))
        
        if self.C_P == 1: ##Call
            #print("PDE Call")
            Arr_OptVal[:,-1] = self.st_max-self.K*np.exp(-self.r*(self.T_mat-vec_T))
            Arr_OptVal[:,0] = np.zeros(vec_T.shape)
            Arr_OptVal[-1,:] = np.maximum(vec_S-self.K,0) #...
        else:
            Arr_OptVal[:,-1] = np.zeros(vec_T.shape)
            Arr_OptVal[:,0] = np.full(vec_T.shape,self.K)
            Arr_OptVal[-1,:] = np.maximum(self.K - vec_S,0)
            
        return Arr_OptVal

    def compute_abc(self,vec_T,vec_S, dt, dS):
        
        #print("(PDE) Creating transition Matrix")
        # print(vec_T)
        # print(vec_S)
        
        vec_S = vec_S[1:-1]
        
        a = (self.r/2)*(vec_S/dS)*dt -0.5*self.sigma**2 * dt *((vec_S/dS)**2)
        b = 1 + (self.sigma**2) * ((vec_S/dS)**2) * dt + self.r * dt
        c = -(self.r/2)*(vec_S/dS)*dt -0.5*self.sigma**2 * dt *((vec_S/dS)**2)
        
        return scipy.sparse.diags( [a[1:],b,c[:-1]],offsets=[-1,0,1])
            
    def price_PDE_implicit(self):
        
        stockprice = self.st0
        dt = self.T_mat/self.time_grid 
        M = int(self.stock_grid)
        
        N = int(self.time_grid)
        
        # Max Stock Price set at return of 8 standard deviation away
 
        dS = (self.st_max-self.st_min)/self.stock_grid
        
        #Initiation
        vec_S = np.linspace(self.st_min,self.st_max,self.stock_grid+1) # M+1 and N+1 so that it starts with 0
        vec_T = np.linspace(0,self.T_mat,self.time_grid+1)
        arr_V = np.zeros((self.time_grid+1,self.stock_grid+1)) #...
        M_Stock = np.zeros((self.time_grid+1,self.stock_grid+1))
        
        V = self.init_grid(vec_T,vec_S)
        # print(V.shape)
        #print(V)
        InvLambda = inv(self.compute_abc(vec_T, vec_S, dt, dS).toarray())
        #print(InvLambda)
        #print("timegrid: ",self.time_grid,"StockGrid:",self.stock_grid)
        print(vec_S[1:M])
        print(vec_S[1:M]-self.K)
        print(np.maximum(vec_S[1:M],vec_S[1:M]-self.K))
        for i in range(self.time_grid,0,-1):
            # print(V[i-1,1:M])
            # print(V[i,1:M])
            # print(V[i,1:M].shape)
            # print(InvLambda.shape)
            
            V[i-1,1:M] = InvLambda.dot(V[i,1:M])
            
            ##American Feature
            if self.Amer:
                if self.C_P:
                    V[i-1,1:M] = np.maximum(V[i-1,1:M],vec_S[1:M]-self.K)
                else:
                    V[i-1,1:M] = np.maximum(V[i-1,1:M],self.K-vec_S[1:M])
                    
        # print(vec_S)
        # print(self.st0)
        Value_Index = np.where(vec_S == self.st0 )[0][0]
        # print(Value_Index)
        OptionValue = V[0,Value_Index]
        # print(V[0,Value_Index])
        return OptionValue

                
    def price_Formula(self):        
        d1 = (np.log(self.st0/self.K) + (self.r + self.sigma**2/2)*self.T_mat) / (self.sigma*np.sqrt(self.T_mat))
        d2 = d1 - self.sigma * np.sqrt(self.T_mat)
            
        if self.C_P: 
          return  self.st0 * Norm_CDF(d1) - self.K * np.exp(-self.r*self.T_mat)* Norm_CDF(d2)
        else:
          return self.K*np.exp(-self.r*self.T_mat)*Norm_CDF(-d2) - self.st0*Norm_CDF(-d1)
                              
@xw.func(ret='float64')    
def BS_OutputFunction(call_put,Eur_Amer,st0,K,vol,T_Mat,r,q,st_min,st_max,stock_grid,time_grid):

    Cls_BS = BS_Inputs(call_put,Eur_Amer,st0,K,vol,T_Mat,r,q,st_min,st_max,stock_grid,time_grid)        
    
    print("Black Scholes Formula: " , Cls_BS.price_Formula())
    print("PDE :", Cls_BS.price_PDE_implicit())
    
    return Cls_BS.price_PDE_implicit(),Cls_BS.price_Formula()

                    
print(BS_OutputFunction(1,1,526.40,600,0.15597,0.501369863013699,0.00176,0,0,2632.0,10,5))                
                