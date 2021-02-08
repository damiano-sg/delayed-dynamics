import scipy as sp
import numpy as np
import numpy.random as rd
import scipy.linalg as LA
import scipy.sparse.linalg as sLA

#CLASS THAT GENERATE THE RANDOM MATRIX

class RMG :
    def __init__(self, S,C,d,mu,sigma,pdf='normal'):
    	#size of the matrix
        self.S = S
        #connectance of the matrix (non-zero entries)
        self.C = C
        #diagonal entries
        self.d = d
        #mean val of the entries
        self.mu = mu
        #st. dev. of the entries
        self.sigma = sigma
        #type of pdf: normal, uniform, multi_norm (multivariate, mu is a vector and sigma a cov matrix)
        self.pdf = pdf
        #the matrix
        self.M = 0
    
    #INITIALIZATIONS OF RANDOM MATRICES
    def uniform (self):
        return rd.uniform(0,1,(self.S,self.S))
    def initialize(self):
        if self.pdf=='normal' :
            return rd.normal(self.mu,self.sigma,(self.S,self.S))
        if self.pdf=='uniform' :
            return rd.uniform(self.mu-np.sqrt(3)*self.sigma,self.mu+np.sqrt(3)*self.sigma,(self.S,self.S)) 
        if self.pdf=='multi_norm' :
            m = rd.multivariate_normal(self.mu, self.sigma,(self.S,self.S)) #mu: 2 and sigma(cov mat): 2x2
            return np.triu(m[:,:,0],1) + np.triu(m[:,:,1],1).T
        #this is a slower method for multivariate when the size S is large
        if self.pdf=='multi_norm_for_cycle' :                    
            m = rd.multivariate_normal(self.mu, self.sigma,int(self.S*(self.S-1)/2))
            mij, mji = m[:,0], m[:,1]
            k=0
            mm = np.zeros((self.S,self.S))
            for i in range(self.S):
                for j in range(i+1,self.S):
                    mm[i,j], mm[j,i] = mij[k], mji[k]
                    k+=1
            return mm
        
    #RANDOM MATRICES GENERATORS
    #PURE RANDOM       
    def ran_matrix (self):
        c = self.uniform() <= self.C
        self.M = self.initialize()*c-self.d*np.eye(self.S)
    #SYMMETRIC
    def sym_matrix (self):
        self.ran_matrix()
        self.M = 0.5*(self.M+(self.M).T)
    #PREDATOR-PREY TYPE
    def pp_matrix (self):
        m = self.initialize()
        m = abs(np.triu(m,0))-abs(np.tril(m,-1))
        t = self.uniform()
        t = np.triu(t,1)+np.triu(t,1).T <= 0.5 
        m = m*(2*t-1) 
        u = self.uniform()
        c = (np.triu(u,0)+np.triu(u,1).T) <= self.C
        self.M = m*c-self.d*np.eye(self.S)
    #MIX TYPE (COOPERATIVE-COMPETITIVE)
    def mix_matrix(self):
        m = abs(self.initialize())
        t = self.uniform()
        t = np.triu(t,1)+np.triu(t,1).T <= 0.5 
        m = m*(2*t-1)
        u = self.uniform()
        c = np.triu(u,0)+np.triu(u,1).T <= self.C
        self.M = m*c-self.d*np.eye(self.S)
    #MULTIVARIATE
    def multi_matrix(self):
        m = self.initialize()
        u = self.uniform()
        c = np.triu(u,0)+np.triu(u,1).T <= self.C
        self.M = m*c-self.d*np.eye(self.S)

    # GENERATE DIRECTLY THE EIGENVALUES USING THEIR DISTRIBUTION (it doesn't work for multivariate distributions)
    def generate_eigenvalues(self, mat):
        #here I'm using a uniform distribution of eigenvalues (OK ONLY AT FIRST APPROX FOR LARGE S -> circular law for non-symm matrices)
        x,y = rd.uniform(-1,1,int(4*self.S)), rd.uniform(-1,1,int(4*self.S))
        evl = 0
        V = self.C*(self.sigma**2 + (1-self.C)*self.mu**2)
        c1 = np.sqrt(self.S * V)
        c2 = - self.d + self.C * self.mu
        a = self.sigma * np.sqrt(self.S * self.C)*(1+(2/np.pi))
        b = self.sigma * np.sqrt(self.S * self.C)*(1-(2/np.pi))
        
        ind = np.sqrt((x)**2+(y)**2) < 1
        if mat =='ran' :
            evl = c1*(x[ind]+ 1.j*y[ind])
        if mat =='sym' :
            evl = np.sqrt(2)*c1*x[ind]
        if mat =='pp' :
            evl = b*x[ind]+ 1.j*a*y[ind]
        if mat =='mix' :
            evl = a*x[ind]+ 1.j*b*y[ind]
            
        evl = np.concatenate((evl[:int(self.S/2)],np.conj(evl[:int(self.S/2)]))) + c2
        return evl
        
        
        
        
