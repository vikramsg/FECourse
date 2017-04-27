import numpy as np
import re

class readGMSH:

    def __init__(self, fileName):
        self.fineName = fileName

        self.nodes    = None
        self.numNodes = 0 
        self.numEle   = 0 
        self.numDOF   = 0 

        self.LM       = None

    def read(self):
        op = open(fileName, 'r')
        lines = op.readlines()
        op.close()

        for it, i in enumerate(lines):
            if re.search('\$Nodes', i) is not None:
                break
        numNodes = int(lines[it + 1])

        mshNodes = np.zeros((numNodes, 2))
        for coun, i in enumerate(lines[it + 2: it + 2 + numNodes]):
            temp     = i.split()
            mshNodes[coun, 0] = float(temp[1])
            mshNodes[coun, 1] = float(temp[2])

        for it, i in enumerate(lines):
            if re.search('\$Elements', i) is not None:
                break
        numElements = int(lines[it + 1])

        numLineElements = 0
        numTriElements  = 0
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            if (int(temp[1]) == 1): numLineElements = numLineElements + 1
            if (int(temp[1]) == 2): numTriElements  = numTriElements  + 1

        lineElements = np.zeros((numLineElements, 2), dtype = int)
        triElements  = np.zeros((numTriElements,  3), dtype = int)
        lineCoun = 0
        triCoun  = 0
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            if (int(temp[1]) == 1): 
                for j in range(2):
                    lineElements[lineCoun, j]  = int(temp[5 + j])
                lineCoun = lineCoun + 1
            if (int(temp[1]) == 2): 
                for j in range(3):
                    triElements[triCoun, j]    = int(temp[5 + j])
                triCoun = triCoun + 1

        return mshNodes, lineElements, triElements 

    def getConnectivity(self):
        mshNodes, lineElements, triElements = self.read()

        numNodes    = mshNodes.shape[0]

        numEle      = triElements.shape[0]
        numDOFNodes = mshNodes.shape[0] - np.unique(lineElements).shape[0]

        nodes = np.zeros((numEle, 4, 2))

        for i in range(numEle):
            for j in range(3):
                nodes[i, j, 0] = mshNodes[triElements[i, j] - 1, 0]
                nodes[i, j, 1] = mshNodes[triElements[i, j] - 1, 1]
        
        LM = triElements     #For triangle 

        for cn1, i in enumerate(LM):
            for cn2, j in enumerate(i):
                for k in np.unique(lineElements):
                    if (j == k): LM[cn1, cn2] = 0

        arr = np.unique(LM)

        newDOF = np.zeros((numDOFNodes, 2)) #Associating olf DOFs with new DOFs

        cn = 0
        for i in arr:
            if ( i != 0):
                newDOF[cn, 0] = i 
                newDOF[cn, 1] = cn + 1
                cn            = cn + 1

        #Renumber LM with new DOFs
        for cn1, i in enumerate(LM):
            for cn2, j in enumerate(i):
                if (LM[cn1, cn2] != 0):
                    for k in newDOF:
                        if (LM[cn1, cn2] == k[0]):
                            LM[cn1, cn2] = k[1]

        LM = LM - 1

        return numEle, numNodes, numDOFNodes, nodes, LM, newDOF 



class FE:

    def __init__(self, kappa, numEle, numNodes, numDOF, node, LM):
        self.kappa    = kappa

        self.numEle   = numEle
        self.numDOF   = numDOF

        self.numNodes = numNodes

        self.nodes    = nodes
        self.LM       = LM 

    def build(self):
        K = np.zeros((self.numDOF, self.numDOF))
        M = np.zeros((self.numDOF, self.numDOF))

        F = np.zeros((self.numDOF))

        for i in range(numEle):
            det_D = self.getLocalStiffNess(i, K)
#            self.getLocalMass(i, M, det_D)
            self.getLocalLoad(i, F, det_D)

        return K, F

    def getLocalLoad(self, elNo, F, det_D):
        numNodes = 3 #For triangle 
        m = np.zeros((numNodes, numNodes))

        #Triangle mass matrix 
        m[0, 0] = 2 
        m[0, 1] = 1 
        m[0, 2] = 1 

        m[1, 0] = 1 
        m[1, 1] = 2 
        m[1, 2] = 1 

        m[2, 0] = 1 
        m[2, 1] = 1 
        m[2, 2] = 2 

        m = m*(det_D/24.0)

        g = np.ones((numNodes)) #Assuming all loads are 1       
        f = np.zeros((numNodes)) 
        for i in range(numNodes):
            for j in range(numNodes):
                f[i] = f[i] + g[j]*m[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]



    def getLocalMass(self, elNo, M, det_D):
        numNodes = 3 #For triangle 
        m = np.zeros((numNodes, numNodes))

        #Triangle mass matrix 
        m[0, 0] = 2 
        m[0, 1] = 1 
        m[0, 2] = 1 

        m[1, 0] = 1 
        m[1, 1] = 2 
        m[1, 2] = 1 

        m[2, 0] = 1 
        m[2, 1] = 1 
        m[2, 2] = 2 

        m = m*(det_D/24.0)

        for i in range(numNodes):
            ip = LM[elNo, i]
            for j in range(numNodes):
                iq = LM[elNo, j]
                if (ip != -1) and (iq != -1):
                    M[ip, iq] = M[ip, iq] + m[i, j]



    def getLocalStiffNess(self, elNo, K):
        numNodes = 3 #For triangle
        k = np.zeros((numNodes, numNodes)) 

        #Triangle stiffness calculation
        D = np.zeros((numNodes, numNodes)) 

        for m in range(numNodes):
            D[m, 0] = 1.0
            D[m, 1] = self.nodes[elNo, m, 0] 
            D[m, 2] = self.nodes[elNo, m, 1] 

        det_D = np.linalg.det(D)

        A = np.zeros((numNodes, 2)) 

        A[0, 0] = (self.nodes[elNo, 1, 1] - self.nodes[elNo, 2, 1])/det_D
        A[0, 1] = (self.nodes[elNo, 2, 0] - self.nodes[elNo, 1, 0])/det_D

        A[1, 0] = (self.nodes[elNo, 2, 1] - self.nodes[elNo, 0, 1])/det_D
        A[1, 1] = (self.nodes[elNo, 0, 0] - self.nodes[elNo, 2, 0])/det_D

        A[2, 0] = (self.nodes[elNo, 0, 1] - self.nodes[elNo, 1, 1])/det_D
        A[2, 1] = (self.nodes[elNo, 1, 0] - self.nodes[elNo, 0, 0])/det_D

        for m in range(numNodes):
            for n in range(numNodes):
                k[m, n] = self.kappa*0.5*det_D*(A[m, 0]*A[n, 0] + A[m, 1]*A[n, 1])

        for m in range(numNodes):
            ip = LM[elNo, m]
            for n in range(numNodes ):
                iq = LM[elNo, n]
                if (ip != -1) and (iq != -1):
                    K[ip, iq] = K[ip, iq] + k[m, n]

        return det_D


if __name__=="__main__":
    fileName                                    = 'test_tri.msh'
    msh                                         = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM, newDOF = msh.getConnectivity()

    kappa = 1.0
    run   = FE(kappa, numEle, numNodes, numDOF, nodes, LM)
    K, F  = run.build()

    U     = np.linalg.solve(K, F)

    for i, l in enumerate(newDOF):
        print(l[0], U[i])

