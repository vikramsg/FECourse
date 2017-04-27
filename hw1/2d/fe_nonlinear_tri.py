import numpy as np
import re

class readGMSH:

    def __init__(self, fileName):
        self.fineName = fileName

        self.mshNodes = None

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

        self.mshNodes = mshNodes

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

        newDOF = np.zeros((numDOFNodes, 2), dtype = int) #Associating olf DOFs with new DOFs

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

    def buildNonLinear(self, U):
        K = np.zeros((self.numDOF, self.numDOF))

        F = np.zeros((self.numDOF))
        R = np.zeros((self.numDOF))

        G = np.zeros((self.numDOF))

        for i in range(numEle):
            det_D = self.assembleStiffness(i, K)
            self.assembleLoad(i, F, det_D)
            self.assembleG(i, U, G)

        R = F - G

        return R, K

 
    def assembleG(self, elNo, U, G):
        numNodes = 3 #For triangle 

        det_D, k = self.getLocalStiffness(elNo)

        g = np.zeros((numNodes)) #Element g
        for i in range(numNodes):
            for j in range(numNodes):
                jp = LM[elNo, j]
                #Need to multiply with load at the point so U[LM(a, j)]
                #This assumes there is no load at boundary points
                if (jp != -1):
                    g[i] = g[i] + U[jp]*k[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    G[ip] = G[ip] + g[i]


    def assembleLoad(self, elNo, F, det_D):
        numNodes = 3 #For triangle 

        m = self.getLocalMass(elNo, det_D)

        g = np.ones((numNodes)) #Assuming all loads are 1       
        f = np.zeros((numNodes)) 
        for i in range(numNodes):
            for j in range(numNodes):
                f[i] = f[i] + g[j]*m[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]

    def assembleStiffness(self, elNo, K):
        numNodes = 3 #For triangle 

        det_D, k = self.getLocalStiffness(elNo)

        for m in range(numNodes):
            ip = LM[elNo, m]
            for n in range(numNodes ):
                iq = LM[elNo, n]
                if (ip != -1) and (iq != -1):
                    K[ip, iq] = K[ip, iq] + k[m, n]

        return det_D



    def getLocalMass(self, elNo, det_D):
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

        return m


    def getLocalStiffness(self, elNo):
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

        return det_D, k




if __name__=="__main__":
    fileName                                    = 'old/test_tri.msh'
#    fileName                                    = 'proj_test.msh'
    msh                                         = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM, newDOF = msh.getConnectivity()

    tol    = 1e-10 #Tolerance

    kappa  = 1.0
    run    = FE(kappa, numEle, numNodes, numDOF, nodes, LM)

    U      = np.zeros((numDOF)) #Solution vector

    for i in range(10):
        R, K   = run.buildNonLinear(U)
        
        deltaD = np.linalg.solve(K, R)
        
        U      = U + deltaD

        if (np.linalg.norm(deltaD) < tol): break


    for i, l in enumerate(newDOF):
#        print(l[0], U[i])
        print(msh.mshNodes[l[0] - 1][0], msh.mshNodes[l[0] - 1][1], U[i])

