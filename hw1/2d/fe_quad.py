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
        numQuadElements = 0
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            if (int(temp[1]) == 1): numLineElements = numLineElements + 1
            if (int(temp[1]) == 3): numQuadElements = numQuadElements + 1

        lineElements = np.zeros((numLineElements, 2), dtype = int)
        quadElements = np.zeros((numQuadElements, 4), dtype = int)
        lineCoun = 0
        quadCoun = 0
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            if (int(temp[1]) == 1): 
                for j in range(2):
                    lineElements[lineCoun, j]  = int(temp[5 + j])
                lineCoun = lineCoun + 1
            if (int(temp[1]) == 3): 
                for j in range(4):
                    quadElements[quadCoun, j]  = int(temp[5 + j])
                quadCoun = quadCoun + 1

        return mshNodes, lineElements, quadElements 

    def getConnectivity(self):
        mshNodes, lineElements, quadElements = self.read()

        numNodes    = mshNodes.shape[0]

        numEle      = quadElements.shape[0]
        numDOFNodes = mshNodes.shape[0] - np.unique(lineElements).shape[0]

        nodes = np.zeros((numEle, 4, 2))

        for i in range(numEle):
            for j in range(4):
                nodes[i, j, 0] = mshNodes[quadElements[i, j] - 1, 0]
                nodes[i, j, 1] = mshNodes[quadElements[i, j] - 1, 1]
        
        LM = quadElements     #For quad 

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

        return numEle, numNodes, numDOFNodes, nodes, LM 


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
            self.getLocalStiffNess(i, K)
            self.getLocalMass(i, M)
            self.getLocalLoad(i, F)

        return K, F 

    def getLocalLoad(self, elNo, F):
        numNodes = 4 #For quad 
        m = np.zeros((numNodes, numNodes))

        #Crude way to calculate element dimension 
        hx = (np.max(self.nodes[elNo, :, 0]) - np.min(self.nodes[elNo, :, 0]))
        hy = (np.max(self.nodes[elNo, :, 1]) - np.min(self.nodes[elNo, :, 1]))

        #Quad mass matrix 
        m[0, 0] = 4 
        m[0, 1] = 2 
        m[0, 2] = 1 
        m[0, 3] = 2 

        m[1, 0] = 2 
        m[1, 1] = 4 
        m[1, 2] = 2 
        m[1, 3] = 1 

        m[2, 0] = 1 
        m[2, 1] = 2 
        m[2, 2] = 4 
        m[2, 3] = 2 

        m[3, 0] = 2 
        m[3, 1] = 1 
        m[3, 2] = 2 
        m[3, 3] = 4 

        m = m*hx*hy/36.0

        g = np.ones((numNodes)) #Assuming all loads are 1       
        f = np.zeros((numNodes)) 
        for i in range(numNodes):
            for j in range(numNodes):
                f[i] = f[i] + g[j]*m[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]


    def getLocalMass(self, elNo, M):
        numNodes = 4 #For quad 
        m = np.zeros((numNodes, numNodes))

        #Quad mass matrix 
        m[0, 0] = 4 
        m[0, 1] = 2 
        m[0, 2] = 1 
        m[0, 3] = 2 

        m[1, 0] = 2 
        m[1, 1] = 4 
        m[1, 2] = 2 
        m[1, 3] = 1 

        m[2, 0] = 1 
        m[2, 1] = 2 
        m[2, 2] = 4 
        m[2, 3] = 2 

        m[3, 0] = 2 
        m[3, 1] = 1 
        m[3, 2] = 2 
        m[3, 3] = 4 

        m = m/36.0

        for i in range(numNodes):
            ip = LM[elNo, i]
            for j in range(numNodes):
                iq = LM[elNo, j]
                if (ip != -1) and (iq != -1):
                    M[ip, iq] = M[ip, iq] + m[i, j]


    def getLocalStiffNess(self, elNo, K):
        numNodes = 4 #For quad 
        k = np.zeros((numNodes, numNodes)) 

        #Quad stiffness calculation
        alpha = 1.0 #For square quad

        t = np.zeros((numNodes)) 

        t[0] = 1/(3*alpha) + alpha/3.0
        t[1] = alpha/6.0   - 1/(3*alpha) 
        t[2] =-alpha/6.0   - 1/(6*alpha) 
        t[3] = 1/(6*alpha) - alpha/3.0

        for m in range(numNodes):
            k[m, m] = t[0]
            k[0, m] = t[m]

        k[1, 0] = t[1]
        k[1, 2] = t[3]
        k[1, 3] = t[2]

        k[2, 0] = t[2]
        k[2, 1] = t[3]
        k[2, 3] = t[1]

        k[3, 0] = t[3]
        k[3, 1] = t[2]
        k[3, 2] = t[1]

        for m in range(numNodes):
            ip = LM[elNo, m]
            for n in range(numNodes ):
                iq = LM[elNo, n]
                if (ip != -1) and (iq != -1):
                    K[ip, iq] = K[ip, iq] + k[m, n]


if __name__=="__main__":
    fileName                            = 'test_quad.msh'
    msh                                 = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM = msh.getConnectivity()

    kappa = 1.0
    run   = FE(kappa, numEle, numNodes, numDOF, nodes, LM)
    K, F  = run.build()

    U     = np.linalg.solve(K, F)
