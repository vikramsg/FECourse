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


class testMesh:

    def __init__(self):
        self.nodes    = None
        self.numNodes = 0 
        self.numEle   = 0 
        self.numDOF   = 0 

        self.LM       = None

    def getConnectivity(self):
        self.LM = np.zeros((self.numEle, 4), dtype=int) #For quad 

        self.LM[0, 0] = -1 
        self.LM[0, 1] = -1
        self.LM[0, 2] =  1
        self.LM[0, 3] =  0

        self.LM[1, 0] = -1
        self.LM[1, 1] = -1
        self.LM[1, 2] = -1
        self.LM[1, 3] =  1

        self.LM[2, 0] =  0
        self.LM[2, 1] =  1
        self.LM[2, 2] =  3
        self.LM[2, 3] =  2

        self.LM[3, 0] =  1
        self.LM[3, 1] = -1
        self.LM[3, 2] = -1
        self.LM[3, 3] =  3

        self.LM[4, 0] =  2
        self.LM[4, 1] =  3
        self.LM[4, 2] =  5
        self.LM[4, 3] =  4

        self.LM[5, 0] =  3
        self.LM[5, 1] = -1
        self.LM[5, 2] =  6
        self.LM[5, 3] =  5

        self.LM[6, 0] = -1
        self.LM[6, 1] = -1
        self.LM[6, 2] =  7
        self.LM[6, 3] =  6

        self.LM[7, 0] = -1
        self.LM[7, 1] = -1
        self.LM[7, 2] = -1
        self.LM[7, 3] =  7


        return self.LM



    def getMesh(self, numEle):
        self.numEle = numEle

        nodes = np.zeros((numEle, 4, 2))

        nodes[0, 0, 0] = -0.0
        nodes[0, 0, 1] =  0.0
        nodes[0, 1, 0] =  1.00
        nodes[0, 1, 1] =  0.0
        nodes[0, 2, 0] =  1.00 
        nodes[0, 2, 1] =  1.0
        nodes[0, 3, 0] =  0.00 
        nodes[0, 3, 1] =  1.0

        nodes[1, 0, 0] =  1.0
        nodes[1, 0, 1] =  0.0
        nodes[1, 1, 0] =  2.00
        nodes[1, 1, 1] =  0.0
        nodes[1, 2, 0] =  2.00 
        nodes[1, 2, 1] =  1.0
        nodes[1, 3, 0] =  1.00 
        nodes[1, 3, 1] =  1.0

        nodes[2, 0, 0] =  0.0
        nodes[2, 0, 1] =  1.0
        nodes[2, 1, 0] =  1.00
        nodes[2, 1, 1] =  1.0
        nodes[2, 2, 0] =  1.00 
        nodes[2, 2, 1] =  2.0
        nodes[2, 3, 0] =  0.00 
        nodes[2, 3, 1] =  2.0

        nodes[3, 0, 0] =  1.0
        nodes[3, 0, 1] =  1.0
        nodes[3, 1, 0] =  2.00
        nodes[3, 1, 1] =  1.0
        nodes[3, 2, 0] =  2.00 
        nodes[3, 2, 1] =  2.0
        nodes[3, 3, 0] =  1.00 
        nodes[3, 3, 1] =  2.0

        nodes[4, 0, 0] =  0.0
        nodes[4, 0, 1] =  2.0
        nodes[4, 1, 0] =  1.00
        nodes[4, 1, 1] =  2.0
        nodes[4, 2, 0] =  1.00 
        nodes[4, 2, 1] =  3.0
        nodes[4, 3, 0] =  0.00 
        nodes[4, 3, 1] =  3.0

        nodes[5, 0, 0] =  1.0
        nodes[5, 0, 1] =  2.0
        nodes[5, 1, 0] =  2.00
        nodes[5, 1, 1] =  2.0
        nodes[5, 2, 0] =  2.00 
        nodes[5, 2, 1] =  3.0
        nodes[5, 3, 0] =  1.00 
        nodes[5, 3, 1] =  3.0

        nodes[6, 0, 0] =  2.0
        nodes[6, 0, 1] =  2.0
        nodes[6, 1, 0] =  3.00
        nodes[6, 1, 1] =  2.0
        nodes[6, 2, 0] =  3.00 
        nodes[6, 2, 1] =  3.0
        nodes[6, 3, 0] =  2.00 
        nodes[6, 3, 1] =  3.0

        nodes[7, 0, 0] =  3.0
        nodes[7, 0, 1] =  2.0
        nodes[7, 1, 0] =  4.00
        nodes[7, 1, 1] =  2.0
        nodes[7, 2, 0] =  4.00 
        nodes[7, 2, 1] =  3.0
        nodes[7, 3, 0] =  3.00 
        nodes[7, 3, 1] =  3.0

        numNodes = 16

        numDOF   = 8

        self.numNodes = numNodes
        self.nodes    = nodes

        self.numDOF   = numDOF  

        return nodes, numNodes, numDOF

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

        for i in range(numEle):
            self.getLocalStiffNess(i, K)
            self.getLocalMass(i, M)

        print(M)

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


    def getLocalLoad(self):
        l = np.zeros((2)) #For 1D

if __name__=="__main__":
    fileName                            = 'test.msh'
    msh                                 = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM = msh.getConnectivity()

#    numEle = 8
#
#    msh                     = testMesh()
#    nodes, numNodes, numDOF = msh.getMesh(numEle)
#    LM                      = msh.getConnectivity()
#
    kappa = 1.0
    run = FE(kappa, numEle, numNodes, numDOF, nodes, LM)
    run.build()
