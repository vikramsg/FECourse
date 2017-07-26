import numpy as np
import re
import time

class readGMSH:

    def __init__(self, fileName):
        '''
        Initialize variables
        '''
        self.fileName = fileName

        self.mshNodes = None

        self.totalNumEle   = 0 

        self.nodes    = None
        self.numNodes = 0 
        self.numEle   = 0 
        self.numDOF   = 0 

        self.LM       = None

    def read(self):
        '''
        Read in mesh
        '''
        op = open(fileName, 'r')
        lines = op.readlines()
        op.close()

        for it, i in enumerate(lines):
            if re.search('\$Nodes', i) is not None:
                break
        numNodes = int(lines[it + 1])

        #Read in nodes
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

        self.totalNumEle = numElements

        numPointElements = 0 #Boundary nodes
        numLineElements  = 0 #Boundary lines
        numTriElements   = 0
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            if (int(temp[1]) ==15): numPointElements = numPointElements + 1
            if (int(temp[1]) == 1): numLineElements  = numLineElements  + 1
            if (int(temp[1]) == 2): numTriElements   = numTriElements   + 1

        pointElements = np.zeros((numPointElements, 2), dtype = int)
        lineElements  = np.zeros((numLineElements,  3), dtype = int)
        triElements   = np.zeros((numTriElements,   4), dtype = int)
        pointCoun = 0
        lineCoun  = 0
        triCoun   = 0
        #Store element connectivity for boundary and interior elements
        for coun, i in enumerate(lines[it + 2: it + 2 + numElements]):
            temp     = i.split()
            numTags = int(temp[2])
            physTag = int(temp[3])
            subTag  = 3 + numTags  # Element in line from where element connectivity starts
            if (int(temp[1]) == 15): 
                numNodes = 1
                for j in range(numNodes):
                    pointElements[pointCoun, j]    = int(temp[subTag + j])
                pointElements[pointCoun, numNodes] = physTag 
                pointCoun = pointCoun + 1
            if (int(temp[1]) == 1): 
                numNodes = 2
                for j in range(numNodes):
                    lineElements[lineCoun, j]     = int(temp[subTag + j])
                lineElements[lineCoun, numNodes]  = physTag 
                lineCoun = lineCoun + 1
            if (int(temp[1]) == 2): 
                numNodes = 3
                for j in range(numNodes):
                    triElements[triCoun, j]       = int(temp[subTag + j])
                triElements[triCoun, numNodes]    = physTag 
                triCoun = triCoun + 1

        return mshNodes, pointElements, lineElements, triElements 

    def getConnectivity(self):
        '''
        Create LM matrix
        '''
        mshNodes, pointElements, lineElements, triElements = self.read()

        numNodes    = mshNodes.shape[0]

        numEle      = triElements.shape[0]
        numDOFNodes = mshNodes.shape[0] - np.unique(np.delete(pointElements, -1, 1)).shape[0] # Only point nodes are closed

        #For each element, get the x, y values of nodes
        nodes = np.zeros((numEle, 3, 2))
        IEN   = np.zeros((numEle, 3), dtype = int)
        for i in range(numEle):
            for j in range(3):
                nodes[i, j, 0] = mshNodes[triElements[i, j] - 1, 0]
                nodes[i, j, 1] = mshNodes[triElements[i, j] - 1, 1]

                IEN[i, j]      = triElements[i, j] - 1

        LM      = triElements     #For triangle 

        #Associate boundary nodes with 0
        for cn1, i in enumerate(LM):
            for cn2, j in enumerate(i[:-1]):
                for k in np.unique(np.delete(pointElements, -1, 1)): # Only point nodes are closed
                    if (j == k): LM[cn1, cn2] = 0
        
        arr = LM
        arr = np.delete(arr, -1, 1) # Remove physical tags
        arr = np.unique(arr)        # Get unique node numbers

        boundLM = lineElements    #For line elements to implement boundary conditions 

        newDOF = np.zeros((numDOFNodes, 2), dtype = int) #Associating olf DOFs with new DOFs

        cn = 0
        for i in arr:
            if ( i != 0):
                newDOF[cn, 0] = i 
                newDOF[cn, 1] = cn + 1
                cn            = cn + 1

        #Renumber LM with new DOFs
        for cn1, i in enumerate(LM):
            for cn2, j in enumerate(i[:-1]):
                if (LM[cn1, cn2] != 0):
                    for k in newDOF:
                        if (LM[cn1, cn2] == k[0]):
                            LM[cn1, cn2] = k[1]

        for cn1, i in enumerate(LM):
            for cn2, j in enumerate(i[:-1]):
                LM[cn1, cn2] = LM[cn1, cn2] - 1 #Boundary nodes become -1

        for cn1, i in enumerate(boundLM):
            for cn2, j in enumerate(i[:-1]):
                boundLM[cn1, cn2] = boundLM[cn1, cn2] - 1 #Boundary nodes become -1

        # Array defining connectivity between IEN and boundLM 
        # Number of boundary elements it contains, boundary element number 1, 2, 3
        boundToLM = -1*np.ones((numEle,  4), dtype = int) 

        for cn1, i in enumerate(IEN):
            for cn2, j in enumerate(boundLM):
                isPresent = 0 
                for cn3, k in enumerate(j[:-1]):
                    for cn4, l in enumerate(i):
                        if (k == l):
                            isPresent = isPresent + 1 
                if (isPresent == 2):
                    if (boundToLM[cn1][0] == -1):
                        boundToLM[cn1][0] = 1 
                    else:
                        boundToLM[cn1][0] = boundToLM[cn1][0] + 1
                    for it in range(3):
                        if (boundToLM[cn1][it] == -1):
                            boundToLM[cn1][it] = cn2
                            break
        return numEle, numNodes, numDOFNodes, newDOF, nodes, IEN, LM, boundLM, boundToLM



class FE:
    '''
    Class for building matrices and vectors required for the problem
    '''

    def __init__(self, numEle, numNodes, numDOF, nodes, IEN, LM, boundLM, boundToLM):
        '''
        Initialize variables
        '''
        self.eleNumNodes = 3          #For triangle

        self.kappa       = 1.0 

        self.numEle      = numEle     #Number of elements
        self.numDOF      = numDOF     #Number of open DOF

        self.numNodes    = numNodes   #Total number of nodes in the mesh

        self.nodes       = nodes      #Array containing x, y co-ordinates of nodes in each element
        self.IEN         = IEN        #Array containing global number of nodes in each element 
        self.LM          = LM         #LM array 
        self.boundLM     = boundLM    #boundary LM array 
        self.boundToLM   = boundToLM  #Array defining connectivity between elements and boundary

        self.load_coef   = 0.0        #Default load


    def buildNonLinear(self, U):
        '''
        Build non-linear matrices
        '''
        #Initialize vectors and matrices
        K = np.zeros((self.numDOF, self.numDOF))

        F = np.zeros((self.numDOF))
        R = np.zeros((self.numDOF))

        G = np.zeros((self.numDOF))

        for i in range(self.numEle):
            det_D = self.assembleStiffness(i, U, K) #Assemble tangent stiffness matrix
            self.assembleLoad(i, F, det_D)          #Assemble load vector
#            self.assembleG(i, U, G)                 #Assemble G vector 
#
#        R = F - G

        return R, F, K

    def assembleKappa(self, U):
        '''
        Assemble kappa for each element. For post-processing
        '''
        Kappa = np.zeros(self.numEle)
        for i in range(self.numEle):
            Kappa[i] = self.getKappa(i, U)
        return Kappa 

    def assembleDerivative(self, U):
        '''
        Assemble derivatives for each element. For post-processing
        '''
        U_x = np.zeros(self.numEle)
        U_y = np.zeros(self.numEle)
        for i in range(self.numEle):
            u_x, u_y = self.getDerivative(i, U)
            U_x[i] = u_x
            U_y[i] = u_y

        return U_x, U_y

    def getDerivative(self, elNo, U):
        '''
        Get du/dx and du/dy
        '''
        numNodes = self.eleNumNodes 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        u_x = 0.0; u_y = 0.0;
        for j in range(numNodes):
            jp = self.LM[elNo, j]
            if (jp != -1):
                u_x = u_x + A[j, 1]*U[jp]
                u_y = u_y + A[j, 2]*U[jp]

        return u_x, u_y

    def assembleG(self, elNo, U, G):
        '''
        Assemble G from each element g
        '''
        numNodes = self.eleNumNodes 

        det_D, k = self.getLocalK(elNo, U)

        g = np.zeros((numNodes)) #Element g
        for i in range(numNodes):
            for j in range(numNodes):
                jp = LM[elNo, j]
                #Need to multiply with load at the point so U[LM(a, j)]
                if (jp != -1):
                    g[i] = g[i] + U[jp]*k[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    G[ip] = G[ip] + g[i]

        

    def assembleLoad(self, elNo, F, det_D):
        '''
        Assemble load from each element load
        '''
        numNodes = self.eleNumNodes 

        m = self.getLocalMass(elNo, det_D)

        g = np.ones((numNodes))*self.load_coef #Assuming uniform loading 
        f = np.zeros((numNodes)) 
        for i in range(numNodes):
            for j in range(numNodes):
                f[i] = f[i] + g[j]*m[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]

        q1 =  170 # Outflow load
        q2 = -170 # Inflow load
        # Only do the following for elements that have a Neumann boundary
        if (boundToLM[elNo][0] != -1):
            for i in range(boundToLM[elNo][0]): # For each boundary element
                boundEl   = boundToLM[elNo][1]  # Line element of boundary
                lineNodes = np.zeros((2, 2))
                coun      = 0
                for j in boundLM[boundEl][:-1]:
                    for k in range(numNodes):
                        if (j == IEN[elNo][k]): # Get nodes of boundary element
                            lineNodes[coun, 0] = self.nodes[elNo, k, 0]
                            lineNodes[coun, 1] = self.nodes[elNo, k, 1]
                            coun = coun + 1
                mStar  = self.getLineMass(lineNodes) # Get line element mass matrix

                line_f = np.zeros((numNodes)) 
                coun   = 0
                for j in boundLM[boundEl][:-1]:
                    for k in range(numNodes):
                        if (j == IEN[elNo][k]): # Get load due to Neumann boundary condition
                            if (boundLM[boundEl][-1] == 3): # boundary tag 3 for inflow
                                load = q2
                            elif (boundLM[boundEl][-1] == 4): # boundary tag 4 for outflow 
                                load = q1
                            line_f[k] = mStar[coun, 0]*load + mStar[coun, 1]*load
                            coun = coun + 1

                for i in range(numNodes):
                    ip = LM[elNo, i]
                    if (ip != -1):
                        F[ip] = F[ip] + line_f[i]



    def getLineMass(self, nodes):
        numNodes = nodes.shape[0] 
        assert(numNodes == 2)

        x1 = nodes[0][0]
        y1 = nodes[0][1]

        x2 = nodes[1][0]
        y2 = nodes[1][1]

        h = (x1 - x2)**2 + (y1 - y2)**2
                
        m = np.zeros((2, 2))

        m[0, 0] = 2 
        m[0, 1] = 1 

        m[1, 0] = 1 
        m[1, 1] = 2 

        m = (h/6.0)*m

        return m
 

    def assembleStiffness(self, elNo, U, K):
        '''
        Assemble K from each element k
        '''
        numNodes = self.eleNumNodes 

        det_D, k = self.getLocalStiffness(elNo, U)

        for m in range(numNodes):
            ip = LM[elNo, m]
            for n in range(numNodes ):
                iq = LM[elNo, n]
                if (ip != -1) and (iq != -1):
                    K[ip, iq] = K[ip, iq] + k[m, n]

        return det_D


    def getLocalK(self, elNo, U):
        '''
        Get local linear stiffness
        '''
        numNodes = self.eleNumNodes 
        k = np.zeros((numNodes, numNodes)) 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        kappa = self.getKappa(elNo, U)

        for m in range(numNodes):
            for n in range(numNodes):
                k[m, n] = kappa*0.5*det_D*(A[m, 1]*A[n, 1] + A[m, 2]*A[n, 2])

        return det_D, k



    def getLocalMass(self, elNo, det_D):
        '''
        Get generic mass matrix for element
        '''
        numNodes = self.eleNumNodes 
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


    def getTriShapeFnCoeff(self, elNo):
        '''
        Get the A matrix for triangle linear shape function
        '''
        numNodes = self.eleNumNodes 

        #Triangle stiffness calculation
        D = np.zeros((numNodes, numNodes)) 

        for m in range(numNodes):
            D[m, 0] = 1.0
            D[m, 1] = self.nodes[elNo, m, 0] 
            D[m, 2] = self.nodes[elNo, m, 1] 

        det_D = np.linalg.det(D)

        A = np.zeros((numNodes, numNodes)) 
        A[0, 0] = (self.nodes[elNo, 1, 0]*self.nodes[elNo, 2, 1] - self.nodes[elNo, 1, 1]*self.nodes[elNo, 2, 0])/det_D
        A[0, 1] = (self.nodes[elNo, 1, 1] - self.nodes[elNo, 2, 1])/det_D
        A[0, 2] = (self.nodes[elNo, 2, 0] - self.nodes[elNo, 1, 0])/det_D

        A[1, 0] = (self.nodes[elNo, 2, 0]*self.nodes[elNo, 0, 1] - self.nodes[elNo, 2, 1]*self.nodes[elNo, 0, 0])/det_D
        A[1, 1] = (self.nodes[elNo, 2, 1] - self.nodes[elNo, 0, 1])/det_D
        A[1, 2] = (self.nodes[elNo, 0, 0] - self.nodes[elNo, 2, 0])/det_D

        A[2, 0] = (self.nodes[elNo, 0, 0]*self.nodes[elNo, 1, 1] - self.nodes[elNo, 0, 1]*self.nodes[elNo, 1, 0])/det_D
        A[2, 1] = (self.nodes[elNo, 0, 1] - self.nodes[elNo, 1, 1])/det_D
        A[2, 2] = (self.nodes[elNo, 1, 0] - self.nodes[elNo, 0, 0])/det_D

        return det_D, A


    def getKappa(self, elNo, U):
        kappa = self.kappa

        return kappa


    def getLocalStiffness(self, elNo, U):
        '''
        At each element, calculate the non-linear part of the tangent stiffness
        and add to the linear part
        '''
        numNodes = self.eleNumNodes

        k = np.zeros((numNodes, numNodes)) 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        kappa = self.getKappa(elNo, U)

        u_x, u_y = self.getDerivative(elNo, U)
        
        det_D, k = self.getLocalK(elNo, U) #Get linear part of k

        return det_D, k


    def project(self, fn, U):
        '''
        Project a function on the solution vector
        '''
        numNodes = 3 #For triangle

        for i in range(self.numEle):
            for j in range(numNodes):
                x     = self.nodes[i, j, :]
                jp    = self.LM[i, j]
                if (jp != -1):
                    U[jp] = fn(x)


def fn(x):
    u = np.cos(0.5*x[1])

    return u

class Post_process:

    def __init__(self, rdFile, wrFile, pointSize):
        '''
        Initialize files to read and write
        pointSize is the the number of points in the mesh
        '''
        self.rdFile = rdFile
        self.wrFile = wrFile

        self.pointSize = pointSize

        rd = open(self.rdFile, 'r')

        wr = open(self.wrFile, 'w')
        
        for i in rd:
            wr.write(i)

        wr.close()

        rd.close()

    def write(self, DOF_connectivity, data, variableName, dataType, writeHeader ):
        '''
        Write data to the post processing file
        boundary points will have 0 written
        '''
        wr = open(self.wrFile, 'a')

        size = data.shape[0]

        if (writeHeader == True):
            if (dataType == 'Point'):
                wr.write('\nPOINT_DATA ' + str(self.pointSize) + '\n')
            elif (dataType == 'Cell'):
                wr.write('\nCELL_DATA ' + str(size) + '\n')
                
        wr.write('\nSCALARS '+ variableName +' float 1 \nLOOKUP_TABLE default \n \n')

        counter = 0

        if (dataType == 'Point'):
            for i in range(self.pointSize):
                '''
                DOF_connectivity is sorted in ascending order
                So we go through the connectivity, get the first free DOF put in the value
                and so on. For closed DOF we simply print 0
                '''
                if (DOF_connectivity[counter, 0] - 1 == i ): #i starts from 0
                    index = DOF_connectivity[counter, 1] - 1
                    st = str(data[index] ) + '\n'
                    counter = counter + 1
                else:
                    st = '0 \n'
    
                wr.write(st)

        elif (dataType == 'Cell'):
            for i in range(size):
                st = str(data[i] ) + '\n'
                wr.write(st)

        wr.close()



if __name__=="__main__":
    fileName                                                             = 'ref_proj.msh'
    msh                                                                  = readGMSH(fileName)
    numEle, numNodes, numDOF, newDOF, nodes, IEN, LM, boundLM, boundToLM = msh.getConnectivity()

    tol    = 1e-4 #Tolerance

    run    = FE(numEle, numNodes, numDOF, nodes, IEN, LM, boundLM, boundToLM)

    U      = np.zeros((numDOF)) #Solution vector initialized to 0

    startTime = time.time() # To measure time required for convergence

    R, F, K   = run.buildNonLinear(U)

    U = np.linalg.solve(K, F)

    U_x, U_y = run.assembleDerivative(U)
    resized_U_x  = np.zeros(msh.totalNumEle); resized_U_x[msh.totalNumEle - numEle:]  = U_x
    resized_U_y  = np.zeros(msh.totalNumEle); resized_U_y[msh.totalNumEle - numEle:]  = U_y

    postFile                        = 'ref_proj.vtk'
    outFile                         = 'post_tri.vtk'

    post = Post_process(postFile, outFile, numNodes)
    post.write(newDOF, U, 'u', dataType = 'Point', writeHeader = True)
    post.write(newDOF, resized_U_x,  'u_x',   dataType = 'Cell', writeHeader = True)
    post.write(newDOF, resized_U_y,  'u_y',   dataType = 'Cell', writeHeader = False)

