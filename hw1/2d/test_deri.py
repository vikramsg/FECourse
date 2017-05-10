import numpy as np
import re


class Mesh:

    def __init__(self):
        self.nodes    = None
        self.numNodes = 0 
        self.numEle   = 0 

        self.LM       = None

    def getConnectivity(self):
        self.numEle = 4
        self.getMesh(self.numEle)

        self.LM = np.zeros((self.numEle, 3)) #For triangle

        self.LM[0, 0] = 0
        self.LM[0, 1] = 1
        self.LM[0, 2] = 2

        self.LM[1, 0] = 4
        self.LM[1, 1] = 2
        self.LM[1, 2] = 1

        self.LM[2, 0] = 1
        self.LM[2, 1] = 3
        self.LM[2, 2] = 4

        self.LM[3, 0] = 2
        self.LM[3, 1] = 4
        self.LM[3, 2] = 5

        newDOF = np.zeros((self.numNodes, 2), dtype = int) #Associating olf DOFs with new DOFs
        for i, j in enumerate(newDOF):
            newDOF[i, 0] = i
            newDOF[i, 1] = i

        return self.numEle, self.numNodes, self.numNodes, self.nodes, self.LM, newDOF



    def getMesh(self, numEle):
        self.numEle = numEle

        nodes = np.zeros((numEle, 3, 2))

        nodes[0, 0, 0] = -1.5
        nodes[0, 0, 1] =  0.0
        nodes[0, 1, 0] = -0.75
        nodes[0, 1, 1] =  0.0
        nodes[0, 2, 0] = -0.75 
        nodes[0, 2, 1] =  0.5

        nodes[1, 0, 0] = -0.0
        nodes[1, 0, 1] =  0.5
        nodes[1, 1, 0] = -0.75
        nodes[1, 1, 1] =  0.5
        nodes[1, 2, 0] = -0.75 
        nodes[1, 2, 1] =  0.0

        nodes[2, 0, 0] = -0.75
        nodes[2, 0, 1] =  0.0
        nodes[2, 1, 0] = -0.00
        nodes[2, 1, 1] =  0.0
        nodes[2, 2, 0] = -0.00 
        nodes[2, 2, 1] =  0.5

        nodes[3, 0, 0] = -0.75
        nodes[3, 0, 1] =  0.5
        nodes[3, 1, 0] = -0.00
        nodes[3, 1, 1] =  0.5
        nodes[3, 2, 0] = -0.00 
        nodes[3, 2, 1] =  1.0

        numNodes = 6

        self.numNodes = numNodes
        self.nodes    = nodes

        return nodes, numNodes


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

        numLineElements = 0 #Boundary nodes
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

        #For each element, get the x, y values of nodes
        nodes = np.zeros((numEle, 3, 2))
        for i in range(numEle):
            for j in range(3):
                nodes[i, j, 0] = mshNodes[triElements[i, j] - 1, 0]
                nodes[i, j, 1] = mshNodes[triElements[i, j] - 1, 1]
        
        LM = triElements     #For triangle #It is not necessarily numbered from 0

#        print(LM)

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

        LM = LM - 1 #Boundary nodes become -1

        return numEle, numNodes, numDOFNodes, nodes, LM, newDOF 



class FE:

    def __init__(self, kappa, case, numEle, numNodes, numDOF, node, LM):
        self.eleNumNodes = 3 #For triangle

        self.kappa    = kappa

        self.numEle   = numEle
        self.numDOF   = numDOF

        self.numNodes = numNodes

        self.nodes    = nodes
        self.LM       = LM 

        self.case     = case 

    def buildNonLinear(self, U):
        K = np.zeros((self.numDOF, self.numDOF))

        F = np.zeros((self.numDOF))
        R = np.zeros((self.numDOF))

        G = np.zeros((self.numDOF))

        for i in range(self.numEle):
            det_D = self.assembleStiffness(i, U, K)
            self.assembleLoad(i, F, det_D) #Assuming uniform load of f = 1
            self.assembleG(i, U, G)

        R = F - G

        return R, K

 
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
        numNodes = self.eleNumNodes 

        m = self.getLocalMass(elNo, det_D)

        g = 0.30*np.ones((numNodes)) #Assuming all loads are 1       
        f = np.zeros((numNodes)) 
        for i in range(numNodes):
            for j in range(numNodes):
                f[i] = f[i] + g[j]*m[i, j]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]


    def assembleStiffness(self, elNo, U, K):
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
        numNodes = self.eleNumNodes 
        k = np.zeros((numNodes, numNodes)) 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        kappa = self.getKappa(elNo, U)

        for m in range(numNodes):
            for n in range(numNodes):
                k[m, n] = kappa*0.5*det_D*(A[m, 1]*A[n, 1] + A[m, 2]*A[n, 2])

        return det_D, k



    def getLocalMass(self, elNo, det_D):
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
        u_x, u_y = self.getDerivative(elNo, U)

        '''
        1. kappa = kappa_0
        2. kappa = kappa_0/(1 + u_x**2 + u_y**2)
        3. kappa = kappa_0*(A*|grad u|**4 + B*|grad u|**2 + 1)
        '''

        if (self.case == 1):
            kappa = self.kappa
        elif (self.case == 2):
            kappa = self.kappa/(1 + u_x**2 + u_y**2)
        elif (self.case == 3):
            A =  0.18
            B = -0.82
            l1_gradu = np.sqrt(u_x**2 + u_y**2)
            kappa    = self.kappa*(A*l1_gradu**4 + B*l1_gradu**2 + 1)

        return kappa


    def getLocalStiffness(self, elNo, U):
        numNodes = self.eleNumNodes

        k = np.zeros((numNodes, numNodes)) 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        kappa = self.getKappa(elNo, U)

        u_x, u_y = self.getDerivative(elNo, U)
        
        det_D, k = self.getLocalK(elNo, U)

        for m in range(numNodes):
            for n in range(numNodes):
                if (self.case == 2):
                    dk_dux  = -2*self.kappa*u_x/(1 + u_x**2 + u_y**2) #dk/du_x
                    dk_duy  = -2*self.kappa*u_y/(1 + u_x**2 + u_y**2)
                    dux_ddb = A[n, 1] #du_x/d d_b = A_b2
                    duy_ddb = A[n, 2] #du_x/d d_b = A_b3
                    for p in range(numNodes):
                        jp      = self.LM[elNo, p]
                        if (jp != -1): #U = 0 for boundary
                            d_c     = U[jp]
                            k_ac    = 0.5*det_D*(A[m, 1]*A[p, 1] + A[m, 2]*A[p, 2])
                            k[m, n] = k[m, n] + d_c*k_ac*(dk_dux*dux_ddb + dk_duy*duy_ddb)
                            
                elif (self.case == 3):
                    A_k =  0.18
                    B_k = -0.82

                    dk_dux  = self.kappa*(2*B_k*u_x + A_k*(4*u_x**3 + 4*u_x*u_y**2)) #dk/du_x
                    dk_duy  = self.kappa*(2*B_k*u_y + A_k*(4*u_y**3 + 4*u_y*u_x**2)) #dk/du_y
                    dux_ddb = A[n, 1] #du_x/d d_b = A_b2
                    duy_ddb = A[n, 2] #du_x/d d_b = A_b3
                    for p in range(numNodes):
                        jp      = self.LM[elNo, p]
                        if (jp != -1): #U = 0 for boundary
                            d_c     = U[jp]
                            k_ac    = 0.5*det_D*(A[m, 1]*A[p, 1] + A[m, 2]*A[p, 2])
                            k[m, n] = k[m, n] + d_c*k_ac*(dk_dux*dux_ddb + dk_duy*duy_ddb)

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

    def __init__(self, fileName):
        self.fileName = fileName

    def write(self, outFile, DOF_connectivity, data, pointSize):
        '''
        Write data to the post processing file
        pointSize is the the number of points in the mesh
        boundary points will have 0 written
        '''
        rd = open(self.fileName, 'r')
        wr = open(outFile, 'w')

        for i in rd:
            wr.write(i)

        size = data.shape[0]

        wr.write('\nPOINT_DATA ' + str(pointSize) + '\n')

        wr.write('\nSCALARS data float 1 \nLOOKUP_TABLE default \n \n')

        counter = 0
        for i in range(pointSize):
            '''
            DOF_connectivity is sorted in ascending order
            So we go through the connectivity, get the first free DOF put in the value
            and so on. For closed DOF we simply print 0
            '''
            if (DOF_connectivity[counter, 0] - 1 == i ): #i starts from 0
                index = DOF_connectivity[counter, 1] - 1
                st = str(U[index] ) + '\n'
                counter = counter + 1
            else:
                st = '0 \n'

            wr.write(st)


        rd.close()
        wr.close()



if __name__=="__main__":
#    fileName                                    = 'old/test_tri.msh'
    fileName                                    = 'proj_test.msh'
    msh                                         = readGMSH(fileName)
#    msh                                         = Mesh()
    numEle, numNodes, numDOF, nodes, LM, newDOF = msh.getConnectivity()

    tol    = 1e-10 #Tolerance

    kappa  = 1.0
    case   = 3
    '''
    1. kappa = kappa_0
    2. kappa = kappa_0/(1 + u_x**2 + u_y**2)
    3. kappa = kappa_0*(A*|grad u|**4 + B*|grad u|**2 + 1)
    '''

    run    = FE(kappa, case, numEle, numNodes, numDOF, nodes, LM)

    U      = np.zeros((numDOF)) #Solution vector
    
    for i in range(50):
        R, K   = run.buildNonLinear(U)
        
        deltaD = np.linalg.solve(K, R)
        
        U      = U + deltaD

        normD  = np.linalg.norm(deltaD)

        if (normD < tol):# FIXME Not all conditions have been included
            break

        print("iteration ", i, ", norm ", normD)

    for i, l in enumerate(newDOF):
#        print(l[0], l[1])
        print(l[0], U[i])
#        print(msh.mshNodes[l[0] - 1][0], msh.mshNodes[l[0] - 1][1], U[i])

    
#    postFile                                    = 'old/test_tri.vtk'
    postFile                                    = 'proj_test.vtk'
    outFile                                     = 'post_tri.vtk'

    post = Post_process(postFile)
    post.write(outFile, newDOF, U, numNodes)
