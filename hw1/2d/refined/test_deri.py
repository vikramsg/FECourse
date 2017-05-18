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
        #Store element connectivity for boundary and interior elements
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
        '''
        Create LM matrix
        '''
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
        
        LM = triElements     #For triangle 

        #Associate boundary nodes with 0
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

    def __init__(self, kappa, case, numEle, numNodes, numDOF, nodes, LM):
        self.eleNumNodes = 3 #For triangle

        self.kappa       = kappa

        self.numEle      = numEle
        self.numDOF      = numDOF

        self.numNodes    = numNodes

        self.nodes       = nodes
        self.LM          = LM 

        self.case        = case 

        self.load_coef   = 1.0


    def buildNonLinear(self, U, load_coef):
        '''
        Build non-linear matrices
        '''
        self.load_coef = load_coef

        K = np.zeros((self.numDOF, self.numDOF))

        F = np.zeros((self.numDOF))
        R = np.zeros((self.numDOF))

        G = np.zeros((self.numDOF))

        for i in range(self.numEle):
            det_D = self.assembleStiffness(i, U, K)
            self.assembleLoad(i, F, det_D) #Assuming uniform load of f = 1
            self.assembleG(i, U, G)

        R = F - G

        return R, F, K

    def assembleKappa(self, U):
        Kappa = np.zeros(self.numEle)
        for i in range(self.numEle):
            Kappa[i] = self.getKappa(i, U)
        return Kappa 

    def assembleDerivative(self, U):
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
        Assemble G for each element
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
        Assemble load for each element
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


    def assembleStiffness(self, elNo, U, K):
        '''
        Assemble K for each element
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
            kappa    = self.kappa*(A*(l1_gradu**4) + B*(l1_gradu**2) + 1)

        return kappa


    def getLocalStiffness(self, elNo, U):
        numNodes = self.eleNumNodes

        k = np.zeros((numNodes, numNodes)) 

        det_D, A = self.getTriShapeFnCoeff(elNo)

        kappa = self.getKappa(elNo, U)

        u_x, u_y = self.getDerivative(elNo, U)
        
        det_D, k = self.getLocalK(elNo, U) #Get linear part of k

        for m in range(numNodes):
            for n in range(numNodes):
                if (self.case == 2):
                    dk_dux  = -2*self.kappa*u_x/(1 + u_x**2 + u_y**2)**2 #dk/du_x
                    dk_duy  = -2*self.kappa*u_y/(1 + u_x**2 + u_y**2)**2
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

    def __init__(self, rdFile, wrFile, pointSize):
        '''
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
    fileName                                    = 'ref_proj.msh'
    msh                                         = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM, newDOF = msh.getConnectivity()

    tol    = 1e-8 #Tolerance

    kappa        = 1.0 #kappa_0
    kappa_case   = 3
    '''
    1. kappa = kappa_0
    2. kappa = kappa_0/(1 + u_x**2 + u_y**2)
    3. kappa = kappa_0*(A*|grad u|**4 + B*|grad u|**2 + 1)
    '''

    run    = FE(kappa, kappa_case, numEle, numNodes, numDOF, nodes, LM)

    U      = np.zeros((numDOF)) #Solution vector initialized to 0

    newton_case = 1
    '''
    1. Newton 
    2. Modified Newton 
    '''

    inc_case    = 1
    '''
    1. No incremental loading
    2. Incremental loading with 10 steps
    '''
    if (inc_case == 2):
        num_increments = 10
    elif (inc_case == 1):
        num_increments = 1 

    load = 0.1940 #f_0

    if (inc_case == 1):
        load_coef   = load #For incremental loading
    elif (inc_case == 2):
        load_coef   = load/num_increments #For incremental loading

    u_vs_load = np.zeros((num_increments, 2)) #For plotting U at a particular node vs load

    iterR = [] #List containing norm of R for each iteration

    startTime = time.time() # To measure time required for convergence
    for j in range(1, num_increments + 1):
        load_coef = j*(load/num_increments)
            
        print("load increment iteration ", j, ", load_coef ", load_coef, " norm(U) ", np.linalg.norm(U))

        I_rc = 3 #Convergence flag

        R, F, K   = run.buildNonLinear(U, load_coef)
        for i in range(500):
            deltaD = np.linalg.solve(K, R)
            
            U      = U + deltaD
            
            normU  = np.linalg.norm(U)

            normD  = np.linalg.norm(deltaD)

            if (normD/normU < tol) :
                '''
                1. increment should be within tolerance
                '''
                I_rc = 2
    
            iterR.append(np.linalg.norm(R))

            if (newton_case == 1):
                R, F, K       = run.buildNonLinear(U, load_coef)
            elif (newton_case == 2):
                R, F, K_temp  = run.buildNonLinear(U, load_coef) #Use original K

            normR  = np.linalg.norm(R)
            normF  = np.linalg.norm(F)

            print("iteration ", i, ", normD ", normD, " normR ", normR, " normF ", normF )

            if (normR/normF < tol) and (I_rc == 2):
                '''
                1. increment should be within tolerance
                '''
                break
            elif (normR/normF < tol) and (I_rc != 2):
                I_rc = 1
            elif (normR/normF > tol) and (I_rc != 2):
                I_rc = 3

        stopTime = time.time()

        u_vs_load[j - 1, 0] = load_coef 
        u_vs_load[j - 1, 1] = U[LM[900, 1]] #Store u at node 1 of element 900

    wrFile = open('node_u_vs_load_' + str(kappa_case) + '_' + str(newton_case) + '_' + str(inc_case) +'.dat', 'w')
    st = 'kappa_case ' + str(kappa_case) + ' Newton_method ' + str(newton_case) + ' loading ' + str(inc_case) + '\n' 
    wrFile.write(st)

    for i in u_vs_load:
        '''
        Write u at particular node vs f_0 with each incremental load to file
        '''
        st = str(i[0]) + '\t' + str(i[1]) + '\n'
        wrFile.write(st)
    wrFile.close()

    wrFile = open('normR_vs_iteration_' + str(kappa_case) + '_' + str(newton_case) + '_' + str(inc_case) +'.dat', 'w')
    st = 'kappa_case ' + str(kappa_case) + ' Newton_method ' + str(newton_case) + ' loading ' + str(inc_case) + '\n' 
    st = st  + 'stopTime - startTime ' + str(stopTime - startTime) + '\n'
    wrFile.write(st)

    for coun, i in enumerate(iterR):
        '''
        Write norm of R vs iterations 
        '''
        st = str(coun) + '\t' + str(i) + '\n'
        wrFile.write(st)
    wrFile.close()

    postFile                        = 'ref_proj.vtk'
    outFile                         = 'post_tri_' +str(kappa_case) + '_' + str(newton_case) + '_' + str(inc_case) + '.vtk'

    U_x, U_y = run.assembleDerivative(U)
    Kappa    = run.assembleKappa(U)

    resized_U_x  = np.zeros(msh.totalNumEle); resized_U_x[msh.totalNumEle - numEle:]  = U_x
    resized_U_y  = np.zeros(msh.totalNumEle); resized_U_y[msh.totalNumEle - numEle:]  = U_y
    resizedKappa = np.zeros(msh.totalNumEle); resizedKappa[msh.totalNumEle - numEle:] = Kappa 

    post = Post_process(postFile, outFile, numNodes)
    post.write(newDOF, U, 'u', dataType = 'Point', writeHeader = True)
    post.write(newDOF, resized_U_x,  'u_x',   dataType = 'Cell', writeHeader = True)
    post.write(newDOF, resized_U_y,  'u_y',   dataType = 'Cell', writeHeader = False)
    post.write(newDOF, resizedKappa, 'kappa', dataType = 'Cell', writeHeader = False)
