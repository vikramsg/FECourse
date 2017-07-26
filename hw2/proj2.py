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
    '''
    Class for building matrices and vectors required for the problem
    '''

    def __init__(self, kappa, rho_c, numEle, numNodes, numDOF, nodes, LM):
        '''
        Initialize variables
        '''
        self.eleNumNodes = 3        #For triangle

        alpha            = None
        dt               = None

        self.kappa       = kappa    #k_0
        self.rho_c       = rho_c    #rho*c

        self.numEle      = numEle   #Number of elements
        self.numDOF      = numDOF   #Number of open DOF

        self.numNodes    = numNodes #Total number of nodes in the mesh

        self.nodes       = nodes    #Array containing x, y co-ordinates of nodes in each element
        self.LM          = LM       #LM array 

        self.load_case   = None     #Default load

        self.pTime       = None     #Time at the present time step

    def problem_param(self, alpha, dt, load_case):
        self.alpha          = alpha
        self.dt             = dt
        self.load_case      = load_case


    def buildMStar(self):
        '''
        Build Mstar 
        '''
        tol   = 10e-8
        if (np.abs(alpha - 0) > tol):
            Mstar = np.zeros((self.numDOF, self.numDOF))
        else:
            Mstar = np.zeros((self.numDOF))

        for i in range(self.numEle):
            det_D, k = self.getLocalK(i, U)      #Assemble tangent stiffness matrix
            self.assembleMStar(i, det_D, U, Mstar)

        return Mstar


    def buildGTFE(self, U, pTime):
        '''
        Build Mstar and F 
        '''
        self.pTime          = pTime
        
        F     = np.zeros((self.numDOF))

        alpha = self.alpha

        for i in range(self.numEle):
            det_D, k = self.getLocalK(i, U)      #Assemble tangent stiffness matrix
            self.assembleLoad(i, det_D, U, F)    #Assemble load vector

        return F 


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



    def assembleMStar(self, elNo, det_D, U, M):
        '''
        Assemble M from each element g
        '''
        numNodes = self.eleNumNodes 
    
        alpha    = self.alpha
        dt       = self.dt

        rho_c    = self.rho_c

        if (np.abs(alpha - 0) > tol):
            m        = self.getLocalMass(elNo, det_D)
        else:
            m        = self.getLumpedMass(elNo, det_D)

        det_D, k = self.getLocalK(elNo, U)

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (np.abs(alpha - 0) > tol):
                for j in range(numNodes ):
                    iq = LM[elNo, j]
                    if (ip != -1) and (iq != -1):
                        M[ip, iq] = M[ip, iq] + (rho_c*m[i, j] + alpha*dt*k[i, j])
            else:
                if (ip != -1):
                    M[ip] = M[ip] + rho_c*m[i] 



    def assembleLoad(self, elNo, det_D, U, F):
        '''
        Assemble load from each element load
        '''
        numNodes = self.eleNumNodes 

        alpha    = self.alpha
        t        = self.pTime
        dt       = self.dt

        rho_c    = self.rho_c

        tol   = 10e-8

        '''
        If alpha == 0, use lumped mass matrix 
        '''
        if (np.abs(alpha - 0) > tol):
            m        = self.getLocalMass(elNo, det_D)
        else:
            m        = self.getLumpedMass(elNo, det_D)

        det_D, k = self.getLocalK(elNo, U)

        '''
        Create element load vector
        '''
        f = np.zeros(numNodes)
        g = np.zeros(numNodes)
        for i in range(numNodes):
            g[i] = self.getGTLoad(t, dt)
        if (np.abs(alpha - 0) > tol):
            f = np.dot(m, g) # f = f_b*m_ab 
        else:
            for i in range(numNodes):
                f[i] = g[i]*m[i]

        for i in range(numNodes):
            for j in range(numNodes):
                jp = LM[elNo, j]
                if (jp != -1):
                    if (np.abs(alpha - 0) > tol):
                        f[i] = f[i] + (rho_c*m[i, j] - (1 - alpha)*dt*k[i, j]) * U[jp]
                    else:
                        if (i == j):
                            m_mass = m[i]
                        else:
                            m_mass = 0 
                        f[i] = f[i] + (rho_c*m_mass - (1 - alpha)*dt*k[i, j]) * U[jp]

        for i in range(numNodes):
            ip = LM[elNo, i]
            if (ip != -1):
                    F[ip] = F[ip] + f[i]


    def getGTLoad(self, t, dt):
        '''
        dt*(alpha*F_{n + 1} + (1 - alpha)*F_n)
        '''
        alpha = self.alpha

        f_np1 = self.getUnsteadyLoad(t + dt)
        f_n   = self.getUnsteadyLoad(t)

        return dt*(alpha*f_np1 + (1 - alpha)*f_n)

    def getUnsteadyLoad(self, t):
        load_case = self.load_case 
        if self.load_case == 1: 
            f_0 = 2.6*(10**4)
            T   = 300.0
            f   = np.max([0, f_0*np.sin(2*np.pi*t/T)])
        else:
            f = 0

        return f


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


    def getLumpedMass(self, elNo, det_D):
        '''
        Get generic mass matrix for element
        '''
        numNodes = self.eleNumNodes 
        m = np.zeros((numNodes))

        #Triangle Lumped mass matrix 
        m[0] = 4 
        m[1] = 4 
        m[2] = 4 

        m = m*(det_D/24.0)

        return m

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


    def project(self, fn, U, init_case):
        '''
        Project a function on the solution vector
        '''
        numNodes = 3 #For triangle

        for i in range(self.numEle):
            for j in range(numNodes):
                x     = self.nodes[i, j, :]
                jp    = self.LM[i, j]
                if (jp != -1):
                    U[jp] = fn(x, init_case)


def fn(x, init_case):
    if init_case == 0:
        u = 100.0 
    else:
        u = 0.0

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
    fileName                                    = 'ref_proj.msh'
    msh                                         = readGMSH(fileName)
    numEle, numNodes, numDOF, nodes, LM, newDOF = msh.getConnectivity()


    tol    = 1e-8 #Tolerance

    kappa  = 250.0 #kappa_0
    rho_c  = 1.76*10**(6)
    run    = FE(kappa, rho_c, numEle, numNodes, numDOF, nodes, LM)

    load_case = 1
    '''
    0: No load
    1: Unsteady load f = max{0, f_0.sin(2.pi.t/T)}
    '''

    alpha =  1.00
    dt    = 15.00 
    run.problem_param(alpha, dt, load_case) #Set problem parameters

    U      = np.zeros((numDOF))    #Solution vector initialized to 0

    init_case = 1
    '''
    0: U_0 = 100
    1: U_0 = 0
    '''
    run.project(fn, U, init_case)

    startTime = time.time() # To measure time required for convergence

    pTime = 0    #Physical time
    U0    = U

    base_name =  "data_at_b_gf"# Write data at point B and a point near GF
    load_str  =  "_load_case_" + str(load_case) 
    init_str  =  "_init_case_" + str(init_case) 
    alph_str  =  "_alpha_" + str(alpha) 
    u_b_file  = base_name + load_str + init_str + alph_str + '.dat' 
    wt_b_file = open(u_b_file, 'w')
       
    Mstar = run.buildMStar() # Mstar is constant for a given alpha and dt 
    if (np.abs(alpha - 0) > tol):
        Minv = np.linalg.inv(Mstar) # Needs to be inverted only once

    for i in range(3000):
        F = run.buildGTFE(U0, pTime)

        if (np.abs(alpha - 0) > tol):
            U1 = np.dot(Minv, F)
        else:
            U1 = F/Mstar

        U0 = U1

        maxU  = np.max(U0)

        pTime  = pTime + dt

        print("iteration ", i, ", maxU ", maxU, " time ", pTime)

        st_b = str(pTime) + '\t' + str(U1[1]) + "\t" + str(U1[505]) + '\n' #Value of u at B
        wt_b_file.write(st_b)

#        if (i == 5) or (i == 60):
#            postFile  = 'ref_proj.vtk'
#            time_str  =  "_time_" + str(i*dt) 
#            outFile   = 'post_tri'+ load_str + init_str + alph_str + time_str + '.vtk' 
#            
#            post = Post_process(postFile, outFile, numNodes)
#            post.write(newDOF, U1, 'u', dataType = 'Point', writeHeader = True)

    stopTime = time.time()

    wt_b_file.close()

