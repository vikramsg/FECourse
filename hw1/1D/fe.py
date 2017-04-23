import numpy as np

class Mesh:

    def getMesh(self, numEle, x_l, x_r):
        dx = (x_r - x_l)/numEle

        numNodes = numEle + 1

        nodes = np.zeros((numEle, 2))

        for i in range(0, numEle):
            nodes[i, 0] = i*dx
            nodes[i, 1] = i*dx + dx

        return nodes, numNodes

class FE:

    def __init__(self, k_coeff, numEle, numNodes, node):
        self.k_coeff = k_coeff

        self.numEle   = numEle
        self.numNodes = numNodes
        self.nodes    = nodes

    def build(self):
        K = np.zeros((self.numNodes, self.numNodes))

        for i in range(numEle):
            self.getLocalStiffNess(i, K)

        print(K)

    def getLocalStiffNess(self, i, K):
        numNodes = 2
        k = np.zeros((numNodes, numNodes)) #For 1D

        h = self.nodes[i, 1] - self.nodes[i, 0]

        k[0, 0] = (k_coeff/h)
        k[0, 1] =-(k_coeff/h)
        k[1, 0] =-(k_coeff/h)
        k[1, 1] = (k_coeff/h)

        for m in range(i, i + 2):
            for n in range(i, i + 2):
                K[m, n] = K[m, n] + k[m - i, n - i]

    def getLocalLoad(self):
        l = np.zeros((2)) #For 1D

if __name__=="__main__":
    numEle = 5
    x_l    = 0.0
    x_r    = 1.0

    msh             = Mesh()
    nodes, numNodes = msh.getMesh(numEle, x_l, x_r)

    k_coeff = 0.1
    run = FE(k_coeff, numEle, numNodes, nodes)
    run.build()
