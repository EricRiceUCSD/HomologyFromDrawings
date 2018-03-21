import numpy as np
from matplotlib import pyplot as plt
from math import *
from itertools import permutations
from PIL import Image as im

def getValsByCoord(png):
    """Given a black and white image png given as a string
    of the filename, gathers the value of the pixel at
    position (x, y) in png as 255 (white) or 0 (black)
    and returns a 2D array of these values indexed as
    (x, y) = [y][x] (following row, column convention)."""

    canvas = im.open(png)
    canvas = canvas.convert(mode="L") # converts to black and white
    vals = np.array(canvas.getdata()) # 1D array with values of pixels left to right, row by row 
    vals = np.split(vals, canvas.getbbox()[3]) # converts into 2D array
                                               # unfortunately, (x, y) = vals[y][x]

    return vals


def getReducedVertexSet(vals, b):
    """Given a rectangular 2D array of pixel values vals
    and an integer b, decomposes vals into b x b
    blocks, ignoring any remainder.  Returns a list S
    of block coordinates of blocks which contain at
    least one black pixel."""

    # Checks if vals is rectagular
    check = True
    for y in range(0, len(vals)):
        if len(vals[0]) != len(vals[y]):
            check = False
    if check == False:
        print("In function getVertexSet, input array was not rectangular.")
        return None

    # Checks if b is valid
    if b > len(vals) or b > len(vals[0]):
        print("Block size " + str(b) + " is too large.")
    
    h = len(vals) // b      # height in block coords
    w = len(vals[0]) // b   # width in block coords

    S = []
    for y in range(0, h):
        for x in range(0, w):
            for s in range(0, b):
                if ((x, y),) in S:
                    break
                else:
                    for t in range(0, b):
                        if vals[y*b + s][x*b + t] == 0:
                            S += [((x, y),)]  # tuple of tuple for consistency
                            break             # with higher dimensions

    return S


def dist(p, q):
    """Given two points p, q in R^n, returns
    the Euclidean distance between them."""

    if len(p) != len(q):
        print("In dist(p, q), p and q are not the same length.")
        return None
    
    normsq = 0
    for i in range(0, len(p)):
        normsq += (p[i] - q[i])**2

    return sqrt(normsq)


class CechComplex():
    """Given a vertex set S of points in R^n (where each point
    is contained in a tuple of length 1) and a radius r,
    returns the Cech complex K of S where:
    - K.complex is a list of lists
    - K.complex[0] = S is the collection of 0-simplices
    - K.complex[1] = all pairs of points in S within Euclidean
             distance 2r from each other; these are
             the 1-simplices
    - K.complex[2] = all triples of points in S within Euclidean
             distance 2r from each other; these are
             the 2-simplices, etc."""

    complex = None
    radius = None
    
    def __init__(self, S, r):
        K = [S, [], []]
        p = 0
        for j in range(0, len(K[p])):
            for k in range(j + 1, len(S)):  # Avoids double checking by not repeating dist(p, q) via dist(q, p)
                if dist(S[j][0], S[k][0]) <= 2*r:
                    K[1] += [(S[j][0], S[k][0])]

        p = 1          
        while len(K[p]) != 0:
            for j in range(0, len(K[p])):
                for k in range(0, len(S)):

                    # Checks if the vertex S[k][0] is within 2*r of all of the
                    # vertices of the simplex K[p][j]
                    close = True
                    for b in range(0, len(K[p][j])):
                        if dist(S[k][0], K[p][j][b]) > 2*r:
                            close = False
                    if S[k][0] not in K[p][j] and close == True:
                        
                        # Checking for present tuples which are permutations 
                        check = True
                        perms = list(permutations((S[k][0],) + K[p][j]))
                        for t in range(0, len(perms)):
                            if perms[t] in K[p + 1]:
                                check = False

                        if check == True:
                            K[p + 1] += [(S[k][0],) + K[p][j]]
    
            p += 1
            K += [[]]
            
        self.complex = K
        self.radius = r

    def bdryMatrix(self):
        """Generates a list of the boundary matrices of the complex
        in Smith normal form, where the index of the list matches
        the dimension."""

        n = [len(self.complex[0])] # ranks of chain groups, n_p in text
        bdry = [None]
        for p in range(1, len(self.complex)):
            n += [len(self.complex[p])]
            bdry += [np.zeros((n[p - 1], n[p]), dtype=int)]
            for i in range(0, n[p - 1]):
                for j in range(0, n[p]):
                    face = True
                    for k in range(0, len(self.complex[p - 1][i])):
                        if self.complex[p - 1][i][k] not in self.complex[p][j]:
                            face = False
                    if face:
                        bdry[p][i][j] = 1

        # Reduction to Smith normal form
        k = None
        ell = None
        for p in range(1, len(self.complex)):
            for x in range(0, min(n[p - 1], n[p])):
                for i in range(x, n[p - 1]):
                    for j in range(x, n[p]):
                        if bdry[p][i][j] == 1:
                            k = i
                            ell = j
                        
                if k != None and ell != None:
                    bdry[p][[x,k]] = bdry[p][[k,x]]          # swaps rows k and x
                    bdry[p][:,[x,ell]] = bdry[p][:,[ell,x]]  # swaps columns ell and x

                    for i in range(x + 1, n[p - 1]):
                        if bdry[p][i][x] == 1:
                            for m in range(0, n[p]):
                                bdry[p][i][m] = (bdry[p][i][m] + bdry[p][x][m]) % 2
                    for j in range(x + 1, n[p]):
                        if bdry[p][x][j] == 1:
                            for m in range(0, n[p - 1]):
                                bdry[p][m][j] = (bdry[p][m][j] + bdry[p][m][x]) % 2
            
        return bdry

    def betti(self):
        """Calculates the Betti numbers of the complex, where
        the output is a list whose index matches the dimension."""

        bdry = self.bdryMatrix()
        n = [] # ranks of chain groups, n_p in text
        for p in range(0, len(self.complex)):
            n += [len(self.complex[p])]

        # Building b_{p - 1} and z_p (the ranks of the boundary and
        # cycle groups) for all p that we have a boundary matrix.
        b = []
        z = [n[0]]
        for p in range(1, len(self.complex)):
            index = -1
            while index + 1 < min(n[p - 1], n[p]):
                if bdry[p][index + 1][index + 1] == 1:
                    index += 1
                else:
                    break

            b += [index + 1]   # Indices are offset from ranks by 1
            z += [n[p] - b[p - 1]]
        b += [0]   # Highest dimension boundary group is trivial

        beta = []
        for p in range(0, len(self.complex)):
            beta += [z[p] - b[p]]
        
        return beta

    def graph(self):
        S = self.complex[0]
        x = []
        y = []
        for k in range(0, len(S)):
            x += [S[k][0][0]]
            y += [-S[k][0][1]]

        plt.scatter(x, y)
        plt.xlim(min(x) - 1, max(x) + 1)
        plt.ylim(min(y) - 1, max(y) + 1)
        plt.show()


def splitByComponent(K):
    """Given a simplicial complex K, returns a list
    of complexes C of the connected components of K."""

    beta = K.betti()
    r = K.radius
    C = []

    # Identifies the connected components of K
    if beta[0] > 1:

        # Initializing a list of sets of vertices, which we will
        # manipulate so that two vertices belong in the same
        # set if and only if they are part of the same component.
        V = []
        for i in range(0, len(K.complex[0])):
            V += [set(K.complex[0][i])]

        for i in range(0, len(K.complex[0])):
            for j in range(0, len(K.complex[0])):
                find_i = None
                find_j = None
                for k in range(0, len(V)):
                    if K.complex[0][i][0] in V[k]:
                        find_i = k
                    if K.complex[0][j][0] in V[k]:
                        find_j = k
                if ((K.complex[0][i][0], K.complex[0][j][0]) in K.complex[1]\
                  or (K.complex[0][j][0], K.complex[0][i][0]) in K.complex[1])\
                  and find_i != find_j:
                    V[find_i] = V[find_i].union(V[find_j])
                    del V[find_j]

        # Checking if separation was successful
        verts = []
        for i in range(0, len(V)):
            for j in range(0, len(V[i])):
                verts += [(list(V[i])[j],)]
        if len(V) != beta[0] or set(verts) != set(K.complex[0]):
            print("splitByComponent(K) was not successful.")
            return None
                    
        # Creates a new complex for each component (given as
        # a set of vertices)
        S = []
        for i in range(0, len(V)):
            S += [[]]
            V[i] = list(V[i])
            for j in range(0, len(V[i])):
                S[i] += [(V[i][j],)]   # Putting entries of V[i] in proper form
            C += [CechComplex(S[i], r)]

        # Generates a list order of indices so that
        # C[order[0]] is the leftmost complex (based on
        # the smallest x-coordinate of a vertex) and
        # C[order[i]] is further left than C[order[j]]
        # if and only if i < j.   
        mins = []   # list to keep track of smallest x coordinates
        for i in range(0, len(C)):
            leftest = C[i].complex[0][0][0][0]
            for k in range(0, len(C[i].complex[0])):
                if C[i].complex[0][k][0][0] < leftest:
                    leftest = C[i].complex[0][k][0][0]
            mins += [leftest]

        order = []
        unordered = C
        for i in range(0, len(C)):
            a = min(mins)
            order += [mins.index(a)]
            mins[mins.index(a)] = max(mins) + 1
            C[i] = unordered[order[i]]

    else:
        C = [K]

    return C

# --- Tests --- #
    
##vals = getValsByCoord("CanvasPNG.png")
##
##b = 20
##r = 0.8
##S = getReducedVertexSet(vals, b)
##S = [((0, 1),), ((0, 0),), ((1, 0),)]                                        # triangle
##S = [((1, 0),), ((sqrt(2)/2, sqrt(2)/2),), ((0, 1),),\                       # empty circle
##     ((-sqrt(2)/2, sqrt(2)/2),), ((-1, 0),), ((-sqrt(2)/2, -sqrt(2)/2),),\
##     ((0, -1),), ((sqrt(2)/2, -sqrt(2)/2),)]
##S = [((0, 0, 0),), ((1, 0, 0),), ((0, 1, 0),), ((0, 0, 1),)]                 # tetrahedron
##
##K = CechComplex(S, r)
##C = splitByComponent(K)
