from __future__ import print_function
import os
import math
import sys
import optparse    # for option sorting
import uuid
import random
import numpy as np
import numpy.linalg as linalg
from sklearn.manifold import MDS
import time
from scipy import optimize
from scipy.optimize import minimize
from decimal import *

class dStruct():
        """
        the class with distance
        """
        def __init__(self, aa, rr, output, ss, no, L, ctype):
                #constructor, initializing values: target rr file to contact map, ss, aa, and other values
                self.aa = aa
                self.rr = rr
                self.ss = ss
                self.no = no
                self.L = L
                self.ctype = ctype
                self.output = output                
                _, aa_file = os.path.split(aa)
                aa_file = aa_file.split('.')[0]
                o, _ = os.path.split(output)
                _, rrfile = os.path.split(rr)
                _, ssfile = os.path.split(ss)
                
                self.aa_file = aa_file 

                


                print ('Target             : ' + aa_file)

                self.target = rr
                fr = open(self.target, 'r')
                lines = fr.readlines()
                faa = open(aa, 'r')
                falines = faa.readlines()
                aafasta = falines[1].strip()
                Sequence = aafasta ## new
                n = len(aafasta) #taking the number of residues
                pairList = []
                self.CM = [[0 for x in range(n)] for y in range(n)]
                self.gdfCM = [[0 for x in range(n)] for y in range(n)]
                self.pair_u = {} # new
                self.pair_l = {} # new
                maxt = 0 #new

                l = 0
                for line in lines[1:]:
                        line = line.strip().split()
                        
                        if(l > n * float(L)):

                                break
                        self.CM[int(line[0]) - 1][int(line[1]) - 1] = self.CM[int(line[1]) - 1][int(line[0]) - 1] = 1
                        l += 1
                        pairList.append((int(line[0]) - 1, int(line[1]) - 1))

                        self.pair_u[(int(line[0]) - 1, int(line[1]) - 1)] = float(line[3]) # new 
                        self.pair_u[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[3]) # new 

                        self.pair_l[(int(line[0]) - 1, int(line[1]) - 1)] = float(line[2]) # new
                        self.pair_l[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[2]) # new

                        if (maxt < float(line[3])):
                                maxt = float(line[3])

                pairlist = pairList.reverse()

                for i in range(len(self.CM)):
                        for j in range(i + 1):
                                if(i - j < 3):
                                        self.CM[i][j] = self.CM[j][i] = 1


                #line1 = lines[1] #new
                #line1 = line1.strip().split() #new
                self.t = maxt #new particularly for dmpfold
                #self.t = float(line1[3]) #considering t as thresh hold, given in each line at column 4

                print ('Length             : ' + str(n))
                print ('Threshold          : ' + str(self.t) +' (Angstrom)')
                if (L == '10000'):
                        print ('Cutoff xL          : All')
                else:
                        print ('Cutoff xL          : ' + L)
                print ('Contact type       : ' + ctype)
                self.T = math.floor(self.t)  #intially taking floor of thresh hold divided by distance between 2 consecutive residues need some calculations based on CM
                fr.close()
                fletter = {'A':'ALA', 'B':'ASX' , 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN',
                        'P':'PRO', 'Q':'GLN', 'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR', 'Z':'GLX', 'U':'   ', 'X':'   ', '*':'   ', '-':'   '} #from https://zhanglab.ccmb.med.umich.edu/FASTA/
                self.res_name = []
                self.fasta = aa
                fs = open(self.fasta, 'r')
                flines = fs.readlines()
                #print ('Sequence           : ' + Sequence))
                for line in flines[1:]:
                        line = line.strip()
                        for i in line:

                                self.res_name.append(fletter[i])
                fs.close()
                
                self.ss = ss

                fss = open(self.ss, 'r')
                self.secondary = fss.read()
                self.secondary = self.secondary.split() #secondary structure is now in self.secondary[0]
                fss.close()
                #print ('Secondary structure: ' + str(self.secondary[0]))
                self.helices = []
                self.sheets = []
                i = 1
                while(i < len(self.secondary[0])):

                        if self.secondary[0][i] == 'H':

                                j = i
                                while (self.secondary[0][j] == 'H' and j < len(self.secondary[0]) - 1):
                                        j = j + 1
                                if (j - i >= 3):
                                        self.helices.append((i, j))
                                        i = j
                        i = i + 1
                i = 1
                while(i < len(self.secondary[0])):

                        if self.secondary[0][i] == 'E':

                                j = i
                                while (self.secondary[0][j] == 'E' and j < len(self.secondary[0]) - 1):
                                        j = j + 1
                                if (j - i >= 3):
                                        self.sheets.append((i, j-1))
                                        i = j
                        i = i + 1

                #print ('helix set: '+ str(self.helices))
                #print ('strand set: ' + str(self.sheets))

                for i in range(len(self.CM)):
                        for j in range(len(self.CM)):
                                self.gdfCM[i][j] = self.CM[i][j]


        def __str__(self):
                return self.get_name()
        def get_name(self):
                return self.__class__.__name__
        def set_target(self):
                self.target = str(target)

        #this method returns partially random distance between two residues based on empirical knowledge about protein
        def distance(self, t, i, j):
                m = max(i , j)
                if i == j:
                        return 0
                elif abs(i - j) == 1:
                        return 3.8
                elif abs(i - j) == 2: # P-SEA d2


                        if(self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H' and self.secondary[0][m-2] == 'H'):
                                return 5.5 + random.uniform(-0.5 , 0.5) # P-SEA
                        elif(self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E' and self.secondary[0][m-2] == 'E'):
                                return 6.7 + random.uniform(-0.6 , 0.6) # P-SEA

                        else:


                                return 6 + random.uniform(-1.5 , 1.5) #COMAR

                elif abs(i - j) == 3: # P-SEA d3

                        if(self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H' and self.secondary[0][m-2] == 'H' and self.secondary[0][m-3] == 'H'):
                                return 5.3 + random.uniform(-0.5 , 0.5) # P-SEA
                                self.gdfCM[i][j] = self.gdfCM[j][i] = 1
                                self.CM[i][j] = self.CM[j][i] = 1
                        elif(self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E' and self.secondary[0][m-2] == 'E' and self.secondary[0][m-3] == 'E'):
                                return 9.9 + random.uniform(-0.9 , 0.9) # P-SEA
                                self.CM[i][j] = self.CM[j][i] = 0
                        else:
                                self.gdfCM[i][j] = self.gdfCM[j][i] = 1
                                return 6.45 + random.uniform(-1, 1) 

                                self.CM[i][j] = self.CM[j][i] = 1

                elif abs(i - j) > 3:

                        if(abs(i - j) == 4 and self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H' and self.secondary[0][m-2] == 'H' and self.secondary[0][m-3] == 'H' and self.secondary[0][m-4] == 'H'):#d4
                                return 6.4 + random.uniform(-0.6 , 0.6) # P-SEA
                                self.gdfCM[i][j] = self.gdfCM[j][i] = 1
                                self.CM[i][j] = self.CM[j][i] = 1
                        elif(abs(i - j) == 4 and self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E' and self.secondary[0][m-2] == 'E' and self.secondary[0][m-3] == 'E' and self.secondary[0][m-4] == 'E'):
                                return 12.4 + random.uniform(-0.1 , 0.1) # P-SEA
                                self.CM[i][j] = self.CM[j][i] = 1 #new
                        else:

                                
                                return random.gauss(6.64, 1.36) 
                                """
                                if (self.CM[i][j] == 1): #new
                                        mid = (self.pair_l[(i,j)] + self.pair_u[(i,j)])/2
                                        std = mid - self.pair_l[(i,j)]
                                        return random.gauss(mid, std)
                                        self.gdfCM[i][j] = 1
                                elif abs(i - j) <= 5:
                                        return random.gauss(6.64, 1.36) #new
                                        self.gdfCM[i][j] = 1 #new
                                """

        #following two methods are for alpha helix optimization using lbfgs
        def lbfgsOpt(self, pos, D):
                D2 = self.coordToDist(pos) #the D2 is initial distrance matrix from pos
                newpos = []
                fu = 0

                for i in range(1, len(pos)):
                        newpoint = pos[i]
                        cons = ({'type':'eq', 'fun': lambda x: (math.sqrt((x[0] - pos[i-1][0]) ** 2 + (x[1] - pos[i-1][0]) ** 2 + (x[2] - pos[i-1][0]) ** 2) - 3.8) ** 2})
                        res = minimize(self.lbfgsOptFunction, (newpoint[0], newpoint[1], newpoint[2]), args = (pos, i, D), method = 'L-BFGS-B', tol=1e-6, options={'maxiter':3000}) #, constraints=cons)
                        pos[i][0] = res.x[0]
                        pos[i][1] = res.x[1]
                        pos[i][2] = res.x[2]
                        fu += float(self.lbfgsOptFunction(newpoint, pos, i, D))

                return pos, fu


        def lbfgsOptFunction(self, newpoint, pos, i, D):
                f2 = 0
                sd2 = 0.5
                for s,e in self.helices:

                        if (i >= s and i < e): #that means i th residue is helics

                                for k in range(s, e):
                                        f2 += ((math.sqrt((newpoint[0] - pos[k][0]) ** 2 + (newpoint[1] - pos[k][1]) ** 2 + (newpoint[2] - pos[k][2]) ** 2) - D[k][i]) ** 2 / (sd2 ** 2))
                return f2

        #this method initally generates a sparse matrix; then by using shortest path algorithm, it obtains proximity map which is refined by secondary structure angular space to distance conversion procedure
        def guess_dist(self):
                n = len(self.CM)

                for i in range(len(self.gdfCM)):
                        for j in range(i+2, len(self.gdfCM)):
                                if (self.gdfCM[i][j] == 1):
                                        k = j
                                        l = j
                                        m = j
                                        while (k != i+2):
                                                if (self.gdfCM[i][k-1] == 0):
                                                        self.gdfCM[i][k-1] = self.gdfCM[i][k] + 1
                                                        k -= 1
                                                elif (self.gdfCM[i][k-1] > self.gdfCM[i][k]):
                                                        self.gdfCM[i][k-1] = self.gdfCM[i][k] + 1
                                                        k -= 1
                                                else:
                                                        break
                                        while (l != len(self.gdfCM)-1):

                                                if (self.gdfCM[i][l+1] == 0):
                                                        self.gdfCM[i][l+1] = self.gdfCM[i][l] + 1
                                                        l += 1
                                                elif (self.gdfCM[i][l+1] > self.gdfCM[i][l]):
                                                        self.gdfCM[i][l+1] = self.gdfCM[i][l] + 1
                                                        l += 1
                                                else:
                                                        break


                D = [[0 for x in range(n)] for y in range(n)]
                D2 = [[0 for x in range(n)] for y in range(n)]
                for i in range (n):
                        for j in range(i+1, n):
                                if self.CM[i][j] == 1:
                                        D[i][j] = self.distance(float(self.t), i, j)
                                        D[j][i] = D[i][j]
                                else:
                                        if (abs(i - j) < 6):
                                                D[i][j] = self.distance(float(self.t), i, j)
                                                D[j][i] = D[i][j]


                                        else: 
                                                D[i][j] = self.gdfCM[i][j] * 5.72 #gdfuz3D #new

                                D[j][i] = D[i][j]
                                D2[i][j] = D[i][j]
                                D2[j][i] = D2[i][j]

                D = self.ChunksImproveDist(D)
                path1 = self.floyd_shortest_path(D)
                return path1

        #shortest path algorithm for upper bound on predicted distance
        def floyd_shortest_path(self,G):
                n = len(G)
                path = G
                for k in range(n):
                        for i in range(n):
                                for j in range(i):
                                        path[i][j] = path[j][i] = min(path[i][j] , path[i][k]+ path[k][j])
                return path

        #cost checking for chirality (lower cost indicates correct chirality)   
        def cAlphaBetaCost(self, pos):
                chiralAlpha = [[0, 0, 0], [0, 0, 0], [0, 0, 0]] #the matrix |xi - xi-1 yi - yi-1 zi - zi-1| for alpha helix
                chiralBeta = [[0, 0, 0], [0, 0, 0], [0, 0, 0]] #the matrix |xi - xi-1 yi - yi-1 zi - zi-1| for beta sheet
                chiDist = [] 
                alphaCost = 0
                betaCost = 0

                idealAlpha = 0.778
                idealBeta = -0.043
                weightAlpha = 1.0
                weightBeta = 0.1
                chiDist.append(math.sqrt((pos[0][0] - pos[1][0]) * (pos[0][0] - pos[1][0]) + (pos[0][1] - pos[1][1]) * (pos[0][1] - pos[1][1]) + (pos[0][2] - pos[1][2]) * (pos[0][2] - pos[1][2]))) #d(1,0) chiDist[0]
                chiDist.append(math.sqrt((pos[1][0] - pos[2][0]) * (pos[1][0] - pos[2][0]) + (pos[1][1] - pos[2][1]) * (pos[1][1] - pos[2][1]) + (pos[1][2] - pos[2][2]) * (pos[1][2] - pos[2][2]))) #d(2,1) chiDist[1]
                for i in range(3, len(pos)):

                        chiDist.append(math.sqrt((pos[i-1][0] - pos[i][0]) * (pos[i-1][0] - pos[i][0]) + (pos[i-1][1] - pos[i][1]) * (pos[i-1][1] - pos[i][1]) + (pos[i-1][2] - pos[i][2]) * (pos[i-1][2] - pos[i][2]))) #d(3,2)-...chiDist[2]-... 

                        if(self.secondary[0][i] == 'E' and self.secondary[0][i-1] == 'E' and self.secondary[0][i-2] == 'E' and self.secondary[0][i-3] == 'E'): #beta sheet

                                chiralBeta[0][0] = pos[i-2][0] - pos[i-3][0] # Xi+1 - Xi
                                chiralBeta[0][1] = pos[i-2][1] - pos[i-3][1] # Yi+1 - Yi
                                chiralBeta[0][2] = pos[i-2][2] - pos[i-3][2] # Zi+1 -Zi

                                chiralBeta[1][0] = pos[i-1][0] - pos[i-2][0] # Xi+2 - Xi+1
                                chiralBeta[1][1] = pos[i-1][1] - pos[i-2][1] # Yi+2 - Yi+1
                                chiralBeta[1][2] = pos[i-1][2] - pos[i-2][2] # Zi+2 - Zi+1

                                chiralBeta[2][0] = pos[i][0] - pos[i-1][0] # Xi+3 - Xi+2
                                chiralBeta[2][1] = pos[i][1] - pos[i-1][1] # Yi+3 - Yi+2
                                chiralBeta[2][2] = pos[i][2] - pos[i-1][2] # Zi+3 - Zi+2
                                chiB = (((np.linalg.det(chiralBeta)))  * (1.0 / (chiDist[i-1] * chiDist[i-2] * chiDist[i-3]))) #(1/d(i)d(i-1)d(i-2))

                                betaCost = betaCost + (2 / (1 + math.exp(-10 * (chiB - idealBeta) * (chiB - idealBeta)))) - 1  #O Lund et al equation 3, s = 10



                        if(self.secondary[0][i] == 'H' and self.secondary[0][i-1] == 'H' and self.secondary[0][i-2] == 'H' and self.secondary[0][i-3] == 'H'):


                                chiralAlpha[0][0] = pos[i-2][0] - pos[i-3][0] # Xi+1 - Xi
                                chiralAlpha[0][1] = pos[i-2][1] - pos[i-3][1] # Yi+1 - Yi
                                chiralAlpha[0][2] = pos[i-2][2] - pos[i-3][2] # Zi+1 -Zi

                                chiralAlpha[1][0] = pos[i-1][0] - pos[i-2][0] # Xi+2 - Xi+1
                                chiralAlpha[1][1] = pos[i-1][1] - pos[i-2][1] # Yi+2 - Yi+1
                                chiralAlpha[1][2] = pos[i-1][2] - pos[i-2][2] # Zi+2 - Zi+1

                                chiralAlpha[2][0] = pos[i][0] - pos[i-1][0] # Xi+3 - Xi+2
                                chiralAlpha[2][1] = pos[i][1] - pos[i-1][1] # Yi+3 - Yi+2
                                chiralAlpha[2][2] = pos[i][2] - pos[i-1][2] # Zi+3 - Zi+2
                                chiA = (((np.linalg.det(chiralAlpha)))  * (1.0 / (chiDist[i-1] * chiDist[i-2] * chiDist[i-3])))

                                alphaCost = alphaCost + (2 / (1 + math.exp(-10 * (chiA - idealAlpha) * (chiA - idealAlpha)))) - 1  # O Lund et al equation 3, s = 10

                return (betaCost + alphaCost) 

        #this method reverses the 3D coordinates to its y axis   
        def mirror(self, pos):

                posy = [[0 for x in range(3)] for y in range (len(pos))]
                for i in range(len(pos)):

                        posy[i][0] = pos[i][0]
                        posy[i][1] = -pos[i][1]
                        posy[i][2] = pos[i][2]


                return posy
        #uses python sklearn library for 3D multidimensional scaling, get initial (20 or more) 3D protein residue coordinates, selects best 5, refines using contact information, ss information and empirical knowledge about protein and generates 5 models as input to the modeller
        def mds3(self):

                D = self.guess_dist()
                D = self.ChunksImproveDist(D)
                mds = MDS(n_components=3, max_iter=300000, eps=1e-9, random_state = np.random.RandomState(seed = int(self.no)),
                   dissimilarity="precomputed", n_jobs=1)

                pos = mds.fit(D).embedding_
                dist = prevDist = 0
                posd = []
                #print ('Generating 20 models ...')
                print ('\nStarting Stage 1 of 3 ...')
                for i in range(20):
                        random.seed(i+int(self.no))
                        D = self.guess_dist()
                        D = self.ChunksImproveDist(D)
                        seed1 = np.random.RandomState(seed = int(i+int(self.no)))
                        mds = MDS(n_components=3, max_iter=300000, eps=1e-9, random_state = seed1,
                                dissimilarity="precomputed", n_jobs=1)
                        posd.append(mds.fit(D).embedding_)
                pos = self.correctCM(pos)

                posc = self.correctCM(pos)

                print('Done! \n')
                print('Starting Stage 2 of 3 ...')
                eps = [x * 0.0085 for x in range(0, 100)]
                eps1 = eps
                #print ('Initial refining 20 models ...')
                for i in range(len(posd)):
                        #as the proximity map is generated using shortest path which is basically upper bound of distances, so, making models compact
                        xc, yc, zc = self.calculateCentroid(posd[i])
                        Rg = 2.2 * (len(posd[i]) ** 0.38)

                        centroid = [xc, yc, zc]
                        cRg = self.calculateRg(posd[i], [xc, yc, zc])

                        posd[i] = self.perturbCentroid(posd[i], [xc, yc, zc], .01)
                        cRg = self.calculateRg(posd[i], [xc, yc, zc])

                        for k in reversed(eps1):
                                if (self.consistant(posd[i]) and k < 0.45):
                                        break

                                posd[i] = self.perturb(posd[i], k)
                                posd[i] = self.perturbCorrect(posd[i], k)
                                posd[i] = self.correctCM(posd[i])


                for i in range(len(posd)):

                        posd[i] = self.correctConSC(posd[i])
                        posd[i] = self.correctClose(posd[i])
                        posd[i] = self.perturbCorrect(posd[i], 0.01)
                        posd[i] = self.correctCM(posd[i])
                        for k in range(25):
                                posd[i] = self.perturbCorrect(posd[i], 0.1)
                                posd[i] = self.correctCM(posd[i])
                for i in range(len(posd)):
                        prevf = float('inf')
                        minf = float('inf')
                        for j in range(5): #for i in range(2 * len(self.helices)):
                                posd[i], f = self.lbfgsOpt(posd[i], D)

                                posd[i] = self.correctCM(posd[i])
                for i in range(len(posd)):


                        c, nc, l = self.cmError(posd[i])

                        cm = c + nc
                        if (cm > 30):
                                for j in range(100):
                                        if (cm < 20):
                                                break
                                        posd[i] = self.perturbCorrect(posd[i], 0.1)
                                        posd[i] = self.correctCM(posd[i])
                                        c, nc, l = self.cmError(posd[i])
                                        cm = c + nc

                #changed to filter best 5 out of 20 based on cm score
                best = []
                for b in range(len(posd)):
                        c, nc, l = self.cmError(posd[b])
                        cmScore = c + nc + len(l)
                        best.append({'pos' : posd[b], 'score' :cmScore})
                sbest = sorted(best, key = lambda i : i['score'])
                #print ('Generated initial top 5 models! ')
                print('Done! \n')
                print('Starting Stage 3 of 3 ...')
                for i in range(5):
                        c = self.cAlphaBetaCost(sbest[i]['pos'])
                        cM = self.cAlphaBetaCost(self.mirror(sbest[i]['pos']))

                        _, _, l = self.cmError(sbest[i]['pos'])
                        #print ('contact map mismatch pairs for model ' + str(i) + ': ' + str(l))

                        if (c < cM):
                                self.toPDB(sbest[i]['pos'], 'bestmodel' + str(i+1) + '.pdb')

                        else:
                                self.toPDB(self.mirror(sbest[i]['pos']), 'bestmodel' + str(i+1) + '.pdb')

        #calculates the center of protein
        def calculateCentroid(self, pos):
                xc = 0
                yc = 0
                zc = 0
                for i in range(len(pos)):
                        #print (pos[i][0])
                        xc += float(pos[i][0])
                        yc += float(pos[i][1])
                        zc += float(pos[i][2])
                xc = xc / len(pos)
                yc = yc / len(pos)
                zc = zc / len(pos)
                return xc, yc, zc

        #getting initial ca trace in pdb format
        def toPDB(self, pos, suffix):
                f = open(self.output + '/' + self.aa_file + suffix, 'w') #write pdb in this file
                atm_srl = 7 #atom serial number initially assumed CA is at 7
                res_seq = 1 #residue sequence number, initially 1
                for i in range(len(pos)):
                        f.write("ATOM  ") #1-6
                        f.write(" " * (5 - len(str(atm_srl))) + str(atm_srl)) #7-11 right alignment
                        f.write("  CA ") #12-16 left alignment* to do: adjust for other atoms
                        f.write(" ") #17 to do next
                        f.write(self.res_name[i]) #18-20 right alignment
                        f.write(" ")
                        #f.write("A") #22 chain identifier, initaially none because this pdb is the input of modeller, where chain is not in consideration in the pipeline
                        f.write(" ") #change log: added ' ' (space) instead of chain identifier assuming single chaing pdb
                        f.write(" " * (4 - len(str(res_seq))) + str(res_seq)) #23-26 right alignment
                        f.write(" ") #to do next
                        f.write("   ")
                        f.write(" " * (8 - len(str(format(pos[i][0],'.3f')))) + str(format(pos[i][0],'.3f'))) #31-38 x coordinate right alignment
                        f.write(" " * (8 - len(str(format(pos[i][1],'.3f')))) + str(format(pos[i][1],'.3f'))) #39-46 y coordinate right alignment
                        f.write(" " * (8 - len(str(format(pos[i][2],'.3f')))) + str(format(pos[i][2],'.3f'))) #47-54 z coordinate right alignment
                        f.write("\n")
                        atm_srl = atm_srl + 7
                        res_seq = res_seq + 1
                f.close()

        #difference between two distance map (kept for further usage)
        def STRESS(self, D1, D2):
                s = 0
                for i in range(len(D1)):
                        for j in range(i):
                                s = s + ((D1[i][j] - D2[i][j]) * (D1[i][j] - D2[i][j]))
                #print ("stress is " + str(s))
                return s

        #returns stress between two pbds
        def comparepdb(self, pdb1, pdb2):
                fpdb1 = open(pdb1, 'r')
                fpdb2 = open(pdb2, 'r')
                pdb1lines = fpdb1.readlines()
                pdb2lines = fpdb2.readlines()
                pos1 = []
                pos2 = []
                for line in pdb1lines:
                        if(line[:4] == "ATOM" and line[12:16].strip() == "CA"):
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                pos1.append([x, y, z])
                for line in pdb2lines:
                        if(line[:4] == "ATOM" and line[12:16].strip() == "CA"):
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                pos2.append([x, y, z])

                n1 = len(pos1)
                n2 = len(pos2)
                if (n1 != n2):
                        return -1
                D1 = [[0 for x in range(n1)] for y in range(n1)]
                D2 = [[0 for x in range(n2)] for y in range(n2)]

                for i in range(len(pos1)):
                        for j in range(len(pos1)):
                                d = math.sqrt((pos1[i][0] - pos1[j][0]) ** 2 + (pos1[i][1] - pos1[j][1]) ** 2 + (pos1[i][2] - pos1[j][2]) ** 2)
                                D1[i][j] = d
                for i in range(len(pos2)):
                        for j in range(len(pos2)):
                                d = math.sqrt((pos2[i][0] - pos2[j][0]) ** 2 + (pos2[i][1] - pos2[j][1]) ** 2 + (pos2[i][2] - pos2[j][2]) ** 2)
                                D2[i][j] = d
                s = self.STRESS(D1, D2)/len(pos1)

                fpdb1.close()
                fpdb2.close()
                return s



        #takes 3D coordinates and returns distance matrix (kept for further usage)                      
        def coordToDist(self, pos):

                Dist = [[0 for x in range(len(pos))] for y in range(len(pos))]
                for i in range(len(pos)):
                        for j in range(i):
                                Dist[i][j] = Dist[j][i] = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                return Dist

        #correcting the corrdinates according to contact map 
        def correctCM(self, pos):
                for i in range(len(self.CM)):
                        #print (pos) #new
                        maxu = 0

                        F = [0.0, 0.0, 0.0] #force vector
                        D0 = float('inf') #minimum distance within true no-contact cm[i,j] = 0
                        D1 = 0.0 #maximum distance within true contact cm[i,j] = 1
                        ri = 0 #radius of the sphere centered to pos[i], pos[i] may be relocated to the surface of this sphere
                        for j in range(len(self.CM)):
                                if (abs(i-j) <= 5):
                                        continue
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))


                                if (self.CM[i][j] == 1 and dij > self.pair_u[(i,j)]): #new
                                        #r = dij - self.t
                                        #changing the vector as inconsistancy found
                                        F[0] = F[0] - (pos[i][0] - pos[j][0]) / dij
                                        F[1] = F[1] - (pos[i][1] - pos[j][1]) / dij
                                        F[2] = F[2] - (pos[i][2] - pos[j][2]) / dij


                                elif (self.CM[i][j] == 0 and dij <= self.t):
                                        #r = self.t - dij
                                        #changing the vector as inconsistancy found                                             
                                        F[0] = F[0] + (pos[i][0] - pos[j][0]) / dij
                                        F[1] = F[1] + (pos[i][1] - pos[j][1]) / dij
                                        F[2] = F[2] + (pos[i][2] - pos[j][2]) / dij
                                else:
                                        if (self.CM[i][j] == 0):
                                                D0 = min(D0, dij)
                                        else:
                                                D1 = max(D1, dij)
                                                maxu = max(maxu, self.pair_u[(i,j)])   
                        ri = min(D0 - self.t, maxu - D1)
                        Fmod = math.sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2])
                        if (Fmod == 0):
                                Fmod = 1
                        pos[i][0] = pos[i][0] + F[0] * (ri / Fmod)
                        pos[i][1] = pos[i][1] + F[1] * (ri / Fmod)
                        pos[i][2] = pos[i][2] + F[2] * (ri / Fmod)
                return pos

        #correcting distance of two consecutive residues
        def correctConSC(self, pos):
                #this program corrects consecutive residues distance, if di,i-1 < 3.8 then corrects
                for i in range(1, len(self.CM)):
                        F = [0, 0, 0]
                        ri = 0
                        for j in range(i-1, i, 2): #consider only previous residue
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                if(dij < 3.8):
                                        #move far

                                        F[0] = F[0] + (pos[i][0] - pos[j][0]) / dij
                                        F[1] = F[1] + (pos[i][1] - pos[j][1]) / dij
                                        F[2] = F[2] + (pos[i][2] - pos[j][2]) / dij
                                        ri = (3.8 - dij)
                                elif(dij > 3.8):

                                        #move close

                                        F[0] = F[0] - (pos[i][0] - pos[j][0]) / dij
                                        F[1] = F[1] - (pos[i][1] - pos[j][1]) / dij
                                        F[2] = F[2] - (pos[i][2] - pos[j][2]) / dij
                                        ri = (dij - 3.8)



                        Fmod = math.sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2])
                        if (Fmod == 0):
                                Fmod = 1


                        pos[i][0] = pos[i][0] + F[0] * (ri / Fmod)
                        pos[i][1] = pos[i][1] + F[1] * (ri / Fmod)
                        pos[i][2] = pos[i][2] + F[2] * (ri / Fmod)


                return pos

        #correcting the distance of two close residues
        def correctClose(self, pos):
                for i in range(len(self.CM)):
                        F = [0, 0, 0]
                        ri = 9999

                        for j in range(len(self.CM)):
                                if (j == i):
                                        continue
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                if (dij < 3.5):
                                        F[0] = F[0] + (pos[i][0] - pos[j][0]) / dij
                                        F[1] = F[1] + (pos[i][1] - pos[j][1]) / dij
                                        F[2] = F[2] + (pos[i][2] - pos[j][2]) / dij
                                        ri = min((3.5 - dij), ri)

                        if (ri == 9999):
                                ri = 0
                        Fmod = math.sqrt(F[0] * F[0] + F[1] * F[1] + F[2] * F[2])
                        if (Fmod == 0):
                                Fmod = 1

                        pos[i][0] = pos[i][0] + F[0] * (ri / Fmod)
                        pos[i][1] = pos[i][1] + F[1] * (ri / Fmod)
                        pos[i][2] = pos[i][2] + F[2] * (ri / Fmod)
                return pos



        #the following two methods are for making the pdb coordinates compact
        def perturb(self, pos, eps):
                for i in range(len(self.CM)):

                        for j in range(len(self.CM)):

                                if (abs(i-j) < 5): #new
                                        uu = self.t #new
                                elif (self.CM[i][j] == 1):
                                        uu = self.pair_u[(i,j)]

					
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                if (self.CM[i][j] == 1 and uu - eps < dij <= uu):
                                        #close to thresh hold, make i,j closer
                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * eps / 1.0
                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * eps / 1.0
                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * eps / 1.0


                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * eps / 1.0
                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * eps / 1.0
                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * eps / 1.0
                                elif(self.t < dij < self.t + eps and self.CM[i][j] == 0):
                                        #close to thresh hold, make i,j farthur 
                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * eps / 1.0
                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * eps / 1.0
                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * eps / 1.0


                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * eps / 1.0
                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * eps / 1.0
                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * eps / 1.0

                return pos

        def perturbCentroid(self, pos, centroid, eps):
                for i in range(len(pos)):
                        dij = math.sqrt((pos[i][0] - centroid[0]) ** 2 + (pos[i][1] - centroid[1]) ** 2 + ( pos[i][2] - centroid[2]) ** 2)
                        # make closer
                        pos[i][0] = pos[i][0] - (pos[i][0] - centroid[0]) * dij * eps
                        pos[i][1] = pos[i][1] - (pos[i][1] - centroid[1]) * dij * eps
                        pos[i][2] = pos[i][2] - (pos[i][2] - centroid[2]) * dij * eps
                return pos

        #calculating Rg given 3D coordinates
        def calculateRg(self, pos, centroid):
                dij = 0
                for i in range(len(pos)):
                        dij += ((pos[i][0] - centroid[0]) ** 2 + (pos[i][1] - centroid[1]) ** 2 + ( pos[i][2] - centroid[2]) ** 2)
                dij = math.sqrt(dij / (len(pos)))
                #print (len(pos))
                return dij
        #correcting coordinates based on secondary structure information                                
        def perturbCorrect(self, pos, eps):
                for i in range(len(self.CM)):
                        for j in range(len(self.CM)):
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                m = max(i, j)
                                if (abs(i - j) == 2):

                                        if (self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H'):
                                                if (dij > 6.0):
                                                        #move close

                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 6.0) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 6.0) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 6.0) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 6.0) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 6.0) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 6.0) / 2.0


                                                elif (dij < 5.0):
                                                        #move far


                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (5.0 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (5.0 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (5.0 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (5.0 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (5.0 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (5.0 - dij) / 2.0




                                        if (self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E'):
                                                if (dij > 7.3):
                                                        #move close

                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 7.3) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 7.3) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 7.3) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 7.3) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 7.3) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 7.3) / 2.0


                                                if (dij < 6.1):
                                                        #move far

                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (6.1 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (6.1 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (6.1 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (6.1 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (6.1 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (6.1 - dij) / 2.0




                                if (abs(i - j) == 3):
                                        if (self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H' and self.secondary[0][m-2] == 'H'):
                                                if (dij > 5.8):
                                                        #move close
                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 5.8) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 5.8) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 5.8) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 5.8) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 5.8) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 5.8) / 2.0


                                                if (dij < 4.8):
                                                        #move far
                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (4.7 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (4.7 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (4.7 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (4.7 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (4.7 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (4.7 - dij) / 2.0



                                        if (self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E' and self.secondary[0][m-2] == 'E'):
                                                if (dij > 10.8):
                                                        #move close
                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 10.8) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 10.8) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 10.8) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 10.8) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 10.8) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 10.8) / 2.0


                                                if (dij < 9):
                                                        #move far
                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (9 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (9 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (9 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (9 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (9 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (9 - dij) / 2.0
                                if (abs(i - j) == 4):
                                        if (self.secondary[0][m] == 'H' and self.secondary[0][m-1] == 'H' and self.secondary[0][m-2] == 'H' and self.secondary[0][m-3] == 'H'):
                                                if (dij > 7.0):
                                                        #move close
                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 7.0) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 7.0) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 7.0) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 7.0) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 7.0) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 7.0) / 2.0
                                                if (dij < 5.8):
                                                        #move far
                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (5.8 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (5.8 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (5.8 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (5.8 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (5.8 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (5.8 - dij) / 2.0


                                        if (self.secondary[0][m] == 'E' and self.secondary[0][m-1] == 'E' and self.secondary[0][m-2] == 'E' and self.secondary[0][m-3] == 'E'):

                                                if (dij > 12.5):

                                                        #move close
                                                        pos[i][0] = pos[i][0] - (pos[i][0] - pos[j][0]) / dij * (dij - 12.5) / 2.0
                                                        pos[i][1] = pos[i][1] - (pos[i][1] - pos[j][1]) / dij * (dij - 12.5) / 2.0
                                                        pos[i][2] = pos[i][2] - (pos[i][2] - pos[j][2]) / dij * (dij - 12.5) / 2.0

                                                        pos[j][0] = pos[j][0] - (pos[j][0] - pos[i][0]) / dij * (dij - 12.5) / 2.0
                                                        pos[j][1] = pos[j][1] - (pos[j][1] - pos[i][1]) / dij * (dij - 12.5) / 2.0
                                                        pos[j][2] = pos[j][2] - (pos[j][2] - pos[i][2]) / dij * (dij - 12.5) / 2.0
                                                if (dij < 12.2):

                                                        #move far
                                                        pos[i][0] = pos[i][0] + (pos[i][0] - pos[j][0]) / dij * (12.2 - dij) / 2.0
                                                        pos[i][1] = pos[i][1] + (pos[i][1] - pos[j][1]) / dij * (12.2 - dij) / 2.0
                                                        pos[i][2] = pos[i][2] + (pos[i][2] - pos[j][2]) / dij * (12.2 - dij) / 2.0

                                                        pos[j][0] = pos[j][0] + (pos[j][0] - pos[i][0]) / dij * (12.2 - dij) / 2.0
                                                        pos[j][1] = pos[j][1] + (pos[j][1] - pos[i][1]) / dij * (12.2 - dij) / 2.0
                                                        pos[j][2] = pos[j][2] + (pos[j][2] - pos[i][2]) / dij * (12.2 - dij) / 2.0

                return pos

        #checking if the 3D structure is consistant with contact map
        def consistant(self, pos):
                for i in range(len(self.CM)):
                        for j in range(i):
                                if (i - j <= 5):
                                        continue
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                if(self.CM[i][j] == 1 and dij > self.pair_u[(i,j)]):
                                        return False
                                elif(self.CM[i][j] == 0 and dij <= self.t):
                                        return False
                return True

        #returns the contact error, no contact error and list of pairs mismatched with given contact map
        def cmError(self, pos):
                e1 = 0
                e2 = 0
                elist = []
                for i in range(len(self.CM)):
                        for j in range(i):
                                if (i - j <= 5):
                                        continue
                                dij = math.sqrt((pos[i][0] - pos[j][0]) * (pos[i][0] - pos[j][0]) + (pos[i][1] - pos[j][1]) * (pos[i][1] - pos[j][1]) + (pos[i][2] - pos[j][2]) * (pos[i][2] - pos[j][2]))
                                if(self.CM[i][j] == 1 and dij > self.pair_u[(i,j)]):
                                        e1 = e1 + ((dij - self.pair_u[(i,j)]) ** 2)
                                        elist.append((i,j))
                                elif(self.CM[i][j] == 0 and dij <= self.t):
                                        e2 = e2 + ((dij - self.t) ** 2)
                                        elist.append((i,j))
                return e1, e2, elist


############# The following 4 functions are takeing the ss information, generating restraints in angular space, converting the restraints to distance which is finally used in initial distance map refinement
##Credit: Dr. Debswapna Bhattacharya##
        #improves distance map based on given secondary structure
        def ChunksImproveDist(self, D):
                #this method takes initial distance map D and imporves based on knwon alpha helix or beta distance info
                helixAlphaI = math.radians(50) #ideal mean from p-sea
                helixTauI = math.radians(89) #ideal mean from p-sea
                betaSheetAlpha = math.radians(-170)
                betaSheetTau = math.radians(124)

                for i,j in self.helices:

                        if (abs(j - i) < 3 or i-2 < 0):
                                continue

                        segment = []
                        segment.append([0, -3.8, 0])
                        segment.append([0, 0, 0])
                        segment.append([3.8, 0, 0])

                        for k in range(j-i):
                                segpos = self.setCoordinate(segment, helixAlphaI, helixTauI)

                                for l in range(k):

                                        d = math.sqrt((segpos[k][0] - segpos[l][0]) ** 2 + (segpos[k][1] - segpos[l][1]) ** 2 + (segpos[k][2] - segpos[l][2]) ** 2)

                                        D[k+i][l+i] = D[l+i][k+i] = d


                for i,j in self.sheets:

                        if (abs(j - i) < 3 or i-2 < 0):
                                continue

                        segment = []
                        segment.append([0, -3.8, 0])
                        segment.append([0, 0, 0])
                        segment.append([3.8, 0, 0])

                        for k in range(j-i):
                                segpos = self.setCoordinate(segment, betaSheetAlpha, betaSheetTau)

                                for l in range(k):

                                        d = math.sqrt((segpos[k][0] - segpos[l][0]) ** 2 + (segpos[k][1] - segpos[l][1]) ** 2 + (segpos[k][2] - segpos[l][2]) ** 2)

                                        D[k+i][l+i] = D[l+i][k+i] = d

                return D

        #takes pos (3d cordinate list), alpha and tau angle and returns pos with an appended corrdinate(x, y, z)
        def setCoordinate(self, pos, alpha, tau):
                x = 0
                y = 1
                z = 2

                normw = 3.8
                c2 = pos[len(pos) - 1]
                c1 = pos[len(pos) - 2]
                c0 = pos[len(pos) - 3]


                u1 = c1[x] - c0[x]
                u2 = c1[y] - c0[y]
                u3 = c1[z] - c0[z]

                v1 = c2[x] - c1[x]
                v2 = c2[y] - c1[y]
                v3 = c2[z] - c1[z]


                norm = math.sqrt((v1 * v1) + (v2 * v2) + (v3 * v3))

                v1 /= norm
                v2 /= norm
                v3 /= norm


                pvuv1 = u2 * v3 - u3 * v2
                pvuv2 = u3 * v1 - u1 * v3
                pvuv3 = u1 * v2 - u2 * v1

                norm = math.sqrt((pvuv1 * pvuv1) + (pvuv2 * pvuv2) + (pvuv3 * pvuv3))
                pvuv1 /= norm
                pvuv2 /= norm
                pvuv3 /= norm

                pvvuv1 = v2 * pvuv3 - v3 * pvuv2
                pvvuv2 = v3 * pvuv1 - v1 * pvuv3
                pvvuv3 = v1 * pvuv2 - v2 * pvuv1

                norm = math.sqrt((pvvuv1 * pvvuv1) + (pvvuv2 * pvvuv2) + (pvvuv3 * pvvuv3))
                pvvuv1 /= norm
                pvvuv2 /= norm
                pvvuv3 /= norm

                nca = math.cos(alpha)
                nsa = math.sin(alpha)
                nct = math.tan(tau - (3.1416/2))

                u1 = nca * (-pvvuv1) + nsa * pvuv1 + v1 * nct
                u2 = nca * (-pvvuv2) + nsa * pvuv2 + v2 * nct
                u3 = nca * (-pvvuv3) + nsa * pvuv3 + v3 * nct

                norm = math.sqrt((u1 * u1) + (u2 * u2) + (u3 * u3))
                u1 = u1 * normw/norm
                u2 = u2 * normw/norm
                u3 = u3 * normw/norm

                newx = u1 + c2[x]
                newy = u2 + c2[y]
                newz = u3 + c2[z]

                X = []
                X.append(newx)
                X.append(newy)
                X.append(newz)

                pos.append(X)

                return pos

        def getAngle(self, x1, x2, x3, y1, y2, y3, z1, z2, z3):
                acc = 0.0
                d1 = 0
                d2 = 0
                acc = (x2 - x1) * (x2 - x3) + (y2 - y1) * (y2 - y3) + (z2 - z1) * (z2 - z3)
                d1 = math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
                d2 = math.sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
                if (d1 == 0 or d2 == 0):
                        return 0
                acc = acc / (d1 * d2);
                if (acc > 1.0):
                        acc = 1.0
                elif (acc < -1.0):
                        acc = -1.0;
                acc = math.acos(acc);
                return acc;


        def getDihedral(self, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4):

                #getDifferences
                #d1 = difference(p1, p2)
                qx = x2 - x1;
                qy = y2 - y1;
                qz = z2 - z1;

                #d2 = difference(p3, p2)
                rx = x2 - x3;
                ry = y2 - y3;
                rz = z2 - z3;

                #d3 = difference(p4, p3)
                sx = x3 - x4;
                sy = y3 - y4;
                sz = z3 - z4;

                #get cross product#
                #cp1 = crossProduct(d1, d2)
                tx = qy * rz - qz * ry;
                ty = qz * rx - qx * rz;
                tz = qx * ry - qy * rx;

                #cp2 = crossProduct(d3, d1)

                ux = sy * rz - sz * ry;
                uy = sz * rx - sx * rz;
                uz = sx * ry - sy * rx;

                #cp3 = crossProduct(cp2, cp1)
                vx = uy * tz - uz * ty;
                vy = uz * tx - ux * tz;
                vz = ux * ty - uy * tx;

                #get dot product#
                #dot = getDotProduct(cp3, d2)
                w = vx * rx + vy * ry + vz * rz;

                #get angle#
                zx = zy = zz = 0.0;
                accTau = 0.0;
                accTau = (zx - tx) * (zx - ux) + (zy - ty) * (zy - uy) + (zz - tz) * (zz - uz);
                d1Tau = math.sqrt((tx - zx) * (tx - zx) + (ty - zy) * (ty - zy) + (tz - zz) * (tz - zz));
                d2Tau = math.sqrt((ux - zx) * (ux - zx) + (uy - zy) * (uy - zy) + (uz - zz) * (uz - zz));

                if (d1Tau == 0 or d2Tau == 0):
                        return 0

                accTau = accTau / (d1Tau * d2Tau);

                if (accTau > 1.0):
                        accTau = 1.0;
                elif (accTau < -1.0):
                        accTau = -1.0;
                accTau = math.acos(accTau);
                #end of getAngle#
                if (w < 0):
                        accTau = -accTau;
                return accTau

        #showing angle (for future use)                                                                                                                                                                                                                                  
        def showAngles(self, pos):
                x = 0
                y = 1
                z = 2
                for i in range(1, len(pos) - 2):
                        print ('res no ' + str(i) + ', secondary: ' + str(self.secondary[0][i]))
                        print ('dihedral angle : ' + str(180 / 3.1416 * self.getDihedral(pos[i-1][x], pos[i][x], pos[i+1][x], pos[i+2][x], pos[i-1][y], pos[i][y], pos[i+1][y], pos[i+2][y], pos[i-1][z], pos[i][z], pos[i+1][z], pos[i+2][z])))
                        print ('tau angle : ' +str(180 / 3.1416 * self.getAngle(pos[i-1][x], pos[i][x], pos[i+1][x], pos[i-1][y], pos[i][y], pos[i+1][y], pos[i-1][z], pos[i][z], pos[i+1][z])))


        #making alignment file for modeller 
        def same_ali(self, aa_file, model, output):

                f = open(aa_file, 'r')
                fo = open(output, 'w')
                lines = f.readlines()
                aa = lines[1].strip()
                fo.write('>P1;query\n')
                fo.write('sequence:query:1 :  :   :  : q1 : q2 : q3 :\n')
                fo.write(aa)
                fo.write('*')
                fo.write('\n')
                fo.write('>P1;')
                fo.write(model)
                fo.write('\n')
                fo.write('structureX:')
                fo.write(model)
                fo.write('  :1 :  :   :  : t1 : t2 : t3 :\n')
                fo.write(aa)
                fo.write('*')
                f.close()
                fo.close()

          #making multiple alignment for modeller
        def multiali(self, aa_file, inputmodel, output, totalmodels):

                f = open(aa_file, 'r')
                fo = open(output, 'w')
                lines = f.readlines()
                aa = lines[1].strip()
                fo.write('>P1;query\n')
                fo.write('sequence:query:1 :  :   :  : q1 : q2 : q3 :\n')
                fo.write(aa)
                fo.write('*')
                fo.write('\n')

                for i in range(1, totalmodels + 1):
                        fo.write('>P1; ' + inputmodel + str(i) + '.B99990001')
                        fo.write('\n')
                        fo.write('structureX:')
                        fo.write(inputmodel + str(i) + '.B99990001')
                        fo.write('  :1 :  :   :  : t1 : t2 : t3 :\n')
                        fo.write(aa)
                        fo.write('*')
                        fo.write('\n')
                f.close()
                fo.close()



        #format final pdb
        def format_casp_pdb(self, inpdb, outpdb):
                with open(inpdb, "r") as rFile:
                        lines = []
                        for line in rFile:
                                if(len(line.split())>0 and line.split()[0] == "ATOM"):
                                        lines.append(line)

                outFile = open(outpdb, "w")
                #outFile.write("PFRMAT TS\n")
                #outFile.write("TARGET "+targetName+"\n")
                #outFile.write("AUTHOR 0000-0000-0000\n")
                #outFile.write("REMARK Residue error estimate is CA-CA distance in Angstroms\n")
                #outFile.write("METHOD iCon3D"+"\n")
                #outFile.write("MODEL  1\n")
                #outFile.write("PARENT "+"N/A\n")

                #starting residue sequence number##
                startResNo = lines[0][22:(22+4)]
                #print(startResNo)

                counter=0
                prev_residueSeqNum=startResNo

                for k in range (len(lines)):

                        atom = lines[k][0:(0+4)]
                        atomSerial = lines[k][6:(6+5)]
                        atomName = lines[k][12:(12+4)]
                        residueName = lines[k][17:(17+3)]
                        residueSeqNum = lines[k][22:(22+4)]
                        xCoord = lines[k][30:(30+8)]
                        yCoord = lines[k][38:(38+8)]
                        zCoord = lines[k][46:(46+8)]
                        occupency = lines[k][54:(54+6)]
                        tempFact = lines[k][60:(60+6)]


                        if(prev_residueSeqNum!=residueSeqNum):
                                counter=counter+1
                        #tempFact = (Decimal(ei[counter]).quantize(TWOPLACES))
                        prev_residueSeqNum=residueSeqNum

                        outFile.write(atom+"  "+'%5s'%atomSerial+" "+'%4s'%atomName+" "+'%3s'%residueName+"  "+
                                  '%4s'%residueSeqNum+"    "+'%8s'%xCoord+'%8s'%yCoord+'%8s'%zCoord+'%6s'%occupency+'%6s'%tempFact+'%15s'%"\n")

                outFile.write("TER\n")
                outFile.write("END")
                outFile.close()


        #generating multialignment modeller script
        def generateModellerMulti(self):
                f = open('model_multi.py', 'w')
                f.write('from modeller import * \n')
                f.write('from modeller.automodel import * \n')
                f.write('import os,sys,optparse\n')

                f.write('import math\n')
                f.write('parser=optparse.OptionParser()\n')

                f.write('parser.add_option(\'-f\', dest=\'alignment_file\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'alignment file\')\n')
                f.write('parser.add_option(\'--target\', dest=\'target\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'name of target\') \n')

                f.write('(options,args) = parser.parse_args()\n')

                f.write('alignment_file = options.alignment_file\n')
                f.write('target = options.target\n')
                f.write('log.verbose() \n')
                f.write('env = environ() \n')
                f.write('a = automodel(env,  \n')
                f.write('       alnfile=alignment_file, \n')
                f.write('       knowns= (target + str(\'1.B99990001\') , target + str(\'2.B99990001\'), target + str(\'3.B99990001\'), target + str(\'4.B99990001\'), target + str(\'5.B99990001\')), \n')#target' + str('3.B99990001.pdb') + ', target' + str('4.B99990001.pdb') + ', target' + str('5.B99990001.pdb') + ')\n')
                f.write('       sequence = \'query\')\n')
                f.write('a.starting_model = 1\n')
                f.write('a.ending_model = 5 \n')
                f.write('a.make()')

        #generating modeller script to run modeller
        def generateModellerScript(self):



                f = open('model_ss_cm.py', 'w')

                f.write('from modeller import *\n')

                f.write('from modeller.automodel import *\n')

                f.write('import os,sys,optparse\n')

                f.write('import math\n')
                f.write('parser=optparse.OptionParser()\n')

                f.write('parser.add_option(\'-f\', dest=\'alignment_file\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'alignment file\')\n')
                f.write('parser.add_option(\'--ss\', dest=\'ss\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'secondary structure file\') \n')

                f.write('parser.add_option(\'--rr\', dest=\'rr\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'contact map in standard rr format\') \n')

                f.write('parser.add_option(\'--target\', dest=\'target\', \n')

                f.write('        default= \'\',    #default empty!\' \n')

                f.write('        help= \'name of target\') \n')

                f.write('parser.add_option(\'--template\', dest=\'template\', \n')

                f.write('        default= \'\',    #default empty! \n')

                f.write('        help= \'name of template\') \n')
                f.write('parser.add_option(\'--ctype\', dest=\'ctype\', \n')

                f.write('        default= \'ca\',    #default empty! \n')

                f.write('        help= \'ca or cb. default ca\') \n')

                f.write('(options,args) = parser.parse_args()\n')


                f.write('alignment_file = options.alignment_file\n')

                f.write('target = options.target\n')

                f.write('template = options.template\n')

                f.write('ss = options.ss\n')

                f.write('rr = options.rr\n')

                f.write('ctype = options.ctype\n')

                f.write('class MyModel(automodel):\n')
                f.write('       def special_restraints(self, aln): \n')
                f.write('               rsr = self.restraints \n')
                f.write('               at = self.atoms \n')

                f.write('               fss = open(ss, \'r\') \n')
                f.write('               self.secondary = fss.read() \n')
                f.write('               self.secondary = self.secondary.split() #secondary structure is now in self.secondary[0] \n')
                f.write('               fss.close() \n')
                f.write('               print (self.secondary[0]) \n')
                f.write('               frr = open(rr, \'r\') \n')
                f.write('               rrlines = frr.readlines() \n')
                f.write('               faln = open(alignment_file, \'r\') \n')
                f.write('               alnlines = faln.readlines() \n')
                f.write('               n = len(alnlines[2]) - 2 \n')
                f.write('               pair_u = {} \n')
                f.write('               pair_l = {} \n')

                f.write('               self.CM = [[0 for x in range(n)] for y in range(n)] \n')

                f.write('               for line in rrlines[1:]: \n')

                f.write('                       line = line.strip().split() \n')
                f.write('                       self.CM[int(line[0]) - 1][int(line[1]) - 1] = self.CM[int(line[1]) - 1][int(line[0]) - 1] = 1   \n')


                f.write('                       if (float(line[4]) <= 0.85): \n')
                f.write('                               self.CM[int(line[0]) - 1][int(line[1]) - 1] = self.CM[int(line[1]) - 1][int(line[0]) - 1] = -1 \n')

                f.write('                       pair_l[(int(line[0]) - 1, int(line[1]) - 1)] = pair_l[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[2]) \n')
                f.write('                       pair_u[(int(line[0]) - 1, int(line[1]) - 1)] = pair_u[(int(line[1]) - 1, int(line[0]) - 1)] = float(line[3]) \n')

                f.write('               mdl = template + \'.B99990001.pdb\' \n')
                f.write('               if os.path.exists(mdl): \n')
                f.write('                       print (\'changing\') \n')
                f.write('                       os.system(\'mv \' + template + \'.B99990001.pdb \' + template + \'.pdb\') \n')
                f.write('               fpdb = open(template + \'.pdb\', \'r\') \n')
                f.write('               fpdblines = fpdb.readlines() \n')
                f.write('               pos = [] \n')
                f.write('               if(ctype == \'cb\' or ctype == \'CB\'): \n')
                f.write('                        ct = \"CB\" \n')
                f.write('               else: \n')
                f.write('                        ct = \"CA\" \n')
                f.write('               for line in fpdblines: \n')
                f.write('                        if(line[:4] == \"ATOM\" and line[12:16].strip() == ct): \n')
                f.write('                                x = float(line[30:38].strip()) \n')
                f.write('                                y = float(line[38:46].strip()) \n')
                f.write('                                z = float(line[46:54].strip()) \n')
                f.write('                                pos.append([x, y, z]) \n')

                f.write('               for i in range(len(pos)): \n')
                f.write('                       for j in range(i+6, len(pos)): \n')
                f.write('                               if (abs(i - j) < 6): \n')
                f.write('                                       continue \n')
                f.write('                               d = math.sqrt((pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2 + (pos[i][2] - pos[j][2]) ** 2) \n')
                f.write('                               if (self.CM[i][j] == 1 and d > pair_u[(i,j)]): \n')
                f.write('                                       firstres = ct + \':\' + str(i+1) \n')
                f.write('                                       secondres = ct + \':\' + str(j+1) \n')
                f.write('                                       if(alnlines[2][i] == \'G\'): \n')
                f.write('                                               firstres = \'CA:\' + str(i+1) \n')

                f.write('                                       if(alnlines[2][j] == \'G\'): \n')
                f.write('                                               secondres = \'CA:\' + str(i+1) \n')


                f.write('                                       rsr.add(forms.lower_bound(group=physical.xy_distance, \n')
                f.write('                                              feature=features.distance(at[firstres], \n')
                f.write('                                                                at[secondres]), \n')
                f.write('                                              mean=3.6, stdev=.1)) \n') #can change mean and stdev


                f.write('                                       rsr.add(forms.upper_bound(group=physical.xy_distance, \n')
                f.write('                                                       feature=features.distance(at[firstres], \n')
                f.write('                                                                at[secondres]), \n')
                f.write('                                                      mean=pair_u[(i,j)], stdev=0.1)) \n') #can change mean and stdev

                f.write('                               elif (self.CM[i][j] == 0 and d < '+str(self.t)+'): \n')
                f.write('                                        firstres = ct + \':\' + str(i+1) \n')
                f.write('                                        secondres = ct + \':\' + str(j+1) \n')

                f.write('                                        if(alnlines[2][i] == \'G\'): \n')
                f.write('                                               firstres = \'CA:\' + str(i+1) \n')

                f.write('                                        if(alnlines[2][j] == \'G\'): \n')
                f.write('                                               secondres = \'CA:\' + str(i+1) \n')

                f.write('                                        rsr.add(forms.lower_bound(group=physical.xy_distance, \n')
                f.write('                                               feature=features.distance(at[firstres], \n')
                f.write('                                                                 at[secondres]), \n')
                f.write('                                               mean='+str(self.t-0.1)+', stdev=0.1)) \n') #can change mean and stdev



                f.write('               self.helices = [] \n')
                f.write('               self.sheets = [] \n')
                f.write('               i = 1 \n')
                f.write('               while(i < len(self.secondary[0])): \n')

                f.write('                       if self.secondary[0][i] == \'H\': \n')

                f.write('                               j = i \n')
                f.write('                               while (self.secondary[0][j] == \'H\' and j < len(self.secondary[0]) - 1): \n')
                f.write('                                       j = j + 1 \n')
                f.write('                               if (j - i >= 3): \n')
                f.write('                                       self.helices.append((i, j)) \n')
                f.write('                                       i = j \n')
                f.write('                       i = i + 1 \n')
                f.write('               i = 1 \n')
                f.write('               while(i < len(self.secondary[0])): \n')

                f.write('                       if self.secondary[0][i] == \'E\': \n')

                f.write('                               j = i \n')
                f.write('                               while (self.secondary[0][j] == \'E\' and j < len(self.secondary[0]) - 1): \n')
                f.write('                                       j = j + 1 \n')
                f.write('                               if (j - i >= 3): \n')
                f.write('                                       self.sheets.append((i, j-1)) \n')
                f.write('                                       i = j \n')
                f.write('                       i = i + 1 \n')
                f.write('               for (i,j) in self.helices: \n')
                f.write('                       start = str(i) + \':\' \n')
                f.write('                       end = str(j) + \':\' \n')
                f.write('                       rsr.add(secondary_structure.alpha(self.residue_range(start, end))) \n')
                f.write('                       print (start) \n')
                f.write('                       print (end) \n')
                f.write('               for (i,j) in self.sheets: \n')
                f.write('                       start = str(i) + \':\' \n')
                f.write('                       end = str(j) + \':\' \n')
                f.write('                       rsr.add(secondary_structure.strand(self.residue_range(start, end))) \n')
                f.write('                       print (start) \n')
                f.write('                       print (end) \n')
                f.write('log.verbose() \n')

                #env = environ(rand_seed=-4000) #for making variable model each time
                f.write('env = environ() #consistant models each time \n')

                f.write('a = MyModel(env,alnfile=alignment_file, \n')

                f.write('                knowns=template,sequence=target, \n')

                f.write('                assess_methods=(assess.DOPE)) \n')


                f.write('a.starting_model= 1 \n')

                f.write('a.ending_model = 1 \n')


                f.write('a.make() \n')
                f.write('mdl = template + \'.B99990001.pdb\' \n')
                f.write('if not os.path.exists(mdl): \n')

                f.write('       os.system(\'mv \' + template + \'.pdb \' + template + \'.B99990001.pdb\') \n')



                f.close()

