#!/usr/bin/python

#  DConStruct: Hybridized distance- and contact-based hierarchical protein folding
#
#  Copyright (C) Bhattacharya Laboratory 2019
#
#  DConStruct is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  DConStruct is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with DConStruct.  If not, see <http://www.gnu.org/licenses/>.
#
############################################################################

import os
import sys
import optparse    # for option sorting
import time
import numpy as np
import random
from decimal import *
import lib.dStruct as ds
import lib.cStruct as cs
os.system('export OPENBLAS_NUM_THREADS=1')
# take input arguments
parser = optparse.OptionParser()
parser.add_option('-r', dest='rr',
        default = '',    # default empty
        help = 'rr file in CASP format containing the contact map (mandatory)')
parser.add_option('-a', dest='aa',
        default = '',    # default empty
        help = 'fasta file containing the amino acid sequence (mandatory)')
parser.add_option('-s', dest='ss',
        default = '',    # default empty
        help = 'secondary structure file (mandatory)')
parser.add_option('-m', dest='M',
        default = '',    # default empty
        help = 'MODELLER program path that contains modpy.sh script (mandatory)')
parser.add_option('-o', dest='output',
        default = '',    # default empty
        help = 'existing output directory path (mandatory)')
parser.add_option('-n', dest='no',
        default = '7',    # default 7
        help = 'positive integer to be used as seed (optional); default 7')
parser.add_option('-c', dest='ctype',
        default = 'ca',    # default ca 
        help = 'contact type ca or cb (optional); default ca')
parser.add_option('-x', dest='L',
        default = '10000',    # default all
        help = 'top xL contacts, where L is the sequence length (optional); default all')


(options,args) = parser.parse_args()
rr = options.rr
aa = options.aa
output = options.output
ss = options.ss
no = options.no
ctype = options.ctype
L = options.L
M = options.M


#header

print("\n*************************************************************************")
print("*                            DConStruct                                 *")
print("*  Hybridized distance- and contact-based hierarchical protein folding  *")
print("*  For comments, please email to bhattacharyad@auburn.edu               *")
print("*************************************************************************\n")

def print_usage():
        print("\nUsage: DConStruct.py [options]\n")

        print("Options:")
        print(" -h, --help  show this help message and exit")
        print(" -r RR       rr file in CASP format containing the contact map (mandatory)")
        print(" -a AA       fasta file containing the amino acid sequence (mandatory)")
        print(" -s SS       secondary structure file (mandatory)")
        print(" -m M        MODELLER program path that contains modpy.sh script")
        print("             (mandatory)")
        print(" -o OUTPUT   existing output directory path (mandatory)")
        print(" -n NO       positive integer to be used as seed (optional); default 7")
        print(" -c CTYPE    contact type ca or cb (optional); default ca")
        print(" -x L        top xL contacts, where L is the sequence length (optional);")
        print("             default all\n")


try:
        check_r = int(no)
except ValueError:
        print ('Error! -n must be positive integer. Exiting ...')
        print_usage()
        sys.exit()

if (float(no) <= 0):
        print ('Error! -n must be positive integer. Exiting ...')
        print_usage()
        sys.exit()

seed = np.random.RandomState(seed = int(no)) #giving seed for sklearn MDS state
random.seed(1000 + int(no)) #giving a seed to generate same random values always


#basic input check
if (rr == ''):
        print ('Error! RR file must be provided. Exiting ...')
        print_usage()
        sys.exit()
if (aa == ''):
        print ('Error! Fasta file must be provided. Exiting ...')
        print_usage()
        sys.exit()
if (ss == ''):
        print ('Error! Secondary structure file must be provided. Exiting ...')
        print_usage()
        sys.exit()
if (output == ''):
        print ('Error! Output directory must be provided. Exiting ...')
        print_usage()
        sys.exit()

#existence check
if not os.path.exists(rr):
        print ('Error! No such rr file found. Exiting ...')
        print_usage()
        sys.exit()
if not os.path.exists(aa):
        print ('Error! No such fasta file found. Exiting ...')
        print_usage()
        sys.exit()
if not os.path.exists(ss):
        print ('Error! No such secondary structure file found. Exiting ...')
        print_usage()
        sys.exit()
if not os.path.exists(output):
        print ('Error! No such output directory found. Exiting ...')
        print_usage()
        sys.exit()

#rr, fasta and ss sanity check
frr = open(rr, 'r')
frlines = frr.readlines()
faa = open(aa, 'r')
falines = faa.readlines()
fss = open(ss, 'r')
fslines = fss. readlines()

if (len(falines) != 2):
        print ('Incompatible fasta file. File format should be - first line: fasta name, second line: sequence. Exiting ...')
        print_usage()
        sys.exit()
if (len(fslines) != 1):
        print ('Incompatible ss file. File format should be - first line: ss sequence. Exiting ...')
        print_usage()
        sys.exit()
rrfasta = frlines[0].strip()
aafasta = falines[1].strip()
ssfasta = fslines[0].strip()

if (len(rrfasta) != len(aafasta) or len(ssfasta) != len(aafasta)):
        print ('RR file, fasta file and ss file length mismatch. Exiting ...')
        print_usage()
        sys.exit()
for i in ssfasta:
        if (i == 'H' or i == 'C' or i == 'E'):
                pass
        else:
                print ('Invalid secondary structure. Exiting ...')
                print_usage()
                sys.exit()

frr.close()
faa.close()
fss.close()

#reading paths 
_, aa_file = os.path.split(aa)
aa_file = aa_file.split('.')[0]
o, _ = os.path.split(output)
_, rrfile = os.path.split(rr)
_, ssfile = os.path.split(ss)

if (M == ''):
        print ('Error! MODELLER program full path must be provided. Exiting ...')
        print_usage()
        sys.exit()

modeller = os.path.abspath(M + '/modpy.sh')
#modeller = M + '/modpy.sh'


if (not os.path.exists(modeller)):
        print ('Error! MODELLER not found. Exiting ...')
        print_usage()
        sys.exit()

if (float(L) <= 0):
        print ('Error! -x must be > 0. Exiting ...')
        print_usage()
        sys.exit()

if (ctype != 'ca' and ctype != 'cb'):
        print ('Error! Contact type must be \'ca\' or \'cb\'. Exiting ...')
        #print(ctype)
        print_usage()
        sys.exit()



def main():
        start = time. time()
        map_type = "contact"

        f = open(rr,'r')
        flines = f.readlines()
        l = flines[1].strip().split()
        fl = float(l[2])
        fu = float(l[3])
        for line in flines[2:]:
                lc = line.strip().split()
                if (fl != float(lc[2]) or fu != float(lc[3])):
                        map_type = "distant"
                        break
        if (map_type == "contact"):
                print ("contact map")
                r1 = cs.cStruct(aa, rr, output, ss, no, L, ctype)
        else:
                print ("distance map")
                r1 = ds.dStruct(aa, rr, output, ss, no, L, ctype)
        #r1 = iCon3D() #initialization
        r1.mds3() #initial models 

        for i in range(1, 6):
                file_path = output + '/' + aa_file + 'bestmodel' + str(i) + '.pdb'
                while not os.path.exists(file_path):
                        time.sleep(1)
        
                
                

        #creating temp directory for modeller refinement
        os.system('mkdir ' + output + '/' + aa_file +'temp')
        os.system('cp ' + rr + ' ' + output + '/' + aa_file + 'temp/')
        os.system('cp ' + ss + ' ' + output + '/' + aa_file + 'temp/')

        #alignment
        for i in range(1, 6):
                r1.same_ali(aa, aa_file + 'bestmodel' + str(i),   output + '/' + aa_file +'temp/' + aa_file + str(i) + '.ali')

        for i in range(1, 6):
                file_path =   output + '/' + aa_file +'temp/'+aa_file + str(i) + '.ali'
                while not os.path.exists(file_path):
                        time.sleep(1)
        for i in range(1, 6):
                os.system('cp ' + output + '/' + aa_file + 'bestmodel' + str(i) + '.pdb ' + output + '/' + aa_file +'temp/')
        for i in range(1, 6):
                file_path =  output + '/' + aa_file +'temp/' + aa_file + 'bestmodel' + str(i) + '.pdb'
                while not os.path.exists(file_path):
                        time.sleep(1)
        save_path = os.getcwd()
        os.chdir( output + '/' + aa_file +'temp/')
        r1.generateModellerScript()
        time.sleep(1) #wait to generate the script
        

        #print ('Running MODELLER ...')
        for i in range(1, 6):
                for lp in range(30):
                        os.system(modeller + ' python model_ss_cm.py -f '  + aa_file + str(i) + '.ali --ss ' + ssfile  + ' --rr ' + rrfile + ' --target ' +  aa_file + 'bestmodel' + str(i) + ' --ctype '+ ctype +' --template '  +  aa_file + 'bestmodel' + str(i) + ' > modeller.log')

                        file_path = aa_file + 'bestmodel' + str(i) + '.B99990001.pdb'
                        while not os.path.exists(file_path):
                                time.sleep(1)


        for i in range(1, 6):
                file_path = aa_file + 'bestmodel' + str(i) + '.B99990001.pdb'
                while not os.path.exists(file_path):
                        time.sleep(1)
        os.chdir(save_path)

        r1.multiali(aa, aa_file + 'bestmodel', output + '/' + aa_file +'temp/' + aa_file + '.ali', 5) #multiple alignment in best 5 models
        os.chdir(output + '/' + aa_file + 'temp/')
        r1.generateModellerMulti()
        time.sleep(1) #wait to generate the script
        os.system(modeller + ' python model_multi.py -f '  + aa_file + '.ali  --target ' +  aa_file + 'bestmodel' + ' > modeller_multi.log')
        for i in range(1, 6):
                file_path = 'query.B9999000' + str(i) + '.pdb'
                while not os.path.exists(file_path):
                        time.sleep(1)

        os.chdir(save_path)
	### cleaning and creating final model

	final_model =  output + '/' + aa_file +'temp/' + 'query.B99990001.pdb'
	r1.format_casp_pdb(final_model, output + '/' + aa_file + '_model1.pdb')
        os.system('rm -r ' + output + '/' + aa_file +'temp')
        for i in range(1, 6):
                os.system('rm ' + output + '/' + aa_file + 'bestmodel' + str(i) + '.pdb')

        end = time. time()
        print('Done!\n')
        print ('Finished! Time taken = ' + str(round((end - start)/60, 1)) + ' Minutes\n')
        print('Final models are saved in the output directory')

if __name__ == "__main__":
        main()

