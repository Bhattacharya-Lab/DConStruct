# iCon3D

Improved Contact-assisted 3D Protein Folding Protocol

## Installation

The following instructions should work in 64-bit Linux system:
- Make sure you have Python with NumPy, SciPy and Scikit-learn installed. iCon3D has been tested on Python 2.7, but it should run on higher versions as well.
- Install [MODELLER](https://salilab.org/modeller), a license key is required. This can be installed using command `conda install modeller -c salilab`. iCon3D has been tested on MODELLER version 9.20.


## Usage

To see the usage of information, run `python iCon3D_V1.0.py -h`

### File formats and parameters
 
- Fasta (-a): The first line contains the description (header) for the sequence and the second line contains the amino acid sequence. For example, see `./examples/input/T0869.fasta`
- Contacts (-r): The first line contains the amino acid sequence followed by list of contact rows using a five-column format similar to CASP RR format. In each contact row, first two columns are the residue pairs in contact, third and fourth columns are lower and upper bounds of their distance (in Angstrom) respectively, and fifth column is a real number indicating the probability of the two residues being in contact. For example, see `./examples/input/T0869.rr`
- Secondary structure (-s): Single line containing a sequence of 3-state secondary structure (i.e. 'H', 'E' and 'C'). For example, see `./examples/input/T0869.ss`
- Output (-o): Output directory path.
- Random seed (-n): Number to be used as random seed. Default is 7.
- Contact type (-c): To define whether the contact is Ca-Ca or Cb-Cb. Put `-c cb` for Cb-Cb contact.
- Contact cutoff (-x): To select top xL contacts, where L is the sequence length of protein. For example, put `-x 8` to select top 8L contacts.
- Modeller path (-m): Modeller program path that contains `modpy.sh`.

### Test iCon3D

We give an example of running iCon3D on CASP12 FM target T0869.

Create output directory `mkdir output/`. 

Run `python iCon3D_V1.0.py -r examples/input/T0869.rr -a examples/input/T0869.fasta -s examples/input/T0869.ss -o output/ -c cb -x 8 -m your/modeller/path`

Top predicted model will be generated at `output/T0869_model1.pdb`.

iCon3D is very fast. On average it takes only several minutes on a single core to run iCon3D. However, the running time depends on the sequence length of the target protein. For longer targets, iCon3D can take a few hours to complete.

## Data

- The list of PDB chains of FRAGFOLD dataset can be found [here](data/FRAGFOLD_150.txt) 
- The list of Target IDs of CASP12 and CASP13 FM targets can be found [here](data/CASP12_13_FM.txt) 
- The list of PDB chains of 510 non-redundant Membrane protein dataset can be found [here](data/Membrane_510.txt) 
- The list of PDB chains of EVfold dataset can be found [here](data/EVfold_15.txt) 
