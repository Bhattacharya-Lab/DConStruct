# DConStruct

<h2>Hybridized distance- and contact-based hierarchical protein folding</h2>

## Installation

Installing DConStruct is very straightforward. The following instructions should work for 64-bit Linux system:

- Make sure you have Python with NumPy, SciPy and Scikit-learn installed. DConStruct has been tested on Python 2.7.5 (numpy version 1.15.2, scipy version 0.12.1 and scikit-learn version 0.19.2), but it should run on higher versions as well.
- Install [MODELLER](https://salilab.org/modeller), a license key is required. This can be installed using command `conda install modeller -c salilab`. DConStruct has been tested on MODELLER version 9.20.

That's it! DConStruct is ready to be used.

## Usage

To see the usage instructions, run `python DConStruct.py -h`

```

*************************************************************************
*                            DConStruct                                 *
*  Hybridized distance- and contact-based hierarchical protein folding  *
*  For comments, please email to bhattacharyad@auburn.edu               *
*************************************************************************

Usage: DConStruct.py [options]

Options:
  -h, --help  show this help message and exit
  -r RR       rr file in CASP format containing the contact map (mandatory)
  -a AA       fasta file containing the amino acid sequence (mandatory)
  -s SS       secondary structure file (mandatory)
  -m M        MODELLER program path that contains modpy.sh script (mandatory)
  -o OUTPUT   existing output directory path (mandatory)
  -n NO       positive integer to be used as seed (optional); default 7
  -c CTYPE    contact type ca or cb (optional); default cb
  -x L        top xL contacts, where L is the sequence length (optional);
              default 8

```

### File formats and parameters
 
- Amino Acid (-a): The first line contains the header of the target protein and the second line contains the amino acid sequence. For example, see `./examples/input/T0968s2.fasta`
- Contacts (-r): The first line contains the amino acid sequence followed by list of contact rows using a five-column format similar to CASP RR format. In each contact row, first two columns are the residue pairs in contact, third and fourth columns are lower and upper bounds of their distance (in Å) respectively, and fifth column is a real number indicating the probability of the two residues being in contact. For example, see `./examples/input/T0968s2.rr`
- Secondary structure (-s): Single line containing a sequence of 3-state secondary structure (i.e. 'H', 'E' and 'C'). For example, see `./examples/input/T0968s2.ss`
- Modeller path (-m): Modeller program path that contains `modpy.sh` script.
- Output (-o): Output directory path. The directory must exist.
- Random seed (-n): Positive integer to be used as random seed.
- Contact type (-c): To define whether the contact is C<sub>α</sub>–C<sub>α</sub> or C<sub>β</sub>–C<sub>β</sub>. Use `-c ca` for C<sub>α</sub>–C<sub>α</sub> contacts and `-c cb` for C<sub>β</sub>–C<sub>β</sub> contacts.
- Contact cutoff (-x): To select top xL contacts, where L is the sequence length of protein. For example, use `-x 8` to select top 8L contacts.


To run DConStruct with predicted distance-based information, we provide a helper script that can generate distance-based 3-class contact (rr file) from distance histogram (distogram) using rawdistpred.current generated by [DMPfold](https://github.com/psipred/DMPfold). The script is available [here](scripts/dmp2rr.py).

### Test DConStruct

We give an example of running DConStruct on CASP13 FM target T0968s2.

Create an output directory `mkdir output/`. 

Run `python DConStruct.py -r examples/input/T0968s2.rr -a examples/input/T0968s2.fasta -s examples/input/T0968s2.ss -o output/ -c cb -x 8 -m your/modeller/path`

Top predicted model will be generated at `output/T0968s2_model1.pdb`. The predicted 3D model is given [here](examples/output/out.pdb) and the output screen should look like [this](examples/output/log).

DConStruct is reasonably fast. Depending on the sequence length of the target protein, DConStruct takes only a few minutes to a few hours to complete.

## Data

- The list of PDB chains of FRAGFOLD dataset can be found [here](data/FRAGFOLD_150.txt) 
- The list of Target IDs of CASP12 and CASP13 FM targets can be found [here](data/CASP12_13_FM.txt) 
- The list of PDB chains of 510 membrane protein dataset can be found [here](data/Membrane_510.txt) 
- The list of PDB chains of EVfold dataset can be found [here](data/EVfold_15.txt) 

## Cite

If you find DConStruct useful, please cite our paper in bioRxiv available as preprint.
