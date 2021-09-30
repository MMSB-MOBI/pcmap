# pcmap : A python module to compute contact map of proteins
pcmap is a PYTHON 3.X library designed to compute pairwise amino acid contacts in protein stuctures. Structures must be provided as PDB coordinates files. Contacts are computed inside a single PDB file or across two PDB files structures. The library can compute one to thousands sets of contacts. Results are produced in JSON format and contacts are encoded in a simple dictionary structure described in the **OUTPUT** section.

## Installation
`pip install pcmap`

or install the package in edition mode, if required

```
git clone pcmap
pip install -e pcmap
```

### Depedencies
pcmap use the [ccmap package](https://github.com/MMSB-MOBI/ccmap) to compute contact maps. This package is a C extension currently available for the following architectures:

* python3.8/OSX.10.14.6
* python3.8/Ubuntu LTS

## Data and testing

The installation folder provides a `data` folder which stores the necessary elements for testing.

```bash
.
+-- README.md
+-- data
    +-- 1A2K_r_u.pdb
    +-- 1A2K_l_u.pdb
    +-- 1A2K_transformations_sample.json
```

Where,
* `1A2K_r_u.pdb` is a single chain protein, later refered as the RECEPTOR.
* `1A2K_r_u.pdb` is a single chain protein, later refered as the LIGAND.
* `1A2K_transformations_sample.json` is a set of transformations described in terms of rotations and translations of the LIGAND protein.

## CLI executable

Once installed, the module can be called as an executable.

### Straight contact map computation on structures

One or several PDB coordinates can be passed to the cli to compute single-body or two-body contact map. Each PDB file defines one body, even if it features many polypeptidic chains.

#### Computing one-body contact map

This will compute the pairwise amino acid contact within the molecule.

##### Single one-body contact map

Just pass the name of the PDB file.

`python -m pcmap dimer data/1A2K_r_u.pdb data/1A2K_l_u.pdb`

#### Many one-body contact maps

Pass a file containing the list of protein as a text file with one PDB file per line:

#### **`sample.lst`**
```txt
data/1A2K_r_u.pdb
data/1A2K_l_u.pdb
```

And pass it to the cli
`python -m pcmap many --structures=sample.lst`

#### Computing many two-body contact map

Again, pass a file containing the list of proteins to compute in tabulated format with two PDB files per line:

#### **`sample_dimer.lst`**

```txt
data/1A2K_r_u.pdb   data/1A2K_l_u.pdb
data/1A2K_r_u.pdb   data/1A2K_l_u.pdb
```

And pass it to the cli
`python -m pcmap many --structures sample_dimer.lst`

### Two-body contact map: applying transformation prior to computation

When dealing with a two body system, it is often convenient to provide the initial conformation of the bodies along with transformations to be applied to generate specific conformations. It is customary to call the first PDB file the RECEPTOR and the second one the LIGAND. The transformation are intended to be applied to LIGAND coordinates, the RECEPTOR remaining unchanged. The transformations are described in terms of rotations and translations of the LIGAND coordinates, using Euler's angles and translation vectors. For mathematical simplicity, LIGAND and RECEPTOR structures can be centered to the origin of the coordinate system.

#### Single two-body contact map
As an example consider the following command:

```shell
python -m pcmap dimer data/1A2K_r_u.pdb data/1A2K_l_u.pdb\
--euler=-1.961,2.066,-2.354 --trans=7.199,16.800,28.799\
--offA=-27.553,-8.229,-80.604 --offB=-67.006,0.11,-77.27
```

1. `1A2K_r_u.pdb` coordinates will be centered onto the origin by the translation `[-27.553,-8.229,-80.604]`.
2. `1A2K_l_u.pdb` coordinates will be centered onto the origin by the translation `[-67.006,0.11,-77.27]`.
3. `1A2K_l_u.pdb` will then be rotated according to `[-1.961,2.066,-2.354]` Euler's angles.
4. `1A2K_l_u.pdb` will then be translated by `[7.199,16.800,28.799]` to obtain the actual conformation.
5. A contact map will be computed across the two structures

###### Obtaining the three dimensional coordinates of the transformed complex

Use the `--apply` flag to generate the corresponding PDB records.
They will be stored into `new_receptor.pdb` and `new_ligand.pdb` files.

#### Many two-body contact maps

When needed, several contact map can be computed by applying a sequence of transformations to the provided RECEPTOR and LIGAND PDB files. Transformations should be described in a JSON file format such as in the following example describing two transformations.

```json
{
    "euler" : [ [-1.96, 2.07, -2.35], [-0.70, 0.95, -0.53] ] ,
    "translation" : [ [7.2, 16.8, 28.8], [21.6, -7.2, -20.4] ],
   "recOffset": [-27.6,-8.2,-80.6],
   "ligOffset": [-67.1,0.1, -77.3],
}
```

Where,

* **euler** references a list of the <span>&alpha;, &beta;, &gamma;</span> Euler angles. Each triplet is transformation specific.
* **translation** references a list of *x,y,z* components of one translation vector. It is transformation specific.
* **recOffset** is the translation vector centering the receptor to coordinates origin. It is common to all transformations.
* **ligOffset** is the translation vector centering the ligand to coordinates origin. It is common to all transformations.

A file example is joined as `data/1A2K_transformations_sample.json`

#### CLI options

##### --apply

Usable on **Single two-body contact map** mode, dumps LIGAND transformed coordinates into  a PDB file named `new_ligand.pdb` (originale RECEPTOR molecule gets written in `new_receptor.pdb`).

##### --dist [default=4.5]

Define the maximal pairwise distance between two heavy atoms to register amino acids contact.

##### --encode [default=False]

If True, contacts are returned as integers. Each integer encoding one pair of atoms/residues positions in contact with this simple formula.

##### --atomic

If True, all atomic contacts are reported.

## PYTHON module

The librairy can be used as a Python module to assemble the pipeline/program of your choice.
First we need to import it.

```python
import pcmap
```

The **pcmap** modules exposes the two following functions:

* *contactMap*
* *contactMapThroughTransform*

### The **contactMap** API

```python
contactMap(proteinA, proteinB=None, **kwargs)
```
 
 The type of the positional parameters controls the function behaviour.

 First parameter can be a PDB file OR a list of PDB files. Second parameter is optional and can also be a PDB file OR a list of PDB files.
 
##### Provided with a PDB file as single parameters:

Compute the internal amino acid contact map of the structure

##### Provided with another PDB file as optional second parameters:

Compute the amino acid contactmap between the two structures

##### Provided with a list of PDB files as single parameter:

Compute the internal amino acid contact maps of each structure indvidually.

##### Provided with a list of PDB files as second parameter:

Compute the amino acid contact map of each pair of structures at identical positions accross the two lists.

#### Examples

```python
# This will compute an internal contact map
c2 = pcmap.contactMap("data/1A2K_r_u.pdb")
# This will compute the contact map across the two structures
c1 = pcmap.contactMap("data/1A2K_r_u.pdb", "data/1A2K_l_u.pdb")
```

### The contactMapThroughTransform API

```python
def contactMapThroughTransform(proteinA,\
    proteinB,\
    eulers, translations,\
    offsetRec,\
    offsetLig,\
    **kwargs):
```

Computes several contact map accross two provided proteins through the applications of provided transormations. Transformations are applied to the SECOND structure.

The *eulers* and *translations* allow to pass even-sized lists of rotation and translation vectors which will be applied to the second structure to generate dimeric conformations.

The *offsetRec* parameter allows to pass a single translation vector to center the first structure barycenter onto the origin.

The *offsetLig* parameter allows to pass a single translation vector to center the second structure barycenter onto the origin.

### module API options

Both **contactMap**and **contactMapThroughTransform** functions share the same set of named parameters:

* *dist*: contact threshold distance[default 4.5]
* *encode*: integer contact encoding[default=False]
* *threadNum*: maximal number of allowed threads[default=8]
* *deserialize*: get results as python dictionary if True, string otherwise [default=True]

## OUTPUT

Amino acid are ranked according to their residue number and chain identifier in the PDB record. Contacts are registred in a one-versus-many fashion: one "root" residue, many "partners" residue, the one residue having the lowest rank. This ensures that parwise contact are registred only once.
The corresponding JSON format is the following:

```json
{'type': 'contactList',
 'data': [ 
        {'root': {'resID': '2 ', 'chainID': 'A'},
            'partners': [ {'resID': '7 ', 'chainID': 'A'},
                        {'resID': '74 ', 'chainID': 'A'}
                    ]
        },
        {'root': {'resID': '9 ', 'chainID': 'A'},
        'partners': [ {'resID': '77 ', 'chainID': 'A'},
                        {'resID': '78 ', 'chainID': 'A'}
                    ]
        }
    ]
}
```

In this example, the residue 2 and 9 of chain A respectively form contacts with residues 7,74 and 77,78 of the same chain.

### Enriched atomic contact map output
In addition to their names and contact distances, the cartesian coordinates of atoms in contacts can be obtained by passing the enrich parameter with a `True` value.

This option is only avaible for the computation of **one single body** contact map. As an example, consider the following call,
```python
contactMap("data/nyxB_monomerB.pdb", atomic=True, enrich=True)
```

which will output (only a sample is shown here):

```json
{'type': 'atomic_rich',
 'data': [
    (('N', 'ALA', '16 ', 'B', 29.392, -64.804, 30.479),
     ('CA', 'ALA', '16 ', 'B', 30.18, -65.518, 29.475),
    1.46),
    (('N', 'ALA', '16 ', 'B', 29.392, -64.804, 30.479),
    ('N', 'PHE', '17 ', 'B', 30.994, -63.28, 28.948),
    2.69),
    (('CA', 'ALA', '16 ', 'B', 30.18, -65.518, 29.475),
    ('N', 'ALA', '16 ', 'B', 29.392, -64.804, 30.479),
    1.46)
    ]
}
```

