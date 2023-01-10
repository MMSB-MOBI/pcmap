# pcmap : A python module to compute contact map of proteins
pcmap is a PYTHON 3.X library designed to compute pairwise amino acid contacts and residues Solvant Accessible Surface Area in protein stuctures. Structures must be provided as PDB coordinates or MDAnalysis trajectory files. Contacts are computed inside a single PDB file or across two PDB files structures. The library can compute one to thousands sets of contacts. Results are produced in JSON format and contacts are encoded in a simple dictionary structure described in the **OUTPUT** section.

## Installation
`pip install pcmap`

or install the package in edition mode, if required

```sh
git clone pcmap
pip install -e pcmap
```

### Dependencies

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

Just pass the name of the PDB file to the `single` command.

`python -m pcmap single data/1A2K_r_u.pdb`

##### Single two-body contact map

Just pass the name of two PDB files to the `dimer` command.

`python -m pcmap single data/1A2K_r_u.pdb`

#### Many one-body contact maps

Create a file containing the list of protein as a text file with one PDB file per line:

#### **`sample.lst`**
```txt
data/1A2K_r_u.pdb
data/1A2K_l_u.pdb
```

And pass it to the `many` command
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

##### --atomic [default=False]

If True, all atomic contacts are reported.

##### --rich [default=False]

If True, add cartesian coordinates to contact map, only compatible with one single body atomic computation

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

which will return (only a sample is shown here):

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

## Computing Solvant Accessible Surface Area
The pcmap package implements a point-count based technique for solvant accessible surface area computation.
Exclusion solvant sphere surface are discretized as Fibonacci grid based on Rodrigo Azevedo Moreira da Silva implementation and surface calculation follows A.Gonzalez [method](https://arxiv.org/pdf/0912.4540.pdf).

A single function can process multiple PDB files or snapshots along a trajectory
The first positional parameter is the source of one of two possible types: a list of paths to PDB files or a MDAnalysis Universe object. The API is designed to easily process All Atoms of CG structures.

##### Options are:

* `npos:Int`, the total number of elements to process
* `step: Int`, the increment between two consecutive processed elements
* `selector:str`, A simple segid selector: a string of the form "segid A or segid B ..."
* `probe:float`, the radius of the water probe in Ang.
* `hres:Boolean`, if set to True increase the Fibonacci grid resolution, trading accuracy for speed
* `vdw_map:Dict`, a map a VDW radii, not defined it will be guessed from the Universe object. For pdb inputs a default all atom map is used.
* `martini3:Boolean`, if set to True the default martini3 VDW radii map is used, overriding the default all atom radii map (usefull for martini3 PDB records).
##### Multithreading options are:

* `chunk_sz:Int` the number of elements to process per task
* `ncpu:Int`, the number of workers

```python
import MDAnalysis as md
from pcmap.sasa import compute_many
```

### Working with MDAnalysis trajectory

```python
ucg = md.Universe('../LSB_data/SPC_L11_CG_2.tpr', '../LSB_data/SPC_L11_CG_1_whole_skip100.xtc')

sasa_res = compute_many(ucg,npos=50,\
        chunk_sz=5, selector="segid seg_0_A", probe=1.91, hres=True)
```

The returned object is a SASA_Results with the following properties:

* `resname:np.array`
* `resID:np.array`
* `chainID:np.array`
* `sasa:np.2Darray`
* `nframe:Int`
* `nresidue:Int`

##### Compute the average and sigma SASA of each residue along the trajectory with `SASA_Results.stats()`
```python
sasa_res.stats()
#Displays
#[(('MET', '1', 'seg_0_A'), 180.49615112304687, 21.59135138757801),
# (('LEU', '2', 'seg_0_A'), 71.09777526855468, 31.050710406745228),
# (('SER', '3', 'seg_0_A'), 70.42540771484374, 8.279207237924606),
# (('LEU', '4', 'seg_0_A'), 106.15115112304687, 9.595877493939815),
# (('ASP', '5', 'seg_0_A'), 32.37488525390625, 15.482186278685225),
#  ...
# ]
```

### Working with PDB files

#### A single PDB file
Results are returned as a simple Dict.
```python
from pcmap.sasa import compute_from_pdb
pdb_path = "./Single_conf/AA.pdb"
n = 5

sasa_res = compute_from_pdb(pdb_path, selector="segid A", probe=1.4)
sasa_res
#Displays
#{'freeASA': [{'resname': 'MET',
#   'resID': '0 ',
#   'chainID': 'A',
#   'SASA': 168.6180419921875,
#   'frac': 0.8308700323104858},
#  {'resname': 'LEU',
#   'resID': '1 ',
#   'chainID': 'A',
#   'SASA': 136.68359375,
#   'frac': 0.8641166090965271},
# ...
# ]}
```

#### Several PDB files

Results are returned as a SASA_Results objects, with a `raw` attributes which stores a simple list of Dict.
```python
from pcmap.sasa import compute_many
pdb_path = "./Single_conf/AA.pdb"
n = 100

sasa_res = compute_many([ pdb_path for _ in range(0, n) ],npos=50,\
           chunk_sz=5, selector="segid A", probe=1.4, hres=True)
sasa_res.raw
#Displays
# [
# [{'freeASA': [{'resname': 'MET',
#     'resID': '0 ',
#     'chainID': 'A',
#     'SASA': 167.74444580078125,
#     'frac': 0.8317462801933289},
#    {'resname': 'LEU',
#     'resID': '1 ',
#     'chainID': 'A',
# ...
# }],
# ...
#]
```