# pcmap : A python module to compute contact map of proteins
BLABLABLA
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


```
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

Pass a file containing the list of protein to compute in the simplest format, one PDB file per line:

#### **`sample.lst`**

```txt
data/1A2K_r_u.pdb
data/1A2K_l_u.pdb
```

And pass it to the cli
`python -m pcmap dimer data/1A2K_r_u.pdb data/1A2K_l_u.pdb`

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

```
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


## PYTHON module

The librairy can be used as a Python module to be assemble in the pipeline/program of your choice.
First we need to import it.

```python
import pcmap
```

The **pcmap** modules exposes the two following functions:

- *contactMap*
- *contactMapThroughTransform*

### The **contactMap** API

```python
contactMap(proteinA, proteinB=None, **kwargs)
```
 
 The type of the positianl parameters controls the function behaviouf
 First parameter can be PDB file OR a list of PDB files. Second parameter can be PDB file OR a list of PDB files.
 
##### Provided with a PDB file as single parameters:

Compute the internal amino acid contact map of the structure

##### Provided with another PDB file as optional second parameters:

Compute the amino acid contactmap between the two structures

##### Provided with a list of PDB files as single parameter:

Compute the all internal amino acid contact map of all structures indvidually.

##### Provided with a list of PDB files as second parameter:

Compute the amino acid contactmap accross pair of structures at identical positions in the two lists.

#### Examples
```python
# This will compute an internal contact map
c2 = pcmap.contactMap("data/1A2K_r_u.pdb")
# This will compute contact across the two structures
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

Both functions share the same set of named parameters:

* *dist*: cotnact threshold distance[default 4.5]
* *encode*: integer contact encoding[default=False]
* *threadNum*: maximal number of allowed threads[default=8]
* *deserialize*: get results as python dictionary if True, string otherwise [default=True]


## OUTPUT