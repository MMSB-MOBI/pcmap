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

### Computing one-body contact map
The pairwise amino acid contact inside a molecule can be obtained with the following call.

#### Single one-body contact map

#### Many one-body contact maps



### Computing two-body contact map

#### Single two-body contact map

#### Straight contact map computation on structures

`python -m pcmap dimer data/1A2K_r_u.pdb data/1A2K_l_u.pdb`

##### Applying transformation prior to contact map computation

```
python -m pcmap dimer data/1A2K_r_u.pdb data/1A2K_l_u.pdb\
--euler=-1.961,2.066,-2.354 --trans=7.199,16.800,28.799\
--offA=-27.553,-8.229,-80.604 --offB=-67.006,0.11,-77.27
```

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

* **euler** reference a list of the <span>&alpha;, &beta;, &gamma;</span> Euler angles. Each triplet is transformation specific.
* **translation** reference a list of *x,y,z* components of one translation vector. It is transformation specific.
* **recOffset** is the translation vector centering the receptor to coordinates origin. It is common to all transformations.
* **ligOffset** is the translation vector centering the ligand to coordinates origin. It is common to all transformations.

A file example is joined as `data/1A2K_transformations_sample.json`
## PYTHON module

## OUTPUT