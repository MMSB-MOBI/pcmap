---
title: 'Pcmap: A Python package for to fast computation of contact map of proteins'
tags:
  - Python
  - bioinformatics
  - protein
  - structure
  - docking
  - modeling
authors:
  - name: Guillaume Launay
    orcid: 0000-0003-0177-8706
    corresponding: true
   # equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Cecile Hilpert
   # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Juliette Martin
    orcid: 000-0002-4787-0885
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
affiliations:
 - name: Laboratory of Biology and Modelling of the Cell, Ecole Normale Supérieure de Lyon, CNRS UMR 5239, Inserm U1293, University Claude Bernard Lyon 1, Lyon, France 
   index: 1   
 - name: Molecular Microbiology and Structural Biochemistry, UMR 5086 CNRS, University Claude Bernard Lyon 1, Lyon, France
   index: 2

date: 09 February 2026
bibliography: paper_assets/pcmap.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Decades of research in structural biology have led to the accumulation of a vast structural knowledge of molecular macromolecules, notably proteins[@burley_updated_2025]. This experimental knowledge served as training corpus for deep learning methods[@jumper_highly_2021], which predictions led to a massive expansion of the structural space of known proteins structures[@fleming_alphafold_2025].
Molecular modeling thrives on this knowledge to achieve critical applications such as the molecular mechanics of biological processes, the identification of molecular causes of diseases or the rational development of drugs.
At the molecular level, biological processes are indeed driven by specific short- and mid-range physical interactions[@dill_molecular_2010]. Drug interaction modes for example correspond to specific chemical interactions between the drug molecule and its target. Likewise, biological function are carried out by specifc modes of association between proteins.
A compact representation of biomolecular structure, like proteins, can be obtained by filtering only the pair of atoms in chemical interactions (ie. physical proximity).
Thereby, full molecular structures can be further compressed in so-called contact matrix, a sparse data structure registering only pair of atoms closer in space than a parameter distance treshold. 
Contact matrix are effective descriptiors of protein folds (when applied to single proteins), protein-protein interactions (when applied to multimeric protein structures or functional molecular motions (when applied to molecular dynamic simulation).
<!-- 
Critical associations between specific amino acids can be identified from the distance matrix.
Applied to molecular dynamics data, contact matrices are typically used to identify relevant functional motions.
-->

The identification of functionaly relevent interaction modes between molecules is paramount in molecuar modeling. In typical modeling pipeline, thousands of molecular complexes are generated and subsequent processing is required to identify the relevant inteaction modes. In such situations, encoding the molecular structures by their contact maps allow for their fast, yet accurate, processing. To answer this needs, the `pcmap` package provides fast computation of amino acid pairwise distances in protein structures.

# Statement of need
A naive approach to detect the relevant distances between pair of atoms in a structure would require the computation of all possible pairwise distances. This makes the problem size quadratic with respect to the total number of atoms in the system.
The `pcmap` package reduces this complexity by projecting atomic coordiantes of protein structure into three dimensional mesh. The parameters of the mesh are chosen such that the set of atomic pairwise distances effectively computed is limited to the atom population within a cell and the ones in direct contact. 

# State of the field     
Alternative efficient implementation of molecular distance matrix softxware exists[@mdanalysis_2016; @mdanalysis_2011; @abraham_gromacs_2015], but they either required the installation of third party software or not suited for the analysis of large batch of structures. The presented Python package is an alternative lightweight, self consistent and yet efficient method for contact map computation, previously applied to a protein-protein associations study[@launay_evaluation_2020].

# Software design

The `pcmap` package is a Python library, which parses molecular structures, perform mutli-thread distance compuations and produces the resulting contact map in JSON format.
To achieve speed performance, protein structures are projected onto a three dimensional mesh managed by the associated CPython extension[@ccmap].
The `pcmap` package can be used on native python data structures or on PDB protein coordinate files[@burley_updated_2025]. The API can be called from user Python code for production purposes or inside jupyter notebook for prototyping and data analysis.
Because, molecular modeling pipeline often involves the processing of structure produced by various softwares, the `pcmap` package features the executables `pcmap-monomer`, `pcmap-dimer` and `pcmap-many` that can be invoked from the terminal.

# Usage and Performances
The `pcmap` modules exposes the two following functions: `contactMap` to compute the contact map of straight protein coordinates and `contactMapThroughTransform` to first transform initial coordinates and then compute their contact map.
Positional parameters can either be single or list of paths to protein coordinate files in PDB format[@burley_updated_2025]. In the following examples, `c1` will store the internal contact map of a single protein while `c2` will store the contact maps of three pairs of structures
```python
from contactMap,contactMapThroughTransform import pcmap

c1 = pcmap.contactMap("data/1A2K_r_u.pdb")
c2 = pcmap.contactMap(
  ["structOne_A.pdb","structTwo_A.pdb", "structThree_A.pdb"],
  ["structOne_B.pdb","structTwo_B.pdb", "structThree_B.pdb"]
  )
```

A variery of inputs can be passed to these functions to control their behaviour, additional documentation can be found on the project page[@pcmap].

The computed contact map stores amino acid ranked according to their residue number and chain identifier in the PDB record. To ensure that contacts are registred only once they are declared with the residue of the lowest rank (aka root) and list of their partners (aka partners). The corresponding JSON format is the following:

```json
{"type": "contactList",
 "data": [
    {"root": {"resID": "69 ", "chainID": "A"},
      "partners": [{"resID": "76 ", "chainID": "A"}]},
    {"root": {"resID": "41 ", "chainID": "B"},
      "partners": [
        {"resID": "72 ", "chainID": "A"},
        {"resID": "73 ", "chainID": "A"}]
      }
  ]
}
```

The C implementation makes it possible for the underlying mesh manupulation functions to release Python Global Interpreter Lock. Hence, "actual" multithreading can be achieved and performances scale decently with the number of workers\autoref{fig:perf}.   

![For this benchmark, up to 50000 protein-protein poses were generated and processed for three coordinate sets of increasing number of atoms: 2000(pdb code: 1GL1) 3500(pdb code: 1F34) 10000(pdb code: 2VIS).\label{fig:perf}](paper_assets/perf.png){ width=20% }

# Conclusion
The `pcmap` Python package computes contact map of proteins was recently released for 3.9 to 3.14 Python under Linux or MacOS operating systems.
The speed performance of the mesh routines underlying the `pcmap` module makes it a promising plateform for the future implementation of additional molecular metrics based on the local enviroment of atoms, such as solvant accessiblity calculations. 

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

<!--
# Acknowledgements

We acknowledge contributions People and support from

-->
# References
