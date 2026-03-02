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
  - name: Cecile Hilpert
   # equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Juliette Martin
    orcid: 000-0002-4787-0885
    # corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1
  - name: Guillaume Launay
    orcid: 0000-0003-0177-8706
    corresponding: true
   # equal-contrib: true
    affiliation: "1" # (Multiple affiliations must be quoted)
  
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

Contact matrices are effective descriptors of protein folds (when applied to single proteins), protein-protein interactions (when applied to multimeric protein structures), and molecular motions (when applied to molecular dynamics simulations). Molecular modeling approaches aim at discovering the biologically relevant protein shape, interactions, and motions. The identification of functional protein interactions, for example, can typically require the generation of thousands of putative molecular complexes, which requires subsequent processing and ranking. In such situations, representing the  molecular structures by their contact maps allows for fast, yet accurate, computation. To address this need, the `pcmap` package provides fast computation of amino acid pairwise distances in protein structures.

# Statement of need

Decades of research in structural biology have led to the accumulation of vast structural knowledge about molecular macromolecules, notably proteins [@burley_updated_2025]. This experimental knowledge has served as a training corpus for deep learning methods [@jumper_highly_2021], whose predictions have led to a massive expansion of the structural space of known protein structures [@fleming_alphafold_2025].
Molecular modeling thrives on this knowledge to achieve critical applications, such as understanding the molecular mechanics of biological processes, identifying the molecular causes of diseases, and rationally developing drugs. At the molecular level, biological processes are driven by specific short- and mid-range physical interactions [@dill_molecular_2010]. Drug interaction modes, for example, correspond to specific chemical interactions between the drug molecule and its target. Similarly, biological functions are carried out by specific modes of association between proteins.
A compact representation of biomolecular structure, such as proteins, can be obtained by filtering only the pairs of atoms involved in chemical interactions (i.e., those in physical proximity). Full molecular structures can be further compressed into so-called contact matrices, a sparse data structure that registers only pairs of atoms closer in space than a parameterized distance threshold.

A naive approach to detecting relevant distances between pairs of atoms in a structure would require computing all possible pairwise distances, making the problem size quadratic with respect to the total number of atoms in the system. The `pcmap` package reduces this complexity by projecting atomic coordinates of protein structures onto a three-dimensional mesh. The mesh parameters are chosen such that the set of atomic pairwise distances effectively computed is limited to the atom population within a cell and those in direct contact.

# State of the field  

Alternative efficient implementations of molecular distance matrix software exist [@mdanalysis_2016; @mdanalysis_2011; @abraham_gromacs_2015], but they either require the installation of third-party software or are not suited for analyzing large batches of structures. The presented Python package is a lightweight, self-contained, and efficient alternative for contact map computation, previously applied to large scale structural studies of proteins[@tam_alphacutter_2023; @launay_evaluation_2020].

# Software design

The `pcmap` package is a Python library that parses molecular structures, performs multi-threaded distance computations, and produces the resulting contact map in JSON format. To achieve high performance, protein structures are projected onto a three-dimensional mesh managed by an associated CPython extension [@ccmap].
The `pcmap` package can be used with native Python data structures or on PDB protein coordinate files [@burley_updated_2025]. The API can be called from user Python code for production purposes or within Jupyter notebooks for prototyping and data analysis. Since molecular modeling pipelines often involve processing structures produced by various software, the `pcmap` package features the executables `pcmap-monomer`, `pcmap-dimer`, and `pcmap-many`, which can be invoked from the terminal.


# Usage and Performances

The `pcmap` module exposes two main functions: `contactMap` to compute the contact map of protein coordinates directly, and `contactMapThroughTransform` to first transform initial coordinates and then compute their contact map. Positional parameters can be either single or lists of paths to protein coordinate files in PDB format [@burley_updated_2025]. In the following examples, `c1` will store the internal contact map of a single protein, while `c2` will store the contact maps at the interface of three pairs of structures:
```python
from contactMap,contactMapThroughTransform import pcmap

c1 = pcmap.contactMap("path/to/structure.pdb")
c2 = pcmap.contactMap(
  ["structOne_A.pdb","structTwo_A.pdb", "structThree_A.pdb"],
  ["structOne_B.pdb","structTwo_B.pdb", "structThree_B.pdb"]
  )
```

A variety of inputs can be passed to these functions to control their behavior; additional documentation can be found on the project page [@pcmap].
The computed contact map stores amino acids ranked according to their residue number and chain identifier in the PDB record. To ensure that contacts are registered only once, they are declared with the residue of the lowest rank (the "root") and a list of its partners. The corresponding JSON format is as follows:

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

The C implementation allows the underlying mesh manipulation functions to release Python's Global Interpreter Lock (GIL). This enables "actual" multithreading, and performance scales well with the number of workers (see \autoref{fig:perf}).

![Benchmarking batches of contact map computations. Amino acid contact were computed at the interface for protein-protein complexes, with a threshold contact distance of 4.5\AA. From 100 to 50000 different protein-protein poses were generated with the ZDOCK software[@pierce_zdock_2014] and processed to compute their interface contact maps. The two initial  protein-protein complexes feature respective sizes of  1974 atoms (pdb code: 1GL1) and 10677 atoms (pdb code: 2VIS).\label{fig:perf}](paper_assets/perf.png){ width=100% }

# Conclusion
The `pcmap` Python package[@pcmap], which computes contact maps of proteins, was recently released for Python versions 3.9 to 3.14 on Linux or macOS operating systems. The high performance of the mesh routines underlying the pcmap module makes it a promising platform for future implementations of additional molecular metrics based on the local environment of atoms, such as solvent accessibility calculations.

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

<!--
# Acknowledgements

We acknowledge contributions People and support from

-->
# References
