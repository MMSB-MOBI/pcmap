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
 - name: LBMC
   index: 1   
 - name: MMSB
   index: 2

date: 09 February 2026
bibliography: paper_assets/pcmap.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Decades of research in Structural biologiy led to the accumulation of vast knowledge structural knowledge with [@burley_updated_2025] further expands with Deep Learing method[@jumper_highly_2021].
Application of Molecular modeling to develop drugs, study the molecular details of biological process.
Many of tehese iological process are drvie by specific 

A major contribtuion to protzein stability and aa tions are short and mid range moelcular intectio,n (electroratic vdw).

A convienrient represnration of rptein strucuter can be made by patiwirese distances between atoms.
The distances can be conbveniatnely be stored into pairwiese distance matrix.
Such matrces ahave proven ti be effective descitoporit of protein fold when applied to single proitein.
Statiublitsing secondary structure elemnt s corresponds to specific regions of the distnfe matrix.
effective desritproi of peotein-protein assoaction surface when distance are compute between atoms of seprated proteins.
Where critival associaiton betwene specif amino acid can be identified frepm the distznce matrix
When applied to molceualrt dynamics data, distancesmatrix are ripycially uised to identify relevant functional motions
In drug design, docking protocols aims at identiyfin the relevant interaction mode betwenn target and prey mocleuls.
Thousands of poisssible associaotnn mldoes needs to be ranked to identify the functironall ones, contact allowing fast still accurate charactreriont of the assoaicain modes. 
The pcmap packages provides fasta computation of amino acid pairwse in protein structure.


# Statement of need
In a first appraoch, the detection of relvenat disntace between pair of atoms in a structure would rerqire computation of all possible pariswise distance making this probkemen size quadratic with resprco to the totacal number of atoms.
The atomic coordiantes of the protein structure are projected into of 3D mesh. The dimensions of the meh cell are chosen such that the set of atomic pairwise distances effectively computed is reduced the atom population within a cell Mesh projection.



Some `stuff` [@launay_evaluation_2020]`
# State of the field     
Alternative efficient implementation of molecular distance matrix softxware exists{}, but they whether required the installation of third party software or not suited for the analysis of large batch of structure.
Was applied to previous study[@launay_evaluation_2020].

# Software design

Front-end Python library, parsing of inputs and the thread mnagnagemet
Mesh calculation are carreid out bu the CPython extension[@ccmap]
The APi is meant to be used on cli, on native python daat structure or on mdanalysis objects[@michaud-agrawal_mdanalysis_2011,@gowers_mdanalysis_2016]. These fucntions can be called from user Python code for production poruporses or indode jupyter notebook for protiptyping or data analysis.
Molecular modeling pipeline often involves the estiamtion of structure produced by various softwares. For this situaition, the pcmap package features executable program from the terminal.


# Usage and Performances
```python
from pcmap import contactMap
contactMap("PATH_TO/PDB_FILE_A", "PATH_TO/PDB_FILE_B")
```

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
   

![Caption for example figure.\label{fig:perf}](paper_assets/perf.png){ width=20% }

and referenced from text using \autoref{fig:perf}.

# AI usage disclosure

No generative AI tools were used in the development of this software, the writing
of this manuscript, or the preparation of supporting materials.

# Acknowledgements

We acknowledge contributions People and support from

# References
