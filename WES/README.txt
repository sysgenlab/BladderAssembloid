This code generates mutational concordance between primary tumors and organoid samples

Requirements
- python (v 2.7.13)
- pandas (v. 0.24.2)
- openpyxl (v 2.6.4)
- xlrd (v 1.2.0)

Installation instruction
- All python packages can be downloaded and installed using pip (https://pypi.org/project/pip/).
- e.g. pip install [package name].
- Installation of all of the required packages should take few minutes.

Code (python)
- "analyze_organoid_tumor_seq_WES.py" to calculate mutational concordances between primary tumors and matching organoid samples. This code also generates variant summary files.
- "annotation_utilities.py" parses gene symbols and correponding uniprot IDs.

Data
- Data folder includes Variant Call Format (VCF) file of primary tumors and organoid samples.