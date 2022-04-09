
Author: Rebecca Weiss

Overview:
get_gene_info.py
Takes two arguments: Host and Gene to extract gene expressed tissue. If no arguments are given, will default to
host "Homo sapiens" and gene "TGM1".

Example:
python3 get_gene_info.py
python3 get_gene_info.py -host "Homo sapiens" -gene AATK
python3 gene_information_query.py -host horse -gene API5

assignment5/config.py and assignment5/my_io.py are submodules of get_gene_info.py

Test scripts for the submodules can be found in:
tests/unit/test_my_io.py
tests/unit/test_config.py

pytest:
pytest --cov-report html --cov --cov-config=.coveragerc
