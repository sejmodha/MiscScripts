# MiscScripts

### UpdateKrakenDatabases.py
This script requires following packages.
- Python 3.0 or higher
- BioPython
- Pandas

Output: The default name for the output database is HumanVirusBacteria.

### reorderHCMVContig.py
This script requires following packages.
- Python 3.0 or higher
- BioPython
- pyfaidx (https://github.com/mdshw5/pyfaidx)
- pymummer (https://github.com/sanger-pathogens/pymummer)
- mummer3 (http://mummer.sourceforge.net/)
```
usage: reorderHCMVContig.py [-h] -q QUERY -r REFERENCE -c COORDINATES -i
                            MINIDENTITY -o OUTPUT [-n HEADER]

This script runs mummer using python

optional arguments:
  -h, --help            show this help message and exit
  -q QUERY, --query QUERY
                        Input query genome fasta file name with one sequence
  -r REFERENCE, --reference REFERENCE
                        Input reference fasta file
  -c COORDINATES, --coordinates COORDINATES
                        Output coordinates file name
  -i MINIDENTITY, --minidentity MINIDENTITY
                        Minimum percent identity for filtering results
  -o OUTPUT, --output OUTPUT
                        Output file in fasta format
  -n HEADER, --header HEADER
                        Header name for output fasta format
```
