# signor2gocam

[![Join the chat at https://gitter.im/geneontology/signor2gocam](https://badges.gitter.im/geneontology/signor2gocam.svg)](https://gitter.im/geneontology/signor2gocam?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Converting Signor pathways to GO-CAM

## Dependencies
[gocamgen](https://github.com/dustine32/gocamgen) and 
[ontobio](https://github.com/biolink/ontobio).
Also, complex info lookup file ("SIGNOR_complexes.csv") is currently hard-coded in and should be downloaded from [here](https://signor.uniroma2.it/downloads.php) before running. Will soon make this a changeable parameter.

## Running
Specify Signor pathway file (-f) as input and an outfile path with -o:
```
python3 pathway_importer.py -f "SIGNOR-G2-M_trans_02_03_18.tsv" -o "outfile.ttl"
```
