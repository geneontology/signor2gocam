# signor2gocam

[![Join the chat at https://gitter.im/geneontology/signor2gocam](https://badges.gitter.im/geneontology/signor2gocam.svg)](https://gitter.im/geneontology/signor2gocam?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Converting SIGNOR pathways to GO-CAM

## Dependencies
[gocamgen](https://github.com/dustine32/gocamgen) and 
[ontobio](https://github.com/biolink/ontobio).

## Running
Specify Signor pathway file `-f` as input, model title `-t`, and an outfile path with `-o`:
```bash
python3 pathway_importer.py -f SIGNOR-G2-M_trans_02_03_18.tsv -t SIGNOR-G2-M_trans -o outfile.ttl
```
Run tests:
```bash
python3 test.py
```

## Download and convert all SIGNOR pathways
```bash
python3 download_all_pathways.py -d downloaded_data
./generate_all_models.sh downloaded_data models
```
