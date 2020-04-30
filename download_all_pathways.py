import requests
import csv
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dest_folder')
args = parser.parse_args()

# Get list of all pathways
pathway_list_url = "https://signor.uniroma2.it/getPathwayData.php?description"
response = requests.get(pathway_list_url)
if response.status_code == 200:
    results = response.content.decode('utf-8').splitlines()
    reader = csv.reader(results, delimiter="\t")
    next(reader)  # skip over headers
    pathway_list = set([r[0] for r in reader])
    pathway_list.remove(' ')  # This is caused by weird separator: "^M    Daniela Posca"
    print(pathway_list)
    # {'SIGNOR-NOTCH_Myogenesis', 'SIGNOR-P38', 'SIGNOR-TA', 'SIGNOR-SAPK-JNK', 'SIGNOR-P38_Myogenesis',
    # 'SIGNOR-M1M2', 'SIGNOR-NOTCH', 'SIGNOR-RMS', 'SIGNOR-Autophagy', 'SIGNOR-PDAP', 'SIGNOR-WNT', 'SIGNOR-IL1R',
    # 'SIGNOR-GCR', 'SIGNOR-Hedgehog', 'SIGNOR-PC', 'SIGNOR-HPP', 'SIGNOR-GBM', 'SIGNOR-DR', 'SIGNOR-AD',
    # 'SIGNOR-MM', 'SIGNOR-WNT_Myogenesis', 'SIGNOR-AC', 'SIGNOR-HT', 'SIGNOR-AMPK', 'SIGNOR-MS', 'SIGNOR-TLR',
    # 'SIGNOR-CRC', 'SIGNOR-AML-MiniPathway', 'SIGNOR-IOA', 'SIGNOR-NFKBC', 'SIGNOR-NS', 'SIGNOR-AML',
    # 'SIGNOR-FapINS', 'SIGNOR-IS', 'SIGNOR-MCAPO', 'SIGNOR-G1-S_trans', 'SIGNOR-PD', 'SIGNOR-TCA', 'SIGNOR-LBC',
    # 'SIGNOR-Myogenesis', 'SIGNOR-G2-M_trans', 'SIGNOR-TGFb', 'SIGNOR-TC', 'SIGNOR-NSCLCN', 'SIGNOR-INSR',
    # 'SIGNOR-EGF', 'SIGNOR-IL6', 'SIGNOR-PI3K-AKT', 'SIGNOR-NFKBNC', 'SIGNOR-FSGS'}
    print(len(pathway_list), "pathways to download")
    print("Writing to", args.dest_folder)

    # Get data for each pathway - e.g. https://signor.uniroma2.it/getPathwayData.php?pathway=SIGNOR-MM&relations=only
    pathway_data_url = "https://signor.uniroma2.it/getPathwayData.php?pathway={}&relations=only"
    for pthwy_id in pathway_list:
        print("Getting", pthwy_id)
        pthwy_resp = requests.get(pathway_data_url.format(pthwy_id))
        if not os.path.exists(args.dest_folder):
            os.makedirs(args.dest_folder)
        dest_filename = f"{os.path.join(args.dest_folder, pthwy_id)}.tsv"
        with open(dest_filename, "w+") as out_f:
            out_f.write(pthwy_resp.content.decode('utf-8'))
    print("done")
else:
    print(response.status_code, "Something other than success occurred:\n", response.content)

