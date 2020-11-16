import requests
import json
import csv
import os 
import sys

# Set up a dictionary for collecting the informatiom form the sample sheet. 
collector = {}

# Take in the sample sheet file from the command line.
inputFile = sys.argv[1]

# Set up a list for...
# Does not appear to be used at all in this script. 
primary_site = []

# Set up a counter for reading each row of the sample sheet file and skip the header. 
rowcount=0

# Open the sample sheet file.
with open(inputFile) as csvfile:
    csvreader = csv.reader(csvfile, delimiter="\t", quotechar='"')
    for row in csvreader:
        rowcount +=1
        if rowcount ==1:
            continue
        else:
           #File ID File Name       Data Category   Data Type       Project ID      Case ID Sample ID       Sample Type
           #558c4654-1be9-4eb2-9322-af29c48035ad    026527d4-007c-4c4e-9bd5-855c44bbe7b0.htseq.counts.gz    Transcriptome Profiling Gene Expression Quantification  TCGA-GBM        TCGA-06-0681    TCGA-06-0681-11A      Solid Tissue Normal
            
            # Save the file ID, project ID, and sample type. 
            fil_id   = row[0]
            proj_id  = row[4]
            sam_type = row[7].replace(" ","-")

            # Update collector dictionary
            
            # Do the project id and sample type exist in the dictionary? If not, create it. 
            if not(proj_id in collector):
                collector[proj_id] = {}
            if not(sam_type in collector[proj_id]):
                collector[proj_id][sam_type] = []
           
            # Add the file Id to the collector dictionary.
            collector[proj_id][sam_type].append(fil_id)
            




# Download everything!
files_endpt = "https://api.gdc.cancer.gov/files"
data_endpt = "https://api.gdc.cancer.gov/data"

# Set the path to the download destination.
path0 = "/mnt/c/Users/caule/OncoMX/survival_dataset/TCGA_UVM_SKCM/"

# Set up the collector dictionary where the structure is project_id : sample_type : file_id.
for proj_id in collector:

    # Create paths to establish destination folders.
    path1 = os.path.join(path0,proj_id)
    if not(os.path.exists(path1)):
        os.mkdir(path1)
    for sam_type in collector[proj_id]:
        path2 = os.path.join(path1,sam_type)
        if not(os.path.exists(path2)):
            os.mkdir(path2)

        # Make a list of all the files.
        fil_list = collector[proj_id][sam_type]

        # Create a file list of only unique names.
        uniq_fil_list = list(set(fil_list))


        print(proj_id, sam_type, len(fil_list), len(uniq_fil_list))


        params = {"ids": uniq_fil_list}

        response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type":"application/json"})

        file_name = path2+"/results.tar.gz"

        with open(file_name, "wb") as output_file:
            output_file.write(response.content)
