#wrapper to run download for both sample UVM and SKCM studies


python3 ../library/get_data_all_samples.py /mnt/c/Users/caule/OncoMX/survival_dataset/normalized_read_counts/TCGA-PRAD/gdc_sample_sheet.2020-11-18.tsv | tee  ../logs/get_data_all_samples.log
