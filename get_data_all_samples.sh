#wrapper to run download for both sample UVM and SKCM studies


python3 library/get_data_all_samples.py /mnt/c/Users/caule/OncoMX/survival_dataset/TCGA_UVM_SKCM/gdc_sample_sheet.2020-11-13.tsv | tee  log/get_data_all_samples.log
