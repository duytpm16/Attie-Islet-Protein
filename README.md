# Attie-Islet-Protein  
  
QTL2 pipeline for the Attie Islet Protein project.  

### Required original files:  
 1. raw protein level text file: "DO_islet_proteomics_non_normalized.txt"  
 2. sample annotation text file: "attie_DO_sample_annot.txt"  
 3. sample ChrM_Y info csv file: "attie_sample_info_ChrM_Y.csv"  
 4. rds dataframe of marker map: "marker_grid_0.02cM_plus.rds"  
  
Be consistent with the ___prefix___ variable. Use the same ___prefix___ in the normalization.R file with other scripts.  
  
  
### 1. Filter and Normalization. 
Filter the raw data as needed. (Ex. remove data with _x_ number of _NAs_)  
Input missing values by PCA and normalize the data by comBat normalization.  
  
Modify the _normalization.R_ file to save 4 output rds files:  
 1. _prefix_samples_annotation.rds_   
 2. _prefix_filtered_raw.rds_       
 3. _prefix_normalized.rds_   
 4. _prefix_rZ_normalized.rds_ 
 
### 2. Gather Data for Scan1. 
 
### 3. Run Scan1. 

### 4. Find LOD Peaks and Allele Effects. 

### 5. Format Data for QTL Viewer

