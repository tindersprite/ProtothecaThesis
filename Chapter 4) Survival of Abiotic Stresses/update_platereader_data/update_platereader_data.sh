#!/bin/bash

# Version	 1.0.0
# Author	 Chris Shave
# Date		 28-01-2021

help(){
echo 
echo '################################################################################' 
echo '################################################################################'
echo
echo SYNOPSIS
echo update_platereader_data.sh
echo
echo DESCRIPTION
echo
echo EXAMPLES
echo
echo '################################################################################'
echo '################################################################################'
echo
}

# FUTURE PLANS
## 

# find all the mapping info and match it to all the raw values from the plate reader
# Create a temporary file with 

# regenerate raw data

# Once raw data is combined with maps, create a normalised OD column

# Save normalised data

# Compile all the normalised data together 

# CODE

cd /Users/u1984259/Documents/Main_Project/Scripts/update_platereader_data


### Remove previous error message, if it exists
if [ -e error_message ]
then rm error_message
fi


### Clear the temporary working directory, if it exists
if [ -e tmp ]
then rm tmp/*
fi


### Make a temporary working directory, if it doesn't exist
if ! [ -e tmp ]
then mkdir tmp
fi


### Define columns that should exist in map_column_names.txt



### Combine data files from experiments, if they were interrupted
Rscript Subscripts/combine_fragmented_data.R

### Find all the mapping info
ls plate_reader_link/Plate_Reader_Data/Mapping_Data/ | grep xlsx | grep -v bak > tmp/cond_maps.txt


### Find all the raw data
ls plate_reader_link/Plate_Reader_Data/Raw_Data/ | cut -f1 -d. > tmp/raw_data_full.txt


### Create a list of files that excludes endpoint data
cat tmp/raw_data_full.txt | grep -v endpoint > tmp/raw_data_filtered.txt

if ! [ -e tmp/raw_data.txt ]
then rm tmp/raw_data.txt 
fi

ln -s raw_data_filtered.txt tmp/raw_data.txt 
# ln -s raw_data_full.txt tmp/raw_data.txt 

### Match raw data files to mapping files. Save as a .csv
echo "date,experiment,plate_number,raw_data_file,mapping_data_file" >  tmp/mapping_pairs.csv
awk -F_ '{print $4 "," $1 "," $3 "," $0 ",NA"}' tmp/raw_data.txt >> tmp/mapping_pairs.csv

Rscript Subscripts/pair_maps_and_data.R

### Test to see if there are any NAs left in the mapping_pairs file
### If NAs exist then that would indicate missing condition map or a broken file name

if [ -e error_message ]
    then cat error_message
    exit 1
fi 

### Use the mapping_pairs file to combine all raw data with the correct mapping file.
### Save all mapped data in a new directory
### Produce graphs that show the mapped data faceted by strain for quality checking

Rscript Subscripts/map_raw_data.R

### Iterate through the mapped data and normalise it
### Save normalised data in a new directory

 Rscript Subscripts/normalise_mapped_data.R

### Compile data from similar experiments into single files

Rscript Subscripts/compile_data.R

### Create an R file that can take the compiled platereader data and make graphs from it

