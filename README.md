# Calculating Cp statistics from raw data

## Python
Python is used to import the data from AWS DynamoDB (for wifi/cellular motes), an SD card (for disconnected motes), or an FTP server (for CR300 wind measurements). This can be done as follows:
```python
python3 ImportDDBData.py table_name  # note - this requires DynamoDBtoCSV to be installed

python3 ImportFTPData.py prev_data_to_append.mat

python3 ImportSDData.py /Volumes/SD_card/path/ table_name 30  # the final argument is num_meas, the FIFO number of measurements cached by the sensor (see Arduino code)
```

## MATLAB
The two scripts in this repository are:
- [Cp_calc_650Cal.m](Cp_calc_650Cal.m)
- [Cp_calc_SN_xlsinput.m](Cp_calc_SN_xlsinput.m)

To run these files, the directories must be updated to point to the .mat data (which can be generated as described in the above section).

The rest of the files are dependent functions.

## Retrieving ASOS data
Weather station data from a nearby airport is used to calculate Cp. This can be done from [this website](https://mesonet.agron.iastate.edu/request/download.phtml?network=CA_ASOS). 
