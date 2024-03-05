import numpy as np
import pandas as pd
import scipy.io as sio
import sys
import re
import os

pd.options.mode.chained_assignment = None  # default='warn'
# pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

def drop_time_duplicates(df, ms_thresh):
    # Consider any neighboring measurements within ms_thresh to be duplicates:
    df['t_rounded'] = df['t'].apply(lambda x: round(x * ms_thresh) / ms_thresh)
    df.drop_duplicates(subset='t_rounded', inplace=True)
    df.dropna(subset='t', inplace=True)
    df.drop(columns='t_rounded', inplace=True)
    # df['t_unix'] = pd.to_datetime(df['t_unix']*1000, unit='ms')
   
    return df

def process_mote_data(filename, freq, num_meas):
    # Remove lines with non-numeric data or wrong number of measurements, written in error to the SD card:
    rf = open(filename, 'r')
    lines = rf.readlines()
    temp_filename = filename[:-4] + '_numeric.txt'
    wf = open(temp_filename, 'w')
    for line in lines:
        if not re.search('[a-zA-Z]', line) and line.count(',') == (num_meas*3 + 3) or (num_meas*2 + 2):
            wf.write(line)
    rf.close()
    wf.close()

    df_fifo = pd.read_csv(temp_filename, header=None, on_bad_lines='skip', low_memory=False)
    os.system('rm ' + temp_filename.replace(' ', r'\ '))

    df_fifo = df_fifo.apply(pd.to_numeric, errors='coerce')

    # Check for erroneous periods with gaps (usually due to sensor failure):
    df_fifo['dt'] = (df_fifo[0]-df_fifo[0].shift()).fillna(0)
    old_leng = df_fifo.shape[0]
    # Exclude any which are not between 0.95*dt and 1.05*dt, where dt is the nominal sample period
    df_fifo = df_fifo[np.logical_and(df_fifo['dt'] < 1.05*num_meas/freq, df_fifo['dt'] > 0.95*num_meas/freq)]
    new_leng = df_fifo.shape[0]
    print('%d erroneously spaced times removed' %(old_leng - new_leng))

    data_type = 'pres'
    if len(df_fifo.columns)==5: # handles the ref mote case (2x pressure, 2x temp, not FIFO)
        df = df_fifo
        df.columns = ['t', 'Pa', 'Pb', 'Ta', 'Tb']
        data_type = 'ref'
    else:
        # num_meas = int((len(df_fifo.columns)-4)/3) # number of FIFO measurements per row

        df_cols = ['t', 'Pa', 'Pb', 'Pc', 'Ta', 'Tb', 'Tc']
        df = pd.DataFrame(columns=df_cols)

        for i in range(num_meas):
            curr_df = pd.concat([df_fifo[0] - (num_meas-1-i)/freq, df_fifo[1+i*3], df_fifo[2+i*3], df_fifo[3+i*3], df_fifo[1+num_meas*3], df_fifo[2+num_meas*3], df_fifo[3+num_meas*3]], axis = 1)
            curr_df.columns = df_cols
            df = pd.concat([df, curr_df])

    df.sort_values(by='t', inplace=True) # sort by time
    df.reset_index(inplace=True)
    df.drop('index', axis=1, inplace=True) # delete old index column

    # Fill missing temperatures:
    # if data_type != 'ref' and 'Ta' in set(df):
    #     df[['Ta', 'Tb', 'Tc']] = df[['Ta', 'Tb', 'Tc']].fillna(method='bfill')
    #     df[['Ta', 'Tb', 'Tc']] = df[['Ta', 'Tb', 'Tc']].fillna(method='ffill')

    # Filter out erroneous times
    df = df[np.logical_and(df['t'] > 1641030199, df['t'] < 1900000000)]

    # Remove erroneous measurements:
    df['Pa'].loc[np.logical_or(df['Pa'] < -7000, df['Pa'] > 7000)] = np.nan
    df['Pb'].loc[np.logical_or(df['Pb'] < -7000, df['Pb'] > 7000)] = np.nan
    if data_type == 'pres':
        df['Pc'].loc[np.logical_or(df['Pc'] < -7000, df['Pc'] > 7000)] = np.nan
        # Add 100000 Pa to measurements (this is subtracted on the Arduino for data compression):
        df[['Pa', 'Pb', 'Pc']] = df[['Pa', 'Pb', 'Pc']] + 100000
    elif data_type == 'ref':
        df[['Pa', 'Pb']] = df[['Pa', 'Pb']] + 100000

    # Exclude nonsensical times:
    df = df[df['t'] > 1577865600]  # should be later than 2020

    # Drop duplicates or NaT:
    df = drop_time_duplicates(df, 20)

    return df

def import_data(path, table_name, num_meas):
    df = process_mote_data(path + table_name + '.txt', 12.5, num_meas)
    sio.savemat(path + table_name + '.mat', {table_name:df.to_dict("list")})
    return True

if __name__ == "__main__":
    import_data(str(sys.argv[1]), str(sys.argv[2]), int(sys.argv[3]))
