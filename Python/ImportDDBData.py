import numpy as np
import pandas as pd
import scipy.io as sio
import base64
import json
import os
import sys

pd.options.mode.chained_assignment = None  # default='warn'

def decode(s):
    try:
        return base64.b64decode(s).decode('utf-8')
    except ValueError as e:
        return e

def read_json(s):
    try:
        return json.loads(s)
    except Exception as e:
        pass

def drop_time_duplicates(df, ms_thresh):
    # Consider any neighboring measurements within ms_thresh to be duplicates:
    df['t_rounded'] = df['t'].apply(lambda x: round(x * ms_thresh) / ms_thresh)
    df.drop_duplicates(subset='t_rounded', inplace=True)
    df.dropna(subset='t', inplace=True)
    df.drop(columns='t_rounded', inplace=True)
    # df['t_unix'] = pd.to_datetime(df['t_unix']*1000, unit='ms')
   
    return df

def download_mote_data(table_name):
    # Requires DynamoDBtoCSV installed in ~/DynamoDBtoCSV/ : https://github.com/edasque/DynamoDBtoCSV
    os.system('/usr/local/bin/node ~/DynamoDBtoCSV/dynamoDBtoCSV.js -t ' + table_name + ' > ' + table_name + '.csv')

def process_mote_data(filename, freq, df_error):
    df_fifo = pd.read_csv(filename, on_bad_lines='skip')
    df_fifo.drop(df_fifo.tail(1).index,inplace=True) # drop last line
    df_fifo = df_fifo.apply(pd.to_numeric, errors='coerce')
    df_fifo.dropna(subset=['t'], inplace=True) # drop lines with NaN in time column

    df_fifo = pd.concat([df_fifo, df_error])

    data_type = 'pres'
    if 'WS' in df_fifo.columns:
        df = df_fifo
        data_type = 'wind'

    elif len(df_fifo.columns)==5: # handles the ref mote case (2x pressure, 2x temp, not FIFO)
        df = df_fifo
        data_type = 'ref'

    elif len(df_fifo.columns)==4: # handles the non-FIFO (no temperatures) case
        df = df_fifo

    elif len(df_fifo.columns)==7: # handles the non-FIFO (with temperatures) case
        df = df_fifo

    else:
        num_meas = int((len(df_fifo.columns)-4)/3) # number of FIFO measurements per row

        df_cols = ['t', 'Pa', 'Pb', 'Pc', 'Ta', 'Tb', 'Tc']
        df = pd.DataFrame(columns=df_cols)

        for i in range(num_meas):
            curr_df = pd.concat([df_fifo['t'] + i/freq, df_fifo['Pa' + str(i)], df_fifo['Pb' + str(i)], df_fifo['Pc' + str(i)], df_fifo['Ta'], df_fifo['Tb'], df_fifo['Tc']], axis = 1)
            curr_df.columns = df_cols
            df = pd.concat([df, curr_df])

    df.sort_values(by='t', inplace=True) # sort by time
    df.reset_index(inplace=True)
    df.drop('index', axis=1, inplace=True) # delete old index column

    # Fill missing temperatures:
    if data_type != 'ref' and 'Ta' in set(df):
        df[['Ta', 'Tb', 'Tc']] = df[['Ta', 'Tb', 'Tc']].fillna(method='bfill')
        df[['Ta', 'Tb', 'Tc']] = df[['Ta', 'Tb', 'Tc']].fillna(method='ffill')

    # Rest of processing depends on the data type:
    if data_type == 'pres':
        # Remove erroneous measurements:
        df['Pa'].loc[np.logical_or(df['Pa'] < -5000, df['Pa'] > 5000)] = np.nan
        df['Pb'].loc[np.logical_or(df['Pb'] < -5000, df['Pb'] > 5000)] = np.nan
        df['Pc'].loc[np.logical_or(df['Pc'] < -5000, df['Pc'] > 5000)] = np.nan
        # Add 100000 Pa to measurements (this is subtracted on the Arduino for data compression):
        df[['Pa', 'Pb', 'Pc']] = df[['Pa', 'Pb', 'Pc']] + 100000

        # Remove time duplicates (by threshold):
        df = drop_time_duplicates(df, 20)
    elif data_type == 'ref':
        df['Pa'].loc[np.logical_or(df['Pa'] < -5000, df['Pa'] > 5000)] = np.nan
        df['Pb'].loc[np.logical_or(df['Pb'] < -5000, df['Pb'] > 5000)] = np.nan
        df[['Pa', 'Pb']] = df[['Pa', 'Pb']] + 100000
        # Drop time duplicates:
        df.drop_duplicates(subset='t', inplace=True)
    elif data_type == 'wind':
        df.drop_duplicates(subset='t', inplace=True)
        # Remove and fill erroneous measurements:
        df['WS'].loc[df['WS'] > 50] = np.nan
        df[['WS', 'WDir']] = df[['WS', 'WDir']].fillna(method='bfill')
    
    # Exclude nonsensical times:
    df = df[df['t'] > 1577865600]  # should be later than 2020

    # Drop rows with NaT and convert UNIX time to datetime:
    df.dropna(subset=['t'], inplace=True) # drop lines with NaN in time column
    # df['t'] = pd.to_datetime(df['t']*1000, unit='ms')

    return df

def import_data(table_name):
    print('Downloading ' + table_name + '...')
    download_mote_data(table_name)
    # error_table = table_name[:-4] + 'error'
    # print('Downloading ' + table_name[:-4] + 'error...')
    # download_mote_data(error_table)
    # try:
    #     df_error = read_error_data(error_table + '.csv')
    # except:
    df_error = pd.DataFrame()

    # os.system('rm ' + error_table + '.csv')
    df = process_mote_data(table_name + '.csv', 12.5, df_error)

    sio.savemat(table_name + '.mat', {table_name:df.to_dict("list")})
    os.system('rm ' + table_name + '.csv')
    # os.system('rm ' + error_table + '.csv')
    return True

def read_error_data(table_name):
    raw = pd.read_csv(table_name, usecols=['payload'])
    raw['idx0'] = raw['payload'].str.find('base64OriginalPayload').astype(int)
    raw['idx1'] = raw['payload'].str.find('clientId').astype(int)
    raw.dropna(inplace=True)
    raw['encoded_message'] = raw.apply(lambda x: x['payload'][x['idx0']+24:x['idx1']-3], axis=1) 
    raw['message'] = raw['encoded_message'].apply(decode)
    raw['message'] = raw['message'].str.replace('nan', '-100000')
    raw['message'] = raw['message'].str.replace('ovf', '-100000')
    raw['json_message'] = raw['message'].apply(read_json)
    df = pd.json_normalize(raw['json_message'])
    
    return df

if __name__ == "__main__":
    import_data(str(sys.argv[1]))
