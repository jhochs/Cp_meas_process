import numpy as np
import pandas as pd
import scipy.io as sio
import os
import sys
from ftplib import FTP
from io import BytesIO

FTP_host = '3.132.186.69'
FTP_user = 'datalog'
FTP_pass = 'windEng0!'

def read_mat(file):
    mat = sio.loadmat(file)
    keys = list(mat.keys())
    mdata = mat[keys[-1]]
    mdtype = mdata.dtype 
    dict_data = {n: mdata[n][0, 0].flatten() for n in mdtype.names}
    df = pd.DataFrame.from_dict(dict_data)

    return df

def import_data():
    # Connect to FTP server and get file listing:
    ftp = FTP(host=FTP_host)
    ftp.login(user=FTP_user, passwd=FTP_pass)
    ftp.cwd('650CalRoof')
    files = ftp.nlst()
    
    first = True
    for file in files:
        # Fetch data:
        flo = BytesIO()
        ftp.retrbinary('RETR ' + file, flo.write)
        flo.seek(0)
        df_cur = pd.read_csv(flo, header=None, skiprows=4)
        df_cur.rename(columns={0:'t', 2:'WDir', 3:'WS'}, inplace=True)
        df_cur = df_cur[['t', 'WS', 'WDir']]
        df_cur['t'] = pd.to_datetime(df_cur['t']).astype(int) / 10**9

        # Append:
        if first:
            df = df_cur
            first = False
        else:
            df = pd.concat([df, df_cur])
        # print('File %s read, %d lines appended; total df size = %d' %(file, df_cur.shape[0], df.shape[0]))

        # Move this file to _processed/ directory
        ftp.rename(fromname='/650CalRoof/' + file, toname='/650CalRoof_processed/' + file)

    return df

if __name__ == "__main__":
    filename = str(sys.argv[1])

    df = read_mat(filename)
    print('Original .mat had %d rows' %(df.shape[0]))
    df_new = import_data()

    # Combine:
    df = pd.concat([df, df_new])
    print('New .mat has %d rows' %(df.shape[0]))

    # Sort and remove duplicates:
    df.sort_values(by='t', inplace=True) # sort by time
    df.drop_duplicates(subset='t', inplace=True)
    df.reset_index(inplace=True, drop=True)

    # Export back to .mat:
    sio.savemat(filename, {'wind':df.to_dict('list')})
