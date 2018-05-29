""" This script converts the SHEBA sounding files into a single CSV file. I have it set
so that only readings below 500 hPa are returned; that can be changed via the variable
'min_pressure'. Also required: set the location of the data to 'drive_loc' and then 
the name you want the output file to have as 'save_prefix'. Still not perfect - there are
a few files with corrupted data when I tried it, with things like $ instead of 0. I've been
able to dodge most of those via the 'find_row_range' function but a few still fail to load.
Soundings that fail to read in have the date printed, everything else is wrapped together and
saved. When a sounding fails to read, a second attempt is made, leaving out the widths argument.
This appears to work for now but is not necessarily very robust."""


import os
import pandas as pd
import numpy as np
from datetime import datetime

drive_loc = 
save_prefix = 
min_pressure = 500

def find_row_range(fname, min_pressure):
    """Scans file to identify beginning, end, and any rows with nonnumerics"""
    filter_char = lambda char: char.isalnum() or char in ['-', ' ', '.', '\n']

    with open(fname) as f:
        first = 0
        last = 0
        lines = f.readlines()
        nrows = len(lines)
        
        for idx in range(len(lines)):
            line = lines[idx]
            if line.split()[0] == '------':
                first = idx + 1
                skip_rows = list(np.arange(first))
            if (idx > first) & (first != 0):
                if float(line.split()[1]) < min_pressure:
                    last = int(idx)
                if not np.all([filter_char(x) for x in line]):
                    skip_rows.append(idx)
                    
        if last == 0:
            last = nrows
            
        if last < nrows:
            skip_rows += list(np.arange(last, nrows + 1))
            
    return skip_rows

                

import os
import pandas as pd
import numpy as np
from datetime import datetime

def find_row_range(fname, min_pressure):
    """Scans file to identify beginning, end, and any rows with nonnumerics"""
    filter_char = lambda char: char.isalnum() or char in ['-', ' ', '.', '\n']

    with open(fname) as f:
        first = 0
        last = 0
        lines = f.readlines()
        nrows = len(lines)
        
        for idx in range(len(lines)):
            line = lines[idx]
            if line.split()[0] == '------':
                first = idx + 1
                skip_rows = list(np.arange(first))
            if (idx > first) & (first != 0):
                if float(line.split()[1]) < min_pressure:
                    last = int(idx)
                if not np.all([filter_char(x) for x in line]):
                    skip_rows.append(idx)
                    
        if last == 0:
            last = nrows
            
        if last < nrows:
            skip_rows += list(np.arange(last, nrows + 1))
            
    return skip_rows

                

def read_sheba_data(fpath, min_pressure=0):
    """Reads text output from the SPARC high resolution
    radiosonde archive. Columns are renamed to be all lower-case
    and be full names rather than abbreviations, and if only data
    below a certain pressure level is desired, min_pressure level 
    can be set to subset the data.
    
    Returns a dataframe as well as a boolean quality-control flag designed
    to catch when the data file has a header but no data."""
    
    colnames = ['date', 'pressure', 'temperature', 'dewpoint', 
                'relative_humidity', 'u_wind', 'v_wind', 'speed',
                'direction', 'ascent_rate', 'longitude', 'latitude', 
                'range', 'azimuth_angle', 'height', 'qp', 'qt',
                'qrh', 'qu', 'qv', 'quv']

    na_values = {'time':9999.0, 'pressure':9999.0, 'temperature':999.0, 'dewpoint':999.0,
             'relative_humidity':9999.0, 'u_wind':999.0, 'v_wind':999.0, 'speed':999.0,
             'direction':999.0, 'ascent_rate':99.0, 'longitude':999.0, 'latitude':999.0,
             'range':999.0, 'azimuth_angle':999.0, 'height':99999.0, 'qp':99.0, 'qt':99.0,
             'qrh':99.0, 'qu':99.0, 'qv':99.0, 'qdz':99.0}
    
    dtype = {'pressure':float, 'temperature':float, 'dewpoint':float,
             'relative_humidity':float, 'u_wind':float, 'v_wind':float, 'speed':float,
             'direction':float, 'ascent_rate':float, 'longitude':float, 'latitude':float,
             'range':float, 'azimuth_angle':float, 'height':float, 'qp':float, 'qt':float,
             'qrh':float, 'qu':float, 'qv':float, 'qdz':float}
    
    widths = np.array([7, 7, 6, 6, 6, 7, 7,
                       6, 6, 6, 9, 8, 6, 6,
                       8, 5, 5, 5, 5, 5, 5])
    
    # Find line 16. Some files have repeated lines in the header
    with open(fpath) as f:
        for i, line in enumerate(f):
            if i == 4:
                date = pd.to_datetime(line[35:57])

    skip_rows = find_row_range(fpath, min_pressure)
    
    try:
        df = pd.read_fwf(fpath, na_values=na_values,
                         skiprows=skip_rows, dtype=dtype, 
                         widths=widths,
                         names=colnames)
        df = df[df['pressure'] >= min_pressure].copy()
        df.date = date
        return df
    except ValueError as e:
        print(date)
        print(e)
        try:
            df = pd.read_fwf(fpath, na_values=na_values,
                             skiprows=skip_rows, dtype=dtype,
                             names=colnames)
            df = df[df['pressure'] >= min_pressure].copy()
            df.date = date
            return df

            
        except ValueError:
            print('Also failed without widths column')
