#### Unzipping the nested 6-second data files ####
# Data downloaded from SPARC ftp and saved to <data_path>/<data_folder>
# This loops through the data folders and unpacks the gzips, and 
# then removes them.

import os
import gzip
import shutil

data_path = 
data_folder = 
pattern = '*.dat.gz'

for root, dirs, files in os.walk(data_path + data_loc):
    for filename in shutil.fnmatch.filter(files, pattern):
        with gzip.open(os.path.join(root, filename), 'rb') as f_in:
            with open(os.path.join(root, filename)[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    ! rm {root + '/*.gz'}
    
#### Processing the data files ####
# Load the data files and save them as csv's.

# Information for making the filenames
data_path = 
read_loc =
save_loc = 
save_prefix = 'sparc-6s-'

stations = ['03020', '03160', '03190', '03198', '03937', '03940', 
            '03948', '03952', '03990', '04102', '04105', '04106', 
            '04830', '04833', '04837', '11501', '11641', '11813', 
            '11818', '12842', '12850', '12919', '12924', '13723', 
            '13841', '13880', '13889', '13897', '13957', '13985', 
            '13995', '13996', '14607', '14684', '14733', '14898', 
            '14918', '14929', '21504', '22010', '22536', '23023', 
            '23047', '23050', '23062', '23066', '23160', '23230', 
            '24011', '24023', '24061', '24127', '24131', '24225', 
            '24232', '25308', '25339', '25501', '25503', '25624', 
            '25713', '26409', '26411', '26510', '26615', '26616', 
            '26617', '27502','40308', '40309', '40504', '40505', 
            '40710', '41415', '53103', '53813', '53819', '53823', 
            '53829', '54762', '61705', '92803', '93734', '93768', 
            '93805', '94008', '94043', '94240', '94703', '94823', 
            '94980', '94982', '94983']

years = [str(x) for x in np.arange(1998,2012)]

# Functions for processing the data files
import pandas as pd
import numpy as np

# Functions for processing text files
def conv_date(fpath):
    """Simple wrapper to extract the date from the SPARC filename"""
    from datetime import datetime
    date = fpath.split('/')[-1].split('-')[1].split('.')[0]
    return datetime.strptime(date,'%Y%m%d%H')

def read_sparc_data(fpath, min_pressure=0):
    """Reads text output from the SPARC high resolution
    radiosonde archive. Columns are renamed to be all lower-case
    and be full names rather than abbreviations, and if only data
    below a certain pressure level is desired, min_pressure level 
    can be set to subset the data.
    
    Returns a dataframe as well as a boolean quality-control flag designed
    to catch when the data file has a header but no data."""
    import pandas as pd
    import numpy as np
    
    colnames = ['date', 'time', 'pressure', 'temperature', 'dewpoint', 
                'relative_humidity', 'u_wind', 
                'v_wind', 'speed', 'direction', 'ascent_rate', 
                'longitude', 'latitude', 'elevation', 
                'azimuth_angle', 'altitude', 'qp', 'qt', 'qrh', 'qu', 'qv', 'qdz']
    units = ['sec', 'mb', 'C', 'C', '%', 'm/s', 'm/s', 'm/s', 'deg', 
             'm/s', 'deg', 'deg', 'deg', 'deg', 'm', 'code', 'code', 
             'code', 'code', 'code', 'code']

    widths = np.array([19,6,6,5,5,5,6,6,5,5,5,8,7,5,5,7,4,4,4,4,4,4])+1
    na_values = {'time':9999, 'pressure':9999, 'temperature':999, 'dewpoint':999,
             'relative_humidity':9999, 'u_wind':9999, 'v_wind':9999, 'speed':999,
             'direction':999, 'ascent_rate':999, 'longitude':9999, 'latitude':999,
             'elevation':999, 'azimuth':999, 'atltitude':99999, 'qp':99, 'qt':99,
             'qrh':99, 'qu':99, 'qv':99, 'qdz':99}
    
    # Find line 16. Some files have repeated lines in the header
    with open(fpath) as f:
        skip = 0
        qc = False # flag to note if any data present
        for line in f:
            if int(line[18:20]) == 16:
                qc = True
                break
            else:
                skip += 1
    
    
        
    df = pd.read_fwf(fpath, widths=widths, na_values=na_values,
                     skiprows=skip, names=colnames)
    
    df = df[df['pressure'] >= min_pressure].copy()
    df.Date = conv_date(fpath)
    return df, qc

# For loop to go through all the stations and make single csv files
for station in stations:
    station_list = []
    for year in years:
        path = os.path.join(drive_loc, read_loc,year)
        if station in os.listdir(path):
            path = os.path.join(drive_loc, read_loc,year, station)
            for file in os.listdir(path)[1:]:
                df, qc = read_sparc_data(os.path.join(path, file))
                if qc == True:
                    station_list.append(df)
    station_df = pd.concat(station_list)
    station_df.to_csv(os.path.join(drive_loc, save_loc, save_prefix + station + '.csv'))
    print(station)
