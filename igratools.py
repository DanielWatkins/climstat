def download_igra(station_list_loc,save_loc):
    """Wrapper for a urllib request to the IGRA-Lite dataset. It requires
    a path to a csv file containing station identifiers for the stations you
    want to download data from. The CSV should have a column 'IGRA2' with the full IGRA2 station codes. An example
    file is 'station_shortlist.csv', which has info for Arctic stations. save_loc is 
    a path to where you want the files to be saved."""
    
    import pandas as pd
    import urllib.request
    
    station_list = pd.read_csv('station_shortlist.csv')
    if saveloc[-1] != '/':
        saveloc = saveloc + '/'
        
    for fname in station_list['IGRA2']:
        data = urllib.request.urlretrieve('ftp://ftp.ncdc.noaa.gov/pub/data/igra/data/data-por/'
               + fname + '-data.txt.zip', saveloc + fname + '-data.txt.zip')
       
    
def get_ID_lines(fname):
    """Each text file has up to a year's worth of radiosonde
    profiles. The record analyzer only can handle a single
    sounding at a time. getIDLines scans through a file and 
    finds the line numbers for sounding headers. fname is the 
    path to the file to analyze. Works for IGRA and for U Wyo 
    soundings."""
    
    fid = open(fname)
    lines = fid.readlines() # returns a list containing all lines of text from fname
    station_id = lines[0].split()[0]
    
    id_loc = []
    for x in range(0,len(lines)):
        # The reason for the try/except thing is that if consists only of white space
        # or has a single entry then the return from the .split() will not be subsettable.
        # In that case, it doesn't matter, since it also won't be a header line
        try:
            if lines[x].split()[0] == station_id:
                id_loc.append(x)
        except IndexError:
            pass
        
    return id_loc
    
def parse_IGRA_header(header_line):
    """Extract information from the header to individual data records in
    the IGRA2 data. Returns a dictionary with date, number of levels, lat,
    and lon."""
    import numpy as np
    import pandas as pd
    
    HEADREC       = header_line[0]
    ID            = header_line[1:12]
    YEAR          = header_line[13:17]
    MONTH         = header_line[18:20]
    DAY           = header_line[21:23]
    HOUR          = header_line[24:26]
    RELTIME       = header_line[27:31]
    NUMLEV        = header_line[32:36]
    P_SRC         = header_line[37:45]
    NP_SRC        = header_line[46:54]
    LAT           = header_line[56:62]
    LON           = header_line[64:71]
    
    date = pd.to_datetime('-'.join([YEAR,MONTH,DAY,HOUR]))
    
    return {'date':date,'numlev':int(NUMLEV),'lat':int(LAT),'lon':int(LON)}

def parse_IGRA_body(begin_idx,lines):
    """Extracts data from single radiosonde launch. Begin_idx is the index of
    the data header supplied by getIDLines(). Returns a dataframe with columns
    'date','p','gph','t','rh'.
    Units: p in hPa, gph in m, t in K, rh in %""" 
    
    import pandas as pd
    import numpy as np
    
    header = parse_IGRA_header(lines[begin_idx])
    date = header['date']
    n = header['numlev']
    PRESS = np.zeros(n)
    GPH = np.zeros(n)
    TEMP = np.zeros(n)
    RH = np.zeros(n)
    
    for j in range(n):
        i = begin_idx + 1 + j
        PRESS[j] = np.float(lines[i][9:15])/100
        GPH[j] = np.int(lines[i][16:21])
        TEMP[j] = np.float(lines[i][22:27])/10+273.15
        RH[j] = np.float(lines[i][28:33])/10
        
    PRESS[PRESS < 0] = np.nan
    GPH[GPH < 0] = np.nan
    TEMP[TEMP < 0] = np.nan
    RH[RH < 0] = np.nan
        
    return pd.DataFrame({'date':np.repeat(date,n),'p':PRESS,'gph':GPH,'t':TEMP,'rh':RH}) 
    
def select_IGRA_date_range(begin_year, end_year, id_loc, lines):
    """Go through the full data file, return only the indices of 
    record headers in id_loc that are within a date range (years).
    begin_year and end_year must be integers. id_loc is the object
    returned by get_ID_loc() and lines a list of all the lines in
    the text file for a particular station. Example usage:
    fname = 'CAN00071043-data.txt'
    id_loc = get_ID_loc(fname)
    
    f = open(fname)
    lines = f.readlines()
    begin_idx, end_idx = select_IGRA_date_range(1990,2005,id_loc,lines)
    """
    
    begin_found = False
    end_found = False
    idx = 0
    
    while ((not begin_found) | (not end_found)) & (idx < len(id_loc)):
        try:
            date, num = parseHeader(lines[id_loc[idx]])
        except IndexError:
            print(len(lines), len(id_loc), idx)
        if not begin_found:
            if date.year >= begin_year:
                begin_idx = idx
                begin_found = True
        else:
            if date.year > end_year:
                end_idx = idx - 1
                end_found = True
                
        idx += 1
    return (begin_idx, end_idx)    
    
    
    
    
    
