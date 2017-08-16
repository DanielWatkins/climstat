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
    
