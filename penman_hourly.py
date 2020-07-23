#!/usr/bin/env python
import os
import shutil
import ftplib
import zipfile
import numpy as np
import pandas as pd

"""This script downloads hourly weather data for a given weather station and calculates the hourly reference crop 
    evapotranspiration according to the FAO 1998 http://www.fao.org/3/X0490E/x0490e00.htm
    
    main('01048'', ('2007-01-01', '2007-05-31'))
    
    known issues:
    - when using nighttime/daytime for estimation of rs/rso, all values at night are zero, so we use 0.4 as constant
    """


def check_availability(dates, stat_id):
    """Checks availability of required station data for calculation of reference evapotranspiration
    returns bool for total availability and string 'solar' or 'sun' dependent on what data is available"""
    # todo: use date to check temporal availability
    try:
        os.mkdir('data/availability')
    except FileExistsError:
        pass
    station_df = {}
    dwd_url = 'opendata.dwd.de'
    dir_hourly = '/climate_environment/CDC/observations_germany/climate/hourly/'
    ftp = ftplib.FTP(dwd_url)
    ftp.login()
    variables = ['air_temperature', 'wind', 'dew_point', 'pressure', 'solar', 'sun']
    radiation = False
    available = True
    for variable in variables:
        try:
            ftp.cwd(dir_hourly + f'{variable}/' + 'recent')
        except ftplib.error_perm:
            ftp.cwd(dir_hourly + f'{variable}/')
        filelst = ftp.nlst()
        filestr = 'Stundenwerte_Beschreibung_Stationen'
        filename = [filename for filename in filelst if filestr in filename][0]
        with open(f'data/availability/{filename}', 'wb') as f:
            ftp.retrbinary("RETR " + filename, f.write)
        fields = ['Stations_id', 'von_datum', 'bis_datum', 'Stationshoehe', 'geoBreite', 'geoLaenge']
        df = pd.read_csv(f'data/availability/{filename}', delim_whitespace=True, usecols=fields)
        df.drop([0], inplace=True)
        df.reset_index(drop=True, inplace=True)
        if variable == 'solar':
            if stat_id in df['Stations_id'].to_list():
                radiation = 'solar'
        elif variable == 'sun' and radiation != 'solar':
            if stat_id in df['Stations_id'].to_list():
                radiation = 'sun'
        elif stat_id not in df['Stations_id'].to_list():
            available = False
            break
        station_df.update({variable: df})
    ftp.quit()
    shutil.rmtree('data/availability')
    return available, radiation


def download_stationdata(stat_id, variable):
    def get_filename():
        filelst = ftp.nlst()
        try:
            fname = [filename for filename in filelst if stat_id in filename][0]
        except IndexError:
            print(f'Error: {variable} data not available for station with station ID {stat_id}')
            raise
        return fname

    dwd_url = 'opendata.dwd.de'
    dir_hourly = '/climate_environment/CDC/observations_germany/climate/hourly/'
    ftp = ftplib.FTP(dwd_url)
    ftp.login()
    try:
        ftp.cwd(dir_hourly + f'{variable}/' + 'recent')
        filename = get_filename()

        with open(f'data/{stat_id}/{variable}/' + filename, 'wb') as f:
            ftp.retrbinary("RETR " + filename, f.write)
        with zipfile.ZipFile(f'data/{stat_id}/{variable}/{filename}', 'r') as zip_ref:
            zip_ref.extractall(f'data/{stat_id}/{variable}/')
        os.remove(f'data/{stat_id}/{variable}/{filename}')
        ftp.cwd(dir_hourly + f'{variable}/' + 'historical')
        contents = ftp.nlst()
        filename = [filename for filename in contents if stat_id in filename][0]
        with open(f'data/{stat_id}/{variable}/' + filename, 'wb') as f:
            ftp.retrbinary("RETR " + filename, f.write)
            ftp.quit()
        with zipfile.ZipFile(f'data/{stat_id}/{variable}/{filename}', 'r') as zip_ref:
            zip_ref.extractall(f'data/{stat_id}/{variable}/')
        os.remove(f'data/{stat_id}/{variable}/{filename}')
    except ftplib.error_perm:
        try:
            ftp.cwd(dir_hourly + f'{variable}/')
            filename = get_filename()
            with open(f'data/{stat_id}/{variable}/' + filename, 'wb') as f:
                ftp.retrbinary("RETR " + filename, f.write)
            with zipfile.ZipFile(f'data/{stat_id}/{variable}/{filename}', 'r') as zip_ref:
                zip_ref.extractall(f'data/{stat_id}/{variable}/')
            os.remove(f'data/{stat_id}/{variable}/{filename}')
            ftp.quit()
        except ftplib.all_errors:
            print(f'{variable} data not available for station with station ID {stat_id}')


class StationData:
    def __init__(self, s_id, daterange):
        self.id = s_id
        self.date = daterange
        self.vars = ['air_temperature', 'wind', 'dew_point', 'pressure', 'solar', 'sun']
        self.df = pd.DataFrame({'datetime': pd.date_range(self.date[0], self.date[1], freq='H')})
        self.df = self.df.set_index('datetime')
        self.pos = 0, 0, 0  # lat [°], lon [°], height [m]

    def download_data(self):
        try:
            os.mkdir(f'data/{self.id}')
        except FileExistsError:
            pass

        for variable in self.vars:

            try:
                os.mkdir(f'data/{self.id}/{variable}')
            except FileExistsError:
                pass
            print(f'Downloading {variable} data...')
            download_stationdata(self.id, variable)

    def read_stationdata(self):

        # definition of DWD column keys and legend
        col_keys = {'wind': ['F'],  # mean wind speed [m/s]
                    'air_temperature': ['TT_TU', 'RF_TU'],  # mean air temperature [°C], relative humidity [%]
                    'dew_point': ['TT'],  # dew point [°C]
                    'pressure': ['P0'],  # air pressure at station height [hPa]
                    'solar': ['ZENIT', 'FG_LBERG'],  # solar zenith angle at mid of interval [°],
                    # hourly sum of solar incoming radiation [J cm-2]
                    'sun': ['SD_SO']}  # hourly sunshine duration [min]

        # add the data for each variable to a column in our station df
        for variable in self.vars:
            # get filenames of data directory
            filelist = []
            for (dirpath, dirnames, filenames) in os.walk(f'data/{self.id}/{variable}/'):
                filelist.extend(filenames)
                break
            data_files = [file for file in filelist if 'produkt' in file]

            # get position of weather station (latitude, longitude, station height)
            geodata = pd.read_csv(f'data/{self.id}/{variable}/Metadaten_Geographie_{self.id}.txt', sep=';',
                                  encoding='ISO-8859-1')
            self.pos = geodata['Geogr.Breite'].iloc[-1], geodata['Geogr.Laenge'].iloc[-1], \
                       geodata['Stationshoehe'].iloc[-1]
            print(f'Reading {variable} data...')

            # the standard case: DWD has one zipfile for historical and one zipfile for recent data
            if len(data_files) == 2:
                df_hist = pd.read_csv(f'data/{self.id}/{variable}/{data_files[0]}', sep=';', encoding='ISO-8859-1')

                # remove whitespaces from column names
                df_hist.columns = df_hist.columns.str.lstrip()

                # remove incosistency in datetime column so format matches
                if df_hist['MESS_DATUM'].dtype != 'int64':
                    df_hist['MESS_DATUM'] = df_hist['MESS_DATUM'].astype(str).str[:-3].astype(np.int64)
                df_hist.index = pd.to_datetime(df_hist['MESS_DATUM'], format='%Y%m%d%H')
                df_recent = pd.read_csv(f'data/{self.id}/{variable}/{data_files[1]}', sep=';', encoding='ISO-8859-1')

                # remove whitespaces from column names
                df_recent.columns = df_recent.columns.str.lstrip()
                df_recent.index = pd.to_datetime(df_recent['MESS_DATUM'], format='%Y%m%d%H')
                df_merged = df_hist.append(df_recent)
                df_merged = df_merged.loc[~df_merged.index.duplicated(keep='first')]
            # the other case: DWD has only one zipfile for historical and recent data
            else:
                df_merged = pd.read_csv(f'data/{self.id}/{variable}/{data_files[0]}', sep=';', encoding='ISO-8859-1')

                # remove whitespaces from column names
                df_merged.columns = df_merged.columns.str.lstrip()

                # when pandas parses the date as string, date contains minutes. In this case we have to change
                # format for matching the date and resample to hours
                if df_merged['MESS_DATUM'].dtype != 'int64':
                    df_merged.index = pd.to_datetime(df_merged['MESS_DATUM'], format='%Y%m%d%H:%M')
                    df_merged = df_merged.resample('1H').mean()
            # add column to object df for given date range
            for col_key in col_keys.get(variable):
                self.df[variable + col_key] = df_merged.loc[self.date[0]: self.date[1], col_key]
            # DWD sets sunshine duration at night time (e.g. 3AM) to NaN, we replace them with zero
            if variable == 'sun':
                self.df[variable + col_keys.get('sun')[0]] = self.df[variable + col_keys.get('sun')[0]].fillna(0)


def lat_to_rad(lat):
    """converts latitude in decimal degrees to radians"""
    return np.pi / 180 * lat


def earth_sun_dist(doy):
    """calculates the inverse relative distance Earth-Sun and the solar declination"""
    d_r = 1 + 0.033 * np.cos((2 * np.pi / 365) * doy)  # inverse relative distance Earth-Sun (Eq. 23)
    delta = 0.0409 * np.sin((2 * np.pi / 365) * doy - 1.39)  # solar declination (Eq.24)
    return d_r, delta


def calc_sunrise(lat, lon, doy, date):
    """calculates time for sunrise in sunset to estimate daytime/nighttime for soil heat flux"""
    suntime_diff = -0.171 * np.sin(0.0337 * doy + 0.465) - 0.1299 * np.sin(0.01787 * doy - 0.168)
    declination = 0.4095 * np.sin(0.016906 * (doy - 80.086))
    lat_rad = lat_to_rad(lat)
    time_diff = 12 * np.arccos((np.sin(-0.0145) - np.sin(lat_rad) * np.sin(declination)) /
                               (np.cos(lat_rad) * np.cos(declination))) / np.pi
    sunrise_local = 12 - time_diff - suntime_diff
    sunset_local = 12 + time_diff - suntime_diff
    sunrise_dec = sunrise_local - lon / 15 + 1
    sunset_dec = sunset_local - lon / 15 + 1

    # convert time in decimals to datetime
    sunrise_dec *= 60
    sunset_dec *= 60
    date = date.replace(hour=0)
    sunrise = date + pd.to_timedelta(sunrise_dec, unit='m')
    sunset = date + pd.to_timedelta(sunset_dec, unit='m')
    return sunrise, sunset


def penman_hourly(station):
    """This function calculates the hourly reference evapotranspiration
    The calculation is performed according to the FAO Penman-Monteith method
    all methods and equations can be looked up at http://www.fao.org/3/X0490E/x0490e08.htm"""
    ets = []
    df = station.df.copy()
    print('Calculating Evapotranspiration...')
    # let's define some constants
    g_sc = 0.0820  # solar constant [MJ m-2 min-1]
    sigma = 2.043 * 10 ** -10  # Stefan-Boltzman constant for hourly intervals [MJ m-2 hour-1]

    for row_index, row in df.iterrows():
        # first we check for missing station data, in this case we set et = np.nan and jump to next row,
        # this can be interpolated later
        if -999 in row.to_list():
            ets.append(np.nan)
            continue
        # read station data for point of time
        date = row_index
        t_hr = row['air_temperatureTT_TU']  # 2m air temperature [°C]
        rh = row['air_temperatureRF_TU']  # 2m relative humidity [%]
        p = row['pressureP0'] / 10  # air pressure at station height [kPa]
        u_2 = row['windF']  # windspeed [m/s]
        r_s = (row['solarFG_LBERG']) / 100  # incoming solar radiation [MJ/m^2]
        w = lat_to_rad((row['solarZENIT']))  # solar zenith angle at mid of interval [RAD]

        # start of et calculation
        # e°(T) saturation vapour pressure at the air temperature T [kPa] (Eq. 11)
        e_o = 0.6108 * np.exp((17.27 * t_hr) / (t_hr + 237.3))
        # average hourly actual vapour pressure [kPa] (Eq. 54)
        e_a = e_o * (rh / 100)

        # calculation of r_n: net radiation at the grass surface [MJ m-2 hour-1]
        # first, we need to calculate extraterrestrial radiation r_a for hourly or shorter periods (Eq. 28)
        # solar time angles at the beginning and end of the period (Eq. 29/30)
        w_1, w_2 = w - np.pi / 24, w + np.pi / 24
        j = date.dayofyear
        # calculate sunrise and sunset time here for later use
        sunr, suns = calc_sunrise(station.pos[0], station.pos[1], j, date)
        # inverse relative distance Earth-Sun and the solar declination, j = DOY
        d_r, delta = earth_sun_dist(j)
        # station latitude [RAD]
        phi = lat_to_rad(station.pos[0])
        # finally: extraterrestrial radiation (Eq. 28)
        r_a = (12 * 60 / np.pi) * g_sc * d_r * ((w_1 - w_2) *
                                                np.sin(phi) * np.sin(delta) +
                                                np.cos(phi) * np.cos(delta) *
                                                (np.sin(w_2) - np.sin(w_1)))
        # r_ns: net shortwave radiation [MJ/m²/h] with albedo of a = 0.23 [-]
        # for hypothetical grass reference crop (Eq. 38)
        r_ns = (1 - 0.23) * r_s

        # Clear-sky solar radiation with a.pos[2] = station height [m] (Eq. 37)
        r_so = (0.75 + 2 * 10 ** -5 * station.pos[2]) * r_a

        # r_s/r_so represents the cloud cover. During daylight we estimate that ratio from measured solar radiation,
        # for nighttime wie assume r_s/r_so = 0.4 for humid climate (see FAO advices for hourly time step)
        # for now this does not work, more testing needed, so: 0.4
        if sunr < date < suns:
            quot_rs = r_s / r_so
        else:
            quot_rs = 0.4
        quot_rs = 0.4
        # r_nl: Net longwave radiation (Eq. 39 with modified Stefan-Boltzman constant)
        r_nl = sigma * ((t_hr + 273.16) ** 4) * (0.34 - 0.14 * np.sqrt(e_a)) * (1.35 * quot_rs - 0.35)

        # Eventually, we calculate the net radiation (Rn) for the penman-monteith-eq
        if r_s == 0:
            r_n = 0
        else:
            r_n = r_ns - r_nl

        # gamma: soil heat flux density [MJ m-2 hour-1] (Eqs. 45 and 46)
        g_day = 0.1 * r_n
        g_night = 0.5 * r_n

        # d: saturation slope vapour pressure curve at Thr [kPa °C-1] (Eq. 13)
        d = 4098 * (0.6018 * np.exp(17.27 * t_hr)) / (t_hr + 237.3) ** 2

        # psychrometric constant [kPa °C-1] (Eq. 8)
        g = 0.665 * 10 ** -3 * p

        # we calculate sunrise and sunset and check if night or day in order to adjust soil heat flux gamma
        if sunr < date < suns:
            gamma = g_day
        else:
            gamma = g_night

        # now that we have all our parameters, we finally calculate crop reference evapotranspiration et
        et = (0.408 * d * (r_n - g) + gamma * (37 / (t_hr + 273)) * u_2 * (e_o * t_hr - e_a)) / (
                d + gamma * (1 + 0.34 * u_2))

        # negative values for reference evapotranspiration may occur, for practical purposes we set the et = 0
        if et < 0:
            et = 0
        ets.append(et)
    df['et'] = ets
    return df


def main(s_id, date):
    available = check_availability(None, s_id)
    print(available)
    if available[0] is False or available[1] != 'solar':
        print(f'Not enough data for station with station id {s_id}')
        return False
    wstation = StationData(s_id, date)
    # comment next line if you have already downloaded station data
    wstation.download_data()
    wstation.read_stationdata()
    df_et = penman_hourly(wstation)
    try:
        os.mkdir('output')
    except FileExistsError:
        pass
    df_et.to_csv(f'output/{s_id}.csv')

    return df_et


ide = '01684'
main(ide, ('2007-04-01', '2007-06-01'))