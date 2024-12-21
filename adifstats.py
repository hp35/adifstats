# -*- coding: utf-8 -*-
"""
For specifications of the ADIF file format, see https://www.adif.org/100/adif_100.htm

Makes use of the following packages:
    pyhamtools       https://github.com/dh1tw/pyhamtools 
    adif_io-0.0.3    https://pypi.org/project/adif-io/
    geopandas  https://geopandas.org/en/stable/gallery/plotting_with_geoplot.html
    geoplot    https://residentmario.github.io/geoplot/index.html

    See also https://github.com/tylert/maidenhead
"""

import adif_io, geopy.distance, math, sys
# import geoplot
import time, datetime
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as image
import numpy as np
import geopandas as gpd
from pyhamtools import LookupLib, Callinfo, locator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from matplotlib import cm
import matplotlib.colors as colors
from astral import LocationInfo
from astral.sun import sun
from pandas.plotting import register_matplotlib_converters
import statistics
import pycountry

register_matplotlib_converters()

class maidenhead:
    longitude_width_degrees=2.0
    latitude_height_degrees=1.0

"""
When extracting the country names from the ADIF file, quite a few of them do
not match the formal names needed for extraction of country codes according to
the ISO 3166-1 standard (see https://en.wikipedia.org/wiki/ISO_3166-1_alpha-2).
As these country codes are needed when searching for country flags (for example
when generating the country histogram of QSO:s), we here take the precaution to
add the correct names as alternatives to the ones as listed in the ADIF format.
Feel free to add alternatives along the way (the following ones are only the
ones I have encountered in my own log files so far, and the list is by no means
exhaustive).
"""
alternativeCountryNames = dict({
    'Czech Republic'          : 'Czechia',
    'Slovak Republic'         : 'Slovakia',
    'Asiatic Turkey'          : 'Turkey',
    'Asiatic Russia'          : 'Russian Federation',
    'European Russia'         : 'Russian Federation',
    'Hawaii'                  : 'United States',
    'Crete'                   : 'Greece',
    'England'                 : 'United Kingdom',
    'Scotland'                : 'United Kingdom',
    'Northern Ireland'        : 'United Kingdom',
    'Wales'                   : 'United Kingdom',
    'Sicily'                  : 'Italy',
    'Sardinia'                : 'Italy',
    'Canary Islands'          : 'Spain',
    'Balearic Islands'        : 'Spain',
    'Madeira Islands'         : 'Portugal',
    'Bosnia-Herzegovina'      : 'Bosnia and Herzegovina',
    'Fed. Rep. of Germany'    : 'Germany',
    'Dodecanese'              : 'Greece',
    'Azores'                  : 'Portugal',
    'Venezuela'               : 'Venezuela, Bolivarian Republic of',
    'European Turkey'         : 'Turkey',
    'Moldova'                 : 'Moldova, Republic of',
    'Svalbard'                : 'Norway',
    'Ceuta & Melilla'         : 'Spain',
    'Corsica'                 : 'France',
    'Republic of Korea'       : 'Korea, Republic of',
    'MARITIME MOBILE'         : None,
    'Taiwan'                  : 'Taiwan, Province of China', # Uh, really?
    'US Virgin Islands'       : 'United States',
    'Republic of Kosovo'      : 'Kosovo',
    'Aland Islands'           : 'Åland Islands',
    'Iran'                    : 'Iran, Islamic Republic of',
    'Trinidad & Tobago'       : 'Trinidad and Tobago',
    'Reunion Island'          : 'France',
    'Shetland Islands'        : 'United Kingdom',
    'Bonaire'                 : 'Bonaire, Sint Eustatius and Saba',
    'Alaska'                  : 'United States',
    'Unknown'                 : None,
    'Falkland Islands'        : 'United Kingdom',
    'West Malaysia'           : 'Malaysia',
    'Bolivia'                 : 'Bolivia, Plurinational State of',
    'South Cook Islands'      : 'Cook Islands',
    'St. Kitts & Nevis'       : 'Saint Kitts and Nevis',
    'African Italy'           : 'Italy',
    'UK Base Areas on Cyprus' : 'Cyprus',
    'Vietnam'                 : 'Viet Nam',
    'Antigua & Barbuda'       : 'Antigua and Barbuda',
    'St. Lucia'               : 'Saint Lucia',
    'Guantanamo Bay'          : 'Cuba',
    'East Malaysia'           : 'Malaysia',
    'Curacao'                 : 'Curaçao'
})

def maidenhead_distance(gridsq1, gridsq2):
    """
    Compute the distance between two locations expressed in the Maidenhead
    grid. https://en.wikipedia.org/wiki/Maidenhead_Locator_System

    Parameters
    ----------
    gridsq1 : str
        The first Maidenhead grid location, expressed as a four- or
        six-character string.
    gridsq2 : str
        The second Maidenhead grid location, expressed as a four- or
        six-character string.

    Returns
    -------
    float
        The shortest distance between the two grid locations, expressed in
        kilometres (km).
    """
    latlong1 = locator.locator_to_latlong(gridsq1)
    latlong2 = locator.locator_to_latlong(gridsq2)
    return geopy.distance.geodesic(latlong1, latlong2).km

def sunriseset(gridsq, timestamp):
    latitude, longitude = locator.locator_to_latlong(gridsq)
    loc = LocationInfo()
    loc.latitude = latitude
    loc.longitude = longitude
    loc.timezone = "UTC"
    date = timestamp.date()
    s = sun(loc.observer, date=date, tzinfo=loc.timezone)
    sunrisetime = s['sunrise']
    sunsettime = s['sunset']
    return sunrisetime, sunsettime

def qsodict(adif_filename):
    qsos, header = adif_io.read_from_file(adif_filename)
    num_qsos = len(qsos)

    """
    Begin with loading the corresponding ISO 3166-1 alpha-2 country codes
    from the pycountry package.
    """
    generateISO3166Table = True
    countrycode_ISO3166 = {}
    if generateISO3166Table:
        f = open("iso3166codes.txt", "w")
    for country in pycountry.countries:
        countrycode_ISO3166[country.name] = country.alpha_2
        if generateISO3166Table:
            f.write("ISO3166-1 code for country '%s' -> '%s'\n"%(country.name,
                    countrycode_ISO3166[country.name]))
    if generateISO3166Table:
        f.close()

    """
    Initiate lookup tables for country identification from call sign.
    """
    my_lookuplib = LookupLib(lookuptype="countryfile")
    cic = Callinfo(my_lookuplib)

    for k in range(num_qsos):
        """Extract fields from the list of individualk QSOs"""
        call = qsos[k]["CALL"]
        gridsq = qsos[k]["GRIDSQUARE"]
#        mode = qsos[k]["MODE"]
#        rst_sent = qsos[k]["RST_SENT"]
#        rst_rcvd = qsos[k]["RST_RCVD"]
#        qso_date = qsos[k]["QSO_DATE"]
#        time_on = qsos[k]["TIME_ON"]
#        qso_date_off = qsos[k]["QSO_DATE_OFF"]
#        time_off = qsos[k]["TIME_OFF"]
#        band = qsos[k]["BAND"]
        freq = float(qsos[k]["FREQ"])
#        station_callsign = qsos[k]["STATION_CALLSIGN"]
        my_gridsq = qsos[k]["MY_GRIDSQUARE"]
#        tx_pwr = qsos[k]["TX_PWR"]
#        comment = qsos[k]["COMMENT"]
#        operator = qsos[k]["OPERATOR"]

        qsos[k]['freq'] = "{:9.4f}".format(freq)
        if len(gridsq) > 0:
            try:
                qsos[k]['latlong'] = locator.locator_to_latlong(gridsq)
            except ValueError:
                print("Unknown Maidenhead grid square '%s' (Ignored in compilation of data)"%gridsq)
                qsos[k]['latlong'] = ""
        else:
            qsos[k]['latlong'] = ""

        """Compute the Maidenhead grid distance from your own QTH to the call"""
        if len(gridsq) > 0:
            distance = maidenhead_distance(my_gridsq, gridsq)
            distance_str = "{:5.1f}".format(distance)
        else:
            distance_str = ""
        qsos[k]['distance'] = distance_str

        """Extract information on country and continent from call sign"""
        try:
            callinfo = cic.get_all(call)
            country = callinfo["country"]
            continent = callinfo["continent"]
        except KeyError:
            print("Unknown key '%s' (Ignored in compilation of data)"%call)
            country = "Unknown"
            continent = "Unknown"

        qsos[k]['country'] = country
        qsos[k]['continent'] = continent

        """
        Extract the two-letter ISO3166-1 country code from the country name.
        Example: https://stackoverflow.com/questions/16253060/\
        how-to-convert-country-names-to-iso-3166-1-alpha-2-values-using-python
        """
        if country in countrycode_ISO3166.keys():
            countrycode = countrycode_ISO3166[country]
        else:
            if country in alternativeCountryNames:
                country = alternativeCountryNames[country]
            else:
                raise ValueError("Failed to find alternative name for '%s' in "
                     "database (unidentifiable ISO 3166-1 code). Please check "
                     "entries of enclosed file iso3166codes.txt for further "
                     "guidance."%country)

        if country == None:
            countrycode = None
        elif country == 'Kosovo':
            countrycode = 'XK'
        else:
            countrycode = countrycode_ISO3166[country]

        qsos[k]['countrycode_ISO3166'] = countrycode

    return qsos, header

def qsodensity(qsos, buckets="maidenhead", txmode="all", freqband="all"):
    numqsos = 0 # Returns the total number of QSO:s matching the specification of band, mode, etc.
    density = {}
    rst_tx_mean = {}
    rst_rx_mean = {}
    rst_num = {}
    distance = {}
    density = {}
    """
    Extract the densities as follows:
        1. Per Maidenhead grid square.
        2. Per country.
    """
    if buckets == "maidenhead":   # Maidenhead grid-wise sampling and analysis
        for k in range(len(qsos)):
            if freqband == qsos[k]["BAND"] or freqband == "all":
                if txmode == qsos[k]["MODE"] or txmode == "all":
                    gridsq = qsos[k]["GRIDSQUARE"]
                    my_gridsq = qsos[k]["MY_GRIDSQUARE"]
                    if len(gridsq) > 0:
                        numqsos += 1
                        distance[gridsq] = maidenhead_distance(my_gridsq, gridsq)
                        if gridsq in density:
                            density[gridsq] += 1
                            if len(qsos[k]["RST_SENT"]) > 0 and len(qsos[k]["RST_RCVD"]) > 0:
                                rst_tx_mean[gridsq] += int(qsos[k]["RST_SENT"])
                                rst_rx_mean[gridsq] += int(qsos[k]["RST_RCVD"])
                                rst_num[gridsq] += 1.0
                        else:
                            density[gridsq] = 1
                            if len(qsos[k]["RST_SENT"]) > 0 and len(qsos[k]["RST_RCVD"]) > 0:
                                rst_tx_mean[gridsq] = int(qsos[k]["RST_SENT"])
                                rst_rx_mean[gridsq] = int(qsos[k]["RST_RCVD"])
                                rst_num[gridsq] = 1.0
        for gridsq in rst_tx_mean:
            rst_tx_mean[gridsq] /= rst_num[gridsq]
            rst_rx_mean[gridsq] /= rst_num[gridsq]
    elif buckets == "countries":  # Country-wise sampling and analysis
        for k in range(len(qsos)):
            if freqband == qsos[k]["BAND"] or freqband == "all":
                if txmode == qsos[k]["MODE"] or txmode == "all":
                    country = qsos[k]["country"]
                    if len(country) > 0:
                        numqsos += 1
                        if country in density:
                            density[country] += 1
                            if len(qsos[k]["RST_SENT"]) > 0 and len(qsos[k]["RST_RCVD"]) > 0:
                                rst_tx_mean[country] += int(qsos[k]["RST_SENT"])
                                rst_rx_mean[country] += int(qsos[k]["RST_RCVD"])
                                rst_num[country] += 1.0
                        else:
                            density[country] = 1
                            if len(qsos[k]["RST_SENT"]) > 0 and len(qsos[k]["RST_RCVD"]) > 0:
                                rst_tx_mean[country] = int(qsos[k]["RST_SENT"])
                                rst_rx_mean[country] = int(qsos[k]["RST_RCVD"])
                                rst_num[country] = 1.0
        for country in rst_tx_mean:
            rst_tx_mean[country] /= rst_num[country]
            rst_rx_mean[country] /= rst_num[country]
    else:
        raise ValueError("Unknown bucket type '%s'"%buckets)
    return density, rst_tx_mean, rst_rx_mean, distance, numqsos

def gridsummary(qsos):
    density, rst_tx_mean, rst_rx_mean, distance, numqsos = qsodensity(qsos, buckets="maidenhead")
    for gridsq in density:
        if (gridsq in rst_tx_mean) and (gridsq in rst_rx_mean):
            print("%8s %4s  %7.2f dB  %7.2f dB  %8.1f km"%(gridsq, density[gridsq], rst_tx_mean[gridsq], rst_rx_mean[gridsq], distance[gridsq]))
    return

def summarize(adif_filename, summary_filename = None):
    """
    Summarize the contents of the supplied ADIF file in plain text to standard
    output.

    Parameters
    ----------
    adif_filename : str
        File name of the ADIF file to summarize.
    """
    qsos, header = qsodict(adif_filename)
    num_qsos = len(qsos)
    num_qsos_num_digits = math.ceil(math.log10(num_qsos))

    if summary_filename == None:
        print("   # Date     Time    Mode  Freq       Callsign  Grid  RST_TX RST_RX  Distance               Country Continent")
    else:
        outfile = open(summary_filename, "w")
        outfile.write("   # Date     Time    Mode  Freq       Callsign   Grid RST_TX RST_RX Distance                Country Continent\n")

    for k in range(num_qsos):
        """Extract fields from the list of individualk QSOs"""
        call = qsos[k]["CALL"]
        gridsq = qsos[k]["GRIDSQUARE"]
        mode = qsos[k]["MODE"]
        rst_sent = qsos[k]["RST_SENT"]
        rst_rcvd = qsos[k]["RST_RCVD"]
        qso_date = qsos[k]["QSO_DATE"]
        time_on = qsos[k]["TIME_ON"]
        qso_date_off = qsos[k]["QSO_DATE_OFF"]
        time_off = qsos[k]["TIME_OFF"]
        band = qsos[k]["BAND"]
        freq = qsos[k]["freq"]
        station_callsign = qsos[k]["STATION_CALLSIGN"]
        my_gridsq = qsos[k]["MY_GRIDSQUARE"]
    #    tx_pwr = qsos[k]["TX_PWR"]
        comment = qsos[k]["COMMENT"]
    #    operator = qsos[k]["OPERATOR"]
        freq = qsos[k]['freq']
        latlong = qsos[k]['latlong']
        distance = qsos[k]['distance']
        country = qsos[k]['country']
        continent = qsos[k]['continent']

        """Summarize all relevant data row-wise"""
        if summary_filename == None:
            print("%*d %8s %6s %4s %6s %11s %5s %5s %6s %10s %22s %9s"
                %(num_qsos_num_digits, k+1, qso_date, time_on, mode, freq,
                  call, gridsq, rst_sent, rst_rcvd, distance, country, continent))
        else:
            outfile.write("%*d %8s %6s %4s %6s %11s %5s %5s %6s %10s %22s %9s\n"
                %(num_qsos_num_digits, k+1, qso_date, time_on, mode, freq,
                  call, gridsq, rst_sent, rst_rcvd, distance, country, continent))

    if summary_filename != None:
        outfile.close()

    """
    geoplot.polyplot(world, figsize=(8, 4))

    # initialize an axis
    fig, ax = plt.subplots(figsize=(8,6))
    # plot map on axis
    countries = gpd.read_file(  
         gpd.datasets.get_path("naturalearth_lowres"))
    countries[countries["name"] == "Australia"].plot(color="lightgrey",
                                                     ax=ax)
    # parse dates for plot's title
    first_month = df["acq_date"].min().strftime("%b %Y")
    last_month = df["acq_date"].max().strftime("%b %Y")
    # plot points
    df.plot(x="longitude", y="latitude", kind="scatter", 
            c="brightness", colormap="YlOrRd", 
            title=f"Fires in Australia {first_month} to {last_month}", 
            ax=ax)
    # add grid
    plt.show()
    """
    return

def figsize_px(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width_px, height_px = bbox.width, bbox.height
    return width_px*fig.dpi, height_px*fig.dpi

def qsomap(qsos, txmode="FT8", freqband="all", sizes="qsodensity", colors="txlevel"):
    density, rst_tx_mean, rst_rx_mean, distance, numqsos = qsodensity(qsos,
                        buckets="maidenhead", txmode=txmode, freqband=freqband)
    num_gridsquares = len(density)

    lat = np.empty(num_gridsquares, dtype=None)
    long = np.empty(num_gridsquares, dtype=None)
    qsodens = np.empty(num_gridsquares, dtype=None)
    rst_tx_mean_gridsq = np.empty(num_gridsquares, dtype=None)
    rst_rx_mean_gridsq = np.empty(num_gridsquares, dtype=None)
    for index, gridsq in enumerate(density):
        """
        Extract latitude and longitude for the Maidenhead grid square,
        taking care to exclude possibly erroneous position codes. At each
        square, extract the following:
            1. Total number of QSO:s at the square.
            2. The mean level in dB which the responding station reported.
            3. The mean level in dB which we received from the counterpart.
        """
        latlong = locator.locator_to_latlong(gridsq)
        if len(latlong) > 0:
            lat[index] = latlong[0]
            long[index] = latlong[1]
            qsodens[index] = density[gridsq]
            rst_tx_mean_gridsq[index] = rst_tx_mean[gridsq]
            rst_rx_mean_gridsq[index] = rst_rx_mean[gridsq]
        else:
            lat[index] = None
            long[index] = None
            qsodens[index] = None
            rst_tx_mean_gridsq[index] = None
            rst_rx_mean_gridsq[index] = None
        # print("%d: %s, %s, %s"%(index, lat[index], long[index], qsodens[index]))

    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    fig, ax = plt.subplots(figsize=(16,10))

    """First plot the world map to get axis limits etc."""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    world.plot(color="lightgrey", ax=ax)

    """Compute a reasonably scaled mapping of marker sizes for scatter plot"""
    ax.grid(b=True, alpha=0.5)
    ax.set_xlim([-198.0, 198.0])
    ax.set_ylim([-80.0, 80.0])

    """
    Define the variable of interest to plot for the colorized data.
    """
#    cmap = cm.viridis
#    cmap = cm.cool
    cmap = cm.rainbow
    for index, gridsq in enumerate(density):
        if colors == "txlevel":
            cvariable = rst_rx_mean_gridsq
            cnorm = mpl.colors.Normalize(vmin=np.min(cvariable), vmax=np.max(cvariable))
        elif colors == "qsodensity":
            cvariable = qsodens
            # cnorm = mpl.colors.Normalize(vmin=np.min(cvariable), vmax=np.max(cvariable))
            cnorm = mpl.colors.LogNorm(vmin=np.min(cvariable), vmax=np.max(cvariable))
        else:
            raise ValueError("Unknown marker color option '%s'"%colors)
        markerColor = cmap(cnorm(cvariable[index]))

        if sizes == "qsodensity":
            if qsodens[index] < 5:
                spotScale = 0.5 + qsodens[index]/10.0
            else:
                spotScale = 1.0 + qsodens[index]/50.0
            markerWidth = maidenhead.longitude_width_degrees*spotScale
            markerHeight = maidenhead.latitude_height_degrees*spotScale
            marker = Ellipse((long[index],lat[index]),
                           width=markerWidth,
                           height=markerHeight, angle=0,
                           linewidth=0, color=markerColor, clip_on=True)
            patch = Rectangle((long[index]-maidenhead.longitude_width_degrees/2.0,
                                       lat[index]-maidenhead.latitude_height_degrees/2.0),
                                      width=maidenhead.longitude_width_degrees,
                                      height=maidenhead.latitude_height_degrees,
                                      angle=0, linewidth=0, clip_on=True, visible=False)
            ax.add_patch(marker)
            ax.add_patch(patch)
            marker.set_clip_path(patch)
        elif sizes == "uniform":
            markerWidth = maidenhead.longitude_width_degrees
            markerHeight = maidenhead.latitude_height_degrees
            marker = Rectangle((long[index],lat[index]),
                           width=markerWidth,
                           height=markerHeight, angle=0,
                           linewidth=0, color=markerColor, clip_on=True)
            ax.add_patch(marker)
        else:
            raise ValueError("Unknown marker size option '%s'"%sizes)

    if freqband == "all":
        ax.set_title("QSO map (%s, all bands, %d QSO:s)"%(txmode, numqsos))
    else:
        ax.set_title("QSO map (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))

    cbartitle = ""
    if colors == "txlevel":
        cbartitle = 'Rx level [dB]'
    elif colors == "rxlevel":
        cbartitle = 'Tx level [dB]'
    elif colors == "qsodensity":
        cbartitle = 'QSO:s per gridsquare'
    else:
        raise ValueError("Unknown marker color option '%s'"%colors)

    fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap), cax=cax,
                 orientation='vertical', label=cbartitle)
    plt.show()
    plt.savefig("qsomap_%s_%s_%s_%s.png"%(sizes, colors, txmode, freqband), dpi=300)
    return

def qsograph(qsos, graphmode="DistanceVsTx", txmode="FT8", freqband="all"):
    mode = []
    rst_sent = np.zeros(len(qsos))
    rst_rcvd = np.zeros(len(qsos))
    qso_date = []
    time_on = []
    qso_date_off = []
    time_off = []
    band = []
    freq = np.zeros(len(qsos))
    station_callsign = []
    my_gridsq = []
    freq = []
    latlong = []
    distance = np.zeros(len(qsos))
    country = []
    continent = []
    qso_time_to_sunrise_hrs = np.zeros(len(qsos))
    qso_time_to_sunset_hrs = np.zeros(len(qsos))
    qso_duration_sec = np.zeros(len(qsos))
    numqsos = 0
    for k, qso in enumerate(qsos):
        if freqband == qso["BAND"] or freqband == "all":
            if txmode == qso["MODE"] or txmode == "all":
                mode.append(qso["MODE"])
                if len(qso["RST_SENT"]) > 0:
                    rst_sent[k] = int(qso["RST_SENT"])
                else:
                    rst_sent[k] = math.nan
                if len(qso["RST_RCVD"]) > 0:
                    rst_rcvd[k] = int(qso["RST_RCVD"])
                else:
                    rst_rcvd[k] = math.nan
                qso_date.append(qso["QSO_DATE"])
                time_on.append(qso["TIME_ON"])
                qso_date_off.append(qso["QSO_DATE_OFF"])
                time_off.append(qso["TIME_OFF"])
                band.append(qso["BAND"])
                freq.append(qso["freq"])
                station_callsign.append(qso["STATION_CALLSIGN"])
                my_gridsq.append(qso["MY_GRIDSQUARE"])
                freq.append(qso['freq'])
                latlong.append(qso['latlong'])
                if len(qso['distance']) > 0:
                    distance[k] = float(qso['distance'])
                else:
                    distance[k] = math.nan
                country.append(qso['country'])
                continent.append(qso['continent'])
                qsodate = datetime.datetime.strptime(qso["QSO_DATE"]+qso["TIME_ON"]+"-UTC","%Y%m%d%H%M%S-%Z")
                qso_on_unixepoch = time.mktime(datetime.datetime.strptime(
                    qso["QSO_DATE"]+qso["TIME_ON"]+"-UTC", "%Y%m%d%H%M%S-%Z").timetuple())
                qso_off_unixepoch = time.mktime(datetime.datetime.strptime(
                    qso["QSO_DATE_OFF"]+qso["TIME_OFF"]+"-UTC", "%Y%m%d%H%M%S-%Z").timetuple())
                qso_duration_sec[k] = qso_off_unixepoch - qso_on_unixepoch
                sunrise, sunset = sunriseset(qso["MY_GRIDSQUARE"], qsodate)
                sunrise_unixepoch = time.mktime(sunrise.timetuple())
                sunset_unixepoch = time.mktime(sunset.timetuple())
                qso_time_to_sunrise_hrs[k] = (qso_on_unixepoch - sunrise_unixepoch)/3600.0
                qso_time_to_sunset_hrs[k] = (qso_on_unixepoch - sunset_unixepoch)/3600.0
                numqsos += 1

    fig, ax = plt.subplots(figsize=(10,8))
    markerStyle = '.'
    if graphmode == "DistanceVsTx":
        ax.plot(distance, rst_sent, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("RST(sent) [dB]")
        ax.set_title("Sent signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "DistanceVsRx":
        ax.plot(distance, rst_rcvd, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("RST(rcvd) [dB]")
        ax.set_title("Received signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "DistanceVsTxMinusRx":
        ax.plot(distance, rst_sent - rst_rcvd, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("RST(sent) - RST(rcvd) [dB]")
        ax.set_title("Sent - received signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "DistanceVsQSODuration":
        ax.plot(distance, qso_duration_sec, marker=markerStyle, linestyle='none')
        median_qso_duration = statistics.median(qso_duration_sec)
        ax.axhline(y=median_qso_duration, color = 'tab:red', label='Median duration of QSO:s = %1.1f s'%median_qso_duration)
        ax.legend(loc='upper right')
        ax.set_xlabel("Distance [km]")
        ax.set_ylabel("QSO duration [s]")
        ax.set_title("QSO duration (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "TimeToSunsetVsTx":
        ax.plot(qso_time_to_sunset_hrs, rst_sent, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Time to sunset [hrs]")
        ax.set_ylabel("RST(sent) [dB]")
        ax.set_title("Sent signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "TimeToSunsetVsRx":
        ax.plot(qso_time_to_sunset_hrs, rst_rcvd, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Time to sunset [hrs]")
        ax.set_ylabel("RST(rcvd) [dB]")
        ax.set_title("Received signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "TimeToSunsetVsTxMinusRx":
        ax.plot(qso_time_to_sunset_hrs, rst_sent - rst_rcvd, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Time to sunset [hrs]")
        ax.set_ylabel("RST(sent) - RST(rcvd) [dB]")
        ax.set_title("Sent - received signal strength (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    elif graphmode == "TimeToSunsetVsDistance":
        ax.plot(qso_time_to_sunset_hrs, distance, marker=markerStyle, linestyle='none')
        ax.set_xlabel("Time to sunset [hrs]")
        ax.set_ylabel("Distance [km]")
        ax.set_title("Distance of QSO (%s @ %s, %d QSO:s)"%(txmode, freqband, numqsos))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("qsograph_%s_%s_%s.png"%(graphmode, txmode, freqband), dpi=300)
    return

def qsohistogram(qsos, txmode="FT8", freqband="all"):
    total_num_qsos = 0
    num_qsos_on_date = {}
    for qso in qsos:
        if freqband == qso["BAND"] or freqband == "all":
            if txmode == qso["MODE"] or txmode == "all":
                if qso["QSO_DATE"] in num_qsos_on_date.keys():
                    num_qsos_on_date[qso["QSO_DATE"]] += 1
                else:
                    num_qsos_on_date[qso["QSO_DATE"]] = 1
                total_num_qsos += 1
    num_qso_days = len(num_qsos_on_date)
    sorted_dates = sorted(num_qsos_on_date.keys())
    date_first = datetime.datetime.strptime(sorted_dates[0], "%Y%m%d").date()
    date_last = datetime.datetime.strptime(sorted_dates[-1], "%Y%m%d").date()
    date_first_str = date_first.strftime("%d %B, %Y")
    date_last_str = date_last.strftime("%d %B, %Y")
    print("Log ranges from %s to %s (mode=%s, band=%s)."%(date_first_str, date_last_str, txmode, freqband))
    print("Number of days with at least one QSO (mode=%s, band=%s): %d"%(txmode, freqband, num_qso_days))
    print("Total number of QSOs (mode=%s, band=%s): %d"%(txmode, freqband, total_num_qsos))
    t = np.zeros(num_qso_days, dtype='datetime64[s]')
    y = np.zeros(num_qso_days, dtype='int')
    ys = np.zeros(num_qso_days, dtype='int')
    for k, date in enumerate(sorted(num_qsos_on_date.keys())):
        t[k] = datetime.datetime.strptime(date, "%Y%m%d").date()
        y[k] = num_qsos_on_date[date]
        if k == 0:
            ys[k] = y[k]
        else:
            ys[k] = ys[k-1] + y[k]

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    ax.bar(t,y)
    ax.xaxis_date()
    ax.set_xlabel("Date")
    ax.set_ylabel("Num QSO:s on date", color=color)
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel("Total number of QSO:s", color=color)
    ax2.plot(t, ys, color=color)
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("QSO:s over time (%s @ all bands, %d QSO:s)"%(txmode, total_num_qsos))
    else:
        ax.set_title("QSO:s over time (%s @ %s, %d QSO:s)"%(txmode, freqband, total_num_qsos))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("qsohistogram_time_%s_%s.png"%(txmode, freqband), dpi=300)

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    median_of_num_qsos_per_day = statistics.median(y)
    max_number_of_qsos_in_a_single_day = np.max(y)
    bins=range(0, max_number_of_qsos_in_a_single_day+1)
    ax.hist(y, bins=bins, rwidth=0.8, align = "left", density=False) # Set density=False to generate counts instead of probability
    ax.axvline(x=median_of_num_qsos_per_day, color = 'tab:red', label='Median = %d QSO:s per day'%median_of_num_qsos_per_day)
    ax.legend(loc='upper right')
    ax.set_xlabel("Number of QSO:s per day")
    ax.set_ylabel("Number of days")
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("Histogram over QSO:s per day (%s @ all bands, %d QSO:s)"%(txmode, total_num_qsos))
    else:
        ax.set_title("Histogram over QSO:s per day (%s @ %s, %d QSO:s)"%(txmode, freqband, total_num_qsos))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("qsohistogram_hist_%s_%s.png"%(txmode, freqband), dpi=300)

    return

def gridhistogram(qsos, txmode="FT8", freqband="all"):
    worked_grid_squares = {}
    num_gridsquares_on_date = {}
    num_new_gridsquares_on_date = {}
    total_num_qsos = 0  # Number of QSO:s with and associated grid square (differs from the total number of QSO:s if counting special signals etc.)
    for qso in qsos:
        if freqband == qso["BAND"] or freqband == "all":
            if txmode == qso["MODE"] or txmode == "all":

                if qso["GRIDSQUARE"] != "": # If the log entry contains a grid square at all
                    
                    if qso["GRIDSQUARE"] in worked_grid_squares:
                        worked_grid_squares[qso["GRIDSQUARE"]] += 1
                    else:
                        worked_grid_squares[qso["GRIDSQUARE"]] = 1
    
                        if qso["QSO_DATE"] in num_new_gridsquares_on_date:
                            num_new_gridsquares_on_date[qso["QSO_DATE"]] += 1
                        else:
                            num_new_gridsquares_on_date[qso["QSO_DATE"]] = 1
    
                    if qso["QSO_DATE"] in num_gridsquares_on_date:
                        num_gridsquares_on_date[qso["QSO_DATE"]] += 1
                    else:
                        num_gridsquares_on_date[qso["QSO_DATE"]] = 1

                    total_num_qsos += 1

    total_num_gridsquares = len(worked_grid_squares)
    sorted_dates = sorted(num_gridsquares_on_date.keys())
    date_first = datetime.datetime.strptime(sorted_dates[0], "%Y%m%d").date()
    date_last = datetime.datetime.strptime(sorted_dates[-1], "%Y%m%d").date()
    date_first_str = date_first.strftime("%d %B, %Y")
    date_last_str = date_last.strftime("%d %B, %Y")
    print("Worked a total of %d unique grid squares between %s and %s."%(total_num_gridsquares, date_first_str, date_last_str))
    print("Number of days with at least one new grid square (mode=%s, band=%s): %d (out of %d days with at least one QSO)."%(txmode, freqband, len(num_new_gridsquares_on_date), len(num_gridsquares_on_date)))
    print("Average number of QSO:s per grid square: %1.2f (mode=%s, band=%s)."%(total_num_qsos/total_num_gridsquares, txmode, freqband))

    t = np.zeros(len(num_new_gridsquares_on_date), dtype='datetime64[s]')
    y = np.zeros(len(num_new_gridsquares_on_date), dtype='int')
    ys = np.zeros(len(num_new_gridsquares_on_date), dtype='int')
    for k, date in enumerate(sorted(num_new_gridsquares_on_date.keys())):
        t[k] = datetime.datetime.strptime(date, "%Y%m%d").date()
        y[k] = num_new_gridsquares_on_date[date]
        if k == 0:
            ys[k] = y[k]
        else:
            ys[k] = ys[k-1] + y[k]

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    ax.bar(t,y, width=0.8, color=color)
    ax.xaxis_date()
    ax.set_xlabel("Date")
    ax.set_ylabel("Number of new grid squares on date", color=color)
    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel("Total number of grid squares", color=color)
    ax2.plot(t, ys, color=color)
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("Number of worked grid squares over time (%s @ all bands, %d QSO:s)"%(txmode, total_num_qsos))
    else:
        ax.set_title("Number of worked grid squares over time (%s @ %s, %d QSO:s)"%(txmode, freqband, total_num_qsos))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("gridhistogram_time_%s_%s.png"%(txmode, freqband), dpi=300)

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    median_of_num_new_gridsquares_per_day = statistics.median(y)
    max_number_of_new_gridsquares_in_a_single_day = np.max(y)
    print(max_number_of_new_gridsquares_in_a_single_day)
    bins=range(0, max_number_of_new_gridsquares_in_a_single_day+1)
    ax.hist(y, bins=bins, rwidth=0.8, align = "left", density=False) # Set density=False to generate counts instead of probability
    ax.axvline(x=median_of_num_new_gridsquares_per_day, color = 'tab:red', label='Median = %d new grid squares per day'%median_of_num_new_gridsquares_per_day)
    ax.legend(loc='upper right')
    ax.set_xlabel("Number of new grid squares per day")
    ax.set_ylabel("Number of days")
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("Histogram over new grid squares per day (%s @ all bands, %d QSO:s)"%(txmode, total_num_qsos))
    else:
        ax.set_title("Histogram over new grid squares per day (%s @ %s, %d QSO:s)"%(txmode, freqband, total_num_qsos))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("gridhistogram_hist_%s_%s.png"%(txmode, freqband), dpi=300)

    return

def countryhistogram(qsos, txmode="FT8", freqband="all", numCountries=30):
    debugFlagLoading = False
    totNumCountries = 0
    totNumQSOs = 0
    worked_countries_all = {}  # Counting the number of QSO:s country-wise
    worked_distance = {}   # Counting the accumulated distance country-wise
    iso3166code = {}
    for qso in qsos:
        if freqband == qso["BAND"] or freqband == "all":
            if txmode == qso["MODE"] or txmode == "all":
                totNumQSOs += 1
                country = qso["country"]

                """Build-up of table (dictionary) of number of QSO:s"""
                if (country not in worked_countries_all) and (country != "Unknown"):
                    worked_countries_all[country] = 1
                    totNumCountries += 1
                    iso3166code[country] = qso["countrycode_ISO3166"]
                    if debugFlagLoading:
                        print("Got countrycode '%s' for new country '%s'"
                              %(iso3166code[country],country))
                elif (country != "Unknown"):
                    worked_countries_all[country] += 1

                """Build-up of table (dictionary) of accumulated distances"""
                if (country not in worked_distance) and (country != "Unknown"):
                    if qso["distance"] != "": # If the log entry contains a nonzero distance at all
                        worked_distance[country] = float(qso["distance"])
#                    iso3166code[country] = qso["countrycode_ISO3166"]
                elif (country != "Unknown"):
                    if qso["distance"] != "": # If the log entry contains a nonzero distance at all
                        worked_distance[country] += float(qso["distance"])

    """
    PART I - HISTOGRAM OF NUMBER OF QSO:S COUNTRY-WISE
    
    Sort the dictionary of worked countries in descending order and keep the
    highest numCountries ones. When creating the sorted dictionary, we rely
    on that dictionaries preserve insertion order since Python 3.7+.
    """
    worked_countries_all = dict(sorted(
        worked_countries_all.items(), key=lambda item: item[1], reverse=True))
    worked_countries = {k: worked_countries_all[k]
                        for k in list(worked_countries_all)[:numCountries]}
    countryNumbers = range(len(worked_countries))  # Indexing along x-axis
    countryNames = worked_countries.keys()
    countryQSOs = worked_countries.values()

    fig = plt.figure("Top-%d countries for number of QSO:s (%s @ %s)"
             %(numCountries, txmode, freqband),figsize=(16,12))
    ax = fig.add_subplot(111)
    plt.xticks(countryNumbers, countryNames, rotation=45, ha='right')
    barplot = ax.bar(countryNumbers, countryQSOs, color='dodgerblue')
    for idx, countryname in enumerate(worked_countries):
        numQSOs = worked_countries[countryname]
        countrycode = iso3166code[countryname]
        flagfilename = 'flags/region-flags-gh-pages/png/%s.png'%countrycode
        if debugFlagLoading:
            print("Loading flag for '%s' -> '%s' (%d QSO:s) from file %s"
                  %(countryname, countrycode, numQSOs, flagfilename))
        x = barplot[idx].get_x()
        barwidth = barplot[idx].get_width()
        barheight = barplot[idx].get_height()
        flagimg = image.imread(flagfilename)
        flagwidth_px = flagimg.shape[1]  # Num columns = width of array in px
        zoomval = 750.0*barwidth/(flagwidth_px*numCountries)
        if debugFlagLoading:
            print("Flag for '%s': Bar width=%1.2f, Flag width=%d px, zoom "
                  "val = %1.2f"%(countryname,barwidth,flagwidth_px,zoomval))
        imagebox = OffsetImage(flagimg, zoom=zoomval)
        imagebox.image.axes = ax
        ab = AnnotationBbox(imagebox,
                            (x+barwidth/2.0, int(barheight)), # Position in graph
                            xybox=(0.0, -20.0), # Position relative top of bar
                            frameon=False, xycoords='data',  boxcoords="offset points", pad=0)
        ax.add_artist(ab)
        ax.text(x+barwidth/2.0, 1.002*barheight,'%d' % int(barheight), ha='center', va='bottom')
        ax.set_ylabel("Number of QSO:s")

    if freqband == "all":
        ax.set_title("Top-%d countries (of %d countries, %s @ all bands, "
             "%d QSO:s)"%(numCountries, totNumCountries, txmode, totNumQSOs))
    else:
        ax.set_title("Top-%d countries (of %d countries, %s @ %s, %d QSO:s)"
             %(numCountries, totNumCountries, txmode, freqband, totNumQSOs))
    plt.xlim(0.0-0.6, numCountries-0.4)
    plt.show()
    plt.savefig("countryhistogram_qso_%s_%s.png"%(txmode, freqband), dpi=300)

    """
    Write a summary of the top accumulated distances to terminal output as
    well as to file.
    """
    f = open("countryhistogram_qso_%s_%s.txt"%(txmode, freqband), "w")
    print("="*80)
    print("Top-%d countries and number of QSO:s for mode %s at %s"%(numCountries,txmode,freqband))
    print("="*80)
    f.write("="*80+"\n")
    f.write("Top-%d countries and number of QSO:s for mode %s at %s\n"%(numCountries,txmode,freqband))
    f.write("="*80+"\n")
    for k, country in enumerate(worked_countries):
        print("#%-4d %30s %5d QSO:s"%(k+1, country,worked_countries[country]))
        f.write("#%-4d %30s %5d QSO:s\n"%(k+1, country,worked_countries[country]))
    print("="*80)
    f.write("="*80+"\n")
    f.close()

    """
    PART II - HISTOGRAM OF ACCUMULATED WORKED DISTANCE COUNTRY-WISE
    
    Sort the dictionary of accumulated distances in descending order and keep
    the highest numCountries ones. Again, when creating the sorted dictionary,
    we rely on that dictionaries preserve insertion order since Python 3.7+.
    """
    worked_distance = dict(sorted(
        worked_distance.items(), key=lambda item: item[1], reverse=True))
    worked_distance = {k: worked_distance[k]
                        for k in list(worked_distance)[:numCountries]}
    countryNumbers = range(len(worked_distance))  # Indexing along x-axis
    countryNames = worked_distance.keys()
    countryDistances = worked_distance.values()

    histLabels = []
    for k, countryName in enumerate(countryNames):
        histLabels.append("%s (%s)"%(countryName,worked_countries_all[countryName]))

    fig = plt.figure("Top-%d countries for distance (%s @ %s)"
             %(numCountries, txmode, freqband), figsize=(16,12))
    ax = fig.add_subplot(111)
    plt.xticks(countryNumbers, histLabels, rotation=45, ha='right')
    barplot = ax.bar(countryNumbers, countryDistances, color='forestgreen')
    ax.set_yscale('log')
    ax.grid(which='both', axis='y', color='#CCCCCC', linestyle='-', alpha=0.5)
    for idx, countryname in enumerate(worked_distance):
        countrycode = iso3166code[countryname]
        flagfilename = 'flags/region-flags-gh-pages/png/%s.png'%countrycode
        x = barplot[idx].get_x()
        barwidth = barplot[idx].get_width()
        barheight = barplot[idx].get_height()
        flagimg = image.imread(flagfilename)
        flagwidth_px = flagimg.shape[1]  # Num columns = width of array in px
        zoomval = 750.0*barwidth/(flagwidth_px*numCountries)
        imagebox = OffsetImage(flagimg, zoom=zoomval)
        imagebox.image.axes = ax
        ab = AnnotationBbox(imagebox,
                            (x+barwidth/2.0, int(barheight)), # Position in graph
                            xybox=(0.0, -20.0), # Position relative top of bar
                            frameon=False, xycoords='data',  boxcoords="offset points", pad=0)
        ax.add_artist(ab)
        ax.text(x+barwidth/2.0, 1.002*barheight,'%d' % int(barheight),
                ha='center', va='bottom', rotation=45)
        numQSOs = worked_countries_all[countryname]
#        ax.text(x+barwidth/2.0, 1.0,'Avg = %d km'%(int(barheight.astype(int))/int(numQSOs)),
#                ha='center', va='bottom', rotation=90)
        ax.set_ylabel("Accumulated distance per country [km]")

    if freqband == "all":
        ax.set_title("Top-%d countries (of %d countries, %s @ all bands, "
             "%d QSO:s)"%(numCountries, totNumCountries, txmode, totNumQSOs))
    else:
        ax.set_title("Top-%d countries (of %d countries, %s @ %s, %d QSO:s)"
             %(numCountries, totNumCountries, txmode, freqband, totNumQSOs))
    plt.xlim(0.0-0.6, numCountries-0.4)
    plt.show()
    plt.savefig("countryhistogram_distance_%s_%s.png"%(txmode, freqband), dpi=300)

    """
    Write a summary of the top accumulated distances to terminal output as
    well as to file.
    """
    f = open("countryhistogram_distance_%s_%s.txt"%(txmode, freqband), "w")
    print("="*80)
    print("Top-%d countries and accumulated distance for mode %s at %s"%(numCountries,txmode,freqband))
    print("="*80)
    f.write("="*80+"\n")
    f.write("Top-%d countries and accumulated distance for mode %s at %s\n"%(numCountries,txmode,freqband))
    f.write("="*80+"\n")
    for k, country in enumerate(worked_distance):
        print("#%-4d %30s %9d km (Avg %5d km, over %d QSO:s)"
              %(k+1, country, worked_distance[country],
                float(worked_distance[country])/float(worked_countries_all[country]),
                worked_countries_all[country]))
        f.write("#%-4d %30s %9d km (Avg %5d km, over %d QSO:s)\n"
              %(k+1, country, worked_distance[country],
                float(worked_distance[country])/float(worked_countries_all[country]),
                worked_countries_all[country]))
    print("="*80)
    f.write("="*80+"\n")
    f.close()

    return

def distancehistogram(qsos, txmode="FT8", freqband="all"):
    t = np.array([], dtype='datetime64[s]')
    distance = np.array([], dtype='float')
    country = np.array([], dtype='str')
    continent = np.array([], dtype='str')
    gridsquare = np.array([], dtype='str')
    call = np.array([], dtype='str')
    for qso in qsos:
        if freqband == qso["BAND"] or freqband == "all":
            if txmode == qso["MODE"] or txmode == "all":
                if qso["distance"] != "": # If the log entry contains a nonzero distance at all
                    t = np.append(t, datetime.datetime.strptime(qso["QSO_DATE"], "%Y%m%d").date())
                    distance = np.append(distance, float(qso["distance"]))
                    country = np.append(country, qso["country"])
                    continent = np.append(continent, qso["continent"])
                    gridsquare = np.append(gridsquare, qso["GRIDSQUARE"])
                    call = np.append(call, qso["CALL"])

    maximum_distance = np.max(distance)
    mean_distance = statistics.mean(distance)
    median_distance = statistics.median(distance)
    p5_percentile_distance = np.percentile(distance, 5.0)
    p95_percentile_distance = np.percentile(distance, 95.0)

    max_index = np.argmax(distance) # Index of array at which maximum distance was logged
    total_num_qsos_with_logged_distance = len(distance)
    sorted_dates = sorted(t)
    date_first_str = sorted_dates[0].strftime("%d %B, %Y")
    date_last_str = sorted_dates[-1].strftime("%d %B, %Y")
    print("Worked a total of %d distance-logged QSO:s (mode=%s, band=%s) between %s and %s."
          %(total_num_qsos_with_logged_distance, txmode, freqband, date_first_str, date_last_str))
    print("Maximum distance: %1.1f km (mode=%s, band=%s)."%(maximum_distance, txmode, freqband))
    print("Mean distance: %1.1f km (mode=%s, band=%s)."%(mean_distance, txmode, freqband))
    print("Median distance: %1.1f km (mode=%s, band=%s)."%(median_distance, txmode, freqband))
    print("5%% distance: %1.1f km (mode=%s, band=%s)."%(p5_percentile_distance, txmode, freqband))
    print("95%% distance: %1.1f km (mode=%s, band=%s)."%(p95_percentile_distance, txmode, freqband))
    print("The maximum distance of %1.1f km (mode=%s, band=%s) was logged on %s (%s, %s, %s, %s)"
          %(distance[max_index], txmode, freqband, t[max_index].strftime("%d %B, %Y"), 
            country[max_index], continent[max_index], gridsquare[max_index], call[max_index]))

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    ax.plot(t, distance, marker='.', linestyle='none', color=color)
    ax.plot(t[max_index], distance[max_index], marker='X', linestyle='none', color='red')
    ax.text(t[max_index], distance[max_index], "  %d km at %s (%s, %s, %s)"
            %(distance[max_index], t[max_index].strftime("%d %B, %Y"),
              country[max_index], gridsquare[max_index], call[max_index]),
            fontdict=None, ha='left', color='red')
    ax.xaxis_date()
    ax.set_xlabel("Date")
    ax.set_ylabel("Distance [km]")
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("Distance over time (%s @ all bands, %d QSO:s)"%(txmode, total_num_qsos_with_logged_distance))
    else:
        ax.set_title("Distance over time (%s @ %s, %d QSO:s)"%(txmode, freqband, total_num_qsos_with_logged_distance))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("distancehistogram_scatter_%s_%s.png"%(txmode, freqband), dpi=300)

    fig, ax = plt.subplots(figsize=(14,8))
    color = 'tab:blue'
    bins=np.linspace(0, math.floor(maximum_distance), 50)
    ax.hist(distance, bins=bins, rwidth=0.8, align = "left", density=False) # Set density=False to generate counts instead of probability
    ax.axvline(x=median_distance, color = 'tab:red', label='Median = %d km'%median_distance)
    ax.axvline(x=p5_percentile_distance, color = 'black', label='5%% percentile = %d km'%p5_percentile_distance)
    ax.axvline(x=p95_percentile_distance, color = 'black', label='95%% percentile = %d km'%p95_percentile_distance)
    ax.legend(loc='upper right')
    ax.set_xlabel("Distance")
    ax.set_ylabel("Frequency")
    plt.autoscale(enable=True, axis='x', tight=True)
    if freqband == "all":
        ax.set_title("Histogram over distances (%s @ all bands, %d QSO:s); max %d km at %s (%s, %s, %s)"
                     %(txmode, total_num_qsos_with_logged_distance,
                       distance[max_index], t[max_index].strftime("%d %B, %Y"),
                       country[max_index], gridsquare[max_index], call[max_index]))
    else:
        ax.set_title("Histogram over distances (%s @ %s, %d QSO:s); max %d km at %s (%s, %s, %s)"
                     %(txmode, freqband, total_num_qsos_with_logged_distance,
                       distance[max_index], t[max_index].strftime("%d %B, %Y"),
                       country[max_index], gridsquare[max_index], call[max_index]))
    ax.grid(b=True, alpha=0.5)
    plt.show()
    plt.savefig("distancehistogram_hist_%s_%s.png"%(txmode, freqband), dpi=300)

    return

def main() -> None:
    if len(sys.argv) > 1:
        adif_filename = sys.argv[1]
    else:
        # adif_filename = "/home/frejon/.local/share/WSJT-X/wsjtx_log.adi"
        adif_filename = "wsjtx_log.adi"
    summarize(adif_filename, "adiflogsummary.txt")

    qsos, header = qsodict(adif_filename)
    for txmode in {"FT8"}:
        gridsummary(qsos)
    #    for freqband in {"15m","20m","40m","all"}:
    #        qsomap(qsos, txmode=txmode, freqband=freqband, sizes="qsodensity", colors="txlevel")
    #        qsomap(qsos, txmode=txmode, freqband=freqband, sizes="uniform", colors="qsodensity")
    #    for freqband in {"15m","20m","40m"}:
    #        qsograph(qsos, graphmode="DistanceVsRx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="DistanceVsTx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="DistanceVsTxMinusRx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="DistanceVsQSODuration", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="TimeToSunsetVsRx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="TimeToSunsetVsTx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="TimeToSunsetVsTxMinusRx", txmode=txmode, freqband=freqband)
    #        qsograph(qsos, graphmode="TimeToSunsetVsDistance", txmode=txmode, freqband=freqband)
    #    for freqband in {"15m","20m","40m","all"}:
    #        qsohistogram(qsos, txmode=txmode, freqband=freqband)
        for freqband in {"15m","20m","40m","all"}:
            gridhistogram(qsos, txmode=txmode, freqband=freqband)
    #    for freqband in {"15m","20m","40m","all"}:
    #        distancehistogram(qsos, txmode=txmode, freqband=freqband)
    #    for freqband in {"15m","20m","40m","all"}:
    #        countryhistogram(qsos, txmode=txmode, freqband=freqband, numCountries=30)

    return

if __name__ == "__main__":
    main()
