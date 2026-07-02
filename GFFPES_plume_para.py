# -*- coding: utf-8 -*-
"""
Simplified CFFEPS processor for the HEMCO thermodynamic plume rise extension.

Reads CFFEPS CSV output and writes NetCDF files containing all variables needed
for any plume rise run mode in a single set of files:

  Always written:
    ECO          : CO emission rate         [kg/m2/s]
    BurnAreaTot  : cumulative burned area   [km2]   (used for mplume column mass)
    Zplume       : CFFEPS plume height      [m]     (max across fires in cell)
    Qo           : accumulated fire energy  [J]     (from CFFEPS CSV; for use_qo mode)
    Qplume       : per-timestep fire energy [J]     (from CFFEPS CSV; for validation)
    TFC          : total fuel consumed      [kg/m2] (for calc_qo mode)
    FireGrowth   : area grown this hour     [km2]   (for calc_qo mode)

  tfc_agg_mode controls how TFC, FireGrowth, Qo, Qplume are aggregated when
  multiple fires share a grid cell in the same hour:

    'weighted'  : TFC = sum(TFC_i * Growth_i) / sum(Growth_i)  (area-weighted average)
                  FireGrowth = sum(Growth_i)
                  Qo = sum(Qo_i),  Qplume = sum(Qplume_i)
                  Plume extension gets H * TFC * FireGrowth * 1e6 * 0.2
                                     = H * sum(TFC_i * Growth_i) * 1e6 * 0.2

    'max_fire'  : TFC, FireGrowth, Qo, Qplume all taken from the single fire
                  with the largest Growth in the cell.

  NOTE: BurnAreaTot uses cumulative Area(t) (the total fire footprint), NOT growth area.
        This is needed for the column mass mplume in the Anderson algorithm.

Based on inprogress_GFFPES_para.py by Robin Stevens.
Simplified and extended by J. Landers.
"""

import datetime as dt
import os
from pathlib import Path
import time
from concurrent.futures import ProcessPoolExecutor

import netCDF4 as nc
import pandas as pd
import numpy as np
import tqdm

# --- PARAMETERS ---

# location of CFFEPS csv input files
filedir = '/data/high_res/CTM/CFFEPS/2019/'
# location for output nc files
outdir = '/data/ctm/HEMCO/CFFEPS/PlumeDataQPlume/'

# NOTE on startdate / enddate:
# These are the UTC forecast dates you want in the OUTPUT files, not the input file dates.
# Each output file for date D is built from the input file named D+1, because CFFEPS files
# are named by the day the fire reports were processed (one day after the fire activity).
# Example: to produce output for April 23-29, set startdate=Apr 23 and enddate=Apr 29.
#          The script will read input files Apr 24 through Apr 30.
startdate = dt.datetime(2019, 1, 1)  # First UTC forecast date to produce output for
enddate   = dt.datetime(2019, 12, 31)  # Last  UTC forecast date to produce output for (inclusive)

dohourly = True  # Should we do hourly output or daily means?
                 # NOTE: hourly mode (True) is required for the plume rise extension,
                 #       since Qo is a per-hour accumulating value.

lat_res = 2.0  # latitude  resolution of output (degrees)
lon_res = 2.5  # longitude resolution of output (degrees)

# Define the boundaries of the output grid (Canada domain).
latmin = 24.
latmax = 71.
lonmin = -170.
lonmax = -52.5

eps = 1e-20  # small number to avoid floating-point exclusion of max boundary

# Aggregation mode for TFC, FireGrowth, Qo, Qplume:
#   'weighted' — area-weighted average TFC; sum FireGrowth, Qo, Qplume across all fires
#   'max_fire' — all four fields taken from the single fire with the largest Growth
#
# Use 'weighted' for the calc_qo agg extension. It sets
#   TFC        = sum(TFC_i * Growth_i) / sum(Growth_i)   (area-weighted mean)
#   FireGrowth = sum(Growth_i)
# so the extension's per-hour energy TFC * FireGrowth = sum(TFC_i * Growth_i), i.e. it sums
# every fire's energy contribution in the cell — consistent with BurnAreaTot and CO, which are
# already summed across all fires. The actual energy calculation stays in the extension
# (TFC * FireGrowth * H * Fire_Eff). 'max_fire' keeps only the single largest-growth fire and
# therefore under-represents energy in multi-fire cells (see qo_reset_condition_bias.md item 11).
tfc_agg_mode = 'weighted'

# --- INITIALIZE ---

ndays = (enddate - startdate).days + 1

# identify the resolution string for the output filename
if lat_res == 1. and lon_res == 1.:
    resstr = '_1x1'
elif lat_res == 2. and lon_res == 2.5:
    resstr = '_2x25'
elif lat_res == 0.1 and lon_res == 0.1:
    resstr = '_01x01'
elif lat_res == 0.5 and lon_res == 0.625:
    resstr = '_05x0625'
else:
    resstr = ''

# make a grid of latitude-longitude values
lats = np.arange(latmin, latmax + eps, lat_res)
lons = np.arange(lonmin, lonmax + eps, lon_res)

nlats = len(lats)
nlons = len(lons)

# Calculate the area of each grid cell in m2.
# Used to convert emissions to per-unit-area.
# 1 degree latitude = 111111 m; 1 degree longitude = 111111 m * cos(lat)
area = 111111 * lat_res * 111111 * lon_res * np.cos(lats * np.pi / 180.)[:, np.newaxis]


# --- PROCESS ONE DAY ---

def process_one_day(the_date: dt.datetime):
    makenc = True

    ntimes = 24 if dohourly else 1

    # Arrays always written
    co_data     = np.zeros([ntimes, nlats, nlons])  # CO emission rate (converted later)
    area_burned = np.zeros([ntimes, nlats, nlons])  # cumulative burned area [km2] -> BurnAreaTot
    zplume_data = np.zeros([ntimes, nlats, nlons])  # CFFEPS plume height [m]

    # All fields always allocated regardless of agg mode
    qo_data            = np.zeros([ntimes, nlats, nlons])  # accumulated fire energy [J]
    qplume_data        = np.zeros([ntimes, nlats, nlons])  # per-timestep fire energy [J]
    tfc_data           = np.zeros([ntimes, nlats, nlons])  # TFC [kg/m2]
    garea_burned       = np.zeros([ntimes, nlats, nlons])  # FireGrowth [km2]
    max_growth_tracker = np.zeros([ntimes, nlats, nlons])  # internal tracker, not written
    if tfc_agg_mode == 'weighted':
        tfc_x_growth   = np.zeros([ntimes, nlats, nlons])  # sum(TFC_i * Growth_i) [kg*km2/m2]
        growth_total   = np.zeros([ntimes, nlats, nlons])  # sum(Growth_i) [km2]

    # Each CFFEPS file is named one day ahead: the file for (the_date + 1) contains
    # the forecasts generated from fires detected on the_date.
    file_date = the_date + dt.timedelta(days=1)

    if makenc:

        fileformat = '2021' if the_date.year < 2022 else '2022'

        datestr = '%04i%02i%02i00' % (file_date.year, file_date.month, file_date.day)
        print("reading " + datestr)

        try:
            dataframe = pd.read_csv(filedir + datestr + "_cffeps_output.csv")
            print(len(dataframe))
        except FileNotFoundError as err:
            print(err)
            makenc = False
        else:

            for i in range(len(dataframe)):
                row = dataframe.loc[[i]]

                # Read lat/lon (column names differ between formats)
                if fileformat == '2021':
                    try:
                        thislat = row['lat']
                        thislon = row[' lon']
                    except KeyError:
                        # Some 2021 files were filled with 2022-format data
                        thislat = row['LAT']
                        thislon = row[' LON']
                        fileformat = '2022'
                elif fileformat == '2022':
                    thislat = row['LAT']
                    thislon = row[' LON']

                assert thislat.shape == (1,)
                assert thislon.shape == (1,)
                thislat = float(thislat.iloc[0])
                thislon = float(thislon.iloc[0])

                # Parse UTC forecast time
                UTC_str = row[' UTC'].iloc[0].strip()
                dt_obj  = pd.to_datetime(UTC_str, format='%Y%m%d_%H:%M')
                UTC_dt  = dt.datetime(dt_obj.year, dt_obj.month, dt_obj.day,
                                      dt_obj.hour, dt_obj.minute)

                in_grid = (thislat > latmin - lat_res / 2. and
                           thislat < latmax + lat_res / 2. and
                           thislon > lonmin - lon_res / 2. and
                           thislon < lonmax + lon_res / 2.)

                in_time = (UTC_dt >= the_date and
                           UTC_dt < the_date + dt.timedelta(hours=24))

                if in_grid and in_time:
                    lon_ind    = np.argmin(abs(lons - thislon))
                    lat_ind    = np.argmin(abs(lats - thislat))
                    time_index = UTC_dt.hour if dohourly else 0

                    if fileformat == '2021':
                        # Sum flaming + smouldering + residual CO [tonnes/hour]
                        for fire_type in 'FSR':
                            co_data[time_index, lat_ind, lon_ind] += float(row[' CO ' + fire_type].iloc[0])

                        # Cumulative area [ha -> km2]: always written as BurnAreaTot
                        area_burned[time_index, lat_ind, lon_ind] += dataframe[' Area(t)'].iloc[i] * 0.01

                        # Take max Zplume across fires in the same grid cell/hour
                        zplume_data[time_index, lat_ind, lon_ind] = max(
                            zplume_data[time_index, lat_ind, lon_ind],
                            float(dataframe[' Zplume'].iloc[i]))

                        growth_km2 = dataframe[' Growth (t)'].iloc[i] * 0.01  # ha -> km2
                        tfc_val    = dataframe[' tfc'].iloc[i]                # [kg/m2]
                        qo_val     = dataframe[' Qo'].iloc[i]                 # [J]
                        qplume_val = dataframe[' Qplume'].iloc[i]             # [J]
                        if tfc_agg_mode == 'weighted':
                            tfc_x_growth[time_index, lat_ind, lon_ind] += tfc_val * growth_km2
                            growth_total[time_index, lat_ind, lon_ind] += growth_km2
                            qo_data[time_index, lat_ind, lon_ind]      += qo_val
                            qplume_data[time_index, lat_ind, lon_ind]  += qplume_val
                        elif tfc_agg_mode == 'max_fire':
                            if growth_km2 > max_growth_tracker[time_index, lat_ind, lon_ind]:
                                max_growth_tracker[time_index, lat_ind, lon_ind] = growth_km2
                                tfc_data[time_index, lat_ind, lon_ind]    = tfc_val
                                garea_burned[time_index, lat_ind, lon_ind] = growth_km2
                                qo_data[time_index, lat_ind, lon_ind]     = qo_val
                                qplume_data[time_index, lat_ind, lon_ind] = qplume_val

                    elif fileformat == '2022':
                        co_data[time_index, lat_ind, lon_ind] += float(row[' ECO '].iloc[0])

                        try:
                            # Cumulative area [ha -> km2]
                            area_burned[time_index, lat_ind, lon_ind] += dataframe[' Area(t)'].iloc[i] * 0.01
                        except (ValueError, KeyError) as err:
                            print('area_burned 2022:', err)

                        # NOTE: tfc, Growth, Qo, Qplume column names in 2022 format not yet verified.
                        # Assuming same names as 2021 format; update if needed.
                        try:
                            growth_km2 = dataframe[' Growth (t)'].iloc[i] * 0.01
                            tfc_val    = dataframe[' tfc'].iloc[i]
                            qo_val     = dataframe[' Qo'].iloc[i]
                            qplume_val = dataframe[' Qplume'].iloc[i]
                            if tfc_agg_mode == 'weighted':
                                tfc_x_growth[time_index, lat_ind, lon_ind] += tfc_val * growth_km2
                                growth_total[time_index, lat_ind, lon_ind] += growth_km2
                                qo_data[time_index, lat_ind, lon_ind]      += qo_val
                                qplume_data[time_index, lat_ind, lon_ind]  += qplume_val
                            elif tfc_agg_mode == 'max_fire':
                                if growth_km2 > max_growth_tracker[time_index, lat_ind, lon_ind]:
                                    max_growth_tracker[time_index, lat_ind, lon_ind] = growth_km2
                                    tfc_data[time_index, lat_ind, lon_ind]    = tfc_val
                                    garea_burned[time_index, lat_ind, lon_ind] = growth_km2
                                    qo_data[time_index, lat_ind, lon_ind]     = qo_val
                                    qplume_data[time_index, lat_ind, lon_ind] = qplume_val
                        except (ValueError, KeyError) as err:
                            print('fire fields 2022:', err)

                elif not in_grid:
                    print('location not in grid: lat %.2f lon %.2f' % (thislat, thislon))

        # --- Post-loop: compute weighted TFC ---
        if tfc_agg_mode == 'weighted':
            # tfc_data = sum(TFC_i * Growth_i) / sum(Growth_i)  [kg/m2]
            # garea_burned = sum(Growth_i)  [km2]
            # Fortran: Qo = H * tfc_data * garea_burned * 1e6 * 0.2 = H * sum(TFC_i*Growth_i) * 1e6 * 0.2
            # np.divide with where= avoids evaluating the division at growth_total<=0
            # cells (np.where would compute it everywhere first -> spurious divide warning).
            tfc_data     = np.divide(tfc_x_growth, growth_total,
                                     out=np.zeros_like(tfc_x_growth), where=growth_total > 0)
            garea_burned = growth_total
            # Qo and Qplume are already summed in the loop

        # --- Unit conversion for CO ---
        # 2021 format: tonnes/hour -> g/s
        if fileformat == '2021':
            co_data *= 1e6 / 3600.

        # g/s -> kg/m2/s
        co_data /= (1e3 * area[np.newaxis, :, :])

        # For daily output, average the emission rate over 24 hours
        if not dohourly:
            co_data /= 24.

        # --- Write NetCDF ---
        datestr = '%04i%02i%02i' % (the_date.year, the_date.month, the_date.day)
        suffix  = '_hourly.nc' if dohourly else '.nc'

        ncpath = outdir + datestr[:4] + '/' + datestr[4:6]
        Path(ncpath).mkdir(parents=True, exist_ok=True)

        ncfilepath = ncpath + '/CFFEPS' + resstr + '_' + datestr + suffix
        if os.path.exists(ncfilepath):
            os.remove(ncfilepath)
        ncfile = nc.Dataset(ncfilepath, 'w')

        # Dimensions
        ncfile.createDimension('time', ntimes)
        ncfile.createDimension('lat',  nlats)
        ncfile.createDimension('lon',  nlons)

        # Time variable
        time_var           = ncfile.createVariable('time', 'f4', ('time',))
        time_var.axis      = 'T'
        time_var.units     = ('hours since %04i-%02i-%02i 00:00:00'
                              % (the_date.year, the_date.month, the_date.day))
        time_var.calendar  = 'standard'
        time_var.long_name = 'time'
        time_var[:]        = np.arange(ntimes)

        # Latitude variable
        lat_var               = ncfile.createVariable('lat', 'f4', ('lat',))
        lat_var.axis          = 'Y'
        lat_var.units         = 'degrees_north'
        lat_var.standard_name = 'Latitude'
        lat_var.long_name     = 'Latitude'
        lat_var[:]            = lats

        # Longitude variable
        lon_var               = ncfile.createVariable('lon', 'f4', ('lon',))
        lon_var.axis          = 'X'
        lon_var.units         = 'degrees_east'
        lon_var.standard_name = 'Longitude'
        lon_var.long_name     = 'Longitude'
        lon_var[:]            = lons

        # CO emissions (always)
        nc_co               = ncfile.createVariable('ECO', 'f4', ('time', 'lat', 'lon'))
        nc_co[:]            = co_data[:]
        nc_co.standard_name = 'carbon monoxide emission rate'
        nc_co.long_name     = 'carbon monoxide emission rate'
        nc_co.units         = 'kg m-2 s-1'

        # Cumulative burned area (always) — feeds mplume column mass in Fortran (eqn 16)
        nc_area               = ncfile.createVariable('BurnAreaTot', 'f4', ('time', 'lat', 'lon'))
        nc_area[:]            = area_burned[:]
        nc_area.standard_name = 'total burned area'
        nc_area.long_name     = 'total burned area'
        nc_area.units         = 'km^2'

        # Plume injection height (always)
        nc_zplume               = ncfile.createVariable('Zplume', 'f4', ('time', 'lat', 'lon'))
        nc_zplume[:]            = zplume_data[:]
        nc_zplume.standard_name = 'plume injection height'
        nc_zplume.long_name     = 'estimated plume injection height'
        nc_zplume.units         = 'm'

        # Qo — accumulated fire energy from CFFEPS CSV (for use_qo mode runs)
        nc_qo               = ncfile.createVariable('Qo', 'f4', ('time', 'lat', 'lon'))
        nc_qo[:]            = qo_data[:]
        nc_qo.standard_name = 'fire energy'
        nc_qo.long_name     = ('accumulated fire energy [%s]' % tfc_agg_mode)
        nc_qo.units         = 'J'
        nc_qo.tfc_agg_mode  = tfc_agg_mode

        # Qplume — per-timestep fire energy from CFFEPS CSV (for validation)
        nc_qplume               = ncfile.createVariable('Qplume', 'f4', ('time', 'lat', 'lon'))
        nc_qplume[:]            = qplume_data[:]
        nc_qplume.standard_name = 'plume energy'
        nc_qplume.long_name     = ('per-timestep fire energy [%s]' % tfc_agg_mode)
        nc_qplume.units         = 'J'
        nc_qplume.tfc_agg_mode  = tfc_agg_mode

        # TFC — for calc_qo mode runs
        nc_tfc               = ncfile.createVariable('TFC', 'f4', ('time', 'lat', 'lon'))
        nc_tfc[:]            = tfc_data[:]
        nc_tfc.standard_name = 'total fuel consumed per unit area'
        nc_tfc.long_name     = ('area-weighted mean TFC [weighted]' if tfc_agg_mode == 'weighted'
                                else 'TFC of largest fire in cell [max_fire]')
        nc_tfc.units         = 'kg m-2'
        nc_tfc.tfc_agg_mode  = tfc_agg_mode

        # FireGrowth — for calc_qo mode runs
        nc_growth               = ncfile.createVariable('FireGrowth', 'f4', ('time', 'lat', 'lon'))
        nc_growth[:]            = garea_burned[:]
        nc_growth.standard_name = 'fire area growth'
        nc_growth.long_name     = ('total area growth this hour [weighted]' if tfc_agg_mode == 'weighted'
                                   else 'area growth of largest fire in cell [max_fire]')
        nc_growth.units         = 'km^2'
        nc_growth.tfc_agg_mode  = tfc_agg_mode

        # Global attributes
        ncfile.setncattr('title',              'CFFEPS fire inputs for plume rise: ' + datestr)
        ncfile.setncattr('contact',            'Julian Landers')
        ncfile.setncattr('conventions',        'COARDS')
        ncfile.setncattr('tfc_agg_mode',       tfc_agg_mode)
        ncfile.setncattr('history',            'File generated on: ' + time.ctime(time.time()))
        ncfile.setncattr('productiondatetime', 'File generated on: ' + time.ctime(time.time()))
        ncfile.setncattr('format',             'NETCDF-4')
        ncfile.setncattr('delta_lon',          lon_res)
        ncfile.setncattr('delta_lat',          lat_res)

        ncfile.close()


if __name__ == "__main__":
    print("Time resolution: " + ["daily", "hourly"][dohourly])
    print("Spatial resolution (degrees): " + str(lat_res) + " x " + str(lon_res))
    print(f"Grid lat[{lats[0]}, {lats[-1]}] lon[{lons[0]}, {lons[-1]}]  (n={len(lats)}x{len(lons)})")
    print(f"tfc_agg_mode={tfc_agg_mode}")

    all_days = [startdate + dt.timedelta(days=i)
                for i in range((enddate - startdate).days + 1)]

    max_workers = min(os.cpu_count() or 1, 6)
    with ProcessPoolExecutor(max_workers=max_workers) as exe:
        list(tqdm.tqdm(exe.map(process_one_day, all_days), total=len(all_days)))
