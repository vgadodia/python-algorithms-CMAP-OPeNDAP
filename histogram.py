# Import Functions
import sys
!{sys.executable} -m pip install netCDF4
!{sys.executable} -m pip install xarray
import opedia
import netCDF4
import os
import numpy as np
import pandas as pd
import xarray as xr
import db
import subset
import common as com
from datetime import datetime, timedelta
import time
from bokeh.io import output_notebook
from math import pi
from bokeh.plotting import figure, show, output_file
from bokeh.layouts import column
from bokeh.models import DatetimeTickFormatter
from bokeh.palettes import all_palettes
from bokeh.models import HoverTool
from bokeh.embed import components
import jupyterInline as jup
if jup.jupytered():
    from tqdm import tqdm_notebook as tqdm
else:
    from tqdm import tqdm

# Histogram Plotting Function
def xarrayPlotDist(tables, variables, startDate, endDate, lat1, lat2, lon1, lon2, depth1, depth2, fname, exportDataFlag, marker='-', msize=30, clr='purple'):
    p = []
    lw = 2
    w = 800
    h = 400
    TOOLS = 'pan,wheel_zoom,zoom_in,zoom_out,box_zoom, undo,redo,reset,tap,save,box_select,poly_select,lasso_select'
    for i in tqdm(range(len(tables)), desc='overall'):
        
        unit = tables[i].variables[variables[i]].attrs['units']
        
        varData = tables[i].sel(TIME = slice(startDate, endDate), LAT_C = slice(lat1, lat2), LON_C = slice(lon1, lon2), DEP_C = slice(depth1, depth2))
        varData = varData.variables[variables[i]][:].values
        varData[varData < 0] = float('NaN')
        
        varData = varData.tolist()[0][0]
        elems = []
        for list in varData:
            for item in list:
                elems.append(item)
        varData = pd.Series(elems).dropna()

        if exportDataFlag:
            exportData(y, tables[i], variables[i], startDate, endDate, lat1, lat2, lon1, lon2, depth1, depth2)
        try:    
            carData = varData[~np.isnan(varData)]     # remove nans
        except Exception as e:
            continue   
            
        hist, edges = np.histogram(varData, density=True, bins=50)
        
        output_notebook()
        
        p1 = figure(tools=TOOLS, toolbar_location="above", plot_width=w, plot_height=h)
        p1.yaxis.axis_label = 'Density'
        p1.xaxis.axis_label = variables[i] + ' (' +str(unit) + ')'
        leg = variables[i]
        fill_alpha = 0.4   
        cr = p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:], fill_color="dodgerblue", line_color=None, hover_fill_color="firebrick", fill_alpha=fill_alpha, hover_alpha=0.7, hover_line_color="white", legend=leg)
        p1.add_tools(HoverTool(tooltips=None, renderers=[cr], mode='mouse'))
        p.append(p1)
    dirPath = 'embed/'
    if not os.path.exists(dirPath):
        os.makedirs(dirPath)        
    
    #if not inline:      ## if jupyter is not the caller
    #    output_file(dirPath + fname + ".html", title="Histogram")
    #    print('test')
    show(column(p))
    return

# File to read data from
xFile = xr.open_dataset('http://3.88.71.225:80/thredds/dodsC/las/id-a1d60eba44/data_usr_local_tomcat_content_cbiomes_20190510_20_darwin_v0.2_cs510_darwin_v0.2_cs510_nutrients.nc.jnl')

# Testing Space
tables = [xFile]           # see catalog.csv  for the complete list of tables and variable names
variables = ['O2']       # see catalog.csv  for the complete list of tables and variable names
startDate = '2010-04-30'
endDate = '2010-04-30'
lat1, lat2 = -80, 80
lon1, lon2 = -80, 80
depth1, depth2 = 0, 10
fname = 'DEP'
exportDataFlag = False      # True if you you want to download data

xarrayPlotDist(tables, variables, startDate, endDate, lat1, lat2, lon1, lon2, depth1, depth2, fname, exportDataFlag)
