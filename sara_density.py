from datetime import datetime, timedelta
import logging, sys, os

import geopandas as gpd
import numpy as np
from scipy import stats
from skimage import measure
from shapely.geometry import Polygon, LineString
from shapely.geometry.collection import GeometryCollection
import psycopg2
from sqlalchemy import create_engine
import utils
from obspy import UTCDateTime
import json
#logging.basicConfig(stream=sys.stdout, level=logging.INFO,
#                    format='%(asctime)s %(levelname)s %(message)s')

conn=utils.connectDB()

extent = (-3981126, -3977756, -1079396, -1076120)
gridStepL = 10
xmin = extent[0]  # np.min(x)
xmax = extent[1]  # np.max(x)
ymin = extent[2]  # np.min(y)
ymax = extent[3]  # np.max(y)

X, Y = np.mgrid[xmin:xmax:gridStepL, ymin:ymax:gridStepL]
stpx = (xmax - xmin) / len(X)
stpy = (ymax - ymin) / len(Y[0])

tb = 'detections'
schema = 'sara4_test'

def getSaraRaw(ts, te,config):
    maoi_bookmark = 'MAOI_buffer_100m'

    sql="SELECT d.*, d.geometry AS geom " \
        "FROM sara4_test.detections  as d " \
        "JOIN sara4_test.bookmarks b ON b.name = '"+maoi_bookmark+"' " \
        "WHERE d.note = '"+config['den_in_type']+"' " \
        " AND d.depth>="+config['den_min_depth']+" "\
        " AND d.depth<="+config['den_max_depth']+" "\
        " AND d.source_ampl>="+config['den_min_source_ampl']+" "\   
        " AND d.utc_time > '"+ts.isoformat()+"' " \
        " AND d.utc_time < '"+te.isoformat()+"' " \
        " AND ST_Within(d.geometry, b.geom);"


    dfo = gpd.read_postgis(sql, conn)
    return dfo


def sara_density(ts, te, X, Y, stpx, stpy,config,den_min, den_max, lines):
    in_type=config['den_in_type']
    print(f'ATH Density calculating: {te.isoformat()}, type: {in_type}')
    dfo = getSaraRaw(ts, te,config)
    nEv = len(dfo)
    maxMag = dfo['source_ampl'].max()

    te_ = []
    nEv_ = []
    maxMag_ = []
    dMax_ = []
    d_line_ = []
    poly_ = []

    if nEv > 3:
        df = dfo.to_crs(epsg=3857)
        x = df['geom'].x
        y = df['geom'].y

        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        Zn = Z * nEv
        dMax = np.max(Zn)
        for d_line in np.linspace(den_min, den_max, lines):
            cc = []
            poly = []
            cc.append(measure.find_contours(Zn, d_line))
            te_.append(te.isoformat())
            dMax_.append(dMax)
            nEv_.append(nEv)
            maxMag_.append(maxMag)
            d_line_.append(d_line)
            if len(cc) > 0:
                for cll in cc:
                    for c in cll:
                        lp = [[X[0][0] + pp[0] * stpx, Y[0][0] + pp[1] * stpy]
                              for pp in c]
                        poly.append(LineString(lp))
            poly_.append(GeometryCollection(poly))
    else:
        te_.append(te.isoformat())
        nEv_.append(nEv)
        maxMag_.append(maxMag)
        dMax_.append(0)
        poly_.append(GeometryCollection([]))
        d_line_.append(0)

    mm = {'utc_time': te_, 'max_density': dMax_, 'rate': nEv_, 'max_mag': maxMag_,
          'geometry': poly_, 'density': d_line_,'out_type':config['den_out_type'],
          'in_type':config['den_in_type'],'pr':'xy'}#,'athena_ids':athena_ids}
    gdf = gpd.GeoDataFrame(mm, crs="EPSG:3857")
    gdf = gdf.to_crs(epsg=4326)
    gdf.to_postgis("density", conn, schema, 'append')


t=UTCDateTime('2023-08-01T00:00')
te=UTCDateTime('2023-12-31T00:00')


ts=UTCDateTime(sys.argv[2])
te=UTCDateTime(sys.argv[3])
conf_file_name=sys.argv[4]
with open(conf_file_name, "r") as f:
    config=json.load(f)
config_json=json.dumps(config)
note=config['den_name']
wnd=config['den_wnd']#3600
sft=config['den_shift']


t=ts
while (t<te):
    try:
        sara_density(t,t+wnd,X,Y,stpx,stpy,config,0,1,1)#min_depth,max_depth,0,1,1,in_type)#'SupAgoDec23')
    except Exception as e:
        print(str(e))
    t=t+sft
print('pippo')