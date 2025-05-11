import traceback
import sys
import time
import os.path
import threading

import obspy
from scipy.optimize import curve_fit
import boto3
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate
from scipy.ndimage import gaussian_filter1d
from scipy import signal
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client
from obspy import read_inventory,Trace,Stream
from obspy.geodetics import  gps2dist_azimuth
from obspy.core.inventory import Channel
import geopandas
from sqlalchemy import create_engine
from shapely import Point
import json

import numpy as np



#net={"st_ch": [["BR.ESM02..HHE", "BR.ESM02..HHN", "BR.ESM02..HHZ"], ["BR.ESM03..HHE", "BR.ESM03..HHN", "BR.ESM03..HHZ"], ["BR.ESM05..HHE", "BR.ESM05..HHN", "BR.ESM05..HHZ"], ["BR.ESM07..HHE", "BR.ESM07..HHN", "BR.ESM07..HHZ"], ["SC.MAC01.00.HNE", "SC.MAC01.00.HNN", "SC.MAC01.00.HNZ"], ["SC.MAC02.00.HNE", "SC.MAC02.00.HNN", "SC.MAC02.00.HNZ"], ["SC.MAC03.00.HNE", "SC.MAC03.00.HNN", "SC.MAC03.00.HNZ"], ["SC.MAC05.00.HNE", "SC.MAC05.00.HNN", "SC.MAC05.00.HNZ"], ["SC.MAC06.00.HNE", "SC.MAC06.00.HNN", "SC.MAC06.00.HNZ"], ["SC.MAC07.00.HNE", "SC.MAC07.00.HNN", "SC.MAC07.00.HNZ"], ["SC.MAC08.00.HNE", "SC.MAC08.00.HNN", "SC.MAC08.00.HNZ"], ["SC.MAC09.00.HNE", "SC.MAC09.00.HNN", "SC.MAC09.00.HNZ"], ["SC.MAC10.00.HNE", "SC.MAC10.00.HNN", "SC.MAC10.00.HNZ"], ["SC.MAC16.00.HNE", "SC.MAC16.00.HNN", "SC.MAC16.00.HNZ"], ["LK.BRK0..EHE", "LK.BRK0..EHN", "LK.BRK0..EHZ"], ["LK.BRK1..EHE", "LK.BRK1..EHN", "LK.BRK1..EHZ"], ["LK.BRK2..EHE", "LK.BRK2..EHN", "LK.BRK2..EHZ"], ["LK.BRK3..EHE", "LK.BRK3..EHN", "LK.BRK3..EHZ"], ["LK.BRK4..EHE", "LK.BRK4..EHN", "LK.BRK4..EHZ"],["BR.ESM01..HH1", "BR.ESM01..HH2", "BR.ESM01..HHZ"], ["BR.ESM08..HH1", "BR.ESM08..HH2", "BR.ESM08..HHZ"], ["BR.ESM10..HH1", "BR.ESM10..HH2", "BR.ESM10..HHZ"], ["SC.MAC04.01.HNE", "SC.MAC04.01.HNN", "SC.MAC04.01.HNZ"], ["SC.MAC11.01.HNE", "SC.MAC11.01.HNN", "SC.MAC11.01.HNZ"], ["SC.MAC12.01.HNE", "SC.MAC12.01.HNN", "SC.MAC12.01.HNZ"], ["SC.MAC13.01.HNE", "SC.MAC13.01.HNN", "SC.MAC13.01.HNZ"],  ["SC.MAC15.01.HNE", "SC.MAC15.01.HNN", "SC.MAC15.01.HNZ"],["BR.ESM04..GP1", "BR.ESM04..GP2", "BR.ESM04..GPZ"],["BR.ESM06..GP1", "BR.ESM06..GP2", "BR.ESM06..GPZ"],["BR.ESM09..GP1", "BR.ESM09..GP2", "BR.ESM09..GPZ"]], "st": ["BR.ESM02", "BR.ESM03", "BR.ESM05", "BR.ESM07", "SC.MAC01", "SC.MAC02", "SC.MAC03", "SC.MAC05", "SC.MAC06", "SC.MAC07", "SC.MAC08", "SC.MAC09", "SC.MAC10", "SC.MAC16", "LK.BRK0", "LK.BRK1", "LK.BRK2", "LK.BRK3", "LK.BRK4","BR.ESM01", "BR.ESM08", "BR.ESM10", "SC.MAC04", "SC.MAC11", "SC.MAC12", "SC.MAC13",  "SC.MAC15","BR.ESM04","BR.ESM06","BR.ESM09"]}
s3 = boto3.client('s3')

inventory = read_inventory('/mnt/seed/stations.xml')
base_path = os.environ.get('SEED_PATH', '/mnt/seed/')
bucket = os.environ.get('S3_BUCKET','geoapp-seed-data')

cl = Client(base_path)
tb='detections'
schema='sara4_test'

def connectDB():
    db_connection_url = "postgresql://postgres:wave*worm@88.99.137.51:5432/maceio_tests"
    engine = create_engine(db_connection_url)
    return engine

def toPath(nslc, t):
    # nslc = ['BR','ESM10','','HHZ']
    # path = 2022/BR/ESM10/HHZ.D/BR.ESM10..HHZ.D.2022.146
    return '%d/%s/%s/%s.D/%s.%s.%s.%s.D.%s.%03d' % (t.year, nslc[0], nslc[1], nslc[3], nslc[0], nslc[1], nslc[2], nslc[3], t.year, t.julday)

def remove_old_1d(nslc,tts):
    ts=tts-24*3600
    path = toPath(nslc, ts)
    fullpath = '%s%s' % (base_path, path)
    if os.path.exists(fullpath):
        os.remove(fullpath)

def checkFile(nslc, ts, te,remove=True):
    remove_old_1d(nslc,ts)
    t=ts
    paths=[]
    while t<te+3600*24:
        paths.append(toPath(nslc,t))
        t+=3600*24
    #paths = set((toPath(nslc, ts), toPath(nslc, te)))
    paths=set(paths)
    res = []
    for p in paths:
        fullpath = '%s%s' % (base_path, p)
        if not os.path.exists(fullpath):
            print('path not found on EFS')
            try:
                try:
                    os.makedirs('/'.join(fullpath.split('/')[:-1]),mode=0o755,exist_ok=True)
                finally:
                    s3.download_file(bucket, p, fullpath)
            except Exception as e:
                print('PATH: ' + fullpath)
                print(traceback.format_exc())
        res.append(fullpath)
    return res



def dist_xyz(seed_id,lat,lon,z):
    s1_coord=inventory.get_coordinates(seed_id)
    d=gps2dist_azimuth(s1_coord['latitude'],s1_coord['longitude'],lat,lon)
    dz=-(s1_coord['elevation']+s1_coord['local_depth'])-z
    return np.sqrt(d[0]**2+dz**2)

def afr(sensor_n,lat,lon,depth,rList):
    r=[]
    for i in sensor_n:
        jj=int(i)
        sr=rList[jj]
        s=sr.split('_')[0]
        s1=sr.split('_')[1]
        try:
            d=dist_xyz(s,lat,lon,depth)
            d1 = dist_xyz(s1, lat, lon, depth)
            r.append(d1/d)
        except:
            print('err afr')
    #print(r)
    return r

def xcorr(x, y, scale='coeff'):
    corr=np.correlate(x,y,mode='full')
    lags = np.arange(-(x.size - 1), x.size)
    if scale == 'biased':
        corr = corr / x.size
    elif scale == 'unbiased':
        corr /= (x.size - abs(lags))
    elif scale == 'coeff':
        corr /= np.sqrt(np.dot(x, x) * np.dot(y, y))
    return corr, lags

def elab_trace(trname,tts,tte,fmin,fmax,lcoffs=0.5):

    nslc = trname.split('.')
    ts = UTCDateTime(tts) - lcoffs * 3600
    te = UTCDateTime(tte) + lcoffs* 3600
    checkFile(nslc, ts, te)

    trace = cl.get_waveforms(
        nslc[0], nslc[1], nslc[2], nslc[3], ts, te).merge(1, 'interpolate')
    trace.remove_response(inventory)
    trace.resample(50, window='hann')
    trace.filter('bandpass', freqmin=fmin, freqmax=fmax, corners=2, zerophase=True)
    trace.trim(tts, tte, pad=True, fill_value=0)
    return trace[0]

def getRawData3D(nslc,tx,fmin,fmax,wnd=30):
    data=[]
    for xx in nslc:
        try:
            data+=elab_trace(xx, tx,tx+wnd, fmin,fmax).data**2
        except:
            s=elab_trace(xx, tx,tx+wnd, fmin,fmax)
            data=s.data**2
    data=np.sqrt(data)
    r=Trace(data)
    r.stats['starttime'] = s.stats['starttime']
    r.stats['sampling_rate'] = s.stats['sampling_rate']
    r.stats['network'] = s.stats['network']
    r.stats['station'] = s.stats['station']
    r.stats['location'] = s.stats['location']
    r.stats['channel'] = '3DE'
    return r

def hilb(id,tx,fmin,fmax,wnd):
    #s.filter('bandpass', freqmin=2, freqmax=25)

    s= getRawData3D(id, tx, fmin, fmax, wnd)
    data=s.data
    t=s.times()
    fs=s.stats['sampling_rate']

    analytic_signal = hilbert(data)
    amplitude_envelope = np.abs(analytic_signal)
    r=Trace(amplitude_envelope)
    r.stats['starttime'] = s.stats['starttime']
    r.stats['sampling_rate'] = fs
    r.stats['network'] = s.stats['network']
    r.stats['station'] = s.stats['station']
    r.stats['location'] = s.stats['location']
    r.stats['channel'] = 'HLB'
    return r

def gauss(id,tx,sigma,fmin,fmax,wnd):
    newid = id[0].rsplit('.', 1)[0] + '.GSS'
    try:

        f=get_filename(newid)
        r=obspy.read(f+".mseed")[0]
        if (r.stats['starttime']<=tx)and(r.stats['endtime']>=(tx+wnd)):
            r.trim(tx,tx+wnd)
            print(newid + "  found")
        else:
            raise ValueError("No GSS file")
    except:
        try:
            print(newid+" GSS no found")
            s= hilb(id,tx,fmin,fmax,wnd)
            r=Trace(gaussian_filter1d(s.data, sigma))
            r.stats['starttime'] = s.stats['starttime']
            r.stats['sampling_rate'] = s.stats['sampling_rate']
            r.stats['network'] = s.stats['network']
            r.stats['station'] = s.stats['station']
            r.stats['location'] = s.stats['location']
            r.stats['channel'] = 'GSS'
            filename=get_filename(r.id)
            r.write(filename+'.mseed')
        except:
            r=[]

    add_ch(id[0], 'GSS')
    return r

def add_ch(nslc,new_channel_code):

    # Parse the NSLC string into its components
    parts = nslc.split('.')
    network_code = parts[0]  # e.g., "BR"
    station_code = parts[1]  # e.g., "ESM02"
    location_code = parts[2]  # may be an empty string, e.g., ""
    channel_code = parts[3]  # e.g., "HHZ"

    found_station = None
    found_network = None
    found_channel = None
    for network in inventory.networks:
        if network.code == network_code:
            for station in network.stations:
                if station.code == station_code:
                    for channel in station:
                        if channel.code == channel_code:
                            found_channel = channel
                            found_station = station
                            found_network = network
                            break

    newch=found_channel
    newch.code=new_channel_code

    found_station.channels.append(newch)


def trARatio(r,corr_thr,ampl_thr):
    a = [rx.id for rx in r]
    c = []#np.empty((len(a), len(a)), dtype=object)
    s={}
    v=np.zeros([len(a),len(a)])

    for i in range(0,len(a)):
        for j in range(i+1,len(a)):
            found = False
            for k in range(0,len(a)):
                if (v[k,i]>0)and(v[k,j]>0):
                    found=True
                    break
            if  not found:
                maxsl0 = np.max(r[i].data)
                maxsl1 = np.max(r[j].data)
                if ((maxsl1 >= ampl_thr) &
                    (maxsl0 >= ampl_thr)):
                        [cc, l] = xcorr(signal.detrend(r[i].data), signal.detrend(r[j].data))
                        corr=np.max(cc)
                        if(corr>corr_thr):
                            c.append(
                                {'id':a[i]+ '_' + a[j],
                               'corr':corr,
                               'ratio':maxsl0/maxsl1,
                               'times':r[0].stats['starttime'],
                               'max0':maxsl0,
                               'max1':maxsl1})
                            v[i,j]=1
                            s[a[i]]=maxsl0;
                            s[a[j]]=maxsl1;

    return c,s

def get_filename(id):
    filename = id +"_"+ str(config['sigma']) +"_"+str(config['fmin']) +"_"+ str(config['fmax'])
    return filename

trr=[]

stzs_mid=[["BR.ESM01..HH1", "BR.ESM01..HH2", "BR.ESM01..HHZ"],
      ["BR.ESM08..HH1", "BR.ESM08..HH2", "BR.ESM08..HHZ"],
      ["BR.ESM10..HH1", "BR.ESM10..HH2", "BR.ESM10..HHZ"],
      ["SC.MAC04.01.HNE", "SC.MAC04.01.HNN", "SC.MAC04.01.HNZ"],
      ["SC.MAC11.01.HNE", "SC.MAC11.01.HNN", "SC.MAC11.01.HNZ"],
      ["SC.MAC12.01.HNE", "SC.MAC12.01.HNN", "SC.MAC12.01.HNZ"],
      ["SC.MAC13.01.HNE", "SC.MAC13.01.HNN", "SC.MAC13.01.HNZ"],
      ["SC.MAC15.01.HNE", "SC.MAC15.01.HNN", "SC.MAC15.01.HNZ"]]

stzs_sup= [["BR.ESM02..HHE", "BR.ESM02..HHN", "BR.ESM02..HHZ"],
       ["BR.ESM03..HHE", "BR.ESM03..HHN", "BR.ESM03..HHZ"],
       ["BR.ESM05..HHE", "BR.ESM05..HHN", "BR.ESM05..HHZ"],
       ["BR.ESM07..HHE", "BR.ESM07..HHN", "BR.ESM07..HHZ"],
       ["SC.MAC01.00.HNE", "SC.MAC01.00.HNN", "SC.MAC01.00.HNZ"],
       ["SC.MAC02.00.HNE", "SC.MAC02.00.HNN", "SC.MAC02.00.HNZ"],
       ["SC.MAC03.00.HNE", "SC.MAC03.00.HNN", "SC.MAC03.00.HNZ"],
       ["SC.MAC05.00.HNE", "SC.MAC05.00.HNN", "SC.MAC05.00.HNZ"],
       ["SC.MAC06.00.HNE", "SC.MAC06.00.HNN", "SC.MAC06.00.HNZ"],
       ["SC.MAC07.00.HNE", "SC.MAC07.00.HNN", "SC.MAC07.00.HNZ"],
       ["SC.MAC08.00.HNE", "SC.MAC08.00.HNN", "SC.MAC08.00.HNZ"],
       ["SC.MAC09.00.HNE", "SC.MAC09.00.HNN", "SC.MAC09.00.HNZ"],
       ["SC.MAC10.00.HNE", "SC.MAC10.00.HNN", "SC.MAC10.00.HNZ"],
       ["SC.MAC16.00.HNE", "SC.MAC16.00.HNN", "SC.MAC16.00.HNZ"],
       ["LK.BRK0..EHE", "LK.BRK0..EHN", "LK.BRK0..EHZ"],
       ["LK.BRK1..EHE", "LK.BRK1..EHN", "LK.BRK1..EHZ"],
       ["LK.BRK2..EHE", "LK.BRK2..EHN", "LK.BRK2..EHZ"],
       ["LK.BRK3..EHE", "LK.BRK3..EHN", "LK.BRK3..EHZ"],
       ["LK.BRK4..EHE", "LK.BRK4..EHN", "LK.BRK4..EHZ"]]




#config_json=json.dumps(config)
#with open("config.json", "w") as f:
#    json.dump(config, f, indent=4)
try:
    conf_file_name=sys.argv[1]
except:
    conf_file_name='config.json'

ts=UTCDateTime(sys.argv[2])
te=UTCDateTime(sys.argv[3])

with open(conf_file_name, "r") as f:
    config=json.load(f)
config_json=json.dumps(config)
note=config['name']

with connectDB().connect() as con:
    sql = "delete from "+schema+"."+tb+" where note='"+note+"'"
    con.execute(sql)
con.close()

t=ts

if config['stz_depth']=='sup':
    stzs=stzs_sup
if config['stz_depth']=='mid':
    stzs=stzs_mid

while t<te:
    inventory = read_inventory('/mnt/seed/stations.xml')
    stz_ph_s=Stream()
    for stz in stzs:
        phg=gauss(stz,t,sigma=config['sigma'],fmin=config['fmin'],fmax=config['fmax'],wnd=3600)
        if len(phg)>0:
            stz_ph_s.append(phg)
    for stz_ph_sx in stz_ph_s.slide(window_length=config['wnd'],step=config['shift']):
        x,s=trARatio(stz_ph_sx,corr_thr=config['corr_thr'],ampl_thr=config['ampl_thr'])
        ids=[cx['id'] for cx in x]
        ratios=[cx['ratio'] for cx in x]
        if len(ratios)>config['min_s']:
            print(ratios)
            def afrl(sensor_n,lat,lon,depth):
                return afr(sensor_n,lat,lon,depth,ids)
            try:
                a = curve_fit(afrl,np.arange(0,len(ids),dtype=float) , ratios, bounds=([-9.65, -35.76, -1500], [-9.62, -35.73, 0]),
                              method=config['method'],tr_solver=config['tr_solver'], full_output=True,
                              maxfev=config['maxfev'], ftol=config['ftol'], loss=config['loss'])
                ds = [dist_xyz(kii, a[0][0], a[0][1], a[0][2]) for kii in list(s.keys())]
                ampls = [s[kii] for kii in list(s.keys())]
                source_ampl = np.mean(np.asarray(ampls) * np.asarray(ds))
                perr = np.sqrt(np.diag(a[1]))

                pdb = {'utc_time': x[0]['times'],
                       'geometry': Point(a[0][1], a[0][0]),
                       'depth': a[0][2],
                       'lat': a[0][0],
                       'lon': a[0][1],
                       'note': note,
                       'config':config_json,
                       'err_lat': perr[0],
                       'err_lon': perr[1],
                       'err_depth': perr[2],
                       'n_st': len(ampls),
                       'sts': '['+','.join(s.keys())+']',
                       'ampl':'['+','.join([str(s[k]) for k in s.keys()])+']',
                       'misfit': 0,
                       'source_ampl': source_ampl,
                       'b': 0}
                print(pdb)
                gdf = geopandas.GeoDataFrame([pdb], crs="EPSG:4326")
                #gdf["geometry"] = gdf["geometry"].apply(lambda geom: geom.wkt)
                gdf.to_postgis('detections',connectDB(), 'sara4_test', 'append')
            except Exception as e:
                print(e)
    t+=3600



