import traceback
import sys
import time
import os.path
import threading

import boto3
from scipy.signal import hilbert
from obspy.signal.cross_correlation import correlate
from scipy.ndimage import gaussian_filter1d
from scipy import signal
from obspy import UTCDateTime
from obspy.clients.filesystem.sds import Client
from obspy import read_inventory,Trace,Stream
import numpy as np



#net={"st_ch": [["BR.ESM02..HHE", "BR.ESM02..HHN", "BR.ESM02..HHZ"], ["BR.ESM03..HHE", "BR.ESM03..HHN", "BR.ESM03..HHZ"], ["BR.ESM05..HHE", "BR.ESM05..HHN", "BR.ESM05..HHZ"], ["BR.ESM07..HHE", "BR.ESM07..HHN", "BR.ESM07..HHZ"], ["SC.MAC01.00.HNE", "SC.MAC01.00.HNN", "SC.MAC01.00.HNZ"], ["SC.MAC02.00.HNE", "SC.MAC02.00.HNN", "SC.MAC02.00.HNZ"], ["SC.MAC03.00.HNE", "SC.MAC03.00.HNN", "SC.MAC03.00.HNZ"], ["SC.MAC05.00.HNE", "SC.MAC05.00.HNN", "SC.MAC05.00.HNZ"], ["SC.MAC06.00.HNE", "SC.MAC06.00.HNN", "SC.MAC06.00.HNZ"], ["SC.MAC07.00.HNE", "SC.MAC07.00.HNN", "SC.MAC07.00.HNZ"], ["SC.MAC08.00.HNE", "SC.MAC08.00.HNN", "SC.MAC08.00.HNZ"], ["SC.MAC09.00.HNE", "SC.MAC09.00.HNN", "SC.MAC09.00.HNZ"], ["SC.MAC10.00.HNE", "SC.MAC10.00.HNN", "SC.MAC10.00.HNZ"], ["SC.MAC16.00.HNE", "SC.MAC16.00.HNN", "SC.MAC16.00.HNZ"], ["LK.BRK0..EHE", "LK.BRK0..EHN", "LK.BRK0..EHZ"], ["LK.BRK1..EHE", "LK.BRK1..EHN", "LK.BRK1..EHZ"], ["LK.BRK2..EHE", "LK.BRK2..EHN", "LK.BRK2..EHZ"], ["LK.BRK3..EHE", "LK.BRK3..EHN", "LK.BRK3..EHZ"], ["LK.BRK4..EHE", "LK.BRK4..EHN", "LK.BRK4..EHZ"],["BR.ESM01..HH1", "BR.ESM01..HH2", "BR.ESM01..HHZ"], ["BR.ESM08..HH1", "BR.ESM08..HH2", "BR.ESM08..HHZ"], ["BR.ESM10..HH1", "BR.ESM10..HH2", "BR.ESM10..HHZ"], ["SC.MAC04.01.HNE", "SC.MAC04.01.HNN", "SC.MAC04.01.HNZ"], ["SC.MAC11.01.HNE", "SC.MAC11.01.HNN", "SC.MAC11.01.HNZ"], ["SC.MAC12.01.HNE", "SC.MAC12.01.HNN", "SC.MAC12.01.HNZ"], ["SC.MAC13.01.HNE", "SC.MAC13.01.HNN", "SC.MAC13.01.HNZ"],  ["SC.MAC15.01.HNE", "SC.MAC15.01.HNN", "SC.MAC15.01.HNZ"],["BR.ESM04..GP1", "BR.ESM04..GP2", "BR.ESM04..GPZ"],["BR.ESM06..GP1", "BR.ESM06..GP2", "BR.ESM06..GPZ"],["BR.ESM09..GP1", "BR.ESM09..GP2", "BR.ESM09..GPZ"]], "st": ["BR.ESM02", "BR.ESM03", "BR.ESM05", "BR.ESM07", "SC.MAC01", "SC.MAC02", "SC.MAC03", "SC.MAC05", "SC.MAC06", "SC.MAC07", "SC.MAC08", "SC.MAC09", "SC.MAC10", "SC.MAC16", "LK.BRK0", "LK.BRK1", "LK.BRK2", "LK.BRK3", "LK.BRK4","BR.ESM01", "BR.ESM08", "BR.ESM10", "SC.MAC04", "SC.MAC11", "SC.MAC12", "SC.MAC13",  "SC.MAC15","BR.ESM04","BR.ESM06","BR.ESM09"]}
s3 = boto3.client('s3')

inventory = read_inventory('/mnt/seed/stations.xml')
base_path = os.environ.get('SEED_PATH', '/mnt/seed/')
bucket = os.environ.get('S3_BUCKET','geoapp-seed-data')

cl = Client(base_path)


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
    paths = set((toPath(nslc, ts), toPath(nslc, te)))
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

def hilb(s):
    #s.filter('bandpass', freqmin=2, freqmax=25)
    s.resample(200,window='hann')
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

def gauss(s,sigma=30):
    r=Trace(gaussian_filter1d(s.data, sigma))
    r.stats['starttime'] = s.stats['starttime']
    r.stats['sampling_rate'] = s.stats['sampling_rate']
    r.stats['network'] = s.stats['network']
    r.stats['station'] = s.stats['station']
    r.stats['location'] = s.stats['location']
    r.stats['channel'] = 'GSS'
    return r



def trARatio(r,corr_thr,ampl_thr):
    a = [rx.id for rx in r]
    c = np.empty((len(a), len(a)), dtype=object)
    for i in range(0,len(a)):
        for j in range(i+1,len(a)):
    #        for item in zip(r[i].slide(cWnd,cSft),r[j].slide(cWnd,cSft)):
            [cc, l] = xcorr(signal.detrend(r[i].data), signal.detrend(r[j].data))
            corr=np.max(cc)
            maxsl0=np.max(r[i].data)
            maxsl1=np.max(r[j].data)
            if(
                    (maxsl1>=ampl_thr)&
                    (maxsl0 >= ampl_thr)&
                    (corr>corr_thr)
                ):
                c[i,j] = {'id':a[i]+ '_' + a[j],
                   'corr':corr,
                   'ratio':maxsl0/maxsl1,
                   'times':r[0].stats['starttime'],
                   'max0':maxsl0,
                   'max1':maxsl1}
    # print('x')
    return c

trr=[]
#p=elab_trace('BR.ESM02..HHE',UTCDateTime('2024-03-02T00:00'),
#             UTCDateTime('2024-03-02T01:00'),1,10)
stzs=[['BR.ESM02..HHE','BR.ESM02..HHN','BR.ESM02..HHZ'],['BR.ESM03..HHE','BR.ESM03..HHN','BR.ESM03..HHZ']]
stz_ph_s=Stream()
for stz in stzs:
    p=getRawData3D(stz,UTCDateTime('2024-03-02T00:00'),fmin=1,fmax=10,wnd=30)
    ph=hilb(p)
    phg=gauss(ph,sigma=30)
    stz_ph_s.append(phg)
    x=trARatio(stz_ph_s,corr_thr=0.5,ampl_thr=1e-5)


print(p)


