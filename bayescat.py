#!/usr/bin/env python

#stdlib
import os.path
import sys
from datetime import datetime,timedelta
import csv
import time
import glob
import sqlite3
import argparse
import socket
import subprocess

#local
from libcomcat.comcat import getPhaseData,getEventData
from neicio.cmdoutput import getCommandOutput
from neicutil import text

#third party
from obspy.fdsn import Client
from obspy import UTCDateTime
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from obspy.core.util.geodetics import gps2DistAzimuth

DEFAULT_RADIUS = 15
DEFAULT_START = datetime(1981,1,31,0,0,0)
DEFAULT_END = datetime.utcnow()
NWEEKS = 12
BAYESLOC = os.path.join(os.path.expanduser('~'),'bayesloc')

CWBHOST = 'cwbpub.cr.usgs.gov'
CWBPORT = 2052

CONFIG = """// The arrival data and station information:
      Data:
      {
        arrival_file = "EVENTFOLDER/arrival.dat";     // the file with the arrival data
        station_file = "EVENTFOLDER/station.dat";     // the file with the station data
      };
      // The specification of the travel-time look-up tables:
      TravelTimeModel:
      {
        tt_prefix = "BAYESLOC/ttimes/ak135."; 
        // Phase labels are mapped to look-up tables by order (Pn to iasp91.Pn, etc.)
};"""
# // parameters specifying the event-origin model and sampling procedure:
# Origin:
# {
#   // parameters specifying the prior model:
#  Prior:
#   {
#     // an input file specifying event origin priors (e.g., known events, if any):
#     file = "EVENTFOLDER/prior.dat";
#     // if not given (""), flat prior assumed for lat-lon-depth, but for time:
#     //def_time_mean_type = "first_station";    // how to set default time mean
#     //def_time_sd = 1.0e12;                    // the default origin time sd
#   };
# };"""


PHASELIST = ['P', 'Pg', 'Pn', 'pP', 'S', 'Sg', 'Sn', 'sP']

MAGARRAY = np.array([0,2.0,2.5,3.0,3.5, 4.0, 4.5, 9.9])
PICKDIST = np.array([0,3.0,6.0,8.0,10.0,12.0,18.0,90.0])

TABLES = {'station':
          {'id':'integer primary key',
           'nscl':'str',
           'lat':'float',
           'lon':'float',
           'elev':'float'},
          'event':
           {'id':'integer primary key',
            'code':'str',
            'lat':'float',
            'lon':'float',
            'depth':'float',
            'time':'datetime',
            'mag':'float',
            'rlat':'float',
            'rlon':'float',
            'rdepth':'float',
            'rtime':'datetime'},
          'phase':
           {'sid':'int',
            'eid':'int',
            'phasetype':'str',
            'time':'float'}}

BAYESDIR = os.path.join(os.path.expanduser('~'),'bayesloc')
BAYESDB = 'bayes.db'
BAYESBIN = 'bayesloc'

TIMEFMT = '%Y-%m-%dT%H:%M:%S'
DATEFMT = '%Y-%m-%d'

FOLDER_ERROR = '''You must create a folder called %s.  Under it, you must have the following contents:
  bin/bayesloc (BayesLoc executable)
  ttimes/ak135.* (Travel times files)''' % BAYESDIR

def maketime(timestring):
    outtime = None
    try:
        outtime = datetime.strptime(timestring,TIMEFMT)
    except:
        try:
            outtime = datetime.strptime(timestring,DATEFMT)
        except:
            raise Exception,'Could not parse time or date from %s' % timestring
    return outtime
  
def createTables(db,cursor):
    for table in TABLES.keys():
        sql = 'CREATE TABLE %s (' % table
        nuggets = []
        for column,ctype in TABLES[table].iteritems():
            nuggets.append('%s %s' % (column,ctype))
        sql += ','.join(nuggets) + ')'
        cursor.execute(sql)
        db.commit()

def getMapLines(dmin,dmax):
    NLINES = 4
    drange = dmax-dmin
    if drange > 4:
        near = 1
    else:
        if drange >= 0.5:
            near = 0.25
        else:
            near = 0.125
    inc = text.roundToNearest(drange/NLINES,near)
    if inc == 0:
        near = pow(10,round(np.log10(drange))) #make the increment the closest power of 10
        inc = text.ceilToNearest(drange/NLINES,near)
        newdmin = text.floorToNearest(dmin,near)
        newdmax = text.ceilToNearest(dmax,near)
    else:
        newdmin = text.ceilToNearest(dmin,near)
        newdmax = text.floorToNearest(dmax,near)
    darray = np.arange(newdmin,newdmax,inc)
    return darray

def getStationCoord(nscl):
    network,station,channel,location = nscl.split('.')
    scode = 'FDSN.%s.%s' % (network,station)
    req = '-a %s -c c \n' % scode
    pad = chr(0) * (80 - len(req))
    req = str(req + pad)
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM,0)
        s.connect((CWBHOST,CWBPORT))
        s.send(req)
        response = s.recv(10241)
        s.close()
    except Exception,msg:
        try:
            time.sleep(2)
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM,0)
            s.connect((CWBHOST,CWBPORT))
            s.send(req)
            response = s.recv(10241)
            s.close()
        except:
            return (None,None,None)
    try:
        parts = response.split('\n')
        parts = parts[0].split(':')
        parts = parts[1].split()
        lat = float(parts[0])
        lon = float(parts[1])
        elev = float(parts[2])/1000.0
    except Exception,msg:
        return (None,None,None)
    return (lat,lon,elev)

def makeMap(eventlist,lat,lon,eventfolder):
    xmin = ymin = 1e9
    xmax = ymax = -1e9
    ymin = min([min(f['lat'],f['rlat']) for f in eventlist])
    ymax = max([max(f['lat'],f['rlat']) for f in eventlist])
    xmin = min([min(f['lon'],f['rlon']) for f in eventlist])
    xmax = max([max(f['lon'],f['rlon']) for f in eventlist])

    #often the map is too zoomed in to these events, lacking geographic context
    #increase the x/y range of the map so that our events occupy the middle third of it in both dimensions
    x_range = xmax-xmin
    y_range = ymax-ymin
    newxrange = x_range * 3
    newyrange = y_range * 3
    xmin = xmin - (newxrange-x_range)/2.0
    xmax = xmax + (newxrange-x_range)/2.0
    ymin = ymin - (newyrange-y_range)/2.0
    ymax = ymax + (newyrange-y_range)/2.0
    
    #Map the input events against outputs
    clat = ymin + (ymax-ymin)/2.0
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_axes([0,0,1.0,1.0])
    bmap = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,
                   urcrnrlon=xmax,urcrnrlat=ymax,
                   resolution='h',projection='merc',lat_ts=clat,ax=ax)
    bmap.drawcoastlines()
    bmap.drawstates()
    bmap.drawcountries()    

    mer = getMapLines(xmin,xmax)
    par = getMapLines(ymin,ymax)
    
    xmap_range = bmap.xmax-bmap.xmin
    ymap_range = bmap.ymax-bmap.ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    bmap.drawmeridians(mer,labels=[0,0,1,0],fontsize=8,
                       linewidth=0.5,color='cyan',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
    bmap.drawparallels(par,labels=[0,1,0,0],fontsize=8,
                             linewidth=0.5,color='cyan',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
    
    bmap.drawmapboundary(color='k',linewidth=2.0)
    for event in eventlist:
        lat1 = event['lat']
        lon1 = event['lon']
        lat2 = event['rlat']
        lon2 = event['rlon']
        bx1,by1 = bmap(lon1,lat1)
        bx2,by2 = bmap(lon2,lat2)
        dx2 = (bx2-bx1)*(bx2-bx1)
        dy2 = (by2-by1)*(by2-by1)
        dist = np.sqrt(dx2 + dy2)
        bmap.plot(bx2,by2,'r.')
        bmap.plot([bx1,bx2],[by1,by2],'b')

    #plot the event in question with a *
    bmap.plot(lon,lat,'c*')
        
    mapfile = os.path.join(eventfolder,'shiftmap.pdf')
    plt.savefig(mapfile)
    plt.close()
    return mapfile

def makeInputMaps(eventdict,stationdict,eventfolder):
    maplist = []
    for eventid,eventinfo in eventdict.iteritems():
        xmin = ymin = 1e6
        xmax = ymax = -1e6
        ecode,elat,elon,edep,etime,emag = eventinfo
        if elat < ymin:
            ymin= elat
        if elat > ymax:
            ymax= elat
        if elon < xmin:
            xmin= elon
        if elon > xmax:
            xmax= elon
        px = []
        py = []
        for stationid,stationtuple in stationdict.iteritems():
            plat,plon,phase = stationtuple
            px.append(plon)
            py.append(plat)
            if plat < ymin:
                ymin= plat
            if plat > ymax:
                ymax= plat
            if plon < xmin:
                xmin= plon
            if plon > xmax:
                xmax= plon

        #Map the input events against outputs
        clat = ymin + (ymax-ymin)/2.0
        clon = xmin + (xmax-xmin)/2.0
        try:
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_axes([0,0,1.0,1.0])
            bmap = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,
                           urcrnrlon=xmax,urcrnrlat=ymax,
                           resolution='c',projection='tmerc',
                           lat_0=clat,lon_0=clon,ax=ax)
        except Exception,msg:
            pass
        bmap.drawcoastlines()
        bmap.drawstates()
        bmap.drawcountries()    
        
        mer = getMapLines(xmin,xmax)
        par = getMapLines(ymin,ymax)
        
        xmap_range = bmap.xmax-bmap.xmin
        ymap_range = bmap.ymax-bmap.ymin
        xoff = -0.09*(xmap_range)
        yoff = -0.04*(ymap_range)
        
        bmap.drawmeridians(mer,labels=[0,0,1,0],fontsize=8,
                           linewidth=0.5,color='cyan',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
        bmap.drawparallels(par,labels=[0,1,0,0],fontsize=8,
                           linewidth=0.5,color='cyan',yoffset=yoff,xoffset=xoff,dashes=[1,0.01])
        bmap.drawmapboundary(color='k',linewidth=2.0)
        bmap.plot(elon,elat,'rx',latlon=True)
        bmap.plot(px,py,'b.',latlon=True)
        mapfile = os.path.join(eventfolder,'%s_event.pdf' % ecode)
        plt.title('Event %s M%.1f' % (ecode,emag))
        plt.savefig(mapfile)
        plt.close()
        maplist.append(mapfile)

    return maplist

def getEventPriors(eventlist,cursor):
    sql = 'SELECT id,rlat,rlon,rdepth,rtime FROM event WHERE code="%s"'
    priors = []
    for event in eventlist:
        query = sql % event['id']
        cursor.execute(query)
        row = cursor.fetchone()
        if row is not None and row[1] is not None:
            eid = row[0]
            rlat = row[1]
            rlon = row[2]
            rdepth = row[3]
            rtime = UTCDateTime(row[4]).timestamp
            priors.append((eid,rlat,rlon,rdepth,rtime))
    return priors

def getProcessedData(eventlist,db,cursor,ndays=0):
    otime = (datetime.utcnow() - timedelta(ndays)).strftime('%Y-%m-%d %H:%M:%S')
    equery = 'SELECT id FROM event where code="%s" and time < "%s"'
    stations = {}
    arrivals = []
    culled = []
    for event in eventlist:
        cursor.execute(equery % (event['id'],otime))
        row = cursor.fetchone()
        if row is None:
            culled.append(event['id'])
            continue
        #get the arrival time stuff
        eventid = row[0]
        pquery = 'SELECT sid,phasetype,time FROM phase WHERE eid=%i' % eventid
        cursor.execute(pquery)
        prows = cursor.fetchall()
        for prow in prows:
            sid = prow[0]
            phase = prow[1]
            time = UTCDateTime(prow[2]).timestamp
            squery = 'SELECT nscl,lat,lon,elev FROM station WHERE id=%i' % sid
            cursor.execute(squery)
            srow = cursor.fetchone()
            scode = srow[0]
            lat = srow[1]
            lon = srow[2]
            elev = srow[3]
            stations[scode] = (lat,lon,elev)
            arrivals.append((eventid,scode,phase,time))
    return (stations,arrivals,culled)
                            
def insertPhaseData(phasedata,db,cursor):
    missing_stations = []
    stations = {}
    arrivals = []
    #remember that the event may already be here, and we may be just updating it
    equery = 'SELECT id FROM event WHERE code="%s"' % phasedata.eventcode
    cursor.execute(equery)
    row = cursor.fetchone()
    eid = None
    if row is not None:
        eid = row[0]
    lat = phasedata.origins[0]['lat']
    lon = phasedata.origins[0]['lon']
    depth = phasedata.origins[0]['depth']
    time = phasedata.origins[0]['time']
    mag = phasedata.magnitudes[0]['magnitude']
    if eid is None:
        efmt = 'INSERT INTO event (code,lat,lon,depth,time,mag) VALUES ("%s",%.4f,%.4f,%.1f,"%s",%.1f)'
        equery = efmt % (phasedata.eventcode,lat,lon,depth,time,mag)
    else:
        efmt = 'UPDATE event set lat=%.4f,lon=%.4f,depth=%.1f,time="%s",mag=%.1f WHERE id=%i'
        equery = efmt % (lat,lon,depth,time,mag,eid)
    cursor.execute(equery)
    db.commit()
    iquery = 'SELECT id FROM event WHERE code="%s"' % phasedata.eventcode
    cursor.execute(iquery)
    row = cursor.fetchone()
    eventid = row[0]
    mdiff = mag - MAGARRAY
    mdiff[mdiff < 0] = np.nan
    imag = mdiff.argmin()
    pdist = PICKDIST[imag]
    for phase in phasedata.phases:
        nscl = phase['nscl']
        distance = phase['distance']
        if distance > pdist:
            continue
        squery = 'SELECT id,lat,lon,elev FROM station WHERE nscl = "%s"' % nscl
        cursor.execute(squery)
        row = cursor.fetchone()
        if row is None:
            slat,slon,elev = getStationCoord(nscl)
            if slat is None:
                missing_stations.append(nscl)
                continue
            sfmt = 'INSERT INTO station (nscl,lat,lon,elev) VALUES ("%s",%.4f,%.4f,%.1f)'
            squery = sfmt % (nscl,slat,slon,elev)
            cursor.execute(squery)
            db.commit()
            squery = 'SELECT id FROM station WHERE nscl="%s"' % nscl
            cursor.execute(squery)
            sid = cursor.fetchone()[0]
        else:
            slat,slon,elev = row[1:]
            sid = row[0]
        stations[nscl] = (slat,slon,elev)
        #now insert the arrival stuff into the phase table
        ptime = phase['time']
        phase = phase['phasetype']
        if phase not in PHASELIST:
            continue
        pfmt = 'INSERT INTO phase (sid,eid,phasetype,time) VALUES (%i,%i,"%s","%s")'
        pquery = pfmt % (sid,eventid,phase,time)
        cursor.execute(pquery)
        db.commit()
        arrivals.append((eventid,nscl,phase,UTCDateTime(ptime).timestamp))
    return (stations,arrivals,missing_stations)

def getStats(cursor):
    query = 'SELECT count(*) FROM event'
    cursor.execute(query)
    nevents = cursor.fetchone()[0]
    query = 'SELECT count(*) FROM station'
    cursor.execute(query)
    nstations = cursor.fetchone()[0]
    query = 'SELECT count(*) FROM phase'
    cursor.execute(query)
    narrivals = cursor.fetchone()[0]
    return (nevents,nstations,narrivals)

def deleteEvents(db,cursor,eventcodes):
    nevents = 0
    for eventcode in eventcodes:
        query = 'SELECT id FROM event WHERE code="%s"' % eventcode
        cursor.execute(query)
        row = cursor.fetchone()
        if row is None:
            continue
        eid = row[0]
        query2 = 'DELETE FROM phase WHERE eid=%i' % eid
        cursor.execute(query2)
        db.commit()
        query3 = 'DELETE FROM event WHERE id=%i' % eid
        cursor.execute(query3)
        db.commit()
        nevents += 1
    return nevents

def main(args):
    eventid = args.id
    radius = args.radius
    #does the bayesloc folder exist?
    if not os.path.isdir(BAYESDIR):
        print FOLDER_ERROR
        sys.exit(1)
    bayesbin = os.path.join(BAYESDIR,'bin',BAYESBIN)
    ttimes = glob.glob(os.path.join(BAYESDIR,'ttimes','ak135.*'))
    if not os.path.isfile(bayesbin):
        print FOLDER_ERROR
        sys.exit(1)
    if not len(ttimes):
        print FOLDER_ERROR
        sys.exit(1)
    bayesdb = os.path.join(BAYESDIR,BAYESDB)
    # if startOver and os.path.isfile(bayesdb):
    #     os.remove(bayesdb)
    #does the database exist - if not, create it
    if not os.path.isfile(bayesdb):
        db = sqlite3.connect(bayesdb)
        cursor = db.cursor()
        createTables(db,cursor)
    else:
        db = sqlite3.connect(bayesdb)
        cursor = db.cursor()

    #Delete selected list of events
    if args.delete:
        nevents = deleteEvents(db,cursor,args.delete)
        print '%i events deleted from the database.' % nevents
        sys.exit(0)
        
    #Return some stats about the current database
    if args.stats:
        nevents,nstations,narrivals = getStats(cursor)
        print 'Your database contains information about:'
        print '\t%i events' % nevents
        print '\t%i stations' % nstations
        print '\t%i picks' % narrivals
        sys.exit(0)
        
    eventinfo = getPhaseData(eventid=eventid)
    if not len(eventinfo):
        print 'Could not find event %s in ComCat.  Returning.'
        sys.exit(1)

    #get the information about the input event
    eventinfo = eventinfo[0]
    eventlat = eventinfo.origins[0]['lat']
    eventlon = eventinfo.origins[0]['lon']
    eventtime = eventinfo.origins[0]['time']
    if eventtime < args.begindate or eventtime > args.enddate:
        fmt = 'Event %s (%s) is outside the time bounds you specified. %s to %s.  Exiting.' 
        print fmt % (eventinfo.eventcode,eventtime,args.begindate,args.enddate)
        sys.exit(1)

    tnow = datetime.utcnow()
    eventfolder = os.path.join(BAYESDIR,'events',eventid)
    if not os.path.isdir(eventfolder):
        os.makedirs(eventfolder)

    eventlist1 = getEventData(radius=(eventlat,eventlon,0,radius),
                             starttime=args.begindate,
                             endtime=args.enddate,catalog='pde')
    eventlist2 = getEventData(radius=(eventlat,eventlon,0,radius),
                              starttime=args.begindate,
                              endtime=args.enddate,catalog='us')
    eventlist = eventlist1 + eventlist2

    if args.count:
        fmt = 'There are %i events inside %.1f km radius around event %s (%.4f,%.4f)'
        print fmt % (len(eventlist),radius,eventid,eventlat,eventlon)
        sys.exit(0)
    
    #check to see if event has already been located - if so, stop, unless we're being forced
    if not args.force:
        sql = 'SELECT id,code,rlat,rlon,rdepth,rtime FROM event WHERE code="%s"' % eventid
        cursor.execute(sql)
        row = cursor.fetchone()
        if row is not None and row[2] is not None:
            print 'Event %s is already in the database.  Stopping.' % eventid
            sys.exit(0)
    
    priors = getEventPriors(eventlist,cursor)
    stations,arrivals,newevents = getProcessedData(eventlist,db,cursor,ndays=NWEEKS*7)
    fmt = 'In database: %i stations, %i arrivals.  %i events not in db.'
    #print fmt % (len(stations),len(arrivals),len(newevents))
    missing_stations = []
    for event in newevents:
        phasedata = getPhaseData(eventid=event)
        if phasedata is None:
            continue
        if not len(phasedata[0].magnitudes):
            continue
        newstations,newarrivals,ms = insertPhaseData(phasedata[0],db,cursor)
        stations = dict(stations.items() + newstations.items())
        arrivals += newarrivals
        missing_stations += ms

    print 'After searching online:'
    fmt = 'In database: %i stations, %i arrivals.  %i missing stations.'
    print fmt % (len(stations),len(arrivals),len(missing_stations))
    stafile = 'station.dat'
    stationfile = os.path.join(eventfolder,stafile)
    f = open(stationfile,'wt')
    f.write('sta_id lat lon elev\n')
    for stationcode,stationvals in stations.iteritems():
        slat,slon,elev = stationvals
        f.write('%s %.4f %.4f %.3f\n' % (stationcode,slat,slon,elev))
    f.close()

    arrfile = 'arrival.dat'
    arrivalfile = os.path.join(eventfolder,arrfile)
    f = open(arrivalfile,'wt')
    f.write('ev_id sta_id phase time\n')
    for arrival in arrivals:
        eid,scode,phase,time = arrival
        f.write('%i %s %s %.3f\n' % (eid,scode,phase,time))
    f.close()

    prifile = 'prior.dat' #??
    priorfile = os.path.join(eventfolder,prifile)
    f = open(priorfile,'wt')
    f.write('ev_id lat_mean lon_mean dist_sd depth_mean depth_sd time_mean time_sd\n')
    for prior in priors:
        evid,plat,plon,pdepth,ptime = prior
        f.write('%i %.4f %.4f 0.0 %.1f 0.0 %.3f 0.0\n' % (evid,plat,plon,pdepth,ptime))
    f.close()

    #write the config file
    configfile = os.path.join(eventfolder,'bayesloc.cfg')
    config = CONFIG.replace('BAYESLOC',BAYESLOC)
    config = config.replace('EVENTFOLDER',eventfolder)
    fcfg = open(configfile,'wt')
    fcfg.write(config)
    fcfg.close()

    #Run the BayesLoc program
    #change to the eventfolder
    cwd = os.getcwd()
    os.chdir(eventfolder)
    bayesbin = os.path.join(BAYESLOC,'bin','bayesloc')
    cmd = '%s %s' % (bayesbin,configfile)
    print 'Running command %s...' % cmd
    t1 = datetime.now()
    # process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    # for c in iter(lambda: process.stdout.read(1), ''):
    #     sys.stderr.write(c)
    res,stdout,stderr = getCommandOutput(cmd)
    t2 = datetime.now()
    if not res:
        print 'BayesLoc command "%s" failed.  \n%s\n%s.' % (cmd,stdout,stderr)
        sys.exit(1)
    else:
        dt = ((t2-t1).seconds)/60.0
        print 'BayesLoc command was successful - took %.1f minutes.' % dt
    os.chdir(cwd)

    resultfile = os.path.join(eventfolder,'output','origins_ned_stats.out')
    f = open(resultfile,'rt')
    f.readline()
    eventlist = []
    fieldlist = ['lat','lon','depth','time','rlat','rlon','rdepth','rtime','mag']
    for line in f.readlines():
        parts = line.split()
        eid = int(parts[0])
        lat = float(parts[1])
        lon = float(parts[2])
        depth = float(parts[3])
        time = UTCDateTime(float(parts[4])).datetime
        efmt = 'UPDATE event set rlat=%.4f,rlon=%.4f,rdepth=%.1f,rtime="%s" WHERE id=%i'
        equery = efmt % (lat,lon,depth,time,eid)
        cursor.execute(equery)
        db.commit()
        query = 'SELECT %s FROM event WHERE id=%i' % (','.join(fieldlist),eid)
        cursor.execute(query)
        row = cursor.fetchone()
        eventlist.append(dict(zip(fieldlist,row)))
    f.close()

    #make a map of all the relocated events
    fname = makeMap(eventlist,eventlat,eventlon,eventfolder)
    print 'Relocated events: %s' % fname
    
    #tell the user what happened with the relocation
    fmt = 'SELECT lat,lon,depth,time,rlat,rlon,rdepth,rtime FROM event WHERE code="%s"'
    query = fmt % (eventid)
    cursor.execute(query)
    row = cursor.fetchone()
    lat,lon,depth,time,rlat,rlon,rdepth,rtime = row
    time = UTCDateTime(time).datetime
    rtime = UTCDateTime(rtime).datetime
    if rtime >= time:
        dt = (rtime-time).seconds + ((rtime-time).microseconds)/float(1e6)
    else:
        dt = (time-rtime).seconds + ((time-rtime).microseconds)/float(1e6)
    dd,az1,az2 = gps2DistAzimuth(lat,lon,rlat,rlon)
    dd /= 1000.0
    print 'Event moved from:'
    print '%s (%.4f,%.4f) %.1f km' % (time.strftime('%Y-%m-%d %H:%M:%S'),lat,lon,depth)
    print '%s (%.4f,%.4f) %.1f km' % (rtime.strftime('%Y-%m-%d %H:%M:%S'),rlat,rlon,rdepth)
    print '%.1f km (%.1f degrees), %.1f seconds' % (dd,az1,dt)
    cursor.close()
    db.close()

if __name__ == '__main__':
    desc = '''Use BayesLoc to help automate the creation of a relocated earthquake catalog.
    ''' 
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i','--id', 
                        help='The NEIC event ID to relocate')
    parser.add_argument('-r','--radius', default=DEFAULT_RADIUS,type=float,
                        help='The radial search distance.')
    parser.add_argument('-s','--stats', action='store_true',
                        help='Display database statistics (number of events, stations, etc.)')
    parser.add_argument('-c','--count', action='store_true',
                        help='Display number of events in search radius.')
    parser.add_argument('-f','--force', action='store_true',
                        help='Force relocation of an already relocated event.')
    parser.add_argument('-d','--delete', nargs='*',
                        help='Delete event(s) from database.')
    parser.add_argument('-b','--begindate', default=DEFAULT_START,type=maketime,
                        help='Start time for search (defaults to %s).  YYYY-mm-dd or YYYY-mm-ddTHH:MM:SS' % DEFAULT_START)
    parser.add_argument('-e','--enddate', default=DEFAULT_END,type=maketime,
                        help='End time for search (defaults to now).  YYYY-mm-dd or YYYY-mm-ddTHH:MM:SS')
    
    pargs = parser.parse_args()
    main(pargs)
    

