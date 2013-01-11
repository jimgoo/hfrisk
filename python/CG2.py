'''
#===========================================================================
# Cladogenesis 2.0
# By Jimmie Goode
#===========================================================================
# 2012-12-12: Created by JG
#===========================================================================
'''

# general
import numpy as np
from datetime import datetime, date, time, timedelta
import pickle, logging, os

#HDF5
import tables, h5py

# for third fridays
import dateutil.rrule as dr
import dateutil.relativedelta as drel

# for printing objects
#from pprint import pprint

# for grouping list of objects
#from collections import defaultdict

import time as time2


'''
Class for displaying output to user and logging.
'''
class MSG:

    ## Setup the logger
    #logging.basicConfig(filename='log_1.log', level=logging.INFO)

    @staticmethod
    def info(str):
        print str

'''
Class for IO helping
'''

class IO:

    ''' Load a dictionary from pickle file. '''
    @staticmethod
    def loadDict(fileName):
        with open(fileName, 'r') as handle:
            return pickle.loads(handle.read())

    ''' Save a dictionary to pickle file. '''
    @staticmethod
    def saveDict(myDict, fileName):
        with open(fileName, 'w+') as handle:
            pickle.dump(myDict, handle)


'''
Dates and Times (DnT) utility class.
'''
class DnT:

    # directory of exchange calendar files
    CALENDAR_DIR = os.getenv("CG_CALENDARS", "../calendars")
    
    def __init__(self):
        pass
    
    '''
    Get NYSE business days within range [tBeg,tEnd].
    Input dates tBeg and tEnd are numeric yyyymmdd format.
    '''
    @staticmethod
    def nyseDates(tBeg, tEnd):
        vnBD = np.genfromtxt(DnT.CALENDAR_DIR + "/dates_nyse.csv", delimiter=',', dtype="int")
        vnBD = vnBD[np.nonzero(np.logical_and(vnBD >= tBeg, vnBD <= tEnd))[0]]
        return vnBD

    ''' Input date is Python datetime. '''
    @staticmethod
    def isThirdFriday(nDate):
        rr = dr.rrule(dr.MONTHLY, byweekday=drel.FR(3), dtstart=nDate, count=1)
        if rr[0] == nDate:
            return True
        else:
            return False

    ''' Split datetime into integer date and times. '''
    @staticmethod
    def datetimeSplit(dt):
        #dtStr = dt.strftime('%Y%m%d%H%M%S%f')
        return [int(dt.strftime('%Y%m%d')), int(dt.strftime('%H%M%S'))]
        

class ColumnMapper:

    # directory of column map files
    COLMAP_DIR = os.getenv("CG_COLMAP", "../data/colmaps")
    
    dictMap = {}
    mapName = ''
    fileName = ''

    INT      = 'INT'
    FLOAT    = 'FLOAT'
    DATE     = 'DATE'
    STRING   = 'STRING'
    STRING30 = 'STRING30'
    BOOL     = 'BOOL'

    def __init__(self, mapName):
        self.dictMap = {}
        self.mapName = mapName
        self.fileName = self.COLMAP_DIR + "/colmap_" + self.mapName + ".pickle"
    
    def load(self):
        self.dictMap = {}
        #with open(self.fileName, 'r') as handle:
        #    self.dictMap = pickle.loads(handle.read())
        IO.loadDict(self.fileName)

    def save(self):
        # with open(self.fileName, 'w+') as handle:
        #     pickle.dump(self.dictMap, handle)
        IO.saveDict(self.dictMap, self.fileName)
        
    def initFromObj(self, obj):
        self.dictMap = {}
        keys = vars(obj).keys()
        keys.sort()
        idx = 0
        for key in keys:
            self.dictMap[key] = [self.FLOAT, idx] # [type,index]
            idx = idx + 1

    def addColumn(self, colName, colType):

        # check if column is already there
        if colName in self.dictMap.keys():
            raise Exception("Column is already present in column index map.")

        if len(self.dictMap.keys()) == 0:
            self.dictMap[colName] = [colType, 0]
            return
        
        # get max column index
        maxIdx = 0
        for key in self.dictMap.keys():
            maxIdx = max(maxIdx, self.dictMap[key][1])

        # add new column with incremented index
        self.dictMap[colName] = [colType, maxIdx + 1]

    def updateType(self, colName, newColType):
        self.dictMap[colName][0] = newColType
    
    def makePytableDesc(self):
        desc = {}
        idx = 1

        # sort by index
        cols = sorted(self.dictMap, key=lambda key: self.dictMap[key][1])
        
        for col in cols:
            t = self.dictMap[col][0]    
            if t == self.INT:
                desc[col] = tables.IntCol(pos=idx)
            elif t == self.FLOAT:
                desc[col] = tables.FloatCol(pos=idx)
            elif t == self.DATE:
                desc[col] = tables.IntCol(pos=idx)
            elif t == self.STRING:
                desc[col] = tables.StringCol(10, pos=idx)
            elif t == self.STRING30:
                desc[col] = tables.StringCol(30, pos=idx)
            elif t == self.BOOL:
                desc[col] = tables.BoolCol(pos=idx)
            else:
                raise Excpetion("Unknown type: " + t + ", for column: " + col)
            idx += 1
        return desc

    def makeH5pyDesc(self):
        desc = {}

        # sort by index
        cols = sorted(self.dictMap, key=lambda key: self.dictMap[key][1])
        
        for col in cols:
            t = self.dictMap[col][0]    
            if t == self.INT:
                desc[col] = 'i4'
            elif t == self.FLOAT:
                desc[col] = 'f4'
            elif t == self.DATE:
                raise Exception("Datetime type not supported in HDF5.")
            elif t == self.STRING:
                desc[col] = 'S10'
            elif t == self.STRING30:
                desc[col] = 'S30'
            elif t == self.BOOL:
                desc[col] = 'b'
            else:
                raise Excpetion("Unknown type: " + t + ", for column: " + col)
        return desc

'''
Enumeration of option contract types.
'''
class OptionType:  
    CALL = 1
    PUT  = 2

    @staticmethod
    def toStr(opttype):
        if opttype == 1:
            return "CALL"
        elif opttype == 2:
            return "PUT"
        else:
            return "NULL"

'''
Container for option properties.
'''
class Option:
    symu = ''
    lastu = np.nan
    exchange = ''
    symo = ''
    optext = ''
    opttype = np.nan
    #expiration = np.nan
    #date = np.nan
    strike = np.nan
    last = np.nan
    bid = np.nan
    ask = np.nan
    #midPrice = np.nan
    vol = np.nan
    #voldoll = np.nan
    oi = np.nan
    iv = np.nan
    delta = np.nan
    gamma = np.nan
    theta = np.nan
    vega = np.nan
    aka = ''
    
    #moneyness = np.nan
    days = np.nan
    itm = np.nan
    fri3 = np.nan 
    thv = np.nan  # <TODO>

    # expiry date
    #expDT = np.nan
    expDate = np.nan
    expTime = np.nan
    expYear = np.nan
    expMonth = np.nan
    expDay = np.nan

    # data observation date
    #dataDT = np.nan
    dataDate = np.nan
    dataTime = np.nan
    
    def __init__(self):
        pass

    @staticmethod
    def getColMap():

        m = ColumnMapper('options_DN')

        m.addColumn('symu', ColumnMapper.STRING)
        m.addColumn('lastu', ColumnMapper.FLOAT)
        m.addColumn('exchange', ColumnMapper.STRING)
        m.addColumn('symo', ColumnMapper.STRING30)
        m.addColumn('optext', ColumnMapper.STRING)
        m.addColumn('opttype', ColumnMapper.INT)
        m.addColumn('strike', ColumnMapper.FLOAT)
        m.addColumn('last', ColumnMapper.FLOAT)
        m.addColumn('bid', ColumnMapper.FLOAT)
        m.addColumn('ask', ColumnMapper.FLOAT)
        m.addColumn('vol', ColumnMapper.INT)
        m.addColumn('oi', ColumnMapper.INT)
        m.addColumn('delta', ColumnMapper.FLOAT)
        m.addColumn('gamma', ColumnMapper.FLOAT)
        m.addColumn('theta', ColumnMapper.FLOAT)
        m.addColumn('vega', ColumnMapper.FLOAT)
        m.addColumn('aka', ColumnMapper.STRING)

        m.addColumn('days', ColumnMapper.INT)
        m.addColumn('itm', ColumnMapper.BOOL)
        m.addColumn('fri3', ColumnMapper.BOOL)
        m.addColumn('thv', ColumnMapper.FLOAT)

        #m.addColumn('expDT', ColumnMapper.DATE)
        m.addColumn('expDate', ColumnMapper.INT)
        m.addColumn('expTime', ColumnMapper.INT)
        m.addColumn('expYear', ColumnMapper.INT)
        m.addColumn('expMonth', ColumnMapper.INT)
        m.addColumn('expDay', ColumnMapper.INT)

        #m.addColumn('dataDT', ColumnMapper.DATE)
        m.addColumn('dataDate', ColumnMapper.INT)
        m.addColumn('dataTime', ColumnMapper.INT)

        return m

        
    def update_itm(self):
        self.itm = np.nan
        if self.opttype == OptionType.CALL:
            if self.lastu >= self.strike:
                self.itm = True
            else:
                self.itm = False
        else:
            if self.lastu <= self.strike:
                self.itm = True
            else:
                self.itm = False

    ''' Create option from Delta Netural CSV record '''
    def init_DN_csv(self, rec):
        
        self.symu = rec['symu']
        self.lastu = rec['lastu']
        self.exchange = rec['exchange']

        symo = rec['symo']
        if len(symo) > 30:
            raise Exception("Symo will be truncated.")
        self.symo = symo
        
        self.optext = rec['optext']

        ## get option type
        if rec['opttype'] == 'call':
            self.opttype = OptionType.CALL
        elif rec['opttype'] == 'put':
            self.opttype = OptionType.PUT
        else:
            raise Exception("Unknown option type: " + rec['opttype'])

        ## format expiration date, ignoring DN provided time " 04:00:00 PM"
        dt = datetime.strptime(rec['expiration'].split(' ')[0], '%m/%d/%Y')
        
        ## set saturday expiry dates to previous Friday
        if dt.weekday() == 5:
            dt = dt - timedelta(days=1)

        ## add time to yield the exchange closing time (4:00 EST)
        expDT = dt + timedelta(hours=16, minutes=0)
        [self.expDate, self.expTime] = DnT.datetimeSplit(dt)
        self.expYear = dt.year
        self.expMonth = dt.month
        self.expDay = dt.day

        ## format observation date, ignoring DN provided time " 04:00:00 PM"
        dt = datetime.strptime(rec['date'].split(' ')[0], '%m/%d/%Y')
        dataDT = dt + timedelta(hours=16, minutes=0)
        [self.dataDate,self.dataTime] = DnT.datetimeSplit(dt)

        self.strike = rec['strike']
        self.last = rec['last']
        self.bid = rec['bid']
        self.ask = rec['ask']
        #self.midPrice = (self.bid + self.ask)/2.0
        self.vol = rec['vol']
        #self.voldoll = self.last * self.vol * 100.0 
        self.oi = rec['oi']
        self.iv = rec['iv']
        self.delta = rec['delta']
        self.gamma = rec['gamma'] / 100.0
        self.theta = rec['theta'] / 100.0
        self.vega = rec['vega'] / 100.0
        self.aka = rec['aka']
        #self.moneyness = 100.0*(self.strike - self.lastu)/self.lastu

        ## days to maturity
        if dataDT < expDT:
            self.days = (expDT - dataDT).days
        else:
            self.days = 0

        # set the ITM indicator
        self.update_itm()

        # Set third friday indicator
        self.fri3 = DnT.isThirdFriday(expDT)


    '''
    Return a list of the fields in "cols": (col_1, col_2, ...)
        Use the sorted output of the initMap() method:
            cols = sorted(mydict, key=lambda key: mydict[key][1])
    '''
    def toList(self, cols):
        s = '('
        for col in cols:
            s += "self." + col + ","
        s = s[0:-1] + ')'
        return eval(s)

# class OptionChain:
#     symu = ''
#     options = []
    
#     def __init__(self, symu, options):
#         self.symu = symu
#         self.options = options

#     # def assignExpiryNumbers(self):
#     #     # <TODO>

# class OptionChainUtil:

#     ''' Group individual options into OptionChains '''
#     @staticmethod
#     def makeOptionChains(optObjs):
        
#         # group options by underlier
#         groups = defaultdict(list)
#         for obj in optObjs:
#             groups[obj.symu].append(obj)
#         grouped = groups.values()

#         chains = []
#         for grp in grouped:
#             chain = OptionChain(grp[0].symu, grp)
#             chains.append(chain)
#         return chains

"""
This serves as the general database class.
"""
class HDF5:

    DATA_DIR = os.getenv("CG_HDF", "../HDF")
    
    def __init__(self):
        pass

    def makeDir():
        if os.path.exists(self.DATA_DIR) == False:
            os.mkdir(self.DATA_DIR)

    def makeOptionFilename(self, source, symu):
        sDir = self.DATA_DIR + "/" + source
        if os.path.exists(sDir) == False:
            os.mkdir(sDir)
        return sDir + "/opt_" + symu + ".hdf5"

    def getOptionHDF(self, source, symu):
        fname = self.makeOptionFilename(source, symu)
        hdfFile = h5py.File(fname, 'a')
        return hdfFile

    '''
    Add list of option objects for one underlier to it's HDF5 file.
    '''
    def loadOptions(self, source, symu, optObjs, colMap):
        
        #fname = self.makeOptionFilename(source, symu)
        #hdfFile = h5py.File(fname, 'a')
        hdfFile = self.getOptionHDF(source, symu)

        try:
            if source not in hdfFile.keys():
                grp = hdfFile.create_group(source)
            else:
                grp = hdfFile.get("/" + source)

            iNew = len(optObjs) # number of new options

            newData = {} # empty dictionary for new option data

            # initiliaze empty arrays for new data
            for col in colMap:
                sType = colMap[col]
                newData[col] = np.empty(iNew, dtype=sType)

            # fill arrays with new data
            idx = 0
            for opt in optObjs:
                if symu != getattr(opt, 'symu'):
                    raise Exception("Mismatching underlier symbols.")
                for col in colMap:
                    newData[col][idx] = getattr(opt, col)
                idx += 1
                
            # Create and/or expand datasets
            for col in colMap:
                sType = colMap[col]

                if col not in grp.keys():
                    # If no column with this name exists, make a new dataset of size iNew.
                    # In this case, newData is very large in size.
                    dset = grp.create_dataset(col, (iNew,), maxshape=(None,),
                                              dtype=sType, data=newData[col])
                else:
                    # If column exists, resize it and append new data.
                    # In this case, newData should be small (ie. when updating DB with new data).
                    dset = grp.get(col)
                    [iRows,] = dset.shape
                    
                    # check that there is only one column
                    #if iCols != 1:
                    #    raise Exception("Dataset should only have one column but has " + str(iCols))

                    # resize
                    iNewSize = iRows + iNew
                    dset.resize((iNewSize,))
                    
                    # set new elements (this should not bring another copy of newData[col] into memory)
                    dset[iRows:iNewSize] = newData[col]
                    
        finally:
            hdfFile.close()

    ''' Verify that all fields have the same length. '''
    def checkOptionLengths(self, source, symu):
        hdfFile = self.getOptionHDF(source, symu)
        try:
            grp = hdfFile[source]
            viLen = []
            for col in grp.keys():
                viLen.append(grp[col].shape[0])

            if len(np.unique(viLen)) != 1:
                return [False, viLen]
            else:
                return [True, viLen]
        finally:
            hdfFile.close()

    
            
class DeltaLoader:

    # directory of extracted data
    dn_dir = os.getenv("CG_DN", "../dn")

    # shorthand name of data source
    sourceName = "DN"
    
    # define CSV data types
    type_stk  = [('symu','S10'),('date','S12'),('op','f'),
                 ('hi','f'),('lo','f'),('cl','f'),('vol','i')]

    type_stat = [('symu','S10'),('date','S12'),('civ','f'),
             ('piv','f'),('meaniv','f'),('cvol','i'),
             ('pvol','i'),('coi','i'),('poi','i')]

    type_opt  = [('symu','S10'),('lastu','f'),('exchange','S10'),
             ('symo','S30'),('optext','S10'),('opttype','S10'),
             ('expiration','S30'),('date','S30'),
             ('strike','f'),('last','f'),('bid','f'),('ask','f'),
             ('vol','i'),('oi','i'),('iv','f'),
             ('delta','f'),('gamma','f'),('theta','f'),('vega','f'),
             ('aka','S30')]

    def getRawCsvFilename(self, sDay, sType):
        return self.dn_dir + "/extracted/triple_" + sDay + "/" + sType + "_" + sDay + ".csv"

    def getShardFilename(self, sSymu, runID):
        sDir = self.dn_dir + "/sharded/options/" + runID
        if os.path.exists(sDir) == False:
            os.mkdir(sDir)
        return sDir + "/EQ_" + sSymu + ".csv"
    
    def parseOptionFile(self, fname):
        
        MSG.info("Parsing option file: " + fname)
        tic = time2.time()
            
        opt = np.genfromtxt(fname, delimiter=",", skip_header=0, dtype=self.type_opt)

        toc = time2.time() - tic
        MSG.info("Elapsed time: " + str(toc))
        
        return opt

    def convertOptionParse(self, opts, iMaxRecs):
        # sort options by undelier
        # opts = np.sort(opts, order=['lastu'])

        MSG.info("Converting option parse to objects...")
        tic = time2.time()
        
        objs = []
        
        count = 0
        for rec in opts:
            opt = Option()
            opt.init_DN_csv(rec)
            objs.append(opt)
            count += 1
            if iMaxRecs > 0 and count > iMaxRecs:
                MSG.info("Stopped parsing due to iMaxRecs being exceeded.")
                break
            
        toc = time2.time() - tic
        MSG.info("Elapsed time: " + str(toc))
        
        return objs

    '''
    Split option CSV file into multiple files, one for each underlier.
    If shard files already exist, they are appended to. Note this can result in duplicates!
    '''
    def shardOptionFile(self, sDay, runID):

        tic = time2.time()
        lastSymu = ''
        count = 0
        
        with open(self.getRawCsvFilename(sDay, "options"), "r") as csvfile:
            for line in csvfile:
                idx = line.find(",")
                symu = line[0:idx]

                if symu != lastSymu:
                    if count != 0:
                        fid.close() # close existing open file

                    #print "opening file for " + symu
                    fid = open(self.getShardFilename(symu, runID), "a") # open symu file for appending

                lastSymu = symu
                count += 1
                fid.write(line)

        fid.close()
        toc = time2.time() - tic
        MSG.info("Parse count:  " + str(count))
        MSG.info("Elapsed time: " + str(toc))
        # 521,676 options
        # 3810 underliers

        return count

    def loadOptionShard(self, sSymu, runID, colMapH5):

        sfile = self.getShardFilename(sSymu, runID)
        
        arr = self.parseOptionFile(sfile)
        
        if len(np.unique(arr['symu'])) != 1:
            raise Exception("Multiple underlier symbols present in sharded options file.")

        if any(arr['symu'] != sSymu):
            raise Exception("Mismatching symbol in options file: " + sfile)

        opts = self.convertOptionParse(arr, -1)

        hdf5 = HDF5()
        hdf5.loadOptions(self.sourceName, sSymu, opts, colMapH5)

        return [opts,arr]

