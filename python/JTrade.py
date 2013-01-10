'''
Backtesting/Trading Classes
'''

import numpy as np
from datetime import datetime, date, time
import logging
#import scipy.io as sio

logging.basicConfig(\
        filename='log_1.log', \
        level=logging.DEBUG)

DT_FMT = "%Y-%m-%d %H:%M:%S"

class ContractType:  
    STOCK  = 1
    OPTION = 2
    FUTURE = 3
    FOREX  = 4    

class Contract:
    iID = np.nan
    sID = ''
    contractType = np.nan

    def __init__(self, iID, sID, contractType):
        self.iID = iID
        self.sID = sID
        self.contractType = contractType

    def toStr(self):
        return "[Contract]" \
            + " iID = " + str(self.iID) + ", sID = " + self.sID \
            + ", contractType = " + str(self.contractType)
    
    def printObj(self):
        print self.toStr()
 
class Quote:
    quoteTime = np.nan
    bid = np.nan
    ask = np.nan
    last = np.nan
    bidSize = np.nan
    askSize = np.nan

    def __init__(self, quoteTime, bid, ask, last, bidSize, askSize):
        self.quoteTime = quoteTime
        self.bid = bid
        self.ask = ask
        self.last = last
        self.bidSize = bidSize
        self.askSize = askSize
    
    def toStr(self):
        return "[Quote]" \
            + " quoteTime = " + self.quoteTime.strftime(DT_FMT) \
            + ", bid = " + str(self.bid) \
            + ", ask = " + str(self.ask) \
            + ", last = " + str(self.last) \
            + ", bidSize = " + str(self.bidSize) \
            + ", askSize = " + str(self.askSize)

    def printObj(self):
        print self.toStr();
        
class QuoteSource:
    def __init__(self):
        #print "QuoteSource Init"
        pass
    def connect(self):
        raise Exception("NotImplementedException")
    def disconnect(self):
        raise Exception("NotImplementedException")
    def getQuote(self, exchangeTime, contract):
        raise Exception("NotImplementedException")
    def getBars(self, tBeg, tEnd):
        raise Exception("NotImplementedException")

class QuoteSourceCSI(QuoteSource):

    ddir = "/Users/jimmiegoode/Documents/Emery/DATA_CSI/"
    
    def __init__(self):
        QuoteSource.__init__(self)
        logging.info("QuoteSourceCSI Init")
    def connect(self):
        logging.info("Connecting...")
    def disconnect(self):
        logging.info("Disconnecting...")
    def getQuote(self, exchangeTime, contract):
        logging.info("Getting quote...")
    def getBars(self, tBeg, tEnd):
        logging.info("Getting bars...")

class OrderTypes:
    STOP  = 1
    LIMIT = 2
    MOO   = 3
    MOC   = 4
        
class OrderStatus:
    WORKING  = 1
    FILLED   = 2
    CANCELED = 3
    REJECTED = 4

class OrderSides:
    BUY  = 1
    SELL = 2

class Order:
    oid_temp   = np.nan  # Temporary/Broker Order ID
    oid_perm   = np.nan  # Permanent Order ID
    oca_group  = ''      # OCA group name
    pos_group  = np.nan  # ID within a multi-leg position
    pos_id     = np.nan  # ID for position/round-trip
    strat_id   = ''      # ID for strategy
    side       = np.nan
    type       = np.nan
    status     = np.nan
    qty        = np.nan
    lmtPrice   = np.nan
    lmtPriceAux = np.nan

    def __init__(self):
        pass
    
