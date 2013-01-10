import CG2
import timeit
import cProfile
import numpy as np
import time
import h5py

runID = 'initial_load'
symu = 'A'
source = 'DN'

d = CG2.DeltaLoader()

days = ['20121025', '20121026']

iOpt = 0

for sDay in days:
    iOpt += d.shardOptionFile(sDay, runID)

colMap = CG2.Option.getColMap()
colMapH5 = colMap.makeH5pyDesc()

# hd = CG2.HDF5()
# hd.loadOptions(source, symu, opts, colMapH5)

[opts,arr] = d.loadOptionShard(symu, runID, colMapH5)

#================================================================================
# Update routine
#================================================================================

sDay = '20121031'
runID = sDay # runID is sDay
d.shardOptionFile(sDay, runID) 
[opts,arr] = d.loadOptionShard(symu, runID, colMapH5) # shard is just options for one day


#================================================================================
# check option counts

i = 0
for sDay in ['20121025', '20121026', '20121031']:
    csv = d.parseOptionFile(d.getRawCsvFilename(sDay, 'options'))
    idx = np.nonzero(csv['symu'] == symu)
    i += len(idx[0])
    print i

#================================================================================
# query test
    
#np.nonzero(f['DN']['opttype'][:] == 1 and f['DN']['vol'][:] > 100)


#d.shardOptionFile(sDay, runID)
#sfile = d.getShardFilename('AAPL', runID)
#optObjs = d.convertOptionParse(optCsv, -1)


# # #cProfile.run('d.convertOptionParse(optCsv)')
# # chains = CG2.OptionChainUtil.makeOptionChains(optObjs)
# # opt = optObjs[0]

# # #colMapObj = opt.initMap()
# # colMapObj = CG2.ColumnMapper('options_DN')
# # colMapObj.load()

# # colMap = colMapObj.dictMap
# # cols = sorted(colMap, key=lambda key: colMap[key][1]) # sort by index
# # print cols

# reload(CG2)
# map = CG2.ColumnMapper()
# obj = optObjs[0]
# map.initFromObj(obj)
# map.save("options_dn")



#==================== PyTables =======================

# chain = chains[0]
# chainList = []
# for opt in chain.options:
#     chainList.append(opt.toList(cols))
# tblDesc = colMapObj.makePytableDesc(cols)
# import tables
# h5 = tables.openFile('db1.h5', mode='a', title='Describe file here.')
# # group for symu
# group = h5.createGroup("/", chain.symu, 'Describe group here.')
# # table for symu's options
# table = h5.createTable(group, 'options_DN_1', tblDesc, 'Describe table here.')
# table.append(chainList)
# table.flush()  # flush data in the table
# h5.flush()  # flush all pending data
# h5.close()



#
# type_stk  = [('symu','S10'),('date','S12'),('op','f'),
#                  ('hi','f'),('lo','f'),('cl','f'),('vol','i')]

# stk = np.genfromtxt("../dn/extracted/stockquotes_" + str(tBeg) + ".csv", delimiter=",", skip_header=0, dtype=type_stk)

# ds = optCsv['expiration']
# dt = datetime.strptime(ds, '%m/%d/%Y')

# for c in chains:
#     if c.symu == 'AAPL':
#         break
