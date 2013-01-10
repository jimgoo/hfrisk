

# # binName = 'stockquotes.bin'
# # fout = open(binName, 'wb')
# # mnStk.tofile(fout)
# # fout.close()

# # m = np.memmap(binName, dtype=type_stks, mode='r')



## These should all be added as per Java code.

# attrs.put("midprice", midPrice); #ok
# attrs.put("voldoll", voldoll); #ok
# attrs.put("expnum", expnum); <TODO in chains>
# attrs.put("strikenum", strikeNum); <TODO in chains>
# attrs.put("opttype", optTypeStr); #ok

# attrs.put("mo", mo); #ok
# attrs.put("expYear", expYear); #ok
# attrs.put("dayOfYear", dayOfYear); #ok
# attrs.put("dayOfMonth", dayOfMonth); #ok

# attrs.put("days", days); #ok
# attrs.put("moneyness", moneyness); #ok
# attrs.put("itm", itm);  #ok
# attrs.put("fri3", fri3); #ok
# attrs.put("thv", thv); <TODO>
# attrs.put("pou", pou); <TODO>


        map = ColumnMapper('options_DN')
        map.initFromObj(self)
        map.updateType('aka', ColumnMapper.STRING)
        map.updateType('exchange', ColumnMapper.STRING)
        map.updateType('optext', ColumnMapper.STRING)
        map.updateType('symo', ColumnMapper.STRING)
        map.updateType('symu', ColumnMapper.STRING)

        map.updateType('dataDate', ColumnMapper.INT)
        map.updateType('dataTime', ColumnMapper.INT)
        map.updateType('days', ColumnMapper.INT)
        map.updateType('expDate', ColumnMapper.INT)
        map.updateType('expTime', ColumnMapper.INT)
        map.updateType('expDay', ColumnMapper.INT)
        map.updateType('expMonth', ColumnMapper.INT)
        map.updateType('expYear', ColumnMapper.INT)
        map.updateType('oi', ColumnMapper.INT)
        map.updateType('opttype', ColumnMapper.INT)
        map.updateType('vol', ColumnMapper.INT)

        map.updateType('dataDT', ColumnMapper.DATE)
        map.updateType('expDT', ColumnMapper.INT)

        map.updateType('fri3', ColumnMapper.BOOL)
        map.updateType('itm', ColumnMapper.BOOL)

            # def load(self, mapName):
    #     self.dictMap = {}
    #     csvMap = np.genfromtxt(self.COLMAP_DIR + "/colmap_" + mapName + ".csv",
    #                            delimiter=",",
    #                            dtype=[('col','S10'),('idx','i')])
    #     for rec in csvMap:
    #         self.dictMap[rec[0]] = rec[1]

    # def save(self, mapName):
    #     fout = file(self.COLMAP_DIR + "/colmap_" + mapName + ".csv", "w+")
    #     for key in self.dictMap.keys():
    #         fout.write(key + "," + str(self.dictMap[key]) + "\n")
    #     fout.close()
    
