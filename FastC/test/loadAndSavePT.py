import Fast.Utils as FUtils
t,tc,NP,pDict,graphID,graphIBCD = FUtils.loadData(split='single',filedir='.')
FUtils.saveData(t,split='multi',filedir='DATA')
t,tc,NP,pDict,graphID,graphIBCD = FUtils.loadData(split='multi',filedir='.')
