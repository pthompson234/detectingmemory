source('sim_functions_20200326.R')

for (ID in c('GF1004', 'GM1046')) {
  # for bears with four total years we cut into year 2 (with year 1 as training) and year 4 (with year 3 as training)
  data.raw = read.SSF.raw.data(paste0('inputdata/', ID))
  data.cut.1 = cut.raw.data(data.raw, 0, 1/3)
  print(min(data.cut.1$strata.controls))
  print(max(data.cut.1$strata.controls))
  data.cut.2 = cut.raw.data(data.raw, 2/3, 1, train.begin = max(data.cut.1$strata.controls))
  print(min(data.cut.2$strata.controls))
  print(max(data.cut.2$strata.controls))
  write.SSF.raw.data(data.cut.1, paste0('inputdata_subset/', ID, '_set1'))
  write.SSF.raw.data(data.cut.2, paste0('inputdata_subset/', ID, '_set2'))
}

for (ID in c('GF1008', 'GF1016', 'GF1041', 'GF1086', 'GF1107', 'GF1130')) {
  # for bears with three total years we cut into year 3 (with year 2 as training) and year 2 (with year 1 as training)
  data.raw = read.SSF.raw.data(paste0('inputdata/', ID))
  data.cut.1 = cut.raw.data(data.raw, 0, 1/2)
  print(min(data.cut.1$strata.controls))
  print(max(data.cut.1$strata.controls))
  data.cut.2 = cut.raw.data(data.raw, 1/2, 1, train.begin = min(data.cut.1$strata.controls))
  print(min(data.cut.2$strata.controls))
  print(max(data.cut.2$strata.controls))
  write.SSF.raw.data(data.cut.1, paste0('inputdata_subset/', ID, '_set1'))
  write.SSF.raw.data(data.cut.2, paste0('inputdata_subset/', ID, '_set2'))
}
