# interrogating data results

# check results
files = list.files('derived/', pattern = '\\.csv', full.names = T)
df = lapply(files, function(x) read.csv(x))
names(df) <- basename(files)
lapply(df, head)

KO_SALIvsSALI = read.csv("derived/201125_KO_SALIvsSALI.csv", row.names = NULL)


