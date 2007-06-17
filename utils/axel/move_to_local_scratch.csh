du -sm *
mkdatadir_scratch
gt data1


move-VAR-to-data1


mv data data1
mkdatadir_scratch
mv data1/* data
du -sm *
