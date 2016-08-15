# Zip Supporting Information:

#  Supporting Info 1:

# Supporting Info 2 - Conversion Tables
zip(zipfile = 'S2_File.zip',
    files = c('data/input/relation_tables/fullpaleon_conversion_v0.4.csv',
              'data/input/relation_tables/pft.trans.table.csv',
              'data/input/relation_tables/FIA_conversion_v0.2.csv',
              'data/input/relation_tables/plss.pft.conversion_v0.1-1.csv'),
    flags = "-u")

# Supporting Info 3 - Output Files:
zip(zipfile = 'S3_File.zip',
    files = c('data/output/wiki_outputs/*.*'), flags = '-uD')

# Supporting Info 3 - Correction Factors

zip(zipfile = "S4_File.zip",
    files = 'data/input/relation_tables/cogbill_corrections.csv',
    flags = '-u')
