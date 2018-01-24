import process_dropping

jobid = 1
landscapeFile = 'wawhab800.tif'
stFile = 'wawst800.tif'
restorationFile = 'wawzoo800.tif'
R= 1000.0
dispersal = 5.0
outputFile = 'waw_dropped.tif'
outputdata = 'waw_dropped.csv'
previewOnly = False
flowOnly = False
jobid  = 0

success = process_dropping.doDropping(landscapeFile = landscapeFile,    stFile = stFile,    restorationFile = restorationFile,    R = R,   dispersal = dispersal,    outputFile =outputFile , outputdata=outputdata, connection=None, previewOnly= previewOnly , flowOnly =flowOnly, jobid=jobid )
				
				