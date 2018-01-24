#!/usr/bin/env python
import scipy as sp
import numpy as np
import landscape as ls
import condatisHDF as chdf
import condatisNE as cne
import condatiscore as cc
#import condatisplot as cp
import logging
import tempfile
import os
import math
import multiprocessing
import make_thumb
import traceback
import log_all
import serverio
    
# Previous versions allowed voltages to be scaled by a constant factor
scaleVoltage = 1

def createHabitat(h5,landscapeFile,STFile):
    if h5.root.scenarios.__contains__("root"):
        h5.removeNode("/scenarios/root",recursive=True)
    
    log_all.log("get landscape")    
    l=ls.Landscape.fromGIS(landscapeFile,makeone=False)
 
   
    l.scalev(scaleVoltage)
    maxVoltage = max(l.v)
	
	# Assume that where inputs are above 1, it must be a percentage
    if maxVoltage > 1.1:
        l.scalev(1.0/100.0)
        
    c=chdf.CondatisCoreHDF.fromLandscape(h5,'root',l)
    log_all.log("get st")    
    st=ls.Landscape.fromGIS(STFile)
    
    log_all.log("add source")
    c.addSource(st.level(1))
    log_all.log("add target")
    c.addTarget(st.level(2))
    log_all.log("gen comb")
    c.generateCombinedHabitat()
    return c
    
    # if passed, the connection object provides functions that allow communication with a server
	# previewOnly loads in the landscapes to give some info on the number of cells etc.
	# the jobid may be required to identify this instance to the server
	# thumbFile may be an image file for a thumbnail, which will be passed to the server as html

def drop(c,restorationLandscape, N=None, loopParamOne = .85, loopScale = 1, loopParam = 1, connection=None,outputdata="outdata.csv",previewOnly = False,jobid = 1, thumbFile = None):
    
    if not N:
        N=c.habitatL().len()

    if N == 0:
       raise ValueError("Error: empty habitat")
 
       
    # Make a copy of the input condatis object
    cl=c.habitatL()
    
    c2=cne.CondatisCoreNE(cl)
    c2.addSource(c.sourceL())
    c2.addTarget(c.targetL())
    c2.setParams(c.R(),c.dispersal())

    # Make an empty landscape object
    result=ls.Landscape()
    result.attribs=c.attribs

    result.firstFlow = None
    result.lastFlow = None
    restorationLandscape.scalev(scaleVoltage)

	# Assume that where inputs are above 1, it must be a percentage    
    maxVoltage = max(restorationLandscape.v)
    if maxVoltage > 1.1:
        restorationLandscape.scalev(1.0/100.0)
    
    restorationLandscape.removeDuplicates(c2.habitatL())
    c2.modifyHabitat(c2.habitatL().append(restorationLandscape))
    
    
    N=restorationLandscape.len();
    numCells = c2.x().size
    
    minx = min(c2.x())
    maxx = max(c2.x())
    miny = min(c2.y())
    maxy = max(c2.y())

    maxVal =  max(c2.habitatL().v)
    sumVal = sum(c2.habitatL().v)
    if (maxVal > 1.1):
       maxVal = maxVal / 100.0
       sumVal = sumVal / 100.0
    sumVal = sumVal / numCells

	# a very loose estimate of time to process
    estimated_time = int (2 + (numCells*numCells) * (N) / 4e9)
	
    timeScale = "minutes"
    
    if estimated_time > 120:
        estimated_time = estimated_time / 60
        timeScale = "hours"
        if estimated_time > 48:
            estimated_time = estimated_time / 24
            timeScale = "days"
    
    img_tag=""
	
	# html is created to make a nice 'preview' image
	# if present, a thumbnail file is encoded into the html
    if thumbFile is not None:
        data_uri = open(thumbFile, 'rb').read().encode('base64').replace('\n', '')
        img_tag = '<img src="data:image/png;base64,%s">' % data_uri


    previewText = "<table  cellpadding='6'><tr><td>Number of cells </td><td>  "+ str(N) + "</td> <td rowspan='8'>"\
    ""+img_tag+"</td></tr>" \
    "<tr><td> Source size </td><td>" + str(c2.sourceL().x.size)+" cells </td></tr>" \
    "<tr><td> Target size </td><td>" + str(c2.targetL().x.size)+" cells </td></tr>" \
    "<tr><td> Estimated time </td><td>" +str(estimated_time)+" "+timeScale+" </td></tr>" \
    "<tr><td> East-West Range </td><td>" + str(minx) +" - "+str(maxx)+ "</td></tr>" \
    "<tr><td> North-South Range </td><td>" +str(miny) + " - " + str(maxy) +  "</td></tr>" \
    "<tr><td> Greatest cell coverage </td><td>" + str(maxVal*100) +"% </td></tr>" \
    "<tr><td> Average cell coverage </td><td>" + str(sumVal * 100) +"% </td></tr></table>"

	
    if previewOnly:
        if connection is not None:
           if connection.previewReport is not None:
              connection.previewReport(previewText)
        return result
    
	
    restorationLandscapeSet = ls.Landscape.getKeySet(restorationLandscape)
    cpuIndex = multiprocessing.cpu_count()* 150

    if N > cpuIndex:
       loopParam = loopParam + (int) (math.sqrt(N-cpuIndex)/20)
       

    loopParam = loopParam * loopScale   + (1* (1-loopScale))
 
    outputDataFile = open(outputdata,'w')
    outputDataFile.write("i,Speed,x,y,flow\n")
    
    # Dropping loop
    i = 1
    for k in range(N-1):
      
        c2.calc()
        sp=c2.speed()
        # checks against the feasible habitat
        nodeFlow = c2.nodeFlowL()

		# storing the first flow for later reporting
        if i == 1:
            result.firstFlow = nodeFlow
        maxVoltage = max(c2.habitatL().v)
        if maxVoltage > 1.1:
            print "out of range"
            raise ValueError("Failed with an out of range voltage  "+str(maxVoltage))
        
        
        # Get (feasible) node with smallest flow
        ignoreSet=set()
        smallestList=[]
        
        if connection is not None:
            if connection.progressReport is not None:
                connection.progressReport((i*100.0)/N)

		# work out how many cells to drop in this iteration 
				
        loopSize = int (((N-i*1.0)/N) * loopParam ) + 1

		# don't drop more than we have
        if loopSize > N-i:
            loopSize = N -i + 1
        
		# e.g., only drop individually for most significant cells
        if i > (N*loopParamOne):   
            loopSize = 1


           
        for j in range(loopSize):
            smallest=nodeFlow.argminBoth(restorationLandscapeSet,ignoreSet)
            smallestList.append(smallest)

        # Create landscape object from the smallest flow
        # set value to i and append to result
            l=c2.habitatL()[smallest]
            ignoreSet.add( ls.Landscape.make_location_key(l.x,l.y) )
        
            l.v=i
            result=result.append(l)

        land=c2.habitatL()
        
        smallestList.sort(reverse=True)

        # Remove the nodes with the smallest flow

        for smallest in smallestList:
            l=c2.habitatL()[smallest]
            nodeFlow = c2.nodeFlow()[smallest]
            log_all.log("%i out of %i, Speed: %e  x,y = %f,%f, i= %e, %i" % (i,N,sp,l.x,l.y,c2.nodeFlow()[smallest],smallest))
            logging.info("%i out of %i, Speed: %e  x,y = %f,%f, i= %e, %i" % (i,N,sp,l.x,l.y,c2.nodeFlow()[smallest],smallest))
            outputDataFile.write("%i,%e,%f,%f,%e,%i\n" % (i,sp,l.x,l.y,nodeFlow,smallest))
            c2.modifyHabitat(land.delete(smallest))
            i = i + 1
            if connection is not None:
                if connection.progressReport is not None:
                    connection.progressReport((i*100.0)/N)

        if i > N:
            break

    outputDataFile.close()
    if connection is not None:
        if connection.progressReport is not None:
           connection.progressReport(100)

    c2.calc()            
    result.lastFlow = c2.nodeFlowL()  
    return result

    
    
    
    
    
    
    
    
# very similar code to dropping, but only computes the flow

    
def flow(c,restorationLandscape, N=None, loopParamOne = .85, loopScale = 1, loopParam = 1, connection=None,outputdata="outdata.csv",previewOnly = False,normalise = None, jobid=0, thumbFile=None):
    
    if not N:
        N=c.habitatL().len()

    # Make a copy of the input condatis object
    cl=c.habitatL()
    
    c2=cne.CondatisCoreNE(cl)
    if N == 0:
       raise ValueError("Error: empty habitat")

       
          
    c2.addSource(c.sourceL())
    c2.addTarget(c.targetL())
    c2.setParams(c.R(),c.dispersal())

    # Make an empty landscape object
    result=ls.Landscape()
    result.attribs=c.attribs

       
    N=c2.habitatL().len();
    numCells = c2.x().size

    
    minx = min(c2.x())
    maxx = max(c2.x())
    miny = min(c2.y())
    maxy = max(c2.y())
    
    estimated_time = int (2 + (numCells) * (N) / 4e8)
    maxVal = max(c2.habitatL().v)
    sumVal = sum(c2.habitatL().v)
    if (maxVal > 1.1):
       maxVal = maxVal / 100.0
       sumVal = sumVal / 100.0
    sumVal = sumVal / numCells
    
    img_tag=""
    if thumbFile is not None:
        data_uri = open(thumbFile, 'rb').read().encode('base64').replace('\n', '')
        img_tag = '<img src="data:image/png;base64,%s">' % data_uri


    previewText = "<table  cellpadding='6'><tr><td>Number of cells </td><td>  "+ str(N) + "</td> <td rowspan='8'>"\
    ""+img_tag+"</td></tr>" \
    "<tr><td> Source size </td><td>" + str(c2.sourceL().x.size)+" cells </td></tr>" \
    "<tr><td> Target size </td><td>" + str(c2.targetL().x.size)+" cells </td></tr>" \
    "<tr><td> Estimated time </td><td>" +str(estimated_time)+" minutes </td></tr>" \
    "<tr><td> East-West Range </td><td>" + str(minx) +" - "+str(maxx)+ "</td></tr>" \
    "<tr><td> North-South Range </td><td>" +str(miny) + " - " + str(maxy) +  "</td></tr>" \
    "<tr><td> Greatest cell coverage </td><td>" + str(maxVal*100) +"% </td></tr>" \
    "<tr><td> Average cell coverage </td><td>" + str(sumVal * 100) +"% </td></tr></table>"

    if previewOnly:
        if connection is not None:
           if connection.previewReport is not None:
              connection.previewReport(previewText)
        return result
      

    c2.calc()
    sp=c2.speed()


    nodeFlow = c2.nodeFlowL()
    outputDataFile = open(outputdata,'w')
    
    outputDataFile.write("Speed =, %f\n" % sp)
    outputDataFile.write("x,y,flow\n")
    
    for l in nodeFlow:
       outputDataFile.write("%f,%f,%e\n" % (l.x,l.y,l.v))
        
    if normalise is not None:
        maxVal = max(nodeFlow.v)
        nodeFlow.v = nodeFlow.v * normalise /maxVal
        
    outputDataFile.close()
         
    return nodeFlow


# external entry point    
    
def doDropping(landscapeFile = "wawhab800.tif",
    stFile = "wawst800.tif",
    restorationFile = "wawzoo800.tif",
    R = 1000.0,
    dispersal = 5.0,
    outputFile ="outfile.tif",outputdata="outdata.csv",connection = None, previewOnly = False, flowOnly=False,
    jobid=1
    ):

    log_all.log("Condatis doDropping called")                      

    
    
    thumbFile = None
    
     
    if previewOnly:
       thumbFile = landscapeFile+"_thumb.png"
       log_all.log("Making thumbnail " + thumbFile)        
       try:
         if flowOnly:
            make_thumb.makeThumb(thumbFile, landscapeFile, stFile)
            serverio.upload_thumb(thumbFile,jobid)
         else:
            make_thumb.makeThumb(thumbFile, landscapeFile, stFile, restorationFile )
            serverio.upload_thumb(thumbFile,jobid)
       except:
          print "Thumbnail error"
  

    caughtException = None
    success = False

    log_all.log("Starting process " + str(previewOnly))  
    
    folder = tempfile.mkdtemp()
    filename = os.path.join(folder,"data.hdf5")

    log_all.log("Make project "+filename)  
    h5=chdf.makeProject(filename,'root')
    log_all.log("create habitat "+landscapeFile)  
    k=createHabitat(h5,landscapeFile,stFile)
        
    log_all.log("created habitat "+landscapeFile)        
    try:

        log_all.log("Check size")  
   
        if k.sourceL().x.size == 0:
           raise ValueError("Error: empty Source")
       
        if k.targetL().x.size == 0:
           raise ValueError("Error: empty Target")     
  
        
        log_all.log("Condatis Habitat created " + str(k.habitatL().len()) +" src:  "+ str(k.sourceL().x.size) +" tgt: "+str(k.targetL().x.size) )
         
        
        k.setParams(R,dispersal)
      
        if flowOnly:

            log_all.log("Call Compute flow")            
            r=flow(k,None,connection=connection,outputdata=outputdata,previewOnly=previewOnly,jobid=jobid,thumbFile=thumbFile)
            log_all.log("Call Restore to bitmap")            
            r.restoreToBitmap()
            log_all.log("Call Export")            
            r.export(outputFile)
        else:
            restorationLandscape =ls.Landscape.fromGIS(restorationFile,makeone=False)

            maxVoltage = max(restorationLandscape.v)
            minVoltage = min(restorationLandscape.v)
            if maxVoltage > 1.1 or minVoltage < -0.1:
                restorationLandscape.scalev(1.0/100.0)

            r=drop(k,restorationLandscape,connection=connection,outputdata=outputdata,previewOnly=previewOnly,jobid=jobid,thumbFile=thumbFile )
            r.restoreToBitmap()
            r.export(outputFile)
            
            if r.lastFlow is not None:
                r.lastFlow.restoreToBitmap()
                #don't forget to amend the filenames when copying to server - see condatis.py
                r.lastFlow.export(outputFile[:-4]+"_endflow.tif")
                r.firstFlow.restoreToBitmap()
                r.firstFlow.export(outputFile[:-4]+"_startflow.tif")
            
        success = True
    except Exception, e:
        log_all.log("exception " + str(e)+" "+str(traceback.format_exc()))            
        caughtException = e
        print str(e)
    finally:
        h5.close()

        os.remove(filename)
        os.rmdir(folder)
        if caughtException is not None:
            raise caughtException
        return success

if __name__=="__main__":
    doDropping(previewOnly = False, flowOnly = False)

    
