import sys
import landscape as ls
import numpy as np
import scipy.misc
from PIL import ImageOps
from PIL import Image     
from PIL import ImageEnhance 
import PIL.ImageFile
import ftplib
from paramiko import SSHClient
import paramiko
import log_all

def makeImage(file, filterV =  None):
    land =ls.Landscape.fromGIS(file,makeone=False)
    land.restoreToBitmap()
 
    if filterV is not None:
        land.v = land.v == filterV
    image = land.image()
    image = np.transpose(image)

    return scipy.misc.toimage(image)


    
def makeThumb(outfile, habfile, stfile, zoofile = None):
    him = makeImage(habfile)
    sim = makeImage(stfile,1)
    tim = makeImage(stfile,2)
     
   
    size = 500, 400

 
    him = ImageOps.colorize(him, (0, 0, 0, 0), "#CCCCCC") 
    him.thumbnail(size, Image.ANTIALIAS)

  
    sim = ImageOps.colorize(sim, (0, 0, 0, 0), "#4040FF") 
    sim.thumbnail(size, Image.ANTIALIAS)

 
    tim = ImageOps.colorize(tim, (0, 0, 0, 0), "#FF3030") 
    tim.thumbnail(size, Image.ANTIALIAS)
    stim = Image.blend(sim, tim, 0.5) 
    
    enhancer = ImageEnhance.Brightness(stim)
    stim = enhancer.enhance(2)
     

     
    if zoofile is not None:
        zim = makeImage(zoofile)
        zim = ImageOps.colorize(zim, (0, 0, 0, 0), "#888822") 
        zim.thumbnail(size, Image.ANTIALIAS)
  
    if zoofile is None:  
        new_img = Image.blend(him, stim, 0.5)   
    else:
        new_img = Image.blend(zim, him, 0.5)      
        enhancer = ImageEnhance.Brightness(new_img)
        new_img = enhancer.enhance(2)  
        new_img = Image.blend(new_img, stim, 0.5)        
    
    enhancer = ImageEnhance.Brightness(new_img)
    
    new_img = enhancer.enhance(2)
    
    
    new_img.convert('RGBA')
    
    new_img.save(outfile, "PNG")
    
    log_all.log("Writing thumbnail " + outfile)     
   

  
    return

        
  