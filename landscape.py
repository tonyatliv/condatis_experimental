import numpy as np
from osgeo import gdal
import logging
import copy
from osgeo import osr
import gisattribs as ga
import traceback
"""
Landscape class.
"""

def biggestImage(l1,l2):
    t1=l1.imageSize()
    t2=l2.imageSize()
    tx1,ty1=t1
    tx2,ty2=t2
    xr=max(tx1,tx2)
    yr=max(ty1,ty2)
    return (xr,yr)
    

class Landscape(object):
    def __init__(self):
        self.x,self.y,self.v=[np.array([],np.int32) for i in range(3)]
        self.gt=None
        self.fname=None
        self.attribs=None
        self.name="Landscape"
#        self.viewSize=None

    @classmethod
    def fromVecs(cls,x,y,v=None,attribs=None):
        """
        Create a Landscape object from vectors (1D arrays).
        The arrays should have the same length.

        Args: 
        x: (1D numpy array), x position of point.
        y: (1D numpy array), y position of point.
        v: (1D numpy array, Optional), Value at the location.(1D numpy array).
           Will be set to 1 if omitted.

        Returns: 
        Will return a Landscape object if successful.
        """     
        
   
        c=cls()
        c.x=x
        c.y=y
        if v is None:
            c.v=np.ones(x.size)
        else:
            c.v=v
        c.attribs=attribs
        return c

    @classmethod
    def fromScalers(cls,x,y,v=None,attribs=None):
        """
        Create a Landscape object from scaler arguments.

        Args: 
        x: (number), x position of point.
        y: (number), y position of point.
        v: (number, Optional), Value at the location.(1D numpy array). 
           Will be set to 1 if omitted.

        Returns: 
        Will return a Landscape object if successful.
        """         
        c=cls()
        c.x=np.array([x])
        c.y=np.array([y])
        if v==None:
            c.v=np.array([1.0])
        else:
            c.v=np.array([v])
        c.attribs=attribs
        return c

        
    @classmethod
    def fromGIS_Old(cls,filename):
        """
        Create a Landscape object from a GIS file. Uses GDAL.

        Args:
        filename: (string), Name of the file to open.

        Returns:
        Will return a Landscape object if successful.
        """
        gd=gdal.Open(filename)
        h=np.array(gd.GetRasterBand(1).ReadAsArray())
        hab=np.where(h>0)
        y=hab[0]
        x=hab[1]
        a=h[y,x]*1.0
        c=cls.fromVecs(x,y,a)

    @classmethod
    def fromGIS(cls,filename, makeone = False):
        """
        Create a Landscape object from a GIS file. Uses GDAL.

        Args:
        filename: (string), Name of the file to open.

        Returns:
        Will return a Landscape object if successful.
        """
        gd=gdal.Open(filename)
        h=np.array(gd.GetRasterBand(1).ReadAsArray())
        
        np.set_printoptions(threshold=np.nan)
        
        hab=np.where(h>0)
        dataset = gd
        geotransform = dataset.GetGeoTransform()
        
        y=hab[0]
        x=hab[1]
        a=h[y,x]*1.0 
        
        x = x * geotransform[1]  
        x = x + geotransform[0] 
        
        y = y * geotransform[5] 
        y = y + geotransform[3]
        
        x = x / 1000
        y = y / 1000
        
        y = y - 0.5
        x = x + 0.5
 
        if makeone:
            a = np.ones(a.size)
 
        
        attr=ga.GISAttribs.fromFile(filename)
      
        c=cls.fromVecs(x,y,a,attribs=attr)
        return c
        
    def restoreToBitmap(self):    
        geotransform = self.attribs.geotransform
        self.y = self.y + 0.5
        self.x = self.x - 0.5
        
        self.x = self.x * 1000
        self.y = self.y * 1000

        self.x = self.x - geotransform[0] 
        self.y = self.y - geotransform[3]
        
        self.x = self.x / geotransform[1]  
        self.y = self.y / geotransform[5] 
        
     
        
        
        
    @classmethod
    def fromImage(cls,h, mincut=0.0):
        """
        Create a Landscape object from a 2D numpy array.

        Args:
        h: (2D numpy array), 2D array to extract locations from.  
        mincut: (number, optional), Threshold over which data will be 
                extracted from the array. Defaults to zero if omitted.

        Returns:
        Will return a Landscape object if successful.
        """
        hab=np.where(h>0.0)
        y=hab[0]
        x=hab[1]
        a=h[y,x]
        return cls.fromVecs(x,y,a)

    @classmethod
    def fromNode(cls,node):
        """
        *** NOT IMPLEMENTED ***
        Create a Landscape object from a node in an HDF5 file. Uses pytables.

        Args:
        node: (node in an HDF5 file), The node where the data will be read from.
        
        Returns:
        Will return a Landscape object if successful.
        """
        pass

    def isin_(self,xnew,ynew,x,y):
        for i in range(x.size):
            if (xnew == x[i] and ynew == y[i]):
                return True
        return False

    def removeDups_(self,x,y,ap,xnew,ynew,apnew):
        xr=[]
        yr=[]
        apr=[]
        for i in range(xnew.size):
            if not self.isin_(xnew[i],ynew[i],x,y):
                xr.append(xnew[i])
                yr.append(ynew[i])
                apr.append(apnew[i])
        nxr=np.array(xr)
        nyr=np.array(yr)
        napr=np.array(apr)
        return nxr,nyr,napr

    def __repr__(self):
        s="Landscape Object:\n"
        s+="Size: %d" % self.len() + '\n'
        s+="X: "+np.array_str(self.x) + '\n'
        s+="Y: "+np.array_str(self.y) + '\n'
        s+="V: "+np.array_str(self.v)
        return s

    def __gt__(self,v):
        return self.v > v

    def __lt__(self,v):
        return self.v < v

    def __ge__(self,v):
        return self.v >= v

    def __le__(self,v):
        return self.v <= v

    
    def __add__(self,l):
        v=self.v+l.v
        return self.fromVecs(self.x,self.y,v,attribs=self.attribs)

    def __mul__(self,l):
        v=self.v*l.v
        return self.fromVecs(self.x,self.y,v,attribs=self.attribs)
        
    def __add__old(self,l):
        logging.debug("Overloaded add operator. Don't use me! Use 'append()' instead.")
        x=np.append(self.x,l.x)
        y=np.append(self.y,l.y)
        v=np.append(self.v,l.v)
        return self.fromVecs(x,y,v)

    def reloadAttribs(self):
        fn=self.attribs.filename
        self.attribs=ga.GISAttribs(fn)
    
    def __getitem__(self,index):
        x=np.array(self.x[index])
        y=self.y[index]
        v=self.v[index]
        return self.fromVecs(x,y,v)

    def argmax(self):
        """
        Returns the index of the maximum value in the landscape (argmax).

        Args:

        Returns:
        The index of the maximum value in the landscape.
        """
        return np.argmax(self.v)

    def max(self):
        """
        Returns the point in the landscape with the biggest value.

        Args:

        Returns:
        A landscape object with a length of 1 containing the point with the 
        biggest value.
        """
        return self.__getitem__(self.argmax())

    def argmin(self):
        """
        Returns the index of the minimum value in the landscape (argmin).

        Args:

        Returns:
        The index of the minimum value in the landscape.
        """
        return np.argmin(self.v)
    
        
    maxWidth = None
    
    @classmethod
    def make_location_key(cls,x,y):
     # maxWidth is just used to make a numerical key from x,y co-ordinates, (i.e. key = y*maxWidth + x) 
     # and hence slightly faster lookup to check for set membership

        return (y*cls.maxWidth)+x
        
    @classmethod	
    def getKeySet(cls,restorationLandscape):
        if cls.maxWidth is None:
            cls.maxWidth = int(max(restorationLandscape.x)*2)
        habitatList = []
        for i in range(restorationLandscape.x.size):
            key = Landscape.make_location_key(restorationLandscape.x[i],restorationLandscape.y[i])
            habitatList.append(key)
   
        return set(habitatList)  
        
 
 
    def argminBoth(self, restorationLandscapeSet, ignoreSet = set()):
        """
        Returns the index of the minimum value in both the landscape and the given restoration landscape set(argmin).

        Args:
        restorationLandscapeSet: (set) A set of keys that describe just the landscape to be restored
        Returns:
        The index of the minimum value in the landscape.
        """

        
   

        lowest = -1
        for i in range(1, self.v.size):
            if lowest < 0 or self.v[i] < self.v[lowest]:
                key = self.make_location_key(self.x[i],self.y[i])
                if key in restorationLandscapeSet:
                    if key not in ignoreSet:
                        lowest = i
     
        if lowest < 0:
             raise ValueError("cannot find a cell to drop which is in the restorationLandscapeSet")     
        return lowest;
        
        
    def min(self):
        """
        Returns the point in the landscape with the smallest value.

        Args:

        Returns:
        A landscape object with a length of 1 containing the point with the 
        smallest value.
        """
        return self.__getitem__(self.argmin())

    def level(self,level):
        """
        Extract points from a landscape that have a single value.

        Args:
        level: (number) The value that you want to extract.

        Returns:
        Returns a Landscape object with the points that have been extracted.
        """
        w=np.where(self.v==level)
        return self.fromVecs(self.x[w],self.y[w],self.v[w])

    def where(self,q):
        return np.where(q)

    def assign(self,l):
        self.x=l.x
        self.y=l.y
        self.v=l.v
        self.attribs=l.attribs
        
    def len(self):
        """
        Returns the number of points in the landscape.

        Args:

        Returns:
        The integer number of points in the landscape.
        """
        return self.x.size

    def xlen(self):
        """
        Returns the number of points in the x array.

        Args:

        Returns:
        The integer number of points in the x array.
        """
        return self.x.size

    def ylen(self):
        """
        Returns the number of points in the y array.

        Args:

        Returns:
        The integer number of points in the y array.
        """
        return self.y.size

    def imageSize(self):
        """
        Returns the size of the 2D array that would be produced from this landscape.

        Args:

        Returns:
        (x,y) The size of the image.
        """
        if self.attribs==None:
            return (np.max(self.x)+1,np.max(self.y)+1)
        else:
            return (self.attribs.xsize,self.attribs.ysize)

        # if self.viewSize==None:
        #     return (np.max(self.x)+1,np.max(self.y)+1)
        # else:
        #     return self.viewSize
        
    def _biggest(self,t1,t2):
        tx1,ty1=t1
        tx2,ty2=t2
        xr=max(tx1,tx2)
        yr=max(ty1,ty2)
        return (xr,yr)
    
    def image(self,size=None):
        """
        Returns an image (2D array) representation of the landscape.

        Args:
        size: (int,int), The size of the image in pixels.

        Returns:
        A 2D numpy array A, where A(x,y)=V
        """
        
 
        if size==None:
            a=np.zeros(self.imageSize())
        else:
            z=self.imageSize()
            z2=self._biggest(size,z)
            a=np.zeros(z2)
        a[self.x.astype(int),self.y.astype(int)]=self.v
        return a

    ####### Obsolete ##########
    def arraySize(self):
        logging.debug("landscape.arraySize(). Don't use me! Use imageSize()")
        return (np.max(self.x)+1,np.max(self.y)+1)

    def array(self,size=None):
        logging.debug("landscape.array(). Don't use me! Use image()")
        if size==None:
            a=np.zeros(self.arraySize())
        else:
            a=np.zeros(size)
        a[self.x.astype(int),self.y.astype(int)]=self.v
        return a
    ###########################

    def dataSize(self):
        """
        Returns the (x,y) size (position of top right corner) in the landscape.

        Args:

        Returns:
        (number,number), The (x,y) size (position of top right corner) in the landscape.
        """
        return (np.max(self.x),np.max(self.y))

    def ewSize(self):
        """
        The maximum x size in the landscape.

        Args:

        returns:
        (Number) The maximum x size in the landscape.
        """
        return self.dataSize()[0]
    
    def nsSize(self):
        """
        The maximum y size in the landscape.

        Args:

        returns:
        (Number), The maximum y size in the landscape.
        """
        return self.dataSize()[1]
        
    def setConstant(self,k):
        """
        Set the v array to a constant.

        Args:
        k: (number), The value that will be set into each element of v
        """
        self.v[:]=k

    def removeDuplicates(self,testl):
        """
        Remove any duplicates found in testl from the landscape.

        Args:
        testl (Landscape), The landscape to test against.

        Returns:
        """
        x,y,v=self.removeDups_(testl.x,testl.y,testl.v,self.x,self.y,self.v)
        self.x=x
        self.y=y
        self.v=v

    def copy(self):
        """
        Return a copy of the landscape.

        Args:
        
        Returns:
        (Landscape) An exact copy of the landscape.
        """
        return copy.deepcopy(self)
        # x=copy.deepcopy(self.x)
        # y=copy.deepcopy(self.y)
        # v=copy.deepcopy(self.v)
        # return self.fromVecs(x,y,v)

    def delete(self,inds):
        """
        Return a new landscape with inds deleted.
        """
#        logging.debug("WARNING: map projection info lost")
        c=self.copy()
        c.x=np.delete(self.x,inds)
        c.y=np.delete(self.y,inds)
        c.v=np.delete(self.v,inds)
        return c
    
    def append(self,l,clip=True):
        """
        Append another landscape object and return the result.
        Does not modify either of the landscapes being appended.
        See also 'grow()'
        """
#        xz,yz=self.imageSize()
        c=self.copy()
        c.x=np.append(self.x,l.x)
        c.y=np.append(self.y,l.y)
        c.v=np.append(self.v,l.v)
        return c

    def grow(self,l,clipsize=None):
        """
        Append a new landscape. This will modify the landscape object.
        Returns nothing. See also 'append()'.
        """
        if not clipsize == None:
            l.clip(clipsize)
        self.x=np.append(self.x,l.x)
        self.y=np.append(self.y,l.y)
        self.v=np.append(self.v,l.v)

    def clip(self,clipsize):
        """
        Remove elements outside bounds.
        Elements with x or y values < 0 also removed.
        """
        xz,yz=clipsize
        w=[(self.x<0) | (self.x>xz) | (self.y<0) | (self.y>yz)]
        w=np.logical_not(w[0])
        self.x=self.x[w]
        self.y=self.y[w]
        self.v=self.v[w]
        
    def scalexy(self,k):
        self.x=self.x*k
        self.y=self.y*k
    
    def scalev(self,k):
        self.v=self.v*k

    def cellArea(self):
        return

    def _cellArea(self):
        return self.attribs.xscale**2
    
    def normaliseHabitat(self,id):
        print "+++++++++++  normaliseHabitst()+++++++++++++++",id
        a=self.v
        ca=self._cellArea()
        # Note case -3 requires no alteration (proportion)
        if id==-2:
            print "Percentage. Divide by 100"
            print "Max area",np.max(a)
            a/=100.0
        if id==-3:
            a*=1.0
        if id==-4:
            print "Area. Divide by cell area (m^2)"
            a/=ca
        if id==-5:
            print "Area. Divide by cell area (km^2)"
            print "Cell area:",ca
            a/=ca/1e6
        if id==-6:
            print "Load as mask"
            a=a*0+1.0
        self.v=a

    def export_old(self,fname):
        writeGeoFile(self.image(),fname,self.attribs.fileName)

    def export(self,newFileName):
        attribs=self.attribs
        a=self.image()
        print a.shape
        NDV=0
        if attribs:
            NDV = attribs.ndv
        if not NDV:
            NDV=0
        xsize = attribs.xsize
        ysize = attribs.ysize
        GeoT = attribs.geotransform
        Projection = osr.SpatialReference()
        Projection.ImportFromWkt(attribs.projectionref)
        dtype=gdal.GDT_Float64
        driver = gdal.GetDriverByName("GTiff")
        DataSet = driver.Create( newFileName, xsize, ysize, 1, dtype )
        DataSet.SetGeoTransform(GeoT)
        DataSet.SetProjection( Projection.ExportToWkt() )
        DataSet.GetRasterBand(1).WriteArray(np.transpose(a).astype(np.float64))
        DataSet.GetRasterBand(1).SetNoDataValue(NDV)
        logging.info("GIS file %s written" % newFileName)

    def exportCSV(self,fname):
        with open(fname, 'w') as f:
            for i in range(self.len()):
                s="%i,%e,%e,%e\n" % (i,self.x[i],self.y[i],self.v[i])
                f.write(s)

    def toLonLat(self):
        """
        Not properly implemented yet.
        Needs finishing off.
        Lon lat not correct and has hardcoded origin.
        Need to figure out how to get pixel coordinates.
        """
        g=self.attribs
        xoff, a, b, yoff, d, e = g.geotransform
        xp = a*self.x + b*self.y + xoff
        yp = d*self.x + e*self.y + yoff

        xp=xp*1e-5
        yp=yp*1e-5
        xp=xp-2.0
        yp=yp+49.0
        l=self.fromVecs(xp,yp,self.v)
        l.attribs=self.attribs
        return l

