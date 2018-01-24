import scipy as sp
from scipy import linalg
import pandas as pd
import numpy as np
import numexpr as ne
import landscape as ls
import autosourcetarget as ast
import logging


logging.root.setLevel(logging.DEBUG)

def diag(v):
    x = v.shape[0]
    y = np.identity(x)
    z = y * v
    return np.identity(v.shape[0])*v

def topNinds_(full,N=10):
    # Get the indices for the largest `num_largest` values.
    num_largest = N
    indices = (-full).argpartition(num_largest, axis=None)[:num_largest]
    xl, yl = np.unravel_index(indices, full.shape)
    return xl,yl
    
class CondatisCore(object):

    map_scale_constant = 1.0

     
        
    def __init__(self,land):
        # Habitat data
        self.x_=land.x
        self.y_=land.y
        self.ap_=land.v
        self.attribs=land.attribs
        
        # Species Parameters
        self.R_=100.0
        self.dispersal_=4.0

        # Map scaling
        if land.attribs == None:
            self.map_scale_=self.map_scale_constant
        else:
            self.map_scale_=land.attribs.geotransform[1]
            
        # Source and Target
        self.sx_=np.empty((0))
        self.sy_=np.empty((0))
        self.tx_=np.empty((0))
        self.ty_=np.empty((0))

        # Temporaries
        self.cin_=np.empty((0))
        self.cout_=np.empty((0))
        self.cfree_=np.empty((0))

        self.cumPowerThreshold=50.0

    # Access
    def R(self):
        return self.R_

    def dispersal(self):
        return self.dispersal_
    
    def x(self):
        return self.x_

    def y(self):
        return self.y_

    def ap(self):
        return self.ap_

    def sx(self):
        return self.sx_

    def sy(self):
        return self.sy_

    def tx(self):
        return self.tx_

    def ty(self):
        return self.ty_

    def hasHabitat(self):
        return len(self.x()) > 0

    def hasSource(self):
        return len(self.sx()) > 0

    def hasTarget(self):
        return len(self.tx()) > 0

    def hasFlow(self):
        f=self.nodeFlow()
        return np.sum(f) > 0

    def hasVoltage(self):
        v=self.nodeVoltage()
        return np.sum(v) > 0

    def landscape(self):
        logging.debug("condatiscore.landscape(): Don't use me. Use habitatL()")
        raise ValueError("condatiscore.landscape(): Don't use me. Use habitatL()")
        return ls.Landscape.fromVecs(self.x(),self.y(),self.ap())

    def imageSize(self):
        return self.habitatL().imageSize()
    
    def habitatL(self):
        l=ls.Landscape.fromVecs(self.x(),self.y(),self.ap())
        l.attribs=self.attribs
        return l
    
    def allHabitat(self):
        x=np.concatenate((self.x(),self.sx(),self.tx()))
        y=np.concatenate((self.y(),self.sy(),self.ty()))
        return ls.Landscape.fromVecs(x,y)

    def sourceL(self):
        l=ls.Landscape.fromVecs(self.sx(),self.sy())
        l.attribs=self.attribs
        return l

    def targetL(self):
        l=ls.Landscape.fromVecs(self.tx(),self.ty())
        l.attribs=self.attribs
        return l
               
    # Setting
    def setParams(self,R,disp):
        self.R_=R
        self.dispersal_=disp
    
    def modifyHabitat(self,land):
        self.x_=land.x
        self.y_=land.y
        self.ap_=land.v

        
    # Source and Target
    def _addSource(self,land):
        self.sx_=land.x
        self.sy_=land.y

    def addSource(self,land):
     #   land.removeDuplicates(self.habitatL())
        temp = self.habitatL()
        temp.removeDuplicates(land)        
        self.modifyHabitat(temp)
        self._addSource(land)

    def _addTarget(self,land):
        self.tx_=land.x
        self.ty_=land.y

    def addTarget(self,land):
        temp = self.habitatL()
        temp.removeDuplicates(land)        
        self.modifyHabitat(temp)
    
     #   land.removeDuplicates(self.habitatL())
        self._addTarget(land)

    # Computation stuff
    def _alpha(self):
        return 2.0/self.dispersal()

    def map_scale(self):
        return self.map_scale_
    
    def _cell(self):
        return self.map_scale()/1000.0

    def _K(self):
        return self.R()*self._alpha()**2/(2.0*np.pi)*self._cell()**4
    
    def _scind(self,v,cell):
        return (v+.5)*self._cell()

    def _calcCin(self):
        x,y,sx,sy=self.x(),self.y(),self.sx(),self.sy()
        ap=self.ap()
        K,cell=self._K(),self._cell()
        logging.debug("Cell: %f" % cell)
#        logging.debug("K: %e" % K)
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(sx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(sy,cell)[:,np.newaxis])**2)
        self.cin_=np.sum(K*ap*np.exp(-alpha*dm),axis=0)

                
    def _calcCout(self):
        x,y,tx,ty=self.x(),self.y(),self.tx(),self.ty()
        ap=self.ap()
        K,cell=self._K(),self._cell()
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(tx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(ty,cell)[:,np.newaxis])**2)
        self.cout_=np.sum(K*ap*np.exp(-alpha*dm),axis=0)
        
    def _cinAll(self):
        x,y,sx,sy=self.x(),self.y(),self.sx(),self.sy()
        ap=self.ap()
        K,cell=self._K(),self._cell()
        logging.debug("Cell: %f" % cell)
#        logging.debug("K: %e" % K)
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(sx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(sy,cell)[:,np.newaxis])**2)
        return K*ap*np.exp(-alpha*dm)
        
    def _coutAll(self):
        x,y,tx,ty=self.x(),self.y(),self.tx(),self.ty()
        ap=self.ap()
        K,cell=self._K(),self._cell()
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(tx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(ty,cell)[:,np.newaxis])**2)
        return K*ap*np.exp(-alpha*dm)
        
    def _calcFree(self):
        x,y,ap=self.x(),self.y(),self.ap()
        K,cell=self._K(),self._cell()
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(x,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(y,cell)[:,np.newaxis])**2)
        apt=ap[:,np.newaxis]
        self.cfree_=K*ap*apt*np.exp(-alpha*dm)
        dd=np.arange(self.cfree_.shape[0])
        self.cfree_[dd,dd]=0
        
    def _calcM0(self):
        self.M0_=diag(self.cin_ + self.cout_ + self.cfree_.sum(axis=0))-self.cfree_

    def _calcTLS(self):
        self.tls_=np.sum(self.cfree_) + np.sum(self.cin_) + np.sum(self.cout_)
        
    def _calcInput(self):
        self._calcCin()
        self._calcCout()
        self._calcFree()
        self._calcM0()
        self._calcTLS()
                        
    def _calcFlow(self):
        V0 = np.linalg.solve(self.M0_,self.cin_)
        I0 = np.sum(V0*self.cout_)
        Iout=V0*self.cout_
        Iin=(1-V0)*self.cin_
        self.Iin_=Iin
        self.Iout_=Iout
        self.cond_=np.sum(Iout)
        cur=self.cfree_*(V0-V0[:,np.newaxis])
        self.flo_=np.sum(np.abs(cur)/2.0,axis=0)
        self.I_=self.flo_+Iout+Iin
        self.V0_=V0
        
    def allx(self):
        return np.append(np.append(self.x(),self.sx),self.tx())

    def ally(self):
        return np.append(np.append(self.y(),self.sy),self.ty())

    def allap(self):
        inap=np.ones(self.sx().size)
        outap=np.ones(self.tx().size)
        return np.append(np.append(self.ap(),inap),outap)

    def allv(self):
        V0=self.nodeVoltage()
        Vin=np.zeros(self.sx().size)
        Vout=np.zeros(self.tx().size)
        allV=np.append(np.append(V0,Vin),Vout)
        return allV

    def _calcPower_wrong(self):
        allV=self.allv()
        self.Vij_=allV-allV[:,np.newaxis]
        #Need tthe resistances for the in and out too

        c=np.append(np.append(self.cfree_,self.cin_),self.cout_)
        print "c.shape",c.shape
        print "allV.shape",allV.shape
        self.I2ij_=self.Vij_*c/2.0
        self.Pij_=self.Vij_*self.I2ij_
        self.PijU=np.triu(self.Pij_)
        self.totalEdgePower=np.sum(self.PijU)
        self.maxEdgePower=np.max(self.PijU)
        self.edgePowerShape=self.PijU.shape
        
        pp=self.PijU
        x=self.allx()
        y=self.ally()
        ap=self.allap()
        w=np.where(pp>self.maxEdgePower/100000.0)
        pps=np.sort(pp[w].flatten())
        wx1=x[w[0]]
        wx2=x[w[1]]
        wy1=y[w[0]]
        wy2=y[w[1]]
        wp=pp[w[0],w[1]]
        
        tn=topNinds_(pp,1000)
        xl=tn[0]
        yl=tn[1]
        self.sigx1=x[xl]
        self.sigy1=y[xl]
        self.sigx2=x[yl]
        self.sigy2=y[yl]
        self.sigpow=pp[tn]
        self.pps=pps

    def _calcPower_a(self):
        self._calcHabLinks()
        self._calc_stLinks()
        asigx1=np.append(np.append(self.sighx1,self.sigsx1),self.sigtx1)
        asigy1=np.append(np.append(self.sighy1,self.sigsy1),self.sigty1)
        asigx2=np.append(np.append(self.sighx2,self.sigsx2),self.sigtx2)
        asigy2=np.append(np.append(self.sighy2,self.sigsy2),self.sigty2)
        asigpow=np.append(np.append(self.sighpow,self.sigspow),self.sigtpow)
        apps=np.append(np.append(self.hpps,self.spps),self.tpps)

        i=np.argsort(asigpow)[2000:2999]
        print "****i****",i
        self.sigx1=asigx1[i]
        self.sigy1=asigy1[i]
        self.sigx2=asigx2[i]
        self.sigy2=asigy2[i]
        self.sigpow=asigpow[i]
        self.pps=apps[i]

    def _calcPower_h(self):
        self._calcHabLinks()
        self._calc_stLinks()
        print self.sighpow.size
        print self.hpps.size
        if self.sighpow.size > 999:
            i=np.argsort(self.sighpow)[0000:999]
        else:
            i=np.argsort(self.sighpow)
        self.sigx1=self.sighx1[i]
        self.sigy1=self.sighy1[i]
        self.sigx2=self.sighx2[i]
        self.sigy2=self.sighy2[i]
        self.sigpow=self.sighpow[i]
        self.pps=self.hpps
        
    def _calcPower_s(self):
        self._calcHabLinks()
        self._calc_stLinks()
        print self.sigspow.size
        print self.spps.size
        if self.sighpow.size > 999:
            i=np.argsort(self.sighpow)[0000:999]
        else:
            i=np.argsort(self.sighpow)
        self.sigx1=self.sigsx1[i]
        self.sigy1=self.sigsy1[i]
        self.sigx2=self.sigsx2[i]
        self.sigy2=self.sigsy2[i]
        self.sigpow=self.sigspow[i]
        self.pps=self.spps
        
    def _calcPower(self):
        self._calcHabLinks()
        #self._calc_stLinks()
#        print self.sigtpow.size
#        print self.tpps.size
        if self.sighpow.size > 999:
            i=np.argsort(self.sighpow)[0000:999]
        else:
            i=np.argsort(self.sighpow)
        print i
        self.sigx1=self.sighx1[i]
        self.sigy1=self.sighy1[i]
        self.sigx2=self.sighx2[i]
        self.sigy2=self.sighy2[i]
        self.sigpow=self.sighpow[i]
        self.pps=self.hpps
        
        
    def _calcHabLinks(self):
        print "_clacHabLinks()"
        V0=self.V0()
        self.Vij_=V0-V0[:,np.newaxis]
        self.I2ij_=self.Vij_*self.cfree_/2.0
        self.Pij_=self.Vij_*self.I2ij_
        self.PijU=np.triu(self.Pij_)
        self.totalEdgePower=np.sum(self.PijU)
        self.maxEdgePower=np.max(self.PijU)
        self.edgePowerShape=self.PijU.shape

        pp=self.PijU
        x=self.x()
        y=self.y()
        ap=self.ap()
        w=np.where(pp>self.maxEdgePower/100000.0)
        pps=np.sort(pp[w].flatten())
        wx1=x[w[0]]
        wx2=x[w[1]]
        wy1=y[w[0]]
        wy2=y[w[1]]
        wp=pp[w[0],w[1]]
                
        tn=topNinds_(pp,1000)
        xl=tn[0]
        yl=tn[1]
        self.sighx1=x[xl]
        self.sighy1=y[xl]
        self.sighx2=x[yl]
        self.sighy2=y[yl]
        self.sighpow=pp[tn]
        self.hpps=pps

    def _calc_stLinks(self):
        V0=self.V0()
        # Voltages between source/hab and target/hab
        Vtg=np.zeros(self.tx().size)
        Vsr=np.ones(self.sx().size)
        VVtg=V0-Vtg[:,np.newaxis]
        VVsr=V0+Vsr[:,np.newaxis]

        print "VVtg.shape",VVtg.shape
        print "VVsr.shape",VVsr.shape

        # Conductances between source/hab and target/hab
        Ctg=self._coutAll()
        Csr=self._cinAll()
        
        print "Ctg.shape",Ctg.shape
        print "Csr.shape",Csr.shape

        # Powers between source/hab and target/hab
        Psr=VVsr*Csr/2.0
        Ptg=VVtg*Ctg/2.0

        pp=Psr
        xs=np.append(self.sx(),self.x())
        ys=np.append(self.sy(),self.y())
        ws=np.where(pp>self.maxEdgePower/100000.0)
        pps=np.sort(pp[ws].flatten())

        tn=topNinds_(pp,1000)
        yl=tn[0]
        xl=tn[1]
        self.sigsx1=xs[xl]
        self.sigsy1=ys[xl]
        self.sigsx2=xs[yl]
        self.sigsy2=ys[yl]
        self.sigspow=pp[tn]
        self.spps=pps

        pp=Ptg
        xt=np.append(self.tx(),self.x())
        yt=np.append(self.ty(),self.y())
        wt=np.where(pp>self.maxEdgePower/100000.0)
        pps=np.sort(pp[wt].flatten())

        tn=topNinds_(pp,1000)
        yl=tn[0]
        xl=tn[1]
        self.sigtx1=xt[xl]
        self.sigty1=yt[xl]
        self.sigtx2=xt[yl]
        self.sigty2=yt[yl]
        self.sigtpow=pp[tn]
        self.tpps=pps

    def _calcPower_orig(self):
        self.Vij_=V0-V0[:,np.newaxis]
        self.I2ij_=self.Vij_*self.cfree_/2.0
        self.Pij_=self.Vij_*self.I2ij_
        self.PijU=np.triu(self.Pij_)
        self.totalEdgePower=np.sum(self.PijU)
        self.maxEdgePower=np.max(self.PijU)
        self.edgePowerShape=self.PijU.shape
        
        pp=self.PijU
        x=self.x()
        y=self.y()
        ap=self.ap()
        w=np.where(pp>self.maxEdgePower/100000.0)
        pps=np.sort(pp[w].flatten())
        wx1=x[w[0]]
        wx2=x[w[1]]
        wy1=y[w[0]]
        wy2=y[w[1]]
        wp=pp[w[0],w[1]]
        
        tn=topNinds_(pp,1000)
        xl=tn[0]
        yl=tn[1]
        self.sigx1=x[xl]
        self.sigy1=y[xl]
        self.sigx2=x[yl]
        self.sigy2=y[yl]
        self.sigpow=pp[tn]
        self.pps=pps
        
    def calc(self):
        self._calcInput()
        self._calcFlow()

    def calcPower(self):
        print "AAAAAAAAAAAAAAAAAAAAAAAAAAAAARRRRGGGGGGGGGGHHHHH"
        self._calcPower()

    def sortedEdgePower(self):
        sigx1=self.sigx1
        sigy1=self.sigy1
        sigx2=self.sigx2
        sigy2=self.sigy2
        sigep=self.sigpow
        asep=np.argsort(sigep)[::-1]

        x1=sigx1[asep]
        y1=sigy1[asep]
        x2=sigx2[asep]
        y2=sigy2[asep]
        ep=sigep[asep]
        return x1,y1,x2,y2,ep
        
    def speed(self):
        return self.cond_

    def time(self):
        co=self.cond_
        if co==0:
            return 0.0
        return 1.0/self.cond_

    def totalLinkStrength(self):
        return self.tls_

    def nodeFlow(self):
        return self.I_

    def nodeFlowL(self):
        return ls.Landscape.fromVecs(self.x(),self.y(),self.nodeFlow(),attribs=self.attribs)

  
        
    def nodeVoltage(self):
        return self.V0_

    def nodeVoltageL(self):
        return ls.Landscape.fromVecs(self.x(),self.y(),self.nodeVoltage(),attribs=self.attribs)

    def mapinfo(self):
        print "Number of cells:",self.habitatL().len()
        print "EW Raster Size:",self.allHabitat().ewSize()
        print "NS Raster Size:",self.allHabitat().nsSize()

    def parameters(self):
        print "Dispersal:",self.dispersal()
        print "R:",self.R()
        
    def solution(self):
        print "Speed: %e" % self.speed()
        print "Time: %10.1e" % self.time()
        print "Resistance to extinction:",0
        print "Total link strength: %10.1f" % self.totalLinkStrength()

    def output(self):
        self.mapinfo()
        print
        self.parameters()
        print
        self.solution()
    
        
