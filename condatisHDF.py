import condatiscore as cc
import condatisNE as cne
import landscape as ls
import autosourcetarget as ast
import logging
import numpy as np
import numexpr as ne
import tables
import gisattribs as ga

def makeProject(projname,title):
    h5=tables.open_file(projname,mode="w",title=title)
    h5.create_group("/",'scenarios','Scenarios')
    h5.root.scenarios._v_attrs.nextIndex=0
    return h5


class CondatisCoreHDF(cc.CondatisCore):
    # Constructors
    def __init__(self,h5):
        self.h5=h5
        self.scenario=None
        self.cumPowerThreshold=50.0
        self.attribs=None
        self.map_x_scale=1000.0

    @classmethod
    def fromLandscape(cls,h5,scname,land):
        c=cls(h5)
        c.makeScenario(scname)
        c.x_,c.y_,c.ap_,c.sx_,c.sy_,c.tx_,c.ty_=[None for i in range(7)]
        c._addHab(land)
        c.attribs=land.attribs
        
        c.generateCombinedHabitat()
        return c

    
 
 
        
    def map_scale(self):
        return self.map_x_scale

 

        # Access
    def getScenario(self,name):
        return self.h5.getNode("/scenarios",name)
    
    def R(self):
        return self.scenario._v_attrs.R

    def dispersal(self):
        return self.scenario._v_attrs.dispersal

    def _cell(self):
        scn=self.scenario
        return (scn._v_attrs.map_x_scale*1.0)/1000.0
    
    def x(self):
        if not self.scenario.__contains__("x"):
            return np.empty(shape=0) 
        return self.x_.read()

    def y(self):
        if not self.scenario.__contains__("y"):
            return np.empty(shape=0) 
        return self.y_.read()

    def ap(self):
        if not self.scenario.__contains__("ap"):
            return np.empty(shape=0) 
        return self.ap_.read()

    def sx(self):
        if not self.scenario.__contains__("or_x"):
            return np.empty(shape=0) 
        return self.scenario.or_x.read()

    def sy(self):
        if not self.scenario.__contains__("or_y"):
            return np.empty(shape=0) 
        return self.scenario.or_y.read()

    def tx(self):
        if not self.scenario.__contains__("tg_x"):
            return np.empty(shape=0) 
        return self.scenario.tg_x.read()

    def ty(self):
        if not self.scenario.__contains__("tg_y"):
            return np.empty(shape=0) 
        return self.scenario.tg_y.read()

    def nodeVoltage(self):
        if not self.scenario.__contains__("V0"):
            logging.debug("No Voltage data")
            return self.x()*0.0
        return self.scenario.V0.read()

    def nodeFlow(self):
        if not self.scenario.__contains__("I"):
            logging.debug("No flow data")
            return self.x()*0.0
        return self.scenario.I.read()

    def speed(self):
        if not self.scenario._v_attrs.__contains__("I0"):
            return 0
        return self.scenario._v_attrs.I0

    def time(self):
        co=self.speed()
        if co==0:
            return 0.0
        return 1.0/self.speed()

    def totalLinkStrength(self):
        if not self.scenario._v_attrs.__contains__("totalLinkStrength"):
            return 0
        return self.scenario._v_attrs.totalLinkStrength

    # Making and deleting scenarios
    def hasScenario(self,name):
        return self.h5.root.scenario.__contains__(name)

    def deleteScenario(self,name):
        node="/scenarios/"+name
        self.h5file.removeNode(node,recursive=True)
        
    def makeScenario(self,name):
        h5=self.h5
        self.scenario=self.h5.create_group("/scenarios",name,"Scenario")
        gname="/scenarios/"+name
      #  h5.createGroup(gname,"input","Landscape data")
        h5.create_group(gname,"combined","Landscape combined with source and sink")
      #  h5.createGroup(gname,"metrics","Population Metrics")
      #  h5.createGroup(gname,"plugins","Plugins")
        self.scenarioDefaults()
        self.scenario._v_attrs.ind=h5.root.scenarios._v_attrs.nextIndex
        h5.root.scenarios._v_attrs.nextIndex+=1
        return self.scenario

    def scenarioDefaults(self,R=100.0,dispersal=4.0):
        scenario=self.scenario
        scenario._v_attrs.R=R
        scenario._v_attrs.dispersal=dispersal
        scenario._v_attrs.mapscale=1.0
        scenario._v_attrs.areascale=0
        scenario._v_attrs.I0=0
            
    def setParams(self,R,dispersal):
        scenario=self.scenario
        scenario._v_attrs.R=R
        scenario._v_attrs.dispersal=dispersal

    def generateCombinedHabitat(self):
        scenario=self.scenario
        h5=self.h5
        x=scenario.x.read()
        y=scenario.y.read()
        ap=scenario.ap.read()
#        cell=scenario.cell.read()
        if scenario.__contains__('or_x') and scenario.__contains__('tg_x'):
            or_x=scenario.or_x.read()
            or_y=scenario.or_y.read()
            tg_x=scenario.tg_x.read()
            tg_y=scenario.tg_y.read()
            xa=np.concatenate((x,or_x,tg_x))
            ya=np.concatenate((y,or_y,tg_y))
        else:
            xa,ya=x,y
#        don't do anything?
#        shape=xa.shape
#        atom=tables.Int32Atom(shape=())
        filters = tables.Filters(complevel=0, complib='zlib')    
        if scenario.combined.__contains__('x'):
            scenario.combined.x.remove()
        if scenario.combined.__contains__('y'):
            scenario.combined.y.remove()
        if scenario.combined.__contains__('V'):
            scenario.combined.V.remove()
        h5.create_array(scenario.combined,'x',xa)
        h5.create_array(scenario.combined,'y',ya)
        h5.create_array(scenario.combined,'V',xa*0)
        
    def modifyHabitat(self,land):
        self._addHab(land)
        
    # Adding data.
    def _addHab(self,land):
        h5=self.h5
        scenario=self.scenario
        if scenario.__contains__('x'):
            scenario.x.remove()
        
        if scenario.__contains__('y'):
            scenario.y.remove()

        if scenario.__contains__('ap'):
            scenario.ap.remove()

        if scenario.__contains__('cell'):
            scenario.cell.remove()

        self.x_=h5.create_array(scenario,'x',land.x)
        self.y_=h5.create_array(scenario,'y',land.y)
        self.ap_=h5.create_array(scenario,'ap',land.v)
        # Can I remove this?
        #h5.createArray(scenario,'cell',np.ones(land.x.size))

        if not land.attribs==None:
            logging.info("Writing attributes to HDF file")
            land.attribs.write(scenario)
            # gt=land.attribs.geotransform
            # scx=np.abs(gt[1]) # Because we get a negative number in the scale
            # scy=np.abs(gt[5]) 
            # orx=gt[0]
            # ory=gt[3]
            # projection=land.attribs.projection
            # projectionName="Not Implemented"
            
            # scenario._v_attrs.filename=land.attribs.filename
            # scenario._v_attrs.datatype=land.attribs.datatype

            # scenario._v_attrs.map_x_size=land.ewSize()+1
            # scenario._v_attrs.map_y_size=land.nsSize()+1
            # scenario._v_attrs.map_x_scale=scx
            # scenario._v_attrs.map_y_scale=scy
            # scenario._v_attrs.map_x_origin=orx
            # scenario._v_attrs.map_y_origin=ory
            # scenario._v_attrs.map_projection=projection
            # scenario._v_attrs.map_projectionName=projectionName
            
    def _addSource(self,land):
        print "add Source hdf5"
        h5=self.h5
        scenario=self.scenario
        if scenario.__contains__('or_x'):
            scenario.or_x.remove()
        if scenario.__contains__('or_y'):
            scenario.or_y.remove()
        self.sx_=h5.create_array(scenario,'or_x',land.x)
        self.sy_=h5.create_array(scenario,'or_y',land.y)

    def _addTarget(self,land):
        h5=self.h5
        scenario=self.scenario
        if scenario.__contains__('tg_x'):
            scenario.tg_x.remove()
        if scenario.__contains__('tg_y'):
            scenario.tg_y.remove()
        self.tx_=h5.create_array(scenario,'tg_x',land.x)
        self.ty_=h5.create_array(scenario,'tg_y',land.y)

    # Calculating
    def _deleteOldCalc(self):
        sc=self.scenario
        if sc.__contains__('V0'):
            sc.V0.remove()
        if sc.__contains__('Vij'):
            sc.Vij.remove()
        if sc.__contains__('M0'):
            sc.M0.remove()
        if sc.__contains__('I'):
            sc.I.remove()
        if sc.__contains__('I2ij'):
            sc.I2ij.remove()
        if sc.__contains__('ipv_in'):
            sc.ipv_in.remove()
        if sc.__contains__('ipv_out'):
            sc.ipv_out.remove()
        if sc.__contains__('ipv_free'):
            sc.ipv_free.remove()

    def _saveCalc(self):
        sc=self.scenario
        h5=self.h5
        self.V0_=h5.create_array(sc,'V0',self.V0_)
        self.I_=h5.create_array(sc,'I',self.I_)
        sc._v_attrs.I0=self.cond_
        sc._v_attrs.totalLinkStrength=self.tls_
        sc.cin_=h5.create_array(sc,"ipv_in",self.cin_)
        sc.cout_=h5.create_array(sc,"ipv_out",self.cout_)
        sc.free_=h5.create_array(sc,"ipv_free",self.cfree_)
        h5.flush()

    def calc(self):
        r=super(CondatisCoreHDF,self).calc()
        self._deleteOldCalc()
        self._saveCalc()

    def _deleteOldPowCalc(self):
        sc=self.scenario
        h5=self.h5
        if sc.__contains__('sig_pow'):
            sc.sig_pow.remove()
        if sc.__contains__('sorted_sig_power'):
            sc.sorted_sig_power.remove()
        if sc.__contains__('sigx1'):
            sc.sigx1.remove()
        if sc.__contains__('sigx2'):
            sc.sigx2.remove()
        if sc.__contains__('sigy1'):
            sc.sigy1.remove()
        if sc.__contains__('sigy2'):
            sc.sigy2.remove()
        if sc.__contains__('edgePower'):
            sc.edgePower.remove()

        h5.flush()

    def _savePowCalc(self):
        sc=self.scenario
        h5=self.h5
        h5.create_array(sc,'edgePower',self.Pij_)
        sc._v_attrs.totalEdgePower=self.totalEdgePower
        sc._v_attrs.maxEdgePower=self.maxEdgePower
        sc._v_attrs.edgePowerShape=self.PijU.shape
        h5.create_array(sc,'sig_pow',self.sigpow)
        h5.create_array(sc,'sorted_sig_power',self.pps)
        h5.create_array(sc,'sigx1',self.sigx1)
        h5.create_array(sc,'sigx2',self.sigx2)
        h5.create_array(sc,'sigy1',self.sigy1)
        h5.create_array(sc,'sigy2',self.sigy2)

    def hasPower(self):
        sc=self.scenario
        return sc._v_attrs.__contains__('maxEdgePower')
    
    def hasCalc(self):
        sc=self.scenario
        return sc._v_attrs.__contains__("I0")

    def calcPower(self):
        super(CondatisCoreHDF,self).calcPower()
        self._deleteOldPowCalc()
        self._savePowCalc()

    # Results access
    def cond(self):
        return self.scenario._v_attrs.I0

    def V0(self):
        return self.V0_.read()

    def Vij(self):
        return self.Vij_.read()
