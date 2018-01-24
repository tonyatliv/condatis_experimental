import condatiscore as cc
import landscape as ls
import autosourcetarget as ast
import logging
import numpy as np
import numexpr as ne
import log_all

class CondatisCoreNE(cc.CondatisCore):
    def __init__(self,land):
        cc.CondatisCore.__init__(self,land)

    def _calcCin(self):
        x,y,sx,sy=self.x(),self.y(),self.sx(),self.sy()
        ap=self.ap()
        K,cell=self._K(),self._cell()
#        logging.debug("Cell: %f" % cell)
#        logging.debug("K: %e" % K)
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(sx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(sy,cell)[:,np.newaxis])**2)
        self.cin_=ne.evaluate("sum(K*ap*exp(-alpha*dm),axis=0)")

    def _calcCout(self):
        x,y,tx,ty=self.x(),self.y(),self.tx(),self.ty()
        ap=self.ap()
        K,cell=self._K(),self._cell()
 #       logging.debug("Cell: %f" % cell)
 #       logging.debug("K: %e" % K)
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(tx,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(ty,cell)[:,np.newaxis])**2)
        self.cout_=ne.evaluate("sum(K*ap*exp(-alpha*dm),axis=0)")

    def _calcFree(self):
        x,y,ap=self.x(),self.y(),self.ap()
        K,cell=self._K(),self._cell()
        alpha=self._alpha()
        dm=np.sqrt((self._scind(x,cell)-self._scind(x,cell)[:,np.newaxis])**2  + (self._scind(y,cell)-self._scind(y,cell)[:,np.newaxis])**2)
        apt=ap[:,np.newaxis]
        self.cfree_=ne.evaluate("K*ap*apt*exp(-alpha*dm)")
        dd=np.arange(self.cfree_.shape[0])
        self.cfree_[dd,dd]=0

    # Same as CondatisCore. Reminder to re-do with NE
    def _calcM0(self):
        self.M0_=cc.diag(self.cin_ + self.cout_ + self.cfree_.sum(axis=0))-self.cfree_

    # Same as CondatisCore. Reminder to re-do with NE
    def _calcTLS(self):
        self.tls_=np.sum(self.cfree_) + np.sum(self.cin_) + np.sum(self.cout_)

