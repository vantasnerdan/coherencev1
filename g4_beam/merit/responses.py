"""Compiled frozen-profile responses using ONLY consumer.csda_nodes geometry.

This is a fit accelerator, not a second stopping/yield model. Final predictions
and reported chi2 are evaluated through forward_yields at fitted parameters.
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline
from g4_beam.forward.consumer import csda_nodes, QE
from g4_beam.forward import model as M
from g4_beam.forward.profile_thermal import transport_stopping


class Response:
    def __init__(self,e,profile):
        self.e=e;layers,self.nads,_=profile
        self.layers=tuple(tuple(float(v) for v in r) for r in layers)
        gx,gw=np.polynomial.hermite_e.hermegauss(5);gw/=gw.sum()
        self.tracks=[(float(2*e*(1+.05/2.3548*g)),float(w*.85)) for g,w in zip(gx,gw)]
        self.tracks += [(float(e*(1+.05/2.3548*g)),float(w*.30)) for g,w in zip(gx,gw)]
        en=[];zz=[];ww=[];incident=[]
        for energy,weight in self.tracks:
            a,z,w,_=csda_nodes(energy,self.layers,transport_stopping)
            a=np.r_[energy,a];z=np.r_[0.,z];w=np.r_[self.nads*1e-31,w]*weight
            en.append(a/2);zz.append(z);ww.append(w);incident.append(np.full(len(a),energy/2))
        self.ec=np.concatenate(en);self.z=np.concatenate(zz);self.w=np.concatenate(ww);self.incident=np.concatenate(incident)
        self.neutral_tables=None

    def ions(self,ub,us,cause,p,ch):
        if cause=='surface':
            value=0.;boundary=p['d']*1e-9
            for energy,weight in self.tracks:
                a,z,w,_=csda_nodes(energy,self.layers,transport_stopping,break_m=boundary)
                u=np.where(z<boundary,us,ub)
                value+=weight*(float(np.sum(w*M.reaction_sigma(np.maximum(a/2,.05),u,ch)))+
                               self.nads*1e-31*float(M.reaction_sigma(max(energy/2,.05),us,ch)))
            return value
        u=ub
        if cause=='screeningrise':u=ub+(us-ub)/(1+np.exp((self.incident-p['E0'])/.5))
        return float(np.sum(self.w*M.reaction_sigma(np.maximum(self.ec,.05),u,ch)))

    def build_neutral_table(self):
        energy=np.geomspace(.5,26.,64)
        u=np.r_[0.,np.geomspace(.005,15.,80)]
        values={ch:[] for ch in ('p','n')}
        for e in energy:
            en,z,w,_=csda_nodes(e,self.layers,transport_stopping)
            ec=np.maximum(np.r_[e,en]/2,.05);weight=np.r_[self.nads*1e-31,w]
            for ch in values:
                y=np.sum(M.reaction_sigma(ec[:,None],u[None,:],ch)*weight[:,None],axis=0)
                values[ch].append(np.log(y))
        self.neutral_tables={ch:RectBivariateSpline(np.log(energy),np.sqrt(u+.05),np.array(vals),kx=3,ky=3) for ch,vals in values.items()}

    def neutral(self,e,u,ch):
        if self.neutral_tables is None:self.build_neutral_table()
        return float(np.exp(self.neutral_tables[ch](np.log(e),np.sqrt(u+.05))[0,0]))
