"""Exact Cox convolution of a finite positive delay measure.

The macro operator consumes accumulated interface fluxes. The underlying ray
ages, boundary locations and fast-window packet Green are unchanged.
"""
import numpy as np


class DelayMeasure:
    def __init__(self,age,weight):
        age=np.asarray(age,float);weight=np.asarray(weight,float)
        order=np.argsort(age);ordered=age[order]
        starts=np.r_[0,np.flatnonzero(np.diff(ordered)!=0)+1] if len(age) else np.zeros(0,dtype=int)
        self.age=ordered[starts];self.weight=np.add.reduceat(weight[order],starts) if len(age) else weight.copy()
        self.mass=np.r_[0.,np.cumsum(self.weight)];self.moment=np.r_[0.,np.cumsum(self.weight*self.age)]
        self.total=float(self.mass[-1]);self.first=float(self.moment[-1]);self.maximum=float(self.age[-1]) if len(age) else 0.
        self._frequency=None
    def _index(self,x):return np.searchsorted(self.age,np.maximum(x,0.),side='right')
    def accumulated(self,x):
        x=np.asarray(x);j=self._index(x)
        return np.where(x>0,x*self.mass[j]-self.moment[j],0.)
    def residence_primitive(self,x):
        x=np.asarray(x);j=self._index(x)
        return np.where(x>0,self.moment[j]+x*(self.total-self.mass[j]),0.)
    def _burst(self,calendar,time,primitive):
        cutoff=max(0.,time-self.maximum);edge=calendar.edges[:-1];state=calendar.state
        lo=int(np.searchsorted(edge,cutoff,side='right')) if time>=self.maximum else 0
        hi=int(np.searchsorted(edge,time,side='right'))
        if getattr(self,'_calendar',None) is not calendar:
            self._calendar=calendar;self._jumps=np.diff(state.astype(float),prepend=0.)
        jumps=self._jumps
        if lo:
            old_state=float(state[lo-1]);old_integral=float(calendar.burst_integral[lo-1]+(cutoff-edge[lo-1])*state[lo-1])
            old=(self.total*(old_integral+old_state*(time-cutoff))-self.first*old_state) if primitive else self.first*old_state
        else:old=0.
        evaluate=self.accumulated if primitive else self.residence_primitive
        return old+float(jumps[lo:hi]@evaluate(time-edge[lo:hi]))
    def _oscillatory(self,calendar,time,primitive):
        w=calendar.omega;phase=calendar.phase;j=int(self._index(time))
        if self._frequency!=w:
            self._prefix_phase=np.r_[0j,np.cumsum(self.weight*np.exp(-1j*w*self.age))]
            self._prefix_difference=np.r_[0j,np.cumsum(self.weight*(-np.expm1(-1j*w*self.age)))]
            self._frequency=w
        if primitive:
            value=np.exp(1j*(w*time+phase))*self._prefix_phase[j]-np.exp(1j*phase)*self.mass[j]
        else:
            value=np.exp(1j*(w*time+phase))*self._prefix_difference[j]+np.exp(1j*phase)*np.expm1(1j*w*time)*(self.total-self.mass[j])
        return .7/w*float(value.imag)
    def departed(self,calendar,time,name):
        if time<=0:return 0.
        result=float(self.accumulated(time))
        if name=='hot_bursts':result=(result+self._oscillatory(calendar,time,True)+9*self._burst(calendar,time,True))/calendar.norm
        return result
    def inventory(self,calendar,time,name):
        if time<=0:return 0.
        result=float(self.residence_primitive(time))
        if name=='hot_bursts':result=(result+self._oscillatory(calendar,time,False)+9*self._burst(calendar,time,False))/calendar.norm
        return result
