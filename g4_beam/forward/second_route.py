"""Independent equations/quadrature route. Deliberately never imports model or fits.
Run from any cwd: python second_route.py. Writes second_route.csv/json.
Adaptive scipy.quad in lab energy, rather than a cumulative log-grid trapezoid.
"""
from pathlib import Path
import csv
import json
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

HERE = Path(__file__).resolve().parent
HOST = {'Pd': (46, 106.42, 5.238, 6.80e22, .159),
        'Ti': (22, 47.867, 4.858, 5.67e22, .142)}
MAT = {'Pd': ('Pd', .1), 'PdD': ('Pd', .7), 'Ti': ('Ti', 1.), 'TiD': ('Ti', 1.7)}
IV = np.array([38,68,94,110,118,123,130,134,136,138,140,145,152,156,160,165,171,176,180,182,192,207,217,222,228,238,251,263,275,283])
GX, GW = np.polynomial.hermite_e.hermegauss(5)
GW /= np.sqrt(2*np.pi)
COMMON_REACTION = False


def cross(e, u, channel='p'):
    a = (55576.,210.54,-.032638,1.4987e-6,1.8181e-10) if channel == 'p' else (53701.,330.27,-.12706,2.9327e-5,-2.5151e-9)
    sf = sum(c*e**i for i,c in enumerate(a))
    if COMMON_REACTION and channel == 'p':
        # Independently evaluate literal Eq17 in SI, then J m² -> keV mb.
        joule_per_keV=1.602176634e-16
        mu=937806.5*joule_per_keV/299792458.**2
        sr_si=np.sqrt(2*986*joule_per_keV)*1.054571817e-34**3*(.040*1.602176634e-19)/(7e-15*mu**1.5*(e*joule_per_keV)**2)
        sr=sr_si/(joule_per_keV*1e-31)
        sf += sr+2*np.sqrt(sf*sr/3)*np.cos(115*np.pi/180)
    return sf / np.sqrt(e*(e+u))*np.exp(-31.397/np.sqrt(e+u))


def nuclear(e, z, mass):
    a = 1+z**.23
    reduced = 32.53*mass*e/(z*(2.01410+mass)*a)
    universal = np.log1p(1.1383*reduced)/(2*(reduced+.01321*reduced**.21226+.19593*np.sqrt(reduced))) if reduced <= 30 else np.log(reduced)/(2*reduced)
    return 8.462*z*2.01410*universal/((2.01410+mass)*a)


def stop(e, material):
    h,x = MAT[material]; z,m,a,_,_ = HOST[h]
    return (a+x*1.254)*np.sqrt(e*1.007276/2.01410)+nuclear(e,z,m)+x*nuclear(e,1,2.01410)


def yield_ion(e0, u, material, channel='p', surface=None, thickness=None):
    if e0 <= .1:
        return 0.
    h,x = MAT[material]; z,m,_,n,dv = HOST[h]
    f = lambda e,uu: x*cross(e/2,uu,channel)/stop(e,material)*1e-9
    if surface is None:
        return quad(lambda e:f(e,u), .1,e0,epsabs=1e-25,epsrel=2e-8)[0]
    # First angular moment: d<cos(theta)>/ds = -n sigma_transport <cos(theta)>.
    # Integrate the resulting projected displacement independently of Target.Q/Lambda.
    gamma = 4*2.01410*m/(2.01410+m)**2
    nh = n/(1+dv*x)
    def depth(e):
        def integrand(t):
            atten = quad(lambda v:2*nuclear(v,z,m)/(gamma*v*stop(v,material)),t,e0,epsabs=1e-10)[0]
            return np.exp(-atten)*1e25/(nh*stop(t,material))
        return quad(integrand,e,e0,epsabs=1e-10,epsrel=2e-7)[0]
    if depth(.1) <= thickness:
        return yield_ion(e0,surface,material,channel)
    ed = brentq(lambda e:depth(e)-thickness,.1,e0,xtol=1e-11)
    return quad(lambda e:f(e,u),.1,ed,epsabs=1e-25)[0]+quad(lambda e:f(e,surface),ed,e0,epsabs=1e-25)[0]


def relative(energies, material, physics, p, meas):
    def yi(e,ch):
        u = p['Ue']
        if physics == 'sigmoid_tt':
            u += (p['Ue_s']-u)/(1+np.exp((e/2-p['E0'])/.5))
        return yield_ion(e,u,material,ch,p.get('Ue_s') if physics == 'layer' else None,p.get('d'))
    def ion(e,ch):
        return sum(w*(.85*yi(2*e*(1+.05/2.3548*g),ch)+.30*yi(e*(1+.05/2.3548*g),ch)) for g,w in zip(GX,GW))
    def plotted(e,ch):
        ic = np.interp(2*e,np.arange(.5,15.01,.5),IV)
        c = (ic/140)**meas.get('b_hi',0) if e >= 2.75 else (ic/138)**meas.get('b_lo',0)
        it = ic/c
        flux = meas.get('f_n',0)*(it if meas.get('delta',0)==0 else 138)
        return (it*ion(e,ch)+flux*yi(meas.get('E_n',5),ch))/ic**(1+meas.get('k',0))
    ref = plotted(6.5,'n'); stitch = plotted(4.5,'n')/plotted(4.5,'p')
    return np.array([(plotted(e,'p')*stitch*meas.get('s',1) if e<=4 else plotted(e,'n'))/ref for e in energies])


def main():
    raw = list(csv.DictReader(l for l in open(HERE/'errbars.csv') if not l.startswith('#')))
    data = {m:np.array([[float(r[k]) for k in ('E_cm','R','sigma_lnR')] for r in raw if r['material']==m and float(r['E_cm'])<6.5]) for m in MAT}
    rows=[]
    for mat,u,cxs,ctt in [('Pd',0,1.92,2.15),('Pd',.8,1.20,1.41),('Ti',.3,1.66,1.88),('TiD',.9,1.10,1.32)]:
        rx = brentq(lambda e:np.log(cross(e,u)/cross(6.5,u))+4*np.log(10),.01,6)
        approx = lambda e:quad(lambda t:cross(t,u)/np.sqrt(t),.0001,e,epsabs=1e-20)[0]
        rt = brentq(lambda e:np.log(approx(e)/approx(6.5))+4*np.log(10),.01,6)
        proper = brentq(lambda e:np.log(yield_ion(2*e,u,mat)/yield_ion(13,u,mat))+4*np.log(10),.06,6)
        rows += [dict(quantity='crossing xs '+mat+' U='+str(u),independent=rx,reference=cxs,tolerance=.005,meaning='exact R=1e-4; earlier reference was a raster row above axis'),dict(quantity='crossing sqrt-stop '+mat+' U='+str(u),independent=rt,reference=ctt,tolerance=.005,meaning='approximate stopping, exact R=1e-4'),dict(quantity='crossing ICRU-ZBL '+mat+' U='+str(u),independent=proper,reference=None,tolerance=.005,meaning='single deuteron, actual material stopping, no beam average')]
    for k in (0,1):
        def resid(v):
            out=[]
            for i,m in enumerate(MAT):
                a=data[m];a=a[a[:,0]>=2.5]
                r=relative(a[:,0],m,'const',{'Ue':v[i]},{'k':k,'s':np.exp(v[4])})
                out.extend((np.log(r)-np.log(a[:,1]))/a[:,2])
            return out
        fit=least_squares(resid,[1.27,2.2,.99,1.34,0] if not k else [.7,1.41,.45,.84,0],bounds=([0,0,0,0,np.log(.3)],[15,15,15,15,np.log(3)]),diff_step=1e-4,ftol=1e-9,xtol=1e-9,gtol=1e-7)
        for m,expected in [('Pd',1.27 if not k else .70),('TiD',1.34 if not k else .84)]:
            rows.append(dict(quantity='bulk '+m+' M'+str(k),independent=float(fit.x[list(MAT).index(m)]),reference=expected,tolerance=.01,meaning='four-material shared proton/neutron stitch; adaptive quadrature'))
    old=json.load(open(HERE/'second_route_reference.json'))
    for cand in ('none','neutral d1','sigmoid_tt','layer wide bounds'):
        rr=min((r for r in old if r['cand']==cand),key=lambda r:r['chi2'])
        d=rr['params']; physics={'none':'const','neutral d1':'const','sigmoid_tt':'sigmoid_tt','layer wide bounds':'layer'}[cand]
        p={'Ue':d['Ue_Pd']}
        for a,b in [('Ues_Pd','Ue_s'),('E0','E0'),('d','d')]:
            if a in d:p[b]=d[a]
        mn=rr['meas'].split()[0]; presets={'M0':(0,0,0),'M1':(1,0,0),'M2':(0,1,0),'M3':(1,1,0),'M4':(0,1,1),'M5':(1,1,1)}
        meas=dict(zip(('k','b_hi','b_lo'),presets[mn]));meas['s']=d['s']
        if cand=='neutral d1':meas.update(f_n=d['f_n'],E_n=d['E_n'],delta=1)
        a=data['Pd'];r=relative(a[:,0],'Pd',physics,p,meas)
        chi=float(np.sum(((np.log(r)-np.log(a[:,1]))/a[:,2])**2))
        rows.append(dict(quantity='Pd frozen chi2 '+cand+' '+mn,independent=chi,reference=rr['per_mat']['Pd'],tolerance=.05,meaning='independent integral at exact saved joint-fit parameters; not independent parameter optimum'))
    for r in rows:
        r['difference']=None if r['reference'] is None else r['independent']-r['reference']
    with open(HERE/'second_route.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    json.dump(rows,open(HERE/'second_route.json','w'),indent=2)
    for r in rows:print(r,flush=True)

def common_comparison():
    """Changed common-reaction consumer only; reuse the unaffected Eq7 table."""
    global COMMON_REACTION
    COMMON_REACTION=True
    raw=list(csv.DictReader(l for l in open(HERE/'errbars.csv') if not l.startswith('#')))
    data={m:np.array([[float(r[k]) for k in ('E_cm','R','sigma_lnR')] for r in raw if r['material']==m and float(r['E_cm'])<6.5]) for m in MAT}
    saved=json.load(open(HERE/'fits.json'))
    rows=[]
    for cand in ('none','neutral d1','sigmoid_tt','layer'):
        r=min((r for r in saved if r['cand']==cand),key=lambda r:r['chi2'])
        d=r['params'];p={'Ue':d['Ue_Pd']}
        for a,b in [('Ues_Pd','Ue_s'),('E0','E0'),('d','d')]:
            if a in d:p[b]=d[a]
        mn=r['meas'].split()[0]
        meas=dict(zip(('k','b_hi','b_lo'),{'M0':(0,0,0),'M1':(1,0,0),'M2':(0,1,0),'M3':(1,1,0),'M4':(0,1,1),'M5':(1,1,1)}[mn]));meas['s']=d['s']
        if cand=='neutral d1':meas.update(f_n=d['f_n'],E_n=d['E_n'],delta=1)
        a=data['Pd'];pred=relative(a[:,0],'Pd',{'none':'const','neutral d1':'const'}.get(cand,cand),p,meas)
        chi=float(np.sum((np.log(pred/a[:,1])/a[:,2])**2))
        rows.append(dict(quantity='common fixed reaction Pd chi2 '+cand+' '+mn,independent=chi,reference=r['per_mat']['Pd'],tolerance=.05))
    for k,mn in [(0,'M0 as stated'),(1,'M1 fig x 1/I')]:
        r=next(r for r in saved if r['cand']=='bulk 2.5-6.5' and r['meas']==mn)
        initial=[r['params']['Ue_'+m] for m in MAT]+[np.log(r['params']['s'])]
        def residual(v):
            out=[]
            for i,m in enumerate(MAT):
                a=data[m];a=a[a[:,0]>=2.5]
                pred=relative(a[:,0],m,'const',{'Ue':v[i]},{'k':k,'s':np.exp(v[4])})
                out.extend(np.log(pred/a[:,1])/a[:,2])
            return out
        optimum=least_squares(residual,initial,bounds=([0,0,0,0,np.log(.3)],[15,15,15,15,np.log(3)]),diff_step=1e-4,ftol=2e-8,xtol=2e-8,gtol=2e-7)
        for mat in ('Pd','TiD'):
            rows.append(dict(quantity='common fixed reaction bulk '+mat+' '+mn,independent=float(optimum.x[list(MAT).index(mat)]),reference=r['params']['Ue_'+mat],tolerance=.01))
    for row in rows:
        row['difference']=row['independent']-row['reference']
        if abs(row['difference'])>row['tolerance']:
            raise AssertionError(row)
    with open(HERE/'second_route_common.csv','w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for row in rows:print(row,flush=True)



if __name__=='__main__':
    import sys
    common_comparison() if '--common' in sys.argv else main()
