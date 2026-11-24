"""Bidirectional finite shared-mode model: identical H in QuTiP and quimb.

Only the even leading-order coupling is physical/exact (zero). The nonzero
calculation uses the explicitly standalone odd witness. Cold, energy-resolved
secular jumps are dressed in this same H: a bare Lindblad jump is not used to
turn virtual nuclear admixture into energy-conserving upconversion.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import csv

import numpy as np
import scipy.sparse as sp
import qutip as qt
import quimb as qu
from scipy.constants import hbar, electron_volt

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
INPUT=ROOT/'g11_cmns/coupling/collective_input.json'
HBAR=hbar/electron_volt*1e12  # eV ps
KAPPA=2/4.136036943628514    # number decay ps^-1, saved amplitude decay4.136ps


def put(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def build(n,cutoff,w,gap,g):
    """H=hw n+Delta sum Pe+sum (g/sqrt(N)) Sy_j i(a†-a).

    One moving D per conditional pair, with the signed Pc and C phase of the
    coupling interface. The optical mode is normalized across N cells. No
    unnormalised Dicke sqrt(N) enhancement is inserted a second time.
    """
    dims=[2]*n+[cutoff]
    a=np.diag(np.sqrt(np.arange(1,cutoff)),1)
    occ=np.diag(np.arange(cutoff,dtype=float))
    pe=np.diag([0.,1.]);sy=np.array([[0,1j],[-1j,0]])
    def place(op,site):return qu.ikron(op,dims,site,sparse=False)
    number=np.asarray(place(occ,n));nuclear=sum(np.asarray(place(pe,j)) for j in range(n))
    momentum=np.asarray(place(1j*(a.T-a),n))
    vertices=[np.asarray(place(sy,j)) for j in range(n)]
    interaction=sum(g/np.sqrt(n)*v@momentum for v in vertices)
    h=w*number+gap*nuclear+interaction
    # Independently assembled tensor-product QuTiP Hamiltonian.
    def qplace(op,site):
        items=[qt.qeye(d) for d in dims];items[site]=qt.Qobj(op);return qt.tensor(items)
    qh=w*qplace(occ,n)+gap*sum(qplace(pe,j) for j in range(n))
    qh += sum(g/np.sqrt(n)*qplace(sy,j)*qplace(1j*(a.T-a),n) for j in range(n))
    if np.max(abs(h-qh.full()))>1e-10:raise AssertionError('Hamiltonians differ')
    return h,qh,number,nuclear,np.asarray(place(a+a.T,n)),vertices,dims


def run(n,occupation,direction,width_ratio,w,gap,g,physical_zero=False):
    cutoff=max(7,occupation+3)
    h,qh,number,nuclear,mode,vertices,dims=build(n,cutoff,w,gap,0. if physical_zero else g)
    eig,u=qu.eigh(h)
    qeig=np.asarray(qh.eigenenergies())
    dim=len(eig)
    e=eig-eig[0]
    # State preparation is a dressed eigenstate, not a sudden bare state with
    # a hidden high-energy spectral tail. Select via its stated bare overlap.
    bare=np.zeros(dim,complex)
    if direction=='up':bare[occupation]=1.
    else:
        for j in range(n):
            bits=[0]*n;bits[j]=1
            index=np.ravel_multi_index(tuple(bits+[occupation]),dims)
            bare[index]=1/np.sqrt(n)
    overlap=abs(u.conj().T@bare)**2
    initial=int(overlap.argmax());p0=np.eye(dim)[initial]
    mats=[u.conj().T@mode@u]+[u.conj().T@v@u for v in vertices]
    qgen=np.zeros((dim,dim));edges=[];jumps=[]
    # Source branching/loss is conditional, not imported from unrelated levels.
    # Gamma and particle rates share width_ratio*KAPPA with chosen .25/.75 split.
    for b in range(dim):
        for a in range(b):
            de=float(e[b]-e[a])
            if de<=1e-8:continue
            if de<4*w:
                rate=KAPPA*abs(mats[0][a,b])**2
                if rate>1e-20:edges.append((a,b,float(rate),'optical',de))
            if de>gap/2 and width_ratio:
                strength=sum(abs(m[a,b])**2 for m in mats[1:])
                for fraction,label in ((.25,'gamma'),(.75,'particle')):
                    rate=width_ratio*KAPPA*fraction*strength
                    if rate>1e-20:edges.append((a,b,float(rate),label,de))
    for a,b,rate,label,de in edges:
        qgen[a,b]+=rate;qgen[b,b]-=rate
        jump=sp.csr_matrix(([np.sqrt(rate)],([a],[b])),shape=(dim,dim))
        jumps.append(qt.Qobj(jump))
    times=np.linspace(0,5/KAPPA,81)
    # Initial rho commutes with H and every energy-resolved jump preserves that
    # subspace. H commutator is exactly zero there; interaction picture avoids
    # resolving 34-MeV phases during picosecond population decay.
    initial_rho=qt.Qobj(sp.csr_matrix(([1.],([initial],[initial])),shape=(dim,dim)))
    pop_ops=[qt.Qobj(sp.csr_matrix(([1.],([j],[j])),shape=(dim,dim))) for j in range(dim)]
    solved=qt.mesolve(qt.qzero(dim),initial_rho,times,jumps,e_ops=pop_ops,
        options={'rtol':1e-10,'atol':1e-13,'nsteps':100000})
    qt_p=np.asarray(solved.expect).T
    # Independent quimb exponential action on the same derived population
    # generator, augmented with residence integrals for exact jump counts.
    aug=sp.bmat([[sp.csr_matrix(qgen),None],[sp.eye(dim),sp.csr_matrix((dim,dim))]],format='csr')
    state=np.r_[p0,np.zeros(dim)]
    evolved=qu.expm_multiply(aug,state,start=0,stop=float(times[-1]),num=len(times))
    qp=np.asarray(evolved)[:,:dim].real;residence=np.asarray(evolved)[-1,dim:].real
    error=float(np.max(abs(qt_p-qp)))
    if error>2e-8:raise AssertionError(error)
    phonon_diag=np.real(np.diag(u.conj().T@number@u))
    virtual_diag=np.real(np.diag(u.conj().T@nuclear@u))
    high=(e>gap/2).astype(float)
    nprob=np.zeros((dim,cutoff))
    for j in range(dim):nprob[j]=np.sum(abs(u[:,j].reshape(*dims))**2,axis=tuple(range(n)))
    initial_energy=float(e[initial]);final_energy=float(qp[-1]@e)
    channel={k:dict(quanta=0.,energy_eV=0.,lines=[]) for k in ('optical','gamma','particle')}
    for a,b,rate,label,de in edges:
        count=float(rate*residence[b]);channel[label]['quanta']+=count;channel[label]['energy_eV']+=count*de
        if count>1e-16:channel[label]['lines'].append(dict(energy_eV=de,expected_quanta=count))
    loss=sum(x['energy_eV'] for x in channel.values())
    balance=final_energy+loss-initial_energy
    # Compare energy with source-scale floating arithmetic, not an unattainable
    # relative target on a near-zero optical difference.
    if abs(balance)>max(2e-6,2e-9*initial_energy):raise AssertionError(balance)
    return dict(N=n,initial_phonons=occupation,direction=direction,
        setup='physical_exact_even_leading_zero' if physical_zero else 'standalone_odd_witness_NOT_physical_DD',
        gap_eV=gap,phonon_eV=w,g_single_cell_eV=0. if physical_zero else g,
        g_each_cell_eV=0. if physical_zero else g/np.sqrt(n),
        signed_mode_extension='one moving D and one stationary pair member, per coupling JSON; every cell /sqrt(N), not a measured foil mode',
        loss_width_over_mode_number_decay=width_ratio,
        loss_status='cold dressed optical bath cutoff4hw; nuclear gamma25% particle75% are specified sensitivities, not measured widths',
        prepared_eigenstate=initial,bare_state_overlap=float(overlap[initial]),
        hilbert_dimension=dim,fock_cutoff=cutoff,time_ps=times.tolist(),
        nuclear_high_band_probability=(qp@high).tolist(),
        bare_nuclear_virtual_occupation=(qp@virtual_diag).tolist(),
        mean_phonons=(qp@phonon_diag).tolist(),
        initial_phonon_spectrum=(p0@nprob).tolist(),final_phonon_spectrum=(qp[-1]@nprob).tolist(),
        maximum_top_fock_probability=float(np.max(qp@nprob[:,-1])),
        initial_energy_eV=initial_energy,final_energy_eV=final_energy,emitted=channel,
        energy_balance_error_eV=balance,qutip_quimb_population_max_abs=error,
        eigenvalue_route_max_abs_eV=float(np.max(abs(eig-qeig))),
        largest_populated_nuclear_energy_band_eV=gap if np.max(qp@high)>1e-8 else 0.,
        emitted_optical_average_quanta_per_s=channel['optical']['quanta']/(times[-1]*1e-12),
        emitted_optical_average_eV_per_s=channel['optical']['energy_eV']/(times[-1]*1e-12),
        rate_semantics='finite prepared toy-state window average; not stationary fusion or measured phonon-production rate')


def main():
    started=time.perf_counter();HERE.mkdir(parents=True,exist_ok=True)
    data=json.loads(INPUT.read_text());inputs=data['collective_model_inputs']
    odd=inputs['standalone_odd_witness'];w=odd['mode']['phonon_eV'];gap=odd['gap_eV']
    g=odd['mode']['conditional_g_times_Pc_0to1_eV'][0]
    results=[]
    # This finite set discriminates N-normalization, occupation and loss; no
    # extrapolation to 10^8-10^9 phonons is performed.
    for n in (1,2,4):
        for occupation in (0,4):
            for direction in ('up','down'):
                for loss in (0.,1.):
                    results.append(run(n,occupation,direction,loss,w,gap,g))
    for direction in ('up','down'):
        results.append(run(2,4,direction,1.,w,gap,0.,physical_zero=True))
    scaling=[]
    for r in results:
        scaling.append({k:r[k] for k in ('N','initial_phonons','direction','setup','g_each_cell_eV',
            'loss_width_over_mode_number_decay','initial_energy_eV','energy_balance_error_eV',
            'qutip_quimb_population_max_abs','maximum_top_fock_probability') } |
            dict(initial_virtual_nuclear_occupation=r['bare_nuclear_virtual_occupation'][0],
                final_nuclear_high_band=r['nuclear_high_band_probability'][-1],
                optical_quanta=r['emitted']['optical']['quanta'],gamma_quanta=r['emitted']['gamma']['quanta'],
                particle_quanta=r['emitted']['particle']['quanta']))
    with (HERE/'scaling.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(scaling[0]));writer.writeheader();writer.writerows(scaling)
    put(HERE/'results.json',dict(runs=results,source_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),runtime_seconds=time.perf_counter()-started,
        versions=dict(qutip=qt.__version__,quimb=qu.__version__,numpy=np.__version__,python=sys.version.split()[0]),
        physical_unknowns=inputs['physical'],
        even_second_order=inputs['standalone_even_path'],
        omissions=['No Coulomb-normalized physical odd continuum or complete physical second-order coupling.',
            'N<=4, occupations0/4, cutoff7; no high-order many-phonon transfer scaling established.',
            'No measured participating mode/active-domain population, matched gamma/particle widths or bath spectrum.']))
    print(json.dumps(dict(runs=len(results),max_population_difference=max(r['qutip_quimb_population_max_abs'] for r in results),runtime_s=time.perf_counter()-started)))


if __name__=='__main__':main()
