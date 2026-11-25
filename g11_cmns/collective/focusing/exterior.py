"""Artificial C1 exteriors of the unchanged accepted on-domain PdD energy."""
import numpy as np
from numba import njit
from scipy.spatial import cKDTree
import engine as E

RMIN = 1.28782409436586
K4 = np.asarray(E.WELL.K4)
W0, S0 = float(E.QE[-1]), float(E.QS[-1])
# W'' just inside the final native Hermite endpoint; stiff is no milder than 4 Jnn.
LAST = len(E.QR)-2
TEND = E.QR[-1]-E.QR[-2]
W_SOFT_CURVATURE = max(2*E.JN, .1)
W_STIFF_CURVATURE = max(6*E.PC[0,LAST]*TEND+2*E.PC[1,LAST], 4*E.JN, 4*W_SOFT_CURVATURE)


@njit(cache=True)
def energy_force(x, force, edges, weights, near, directions, stiff):
    # Same accepted local operator; only the explicit exterior changes below.
    energy = E.local(x,force,False)
    if stiff:
        for i in range(len(x)):
            r2 = np.dot(x[i],x[i]); r = np.sqrt(r2)
            if r<=E.RMAX:
                continue
            z = (x[i,0]**2/r2,x[i,1]**2/r2,x[i,2]**2/r2)
            aq = z[0]*z[1]+z[1]*z[2]+z[2]*z[0]
            ap = z[0]*z[1]*z[2]
            aw = (1-4*aq+9*ap,4*aq-36*ap,27*ap)
            dr = r-E.RMAX
            ed = 1.5*E.RMAX**2*dr**2+E.RMAX*dr**3+.25*dr**4
            sd = 3*E.RMAX**2*dr+3*E.RMAX*dr**2+dr**3
            ee = (K4[0]*ed,K4[1]*ed,K4[2]*ed)
            ss = (K4[0]*sd,K4[1]*sd,K4[2]*sd)
            radial = 0.
            for a in range(3):
                energy += aw[a]*ee[a]
                radial += aw[a]*ss[a]
            deq = 4*(ee[1]-ee[0]); dep = 9*ee[0]-36*ee[1]+27*ee[2]
            for a in range(3):
                dq = 2*x[i,a]/r2*(1-z[a]-2*aq)
                dp = 2*x[i,a]/r2*(z[(a+1)%3]*z[(a+2)%3]-3*ap)
                force[i,a] -= x[i,a]/r*radial+deq*dq+dep*dp
    for b in range(len(edges)):
        i,j = edges[b,0],edges[b,1]
        w = weights[b]
        for a in range(3):
            energy -= w*x[i,a]*x[j,a]
            force[i,a] += w*x[j,a]
            force[j,a] += w*x[i,a]
    curvature = W_STIFF_CURVATURE if stiff else W_SOFT_CURVATURE
    for b in range(len(near)):
        i,j = near[b,0],near[b,1]
        u,v = 0.,0.
        for a in range(3):
            u += directions[b,a]*x[i,a]
            v -= directions[b,a]*x[j,a]
        if u>0 and v>0:
            if u*v<=E.RMAX**2:
                e,d = E.pair_correction(u*v)
            else:
                q = np.sqrt(u*v); dr = q-E.RMAX
                e = W0+S0*dr+.5*curvature*dr*dr-E.JN*u*v
                d = (S0+curvature*dr)/(2*q)-E.JN
            energy += e
            for a in range(3):
                force[i,a] -= d*v*directions[b,a]
                force[j,a] += d*u*directions[b,a]
    return energy


@njit(cache=True)
def evolve(x0,v0,edges,weights,near,directions,bonds,target_bond,ti,tj,dt,steps,stiff,capture):
    """Lift only support stopping; keep a distinct declared model-failure stop."""
    x,v = x0.copy(),v0.copy()
    f,fnew = np.empty_like(x),np.empty_like(x)
    potential = energy_force(x,f,edges,weights,near,directions,stiff)
    initial_kinetic = 0.
    for i in range(len(v)):
        for a in range(3):
            initial_kinetic += v[i,a]**2
    initial = potential+.5*E.MASS*initial_kinetic
    trace = np.empty((steps+1,10))
    # Full simultaneous vectors, never reconstructed from scalar radii.
    states = np.empty((len(capture),len(x),3))
    velocities = np.empty_like(states)
    first_exit,exit_kind,branch_kind = -1.,0,0
    count = 0
    for step in range(steps+1):
        maxsite,minbond,whole,minline,kinetic = 0.,1e30,0.,1e30,0.
        for i in range(len(x)):
            maxsite = max(maxsite,np.sqrt(np.dot(x[i],x[i])))
            for a in range(3):
                kinetic += v[i,a]**2
        for b in range(len(near)):
            i,j = near[b,0],near[b,1]
            rr,vv,line = 0.,0.,0.
            for a in range(3):
                dr = bonds[b,a]+x[j,a]-x[i,a]
                rr += dr*dr
                vv += (v[i,a]-v[j,a])**2
                line += dr*directions[b,a]
            minbond = min(minbond,np.sqrt(rr)); minline = min(minline,line)
            whole = max(whole,.25*E.MASS*vv)
        rr,dv2,rad = 0.,0.,0.
        for a in range(3):
            dr = target_bond[a]+x[tj,a]-x[ti,a]
            rr += dr*dr; dv2 += (v[ti,a]-v[tj,a])**2
            rad += (v[ti,a]-v[tj,a])*dr
        radius = np.sqrt(rr); closing = rad/radius
        total = potential+.5*E.MASS*kinetic
        trace[step] = (step*dt,.25*E.MASS*dv2,.25*E.MASS*closing**2,radius,
                       maxsite,minbond,total,potential,whole,closing)
        count = step+1
        for k in range(len(capture)):
            if step==capture[k]:
                states[k] = x; velocities[k] = v
        if first_exit<0 and (maxsite>E.RMAX or minbond<RMIN):
            first_exit,exit_kind = step*dt,1 if maxsite>E.RMAX else 2
        # These are numerical/model stops, not physical material thresholds.
        if not np.isfinite(total) or abs(total-initial)>.01*initial:
            branch_kind = 3
        elif minbond<.3:
            branch_kind = 1
        elif minline<=0:
            branch_kind = 2
        if branch_kind or step==steps:
            break
        for i in range(len(x)):
            for a in range(3):
                v[i,a] += .5*dt/E.MASS*f[i,a]
                x[i,a] += .5*dt*v[i,a]
                x[i,a] += .5*dt*v[i,a]
        potential = energy_force(x,fnew,edges,weights,near,directions,stiff)
        for i in range(len(x)):
            for a in range(3):
                v[i,a] += .5*dt/E.MASS*fnew[i,a]
        f,fnew = fnew,f
    return x,v,trace[:count],first_exit,exit_kind,branch_kind,states,velocities


def trajectory(lattice,q,v,dt=.25,horizon=2400.,stiff=False,capture=()):
    return evolve(q,v,lattice.edges,lattice.weights,lattice.near,lattice.directions,
                  lattice.bonds,lattice.target_bond,*lattice.target,dt,int(round(horizon/dt)),
                  stiff,np.asarray(capture,dtype=np.int64))


def constants():
    return dict(site_soft='E(.8,n)+S(.8,n)*delta+K*delta^2/2 with unchanged angular weights',
        site_stiff='accepted wide_quartic: soft + K4_ray*(1.5*.8^2*delta^2+.8*delta^3+delta^4/4)',
        pair='W(.8)+Wprime(.8)*delta+curvature*delta^2/2; unchanged q=sqrt(max(n.xi,0)max(-n.xj,0))',
        W0_eV=W0,Wprime_eV_A=S0,W_soft_curvature_eV_A2=float(W_SOFT_CURVATURE),
        W_stiff_curvature_eV_A2=float(W_STIFF_CURVATURE),Jnn_eV_A2=E.JN,
        model_failure_stops='actual NN R<.3 A; reference-axis NN lineage <=0; nonfinite E or >1% energy drift',
        support='all original crossings logged, not stopped; exterior results are outside computed electronic range')


def state_geometry(lattice,x,time_fs):
    radius = np.linalg.norm(x,axis=1)
    imax = int(np.argmax(radius))
    real_bonds = lattice.bonds+x[lattice.near[:,1]]-x[lattice.near[:,0]]
    distances = np.linalg.norm(real_bonds,axis=1)
    bmin = int(np.argmin(distances))
    def atom(i):
        return dict(row=int(i),original_D_index=int(E.BASE['original_indices'][i%32]),
            tile_translation_A=(lattice.positions[i]-E.BASE['positions'][i%32]).tolist(),
            equilibrium_position_A=lattice.positions[i].tolist(), displacement_A=x[i].tolist(),
            radius_A=float(radius[i]),
            normalized_crystal_direction=(x[i]/radius[i]).tolist() if radius[i] else None)
    def pair(i,j,bond,axis):
        vector = bond+x[j]-x[i]
        u,v = max(float(axis@x[i]),0.),max(float(-axis@x[j]),0.)
        return dict(first=atom(i),second=atom(j),reference_bond_A=bond.tolist(),
            actual_bond_vector_A=vector.tolist(),actual_R_A=float(np.linalg.norm(vector)),
            effective_symmetric_q_A=float(np.sqrt(u*v)),inward_projections_A=[u,v])
    wrapped = (lattice.positions+x)%lattice.box
    distance,neighbour = cKDTree(wrapped,boxsize=lattice.box).query(wrapped,k=2)
    self_first = neighbour[:,0]==np.arange(lattice.n)
    closest = np.where(self_first,neighbour[:,1],neighbour[:,0])
    actual_distances = np.where(self_first,distance[:,1],distance[:,0])
    ai = int(np.argmin(actual_distances)); aj = int(closest[ai])
    vector = wrapped[aj]-wrapped[ai]
    vector -= lattice.box*np.rint(vector/lattice.box)
    is_nn = bool(np.any(((lattice.near[:,0]==ai)&(lattice.near[:,1]==aj))|
                        ((lattice.near[:,0]==aj)&(lattice.near[:,1]==ai))))
    any_pair = dict(first=atom(ai),second=atom(aj),actual_bond_vector_A=vector.tolist(),
        actual_R_A=float(np.linalg.norm(vector)),is_original_NN=is_nn,
        scope='minimum over every periodic D pair at this simultaneous snapshot, not a full-time all-pairs bound')
    return dict(time_fs=float(time_fs),largest_site=atom(imax),
        minimum_NN_pair=pair(*lattice.near[bmin],lattice.bonds[bmin],lattice.directions[bmin]),
        target_pair=pair(*lattice.target,lattice.target_bond,lattice.target_axis),
        minimum_any_periodic_D_pair=any_pair)


def campaign(data,traces):
    from continuation import search,save
    if 'fine_N2048_25meV' not in data['cases']:
        raise RuntimeError('ordered fine-step reselection must be complete before exterior runs')
    lattice = E.Lattice(8)
    data['exterior_models'] = constants()
    states_path = E.HERE/'continuation_states.npz'
    states = dict(np.load(states_path)) if states_path.exists() else {}
    for energy in (.0048,.025):
        for stiff in (False,True):
            model = 'stiff' if stiff else 'soft'
            key = f'outside_{model}_{energy*1000:g}meV'
            selected,launches = search(lattice,energy,.25,key,data,traces,
                evolve=lambda q,v:trajectory(lattice,q,v,.25,2400.,stiff))
            row = selected[0]
            q,v = launches[row['candidate']]
            trace = traces[key+'__'+row['candidate']]
            peak = row['peak_index']
            last = len(trace)-1-(0 if row['model_valid_full_horizon'] else 1)
            captures = {peak,max(0,peak-40),min(last,peak+40),
                        int(np.argmax(trace[:last+1,4])),int(np.argmin(trace[:last+1,5]))}
            if row['first_site_or_bond_exit_fs'] is not None:
                captures.add(int(round(row['first_site_or_bond_exit_fs']/.25)))
            captures = sorted(k for k in captures if k<=last)
            result = trajectory(lattice,q,v,.25,2400.,stiff,captures)
            # Replays retain simultaneous displacement/velocity vectors only for
            # the selected peak, neighbourhood, extrema and first original exit.
            states[key+'__steps'] = np.array(captures)
            states[key+'__x'] = result[6]
            states[key+'__v'] = result[7]
            states[key+'__launch_q'] = q
            states[key+'__launch_v'] = v
            geometry = [state_geometry(lattice,x,.25*k) for k,x in zip(captures,result[6])]
            row['simultaneous_geometry'] = geometry
            row['largest_radius_over_all_candidates_A'] = max(
                r['max_site_radius_A_before_exit'] for r in data['cases'][key]['candidates'])
            row['smallest_NN_R_over_all_candidates_A'] = min(
                r['min_bond_A_before_exit'] for r in data['cases'][key]['candidates'])
            row['selected_replay_trace_equal'] = bool(np.array_equal(trace,result[2]))
            extrema = []
            candidates = data['cases'][key]['candidates']
            for event,extreme in (
                ('largest_site_over_search',max(candidates,key=lambda r:r['max_site_radius_A_before_exit'])),
                ('minimum_NN_over_search',min(candidates,key=lambda r:r['min_bond_A_before_exit']))):
                raw = traces[key+'__'+extreme['candidate']]
                stop = len(raw)-(0 if extreme['model_valid_full_horizon'] else 1)
                index = int(np.argmax(raw[:stop,4]) if event.startswith('largest') else np.argmin(raw[:stop,5]))
                if extreme is row and index in captures:
                    slot = captures.index(index)
                    xx,vv = result[6][slot],result[7][slot]
                else:
                    qx,vx = launches[extreme['candidate']]
                    replay = trajectory(lattice,qx,vx,.25,index*.25,stiff,[index])
                    xx,vv = replay[6][0],replay[7][0]
                states[key+'__'+event+'__x'] = xx
                states[key+'__'+event+'__v'] = vv
                extrema.append(dict(event=event,candidate=extreme['candidate'],
                    geometry=state_geometry(lattice,xx,index*.25)))
            data['cases'][key]['search_extrema'] = extrema
            np.savez_compressed(states_path,**states)
            save(data,traces)
    crossing_geometry(data,traces)


def crossing_geometry(data,traces):
    """Retain observed off-ray crossing states, not invented symmetric rays."""
    import patterns as P
    from continuation import save
    lattice = E.Lattice(8)
    path = E.HERE/'continuation_states.npz'
    states = dict(np.load(path))
    for key,case in data['cases'].items():
        if not key.startswith('outside_'):
            continue
        needed = [r for r in case['candidates'] if r['first_site_or_bond_exit_fs'] is not None
                  and 'outer_state_geometry' not in r]
        if not needed:
            continue
        speed = np.sqrt(2*case['energy_per_D_eV']/E.MASS)
        base = list(P.velocity_candidates('direct',lattice.positions,lattice.box,lattice.target,
            lattice.target_axis,speed,lattice.predicted_kernel,focus_times_fs=case['proposed_focus_times_fs']))
        launch = {r['name']:(r['q'],r['v']) for r in base}
        base_rows = [next(r for r in case['candidates'] if r['candidate']==b['name']) for b in base]
        best = max(base_rows,key=lambda r:r['peak_relative_kinetic_eV'])
        q0,v0 = launch[best['candidate']]
        for k,v in enumerate(P.direct_refinements(v0,lattice.target,lattice.target_axis)):
            launch[f'direct_phase_{k}'] = q0,v
        other = next(r for r in base_rows if r is not best)
        vb = launch[other['candidate']][1]
        for c in (-.5,-.25,.25,.5):
            launch[f'direct_mixture_{c:+.2f}'] = q0,P.capped_direct_velocity(
                v0+c*vb,lattice.target,speed,lattice.target_axis)
        largest = case['search_extrema'][0]
        for row in needed:
            raw = traces[key+'__'+row['candidate']]
            index = int(np.argmax(raw[:,4]))
            if row['candidate']==largest['candidate']:
                xx = states[key+'__largest_site_over_search__x']
                vv = states[key+'__largest_site_over_search__v']
            else:
                q,v = launch[row['candidate']]
                replay = trajectory(lattice,q,v,case['dt_fs'],index*case['dt_fs'],
                                    '_stiff_' in key,[index])
                xx,vv = replay[6][0],replay[7][0]
            states[key+'__'+row['candidate']+'__outer_x'] = xx
            states[key+'__'+row['candidate']+'__outer_v'] = vv
            row['outer_state_geometry'] = state_geometry(lattice,xx,index*case['dt_fs'])
    data['geometry_readout'] = 'same-candidate .25fs snapshot acquisition at every crossing candidate radius maximum; no new optimization or timestep refinement'
    np.savez_compressed(path,**states)
    save(data,traces)
