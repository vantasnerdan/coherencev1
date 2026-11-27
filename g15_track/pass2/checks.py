"""Small witnesses at the repaired physical joins, not a repository test suite."""
import numpy as np
from observer import CurvatureObserver,discriminator
from mechanics import ForceHistory


def run(m,previous):
    result=dict(previous)
    result.update(m.positive_loss_check())
    result['curvature_observer']=discriminator()
    joins=[]
    for env in m.D0:
        n=m.AXIS[env];ia,ib=(0,4) if env=='bulk' else (0,3)
        active=np.zeros((3,1)+m.SHAPE,bool);active.reshape(3,-1)[:,:8 if env=='bulk' else 6]=True
        def state(x,y):
            q=np.zeros((1,1,3)+m.SHAPE);q.reshape(1,1,3,-1)[0,0,:,ia]=x*n;q.reshape(1,1,3,-1)[0,0,:,ib]=y*n
            return q
        def energy(x,y):return float(m.lattice_force(state(x,y),env,active[:1])[1].sum())
        u=.24 if env=='bulk' else .12
        c,limit=m.CURVES[env][0][0],m.CURVES[env][1]
        full=energy(u,-u)-energy(u,0)-energy(0,-u)+energy(0,0)
        native=float(m.curve_value(c,np.array(u),limit))
        F=m.lattice_force(state(u,-u),env,active[:1])[0].reshape(3,-1)
        FA=m.lattice_force(state(u,0),env,active[:1])[0].reshape(3,-1)
        FB=m.lattice_force(state(0,-u),env,active[:1])[0].reshape(3,-1)
        derivative=float(-(F[:,ia]-FA[:,ia])@n+(F[:,ib]-FB[:,ib])@n)
        expected=float(m.curve_value(c,np.array(u),limit,1))
        axes=[float(-x*y*m.connected_g(np.array((x-y)/2),env)) for x,y in ((u,0),(0,-u))]
        if abs(full-native)>1e-11 or abs(derivative-expected)>1e-11 or max(abs(z) for z in axes)>1e-15:
            raise ValueError('Native opposed-cut / one-body-axis identity failed')
        q=np.zeros((3,1,3)+m.SHAPE);q[0]=.03*np.sin(np.arange(24)).reshape((1,3)+m.SHAPE);q[1]=-q[0]
        q*=active[:,:,None]
        nsite=8 if env=='bulk' else 6;position=np.zeros((8,3));weights=np.r_[np.ones(nsite),np.zeros(8-nsite)]
        omega=np.array([np.sqrt(m.KG/m.M)])
        observer=CurvatureObserver(m,env,position,weights,np.ones((8,1),complex),omega,np.array([0.]))
        observer.add(0.,.5,q,active)
        measured=float((n@observer.buffer[-1][2][0,0].reshape(3,3)@n).real)
        direction=np.broadcast_to(n.reshape((1,1,3,1,1,1)),q[:1].shape)*active[:1,:,None]
        eps=1e-5
        def stiffness(q0):
            fp=m.lattice_force(q0+eps*direction,env,active[:1])[0]
            fm=m.lattice_force(q0-eps*direction,env,active[:1])[0]
            return -float(np.sum((fp-fm)*direction))/(2*eps)
        finite=(stiffness(q[:1])-stiffness(q[:1]*0))/(m.M*nsite*omega[0]**2)
        error=abs(measured-finite)
        if error>1e-7:raise ValueError('Actual Hamiltonian curvature / observer join failed')
        joins.append(dict(environment=env,native_cut_u_A=u,opposed_connected_eV=full,native_connected_eV=native,
            connected_force_derivative_eV_A=derivative,native_signed_derivative_eV_A=expected,
            one_body_axis_connected_eV=axes,curvature_observer=measured,finite_force_curvature=finite,
            curvature_difference=error,extension='-x*y*g((x-y)/2); smooth even opening is conditional, not native',
            bond='only tagged axial G3 cross term removed before native replacement',local_modes=observer.local_info))
    result['native_and_curvature_joins']=joins
    # The old second moment depended on the repaired connected Hamiltonian.
    pq=np.zeros((2,1,3)+m.SHAPE);pv=np.zeros_like(pq);pq[0,:,0]=.05;pq[1,:,0]=-.05
    pa=np.ones((2,1)+m.SHAPE,bool)
    for _ in range(200):
        pq,pv=m.hamiltonian_step(pq,pv,.5,'bulk',pa);pv-=.5*.1*m.K*pq/m.M
    result.update(even_prior_mean_A=float(pq.mean()),even_prior_second_moment_A2=float(np.mean(pq*pq)))
    s=m.p1.SCENARIOS[0];tr=m.p1.track(2.5,1);z=tr['timeline'][len(tr['timeline'])//2]['depth_A']
    history=ForceHistory(m,2.5,1,np.array([[2.,0.,z]]),'bulk',1.)
    t0=history.clock[0];delay=.74/np.sqrt(2*tr['initial_eV']/m.M)
    before=history.cumulative(np.array([-1.]))
    leading=history.cumulative(np.array([t0]))
    first_at_second=history.motion(np.array([t0+delay]),m.OMEGA,1/(s['tau_opt']*1e15))
    if abs(before).max()!=0 or abs(leading).max()==0:raise ValueError('Causal signed mechanical history failed')
    result['signed_mechanical_clock']=dict(z_A=z,b_A=2.,axial_constituent_delay_fs=delay,
        zero_before_entry_momentum_eV_fs_A=before[0].tolist(),leading_closest_approach_momentum_eV_fs_A=leading[0].tolist(),
        first_q_at_second_A=first_at_second[0][0].tolist(),first_v_at_second_A_fs=first_at_second[1][0].tolist(),
        candidate_transfer_eV=float(history.T[0]),strong_mark=bool(history.strong[0]),
        timing='piecewise accepted projectile z(t), force integral includes leading field; source ledger caps actual signed work')
    result['reuse']='unaffected zero/scalar launch and spatial-gradient/operator witnesses retained; nonzero-prior, native, loss and curvature joins exercised after repair'
    result['consumer_measure_discriminator']=measure_discriminator(m)
    return result

def measure_discriminator(m):
    from observables import ConsumerMeasures,inward_crossing
    def crossing(a,b):
        hit,root=inward_crossing(np.array([[a]],float),np.array([[b]],float),np.array([1.]))
        return int(hit.sum()),float(root[0,0,0])
    through=crossing([-2,0,0],[2,0,0]);tangent=crossing([-2,1,0],[2,1,0])
    endpoint=crossing([2,0,0],[1,0,0]);following=crossing([1,0,0],[.5,0,0])
    reentry=crossing([2,0,0],[-2,0,0])
    if (through[0],tangent[0],endpoint[0],following[0],reentry[0])!=(1,0,1,0,1) or abs(through[1]-.25)>1e-14:
        raise ValueError('Within-step entry/tangency/endpoint convention failed')
    L=3*len(m.BRANCHES);w=np.ones(8);pw=np.ones(1)
    measure=ConsumerMeasures(m,.25,0,8.,w,pw,np.array([2.4]),'bulk')
    velocity=np.zeros((L,8,3));r=np.broadcast_to(np.array([2.88782409436586,0,0]),(L,1,3))
    intervals=[.25,.75,2.5,1.25]
    for dt in intervals:measure.add(dt,velocity,np.zeros(L),r,np.zeros_like(r),np.ones((L,1),bool))
    actual=measure.moments[:,:,:,0].sum(axis=2);expected=8*sum(intervals)
    defect=float(abs(actual-expected).max())
    if defect>1e-12 or abs(measure.moments[:,:,:,1:]).max()!=0:raise ValueError('Physical exposure / exact quiet-energy anchor failed')
    return dict(actual_D_fs=float(actual[0,0]),required_D_fs=expected,exposure_defect_D_fs=defect,
        old_last_step_only_D_fs=8*intervals[-1],zero_source_higher_energy_moments=0.,
        chord_entry_fraction=through[1],tangency_entries=tangent[0],shared_endpoint_entries=endpoint[0]+following[0],
        reentry_entries=reentry[0],crossing_velocity='negative position secant / physical dt at the same interpolant root',
        accuracy='within-step linear-position interpolant, not continuum-exact; endpoint ties assigned once')
