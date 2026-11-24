! Standalone driver linked against NuHamil eb65165; no changes to upstream.
program g11_intrinsic
  use ABodyManager, only: solve_A_body_hamil, get_3_body_jacobi_op_isospin
  use ABodyJacSpaceIso, only: ABodyJacIsoChan
  use ABodyJacOpsIso, only: ABodyJacOpChanIso
  use ClassSys, only: sys
  use LinAlgLib, only: DMat, DVec
  use MyLibrary, only: tjs, geometry_part, red_r_l, init_dtrinomial, fin_dtrinomial, m_proton, m_neutron
  use NuHamilInput, only: SetInputParameters, InputParameters, params
  use OperatorDefinitions, only: SetFrequency
  use Profiler, only: timer
  use ThreeBodyJacobiSpace, only: ThreeBodyJacIsoSpace
  use ThreeBodyJacOpsIso, only: ThreeBodyJacOpIso
  use TwoBodyRelativeSpace, only: TwoBodyRelSpaceIsoHOBasis, TwoBodyRelChanIsoHOBasis, HarmonicOscillator
  use TwoBodyRelOpsIso, only: TwoBodyRelOpIso
  implicit none
  type(sys) :: s
  type(InputParameters) :: input
  type(TwoBodyRelSpaceIsoHOBasis), target :: rel
  type(TwoBodyRelChanIsoHOBasis), pointer :: b, k
  type(HarmonicOscillator), pointer :: ob, ok
  type(TwoBodyRelOpIso) :: pair
  type(ThreeBodyJacIsoSpace), target :: jac
  type(ThreeBodyJacOpIso) :: triple
  type(ABodyJacIsoChan), target :: even, odd
  type(ABodyJacOpChanIso) :: up, down
  type(DVec) :: e_even, e_odd
  type(DMat) :: v_even, v_odd
  integer :: cb, ck, ib, ik, jb, jk, tb, tk, u, pair_b, pair_k, ho0
  real(8) :: mu, control, ref, ho_strength, ho_ref

  call timer%init()
  call SetInputParameters()
  if(params%particle_rank /= 4 .or. params%renorm%val /= 'bare') error stop 'Expected bare A=4 input'
  call SetFrequency(params%hw)
  call init_dtrinomial()
  mu = m_proton*m_neutron/(m_proton+m_neutron)  ! native equal NN HO masses, units MeV/c^2

  input=params
  input%fn_fewbody_wave_function = s%str('even')
  call solve_A_body_hamil(input, 0, +1, 0, 0, e_even, v_even)
  input%fn_fewbody_wave_function = s%str('odd')
  call solve_A_body_hamil(input, 2, -1, 0, 0, e_odd, v_odd)

  ! O=sum_i sigma_i cross (pi_i c) = (2/A) sum_{i<j} (sigma_i-sigma_j) cross (p_ij c).
  ! p_ij=(p_i-p_j)/2. Native pair matrices are reduced in both J and T.
  ! TDM_mag_s is used ONLY as native metadata (J=1, parity=-1, T=0); its built-in matrix is never called.
  call rel%init(params%hw, params%NAmax, params%Jmax2)
  call pair%init(rel, rel, s%str('TDM_mag_s'))
  do ck=1,rel%GetNumberChannels()
    k => rel%GetChannel(ck)
    jk=k%GetJ(); tk=k%GetT()
    do cb=1,rel%GetNumberChannels()
      if(pair%MatCh(cb,ck)%is_zero) cycle
      b => rel%GetChannel(cb)
      jb=b%GetJ(); tb=b%GetT()
      if(tb /= tk) error stop 'Expected isoscalar block'
      do ik=1,k%GetNumberStates()
        ok => k%GetP(ik)
        do ib=1,b%GetNumberStates()
          ob => b%GetP(ib)
          pair%MatCh(cb,ck)%m(ib,ik)=sqrt(real(2*tk+1,8))*pair_reduced(jb,ob,jk,ok)
        end do
      end do
    end do
  end do

  ! Independent stretched pair control, in Cartesian coordinates:
  ! <0p,S=0,J=1,M=1|O_z(12)|0s,S=1,J=1,M=1> = sqrt(mu*hw/2) for A=4.
  pair_b=rel%GetIndex(1,-1,0)
  pair_k=rel%GetIndex(1,+1,0)
  ib=pair%MatCh(pair_b,pair_k)%rel_ch_bra%GetIndex(0,1,0)
  ik=pair%MatCh(pair_b,pair_k)%rel_ch_ket%GetIndex(0,0,1)
  control=pair%MatCh(pair_b,pair_k)%m(ib,ik)*geometry_part(2,2,2,2,0,2)
  ref=sqrt(mu*params%hw/2.d0)
  if(abs(control-ref)>1.d-12*ref) error stop 'Pair tensor sign or normalization'

  call jac%init(params%hw, 2*params%NAmax+3, params%NAmax, params%ramp, params%N3max, params%path_to_tmp_dir)
  call even%init(0,+1,0,params%NAmax,jac,params%path_to_tmp_dir)
  call odd%init(2,-1,0,params%NAmax,jac,params%path_to_tmp_dir)
  triple = get_3_body_jacobi_op_isospin(jac, pair, s%str('TDM_mag_s'))
  call up%init(odd, even, s%str('spin_cross_pc'))
  call down%init(even, odd, s%str('spin_cross_pc'))
  call up%set(triple)
  call down%set(triple)
  ! Native triple includes 3 equivalent pairs. A=4 has 6: multiply by 6/3.
  ! Reduce to the Cartesian z component (q=0) between M=0 states, and T=Tz=0.
  up%m=2.d0*geometry_part(2,2,0,0,0,0)*up%m
  down%m=2.d0*geometry_part(0,2,2,0,0,0)*down%m
  if(maxval(abs(up%m-transpose(down%m)))>1.d-9) error stop 'Projected operator orientation or Hermiticity'

  ! From the unique 0 hbar-omega spin-isospin singlet, the span of all odd
  ! 1 hbar-omega states exhausts ||O_z |0s^4>||^2 = 4 m hw = 8 mu hw.
  ho0=even%GetASIndex(0,1)
  ho_strength=sum(up%m(:,ho0)**2)
  ho_ref=8.d0*mu*params%hw
  if(abs(ho_strength/ho_ref-1.d0)>1.d-9) error stop 'Four-body oscillator closure'

  open(newunit=u,file='g11_matrices.bin',form='unformatted',access='stream',status='replace')
  write(u) size(e_even%v),size(e_odd%v)
  write(u) mu,control,ref, maxval(abs(up%m-transpose(down%m))),ho_strength,ho_ref
  write(u) e_even%v,e_odd%v,v_even%m,v_odd%m,up%m,down%m
  close(u)
  call fin_dtrinomial()
  call timer%fin()

contains
  real(8) function cg_coef(j1,m1,j2,m2,j,m) result(r)
    integer, intent(in) :: j1,m1,j2,m2,j,m  ! twice each quantum number
    r=(-1.d0)**((j1-j2+m)/2)*sqrt(real(j+1,8))*tjs(j1,j2,j,m1,m2,-m)
  end function cg_coef

  real(8) function pair_reduced(jb,b,jk,k) result(r)
    integer, intent(in) :: jb,jk
    type(HarmonicOscillator), intent(in) :: b,k
    integer :: nb,nk,lb,lk,sb,sk,mb,mk,q,ql,qs,mlk,mlb,msb,msk,dn
    real(8) :: qweight, orbital, spin, cb, ck, ang, radial
    r=0.d0
    nb=b%GetN(); nk=k%GetN(); lb=b%GetL(); lk=k%GetL(); sb=b%GetS(); sk=k%GetS()
    if(sb==sk) return
    ! p c = i sqrt(mu*hw) (N_b-N_k) (r/b_osc), coordinate HO convention.
    ! A cross B = -i sqrt(2) [A^1 tensor B^1]^1; i from p cancels.
    dn=2*(nb-nk)+lb-lk
    radial=red_r_l(nb,lb,nk,lk)*real(dn,8)*sqrt(mu*params%hw)*sqrt(2.d0)/2.d0
    if(radial==0.d0) return
    do mk=-jk,jk
      do q=-1,1
        mb=mk+q
        if(abs(mb)>jb) cycle
        qweight=geometry_part(2*jb,2,2*jk,2*mb,2*q,2*mk)
        if(abs(qweight)<1.d-14) cycle
        do msk=-sk,sk
          mlk=mk-msk
          if(abs(mlk)>lk) cycle
          ck=cg_coef(2*lk,2*mlk,2*sk,2*msk,2*jk,2*mk)
          do qs=-1,1
            ql=q-qs
            if(abs(ql)>1) cycle
            msb=msk+qs; mlb=mlk+ql
            if(abs(msb)>sb.or.abs(mlb)>lb) cycle
            spin=2.d0
            if(sb==0) spin=spin*(-1.d0)**qs
            cb=cg_coef(2*lb,2*mlb,2*sb,2*msb,2*jb,2*mb)
            orbital=geometry_part(2*lb,2,2*lk,2*mlb,2*ql,2*mlk)
            ang=cg_coef(2,2*qs,2,2*ql,2,2*q)
            r=r+qweight*ck*cb*spin*orbital*ang
          end do
        end do
      end do
    end do
    ! Sum of squared Wigner-Eckart geometrical weights over mb,mk,q is one.
    r=r*radial
  end function pair_reduced
end program g11_intrinsic
