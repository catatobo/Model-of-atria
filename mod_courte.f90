! This module implements the Courte Model
!
! 08-02-2011: Courte_A_P01 returns currents vector for output (EAH)
!
!------------------------------------------------------------------------------!
module mod_courte
!------------------------------------------------------------------------------!
  use mod_precision
  implicit none
  include 'CRN_parameters.inc'
  !.
  type t_crn
    sequence
    real(rp) :: m, h, j, oa, oi, ua, ui, xr, xs, d, f, fca, nai, ki, cai,  &
                caup, carel, u, ve, w

  end type t_crn
  !.
  public  :: get_parameter_CRN , ic_Courte, Courte_A_P01
  private :: put_param, f_cur, f_comp
  !.
  type, public:: t_prm
    private
    real(rp)  :: p_gna    ! 1	
    real(rp)  :: p_gk1    ! 2  
    real(rp)  :: p_gto    ! 3
    real(rp)  :: p_gkr    ! 4 
    real(rp)  :: p_gks    ! 5 
    real(rp)  :: p_gcal   ! 6 
    real(rp)  :: p_inakm  ! 7
    real(rp)  :: p_inacam ! 8
    real(rp)  :: p_gbca   ! 9
    real(rp)  :: p_gbna   ! 10
    real(rp)  :: p_ipcam  ! 11
    real(rp)  :: p_krel   ! 12
    real(rp)  :: p_iupm   ! 13
    real(rp)  :: p_gach   ! 14
    real(rp)  :: p_ACh    ! 15
  end type t_prm
!------------------------------------------------------------------------------!
!  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&  !
!------------------------------------------------------------------------------!
contains
!------------------------------------------------------------------------------!
function get_parameter_CRN (tcell) result (v_prm)
!------------------------------------------------------------------------------!
  implicit none
  integer(ip), intent(in) :: tcell
  real (rp)               :: v_prm(15)

  select case (tcell)
  case (16)
    v_prm = (/ p_gna, p_gk1, p_gto, p_gkr, p_gks, p_gcal, p_inakm, p_inacam, &
               p_gbca, p_gbna, p_ipcam, p_krel, p_iupm, p_gach, p_ACh /) 
  case default
    write (*,10) tcell; stop
  end select
  return
  10 format ('###.Specified Cell type is not defined:',I3)
end function get_parameter_CRN
!------------------------------------------------------------------------------!
subroutine put_param(v_prm,param)
!------------------------------------------------------------------------------!
  implicit none
  real(rp), intent(in)    :: v_prm(:)
  type (t_prm)            :: param
  param%p_gna      = v_prm(1)  ! 1	
  param%p_gk1      = v_prm(2)  ! 2  
  param%p_gto      = v_prm(3)  ! 3
  param%p_gkr      = v_prm(4)  ! 4 
  param%p_gks      = v_prm(5)  ! 5 
  param%p_gcal     = v_prm(6)  ! 6 
  param%p_inakm    = v_prm(7)  ! 7
  param%p_inacam   = v_prm(8)  ! 8
  param%p_gbca     = v_prm(9)  ! 9
  param%p_gbna     = v_prm(10) ! 10
  param%p_ipcam    = v_prm(11) ! 11
  param%p_krel     = v_prm(12) ! 12
  param%p_iupm     = v_prm(13) ! 13  
  param%p_gach     = v_prm(14) ! 14  
  param%p_ACh      = v_prm(15) ! 15 
  return
end subroutine put_param
!------------------------------------------------------------------------------!
function ic_Courte(ict) result(crn_ic)
!------------------------------------------------------------------------------!
! This function sets the initial conditions for the Model
! tt_m is the main structure which contains the membrane potential, ioninc 
! concentrations, and gate variables
!------------------------------------------------------------------------------!
  implicit none
!
! Original Courte
!
  real(rp), parameter :: mic      = 0.002908;
  real(rp), parameter :: hic      = 0.9649;
  real(rp), parameter :: jic      = 0.9775;
  real(rp), parameter :: oaic     = 0.03043;
  real(rp), parameter :: oiic     = 0.9992;
  real(rp), parameter :: uaic     = 0.004966;
  real(rp), parameter :: uiic     = 0.9986;
  real(rp), parameter :: xric     = 0.00003296;
  real(rp), parameter :: xsic     = 0.01869;
  real(rp), parameter :: dic      = 0.0001367;
  real(rp), parameter :: fic      = 0.9996;
  real(rp), parameter :: fcaic    = 0.7755;
  real(rp), parameter :: naiic    = 11.17;
  real(rp), parameter :: kiic     = 139.0;
  real(rp), parameter :: caiic    = 0.0001013;
  real(rp), parameter :: caupic   = 1.488;
  real(rp), parameter :: carelic  = 1.488;
  real(rp), parameter :: uic      = 0.0;
  real(rp), parameter :: veic     = 1.0;
  real(rp), parameter :: wic      = 0.9992;
  !.
  integer (ip), intent(in)  :: ict
  real(rp)                  :: crn_ic(20)
  
  select case (ict)
  case(0)
    crn_ic = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)   
  case( 16) !.----ATRIA
    crn_ic = (/ mic, hic, jic, oaic, oiic, uaic, uiic, xric, xsic, dic, fic,   &
		  fcaic, naiic, kiic, caiic, caupic, carelic, uic, veic, wic /)
  case default
    crn_ic = (/ mic, hic, jic, oaic, oiic, uaic, uiic, xric, xsic, dic, fic,   &
		  fcaic, naiic, kiic, caiic, caupic, carelic, uic, veic, wic /)
  end select
  !.
  return
end function ic_Courte
!------------------------------------------------------------------------------!
subroutine get_me_struct(str_me,v_me)
!------------------------------------------------------------------------------!
  implicit none
  type (t_crn), intent(in) :: str_me
  real(rp), intent(out)    :: v_me(20)
  v_me(1)  =   str_me%m 
  v_me(2)  =   str_me%h 
  v_me(3)  =   str_me%j
  v_me(4)  =   str_me%oa
  v_me(5)  =   str_me%oi
  v_me(6)  =   str_me%ua
  v_me(7)  =   str_me%ui
  v_me(8)  =   str_me%xr
  v_me(9)  =   str_me%xs
  v_me(10) =   str_me%d
  v_me(11) =   str_me%f
  v_me(12) =   str_me%fca
  v_me(13) =   str_me%nai
  v_me(14) =   str_me%ki
  v_me(15) =   str_me%cai
  v_me(16) =   str_me%caup
  v_me(17) =   str_me%carel
  v_me(18) =   str_me%u
  v_me(19) =   str_me%ve
  v_me(20) =   str_me%w
  return
end subroutine get_me_struct
!------------------------------------------------------------------------------!
subroutine put_me_struct(v_me,str_me)
!------------------------------------------------------------------------------!
  implicit none
  real(rp), intent(in)      :: v_me(:)
  type (t_crn), intent(out) :: str_me
  str_me%m        = v_me(1)
  str_me%h        = v_me(2)
  str_me%j        = v_me(3)
  str_me%oa       = v_me(4)
  str_me%oi       = v_me(5)
  str_me%ua       = v_me(6)
  str_me%ui       = v_me(7)
  str_me%xr       = v_me(8)
  str_me%xs       = v_me(9)
  str_me%d        = v_me(10)
  str_me%f        = v_me(11)
  str_me%fca      = v_me(12)
  str_me%nai      = v_me(13)
  str_me%ki       = v_me(14)
  str_me%cai      = v_me(15)
  str_me%caup     = v_me(16)
  str_me%carel    = v_me(17)
  str_me%u        = v_me(18)
  str_me%ve       = v_me(19)
  str_me%w        = v_me(20)
  return
end subroutine put_me_struct
!------------------------------------------------------------------------------! 
subroutine f_cur(v,nai,ki,cai,prm,ena,ek,INA,IK1,Ito,IKur,IKr,IKs,ICAL,INAK,  &
		   INACA,IBNA,IBCA,IPCA,IREL,ITR,IUP,IUPLEAK,IKACh,m,h,j,oa,oi, &
		   ua,ui,xr,xs,d,f,fca,u,ve,w,carel,caup)
!------------------------------------------------------------------------------!
  implicit none
  real(rp),intent(in)    :: v,nai,ki,cai
  type (t_prm)           :: prm
  real(rp),intent(out)   :: ena,ek,INA,IK1,Ito,IKur,IKr,IKs,ICAL,INAK,INACA,  &
				IBNA,IBCA,IPCA,IREL,ITR,IUP,IUPLEAK,IKACh
  real(rp),intent(inout) :: m,h,j,oa,oi,ua,ui,xr,xs,d,f,fca,u,ve,w,carel,caup
  real(rp)               :: FRT,eca,FVRT,alpa,fnak,gkur,nmr,dnm,bikach_clo,bik1_clo,Dclo

  FRT = (p_R*p_Te)/p_FA;
  ena = FRT*log(p_nao/nai);
  ek  = FRT*log(p_ko/ki);
  eca = 0.5*FRT*log(p_cao/cai);
  FVRT = (p_FA*v)/(p_R*p_Te);
  alpa = (1.0/7.0)*((exp(p_nao/67.3))-1.0);
  fnak = 1.0/(1.0+(0.1245*exp(-0.1*FVRT))+(0.0365*alpa*exp(-1.0*FVRT)));
  gkur = 0.005+(0.05/(1.0+exp((v-15.0)/-13.0)));

  Dclo = 2.0;  !---1-8.7 uM---!
  INA = 100.0*prm%p_gna*(m**3)*h*j*(v-ena);
  bik1_clo = 1.0/(1.0+(1.0/Dclo));
  IK1 =(1.0-bik1_clo)*(2.0*100.0*prm%p_gk1*(v-ek))/(1.0+exp(0.07*(v+80.0))); 
  Ito  = 0.5*100.0*prm%p_gto*(oa**3)*oi*(v-ek); 
  IKur = 0.5*100.0*gkur*(ua**3)*ui*(v-ek);! 
  IKr  = (100.0*prm%p_gkr*xr*(v-ek))/(1.0+exp((v+15.0)/22.4));
  IKs  = 100.0*prm%p_gks*(xs**2)*(v-ek);
  ICAL = 0.3*100.0*prm%p_gcal*d*f*fca*(v-65.0);
  INAK = 100.0*prm%p_inakm*fnak*(1.0/(1.0+((p_kmnai/nai)**1.5)))*(p_ko/(p_ko+p_kmko)); 
  nmr = (((nai**3)*p_cao)*exp(0.35*FVRT))-(((p_nao**3)*cai)*exp((0.35-1.0)*FVRT)); 
  dnm = ((p_kmna**3)+(p_nao**3))*(p_kmca+p_cao)*(1.0+(0.1*exp((0.35-1.0)*FVRT))); 
  INACA = 100.0*prm%p_inacam*(nmr/dnm);
  IBNA = 100.0*prm%p_gbna*(v-ena);
  IBCA = 100.0*prm%p_gbca*(v-eca);
  IPCA = 100.0*prm%p_ipcam*(cai/(0.0005+cai));
  IREL = prm%p_krel*(u**2)*ve*w*(carel-cai);            
  ITR = (caup-carel)/p_tautr;
  IUP = prm%p_iupm/(1.0+(p_kup/cai));
  IUPLEAK = prm%p_iupm*caup/p_caupm;
  bikach_clo = 1.0/(1.0+(0.97/Dclo));
  IKACh = (1.0-bikach_clo)*100.0*prm%p_gach*10.0/(1.0+9.13652/(prm%p_ACh**0.477811))*(0.0517+0.4516/(1.0+exp((v+59.53)/17.18)))*(v-ek);
   
  return
end subroutine f_cur
!------------------------------------------------------------------------------! 
subroutine f_comp(v,dt,INA,IK1,Ito,IKur,IKr,IKs,ICAL,IPCA,INAK,INACA,IBNA,    &
		    IBCA,IREL,ITR,IUP,IUPLEAK,IKACh,prm,cai,nai,ena,ek,ki,caup,  &
		    carel,m,h,j,u,ve,w,d,f,fca,xs,xr,ua,ui,oa,oi)
!------------------------------------------------------------------------------! 
  implicit none
  real(rp),intent(in)    :: v,dt,INA,IK1,Ito,IKur,IKr,IKs,ICAL,IPCA,INAK,    &
				INACA,IBNA,IBCA,IREL,ITR,IUP,IUPLEAK,IKACh
   type (t_prm)          :: prm
  real(rp),intent(inout) :: cai,nai,ena,ek,ki,caup,carel,m,h,j,u,ve,w,d,f,   &
 				fca,xs,xr,ua,ui,oa,oi
  real(rp)               :: am,bm,ah,bh,aj,bj,taum,tauh,tauj,mm,hm,jm,aoa,   &
				aoi,boa,boi,tauoa,tauoi,oam,oim,aua,bua,aui,bui, &
				uam,uim,tauua,tauui,axr,bxr,xrm,tauxr,axs,bxs,   &
				xsm,tauxs,taud,dm,tauf,fm,fcam,alpa,fnak,FVRT,   &
				fn,um,tauve,vem,tauw,wm,aux1,aux2,aux3,aux4,aux5

  fn = (1.0E-12)*p_volrel*IREL-((5.0E-13)/p_FA)*(0.5*ICAL-0.2*INACA);
  if (v == -47.13) then
      am = 3.2;
  else
      am = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
  endif
  bm = 0.08*exp(-v/11.0);
  taum = 1.0/(am+bm); 
  mm = am*taum;
  if (v >= -40.0) then
      ah = 0.0;
      bh = 1.0/(0.13*(1.0+exp((v+10.66)/-11.1)));
      aj = 0.0;
      bj = 0.3*(exp((-2.535E-7)*v))/(1.0+exp(-0.1*(v+32.0)));
  else
      ah = 0.135*exp((v+80.0)/-6.8);
      bh = (3.56*exp(0.079*v))+((3.1E5)*exp(0.35*v));
      aux4 = (v+37.78)/(1.0+exp(0.311*(v+79.23)));
      aj = (-127140.0*exp(0.2444*v)-(3.474E-5)*exp(-0.04391*v))*aux4;
      bj = 0.1212*(exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
  endif
  tauh = 1.0/(ah+bh); 
  hm = ah*tauh;
  tauj = 1.0/(aj+bj); 
  jm = aj*tauj;
  aoa = 0.65/(exp((v+10.0)/-8.5)+exp((v-30.0)/-59.0));
  boa = 0.65/(2.5+exp((v+82.0)/17.0));
  tauoa = (1.0/(aoa+boa))/3.0;
  oam = 1.0/(1.0+exp((v+20.47)/-17.54));
  aoi = 1.0/(18.53+exp((v+113.7)/10.95));
  boi = 1.0/(35.56+exp((v+1.26)/-7.44));
  tauoi = (1.0/(aoi+boi))/3.0;
  oim = 1.0/(1.0+exp((v+43.1)/5.3));
  aua = 0.65/(exp((v+10.0)/-8.5)+exp((v-30.0)/-59.0));
  bua = 0.65/(2.5+exp((v+82.0)/17.0));
  tauua = (1.0/(aua+bua))/3.0;
  uam = 1.0/(1.0+exp((v+30.3)/-9.6));
  aui = 1.0/(21.0+exp((v-185.0)/-28.0));
  bui = 1.0/(exp((v-158.0)/-16.0));
  tauui = (1.0/(aui+bui))/3.0;
  uim = 1.0/(1.0+exp((v-99.45)/27.48));
  axr = 0.0003*(v+14.1)/(1.0-exp((v+14.1)/-5.0));
  bxr = (7.3898E-5)*(v-3.3328)/((exp((v-3.3328)/5.1237))-1.0);
  tauxr = 1.0/(axr+bxr);
  xrm = 1.0/(1.0+exp((v+14.1)/-6.5));
  axs = (4.0E-5)*(v-19.9)/(1.0-exp((v-19.9)/-17.0));   
  bxs = (3.5E-5)*(v-19.9)/(exp((v-19.9)/9.0)-1.0);
  tauxs = 0.5/(axs+bxs);
  xsm = 1.0/((1.0+exp((v-19.9)/-12.7))**0.5);
  taud = (1.0-exp((v+10.0)/-6.24))/(0.035*(v+10.0)*(1.0+exp((v+10.0)/-6.24)));
  dm = 1.0/(1.0+exp((v+10.0)/-8.0));
  tauf = 9.0/(0.0197*exp(-1.0*0.0337*0.0337*((v+10.0)**2))+0.02);
  fm = 1.0/(1.0+exp((v+28.0)/6.9));
  fcam = 1.0/(1.0+(cai/0.00035));
  um = 1.0/(1.0+exp((fn-3.4175E-13)/-13.67E-16));
  tauve = 1.91+(2.09/(1.0+exp((fn-3.4175E-13)/-13.67E-16)));
  vem = 1.0-(1.0/(1.0+exp((fn-6.835E-14)/-13.67E-16)));
  tauw = 6.0*(1.0-exp((v-7.9)/-5.0))/(( 1.0+0.3*exp((v-7.9)/-5.0))*(v-7.9));
  wm = 1.0-(1.0/(1.0+exp((v-40.0)/-17.0)));

  nai = nai + (dt*((-3.0*INAK-3.0*INACA-IBNA-INA)/(p_FA*p_voli)));
  ki = ki + (dt*((2.0*INAK-IK1-Ito-IKur-IKr-IKs)/(p_FA*p_voli))); 
  aux1 = (2.0*INACA-IPCA-ICAL-IBCA)/(2.0*p_FA*p_voli);
  aux2 = ((p_volup*(IUPLEAK-IUP))+(IREL*p_volrel))/p_voli;
  aux5 = p_cmdnm*p_kmcmdn/((cai+p_kmcmdn)**2)
  aux3 = 1.0+(p_trpnm*p_kmtrpn/((cai+p_kmtrpn)**2))+(aux5);
  cai = cai + (dt*((aux1+aux2)/aux3)); 
  caup = caup + (dt*(IUP-IUPLEAK-ITR*p_volrel/p_volup));
  carel = carel + (dt*((ITR-IREL)/(1.0+(p_csqnm*p_kmcsqn)/((carel+p_kmcsqn)**2))));

  m = m + (dt*((mm-m)/taum)); 
  h = h + (dt*((hm-h)/tauh));
  j = j + (dt*((jm-j)/tauj));
  oa = oa + (dt*((oam-oa)/tauoa));
  oi = oi + (dt*((oim-oi)/tauoi));
  ua = ua + (dt*((uam-ua)/tauua));
  ui = ui + (dt*((uim-ui)/tauui));
  xr = xr + (dt*((xrm-xr)/tauxr));
  xs = xs + (dt*((xsm-xs)/tauxs));
  d = d + (dt*((dm-d)/taud)); 	
  f = f + (dt*((fm-f)/tauf));
  fca = fca + (dt*((fcam-fca)/p_taufca));
  u = u + (dt*((um-u)/p_tauu));
  ve = ve + (dt*((vem-ve)/tauve));
  w = w + (dt*((wm-w)/tauw));

  return
end subroutine f_comp
!------------------------------------------------------------------------------! 
subroutine Courte_A_P01 (dt, v, Iion, v_prm, v_crn, v_cr)
!------------------------------------------------------------------------------!
  implicit none
  real(rp),    intent(in)    :: dt, v
  real(rp),    intent(in)    :: v_prm(:)
  real(rp),    intent(out)   :: Iion
  real(rp),    intent(inout) :: v_crn(20)
  real(rp),    intent(out)   :: v_cr(:)
  type(t_crn)                :: courte
  type(t_prm)                :: param
  real(rp)                   :: ek,ena,INA,IK1,Ito,IKur,IKr,IKs,ICAL,IPCA,INAK, &         
  				    INACA,IBNA,IBCA,IREL,ITR,IUP,IUPLEAK,IKACh
  !
  call put_me_struct(v_crn,courte)
  call put_param(v_prm,param)
  !
  call f_cur(v,courte%nai,courte%ki,courte%cai,param,ena,ek,INA,IK1,Ito,IKur,  &
             IKr,IKs,ICAL,INAK,INACA,IBNA,IBCA,IPCA,IREL,ITR,IUP,IUPLEAK,IKACh,  &
             courte%m,courte%h,courte%j,courte%oa,courte%oi,courte%ua,courte%ui,  &
             courte%xr,courte%xs,courte%d,courte%f,courte%fca,courte%u,courte%ve, &
             courte%w,courte%carel,courte%caup)
  
  Iion = INA+IK1+Ito+IKur+IKr+IKs+ICAL+IPCA+INAK+INACA+IBNA+IBCA+IKACh;
  !
  call f_comp(v,dt,INA,IK1,Ito,IKur,IKr,IKs,ICAL,IPCA,INAK,INACA,IBNA,IBCA,IREL,  &
              ITR,IUP,IUPLEAK,IKACh,param,courte%cai,courte%nai,ena,ek,courte%ki,  &
              courte%caup,courte%carel,courte%m,courte%h,courte%j,courte%u,  &
              courte%ve,courte%w,courte%d,courte%f,courte%fca,courte%xs,courte%xr,  &
		courte%ua,courte%ui,courte%oa,courte%oi)
  
  call get_me_struct(courte,v_crn)

  !.
  v_cr(1:20) = (/ ena, INA, ek, IK1, Ito, IKur, IKr, IKs, ICAL, IPCA, INAK,  &
                  INACA, IBNA, IBCA, IREL, ITR, IUP, IUPLEAK, IKACh, Iion /)
  !.
  return
end subroutine Courte_A_P01
!-------------------------------------------------------------------------------
end module mod_courte
!-------------------------------------------------------------------------------




