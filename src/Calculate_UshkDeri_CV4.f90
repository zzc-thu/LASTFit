SUBROUTINE Calculate_UshkDeri_CV4_B(u_inf, v_inf, w_inf, rho_inf, &
& shockv, shockvb, mach_ref, gamma, shockh, shockhb, shockhdxi, &
& shockhdxib, shockhdzeta, shockhdzetab, heta, hetadxi, hetadeta, &
& hetadzeta, wallsxdxi, wallsydxi, wallszdxi, wallsxdzeta, wallsydzeta, &
& wallszdzeta, wallnormalx, wallnormaly, wallnormalz, wallnormalxdxi, &
& wallnormalydxi, wallnormalzdxi, wallnormalxdzeta, wallnormalydzeta, &
& wallnormalzdzeta, cv_4b)
  IMPLICIT NONE
! Input variables ShockH, ShockHdxi, ShockHdzeta
! Output variables: CV_1,CV_2,CV_3,CV_4,CV_5
  REAL(kind=8) :: delta_inf, masn2inf, deltashkuinf
  REAL(kind=8) :: delta_infb, masn2infb, deltashkuinfb
  REAL(kind=8) :: u_inf, v_inf, w_inf, rho_inf
  REAL(kind=8) :: shocknormalx, shocknormaly, shocknormalz
  REAL(kind=8) :: shocknormalxb, shocknormalyb, shocknormalzb
  REAL(kind=8) :: wallnormalx, wallnormaly, wallnormalz
  REAL(kind=8) :: w, rho, mach_ref, gamma
  REAL(kind=8) :: wb, rhob
  REAL(kind=8) :: etax, etay, etaz, nablaeta
  REAL(kind=8) :: etaxb, etayb, etazb, nablaetab
  REAL(kind=8) :: xdxi, ydxi, zdxi, xdzeta, ydzeta, zdzeta, jaco
  REAL(kind=8) :: xdxib, ydxib, zdxib, xdzetab, ydzetab, zdzetab, jacob
  REAL(kind=8) :: xdeta, ydeta, zdeta
  REAL(kind=8) :: xdetab, ydetab, zdetab
  REAL(kind=8) :: heta, hetadxi, hetadeta, hetadzeta
  REAL(kind=8) :: wallsxdxi, wallsydxi, wallszdxi
  REAL(kind=8) :: wallsxdzeta, wallsydzeta, wallszdzeta
  REAL(kind=8) :: wallnormalxdxi, wallnormalydxi, wallnormalzdxi
  REAL(kind=8) :: wallnormalxdzeta, wallnormalydzeta, wallnormalzdzeta
  REAL(kind=8) :: temp
  REAL(kind=8) :: tempb
  REAL(kind=8), INTENT(IN) :: shockh, shockhdxi, shockhdzeta, shockv
  REAL(kind=8) :: shockhb, shockhdxib, shockhdzetab, shockvb
  !REAL(kind=8) :: cv_1, cv_2, cv_3, cv_4, cv_5
  REAL(kind=8) :: cv_4b
  INTRINSIC SQRT
  REAL(kind=8) :: tempb0
  REAL(kind=8) :: tempb1
  REAL(kind=8) :: tempb2
  xdxi = wallsxdxi + (wallnormalxdxi*shockh+wallnormalx*shockhdxi)*heta &
&   + wallnormalx*shockh*hetadxi
  ydxi = wallsydxi + (wallnormalydxi*shockh+wallnormaly*shockhdxi)*heta &
&   + wallnormaly*shockh*hetadxi
  zdxi = wallszdxi + (wallnormalzdxi*shockh+wallnormalz*shockhdxi)*heta &
&   + wallnormalz*shockh*hetadxi
  xdeta = wallnormalx*shockh*hetadeta
  ydeta = wallnormaly*shockh*hetadeta
  zdeta = wallnormalz*shockh*hetadeta
  xdzeta = wallsxdzeta + (wallnormalxdzeta*shockh+wallnormalx*&
&   shockhdzeta)*heta + wallnormalx*shockh*hetadzeta
  ydzeta = wallsydzeta + (wallnormalydzeta*shockh+wallnormaly*&
&   shockhdzeta)*heta + wallnormaly*shockh*hetadzeta
  zdzeta = wallszdzeta + (wallnormalzdzeta*shockh+wallnormalz*&
&   shockhdzeta)*heta + wallnormalz*shockh*hetadzeta
  temp = xdxi*(ydeta*zdzeta-zdeta*ydzeta) + xdeta*(ydzeta*zdxi-ydxi*&
&   zdzeta) + xdzeta*(ydxi*zdeta-ydeta*zdxi)
  jaco = 1.d0/temp
  etax = -(jaco*(ydxi*zdzeta-zdxi*ydzeta))
  etay = jaco*(xdxi*zdzeta-zdxi*xdzeta)
  etaz = -(jaco*(xdxi*ydzeta-ydxi*xdzeta))
  nablaeta = SQRT(etax*etax + etay*etay + etaz*etaz)
! Given ShockNormalX, ShockNormalY, ShockNormalZ 
  shocknormalx = etax/nablaeta
  shocknormaly = etay/nablaeta
  shocknormalz = etaz/nablaeta
! Velocity ahead of the shocks
  delta_inf = u_inf*shocknormalx + v_inf*shocknormaly + w_inf*&
&   shocknormalz - shockv*(wallnormalx*shocknormalx+wallnormaly*&
&   shocknormaly+wallnormalz*shocknormalz)
! Shock normal Mach number^2
  masn2inf = (mach_ref*delta_inf)**2.d0
! Velocity normal to the shock after the shock
  deltashkuinf = 2.d0*delta_inf/(gamma+1.d0)*(1.d0/masn2inf-1.d0)
  w = w_inf + deltashkuinf*shocknormalz
  rho = rho_inf*(gamma+1.d0)*masn2inf/((gamma-1.d0)*masn2inf+2.d0)
  rhob = w*cv_4b
  wb = rho*cv_4b
  tempb2 = rho_inf*(gamma+1.d0)*rhob/((gamma-1.d0)*masn2inf+2.d0)
  deltashkuinfb = shocknormalz*wb
  masn2infb = (1.0-(gamma-1.d0)*masn2inf/((gamma-1.d0)*masn2inf+2.d0))*&
&   tempb2 - delta_inf*2.d0*deltashkuinfb/(masn2inf**2*(gamma+1.d0))
  delta_infb = (1.0/masn2inf-1.d0)*2.d0*deltashkuinfb/(gamma+1.d0) + &
&   mach_ref**2*2.d0*delta_inf*masn2infb
  shockvb = -((wallnormalx*shocknormalx+wallnormaly*shocknormaly+&
&   wallnormalz*shocknormalz)*delta_infb)
  tempb2 = -(shockv*delta_infb)
  shocknormalzb = deltashkuinf*wb + w_inf*delta_infb + wallnormalz*&
&   tempb2
  shocknormalxb = u_inf*delta_infb + wallnormalx*tempb2
  shocknormalyb = v_inf*delta_infb + wallnormaly*tempb2
  nablaetab = -(etaz*shocknormalzb/nablaeta**2) - etay*shocknormalyb/&
&   nablaeta**2 - etax*shocknormalxb/nablaeta**2
  !IF (etax**2 + etay**2 + etaz**2 .EQ. 0.0) THEN
  !  tempb2 = 0.0_8
  !ELSE
    tempb2 = nablaetab/(2.0*SQRT(etax**2+etay**2+etaz**2))
  !END IF
  etazb = shocknormalzb/nablaeta + 2*etaz*tempb2
  etayb = shocknormalyb/nablaeta + 2*etay*tempb2
  etaxb = shocknormalxb/nablaeta + 2*etax*tempb2
  jacob = (xdxi*zdzeta-zdxi*xdzeta)*etayb - (xdxi*ydzeta-ydxi*xdzeta)*&
&   etazb - (ydxi*zdzeta-zdxi*ydzeta)*etaxb
  tempb2 = -(jaco*etazb)
  xdxib = ydzeta*tempb2
  ydzetab = xdxi*tempb2
  ydxib = -(xdzeta*tempb2)
  xdzetab = -(ydxi*tempb2)
  tempb2 = jaco*etayb
  zdzetab = xdxi*tempb2
  zdxib = -(xdzeta*tempb2)
  tempb = -(jacob/temp**2)
  xdxib = xdxib + zdzeta*tempb2 + (ydeta*zdzeta-zdeta*ydzeta)*tempb
  xdzetab = xdzetab + (ydxi*zdeta-ydeta*zdxi)*tempb - zdxi*tempb2
  tempb2 = -(jaco*etaxb)
  ydxib = ydxib + zdzeta*tempb2
  zdxib = zdxib - ydzeta*tempb2
  tempb0 = xdxi*tempb
  xdetab = (ydzeta*zdxi-ydxi*zdzeta)*tempb
  tempb1 = xdeta*tempb
  zdzetab = zdzetab + ydxi*tempb2 + ydeta*tempb0 - ydxi*tempb1
  ydzetab = ydzetab + zdxi*tempb1 - zdxi*tempb2 - zdeta*tempb0
  tempb2 = xdzeta*tempb
  ydxib = ydxib + zdeta*tempb2 - zdzeta*tempb1
  zdetab = ydxi*tempb2 - ydzeta*tempb0
  ydetab = zdzeta*tempb0 - zdxi*tempb2
  zdxib = zdxib + ydzeta*tempb1 - ydeta*tempb2
  shockhb = (wallnormalzdzeta*heta+wallnormalz*hetadzeta)*zdzetab + (&
&   wallnormalydzeta*heta+wallnormaly*hetadzeta)*ydzetab + (&
&   wallnormalxdzeta*heta+wallnormalx*hetadzeta)*xdzetab + wallnormalz*&
&   hetadeta*zdetab + wallnormaly*hetadeta*ydetab + wallnormalx*hetadeta&
&   *xdetab + (wallnormalzdxi*heta+wallnormalz*hetadxi)*zdxib + (&
&   wallnormalydxi*heta+wallnormaly*hetadxi)*ydxib + (wallnormalxdxi*&
&   heta+wallnormalx*hetadxi)*xdxib
  shockhdzetab = wallnormalz*heta*zdzetab + wallnormaly*heta*ydzetab + &
&   wallnormalx*heta*xdzetab
  shockhdxib = wallnormalz*heta*zdxib + wallnormaly*heta*ydxib + &
&   wallnormalx*heta*xdxib
  cv_4b = 0.0_8
END SUBROUTINE Calculate_UshkDeri_CV4_B