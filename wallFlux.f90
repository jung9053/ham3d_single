!>
!
! flux at the wall
!
subroutine wallFlux(nxyz,ql,flx,specRadius,gamma)
!
implicit none
!
real*8, intent(in) :: nxyz(3)
real*8, intent(in) :: ql(5)
real*8, intent(inout) :: flx(5)
real*8, intent(inout) :: specRadius
real*8, intent(in) :: gamma
!
real*8 :: p,u,v,w,a,ri
!
p=(gamma-1)*(ql(5)-0.5*(ql(4)**2+ql(3)**2+ql(2)**2)/ql(1))
ri=1./ql(1)
u=ql(2)*ri
v=ql(3)*ri
w=ql(4)*ri
a=sqrt(gamma*p*ri);
!
specRadius=(sqrt(u*u+v*v+w*w)+a)*sqrt(nxyz(1)**2+nxyz(2)**2+nxyz(3)**2)
!
flx(1)=0
flx(2)=p*nxyz(1)
flx(3)=p*nxyz(2)
flx(4)=p*nxyz(3)
flx(5)=0
!
return
end subroutine wallFlux
!>
!
! Jacobian of wall flux
!
subroutine wallFluxJacobian(nxyz,ql,lmat,rmat,gamma)
!
implicit none
!
real*8, intent(in) :: nxyz(3)
real*8, intent(in) :: ql(5)
real*8, intent(inout) :: lmat(5,5)
real*8, intent(inout) :: rmat(5,5)
real*8, intent(in) :: gamma
!
real*8 :: u,v,w,umag2,gm1
!
u=ql(2)/ql(1)
v=ql(3)/ql(1)
w=ql(4)/ql(1)

umag2=u*u+v*v+w*w
gm1=gamma-1
!
lmat(1,:)=0.
lmat(5,:)=0.
!
! ys
!lmat(2,1) = gm1*(ql(4)/ql(1)-umag2*0.5)*nxyz(1)
lmat(2,1)=gm1*0.5*umag2*nxyz(1)
lmat(2,2)=-gm1*u*nxyz(1)
lmat(2,3)=-gm1*v*nxyz(1)
lmat(2,4)=-gm1*w*nxyz(1)
lmat(2,5)=gm1*nxyz(1)
!
!lmat(3,1) = gm1*(ql(4)/ql(1)-umag2*0.5)*nxyz(2)
lmat(3,1)=gm1*0.5*umag2*nxyz(2)
lmat(3,2)=-gm1*u*nxyz(2)
lmat(3,3)=-gm1*v*nxyz(2)
lmat(3,4)=-gm1*w*nxyz(2)
lmat(3,5)=gm1*nxyz(2)
!
lmat(4,1)=gm1*0.5*umag2*nxyz(3)
lmat(4,2)=-gm1*u*nxyz(3)
lmat(4,3)=-gm1*v*nxyz(3)
lmat(4,4)=-gm1*w*nxyz(3)
lmat(4,5)=gm1*nxyz(3)
!
rmat=0.
!
return
end subroutine wallFluxJacobian
