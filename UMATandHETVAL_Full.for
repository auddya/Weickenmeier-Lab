
! ------------------------------------------------------------------------------------

! ------------------------------------------------------------------------------------
! ------------------------------------------------------------------------------------



!
!     # Changes made by Andreia Cacoilo @Stevens Institute of Technology
!
! **********************************************************************
!
! Copyright 2019 Stephen Connolly
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
! BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!
! Programs tested using Abaqus/Standard v2016 HF6, Microsoft Visual Studio 2010 
!      and Intel Parallel Studio XE 2013.
!
!    UMAT for Ogden model to compute the Cauchy stress (stress), 
!     spatial elasticity modulus (DDSDDE) (in terms of the Jaumann-rate
!     of the Cauchy stress), and strain energy density (SSE) for a 
!     hyperelastic constitutive model in terms of isochoric principal
!     stretches using explicit computation of the eigenvectors.
!
!    User input required in "la_sub" and "kstress"
!
!    Stephen Connolly, December 2019.
!      Tested and implemented in Abaqus 2016 HF6
!      Cannot be used for plane stress
!
! **********************************************************************
!      
!     OPTIONAL ADDITIONAL VARIABLES TO BE DEFINED BY USER
!
! **********************************************************************
!     #ac _ ddsdde = Jacobian matrix of the constitutive model
!     #ac _ stress = stress matrix 6 component 
!     #ac _ sse = specific elastic strain energy 
!     #ac _ sse = optional solution to state variable
!     #ac _ spd = specific plastic dissipation
!     #ac _ ntens = size of the stress/strain component tensor

      subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl, 
     1  ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, 
     2  dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, 
     3  nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel, 
     4  npt, layer, kspt, kstep, kinc)
!
!
      include 'aba_param.inc'
!
!      
      character*(*) cmname

!     #ac _ real(8) statement specifies the variable names to be double
!     #ac _ precision 8-byte real numbers which has 15 digits of accuracy
      real(8) stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3, 3),
     3 dfgrd0(3, 3), dfgrd1(3, 3), time(2)
!     
!     Define arrays for axisymmetric and 3D dummy tangents
      real(8) axi(6,6), k3d(4,4)

      parameter(zero=0d0, one=1d0, two=2d0, three=3d0, four=4d0,
     1          nine=9d0, half=5d-1, third=1d0/3d0, ninth=1d0/9d0)


!      write(*,*) 'Now, I am in UMAT.'


!     #ac _ nshr: number of engineering shear stress components
!      write(*,*) 'I am here.'
      if (cmname.eq.'ZEROSTIFFNESS') then
            call zero_stiffness(dfgrd1,nshr,ntens,stress,ddsdde,sse,statev,
     1                  nstatv,kinc)
      else
            if (nshr.eq.1) then ! ###_ ac _ eq means =
                  call la_sub(dfgrd1,nshr,props,stress,axi,DDSDDE,SSE,
     1                  statev,nstatv,temp,time,dtime,rpl,drplde,drpldt)
            else if (nshr.eq.3) then
                  call la_sub(dfgrd1,nshr,props,stress,DDSDDE,k3d,SSE,
     1                  statev,nstatv,temp,time,dtime,rpl,drplde,drpldt)
            end if
      end if


      return
      end subroutine umat
!
!***********************************************************************
      subroutine zero_stiffness(kdefg,nsh,n,sigma,tangent,energy,stv,
     1      nstv,inc)
!
!    Inputs: 
!     deformation gradient:                     kdefg(3,3)
!     # direct + # shear stress components:     n(integer)
!     # strain components:                      nsh(integer)
!
!    Outputs:
!     Cauchy stress tensor:   sigma(n)
!     Material's tangent:     tangent(n,n)
!     Strain energy density:  energy(real)
!
      implicit none

      real(8), intent(in) :: kdefg(3,3)
      integer, intent(in) :: n, nsh, nstv, inc

      real(8), intent(out) :: sigma(n)
      real(8), intent(out) :: tangent(n,n)
      real(8), intent(out) :: energy

      real(8), intent(inout) :: stv(nstv)

      integer :: k, l
      real(8) :: J, Aimag, afmag, astr
      real(8), dimension(3) :: Ai, af
      real(8), dimension(3,3) :: F, FT, FTinv


!      write(*,*) 'Now, I am in zero_stiffness.'


!    Stress
      do k=1, n
         sigma(k)=0.0d0
      enddo
!    Tangent
      do k=1, n
         do l=1, n
            tangent(k,l)=0.0d0
         enddo
      enddo
!    Energy
      energy=0.0d0

!    Deformation gradient
      F = kdefg
!      write(*,*) 'F:', F
      call bdet(F, nsh, J)
!      write(*,*) 'det:', J
      FT = transpose(F)
      call invertmat33(FT,J,FTinv)


!    Area stretch
!    Initial area vector and magnitude
!      write (*,*) 'stv 1:', stv(1)
      Ai(1) = stv(1)
      Ai(2) = stv(2)
      Ai(3) = stv(3)
      Aimag = stv(4)
!      write(*,*) 'Initial area vector:', Ai
!      write(*,*) 'Initial area:', Aimag
!    Deformed area vector and magnitude
      af = J*matmul(FTinv,Ai)
      afmag = dsqrt(af(1)*af(1) + af(2)*af(2) + af(3)*af(3))
!      write(*,*) 'Deformed area vector:', af
!      write(*,*) 'Deformed area:', afmag
!    Stretch  
      astr = afmag/Aimag
!      write(*,*) 'Area stretch:', astr


!    Saving state variables 
!    Defined in UMAT:
      stv(1)=stv(1)
      stv(2)=stv(2)
      stv(3)=stv(3)
      stv(4)=stv(4)
      stv(5)=astr


!    Triggering diffusion
      stv(6)=0
      if (astr.gt.1.4d0) then
            stv(6)=1
      end if

      
!    Saving state variables 
!    Defined in HETVAL:
      stv(7)=stv(7)



      return
      end subroutine zero_stiffness

!*********************************************************************
      subroutine la_sub(kdefg,knshr,kprops,CST,kEtens,axtens,kSSE,
     1                  sttv,nstv,tmptr,utime,udtime,htsrc,htsrcde,
     2                  htsrcdt)
!
!    This program computes the stress (CST), tangent modulus (kEtens)
!     and strain energy density (SSE) for a hyperelastic material model
!     defined in terms of principal stretches; using analytically 
!     derived expressions.
!
!    Inputs: 
!     deformation gradient:   kdefg(3,3)
!     number of shear terms:  knshr (integer)
!     material properties:    kPROPS(N)
!
!    Outputs:
!     Cauchy stress tensor:   CST(6)
!     3D Elasticity Tensor:   kEtens(6,6)
!     Axi Elasticity Tensor:  axtens(4,4)
!     Specific Elastic SE:    kSSE (real)
!
      implicit none
      real(8) J,tol1,kc,kd,PI,W,U,kSSE,utime(2),htsrc,htsrcde(6),
     1   mxGLx, mxGLy, mxGLz, maxpGL, minpGL, maxpCST, minpCST,udtime,
     2   kprops(12),mu1,mu2,mu3,mu4,kl1,kl2,kl3,kl4,kappa,gc,gnot,
     3   nn(3),nn0(3),lamNN,lamTTmax,lamTT2,lamTTmin,lamNNCST,
     4   tt(3),tt0(3), tt1(3),tt2(3),temp,tmptr,theta,nrm,theg,thegn,
     5   lamTTmax_vec(3),lamTT2_vec(3),lamTTmin_vec(3),htsrcdt,
     6   ratiocrit,lnJe,detfe,detfg,lamTTCSTmax_vec(3),lamTTCST2_vec(3),
     7   lamTTCSTmin,lamTTCSTmin_vec(3),lamTTCSTmax,lamTTCST2,lamTT
      real(8), dimension(3) :: la2,la2bar,la,labar,dWdLa,B_a,
     1   laE, laE2, laCST
      real(8), dimension(3,3) :: kdefg,F,FT,b,na,Ga,Id,d2WLa2, C,
     1   I2A, E, nE, nCST, CSTM,Finv,FinvT,fe, fg, fginv,binv
      real(8), dimension(6) :: bv,n11,n22,n33,Iden,kiso,VSC,CST,n12,
     1   n13,n23,htsrcdes
      real(8), dimension(6,6) :: Idy,Id4,a_vol,n11dy11,n22dy11,n33dy11,
     1   n11dy22,n22dy22,n33dy22,n11dy33,n22dy33,n33dy33,n12dy12,
     2   n13dy13,n23dy23,a_iso,CSTdyId,IddyCST,KETENS
      real(8), dimension(4,4) :: axtens
      integer knshr, a, k1, k2, i, jj, k, l, locmin, locmax,locn, loct,
     1   nstv
      integer xim (6), xjm(6), CSTT(3,3)
      real(8) sttv(nstv)
!      integer, dimension(3) :: locmax(3)

      
!
!    Numerical Parameters
!
!     #ac _ parameter(...) ! defined parameters to use in the formulas below
!
      real(8) zero, one, two, three, four, nine, half, third, ninth
      parameter(zero=0d0, one=1d0, two=2d0, three=3d0, four=4d0,
     1          nine=9d0, half=5d-1, third=1d0/3d0, ninth=1d0/9d0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Define the number of material parameters and their symbolic notation

      data xim/1,2,3,1,1,2/ ! #ac _ xim and xjm are used to convert normal
      data xjm/1,2,3,2,3,3/ ! #ac _ matrix notation to voigt notation
    
!
!    Assign material coefficients
!     #ac _ kprops: array with material property data
!     #ac _ respect this order when you enter material properties in CAE
!
      mu1 = kprops(2)*two
      kl1 = kprops(3) 
      mu2 = kprops(4)*two
      kl2 = kprops(5)
      mu3 = kprops(6)*two
      kl3 = kprops(7)
      mu4 = kprops(8)*two
      kl4 = kprops(9)
      kappa=kprops(10)  ! kappa=E/(3*(1-2v))
      gnot = kprops(11)
      ratiocrit=kprops(12)!kprops(12)

!      write(*,*) 'mu1:', mu1, 'kl1:', kl1, 'mu2:', mu2, 'kl2:', kl2
!      write(*,*) 'mu3:', mu3, 'kl3:', kl3, 'mu4:', mu4, 'kl4:', kl4
!      write(*,*) 'kappa:', kappa, 'gnot:', gnot, 'ratiocrit:', ratiocrit
!      STOP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      write(*,*) 'And now, I am in la_sub.'
!    Define tolerance for L’Hôpital’s rule
      tol1 = 1d-6
!
!    Calculate the determinant of deformation gradient: J(real)
      F = kdefg
!      write(*,*) 'F:', F
      call bdet(F, knshr, J)
!      write(*,*) 'det:', J
!      STOP
!
c...  update atrophy 
      thegn      = sttv(1) + gnot*udtime
      theg       = 1.0d0-thegn
!      write(*,*) 'atrophy:', thegn
!      write(*,*) 'scaling factor:', theg
!      STOP

      do k1=1,3
        do k2=1,3
          fg(k1,k2)=0.0d0
        enddo
      enddo
      fg(1,1) = theg
      fg(2,2) = theg
      fg(3,3) = theg
!
      detfg   = theg*theg*theg
!
      call invertmat33(fg,detfg,fginv)
      fe      = matmul(F,fginv) ! #ac _ Fe(elastic tensor)=F*Fg^{-1}
!
      detfe   = J/detfg
!
      sttv(1) = thegn
!
     

!     #ac _  matI2 is a subroutine which is in utils.f file
      call matI2(I2A, 3) 
!
!    Calculate left Cauchy-Green deformation tensor: b(3,3)=F*(F^T)
      F  = fe
!      write(*,*) 'Fe:', F
      FT = transpose(fe)
      b  = matmul(fe,FT) ! matmul: performs a matrix multiplication
      J  = detfe
!
!    Store b in Voigt notation, bv(6)
!
      bv(1) = b(1,1)
      bv(2) = b(2,2)
      bv(3) = b(3,3)
      bv(4) = b(1,2)
      bv(5) = b(1,3)
      bv(6) = b(2,3)
!
!    Calculate eigenvalues la2(a) and eigenvectors n_a(3,3) of b
!     (Note: Jacobian algorithm destroys upper triangular components 12,13,23)
!
      call DSYEVJ3(b,na,la2)
!
!    Calculate the eigenvalues' square-root: la(3)
!
      la(1) = (la2(1)**half)
      la(2) = (la2(2)**half)
      la(3) = (la2(3)**half)
!
!    Calculate the isochoric eigenvalues: la2bar(3) 
!
      la2bar(1) = la2(1)*(J**(-two/three))
      la2bar(2) = la2(2)*(J**(-two/three))
      la2bar(3) = la2(3)*(J**(-two/three))
!
!    Calculate their square-root: labar(3)
!    
!
      labar(1) = (la2bar(1)**half)
      labar(2) = (la2bar(2)**half)
      labar(3) = (la2bar(3)**half)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!  USER INPUT  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Define uncoupled strain energy density function, its derivatives
!     and any additional variables
!
!    Note: additional variables must be declared if "implicit none" is
!           is used
!
!     Define the isochoric energy, W
!
      W= (((mu1)/(kl1**two))* 
     1   ((labar(1)**(kl1))+(labar(2)**(kl1))+(labar(3)**(kl1))-three))
      if (kprops(1).ge.2) then !#ac _ ge means >=
      W= W + (((mu2)/(kl2**two))* 
     1   ((labar(1)**(kl2))+(labar(2)**(kl2))+(labar(3)**(kl2))-three))
      end if
      if (kprops(1).ge.3) then
      W= W + (((mu3)/(kl3**two))* 
     1   ((labar(1)**(kl3))+(labar(2)**(kl3))+(labar(3)**(kl3))-three))
      end if
      if (kprops(1).ge.4) then
      W= W + (((mu4)/(kl4**two))* 
     1   ((labar(1)**(kl4))+(labar(2)**(kl4))+(labar(3)**(kl4))-three))
      end if
!
!     Define the Volumetric energy, U
!
      U = kappa/four*(J**two-one-two*dlog(J)) !#ac _ dlog = calculates the logarithm
!
!     Combine the additive energy contributions to calculate the
!      total strain energy density
!
      kSSE = W + U
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch labar(k1)
!
      do k1=1,3
      dWdLa(k1) = ((mu1/kl1)*labar(k1)**(kl1-one))
      if (kprops(1).ge.2) then
      dWdLa(k1) = dWdLa(k1) + ((mu2/kl2)*labar(k1)**(kl2-one))
      end if
      if (kprops(1).ge.3) then
      dWdLa(k1) = dWdLa(k1) + ((mu3/kl3)*labar(k1)**(kl3-one))
      end if
      if (kprops(1).ge.4) then
      dWdLa(k1) = dWdLa(k1) + ((mu4/kl4)*labar(k1)**(kl4-one))
      end if
      end do
!
!    Calculate the first derivative of the constitutive model with
!     respect to the isochoric principal stretch labar(k1)
!
      do k1=1,3
       do k2=1,3
        if (k1.eq.k2) then
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a = b
!
         d2WLa2(k1,k2) = ((mu1/kl1)*(kl1-one)*(labar(k1)**(kl1-two)))
         if (kprops(1).ge.2) then
         d2WLa2(k1,k2) = d2WLa2(k1,k2) + ((mu2/kl2)*(kl2-one)*
     1        (labar(k1)**(kl2-two)))
         end if
         if (kprops(1).ge.3) then
         d2WLa2(k1,k2) = d2WLa2(k1,k2) + ((mu3/kl3)*(kl3-one)*
     1        (labar(k1)**(kl3-two)))
         end if
         if (kprops(1).ge.4) then
         d2WLa2(k1,k2) = d2WLa2(k1,k2) + ((mu4/kl4)*(kl4-one)*
     1        (labar(k1)**(kl4-two)))
         end if
!
        else
!
!    2nd derivative of W with respect to Labar(a) and Labar(b), a /= b
!     (a not equal b)            
!
         d2WLa2(k1,k2) = zero
!
        end if
       end do
      end do
!
!    Calculate volumetric stress coefficient, kc=dU/dJ
!
      kc = kappa/two*(one*J-one/J)
!
!    Calculate volumetric elasticity coefficient, kd=d^2U/dJ^2
!
      kd = kappa/two*(one+one/(J**two))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! END OF USER INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!***********************************************************************
!                   CALCULATE THE CAUCHY STRESS
!***********************************************************************
!
!    Calculate the stress coefficients B_a(3)
!
      do k1=1,3
          B_a(k1) = dWdLa(k1)*labar(k1) - third*(dWdLa(1)*labar(1) +
     1         dWdLa(2)*labar(2) + dWdLa(3)*labar(3))
      end do
!
!    Calculate the eigenvector dyadic products (eigenvalue bases)
!
!    n11
      n11(1) = na(1,1)*na(1,1)
      n11(2) = na(2,1)*na(2,1)
      n11(3) = na(3,1)*na(3,1)
      n11(4) = na(1,1)*na(2,1)
      n11(5) = na(1,1)*na(3,1)
      n11(6) = na(2,1)*na(3,1)
!    n22
      n22(1) = na(1,2)*na(1,2)
      n22(2) = na(2,2)*na(2,2)
      n22(3) = na(3,2)*na(3,2)
      n22(4) = na(1,2)*na(2,2)
      n22(5) = na(1,2)*na(3,2)
      n22(6) = na(2,2)*na(3,2)
!    n33
      n33(1) = na(1,3)*na(1,3)
      n33(2) = na(2,3)*na(2,3)
      n33(3) = na(3,3)*na(3,3)
      n33(4) = na(1,3)*na(2,3)
      n33(5) = na(1,3)*na(3,3)
      n33(6) = na(2,3)*na(3,3)
!
!    Define the second-order identitiy tensor, Id(3,3)
!
      do k1=1,3
       do k2=1,3
        if (k1.eq.k2) then
         Id(k1,k2) = 1d0
        else
         Id(k1,k2) = 0d0
        end if
       end do
      end do
!
!    Define the second-order identitiy tensor in Voigt notation, Iden(6)
!
      do k1=1,6
          if (k1.lt.4) then ! #ac _ lt means <
              Iden(k1) = 1d0
          else
              Iden(k1) = 0d0
          end if
      end do
!
!    Calculate the isochoric Kirchhoff stress, kiso(6)
!
      do k1=1,6
          kiso(k1) = (B_a(1)*n11(k1))+(B_a(2)*n22(k1))+(B_a(3)*n33(k1))
      end do
!
!    Calculate the volumetric Kirchhoff Stress, VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the volumetric Kirchhoff Stress, VSC(6)
!
      do k1=1,6
          VSC(k1) = kc*J*Iden(k1)
      end do
!
!    Calculate the Cauchy Stress Tensor, CST(6)
!
      do k1=1,6
       CST(k1) = (one/J)*(kiso(k1) + VSC(k1))
      end do
!
!
!***********************************************************************
!                   CALCULATE THE ELASTICITY MODULI
!***********************************************************************
!
!    Calculate the isochoric elasticity coefficient, Ga(3,3)
!
      do k1=1,3
       do k2=1,3
        Ga(k1,k2) = ((d2WLa2(k1,k2)*labar(k1)*labar(k2))+(dWdLa(k1) 
     &    *Id(k1,k2)*labar(k1))) + 
     &    (ninth*( 
     &    ((d2WLa2(1,1)*la2bar(1))+(dWdLa(1)*labar(1))) + 
     &    ((d2WLa2(2,2)*la2bar(2))+(dWdLa(2)*labar(2))) + 
     &    ((d2WLa2(3,3)*la2bar(3))+(dWdLa(3)*labar(3))) + 
     &    (two*( 
     &    ((d2WLa2(1,2))*labar(1)*labar(2)) + 
     &    ((d2WLa2(1,3))*labar(1)*labar(3)) + 
     &    ((d2WLa2(2,3))*labar(2)*labar(3)))))) - 
     &    (third*( 
     &    ((d2WLa2(k1,1)*labar(k1)*labar(1))+(dWdLa(1) 
     &    *Id(k1,1)*labar(1))) + 
     &    ((d2WLa2(k1,2)*labar(k1)*labar(2))+(dWdLa(2) 
     &    *Id(k1,2)*labar(2))) + 
     &    ((d2WLa2(k1,3)*labar(k1)*labar(3))+(dWdLa(3) 
     &    *Id(k1,3)*labar(3))) + 
     &    ((d2WLa2(1,k2)*labar(1)*labar(k2))+(dWdLa(k2) 
     &    *Id(1,k2)*labar(k2))) + 
     &    ((d2WLa2(2,k2)*labar(2)*labar(k2))+(dWdLa(k2) 
     &    *Id(2,k2)*labar(k2))) + 
     &    ((d2WLa2(3,k2)*labar(3)*labar(k2))+(dWdLa(k2) 
     &    *Id(3,k2)*labar(k2)))))
       end do
      end do
!
!    Calculate the "average" eigenvector dyads
!
!    na12:
      n12(1) = (na(1,1)*na(1,2))
      n12(2) = (na(2,1)*na(2,2))
      n12(3) = (na(3,1)*na(3,2))
      n12(4) = 0.5d0*((na(1,1)*na(2,2)) + (na(1,2)*na(2,1)))
      n12(5) = 0.5d0*((na(1,1)*na(3,2)) + (na(1,2)*na(3,1)))
      n12(6) = 0.5d0*((na(2,1)*na(3,2)) + (na(2,2)*na(3,1)))
!    na13:
      n13(1) = (na(1,1)*na(1,3))
      n13(2) = (na(2,1)*na(2,3))
      n13(3) = (na(3,1)*na(3,3))
      n13(4) = 0.5d0*((na(1,1)*na(2,3)) + (na(1,3)*na(2,1)))
      n13(5) = 0.5d0*((na(1,1)*na(3,3)) + (na(1,3)*na(3,1)))
      n13(6) = 0.5d0*((na(2,1)*na(3,3)) + (na(2,3)*na(3,1)))
!    na23:
      n23(1) = (na(1,2)*na(1,3))
      n23(2) = (na(2,2)*na(2,3))
      n23(3) = (na(3,2)*na(3,3))
      n23(4) = 0.5d0*((na(1,2)*na(2,3)) + (na(1,3)*na(2,2)))
      n23(5) = 0.5d0*((na(1,2)*na(3,3)) + (na(1,3)*na(3,2)))
      n23(6) = 0.5d0*((na(2,2)*na(3,3)) + (na(2,3)*na(3,2)))
!
!    Calculate the fourth-order eigenvector dyadic products
!
      call m_dyad(n11, n11, n11dy11)
      call m_dyad(n22, n11, n22dy11)
      call m_dyad(n33, n11, n33dy11)
      call m_dyad(n11, n22, n11dy22)
      call m_dyad(n22, n22, n22dy22)
      call m_dyad(n33, n22, n33dy22)
      call m_dyad(n11, n33, n11dy33)
      call m_dyad(n22, n33, n22dy33)
      call m_dyad(n33, n33, n33dy33)
!
      call m_dyad(n12, n12, n12dy12)
!
      call m_dyad(n13, n13, n13dy13)
!
      call m_dyad(n23, n23, n23dy23)
!
!    Calculate the dyad of Iden and Iden, Idy(6,6)
!
      call m_dyad(Iden, Iden, Idy)
!
!    Calculate the symmetric dyad of Iden and Iden, Id4(6,6)
!
      call m_symdy(Iden, Iden, Id4)
!
!    Compute the isochoric tangent modulus: a_iso(6,6)
!
!     First calculate terms independent of eigenvalue similarity
!
      do k1=1,6
       do k2=1,6
        a_iso(k1,k2) = 
     &       ((Ga(1,1)*n11dy11(k1,k2)) + (Ga(1,2)*n11dy22(k1,k2)) + 
     &        (Ga(1,3)*n11dy33(k1,k2)) + (Ga(2,1)*n22dy11(k1,k2)) + 
     &        (Ga(2,2)*n22dy22(k1,k2)) + (Ga(2,3)*n22dy33(k1,k2)) + 
     &        (Ga(3,1)*n33dy11(k1,k2)) + (Ga(3,2)*n33dy22(k1,k2)) + 
     &        (Ga(3,3)*n33dy33(k1,k2)))- 
     &        (two*((B_a(1)*n11dy11(k1,k2))+(B_a(2)*n22dy22(k1,k2)) + 
     &        (B_a(3)*n33dy33(k1,k2))))
       end do
      end do
!
!     If la(1)=la(2) AND If la(2)=la(3)
!
      if ((abs(la(1)-la(2))).lt.tol1) then
       if ((abs(la(2)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(2) AND If la(1)=la(3)
!
       else if ((abs(la(1)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(2)
!
       else
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If la(1)=la(3) AND If la(2)=la(3)
!
      else if((abs(la(1)-la(3))).lt.tol1) then
        if ((abs(la(2)-la(3))).lt.tol1) then
        do k1=1,6
         do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2)+ 
     &        ((((la2(1)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(1,2)))*(two*n12dy12(k1,k2))) + 
     &        ((((la2(2)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(2,1)))*(two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
          end do
        end do
!
!     If la(1)=la(3)
!
        else
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((la2(1)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(1,3)))*(two*n13dy13(k1,k2))) + 
     &        ((((la2(3)*(1/la2(1)))*((half*Ga(1,1))-B_a(1))) - 
     &        (half*Ga(3,1)))*(two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
        end do
       end do
       end if
!
!     If la(2)=la(3)
!
      else if((abs(la(2)-la(3))).lt.tol1) then
       do k1=1,6
        do k2=1,6
          a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((la2(2)*(1/la2(3)))*((half*Ga(3,3))-B_a(3))) - 
     &        (half*Ga(2,3)))*(two*n23dy23(k1,k2))) + 
     &        ((((la2(3)*(1/la2(2)))*((half*Ga(2,2))-B_a(2))) - 
     &        (half*Ga(3,2)))*(two*n23dy23(k1,k2)))
        end do
       end do
!
!    If eigenvalues are unique
!
      else
       do k1=1,6
        do k2=1,6
         a_iso(k1,k2) = a_iso(k1,k2) + 
     &        ((((B_a(2)*la2(1))-(B_a(1)*la2(2)))/(la2(2)-la2(1)))*
     &        (two*n12dy12(k1,k2))) + 
     &        ((((B_a(1)*la2(2))-(B_a(2)*la2(1)))/(la2(1)-la2(2)))*
     &        (two*n12dy12(k1,k2))) + 
!
     &        ((((B_a(3)*la2(1))-(B_a(1)*la2(3)))/(la2(3)-la2(1)))*
     &        (two*n13dy13(k1,k2))) + 
     &        ((((B_a(1)*la2(3))-(B_a(3)*la2(1)))/(la2(1)-la2(3)))*
     &        (two*n13dy13(k1,k2))) + 
!
     &        ((((B_a(3)*la2(2))-(B_a(2)*la2(3)))/(la2(3)-la2(2)))*
     &        (two*n23dy23(k1,k2))) + 
     &        ((((B_a(2)*la2(3))-(B_a(3)*la2(2)))/(la2(2)-la2(3)))*
     &        (two*n23dy23(k1,k2)))
!
        end do
       end do
      end if
!
!    Compute the volumetric tangent modulus: a_vol(6,6)
!    
      do k1 = 1, 6
       do k2 = 1, 6
        a_vol(k1,k2) = ((J*(kc+(J*kd))*Idy(k1,k2))) -
     1                  (two*kc*J*Id4(k1,k2))
       end do
      end do
!
!    Calculate symmetric dyadic products (geometric tangent contributions)
!
      call m_symdy(Iden,CST,IddyCST)
      call m_symdy(CST,Iden,CSTdyId)
!
!     Calculate the spatial tangent modulus, kEtens(6,6)
!
      CSTT(1,1)=CST(1)
      CSTT(2,2)=CST(2)
      CSTT(3,3)=CST(3)
      CSTT(2,3)=CST(4)
      CSTT(1,3)=CST(5)
      CSTT(1,2)=CST(6)
      CSTT(2,1)=CSTT(1,2)
      CSTT(3,1)=CSTT(1,3)
      CSTT(3,2)=CSTT(2,3)
      
!     Eigenvalues Cauchy Green Stress !!!!!
      CSTM(1,1)=CST(1)
      CSTM(2,2)=CST(2)
      CSTM(3,3)=CST(3)
      CSTM(2,3)=CST(4)
      CSTM(1,3)=CST(5)
      CSTM(1,2)=CST(6)
      CSTM(2,1)=CSTM(1,2)
      CSTM(3,1)=CSTM(1,3)
      CSTM(3,2)=CSTM(2,3)

      call DSYEVJ3(CSTM, nCST, laCST)

!    #ac _ calculation of max principal Cauchy-Green stress
      locmax=MAXLOC(laCST, DIM=1)
      maxpCST=laCST(locmax)
      
!    #ac _ calculation of min principal Cauchy-Green stress
      locmin=MINLOC(laCST, DIM=1)
      minpCST=laCST(locmin)
      
      do k1=1,6
       do k2=1,6
        i=xim(k1)
       jj=xjm(k1)
        k=xim(k2)
        l=xjm(k2)
        kEtens(k1,k2) = (1/J)*(a_iso(k1,k2) + a_vol(k1,k2)) +
     &        (IddyCST(k1,k2)+CSTdyId(k1,k2)) +
     &        0.5*(I2A(i,k)*CSTT(jj,l) + I2A(i,l)*CSTT(jj,k) +
     &        I2A(jj,k)*CSTT(i,l) + I2A(jj,l)*CSTT(i,k))
                  
       end do
      end do
!

!     If axisymmetric or plane strain, use axtens(4,4)
!
      do k1=1,4
       do k2=1,4
        axtens(k1,k2) = kEtens(k1,k2)
       end do
      end do
!
      return 
      end subroutine la_sub
!***********************************************************************
      subroutine bdet(mat_A, nshr, DET)
!
!    This subroutine calculates the determinant [DET] of a 3x3 
!    matrix [mat_A].
!
      implicit none
      real(8) :: DET, mat_A(3,3)
      integer, intent(IN) :: nshr
!
      DET = mat_A(1,1)*mat_A(2,2)*mat_A(3,3) -
     1      mat_A(1,2)*mat_A(2,1)*mat_A(3,3)
      if (nshr.eq.3) then
          DET = DET + mat_A(1,2)*mat_A(2,3)*mat_A(3,1)
     1          + mat_A(1,3)*mat_A(2,1)*mat_A(3,2)
     2          - mat_A(1,1)*mat_A(2,3)*mat_A(3,2)
     3          - mat_A(1,3)*mat_A(2,2)*mat_A(3,1)
      end if
      return
      end subroutine bdet
!**********************************************************************
      subroutine m_dyad(mat_B, mat_C, dyad_D)
!
!    This subroutine calculates the dyadic product [dyad_D] of two
!    symmetric 2nd-order tensors [mat_B & mat_C] written in Voigt form.
!
      implicit none
      real(8) :: mat_B(6), mat_C(6), dyad_D(6,6)
      integer :: k1,k2
!
      do k1=1, 6
          do k2=1, 6
              dyad_D(k1,k2) = mat_B(k1)*mat_C(k2)
          end do
      end do
      return
      end subroutine m_dyad
!**********************************************************************
      subroutine m_symdy(mat_B, mat_C, symdy_D)
!
!    This subroutine calculates the symmetric dyadic product [symdy_D]
!    of two symmetric 2nd-order tensors [mat_B & mat_C] in Voigt 
!    notation for an isotropic material.
!
      implicit none
      real(8) mat_B(6),mat_C(6),mat_B2(6,6),mat_C2(6,6),symdy_D(6,6)
      integer k1,k2
!
      mat_B2(1,1) = mat_B(1)*mat_C(1)
      mat_B2(1,2) = mat_B(4)*mat_C(4)
      mat_B2(1,3) = mat_B(5)*mat_C(5)
      mat_B2(1,4) = mat_B(1)*mat_C(4)
      mat_B2(1,5) = mat_B(1)*mat_C(5)
      mat_B2(1,6) = mat_B(4)*mat_C(5)
      mat_B2(2,2) = mat_B(2)*mat_C(2)
      mat_B2(2,3) = mat_B(6)*mat_C(6)
      mat_B2(2,4) = mat_B(4)*mat_C(2)
      mat_B2(2,5) = mat_B(4)*mat_C(6)
      mat_B2(2,6) = mat_B(2)*mat_C(6)
      mat_B2(3,3) = mat_B(3)*mat_C(3)
      mat_B2(3,4) = mat_B(5)*mat_C(6)
      mat_B2(3,5) = mat_B(5)*mat_C(3)
      mat_B2(3,6) = mat_B(6)*mat_C(3)
      mat_B2(4,4) = mat_B(1)*mat_C(2)
      mat_B2(4,5) = mat_B(1)*mat_C(6)
      mat_B2(4,6) = mat_B(4)*mat_C(6)
      mat_B2(5,5) = mat_B(1)*mat_C(3)
      mat_B2(5,6) = mat_B(4)*mat_C(3)
      mat_B2(6,6) = mat_B(2)*mat_C(3)
      do k1=1, 6
        do k2=1, k1-1
          mat_B2(k1,k2) = mat_b2(k2,k1)
        end do
      end do
!
      mat_C2(1,1) = mat_B(1)*mat_C(1)
      mat_C2(1,2) = mat_B(4)*mat_C(4)
      mat_C2(1,3) = mat_B(5)*mat_C(5)
      mat_C2(1,4) = mat_B(4)*mat_C(1)
      mat_C2(1,5) = mat_B(5)*mat_C(1)
      mat_C2(1,6) = mat_B(5)*mat_C(4)
      mat_C2(2,2) = mat_B(2)*mat_C(2)
      mat_C2(2,3) = mat_B(6)*mat_C(6)
      mat_C2(2,4) = mat_B(2)*mat_C(4)
      mat_C2(2,5) = mat_B(6)*mat_C(4)
      mat_C2(2,6) = mat_B(6)*mat_C(2)
      mat_C2(3,3) = mat_B(3)*mat_C(3)
      mat_C2(3,4) = mat_B(6)*mat_C(5)
      mat_C2(3,5) = mat_B(3)*mat_C(5)
      mat_C2(3,6) = mat_B(3)*mat_C(6)
      mat_C2(4,4) = mat_B(4)*mat_C(4)
      mat_C2(4,5) = mat_B(5)*mat_C(4)
      mat_C2(4,6) = mat_B(5)*mat_C(2)
      mat_C2(5,5) = mat_B(5)*mat_C(5)
      mat_C2(5,6) = mat_B(5)*mat_C(6)
      mat_C2(6,6) = mat_B(6)*mat_C(6)
      do k1=1, 6
        do k2=1, k1-1
          mat_C2(k1,k2) = mat_C2(k2,k1)
        end do
      end do
!
      do k1=1, 6
          do k2=1, 6
              symdy_D(k1,k2) = (0.5d0*(mat_B2(k1,k2)+mat_C2(k1,k2)))
          end do
      end do
      return
      end subroutine m_symdy
!**********************************************************************
* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrices
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X)) !gt means >
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
      END SUBROUTINE
c***********************************************************************
      subroutine invertmat33(a,deta,ainv)
      implicit none
      real*8   a(3,3),deta
      real*8   detainv
      real*8   ainv(3,3)

      detainv=1.0d0/deta
      ainv(1,1)=(+a(2,2)*a(3,3)-a(2,3)*a(3,2))*detainv
      ainv(1,2)=(-a(1,2)*a(3,3)+a(1,3)*a(3,2))*detainv
      ainv(1,3)=(+a(1,2)*a(2,3)-a(1,3)*a(2,2))*detainv
      ainv(2,1)=(-a(2,1)*a(3,3)+a(2,3)*a(3,1))*detainv
      ainv(2,2)=(+a(1,1)*a(3,3)-a(1,3)*a(3,1))*detainv
      ainv(2,3)=(-a(1,1)*a(2,3)+a(1,3)*a(2,1))*detainv
      ainv(3,1)=(+a(2,1)*a(3,2)-a(2,2)*a(3,1))*detainv
      ainv(3,2)=(-a(1,1)*a(3,2)+a(1,2)*a(3,1))*detainv
      ainv(3,3)=(+a(1,1)*a(2,2)-a(1,2)*a(2,1))*detainv

      return
      end

c***********************************************************************
      subroutine matI2(I2, ndim)
      implicit none
      integer ndim, i
      real*8 I2(ndim,ndim)

      I2 = 0.0
      do i = 1,ndim
            I2(i,i) = 1.0
      enddo

      return
      end
c*******E*N*D****O*F****S*U*B*R*O*U*T*I*N*E*****************************




