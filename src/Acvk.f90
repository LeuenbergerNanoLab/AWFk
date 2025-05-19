


   Program Acvk1
    integer, parameter    :: nv = 6                   ! valence bands
    integer, parameter    :: nc = 6                   ! conduction bands
    integer, parameter    :: nk = 324                 ! k-points 
    integer               :: iv,ic,ik
    real(8)               :: Abs,phi
    real(8)               :: Acvk_re(nv,nc,nk)        ! envelope functions from Yambo re and im parts
    real(8)               :: Acvk_im(nv,nc,nk)
    real(8)               :: kp(3,nk)                 ! k-points in cartesian coordinates
    real(8)               :: abs0                     ! max value of abs|Acvk|
    character(1)          :: fstr
    character(15)         :: name

    call getarg(1,name)                               ! read name of the file
    print *,'read file ',name
    call read_Acvk(name,Acvk_re,Acvk_im,nc,nv,nk)
    call read_k_points(nk,kp)
    call calc_max(Acvk_re,Acvk_im,nc,nv,nk,kp)
    call calc_max_abs(Acvk_re,Acvk_im,nc,nv,nk,abs0,kp)
!    call write_tecplot(kp,Acvk_re,Acvk_im,nc,nv,nk)
    call write_gnuplot(kp,Acvk_re,Acvk_im,nc,nv,nk,abs0)
   end Program Acvk1



   subroutine read_Acvk(name,Acvk_re,Acvk_im,nc,nv,nk)
    character(*)          :: name
    integer               :: nc,nv,nk
    real(8)               :: Acvk_re(nv,nc,nk)
    real(8)               :: Acvk_im(nv,nc,nk)
    open(unit=1,file=name)
     read(1,*)
     do ik = 1,nk
      do ic = 1,nc
       do iv = 1,nv
        read(1,*) Acvk_re(iv,ic,ik),Acvk_im(iv,ic,ik)
       enddo
      enddo
     enddo 
    close(unit=1)
   end subroutine read_Acvk



   subroutine write_gnuplot(kp,Acvk_re,Acvk_im,nc,nv,nk,abs0)
    use QSHEP_MOD
    integer               :: nc,nv,nk
    real(8)               :: Acvk_re(nv,nc,nk)
    real(8)               :: Acvk_im(nv,nc,nk)
    real(8)               :: kp(3,nk)           ! k-points in cartesian coordinates
    character(3)          :: A3
    real(8), allocatable  :: A11r(:),A11i(:)
    real(8)               :: abs0
    integer               :: ik
    real(8)               :: kxp,kyp
    real(8), allocatable  :: kx(:)
    real(8), allocatable  :: ky(:)
    character(1)          :: fstr
    character(7)          :: color
    real(8)               :: A11rs,A11is
    integer               :: I2
    real(8)               :: kpK(2,6)
    integer               :: ip1,jp1
    integer               :: I20
    logical               :: lh
    
! QSHEP
    integer, parameter   :: NR = 10
    integer              :: LCELLr(NR,NR)
    integer              :: LCELLi(NR,NR)
    real(8)              :: XMINr,YMINr,DXr,DYr,RMAXr
    real(8)              :: XMINi,YMINi,DXi,DYi,RMAXi
    integer              :: IER
    real(8)              :: QX,QY                      ! working 
    integer, allocatable :: LNEXTr(:),LNEXTi(:)
    real(8), allocatable :: RSQr(:),RSQi(:)
    real(8), allocatable :: Ashepr(:,:),Ashepi(:,:)
    
    allocate(LNEXTr(nk))                              ! working arrays for QSHEP
    allocate(RSQr(nk))                                ! working arrays for QSHEP
    allocate(Ashepr(5,nk))                            ! working arrays for QSHEP
    allocate(LNEXTi(nk))                              ! working arrays for QSHEP
    allocate(RSQi(nk))                                ! working arrays for QSHEP
    allocate(Ashepi(5,nk))                            ! working arrays for QSHEP
    
    allocate(A11r(nk))
    allocate(A11i(nk))
    allocate(kx(nk))
    allocate(ky(nk))
    
    I20 = 20      ! cut intensity

! hexagonal K-points
       kpK(1:2,1) = (/ 0.296101d0, 0.512862d0/)
       kpK(1:2,2) = (/ 2*0.296101d0, 0.d0/)
       kpK(1:2,3) = (/ 0.296101d0,-0.512862d0/)
       kpK(1:2,4) = (/-0.296101d0,-0.512862d0/)
       kpK(1:2,5) = (/-2*0.296101d0, 0.d0/)
       kpK(1:2,6) = (/-0.296101d0, 0.512862d0/)

    do ik=1,nk
        kx(ik) = kp(1,ik)
        ky(ik) = kp(2,ik)
    enddo
    do iv = 1,nv
     print *,'**** iv=',iv
     do ic = 1,nc
      print *,'ic=',ic
      A3 = 'A'//fstr(iv)//fstr(ic)
      open(unit=2,file=A3//'.gnu')
       write(2,1)
       write(2,2) A3
       write(2,3)
       write(2,4) A3,200,650
       do ik=1,nk
        A11r(ik) = Acvk_re(iv,ic,ik)
        A11i(ik) = Acvk_im(iv,ic,ik)
       enddo

       call QSHEP2(nk,kx,ky,A11r,13,10,NR,LCELLr,LNEXTr,XMINr,YMINr,DXr,DYr,RMAXr,RSQr,Ashepr,IER)
       call QSHEP2(nk,kx,ky,A11i,13,10,NR,LCELLi,LNEXTi,XMINi,YMINi,DXi,DYi,RMAXi,RSQi,Ashepi,IER)

       do ip=1,800
        do jp=1,800
         call convert_p_to_k(ip,jp,kxp,kyp)
         call check_inside_hexagon(kyp,kxp,kpK(1,2),lh)     ! kyp,kxp to orient hexagon horizontally
         if(lh) then
          call QS2GRD(kxp,kyp,nk,kx,ky,A11r,NR,LCELLr,LNEXTr,XMINr,YMINr,DXr,DYr,RMAXr,RSQr,Ashepr,A11rs,QX,QY,IER)
          call QS2GRD(kxp,kyp,nk,kx,ky,A11i,NR,LCELLi,LNEXTi,XMINi,YMINi,DXi,DYi,RMAXi,RSQi,Ashepi,A11is,QX,QY,IER)
          call calc_color(A11rs,A11is,abs0,color,I2)
          if(I2 > I20) then
           write(2,9) dfloat(I2)/255.d0
           write(2,10) ip,jp,color
          endif
         endif
        enddo 
       enddo
       ! plot 6 K-points
       do ik=1,6              ! 6 K-points
        call convert_k_to_p(kpK(1,ik),kpK(2,ik),ip,jp)
        if(ik==1) then
         ip10 = ip
         jp10 = jp
        endif
        color = '#000000'
        write(2,9) 1.d0
        write(2,13) ip,jp,color
        if(ik>=2.and.ik<6) then
         write(2,12) ip,jp,ip1,jp1,color
        elseif(ik==6) then
         write(2,12) ip,jp,ip1,jp1,color
         write(2,12) ip,jp,ip10,jp10,color
        endif         
        ip1 = ip
        jp1 = jp
       enddo
       call write_scale
      write(2,11)
      close(unit=2)
     enddo
    enddo
    deallocate(kx)
    deallocate(ky)
    deallocate(A11r)
    deallocate(A11i)
    deallocate(LNEXTr)
    deallocate(LNEXTi)
    deallocate(RSQr)
    deallocate(RSQi)
    deallocate(Ashepr)
    deallocate(Ashepi)
 1  format('set terminal pdfcairo size 3,3 linewidth 1 font "Times,10"')
 2  format('set output "',A3,'.pdf"'/)
 3  format('set xrange [1:800]'/ &
           'set yrange [1:800]'/ &
           'set size ratio -1'/ &
           'unset key'/ &
           'set notics'/ &
           'unset border'/)
 4  format('set label "',A3,'" at ',I4,',',I4,' right')
 5  format(2I4,2x,A10)
 6  format('plot "+" using (',I3,'):(',I3,'):(',A10,') lc rgb variable pt 7 ps 2,\')
 7  format('     "+" using (',I3,'):(',I3,'):(',A10,') lc rgb variable pt 7 ps 2,\')
 8  format('     "+" using (',I3,'):(',I3,'):(',A10,') lc rgb variable pt 7 ps 2')
 9  format('set style fill transparent solid ',F8.5,' noborder')
10  format('set object circle at ',I4,',',I4,' radius 1 fc rgb "',A7,'"') 
11  format('plot -1 notitle')
12  format('set arrow from ',I4,',',I4,' to ',I4,',',I4,' nohead lw 1 lc rgb "',A7,'"') 
13  format('set object circle at ',I4,',',I4,' radius 3 fc rgb "',A7,'"') 
   end subroutine write_gnuplot



   subroutine write_scale
    real(8), parameter    :: pi=3.141592653589793d0
    integer               :: jc,jl
    character(6)          :: color
    real(8)               :: phi
    integer               :: B1x,B2y,sx,sy
    B1x = 770   ! left bottom corner of the box 21 x 512
    B2y = 150
    sx = 21     ! size of the box
    sy = 512
    write(2,12) B1x,   B2y,    B1x,   B2y+sy    ! plot black box
    write(2,12) B1x,   B2y+sy, B1x+sx,B2y+sy
    write(2,12) B1x+sx,B2y+sy, B1x+sx,B2y
    write(2,12) B1x+sx,B2y,    B1x,   B2y
    write(2,9) B1x-10,B2y+sy,B1x-10,B2y+sy/2,B1x-10,B2y
    jl = B2y + 1
    do jc = -255,255                            ! 511 colors
     phi = dfloat(jc)/255.d0*pi 
     call calc_color_phi(phi,color)             ! define color red-blue for phase
     do ip = 1,20
      write(2,10) B1x+ip,jl,color
     enddo 
     jl = jl + 1
    enddo
 9  format('set label "+{/Symbol p}" at ',I4,',',I4,' right'/ &
           'set label "0"            at ',I4,',',I4,' right'/ &
           'set label "-{/Symbol p}" at ',I4,',',I4,' right')
10  format('set object circle at ',I4,',',I4,' radius 1 fc rgb "#',A6,'"') 
12  format('set arrow from ',I4,',',I4,' to ',I4,',',I4,' nohead lw 1 lc rgb "#000000"'/) 
   end subroutine write_scale



   subroutine write_tecplot(kp,Acvk_re,Acvk_im,nc,nv,nk)
    integer               :: nc,nv,nk
    real(8)               :: Acvk_re(nv,nc,nk)
    real(8)               :: Acvk_im(nv,nc,nk)
    real(8)               :: kp(3,nk)           ! k-points in cartesian coordinates
    character(1)          :: fstr
    do iv = 1,nv
     print *,'**** iv=',iv
     do ic = 1,nc
      print *,'ic=',ic
      open(unit=1,file='A'//fstr(iv)//fstr(ic)//'.dat')
       write(1,1) nk
       do ik=1,nk
        abs = dsqrt(Acvk_re(iv,ic,ik)**2+Acvk_im(iv,ic,ik)**2)
        phi = datan2(Acvk_im(iv,ic,ik),Acvk_re(iv,ic,ik))
        write(1,2) kp(1:2,ik),Acvk_re(iv,ic,ik),Acvk_im(iv,ic,ik),abs,phi
       enddo
      close(unit=1)
     enddo
    enddo
 1  format('VARIABLES = "kx", "ky", "Re", "Im", "Abs", "Phi"'/'ZONE I=',I4,' F=POINT')
 2  format(2F11.6,4E15.5)    
   end subroutine write_tecplot



     character(len=1) function fstr(k)                             !   Convert an integer to character*7
      integer, intent(in) :: k
      write (fstr,'(I1)') k
     end function fstr



   subroutine read_k_points(nk,kp)
    integer         :: nk
    real(8)         :: kp(3,nk)
    integer         :: ik
    print *,'read kpoints.dat'
    open(unit=2,file='kpoints.dat')
     read(2,*)
     read(2,*)
     read(2,*)
     do ik=1,nk
      read(2,*)
      read(2,*)
      read(2,*) kp(1:3,ik)
 !     print *, kp(1:3,ik)
      read(2,*)
      read(2,*)
     enddo
     close(unit=2)
   end subroutine read_k_points



   subroutine calc_max(Acvk_re,Acvk_im,nc,nv,nk,kp)
    integer               :: nc,nv,nk
    real(8)               :: kp(3,nk)
    real(8)               :: Acvk_re(nv,nc,nk)
    real(8)               :: Acvk_im(nv,nc,nk)
    real(8)               :: max,abs,max_re,max_im,abs_re,abs_im
    integer               :: ic,iv,ik
    logical               :: lh
    print 1
    do iv=1,nv
     do ic=1,nc
      max    = 0.d0
      max_re = 0.d0
      max_im = 0.d0
      do ik=1,nk
       call check_inside_hexagon(kp(1,ik),kp(2,ik),2*0.296101d0,lh)    ! check only points inside BZ
       if(lh) then
        abs = dsqrt(Acvk_re(iv,ic,ik)**2+Acvk_im(iv,ic,ik)**2)
        abs_re = dabs(Acvk_re(iv,ic,ik))
        abs_im = dabs(Acvk_im(iv,ic,ik))
        if(max    .lt. abs)    max = abs
        if(max_re .lt. abs_re) max_re = abs_re
        if(max_im .lt. abs_im) max_im = abs_im
       endif
      enddo 
      print 2,iv,ic,max,max_re,max_im
     enddo
    enddo
 1  format('  iv   ic      max_abs        max_re         max_im')
 2  format(2I2,3F25.17)
   end subroutine calc_max



   subroutine calc_max_abs(Acvk_re,Acvk_im,nc,nv,nk,abs0,kp)
    integer               :: nc,nv,nk
    real(8)               :: Acvk_re(nv,nc,nk)
    real(8)               :: Acvk_im(nv,nc,nk)
    real(8)               :: max,abs
    integer               :: ic,iv,ik
    real(8)               :: abs0
    real(8)               :: kp(3,nk)
    logical               :: lh
    max    = 0.d0
    do iv=1,nv
     do ic=1,nc
      do ik=1,nk
       call check_inside_hexagon(kp(1,ik),kp(2,ik),2*0.296101d0,lh)    ! check only points inside BZ
       if(lh) then
        abs = dsqrt(Acvk_re(iv,ic,ik)**2+Acvk_im(iv,ic,ik)**2)
        if(max    .lt. abs)    max = abs
       endif
      enddo 
     enddo
    enddo
    abs0 = max
    print *,'abs0=',abs0
   end subroutine calc_max_abs



   subroutine calc_color(A11r,A11i,abs0,color,I2)
    real(8)            :: A11r,A11i
    character(7)      :: color
    character(6)       :: color1
    character(2)       :: color2
    real(8)            :: abs0
    real(8)            :: abs,phi
    integer            :: I2
    call calc_abs_phi(A11r,A11i,abs,phi)
    call calc_color_phi(phi,color1)              ! define color red-blue for phase
    call calc_color_abs(abs,abs0,color2,I2)      ! define intensity (transparency) for abs value
!    color = '0x'//color2//color1                ! combine transparency + RGB color
    color = '#'//color1                          ! combine transparency + RGB color
   end subroutine calc_color



   subroutine calc_abs_phi(A11r,A11i,abs,phi)
    real(8)             :: A11r,A11i
    real(8)             :: abs,phi
    abs = dsqrt(A11r**2+A11i**2)
    phi = datan2(A11i,A11r)
   end subroutine calc_abs_phi



   subroutine calc_color_phi(phi,color1)               ! define color red-blue for phase
    real(8)            :: phi
    character(6)       :: color1
    character(2)       :: R,G,B
    integer            :: I
    character(2)       :: col
    G = 'FF'
    R = 'FF'
    B = 'FF'
    if(phi>0.0001d0) then
     call convert_phi_to_I(phi,I)                      ! convert phi to intensity of red from 0 to 255
     call convert_I_to_color(I,col)                    ! convert intensity of color to hexadecimal value
     if(I>=0.and.I<=127) then
      R = 'FF'
      G = col
      B = col
     elseif(I>127.and.I<=255) then
      R = col
      G = '00'
      B = '00'
     endif 
    elseif(phi<-0.0001d0) then
     call convert_phi_to_I(-phi,I)                     ! convert phi to intensity of blue from 0 to 255
     call convert_I_to_color(I,col)                    ! convert intensity of color to hexadecimal value
     if(I>=0.and.I<=127) then
      R = col
      G = col
      B = 'FF'
     elseif(I>127.and.I<=255) then
      R = '00'
      G = '00'
      B = col
     endif
    else
     R = 'FF'
     B = 'FF'
    endif
    color1 = R//G//B
   end subroutine calc_color_phi



   subroutine convert_phi_to_I(phi,I)                  ! convert phi to intensity of color from 0 to 255 scale
    real(8)            :: phi
    integer            :: I
    real(8), parameter :: pi=3.141592653589793d0
    I = int((dabs(phi)/pi)*255.d0+0.1d0)
    if(I<0.or.I>255) then
     print *,'convert_phi_to_I   I=',I
     print *,'convert_phi_to_I   phi=',phi
     stop
    endif
   end subroutine convert_phi_to_I



   subroutine convert_I_to_color(I,col)                ! convert intensity of color to hexadecimal value
    integer             :: I,I1
    character(2)        :: col
    if(I>=0.and.I<=127) then
     I1 = 255 - 2*I
     write (col,'(Z2)') I1
    elseif(I>127.and.I<=255) then
!     I1 = int(dabs(255.d0-2.d0*I)/1.83d0)
     I1 = 382 - I
     write (col,'(Z2)') I1
    else
     print *,'convert_I_to_color   I=',I
     print *,'convert_I_to_color   col=',col
     stop
    endif 
    if(col(1:1) == ' ') col(1:1) = '0'
   end subroutine convert_I_to_color



   subroutine calc_color_abs(abs,abs0,color2,I2)        ! define intensity (transparency) for abs value
    real(8)             :: abs,abs0
    character(2)        :: color2
    integer             :: I,I2
    call convert_abs_to_I(abs,abs0,I)                   ! convert abs to intensity of blue from 0 to 255
    I2 = I
    call convert_I_to_color(I,color2)                   ! convert intensity of color to hexadecimal value
   end subroutine calc_color_abs



   subroutine convert_abs_to_I(abs,abs0,I)              ! convert abs to intensity from 0 to 255
    real(8)             :: abs,abs0
    integer             :: I
    I = int((abs/abs0)*255.d0)
    if(I>255) I = 255
   end subroutine convert_abs_to_I



   subroutine convert_p_to_k(ip,jp,kx,ky)               ! convert pixel positions (ip,jp) to (kx,ky)
    integer            :: ip,jp
    real(8)            :: kx,ky,dk
    real(8), parameter :: Size = 0.8d0
    dk = Size/400.d0
    kx = (ip-400)*dk
    ky = (jp-400)*dk
   end subroutine convert_p_to_k



   subroutine convert_k_to_p(kx,ky,ip,jp)               ! convert (kx,ky) coordinates to pixel positions (ip,jp) 
    integer            :: ip,jp
    real(8)            :: kx,ky,dk
    real(8), parameter :: Size = 0.8d0
    dk = Size/400.d0
    ip = kx/dk+400
    jp = ky/dk+400
   end subroutine convert_k_to_p



 subroutine check_inside_hexagon(x,y,R6,lh)          ! check if (x,y) inside hexagone (vertically oriented)
  real(8)        :: x, y                             ! (y,x) - will be horizontally oriented
  real(8)        :: R6                   ! external radius of hexagone
  real(8)        :: R33                  ! internal radius of hexagone 
  real(8)        :: q, r, s, max_abs
  logical        :: lh

  lh = .false.
  R33 = dsqrt(3.d0)/3.d0*R6 

  q = (2.0D0/3.0D0) * (x / R33)
  r = (-1.0D0/3.0D0) * (x / R33) + (dsqrt(3.0D0)/3.0D0) * (y / R33)
  s = -q - r
  max_abs = dmax1(dabs(q), dabs(r), dabs(s))

  if (max_abs <= 1.0D0) then
   lh = .true.             ! x,y is INSIDE the hexagon
  else
   lh = .false.            ! x,y is OUTSIDE the hexagon
  end if

 end subroutine check_inside_hexagon



