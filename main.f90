module dannie
integer(4) i,j,mi,mj,Nv,Nstep,step
real(8) Temp,dr
real(8) randalfa !случайный поворот относительно оси
real(8) randbeta !случайный поворот относительно оси
real(8) randgamma !случайный поворот относительно оси
character(100) Naz, workdir,Nazv, NazF(10)
character(204) Nazv2
real(8) sigmin, sigmax, epsmin, epsmax,MMi
real(8) xm(100),ym(100),zm(100),rm(100),ix(100),iy(100),iz(100),ir(100)
real(8) iep(100),isi(100),sumXm,sumYm,sumZm,sumRm
real(8) iq(100)
real(8) MM(100)
integer(4) Nm
real(8) cxm,cym,czm
real(8) epsi(100),sigma(100)
real(8) alfa(100) !альфа букингема
real(8) alf(30,30)
real(8) rz,rz2,rz6,expp
real(8) bk1(30,30)
real(8) bk2(30,30)
real(8) dbk1(30,30)
real(8) dbk2(30,30)
real(8) xa1(100),xa2(100),ya1(100),ya2(100),za1(100),za2(100)
character(4) ill(100),labelm(100)
real(8) siga(100,100),epsa(100,100)
real(8) dxa,dya,dza
real(8) sintx,sinty,sintz,costx,costy,costz
real(8) xv,yv,zv
real(8) b2dif,b2dif2
real(8) b2int(1000000),rdr(1000000)
integer(4) iri
real(8) b2, Rcut2m,en2
integer(4) ek,sk,tk,sigk,epsk,ak,alfk
real(8) ekt, alfmin,alfmax
real(8) ia(100)
character(20) estr,sstr,astr
character(200) outf,iname
real(8) razna,razno,b2e,pich,abst,otnt,maxra,maxro
real(8) kk1(30,30),kk2(30,30),kk3(30,30),kk4(30,30),siga6
real(8) kkm1(30,30),kkm2(30,30)
integer(4) schet
integer(4) ptip
!
integer(4) bind_kol(30)
integer(4) bind_koli,bind_ci(30),bind1i(30),bind2i(30)
integer(4) bind1ki(30,30),bind2ki(30,30),bind_ki(30),bind_fii(30)
integer(4) bind_cent(30,30)
integer(4) bind1(30,30)
integer(4) bind2(30,30)
integer(4) bind1k(30,30,30)
integer(4) bind2k(30,30,30)
real(8) bind_k(30,30)
real(8) bind_fi(30,30)
real(8) r1_bind,r2_bind,cos_bind
real(8) t_bind_en,bind_en,tors_en
!
integer(4) sample_kol
real(8),allocatable:: xa(:,:),ya(:,:),za(:,:),bind_u(:,:)
real(8),allocatable:: xx(:)
real(8),allocatable:: xmc(:),ymc(:),zmc(:)
real(8),allocatable:: xmc_do(:),ymc_do(:),zmc_do(:),exp_en(:)
real(8) mc_en,mc_en_posle,sum_exp
!!
integer(4) torsi !количество торсионных связей в молекуле
integer(4) fsti(30), seci(30), nroti(100)
integer(4) rotai(100,100), torskoefi(100,5)
integer(4) NMCI,jto,i_s
real(8) torsCosi(30),delta,bind_u_mc
real(8),allocatable:: bind_m(:)
real(8) fmax
real(8) rmax(30,30)
integer(4) zi
real(8) modif
real(8) sign6i
real(8) dekt
real(8) minr,dminr
real(8) pmin1,pmin2
character(100) cdir
character(100) tempst
end module
!**********************************************************************
module generator
 integer(4) n1, n2, n3 !случаные числа
 integer(4) m1,a1,b1,a2,b2,m3,m2 !параметры генератора
 integer(4) max1, max2,max3 !максимумы генераторов
 real(8) randmass(500) !массив случайных чисел
 real(8) outrand
end module
!**********************************************************************
program virial
use dannie
!Ввод исходных данных
namelist /main/ptip,Naz,Tempm,Tempk,kTemp,Nstep,workdir,sigmax,sigmin,epsmin, epsmax,sigk,epsk,alfmin,alfmax,alfk
namelist /molekula/Nmi,iName,ix,iy,iz,ir,iep,ill,ia,isi,iq
open(10,file='input.txt')
 read(10,main)
close(10)
call randomn !задействование случайных чисел
PICH=4*atan(1.0)
modif=5
print *,Pich
workdir=adjustr(workdir) !рабочая папка
 write(Nazv2,'(a,a4,a)') trim(adjustl(workdir)),'/VV/',trim(adjustl(Naz))
!print *, nazv2
!print *,'ok'
!print *,tempk,tempm
open(10,file=Naz) !Nazv2
!read(10,molekula) !считывание данных по молекуле
read(10,'(a)') tempst
read(10,'(i5)') nmi
read(10,'(a)') tempst
read(10,'(a)') iname
do j=1,nmi
    read(10,'(a)') ill(j)
    read(10,'(f20.10)') ix(j)
    read(10,'(f20.10)') iy(j)
    read(10,'(f20.10)') iz(j)
    read(10,'(f20.10)') ir(j)
    read(10,'(f20.10)') iep(j)
    read(10,'(f20.10)') isi(j)
    read(10,'(f20.10)') ia(j)
    read(10,'(f20.10)') iq(j)
enddo

 MM=MMi
 Nm=Nmi
 do j=1,Nm,1
  xm(j)=ix(j)
  ym(j)=iy(j)
  zm(j)=iz(j)
  rm(j)=ir(j)
  alfa(j)=ia(j)
  epsi(j)=iep(j)
  sigma(j)=isi(j)
  labelm(j)=ill(j)
  print *,xm(j),ym(j),zm(j)
  print *,epsi(j),sigma(j)
 enddo
 i=1
 bind_kol(i)=bind_koli
 do j=1,bind_kol(i),1
  bind_cent(i,j)=bind_ci(j) !номер центра
  bind1(i,j)=bind1i(j) !количество центров с одной стороны
  bind2(i,j)=bind2i(j) !количество центров с другой
  print *,bind_cent(i,j),bind1(i,j),bind2(i,j)
  do k=1,bind1(i,j),1
   bind1k(i,j,k)=bind1ki(j,k) !номера поворачиваемых центров
   !первый центр соотвествует центру составляющему угол
   print *,bind1k(i,j,k)
  enddo
  do k=1,bind2(i,j),1
   bind2k(i,j,k)=bind2ki(j,k) !номера поворачиваемых центров
   !первый центр соотвествует центру составляющему угол
   print *, bind2k(i,j,k)
  enddo
  bind_k(i,j)=bind_ki(j)/2.0
  bind_fi(i,j)=bind_fii(j)/180.0*PICh
  print *,bind_fi(i,j),bind_k(i,j)
 enddo
close(10)
!поиск минимума потенциала
do j=1,nm
    print *,'zaryad' ,j, iq(j)
enddo

minr=1000.0
dminr=0.0001
!open(26,file='dpot.txt')
do i=1,100000
    dmr1=float(i)*dminr
    dmr2=float(i+1)*dminr
    pmin1=4.0*epsi(1)*(6.0*(sigma(1)/dmr1)**6-12.0*(sigma(1)/dmr1)**12)/dmr1-iq(1)*iq(2)/dmr1**2
    pmin2=4.0*epsi(1)*(6.0*(sigma(1)/dmr2)**6-12.0*(sigma(1)/dmr2)**12)/dmr2-iq(1)*iq(2)/dmr2**2
    write (26,'(f10.5,a,f10.5,a,f10.5)') dmr1,';',pmin1,';',pmin2
    if (pmin1*pmin2<0.0) then
        print *,dmr1,dmr2, pmin1,pmin2
        minr=(dmr1+dmr2)/2.0
    endif
enddo
!xm(1)=minr*1.0
!close(26)
!
call getcwd(cdir)
print *, cdir
print *, minr, pmin1,';',pmin2 !, sigma(1),epsi(1)
print *,'ptip', ptip
!перенос координат в  центр масс
sumXm=0.0
sumYm=0.0
sumZm=0.0
sumRm=0.0
do j=1,Nm,1
 sumXm=sumXm+xm(j)*rm(j)
 sumYm=sumYm+ym(j)*rm(j)
 sumZm=sumZm+zm(j)*rm(j)
 sumRm=sumRm+rm(j)
enddo
cxm=sumXm/sumRm
cym=sumYm/sumRm
czm=sumZm/sumRm
print *,cxm,cym,czm
!****************************************************************************
!перенос координат атомов в молекуле
do j=1,Nm,1
 xm(j)=xm(j)-cxm
 ym(j)=ym(j)-cym
 zm(j)=zm(j)-czm
enddo
do i=1,Nm,1
 xa1(i)=xm(i)
 ya1(i)=ym(i)
 za1(i)=zm(i)
 xa2(i)=xm(i)
 ya2(i)=ym(i)
 za2(i)=zm(i)
enddo
!
sample_kol=1
nmmax=nm
allocate(xa(sample_kol,nmmax))
allocate(ya(sample_kol,nmmax))
allocate(za(sample_kol,nmmax))
allocate(bind_u(sample_kol,nmmax))
allocate(xx(sample_kol))
allocate(exp_en(sample_kol))
allocate (xmc(nmmax))
allocate (ymc(nmmax))
allocate (zmc(nmmax))
allocate (xmc_do(nmmax))
allocate (ymc_do(nmmax))
allocate (zmc_do(nmmax))
allocate (bind_m(sample_kol))
!111111111111
dr=0.1
schet=0
!меняем параметры потенциала
do ak=1,alfk,1
 if (alfk==1) then
!  alfa(1)=alfmin
!  alfa(2)=alfmin
!  alfa(3)=alfmin
!  alfa(4)=alfmin
!  alfa(5)=alfmin
!  alfa(6)=alfmin
 else
!  alfa(1)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
!  alfa(2)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
!  alfa(3)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
!  alfa(4)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
!  alfa(5)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
!  alfa(6)=alfmin+(alfmax-alfmin)/float(alfk-1)*float(ak-1)
 endif
 do sk=1,sigk,1
  !sigma(1)=-7.0/150.0*alfa(1)+4.77 !6.666-0.137*alfa(1)
  !sigma(2)=sigma(1)
  !--------------------------
  if (sigk==1) then
!   sigma(1)=sigmin
!   sigma(2)=sigmin
!   sigma(3)=sigmin
!   sigma(4)=sigmin
!   sigma(5)=sigmin
!   sigma(6)=sigmin
  else
!   sigma(1)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
!   sigma(2)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
!   sigma(3)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
!   sigma(4)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
!   sigma(5)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
!   sigma(6)=sigmin+(sigmax-sigmin)/float(sigk-1)*float(sk-1)
  endif
  !print *, sigma(1)
 do ek=1,epsk,1
 ! если проверять то надо убрать
 !epsi(2)=607.22136-211.27346*sigma(2)+21.05688*sigma(2)*sigma(2)!(6.0+691.0/750.0)*alfa(1)+32.804 !-35.0*sigma(1)+239.3
 !epsi(2)=epsi(2)
 !-----------
 if (epsk==1) then
!  epsi(1)=epsmin
!  epsi(2)=epsmin
!  epsi(3)=epsmin
!  epsi(4)=epsmin
!  epsi(5)=epsmin
!  epsi(6)=epsmin
 else
!  epsi(1)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
!  epsi(2)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
!  epsi(3)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
!  epsi(4)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
!  epsi(5)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
!  epsi(6)=epsmin+(epsmax-epsmin)/float(epsk-1)*float(ek-1)
 endif
 !-----------
 !print *, epsi(1)
  do i=1,Nm
   do j=1,Nm
    siga(i,j)=(sigma(i)+sigma(j))/2.0
    epsa(i,j)=sqrt(epsi(i)*epsi(j))
    alf(i,j)=sqrt(alfa(i)*alfa(j))
    if (ptip==9) then
     sign6i=6.0**(alf(i,j)/(6.0-alf(i,j))+6.0/(6.0-alf(i,j)))/(alf(i,j)**(6.0/(6.0-&
     &alf(i,j)))*6.0**(alf(i,j)/(6.0-alf(i,j)))-alf(i,j)**(alf(i,j)/(6.0-alf(i,j)))*&
     &6.0**(6.0/(6.0-alf(i,j))))
     epsa(i,j)=sqrt(epsi(i)*epsi(j)*sign6i*sign6i)
     print *,'epsi p9', sign6i, epsa(i,j), alf(i,j)
    endif

    !для букингема
    !bk1(i,j,k,l)=6.0*epsa(i,j,k,l)/(alf(i,j,k,l)-6.0)
    !bk2(i,j,k,l)=-epsa(i,j,k,l)*siga6/(1.0-6.0/alf(i,j,k,l))
    siga6=siga(i,j)**6
    kk1(i,j)=6.0*epsa(i,j)*siga6/alf(i,j)
    kk2(i,j)=-epsa(i,j)*(alf(i,j)+6.0)*siga6/alf(i,j)
    !_______________________
    kk3(i,j)=-6.0*epsa(i,j)*siga6/alf(i,j)
    kk4(i,j)=-epsa(i,j)*siga6
    !для потенциала с непонятным отталкиванием
    kkm1(i,j)=modif*epsa(i,j)*siga(i,j)**modif/alf(i,j)
    kkm2(i,j)=-epsa(i,j)*(alf(i,j)+modif)*siga(i,j)**modif/alf(i,j)
    !
    bk1(i,j)=6.0*epsa(i,j)/(alf(i,j)-6.0)
    bk2(i,j)=-epsa(i,j)*siga(i,j)**6/(1.0-6.0/alf(i,j))
    print *,siga(i,j), epsa(i,j),alf(i,j),kk1(i,j),kk2(i,j)
    !dbk1(i,j,k,l)=6.0*epsa(i,j,k,l)*siga6/(1.0-6.0/alf(i,j,k,l))
    !dbk2(i,j,k,l)=-6.0*epsa(i,j,k,l)/(siga(i,j,k,l)-6.0*siga(i,j,k,l)/alf(i,j,k,l))
   enddo
  enddo
 !Определение обрезания для потенциала букингема
if (ptip==3) then
 do i=1,Nm,1
  do j=1,Nm,1
   fmax=-1000         !обнуление максимума
   do zi=1,1000,1
    rz=zi*siga(i,j)/1000.0
    rz2=1.0/rz/rz
    rz6=rz2*rz2*rz2
    expp=exp(alf(i,j)*(1.0-rz/siga(i,j)))
    ekt=bk1(i,j)*expp+bk2(i,j)*rz6
    !if (ekt>0) then
     !print *, rz, ekt, fmax, rmax(i,j)
    ! pause
    !endif
    if (ekt>fmax) then
     rmax(i,j)=rz
     fmax=ekt
    endif
   enddo
    print *, rmax(i,j)
    !pause
  enddo
 enddo
endif
  !меняем температру
  razno=0.0 !обнуляем разницу
  razna=0.0
  write(estr,'(f10.2)') epsi(1) !к названию файла
  write(sstr,'(f10.3)') sigma(1) !
  write(astr,'(f10.2)') alfa(1)
  write(outf,'(2a)') 'out', trim(adjustl(Naz))      !trim(adjustl(workdir)),'/virK/',trim(adjustl(estr)),'v',trim(adjustl(sstr)),'v',trim(adjustl(astr))
  open(20,file=outf,position='append') !открываем файл на запись кривой
  maxra=0.0
  maxro=0.0
  do tk=1,kTemp,1
   if (kTemp==1) then
    Temp=Tempm
   else
    Temp=Tempm+(Tempk-Tempm)/float(kTemp-1)*float(tk-1)
   endif
   !print *,Temp,Tempm,Tempk
   schet=schet+1
   print *, schet, '  из  ', epsk*sigk*kTemp*alfk
   !print *, Temp
   !Определяем значение 2го вириального коэффициента
   !call set_sample  !сэмпилирование
   call MCinteg
   b2=b2*6.022136736/10.0
   !b2e=109.71-8.4673*10.0**4/Temp-8.1215*10.0**6/Temp/Temp-3.4382*10.0**9/Temp/Temp/Temp
   !b2e=87.776-60061.0/Temp+1.2857*float(10**6)/Temp/Temp-1.0861*float(10**9)/Temp/Temp/Temp
   !b2e=42.859-1.7696*10000.0/Temp+5.2007*100000.0/Temp/Temp-1.6393*100000000.0&
   !&/Temp/Temp/Temp+5.0855*1000000000.0/Temp/Temp/Temp/Temp
   !b2e=4.2859*10-1.7696*10**4/Temp+5.2007*10**5/Temp**2-1.6393*10**8/Temp**3+ &
   !&5.0855*10**9/Temp**4
   !b2e=107.73-8.2548*10**4/Temp+5.2387*10**6/Temp/Temp-1.9764*10**9/Temp/Temp/Temp
   !b2e=4.2859*10-1.7696*10**4/Temp+5.2007*10**5/Temp/Temp-1.6393*10**8/Temp/Temp/Temp+&
   !&5.0855*10**9/Temp/Temp/Temp/Temp
   !b2e=4.4344*10-1.6608*10**4/Temp-3.5430*10**6/Temp**2+2.9832*10**8&
   !&/Temp**3-23.448*10**9/Temp**4
   b2e= 1.0 !4.4697*100-4.7244*10**5/Temp+1.0310*10**8/Temp**2-23.475*10**9&
   !&/Temp**3  !-23.448*10**9/Temp**4
   abst=b2-b2e
   otnt=abs(b2-b2e)/b2e
   write(20,'(f10.4,a1,f18.10,a1,f18.10,a1,f18.10,a1,f18.10)') Temp,';',b2,';',b2e,';',abst,';',otnt
   !print *, Temp,';',b2,';',b2e
   razna=razna+abs(b2-b2e)
   if (maxra<abs(b2-b2e)) then
    maxra=abs(b2-b2e)
   endif
   razno=razno+abs(abs(b2-b2e)/b2e)
   if (maxro<abs(abs(b2-b2e)/b2e)) then
    maxro=abs(abs(b2-b2e)/b2e)
   endif
  enddo
  razna=razna/kTemp
  razno=razno/kTemp
   open(24,file='out.txt', position='append')
    write(24,'(f15.10,a1,f15.10,a1,f15.10,3(a1,f20.10),a1,f20.10)') epsi(1),';',sigma(1),';' &
    & ,alfa(1),';',razna,';',razno,';',maxra,';',maxro
   close(24)

 enddo
enddo
enddo
close(20)
!находим ошибку
end program

!расчет среднего методом монте-карло
subroutine MCinteg
use dannie
open(29, file='prov.txt')
b2=0.0 !обнуляем 2ой вириальный коэффициент
dr=0.002            !измеянем до ? было 0,02
Rcut2=0.01 !начальное значение растояния
do iri=1,3000       !было 3000
 dr=0.005           !было 0,005
 !print * ,siga(1,1)
 if (Rcut2>12.0) then
  dr=0.1
 endif
 Rcut2=Rcut2+dr
 rdr(iri)=dr
 b2dif2=0.0 !обнуляем счетчик на повороты
 en2=0.0
 do step=1,Nstep,1
!поворачиваем одну молекулу
  randalfa=getrand()
  randbeta=getrand()
  randgamma=getrand()
!Находим синусы и косинусы поворачиваемых углов
  sintx=(randalfa-0.5)*2
  costx=sqrt(1.0-sintx*sintx)
  sinty=(randbeta-0.5)*2
  costy=sqrt(1.0-sinty*sinty)
  sintz=(randgamma-0.5)*2
  costz=sqrt(1.0-sintz*sintz)
!Поворачиваем молекулу относительно оси x
  do i=1,Nm,1
   yv=ya1(i)
   zv=za1(i)
   ya1(i)= yv*costx+zv*sintx
   za1(i)= (-yv)*sintx+zv*costx
  enddo
!Поворачиваем молекулу относительно оси y
  do i=1,Nm,1
   xv=xa1(i)
   zv=za1(i)
   xa1(i)= xv*costy+zv*sinty
   za1(i)= (-xv)*sinty+zv*costy
  enddo
! Поворачиваем молекулу относительно оси z
  do i=1,Nm,1
  xv=xa1(i)
   yv=ya1(i)
   xa1(i)= xv*costz+yv*sintz
   ya1(i)= (-xv)*sintz+yv*costz
  enddo
  !
  !
  !поворачиваем другую молекулу
  randalfa=getrand()
  randbeta=getrand()
  randgamma=getrand()
 !Находим синусы и косинусы поворачиваемых углов
  sintx=(randalfa-0.5)*2
  costx=sqrt(1.0-sintx*sintx)
  sinty=(randbeta-0.5)*2
  costy=sqrt(1.0-sinty*sinty)
  sintz=(randgamma-0.5)*2
  costz=sqrt(1.0-sintz*sintz)
 !Поворачиваем молекулу относительно оси x
  do i=1,Nm,1
   yv=ya2(i)
   zv=za2(i)
   ya2(i)= yv*costx+zv*sintx
   za2(i)= (-yv)*sintx+zv*costx
  enddo
 !Поворачиваем молекулу относительно оси y
  do i=1,Nm,1
   xv=xa2(i)
   zv=za2(i)
   xa2(i)= xv*costy+zv*sinty
   za2(i)= (-xv)*sinty+zv*costy
  enddo
 !Поворачиваем молекулу относительно оси z
  do i=1,Nm,1
   xv=xa2(i)
   yv=ya2(i)
   xa2(i)= xv*costz+yv*sintz
   ya2(i)= (-xv)*sintz+yv*costz
  enddo
  b2dif=0.0 !обнуляем счетчтк на одну конфигурацию
  ekt=0.0
  do i=1,Nm,1
   do j=1,Nm,1
    dxa=xa2(j)-xa1(i)+Rcut2
    dya=ya2(j)-ya1(i)
    dza=za2(j)-za1(i)
    !sig2=siga(i,j)
    !sig2=sig2*sig2
    !rc=sig2/Rcut2m !в квадрате
    !rc=rc*rc*rc
    !выражение для вириального коэффициента
    Rcut2m=dxa*dxa+dya*dya+dza*dza
    if (Rcut2m==0.0) then
     ekt=9990000.0+ekt !оч большое число
    else
     if (ptip==1) then
      rz=sqrt(Rcut2m)
      if (rz<0.5) then
       ekt=ekt+999900000.0
      else
       rz2=1.0/Rcut2m
       rz6=rz2*rz2*rz2
!print *, siga(i,j)
       expp=exp(alf(i,j)*(1.0-rz/siga(i,j)))
       !ekt=ekt+(kk1(i,j)*expp+kk2(i,j)*expp+kk3(i,j)+kk4(i,j))*rz6
       dekt=(kk1(i,j)*expp+kk2(i,j))*rz6
       ekt=ekt+dekt
!print *,rz,Rcut2,dekt
!pause
!if (Rcut2>5) then
 !print *, xa2(j),ya2(j),za2(j)
 !print *, xa1(i),ya1(i),za1(i)
 !print *,ekt, dekt,rz, Rcut2
 !pause
!endif
      endif
     endif
     if (ptip==7) then
      rz=sqrt(Rcut2m)
      rz6=(1.0/rz)**modif
      expp=exp(alf(i,j)*(1.0-rz/siga(i,j)))
      !ekt=ekt+(kk1(i,j)*expp+kk2(i,j)*expp+kk3(i,j)+kk4(i,j))*rz6
      ekt=ekt+(kkm1(i,j)*expp+kkm2(i,j))*rz6
     endif
     if (ptip==9) then
      rz=sqrt(Rcut2m)
      ekt=ekt+epsa(i,j)*((siga(i,j)/rz)**alf(i,j)-(siga(i,j)/rz)**6)
     endif
     !----------------------------------------------------------
     if (ptip==2) then
      rz2=siga(i,j)*siga(i,j)/Rcut2m
      rz6=rz2*rz2*rz2
      ekt=4.0*epsa(i,j)*(rz6*rz6-rz6)+ekt
     endif
     !-----------------------------------------------------------
     if (ptip==3) then
      if (Rcut2m>rmax(i,j)*rmax(i,j)) then
       rz=sqrt(Rcut2m)
       rz2=1.0/Rcut2m
       rz6=rz2*rz2*rz2
       expp=exp(alf(i,j)*(1.0-rz/siga(i,j)))
       ekt=ekt+bk1(i,j)*expp+bk2(i,j)*rz6
       else
       ekt=ekt+99999999999.0
      endif
     endif
     if (ptip==6) then
      rz=sqrt(Rcut2m)
      if ((rz>1.8).and.(rz<2.7373)) then
       ekt=ekt-574.163
      else
       rz2=1.0/Rcut2m
       rz6=rz2*rz2*rz2
       expp=exp(alf(i,j)*(1.0-rz/siga(i,j)))
       !ekt=ekt+(kk1(i,j)*expp+kk2(i,j)*expp+kk3(i,j)+kk4(i,j))*rz6
       ekt=ekt+(kk1(i,j)*expp+kk2(i,j))*rz6
      endif
     endif
     if (ptip==8) then
        if (Rcut2m>0.005) then
            rz=sqrt(Rcut2m)
            rz2=siga(i,j)*siga(i,j)/Rcut2m
            rz6=rz2*rz2*rz2
            ekt=4.0*epsa(i,j)*(rz6*rz6-rz6)+ekt+167100.9566*iq(i)*iq(j)/rz
            if ((i==2).and.(j==2)) then
                !print *,rz, Rcut2, ekt, epsa(i,j),siga(i,j)
                !pause
            endif
        endif
     endif
     !----------------------------------------------------------
    endif
   enddo
  enddo
  en2=en2+ekt
  b2dif=b2dif+(exp(-ekt/Temp)-1.0)*Rcut2*Rcut2
  !print *,ekt, temp
  !if ((mod(step,10)==0).and.(Rcut2<30)) then
  ! write(29,'(f30.10,a,f30.10,a,f30.10 )') Rcut2, ';', ekt/float(Nstep),';',b2dif
  !endif
  !суммирование выражения
  b2dif2=b2dif2+b2dif
 enddo
b2int(iri)=(b2dif2/float(Nstep))*(-2.0)*PICH
  !open(
  !if (en2/float(Nstep)<5000) then
   write (29,'(f30.15,a,f30.15)'),rcut2,';',en2/float(Nstep)
  !endif
!print *,b2int(iri)
enddo
do k=1,2999
 b2=b2+(b2int(k)+b2int(k+1))/2.0*rdr(k)
 !print *,b2int(k),b2int(k+1),rdr(k)
enddo
close(29)
!print *,'virkof', b2
end subroutine

!Функция рандомирует случайные числа
!начало рандома
subroutine randomn
use generator
integer(4) i
real(8) sseed1,sseed2,sseed3
integer(4) seed(20)
character(20) a,b,c
real(8) n11,n21,n31
call date_and_time(a,b,c,seed)
sseed1=0.0
sseed2=0.0
sseed3=0.0
!print *,seed
do i=1,20,1 !здесь шарной бред
 sseed1=sseed1+seed(i)
 if (mod(i,2)==0) then
  sseed2=sseed2+seed(i)
 else
  sseed2=sseed2/2.0
 endif
 if (mod(i,3)==0) then
  sseed3=sseed3+seed(i)
 else
  sseed3=sseed3-seed(i)/2.2
 endif
enddo
!print *, sseed1,sseed2,sseed3
!pause
!проверка  тригонометрических функций
n11=abs(cos(sseed1/1000))
n21=abs(sin(sseed2/1000))
n31=abs(cos(sseed3/1000))
!call random_seed(int(sseed))
!call random_number(n11)
!call random_number(n21)
!call random_number(n31)
m1=2147483563
a1=40014
b1=12345
m2=2147483399
a2=40692
b2=54321
max3=2147483647
max1=m1-1
max2=m2-1
n1=ceiling((n11*max1)-1)
n2=ceiling((n21*max2)-1)
n3=ceiling((n31*max3)-1)
!   open(15,file='rand.txt')
 do i=1,500
 call randstart()
 randmass(i)=outrand
!       write(15,'(f30.20)') randmass(i)
 enddo
write (6,'(a10,3i15)') 'zatravka',n1,n2,n3
!   close(15)
end subroutine
!*******************************************************
!Рандомная функция начало
subroutine  randstart()
use generator
n1=abs(mod(a1*n1+b1,m1))
n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
if (float(n3)/float(max3)<0.5) then
 outrand=float(n1)/float(max1)
else
 outrand=float(n2)/float(max2)
endif
end subroutine
!*******************************************************
!Рандомная функция
function getrand()
use generator
integer(4) repick
!n1=abs(mod(a1*n1+b1,m1))
!n2=abs(mod(a2*n2+b2,m2))
n3=abs(n3*1664525+1013904223)
repick=ceiling((float(n3)/float(max3))*500)
getrand=randmass(repick)
call randstart()
randmass(repick)=outrand
!getrand=0.012
end function

subroutine set_sample()
use dannie
t_bind_en=0.0
!первоначальное положение
do i=1,nm,1
 xmc(i)=xm(i)
 ymc(i)=ym(i)
 zmc(i)=zm(i)
enddo
i_s=0
sum_exp=0.0
do while (i_s<=sample_kol)
!подсчет энергии до
 do j=1,nm,1
  xmc_do(j)=xmc(j)
  ymc_do(j)=ymc(j)
  zmc_do(j)=zmc(j)
 enddo
 mc_en=0.0
 call bind_energy()
 !call tors_enengy()
 mc_en=bind_en+tors_en
!поворот по валентным углам
 call bind_ch()
!поворот по торсионным углам
! call tors_ch()
!подсчет энергии после
 call bind_energy()
 !call tors_en()
 mc_en_posle=bind_en+tors_en
!принятие/непринятие перемещения
delta=mc_en_posle-mc_en
if (exp(-delta/Temp)<getrand()) then
 !не принимаем
 do j=1,nm,1
  xmc(j)=xmc_do(j)
  ymc(j)=ymc_do(j)
  zmc(j)=zmc_do(j)
 enddo
else
 i_s=i_s+1
 !принимаем
 do j=1,nm,1
  xa(i_s,j)=xmc(j)
  ya(i_s,j)=xmc(j)
  za(i_s,j)=xmc(j)
 enddo
 !вычисление вероятности
 exp_en(i_s)=exp(-bind_en/Temp)
 bind_m(i_s)=bind_u_mc
 !print*, i_s, bind_m(i_s),exp_en(i_s)
 !pause
 sum_exp=sum_exp+exp_en(i_s)
endif
enddo
do i=1,sample_kol,1
 xx(i)=exp_en(i)/sum_exp
enddo
open (18,file='out-bind')
do i=1,sample_kol,1
 write(18,'(f20.10,a1,f20.10,a1,f20.10)') bind_m(i),';',exp_en(i),';',xx(i)
enddo
close(18)
end subroutine

subroutine bind_energy()
use dannie
integer(4) ib
bind_en=0.0
do ib=1,bind_kol(1),1
!
  xv1=xmc(bind1k(1,ib,1))-xmc(bind_cent(1,ib))
  yv1=ymc(bind1k(1,ib,1))-ymc(bind_cent(1,ib))
  zv1=zmc(bind1k(1,ib,1))-zmc(bind_cent(1,ib))
  xv2=xmc(bind2k(1,ib,1))-xmc(bind_cent(1,ib))
  yv2=ymc(bind2k(1,ib,1))-ymc(bind_cent(1,ib))
  zv2=zmc(bind2k(1,ib,1))-zmc(bind_cent(1,ib))
  !print *, xv1,xv2,yv1,yv2,zv1,zv2
  r1_bind=(xv1*xv1+yv1*yv1+zv1*zv1)
  r2_bind=(xv2*xv2+yv2*yv2+zv2*zv2)
  cos_bind=(xv1*xv2+yv1*yv2+zv1*zv2)/sqrt(r1_bind*r2_bind)
  bind_u_mc=acos(cos_bind)
!
 bind_en=bind_en+bind_k(1,ib)*(bind_u_mc-bind_fi(1,ib))&
 &*(bind_u_mc-bind_fi(1,ib))
enddo
end subroutine

subroutine bind_ch()
use dannie
do i=1,bind_kol(1),1
 if (ymc(bind1k(1,i,1))==0.0) then
  sin1=0.0
  cos1=1.0
 else
  sin1=ymc(bind1k(1,i,1))/sqrt(ymc(bind1k(1,i,1))*ymc(bind1k(1,i,1))&
  &+zmc(bind1k(1,i,1))*zmc(bind1k(1,i,1)))
  cos1=zmc(bind1k(1,i,1))/sqrt(ymc(bind1k(1,i,1))*ymc(bind1k(1,i,1))+&
  &zmc(bind1k(1,i,1))*zmc(bind1k(1,i,1)))
 endif
 do jto=1,nm,1
  yv=-zmc(jto)*sin1+ymc(jto)*cos1
  zv=zmc(jto)*cos1+ymc(jto)*sin1
  ymc(jto)=yv
  zmc(jto)=zv
 enddo

 !после этого центр ориентирован на оси x (y=0,z=0)
 if (zmc(bind1k(1,i,1))==0.0) then
 sin2=0.0
 cos2=1.0
 else
  sin2=zmc(bind1k(1,i,1))/sqrt(xmc(bind1k(1,i,1))*xmc(bind1k(1,i,1))+&
  &zmc(bind1k(1,i,1))*zmc(bind1k(1,i,1)))
  cos2=xmc(bind1k(1,i,1))/sqrt(xmc(bind1k(1,i,1))*xmc(bind1k(1,i,1))+&
  &zmc(bind1k(1,i,1))*zmc(bind1k(1,i,1)))
 endif
 do jto=1,nm,1
  xv=xmc(jto)*cos2+zmc(jto)*sin2
  zv=-xmc(jto)*sin2+zmc(jto)*cos2
  xmc(jto)=xv
  zmc(jto)=zv
 enddo
!передвигаем другую сторону
! после этого (z=0)
 if (zmc(bind2k(1,i,1))==0) then
  sin3=0.0
  cos3=1.0
 else
  sin3=zmc(bind2k(1,i,1))/sqrt(ymc(bind2k(1,i,1))*ymc(bind2k(1,i,1))+&
  &zmc(bind2k(1,i,1))*zmc(bind2k(1,i,1)))
  cos3=ymc(bind2k(1,i,1))/sqrt(ymc(bind2k(1,i,1))*ymc(bind2k(1,i,1))+&
  &zmc(bind2k(1,i,1))*zmc(bind2k(1,i,1)))
 endif
 do jto=1,nm,1
  yv=ymc(jto)*cos3+zmc(jto)*sin3
  zv=-ymc(jto)*sin3+zmc(jto)*cos3
  ymc(jto)=yv
  zmc(jto)=zv
 enddo
! крутим часть молекулы относительно оси z
! bind_rnd=0.0175*(getrand()-0.5)*2.0
 !bind_u(Nbind,ib)=bind_u(Nbind,ib)+bind_rnd !измеряем потом
 !print *,bind_rnd

 sin4=(getrand()-0.5)*0.2
 cos4=sqrt(1.0-sin4*sin4)  !подумать как сделать
 do jto=1,bind1(1,i),1
  yv=ymc(bind1k(1,i,jto))*cos4+xmc(bind1k(1,i,jto))*sin4
  xv=-ymc(bind1k(1,i,jto))*sin4+xmc(bind1k(1,i,jto))*cos4
  xmc(bind1k(1,i,jto))=xv
  ymc(bind1k(1,i,jto))=yv
 enddo
 !крутим обратно
 do jto=1,Nm,1
  yv=ymc(jto)*cos3-zmc(jto)*sin3
  zv=ymc(jto)*sin3+zmc(jto)*cos3
  zmc(jto)=zv
  ymc(jto)=yv
 enddo
 !крутим еще обратно
 do jto=1,nm,1
  xv=xmc(jto)*cos2-zmc(jto)*sin2
  zv=xmc(jto)*sin2+zmc(jto)*cos2
  zmc(jto)=zv
  xmc(jto)=xv
 enddo
 do jto=1,nm,1
  yv=zmc(jto)*sin1+ymc(jto)*cos1
  zv=zmc(jto)*cos1-ymc(jto)*sin1
  ymc(jto)=yv
  zmc(jto)=zv
 enddo
enddo
end subroutine
