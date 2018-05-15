
	program IS

	IMPLICIT REAL*8 (A-H,O-Z), INTEGER(I-N)
	Real(8), ALLOCATABLE :: R(:,:), P(:,:),Pr(:,:),Pm(:,:)
	Real(8), ALLOCATABLE :: AuxFilho(:)
	Integer, ALLOCATABLE ::Mating(:), Mordem(:),Nspecies(:)
	character*30 smigra, sdis

	common /c1/ nparasitas,nlp
	common /c2/ rmu,nG
	OPEN(UNIT=60,FILE='entradaTM.in',STATUS='old')
	Read(60,*) ntime
	Read(60,*) nhostin, sdr
	Read(60,*) ncarry, noff
	Read(60,*) ngenoma, nG
	Read(60,*) rmu
	Read(60,*) dismin,dismax, deltadis
	Read(60,*) tmmin,tmmax, deltatm
	Read(60,*) nlag
	Close(60)	
	
    nlp=ngenoma+3	
	nparasitas=3
	ALLOCATE (R(nhostin,2),P(nparasitas,nlp),Pr(nparasitas,nlp),Pm(nparasitas,nlp))
	ALLOCATE (Mating(nparasitas),Mordem(nparasitas),Nspecies(nparasitas))
	ALLOCATE (AuxFilho(nlp))
	call init_random_seed()
	tmigra=tmmin
	
	OPEN(UNIT=50,FILE='EvolDistance'//TRIM(adjustl(sdis))//'.dat',STATUS='UNKNOWN')
	
  Dotm: Do while (tmigra.le.tmmax)
	dis=dismin
	Dodis: Do while (dis.le.dismax)	
	distance=dis*sdr
	Print*,'distancia=',distance,'migracao=',tmigra
	ngs=0
	nparasitas=3
	nhost=nhostin
	nprmax=nhost*ncarry
	Deallocate(R,P,Pr,Pm,Mating,Mordem,Nspecies)
	ALLOCATE (R(nhost,2),P(nparasitas,nlp),Pr(nparasitas,nlp),Pm(nparasitas,nlp))
	ALLOCATE (Mating(nparasitas),Mordem(nparasitas),Nspecies(nparasitas))
	R=0.0
	P=0.0
	Pr=0.0
	Pm=0.0
	Mating=0
	Mordem=0
	AuxFilho=0
	Do while (ngs.lt.Nint(ngenoma/2.0))
		call random_number(h1)
		nn=1+Int(ngenoma*h1)
		if(P(1,nn).eq.0)then
			P(1,nn)=1
			ngs=ngs+1
		endif
	Enddo	
	Do np=1, nparasitas
		P(np,:)=P(1,:)
	Enddo
	z=Real(ngs)/10.0
	 R(1,1)=z
	ncont=1
	i=1
	Do while(ncont.lt.nhost) 
		ncont=ncont+1
		R(ncont,1)=z+i*distance
		if(ncont.lt.nhost) then
			ncont=ncont+1
			R(ncont,1)=z-i*distance
		endif
		i=i+1
	Enddo
	R(1,2)=nparasitas
	P(:,nlp-2)=1.0
	P(:,nlp-1)=R(1,1)
	P(:,nlp)=Real(1)
!!!!!!
!!!!!!
	Kdis=Int(distance*100.0)
	Kmigra=NInt(tmigra*10000.0)
	write(sdis,*) Kdis
	write(smigra,*) Kmigra				
	

!!!!!
!!!!!
		Doo1:Do n=1, ntime	
			if(nparasitas.eq.0) go to 31
			Deallocate(Nspecies)
			Allocate(Nspecies(nparasitas))
			call species(P,Nspecies,nparasitas,nlp,ng,nespG)
			if(n.eq.1.or.mod(n,nlag).eq.0)then
				nhcol=0
				Do i=1, nparasitas
					nopt=Nint(P(i,nlp))
					write(50,16) distance, n ,P(i,nlp-1),P(i,nlp),R(nopt,1),R(nopt,2), Nspecies(i),Nint( P(i,nlp-2)) 
				Enddo
			endif	
			!Reproduction
			npacu=0
			nauxP=1
			nparasite=0
			npr=(noff+1)*nparasitas
			if(npr.gt.nprmax) npr=nprmax
			Deallocate(Pr)
			Allocate(Pr(npr,nlp))
			ncont=1
			Do nh=1, nhost
				ncontoneh=0
				npacu=npacu+nparasite
			   	nparasite= Nint(R(nh,2))
				if(nparasite.gt.0) then
					Deallocate(Mordem)
					Allocate(Mordem(nparasite))
					call shuffle(Mordem,nparasite,npacu)
					nupais=1
					ndes=nparasite*noff
					Do jj=1, ndes
						npai=Mordem(nupais)
						Deallocate(Mating)
						Allocate(Mating(nparasitas))
						call shuffle(Mating,nparasitas,0)
						call offspring(P,Mating,AuxFilho,npai,nparasite,nfilho)
						if(nfilho.eq.1) then
							Do jjj=1,nlp
								Pr(ncont,jjj)=AuxFilho(jjj)
							Enddo	
							nupais=nupais+1	
							if(nupais.gt.nparasite) nupais=1
							ncont=ncont+1
							ncontoneh=ncontoneh+1
							if(ncontoneh.ge.ncarry) then
							go to 30
							endif
						endif
					Enddo	
					30 continue
				endif
			Enddo
			!Survival and Probing
			npf=ncont-1
			Deallocate(Pm)
			Allocate(Pm(npf,nlp))	
			ncont=0
			Do i=1, npf
				call random_number(h1)
				if(h1.le.tmigra)then
					call random_number(h1)		
					novo=1+int(h1*nhost)		
					z=Pr(i,nlp-1)
					un=R(novo,1)
					f=exp(-((z-un)**2)/(2.0*sdr**2))
					call random_number(h)
					if(h.le.f)then
					    ncont=ncont+1
					    Pm(ncont,:)=Pr(i,:)
					    Pm(ncont,nlp)=Real(novo)	
					endif
				else
					z=Pr(i,nlp-1)
					nh=Nint(Pr(i,nlp))	
					u=R(nh,1)
					f=exp(-((z-u)**2)/(2.0*sdr**2))
					call random_number(h)
					if(h.le.f)then
						ncont=ncont+1
						Pm(ncont,:)=Pr(i,:)
					endif
				endif    
			Enddo
			nparasitas=ncont
			Deallocate(P)
			Allocate(P(nparasitas,nlp))
			!Sort and count
			np=0
			Do nh=1,nhost
				nhe=0
				Do i=1,nparasitas
					nuh=Nint(Pm(i,nlp))
					if (nuh.eq.nh) then
						np=np+1
						nhe=nhe+1
						P(np,:)=Pm(i,:)
					endif
				Enddo
				R(nh,2)=real(nhe)
			Enddo 
		Enddo Doo1
31 continue
	dis=dis+deltadis
	End do dodis
tmigra=tmigra+deltaTM
Enddo Dotm

Close(50)

15	FORMAT(100(1X,f4.2))
16	FORMAT(1x,f15.4,1X,I8,4(1x,f15.4),2(1X,I8))

	END PROGRAM IS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine species(V,Nv,np,nl,nG,nespG)
	IMPLICIT REAL*8 (A-H,O-Z)
  	  Dimension::Nv(np)
	Dimension::V(np,nl)
	Integer, ALLOCATABLE :: NES(:),Nes2(:),Nfreq(:,:)
	
	Nv=0
  	ngrupo=0
	nfim=0
	nGes=Nint(NG*1.0)
	nmax=Nint(V(1,nl-2))
	Do i=2,np
		nt=Nint(V(i,nl-2))
		if(nt.gt.nmax) nmax=nt
	enddo
	ALLOCATE(NES(nmax))

!!!!!
	Do while(nfim.eq.0)
		naux=1
		Do l=1,np
			if(Nv(l).eq.0)naux=0		
		Enddo
		nfim=naux
		Do l=1,np
			if(Nv(l).eq.0)then
				ngrupo=ngrupo+1
				Nv(l)=ngrupo
				go to 22		
			endif
		Enddo
		nespG=0
		Do nng=1,ngrupo
			nfreqG=0
			Do l=1,np
				if(Nv(l).eq.nng) nfreqG=nfreqG+1
			Enddo
			if(nfreqG.gt.10) nespG=nespG+1
		Enddo
		!!!
22      continue
		Do i=1,np
			  if(Nv(i).eq.ngrupo) then	
				Do ii=1,np
					if(Nv(ii).eq.0) then
					ndiff=0
					Do j=1,nl-3
					 n1=Nint(V(i,j))
					 n2=Nint(V(ii,j))
					 ndiff=ndiff+abs(n1-n2)		
					Enddo
					 if(ndiff.le.nGes)then
							Nv(ii)=ngrupo
							ncont=ncont+1
						endif
					endif
				enddo
			endif
		Enddo
	Enddo
	!!!!!
	ALLOCATE(Nfreq(ngrupo,2))
	Nes=0
	mesp=nmax
	Do n=1,nmax
		nfreq=0
		Do i=1,np
			mva=Nint(V(i,nl-2))
			if(mva.eq.n)then		
				Nfreq(nv(i),2)=Nfreq(nv(i),2)+1	
				Nfreq(nv(i),1)=nv(i)	
			endif
		Enddo
		call sortc(Nfreq,ngrupo)
		ngesco=Nfreq(1,1)
		if(Nes(n).eq.0)then
			Do i=1,np
				IF(Nv(i).eq.ngesco)then
					V(i,nl-2)=Float(n)
					Nes(n)=1					
				endif				
			Enddo
		else
			mesp=mesp+1
			Do i=1,np
				if(Nv(i).eq.Nfreq(1,1))then
					V(i,nl-2)=Float(mesp)
				endif
			enddo
		endif
		if(ngrupo.gt.1) then
			do nn=2,ngrupo
					if(Nfreq(nn,2).gt.10)then
						mesp=mesp+1
						Do i=1,np
							if(Nv(i).eq.Nfreq(nn,1))then
								V(i,nl-2)=Float(mesp)
							endif
						enddo
					 Endif
			Enddo
		endif
	enddo

	  
      endsubroutine species

		!!!!!!!!!!!!!!!!!
	subroutine sortc(Nio,nb)
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension:: Nio(nb,2)
	do i=1,nb-1
		ip=nb+1-i
		ipa=nb+1-i-1
		if(Nio(ip,2).gt.Nio(ipa,2)) then
			naux=Nio(ip,2)
			Nio(ip,2)=Nio(ipa,2)
			Nio(ipa,2)=naux
			naux=Nio(ip,1)
			Nio(ip,1)=Nio(ipa,1)
			Nio(ipa,1)=naux
		endif
	enddo	
	endsubroutine sortc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine shuffle (NV,Ni,nacu)
	  IMPLICIT REAL*8 (A-H,O-Z)
      Dimension::NV(Ni)
	  Do n=1,Ni
		NV(n)=nacu+n
	  Enddo
	  Do n=1, Ni
		call random_number(h1)
		call random_number(h2)
		i1=1+int(h1*Ni)
		i2=1+int(h2*Ni)
		vh1=NV(i1)
		NV(i1)=NV(i2)
		NV(i2)=vh1
	 Enddo	  	
     endsubroutine shuffle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine offspring(V,M,Aux,npai,nparasite,nachou)
	  IMPLICIT REAL*8 (A-H,O-Z)
	  common /c1/ nparasitas,nlp
	  common /c2/ rmu,nG
      Dimension::V(nparasitas,nlp)
	  Dimension::M(nparasite)
	  Dimension::Aux(nlp)
	  nachou=0
		Do jm=1,nparasite
			nmae=M(jm)
				kdif=0
				Do j=1,nlp-3
					kdif=kdif+abs(V(npai,j)-V(nmae,j))
				Enddo 
				if(kdif.le.nG) then
					nachou=1
					call random_number(h1)
					ncorte=1+h1*(nlp-3)
					
					Do j=1,nlp-3
						call random_number(h1)
						if(j.le.ncorte) then
							Aux(j)=V(nmae,j)
							if(h1.le.rmu) Aux(j)=abs(V(nmae,j)-1) 			
						else
							Aux(j)=V(npai,j)
							if(h1.le.rmu) Aux(j)=abs(V(npai,j)-1)
						endif
					Enddo 

					Aux(nlp-2)=V(nmae,nlp-2)
           				Aux(nlp-1)=REal(sum(Aux(1:nlp-3)))/10.0
					Aux(nlp)=V(npai,nlp)

				endif
			if(nachou.eq.1) go to 21
		enddo
21 continue
	  
	  	
      endsubroutine offspring


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine sort(za,H, NiH,nb,Nio)
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension::H(NiH), Nio(nb)
	Integer, ALLOCATABLE ::Naux(:)
	Allocate(Naux(nb))

	
	Naux=Nio
	Do i=2,nb
		iha=Nio(i-1)
		ihb=Nio(i)
		da=abs(H(iha)-za)
		db=abs(H(ihb)-za)
		if(db.lt.da)then
			Do j=1,i-1
				ihc=Nio(j)
				dc=abs(H(ihc)-za)
				if(db.lt.dc)then
					Naux(j)=Nio(i)
					Do k=j+1,i
						Naux(k)=Nio(k-1)
					Enddo
					go to 55
				endif
				
			Enddo
		endif
55		continue
		Nio=Naux

	Enddo


	endsubroutine sort


	!!!!!!!!!!!!!!!!!
	subroutine sortb(za,H, NiH,nb,Nio)
	IMPLICIT REAL*8 (A-H,O-Z)
	Dimension::H(NiH), Nio(nb)
	Integer, ALLOCATABLE ::Naux(:)
	Allocate(Naux(nb))
	if(nb.eq.4)then
		Naux=0
		Naux(1)=Nio(1)
		ih=Nio(1)
		d1=abs(H(ih)-za)
		do i=1,nb
			do j=2, nb
				ih=Nio(j)
				d2=abs(H(ih)-za)
				if(d2.lt.d1)then
					
				endif
			enddo
		enddo	
	endif
	endsubroutine sortb


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine newgeneration (VV,VVn,fits,pv,Ni,Nin,rmuv)
	IMPLICIT REAL*8 (A-H,O-Z)
	  
      Dimension::VV(Ni),VVn(Nin),pv(Ni)
	  VVn=0.0d0
		Do ii=1,Nin 
			call random_number(fitf)
			fitf=fitf*fits
			pant=0.0d0	
			Do i=1, Ni					
				if(fitf.gt.pant.and.fitf.le.pv(i)) then 
					ip=i
					pant=pv(i) 
				endif
			Enddo
20              continue	 
			call random_number(hp)
			hp=hp*3.0*rmuv
			f=exp(-(hp**2)/(2*rmuv**2))
			call random_number(ff)
			if(ff.le.f) then
				vr=hp
				else
					go to 20
			endif      					
			call random_number(hp)
			if(hp.le.0.5) vr=-vr
			VVn(ii)=VV(ip)+vr 
		Enddo		
		Ni=Nin
      endsubroutine newgeneration

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine poisson (NV,NiP,NiH)
	  IMPLICIT REAL*8 (A-H,O-Z)
      Dimension::NV(NiP)

	  NV=0
	  Do n=1, NiH
		call random_number(h)
		i=1+int(h*NiP)
		Nv(i)=Nv(i)+1
	  Enddo
	  	
      endsubroutine poisson


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	subroutine output(H,P,NiH,NeH,NiP,NeP,nst,Kep)
	IMPLICIT REAL*8 (A-H,O-Z)
	Real(8), ALLOCATABLE :: Dzr(:),Dza(:) 
	Dimension::H(NiH),P(NiP)

	w=0.02
	zhmin=H(1)
	zhmax=H(1)
	Do i=1, NiH
		if(H(i).gt.zhmax) zhmax=H(i)
		if(H(i).lt.zhmin) zhmin=H(i)
	Enddo
	njaH=2+int((zhmax-zhmin)/w)
	
	

	zpmax=P(1)
	zpmin=P(1)

	Do i=1, NiP
		if(P(i).gt.zpmax) zpmax=P(i)
		if(P(i).lt.zpmin) zpmin=P(i)
	Enddo
	njaP=2+int((zpmax-zpmin)/w)					
	ALLOCATE (Dzr(njaH),Dza(njaP))
	Dzr=0.0d0
	Do nn=1, njaH
		zw=zhmin+w/2.0+(nn-1)*w
		zwi=zw-w/2.0
		zws=zw+w/2.0
		Do i=1, NiH
			if(H(i).ge.zwi.and.H(i).lt.zws) then
				Dzr(nn)=Dzr(nn)+1.0d0 
			endif
		Enddo
		Write(63,11) nst,zw,Dzr(nn),Dfloat(NiH)
	Enddo
	Dza=0.0d0
	Do nn=1, njaP
		zw=zpmin+w/2.0+(nn-1)*w
		zwi=zw-w/2.0
		zws=zw+w/2.0
		Do i=1, NiP
			if(P(i).ge.zwi.and.P(i).lt.zws) then
				Dza(nn)=Dza(nn)+1.0d0 
			endif
		Enddo
		Write(65,11) nst,zw,Dza(nn),Dfloat(NiP)
	Enddo
11	FORMAT(1X,I8,100(1X,f15.8)) 
	endsubroutine output
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE init_random_seed()
        INTEGER :: l, k, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        CALL RANDOM_SEED(size = k)
        ALLOCATE(seed(k))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (l - 1, l = 1, k) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)
	END SUBROUTINE init_random_seed
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
