c		Aaron's nearest neighbour finder
       subroutine fof(x1,n1,nb,nz,ind,bgal,zcut,rd,zd,Nmax, iflag)

       integer n1,ind(Nmax,2),nb,nz
       integer kk,i,j,ic
       
       real*8 x1(n1,4),  bgal(nb), bgal2, zcut(nz), rd(Nmax)
	   real*8 zd(Nmax), radproj, ztemp, zrad

c****   counter  for accumulating close points
        kk=0 
		
c
c If bgal and rscale don't change, set them here
c		
		if(nb.eq.1) then 
		bgal2=bgal(1)**2
		endif
		if(nz.eq.1) then
		ztemp=zcut(1)
		endif
c
c Two for loops for all objects
c
            do  15 i= 1, n1
                  do 10 j =1,n1
				  if(i.eq.j) goto 10
				  
c
c Radial Z scale cut
c
		
			if(nz.ne.1) then
			ztemp=(zcut(i)/2)+(zcut(j)/2)
			endif
			
			zrad = ABS(x1(i,4) - x1(j,4))
			if (zrad.gt.ztemp) goto 10
c
c Projected radius chop
c
			if(nb.ne.1) then
			bgal2 = ((bgal(i)/2) + (bgal(j)/2))**2
			endif
			radproj=0.0
			do 5 ic= 1, 3
                    radproj= radproj + (x1(i,ic) - x1(j,ic))**2
                    if( radproj.gt.bgal2) goto 10
 5             continue

             kk=kk+1

c**** check if there is still space 
             if( kk .gt. Nmax) then 
                iflag= -1
                goto 20 
             else
                ind(kk,1)= i
                ind(kk,2)= j
                rd(kk)= sqrt(radproj)
				zd(kk)= zrad
             endif
       

 10             continue
 15         continue
 
         Nmax=kk  
20       continue


       return
       end