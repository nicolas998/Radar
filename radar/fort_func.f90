module radar_f90

!-----------------------------------------------------------------------
!Variableses globales del tipo de mapa leido
!-----------------------------------------------------------------------
!Globales
real xll,yll !coordenadas de la esquina inferior izquierda
real noData !Valor que representa los no datos
real dx !Largo del mapa leido
real dxP !largo proyectado del mapa, debe ser indicado si no se cononce
integer ncols,nrows !cantidad de columnas y filas del mapa
integer, allocatable :: ObjectsTemp(:,:) !Fila columna e id de los objetos que se encuentran en una imagen
integer, allocatable :: ObjectsNumTemp(:,:) !Lista la cantidad de elementos de cada objeto encontrado
real, allocatable :: LenghtTemp(:) !Array con las longitudes medidas
real, allocatable :: MatrizRadios(:,:) !Matriz Con radios calculados (usado mucho en Steiner 1995)
!funcion para hacer sort
public :: QsortC

!-----------------------------------------------------------------------
!Punto de inicio de funciones del modulo
contains

!-----------------------------------------------------------------------
!Subrutinas para operar con radares como si fueran imagenes
!-----------------------------------------------------------------------

subroutine detect_clouds(image,Grad,kerX,kerY,nc,nf) !detecta las nubes en un dem
	!Variables de entrada
	integer, intent(in) :: nc,nf
	real, intent(in) :: KerX(3,3),KerY(3,3)
	real, intent(in) :: image(nc,nf)
	!variables de salida 
	real, intent(out) :: Grad(nc,nf)
	!f2py intent(in) :: nc,nf,image
	!f2py intent(out) :: Grad
	!Variables locales
	integer i,j
	real DifX(nc,nf),DifY(nc,nf),imageTemp(nc,nf),sumaX,sumaY
	logical mascara(nc,nf)
	!itera por todo el mapa calculando bordes
	imageTemp=image
	where(imageTemp .eq. nodata) imageTemp=0.0
	mascara=.false.
	where(imageTemp .eq. nodata) mascara=.true.
	DifX=0; DifY=0
	do i=2,nc-1
		do j=2,nf-1
			sumaX=0
			sumaY=0
			do ki=1,3
				do kj=1,3
					sumaX=sumaX+imageTemp(i-2+ki,j-2+kj)*KerX(ki,kj)
					sumaY=sumaY+imageTemp(i-2+ki,j-2+kj)*KerY(ki,kj)
				enddo
			enddo
			if (any(mascara(i-1:i+1,j-1:j+1))) then				
				DifX(i,j)=0.0
				DifY(i,j)=0.0
			else	
				DifX(i,j)=sumaX
				DifY(i,j)=sumaY
			endif
		enddo
	enddo
	!Calcula el gradiente de diferencias 
	Grad=sqrt(DifX**2+DifY**2)
end subroutine

subroutine erosion(imageIn,imageOut,kernel,N,nc,nf) !Hace una erosion sobre el elemento
	!Variables de entrada
	integer, intent(in) :: N,nc,nf,kernel
	integer, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf)
	!f2py intent(in) :: N,nc,nf,kernel,imageIn
	!f2py intent(out) :: imageOut
	!Variables locales 
	integer i,j,t,k,imageTemp(nc,nf)
	!Itera La cantidad de veces especificada
	k=kernel/2.0
	imageTemp=imageIn
	imageOut=imageIn
	do t=1,N
		!Itera para toda la imagen binaria
		do i=2,nc-1
			do j=2,nf-1
				!se fija en el kernel
				if (any(imageTemp(i-k:i+k,j-k:j+k) .eq. 0)) then
					imageOut(i,j)=0
				endif
			enddo
		enddo
		imageTemp=imageOut
	enddo
end subroutine

subroutine dilation(imageIn,imageOut,kernel,N,nc,nf) !Hace una erosion sobre el elemento
	!Variables de entrada
	integer, intent(in) :: N,nc,nf,kernel
	integer, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf)
	!f2py intent(in) :: N,nc,nf,kernel,imageIn
	!f2py intent(out) :: imageOut
	!Variables locales 
	integer i,j,t,k,imageTemp(nc,nf)
	!Itera La cantidad de veces especificada
	k=kernel/2.0
	imageTemp=imageIn
	imageOut=imageIn
	do t=1,N
		!Itera para toda la imagen binaria
		do i=2,nc-1
			do j=2,nf-1
				!se fija en el kernel
				if (any(imageTemp(i-k:i+k,j-k:j+k) .eq. 1)) then
					imageOut(i,j)=1
				endif
			enddo
		enddo
		imageTemp=imageOut
	enddo
end subroutine

subroutine opening(imageIn,imageOut,kernel,nc,nf)
	!Variables de entrada
	integer, intent(in) :: nc,nf,kernel
	integer, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf)
	!f2py intent(in) :: nc,nf,kernel,imageIn
	!f2py intent(out) :: imageOut
	!Variables locales 
	integer i,j,t,k,imageTemp(nc,nf)
	!Itera La cantidad de veces especificada
	k=kernel/2.0
	!Itera para hacerlo N veces
	call erosion(imageIn,imageTemp,kernel,1,nc,nf)
	call dilation(imageTemp,imageOut,kernel,1,nc,nf)
end subroutine

subroutine closing(imageIn,imageOut,kernel,nc,nf)
	!Variables de entrada
	integer, intent(in) :: nc,nf,kernel
	integer, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf)
	!f2py intent(in) :: nc,nf,kernel,imageIn
	!f2py intent(out) :: imageOut
	!Variables locales 
	integer i,j,t,k,imageTemp(nc,nf)
	!Itera La cantidad de veces especificada
	k=kernel/2.0
	!Itera para hacerlo N veces
	call dilation(imageIn,imageTemp,kernel,1,nc,nf)
	call erosion(imageTemp,imageOut,kernel,1,nc,nf)
end subroutine

subroutine classify_binary(imageIn,imageOut,nc,nf,Nelem,Npixels)
	!Variables de entrada
	integer, intent(in) :: nc,nf
	integer, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf),Nelem,Npixels
	!f2py intent(in) :: nc,nf,imageIn
	!f2py intent(out) :: imageOut,Nelem,Npixels
	!Variables locales
	integer i,j,cont2,enc,lista(2,nc*nf),ki,kj,x,y
	logical flag
	!Maneja la memoria de la variable temporal de objetos 
	if (allocated(ObjectsTemp)) deallocate(ObjectsTemp)
	allocate(ObjectsTemp(3,nc*nf))	
	if (allocated(ObjectsNumTemp)) deallocate(ObjectsNumTemp)
	allocate(ObjectsNumTemp(1,nc*nf))	
	!Itera sobre la imagen buscando diferentes objetos 
	Nelem=0
	Npixels=0
	lista=0
	do i=1,nc
		do j=1,nf
			!Si encuentra una celda que es uno busca su vecindad para 
			!hallar todos los que pertenecen al evento
			if ((imageIn(i,j).eq.1).and.(imageOut(i,j).eq.0)) then
				Nelem=Nelem+1; Npixels=Npixels+1						
				ObjectsTemp(1,Npixels)=Nelem
				ObjectsTemp(2,Npixels)=i
				ObjectsTemp(3,Npixels)=j
				enc=1; cont2=1
				lista=-999
				lista(1,1)=i
				lista(2,1)=j
				imageOut(i,j)=Nelem				
				!Itera hasta no encontrar mas
				do while (lista(1,cont2) .ne. -999)					
					x=lista(1,cont2)
					y=lista(2,cont2)
					!Busca en la vecindad de los 8
					do ki=1,3
						do kj=1,3
							!Si es 1 lo agrega a la lista
							if ((imageIn(x+ki-2,y+kj-2) .eq. 1) .and. &
							(imageOut(x+ki-2,y+kj-2) .eq. 0) ) then								
								Npixels=Npixels+1
								enc=enc+1								
								lista(1,enc)=x+ki-2
								lista(2,enc)=y+kj-2
								imageOut(x+ki-2,y+kj-2)=Nelem	
								ObjectsTemp(1,Npixels)=Nelem
								ObjectsTemp(2,Npixels)=x+ki-2
								ObjectsTemp(3,Npixels)=y+kj-2								
							endif
						enddo
					enddo
					!actualiza contador interno
					cont2=cont2+1
				enddo
				ObjectsNumTemp(1,Nelem)=enc	
			endif	
		enddo
	enddo
end subroutine

subroutine cut_List_Object(ObjectList,ObjectNum,Npixels,Nelem)
	!varaibles de entrada
	integer, intent(in) :: Npixels,Nelem
	!Variabnles de salida
	integer, intent(out) :: ObjectList(3,Npixels),ObjectNum(1,Nelem)
	!f2py intent(in) :: Npixels,Nelem
	!f2py intent(out) :: ObjectList,ObjectNum
	!Corta 
	if (allocated(ObjectsTemp)) then
		ObjectList(:,:)=ObjectsTemp(:,1:Npixels)
		deallocate(ObjectsTemp)
	else
		print *, 'Error: variable ObjectsTemp no alojada en memoria'
	endif
	!Corta 
	if (allocated(ObjectsNumTemp)) then
		ObjectNum(1,:)=ObjectsNumTemp(1,1:Nelem)
		deallocate(ObjectsNumTemp)
	else
		print *, 'Error: variable ObjectsNumTemp no alojada en memoria'
	endif
end subroutine

subroutine clean_by_size(imageIn,imageOut,nc,nf,ObjectList,ObjectNum,Npixels,Nelem,umbral)
	!Variables de entrada
	integer, intent(in) :: nc,nf,umbral,Npixels,Nelem
	integer, intent(in) :: imageIn(nc,nf),ObjectList(3,Npixels)
	integer, intent(in) :: ObjectNum(Nelem)
	!Variables de salida
	integer, intent(out) :: imageOut(nc,nf)
	!f2py intent(in) :: imageIn,nc,nf,ObjectList,ObjectNum,Npixels,Nelem,umbral
	!f2py intent(out) :: imageOut
	!Variables locales 
	integer i,j
	!Si hay menos elementos de los presentados en el umbral
	!los elimina
	imageOut=imageIn
	do i=1,Nelem
		if (ObjectNum(i).lt.umbral) then
			where(imageOut(:,:).eq.i) imageOut=0
		endif	
	enddo
	where(imageOut.ne.0) imageOut=imageOut/imageOut
end subroutine

subroutine objects_lenght(ObjectList,ObjectNum,DistLenght,MaxLenght,nc,nf,&
	&Npixels,Nelem)
	!Variables de entrada
	integer, intent(in) :: nc,nf,Npixels,Nelem
	integer, intent(in) :: ObjectList(3,Npixels),ObjectNum(Nelem)
	!Variables de salida
	real, intent(out) :: MaxLenght(Nelem),DistLenght(19,Nelem)
	!f2py intent(in) :: nc,nf,Nelem,Npixels,ObjectList,ObjectNum
	!f2py intent(out) :: MaxLenght,DistLenght
	!Variables locales 
	integer elem,p,pos1,pos2,celda1,celda2,cont,cont2,tamano
	real X1,X2,Y1,Y2,lenght,lenghtNueva
	!itera para cada elenmento encontrado dentro de la matriz 
	pos1=1
	cont=0
	do elem=1,Nelem		
		!Guarda las longitudes del elemento
		if (allocated(LenghtTemp)) deallocate(LenghtTemp)
		tamano=cum_sum(ObjectNum(elem))
		allocate(LenghtTemp(tamano))		
		LenghtTemp=0.0
		!Itera para todas las celdas del elemento
		pos2=ObjectNum(elem)+cont
		cont2=1
		lenght=0.0
		do celda1=pos1,pos2-1
			!calcula la posicion X,Y relativa en metros de la celda base			
			X1=ObjectList(2,celda1)*dxp
			Y1=ObjectList(3,celda1)*dxp
			!Itera para comparar con las demas celdas y va guardado en LenghtTemp			
			do celda2=celda1+1,pos2
				X2=ObjectList(2,celda2)*dxp
				Y2=ObjectList(3,celda2)*dxp
				lenght=sqrt((X1-X2)**2+(Y1-Y2)**2)									
				LenghtTemp(cont2)=lenght/1000.0
				cont2=cont2+1
			enddo			
		enddo
		!Guarda la mas larga
		call QsortC(LenghtTemp)	
		MaxLenght(elem)=LenghtTemp(tamano)
		do p=1,19					
			DistLenght(p,elem)=LenghtTemp(ceiling((0.0+p*0.05)*tamano))
		enddo
		!Actualiza contadores
		pos1=ObjectNum(elem)+1+cont
		cont=cont+ObjectNum(elem)
	enddo
end subroutine

subroutine arc_slope(imageIn,slope,nc,nf)
	!Variables de entrada
	integer, intent(in) :: nc,nf
	real, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	real, intent(out) :: slope(nc,nf)
	!f2py intent(in) :: nc,nf,imageIn
	!f2py intent(out) :: slope
	!Variables locales 
	integer i,j
	real imageCopy(nc+2,nf+2)
	real pendX,pendY,k(3,3)
	imageCopy(2:nc,2:nf)=imageIn
	do i=2,nc+1
		do j=2,nf+1
			!Si el punto no es cero caklcula la pendiente 
			if (imageCopy(i,j).ne.0) then
				k=imageCopy(i-1:i+1,j-1:j+1)
				pendX=((k(3,1)+2*k(3,2)+k(3,3))-(k(1,1)+2*k(1,2)+k(1,3)))/(8*dxp)						
				pendY=((k(1,3)+2*k(2,3)+k(3,3))-(k(1,1)+2*k(2,1)+k(3,1)))/(8*dxp)
				slope(i-1,j-1)=sqrt(pendX**2 + pendY**2)	
			endif
		enddo
	enddo
end subroutine

subroutine var2mean(imageIn,ObjectList,meanVar,stdVar,nc,nf,Nelem,Npixels)
	!Variables de entrada 
	integer, intent(in) :: nc, nf, Nelem,Npixels
	integer, intent(in) :: ObjectList(3,Npixels)
	real, intent(in) :: imageIn(nc,nf)
	!Variables de salida
	real, intent(out) :: meanVar(Nelem),stdVar(Nelem)
	!f2py intent(in) :: nc,nf,Nelem,Npixels,ObjectList,imageIn
	!f2py intent(out) :: meanVar
	!Variables locales 	
	integer i
	real Values(Npixels),sumValue,contValue
	!Obtiene los valores de las nubes en todas las celdas de los elementos 
	do i=1,Npixels
		Values(i)=imageIn(ObjectList(2,i),ObjectList(3,i))
	enddo
	!Itera para todos los elementos del objeto 
	do i=1,Nelem
		contValue=count(ObjectList.eq.i)
		sumValue=sum(Values,mask=ObjectList(1,:).eq.i)
		meanVar(i)=sumValue/contValue
		stdVar(i)=sum((Values-meanVar(i))**2,mask=ObjectList(1,:).eq.i)/(contValue)
		stdVar(i)=sqrt(stdVar(i))
	enddo		
end subroutine

!-----------------------------------------------------------------------
!Subrutinas Para calcular dimension fractal
!-----------------------------------------------------------------------

subroutine fractal3d(imageIn,ObjectList,ker,nc,nf,Npixels,a,Fractal)
	!Variables de entrada
	integer, intent(in) :: nc,nf,Npixels,ker,a
	integer, intent(in) :: ObjectList(3,Npixels)
	real, intent(in) :: imageIn(nc,nf)
	!Variables de salida 
	real, intent(out) :: Fractal(nc,nf)
	!f2py intent(in) :: nc,nf,Npixels, ObjectList, imageIn,a
	!f2py intent(out) :: Fractal
	!Variables locales 
	integer i,j,ncn,nfn,halfKer,col,fil
	real MatTemp(nc+ker-1,nf+ker-1)
	real MatKernel(ker,ker)
	!Copia la imagen de reflectividad en la temporal 
	ncn=nc+ker-1; nfn=nf+ker-1
	halfKer=ker/2
	MatTemp(halfKer+1:ncn-halfKer+1,halfKer+1:nfn-halfKer+1)=imageIn(:,:)
	!Itera sobre los pixeles de las nubes encontradas
	Fractal=0.0
	do i=1,Npixels
		!toma la matriz de datos alrededor del pixel seleccionado 
		col=ObjectList(2,i)+halfKer
		fil=ObjectList(3,i)+halfKer
		MatKernel=MatTemp(col-halfKer:col+halfKer-1,fil-halfKer:fil+halfKer-1)
		Fractal(ObjectList(2,i),ObjectList(3,i))=fd(MatKernel,ker,a)
	enddo
end subroutine 

real function fd(Mat,k,a)
	!Defincion de variables
	integer k,a
	real Mat(k,k)
	!Variables locales 
	real NpFin(k), SFin(k)	
	real rp
	integer cont,nr,np,z,j
	!Encuentra la cantidad de multiplos del kernel 	y los multiplos
	cont=0
	do i=1,k
		if (mod(k,i).eq.0) then 
			!contabiliza los que son multiplos
			cont=cont+1			
			!Obtiene el tamano vertical de las celdas 
			rp=i/(1+2*a*std(Mat,k)) 
			!Obtiene la cantidad de celdas verticales 
			if (maxval(Mat) .eq. minval(Mat)) then 
				nr=1.0
			else
				nr=ceiling((maxval(Mat)-minval(Mat))/rp)
			endif
			!Itera la ventana multiplo para obtener la cantidad 
			!de cuadros cubiertos por esta 			
			np=0
			do j=1,k,i
				do z=1,k,i
					np=np+ceiling((maxval(Mat(j:j+i-1,z:z+i-1))-minval(Mat(:,:)))/rp+1)
				enddo
			enddo
			!guarda la cantidad de cuadritos almacenados para esa dimension
			NpFin(cont)=np
			SFin(cont)=i
		endif
	enddo
	!Saca el logaritmo, hace la regresion y con eso el FD
	NpFin(1:cont)=log(NpFin(1:cont))
	SFin(1:cont)=log(SFin(1:cont))
	fd=-1*slope(SFin(1:cont),NpFin(1:cont),cont)
	if (fd .lt. 2.0) fd=2.0	
end function

real function std(Mat,k)
	!Variables de entrada
	integer k
	real Mat(k,k)
	!Variables locales 
	integer i
	real media,sumaDif 
	!calcula el valor medio 
	media=sum(Mat)/(k**2)	
	!Calcula la desviacion 
	sumaDif=0
	do i=1,k
		do j=1,k
			sumaDif=sumaDif+((Mat(i,j)-media)**2)
		enddo
	enddo
	std=sqrt(sumaDif/(k**2))
end function

real function slope(x,y,n)
	!Variables de entrada
	integer n
	real x(n),y(n)	
	!calcula 
	slope=(sum(x*y)-(sum(x)*sum(y))/n)/(sum(x**2)-((sum(x)**2)/n))	
end function

recursive function cum_sum(n) result(res)    
    integer res,n
    if (n .EQ. 0) then
        res = 1
    else
        res = n + cum_sum(n - 1)
    endif
end


!-----------------------------------------------------------------------
!Subrutinas Para clasificar de acuerdo a Steiner 1995
!-----------------------------------------------------------------------
subroutine steiner_genera_radios(radio) !Genera una matriz de radios a partir de un radio dado
	!Variables de entrada
	real, intent(in) :: radio !Ingresa en metros
	!Variables locales 
	integer Tamano, i,j, TamMedio
	!calcula el tamano de la matriz de radios 
	Tamano = ceiling(radio / dxp)+7
	TamMedio = ceiling(Tamano/2.0)
	if (allocated(MatrizRadios)) deallocate(MatrizRadios)
	allocate(MatrizRadios(Tamano,Tamano))
	!Calcula los radios para la matriz 
	MatrizRadios = noData
	do i=1,Tamano
		do j = 1,Tamano
			MatrizRadios(i,j) = sqrt(dxp*(TamMedio-i)**2.0+dxp*(TamMedio-j)**2.0)
		enddo
	enddo
	!mata los radios que superan el umbral 
	where(MatrizRadios .gt. radio) MatrizRadios = noData
end subroutine 

real function steiner1995(MZbg)
	!Variables de entrada
	real MZbg
	!Calcula
	if (MZbg .lt. 0) then 
		steiner1995 = 10
	elseif (MZbg .ge. 0 .and. MZbg .lt. 42.43) then 
		steiner1995 = 10 - (MZbg**2.0)/180.0
	elseif (MZbg .gt. 42.43) then 
		steiner1995 = 0
	endif
end function

real function siriluk2008(MZbg, Zbg, tam)
	!Variables de entrada
	integer tam
	real Zbg(tam,tam), MZbg
	!Variables locales 
	real Zc, P
	!Calcula 
	Zc = minval(Zbg, mask = Zbg .ne. noData)
	P = max(((Zc+2.5)**2.0)/10.0, 140)
	siriluk2008 = 10 - (MZbg**2)/P
end function

subroutine steiner_find_peaks(Ref, ncol, nfil, umbral, radio, metodo, peaks) !Esta funcion encuentra los picos para la clasificacion de conv y estratiforme de Steiner 1995
	!Variables de entrada
	integer, intent(in) :: ncol, nfil, metodo
	real, intent(in) :: Ref(ncol,nfil)
	real, intent(in) :: umbral, radio
	!Variables de salida
	integer, intent(out) :: peaks(ncol,nfil)
	!Variables locales 
	integer i,j
	real, allocatable :: Zbg(:,:)
	real, allocatable :: refLoc(:,:)
	real SumZbg, CountZbg, MeanZbg, DiffZbg
	!inicia los picos en cero 
	peaks = 0
	!Calcula cantidad de celdas e Inicia matriz de fondo para un radio de 11km
	Tamano = ceiling(radio / dxp)+7
	if (mod(Tamano,2) .eq. 0) Tamano = Tamano+1
	TamMedio = floor(Tamano/2.0)
	allocate(Zbg(Tamano,Tamano))
	allocate(refLoc(ncol+Tamano, nfil+Tamano))
	refLoc = nodata
	refLoc(TamMedio:ncol, TamMedio:nfil) = ref
	call steiner_genera_radios(radio)
	! Itera por toda la matriz de reflectividad
	do i=1,ncol
		do j=1,nfil
			!Evalua solo si tiene reflectividad y no es nulo 
			if ((refLoc(i,j) .ne. noData ) .and. (refLoc(i,j) .gt. 0)) then 
				!Evalua por criterio por Intensidad
				if (refLoc(i,j) .gt. umbral) then
					peaks(i,j) = 1
					call steiner_convective_radius()					
				!Si no encuentra picos por ese criterio busca por el criterio de peakness
				else
					!Calcula el Zbg Medio para el radio seleccionado y la diferencia con el punto central
					Zbg = refLoc(i-TamMedio:i+TamMedio, j-TamMedio:j+TamMedio)
					SumZbg = sum(Zbg, mask = MatrizRadios .ne. noData)
					CountZbg = count(MatrizRadios .ne. noData)
					MeanZbg = SumZbg / CountZbg
					DiffZbg = refLoc(i,j) - MeanZbg
					!Calcula el DifZcrit de acuerdo a la metodologia seleccionada
					if (metodo .eq. 1) then
						DiffZcrit = steiner1995(MeanZbg)
					elseif (metodo .eq. 2) then 
						DiffZcrit = siriluk2008(Zbg, Tamano)
					endif
					!Determina si el punto es o no por peakness
					if (DiffZbg .gt. DiffZcrit) then 
						peaks(i,j) = 2
						call steiner_convective_radius()
					endif
				endif
			endif
		enddo
	enddo
	
end subroutine

!-----------------------------------------------------------------------
!Funciones para hacer sort de algo
!Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
!Based on algorithm from Cormen et al., Introduction to Algorithms,
!-----------------------------------------------------------------------
recursive subroutine QsortC(A)
  real, intent(in out), dimension(:) :: A
  integer :: iq
    !f2py intent(inout) :: A
  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine 
subroutine Partition(A, marker) !subrutina utilizada por QsortC para hacer sort
  real, intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1
  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine

end module
