from radar_f90 import *
import numpy as np
import glob
try:
	from mpl_toolkits.basemap import Basemap, addcyclic, shiftgrid, cm
except:
	print 'se desactivan funciones de ploteo con Basemap'
	pass
import pylab as pl
import scipy.ndimage as nd
from netCDF4 import Dataset  
import pickle

Path = __file__
Path = Path[:-9]+'ajuste_multicapaall_77.pkl'
RadProp = []

def __open_pklfiles__(path_pkldata):
    open_pkl = open(path_pkldata, 'rb')
    data = pickle.load(open_pkl)
    open_pkl.close()
    return data

class image_process:
	def __init__(self):
		pass
	#erosiona imagen
	def erosion(self,binaria,kernel=3,N=1):
		'\n'\
		'Descripcion: erosiona la imagen, de acuerdo al tamano del\n'\
		'	kernel, si una celda dentro del kernel movil es 0, la \n'\
		'	celda central se hace igual a cero.  La funcion puede \n'\
		'	operar N veces sobre la imagen binaria \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'binaria : Imagen binaria del radar.\n'\
		'k : Tamano de la ventana para realizar la erosion.\n'\
		'N : Cantidad de veces que se realiza la erosion sobre la imagen.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'eroded : Imagen erosionada.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'eroded=erosion(binImage,5,3).\n'\
		
		return radar_f90.erosion(binaria,kernel,N,radar_f90.ncols,radar_f90.nrows)
	#Dilata la imagen
	def dilation(self,binaria,kernel=3,N=1):
		'\n'\
		'Descripcion: dilata la imagen binaria que sea entregada\n'\
		'	usa el kernel para dilatar, de acuerdo a su tamano \n'\
		'	busca si cualquier celda dentro del kernel es 1, \n'\
		'	las demas celdas se hacen igual a 1, itera N veces \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'binaria : Imagen binaria del radar.\n'\
		'k : Tamano de la ventana para realizar la dilatacion.\n'\
		'N : Cantidad de veces que se realiza la erosion sobre la imagen.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'dilation : Imagen dilatada.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'dilated=dilation(binImage,5,3).\n'\
		
		return radar_f90.dilation(binaria,kernel,N,radar_f90.ncols,radar_f90.nrows)
	#Abre la imagen	
	def opening(self,binaria,kernel=3):
		'\n'\
		'Descripcion: Abre la imagen en terminos de procesamiento\n'\
		'	de imagenes, par alo cual aplica primero una erosion\n'\
		'	seguida de una dilatacion\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'binaria : Imagen binaria del radar.\n'\
		'k : Tamano de la ventana para realizar la erosion.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'opened : Imagen abierta.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'open=openning(binImage,5).\n'\
		
		return self.erosion(self.dilation(binaria,kernel))
	#cierra la imagen
	def closing(self,imageIn,kernel=3):
		'\n'\
		'Descripcion: Cierra la imagen en terminos de procesamiento\n'\
		'	de imagenes, para lo cual aplica primero una dilatacion\n'\
		'	seguida de una erosion\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'binaria : Imagen binaria del radar.\n'\
		'k : Tamano de la ventana para realizar la erosion.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'closed : Imagen cerrada.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'closed=closing(binImage,5).\n'\
		
		return self.dilation(self.erosion(imageIn,kernel))
	#Obtiene bordes de la imagen
	def borders(self,imageIn):
		'\n'\
		'Descripcion: Encuentra los bordes de una imagen binaria\n'\
		'	para lo cual hace lo siguiente: imageIn-erosion(imageIn)\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen binaria del radar.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'borders : Bordes bien definidos de la imagen binaria.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'border=borders(binImage).\n'\
		
		return imageIn-self.erosion(imageIn)

class radar_process:	
	#------------------------------------------------------
	# Subrutinas de lectura y prerpocesamiento inicial
	#------------------------------------------------------
	#Inicia la clase
	def __init__(self,**kwargs):
		'Descripcion: Inicia las variables de radar con las que\n'\
		'	se trabaja: Z, ref, binario, objetos\n'\
		'funcion asigna parametros al modulo radar \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : Inicia las variables vacias.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self : Con las variables iniciadas.\n'\
		
		#Inicia el objeto 
		self.Z=None
		self.ref=None
		self.binario=None
		self.borders=None
		self.objetos=None
		self.image=image_process()
		self.plot=draw_func()					
		#Inicia las propiedades para fortran
		radar_f90.yll =  4.9
		radar_f90.xll = -76.82
		radar_f90.dx = 0.0015
		radar_f90.dxp=125.0
		radar_f90.nodata = -999
		radar_f90.ncols = 1728
		radar_f90.nrows=1728
		
		#si se pasan argumentos se cambia
		for key, value in kwargs.iteritems():      # styles is a regular dictionary
			setattr(radar_f90, key, value)
		
		#Copia las propiedades en la lista de propiedades variable
		RadProp.extend([radar_f90.ncols,
			radar_f90.nrows,
			radar_f90.xll,
			radar_f90.yll,
			radar_f90.dx,
			radar_f90.nodata])
		
	#lee Z y ref
	def read_bin(self,path):	
		'Descripcion: Lee los archivos del radar de reflectividad y los \n'\
		'convierte en Z, en el momento en que lee la imagen, esta \n'\
		'funcion asigna parametros al modulo radar \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'path : Ruta en la que se encuentra el binario de radar.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Z : Matriz con los valores de Z del radar.\n'\
		'ref : Matriz con los valores de reflectividad del radar.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'Z,ref=read_bin(path).\n'\
		
		f = open(path, "r")
		a = np.fromfile(f, dtype = np.uint8)
		f.close()
		reflectividad = (0.5*a)-32.0
		miss_val = (np.where(reflectividad == 95.5))[0]      
		if (len(miss_val) != 0):
			reflectividad[miss_val] = -999.0
		try:
		    ref = np.reshape(reflectividad,(1728,1728))
		except BaseException:
		    print "Binario defectuoso"
		ref[ref==-999]=0.0
		ref=np.flipud(ref)
		ref[ref<5]=0	
		Z = 10.**(ref/10.)
		Z[ref==-999]=-999
		self.Z=Z; self.ref=ref
	def read_netcdf(self,path):
		fid=Dataset(path,'r')		
		ref=fid.variables['DBZ_H'][0][0]
		ref[ref<5]=0
		Z = 10.**(ref/10.)
		Z[ref==-999]=-999
		#Actualiza condiciones para fortran
		radar_f90.nodata = -9999
		radar_f90.ncols = ref.shape[0]
		radar_f90.nrows=ref.shape[1]
		self.Z=Z; self.ref=ref
	#Obtiene el binario
	def detect_clouds(self,k=np.array([[-1,-2,-1],[0,0,0],[1,2,1]]),
		umbral=40):		
		'Descripcion: Obtiene los bordes de una imagen entregada  \n'\
		'puede retornar los bordes y entregar la imagen rellenada \n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen raster con informacion de Z o de ref.\n'\
		'k : forma del filtro que se utiliza para obtener los bordes.\n'\
		'umbral : Umbral utilizado para determinar si un objeto es o no.\n'\
		'	un borde.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'borders : Matriz con los bordes detectados.\n'\
		'binario : Matriz binaria con los bordes rellenos.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'borders,binario=detect_clouds(Z).\n'\
		
		#Obtiene el gradiente, y el binario
		gradiente=radar_f90.detect_clouds(self.Z,k,k.T,radar_f90.ncols,radar_f90.nrows)
		bordes=np.zeros(gradiente.shape)
		bordes[gradiente>umbral]=1
		#Llena el binario y quita el ruido
		binario=nd.binary_fill_holes(bordes)*1
		#Retorna lo obtenido
		self.borders=bordes
		self.binario=binario		
	#Clasifica el binario	
	def classify(self,umbral=100,bordes='yes'):
		'\n'\
		'Descripcion: Limpia la imagen binaria de de acuerdo al \n'\
		'	tamano de los objetos que la componen, todo objeto menor\n'\
		'	al umbral es eliminado y clasifica la imagen\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self.binario : Imagen binaria de radar.\n'\
		'umbral : cantidad minima de pixeles para determinar si borrar o no.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self.classes : los objetos clasificados de la imagen.\n'\
		'self.elements : Vector 3xN con las celdas de los objetos.\n'\
		'self.cant_elem : Vector con la cantidad de celdas de cada objeto.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'.\n'\
	
		#clasifica
		classified,npixels,nelem=radar_f90.classify_binary(self.binario,
			radar_f90.ncols,radar_f90.nrows)
		elem,tam= radar_f90.cut_list_object(nelem,npixels)
		#Itera sobre los elementos buscando nubes menores al umbral
		imageOut = radar_f90.clean_by_size(classified,
			elem,tam,umbral,
			radar_f90.ncols,
			radar_f90.nrows,
			elem.shape[1],
			elem[0,-1])
		#Clasifica las imagenes
		classified,npixels,nelem=radar_f90.classify_binary(imageOut,
			radar_f90.ncols,radar_f90.nrows)
		elem,tam= radar_f90.cut_list_object(nelem,npixels)
		self.classes=classified
		self.elements=elem
		self.cant_elem=tam[0]
		#Clasifica solo bordes
		if bordes=='yes':
			classified,npixels,nelem=radar_f90.classify_binary(
				imageOut-self.image.erosion(imageOut,kernel=5),
				radar_f90.ncols,radar_f90.nrows)
			elem,tam= radar_f90.cut_list_object(nelem,npixels)
			self.border_class=classified
			self.border_elem=elem
			self.border_cant=tam[0]
	#Encuentra distancia de los objetos, distribucion y maxima
	def find_lenght(self):
		'\n'\
		'Descripcion: Determina la longitud de los objetos \n'\
		'	estas longitudes son de extremo a extremo\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'self : la clase debe estar clasificada.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'MaxLenght : Maxima longitud en cada objeto.\n'\
		'DistLenght : Distribucion de los deciles de las distancias.\n'\
		'	en cada objeto .\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'MaxL,DL=find_lenght(binaria).\n'\
			
		#Saca longitudes de todas las nubes
		Distlenght,Maxlenght = radar_f90.objects_lenght(self.border_elem,
			self.border_cant,radar_f90.ncols,radar_f90.nrows,
			self.border_elem.shape[1],self.border_cant.shape[0])
		self.MaxLenght=Maxlenght
		self.DistLenght=Distlenght
	#Geometria basica de los objetos
	def Basics_Geometry(self,coordType='LatLong'):
		'\n'\
		'Descripcion: Calcula el area, el perimetro y el centro\n'\
		'	de cada uno de los objetos encontrados en la clasificada\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen binaria del radar.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Area : Areas de los objetos.\n'\
		'Perim : Perimetros de los objetos.\n'\
		'centroMasa : Coordenadas X,Y de los centros de cada objeto.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'A,P,xy=Basics_Geometry(binImage).\n'\
		
		#calcula area
		self.area=self.cant_elem*radar_f90.dxp**2/1e6
		#calcula centro de masa de cada elemento
		self.centroMasa=[]
		for i in range(len(self.cant_elem)):
			Y=np.median(self.elements[:,np.where(self.elements[0]==i+1)][1][0])
			X=np.median(self.elements[:,np.where(self.elements[0]==i+1)][2][0])
			if coordType=='LatLong':
				X=radar_f90.xll+radar_f90.dx*(X-0.5)
				Y=radar_f90.yll+radar_f90.dx*((radar_f90.nrows-Y)+0.5)
			self.centroMasa.append([X,Y])
		self.centroMasa=np.array(self.centroMasa).T
		#Calcula perimetro como una aproximacion a 1.8 veces la cantidad de 
		#recuadros que conforman un borde		
		self.perim=self.border_cant*radar_f90.dxp*1.8/1e6		
	#Dimension fractal plana
	def Fractal_Dimension_Plain(self):
		'\n'\
		'Descripcion: Estima la dimension fractal del objeto \n'\
		'	analizado como un binario y no como una superficie\n'\
		'	para superficie usar: Fractal_Dimension_Surface\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'ObectNum : Cantidad de pixeles que comoponen a cada objeto.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'FractalD : Dimension fractal plana.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'FD_plana=Fractal_Dimension_Plain(tam).\n'\
				
		if radar_f90.dxp<>0.0:
			return np.log(self.cant_elem)/np.log(radar_f90.dxp)
		else:
			print 'Error: radar.dxp=0.0'
	#Dimension fractal superficie
	def Fractal_Dimension_Surface(self,k=12,a=1):
		'\n'\
		'Descripcion: Determina la dimension fractal superficial\n'\
		'	de los objetos\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen de Z o de reflectividad del radar.\n'\
		'ObjectList : Lista de pixeles que componen los objetos a analizar.\n'\
		'k : Tamano del kernel con el que se va a realizar el calculo.\n'\
		'a : Factor de incremento de la desviacion en el metodo de conteo.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'FD_surf : Matriz con lsa dimensiones Z, indicando la dimension fractal .\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'FD=Fractal_Dimension_Surface(Z,elem,k=10,a=2).\n'\
		
		self.fractal = radar_f90.fractal3d(self.Z,self.elements,k,a,
			radar_f90.ncols,
			radar_f90.nrows,
			self.elements.shape[1])		
	#Pendiente con respecto a los vecinos 
	def slope_arc(self):	
		'\n'\
		'Descripcion: Calcula la pendiente en cada celda como si \n'\
		'	se tratara de un DEM\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen de reflectividad o Z.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'slope : Pendiente calculada.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'.\n'\
	
		return radar_f90.arc_slope(self.Z,radar_f90.ncols,radar_f90.nrows)
	#Obtiene media y desbviacion de una variable
	def var2mean(self,RasterValues):
		'\n'\
		'Descripcion: Toma la imagen clasificada y una varaible\n'\
		'	distribuida a partir de esta segunda calcula la media\n'\
		'	y la desviacion de la variable en cada objeto\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'ObjectList : Lista de 3XN con las celdas que componen los objetos.\n'\
		'RasterValues : Matriz con los valores raster a promediar.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'meanVar : Vector con las medias de cada objeto.\n'\
		'stdVar : Vector con las desviaciones de cada objeto.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'.\n'\
		
		meanvar,stdVar = radar_f90.var2mean(RasterValues,
			self.elements,self.elements[0,-1],
			radar_f90.ncols,
			radar_f90.nrows,
			self.elements.shape[1])
		return meanvar,stdVar
	#Clasifica convectivas y stratiformes
	def Class2ConvStrat_deprecated(self,umbral=None,kernel=3):
		strat=np.copy(self.fractal)
		conv=np.copy(self.fractal)
		strat[strat>0]=1; strat[np.isnan(strat)]=0
		conv[np.isfinite(conv)]=0; conv[np.isnan(conv)]=1
		if umbral<>None:
			conv[self.ref>umbral]=1
		#conv=self.image.dilation(conv,kernel=kernel)
		conv=nd.binary_fill_holes(conv); conv=conv*2
		self.StratConv=strat+conv
		self.StratConv[self.StratConv>2]=2
	#Clasifica convectivas y estratiformes por el metodo de steiner
	def Class2ConvStratiform(self, umbral = 40, radio = 11, metodo = 'yuter',
		ZminSiriluk = 15, a_yuter = 10, b_yuter = 50):
		'\n'\
		'Descripcion: Clasifica entre convectivo y estratiforme \n'\
		'	de acuerdo a las metodologias de Steiner, Siriluk o Yuter\n'\
		'	todas se basan en el mismo principio, sin embargo varian\n'\
		'	algunos parametros.\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'umbral : Cantidad minima de reflectividad para considerar picos.\n'\
		'radio: Radio de busqueda de centros convectivos.\n'\
		'metodo: Metodo de busqueda: yuter, siriluk o steiner.\n'\
		'ZminSiriluk: Valor minimo de Zc del metodo de siriluk.\n'\
		'a_yuter: Valor de a del metodo de yuter.\n'\
		'b_yuter: Valor de b del metodo de yuter.\n'\
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'self.ConvStra : Imagen clasificada en Convectivos, Estratiformes.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'.\n'\
		#Parsea el metodo
		if metodo == 'steiner':
			metNum = 1
		elif metodo == 'yuter':
			metNum = 3
		elif metodo == 'siriluk':
			metNum = 2
		#Invoca la funcion de fortran 
		peaks,self.ConvStra = radar_f90.steiner_find_peaks(self.ref,
			umbral,radio,metNum,ZminSiriluk,a_yuter,b_yuter,
			int(self.ref.shape[0]), int(self.ref.shape[1]))
	
	def save_rain_class(self, ruta):
		gr = Dataset(ruta,'w',format='NETCDF4')
		#Diccionario de propiedades
		Dict = {'ncols':RadProp[0],
		    'nrows': RadProp[1],
		    'xll': RadProp[2],
		    'yll': RadProp[3],
		    'dx': RadProp[4]}
		#Establece tamano de las variables 
		DimNcol = gr.createDimension('ncols',self.ConvStra.shape[0])
		DimNfil = gr.createDimension('nrows',self.ConvStra.shape[1])
		#Crea variables
		ClasStruct = gr.createVariable('Conv_Strat','i4',('ncols','nrows'),zlib=True)
		ClasRain = gr.createVariable('Rain', 'i4', ('ncols','nrows'),zlib=True)
		ClasRainHigh = gr.createVariable('Rhigh', 'i4', ('ncols','nrows'),zlib=True)
		ClasRainLow = gr.createVariable('Rlow', 'i4', ('ncols','nrows'),zlib=True)
		#Asigna valores a las variables
		ClasStruct[:] = self.ConvStra
		#Lluvia normal
		ppt = np.copy(self.ppt['media']) * 1000
		ppt = ppt.astype(float)
		ClasRain[:] = ppt
		#Lluvia alta
		ppt = np.copy(self.ppt['alta']) * 1000
		ppt = ppt.astype(float)
		ClasRainHigh[:] = ppt
		#Lluvia baja
		ppt = np.copy(self.ppt['baja']) * 1000
		ppt = ppt.astype(float)
		ClasRainLow[:] = ppt		
		#Cierra el archivo 
		gr.setncatts(Dict)
		gr.close()
				
	# Genera kernels circulares 
	def CircKernel(self,radio):
		#Centro del kernel y cantidad de datos
		a, b = radio, radio
		n = radio*2+1
		#Obtiene el kernel
		y,x = np.ogrid[-a:n-a, -b:n-b]
		mask = x*x + y*y <= radio*radio
		#Lo vuelve numerico
		array = np.zeros((n, n))
		array[mask] = 1
		return array
	# Convierte reflectividad a lluvia por Sepulveda, 2015
	def DBZ2Rain (self,path=Path):
		str_disdro = '77'
		ajuste_multicapaall = __open_pklfiles__(path)
		
		try:
			aa = len(self.ref)
			ref2 = np.array(self.ref)
		except:
			ref2 = np.array([self.ref])
		
		dbz = ref2
		trunc_dbz = np.copy(dbz)
		trunc_dbz[trunc_dbz >= 40] = 40.0
		
		## Se hace el ajuste para la reflectividad horizontak
		
		trunc_dbz2 = np.copy(dbz)
		trunc = 55. 
		
		c1_all = (1/ajuste_multicapaall['capa1']['mc1'])*trunc_dbz2 - (ajuste_multicapaall['capa1']['bc1']/ajuste_multicapaall['capa1']['mc1'])
		
		if trunc is not None:
			c1_all[c1_all >= trunc] = trunc       
		c2_all = ((10**(c1_all/10.0))/ajuste_multicapaall['capa2']['bc2'])**(1.0/(ajuste_multicapaall['capa2']['mc2']))
		
		self.ppt = {}
		for name, kmin, kmax in zip(['media','baja','alta'],['mc3','Emin_m','Emax_m'],['bc3','Emax_m','Emax_b']):
			#Media		
			c3_all = (ajuste_multicapaall['capa3'][kmin]*c2_all + ajuste_multicapaall['capa3'][kmax])
			c3_all[c1_all <= 5.0] = 0.0
			mask_ajust = np.isnan(c3_all)
			c3_all[mask_ajust == True] = 0 
			c3_all[c3_all < 0] = 0		
			## Obtiene las variables de la lluvia
			self.ppt.update({name:c3_all})
	
class draw_func:
	def __init__(self):
		pass
	
	#Hace un plot elegante de la imagen de radar	
	def plot_radar_elegant(self,imageIn,ruta=None,figsize=(12,12),extra_lat=-0.2,
		extra_long=-0.2,lines_spaces=0.4,mask_value=0.0,xy=None,
		xyColor='red',colorbar=True, texto=None, **kwargs):
		'\n'\
		'Descripcion: Toma cualquier raster con valores del radar\n'\
		'	y hace un plot de este de forma elegante\n'\
		'\n'\
		'Parametros\n'\
		'----------\n'\
		'imageIn : Imagen con valores o binaria del radar\n'\
		'ruta : si es especificada guarda la imagen.\n'\
		'figsize : Tamano de la figura.\n'\
		
		'\n'\
		'Retornos\n'\
		'----------\n'\
		'Plotea la figura y si ruta <> None, la guarda.\n'\
		'\n'\
		'Ejemplo\n'\
		'----------\n'\
		'.\n'\
			
		#Obtiene las latitudes y longitudes
		longs=np.array([radar_f90.xll+(i+0.5)*radar_f90.dx for i in range(radar_f90.ncols)])
		lats=np.array([radar_f90.yll+(i+0.5)*radar_f90.dx for i in range(radar_f90.nrows)])
		X,Y=np.meshgrid(longs,lats)
		Y=Y[::-1]
		#Cambia el tamano de la figura si tiene que hacerlo
		if figsize<>None:
			fig=pl.figure(figsize=figsize)		
		#Genera el lugar de ploteo con el Basemap
		m = Basemap(projection='merc',
			llcrnrlat=lats.min()-extra_lat,
			urcrnrlat=lats.max()+extra_lat,
			llcrnrlon=longs.min()-extra_long,
			urcrnrlon=longs.max()+extra_long,
			resolution='c')
		#Grafica las lineas horizontales y verticales
		m.drawparallels(np.arange(lats.min(),
			lats.max(),lines_spaces),
			labels=[1,0,0,0],
			fmt="%.2f",
			rotation='vertical',
			xoffset=0.1,
			linewidth=0.1)
		m.drawmeridians(np.arange(longs.min(),
			longs.max(),lines_spaces),
			labels=[0,0,1,0], 
			fmt="%.2f",
			yoffset=0.1,
			linewidth=0.1)
		#Cambia zeros por nana
		if mask_value<>None:
			imageIn=np.ma.array(imageIn,mask=imageIn==mask_value)
		#Genera el mapa
		demX,demY=m(X,Y)
		cs=m.contourf(demX, demY, imageIn, 100, **kwargs)
		if colorbar:
			cbar = m.colorbar(cs,location='bottom',pad="5%")	
		#dibuja los circulos		
		XY = []
		for i in [30,60,90,120,250]:			
			x,y = self.draw_circle(m,-75.5276,6.1937,i,c='k',lw=0.5)
			XY.append([x,y])
		#Si hay coordenadas de algo las dibuja
		if xy<>None:
			xc,yc=m(xy[0],xy[1])
			m.plot(xc,yc,color=xyColor,
				#s=30,
				linewidth=1,)
				#edgecolor='black')
		if texto<>None:
			pl.annotate(texto, xy=(0.1, 0.9), xycoords='axes fraction', size=16)
		if ruta<>None:
			pl.savefig(ruta,bbox_inches='tight')
		pl.show()
		#return XY
		
	def draw_circle(self,m, centerlon, centerlat, radius, *args, **kwargs):
	    glon1 = centerlon
	    glat1 = centerlat
	    X = []
	    Y = []
	    for azimuth in range(0, 360):
	        glon2, glat2, baz = self.shoot(glon1, glat1, azimuth, radius)
	        X.append(glon2)
	        Y.append(glat2)
	    X.append(X[0])
	    Y.append(Y[0])
	    #m.plot(X,Y,**kwargs) #Should work, but doesn't...
	    Xn,Yn = m(X,Y)
	    pl.plot(Xn,Yn,**kwargs)
	    return X,Y
	
	def shoot(self,lon, lat, azimuth, maxdist=None):
	    """Shooter Function
	    Original javascript on http://williams.best.vwh.net/gccalc.htm
	    Translated to python by Thomas Lecocq
	    """
	    glat1 = lat * np.pi / 180.
	    glon1 = lon * np.pi / 180.
	    s = maxdist / 1.852
	    faz = azimuth * np.pi / 180.
	 
	    EPS= 0.00000000005
	    if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
	        alert("Only N-S courses are meaningful, starting at a pole!")
	 
	    a=6378.13/1.852
	    f=1/298.257223563
	    r = 1 - f
	    tu = r * np.tan(glat1)
	    sf = np.sin(faz)
	    cf = np.cos(faz)
	    if (cf==0):
	        b=0.
	    else:
	        b=2. * np.arctan2 (tu, cf)
	 
	    cu = 1. / np.sqrt(1 + tu * tu)
	    su = tu * cu
	    sa = cu * sf
	    c2a = 1 - sa * sa
	    x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
	    x = (x - 2.) / x
	    c = 1. - x
	    c = (x * x / 4. + 1.) / c
	    d = (0.375 * x * x - 1.) * x
	    tu = s / (r * a * c)
	    y = tu
	    c = y + 1
	    while (np.abs (y - c) > EPS):
	 
	        sy = np.sin(y)
	        cy = np.cos(y)
	        cz = np.cos(b + y)
	        e = 2. * cz * cz - 1.
	        c = y
	        x = e * cy
	        y = e + e - 1.
	        y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
	              d / 4. - cz) * sy * d + tu
	 
	    b = cu * cy * cf - su * sy
	    c = r * np.sqrt(sa * sa + b * b)
	    d = su * cy + cu * sy * cf
	    glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
	    c = cu * cy - su * sy * cf
	    x = np.arctan2(sy * sf, c)
	    c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
	    d = ((e * cy * c + cz) * sy * c + y) * sa
	    glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
	 
	    baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
	 
	    glon2 *= 180./np.pi
	    glat2 *= 180./np.pi
	    baz *= 180./np.pi
	 
	    return (glon2, glat2, baz)
	
