### VACIOS DE INFROMACION###
rm(list = ls(all = TRUE))
gc()

library(maptools)
library(sp)
library(rasclass)
library(ape)
library(raster)
library(spatstat)
library(rgdal)
library("classInt")


# ############ FUNCIONES############ --------------------------------------
#generar tablas resumen
tabla_resumen=function (campo1,campo_shp, nombre,join_TABLA,shp,rango){
  tabla1=as.data.frame(summary(TABLA[,campo1]))
  tabla2=as.data.frame(matrix(0,nrow(tabla1),2))
  tabla2[,1]=as.character(rownames(tabla1))
  tabla2[,2]=as.numeric(tabla1[,1])
  colnames(tabla2)=c(nombre,"CUENTA")
  
  ## llena  cuenta en shp
  #shp$campo_shp=as.character(shp$campo_shp)
  shp$CUENTA = 0
  for (i in tabla2[,1])
  {shp$CUENTA[shp@data[,campo_shp]==i]=tabla2$CUENTA[tabla2[,1]==i]}
  
  
  ###genera tabla  tabla3 con resumen 
  
  nombres_mun=unique(TABLA[,rango])
  tabla3=merge(tabla2,nombres_mun, by.x=nombre,by.y=join_TABLA) ### para G_F_funtion
  ix=order(tabla3$CUENTA,decreasing=T)
  tabla3=tabla3[ix,]
  tabla3$por=cumsum(tabla3$CUENTA/sum(tabla3$CUENTA))
  resultado=c(shp,tabla3)
  return(resultado)
  
}

# funcion G y F  
G_F_FUNTION=function(x,t){
  SP=ppp(x$longitud,x$latitud,window=ventana)
  plot(SP,axes=T, main=t) 
  G=Gest(SP)
  plot.fv(G,main=t)###Agrupado = G aumenta rápidamente a corta distancia
  ###Uniformidad = G aumenta lentamente hasta una distancia donde la mayoría de los eventos separados, entonces aumenta rápidamente
  Geva=envelope(SP,Gest,nsim=10+nrow(x))
  plot(Geva, main=paste("Gfun",t))
  
  Ffun=Fest(SP)
  #plot.fv(Ffun, main=paste("Ffun",t))
  Feva=envelope(SP,Fest,nsim=10+nrow(x))
  plot(Feva, main=paste("Ffun",t))
  
}

# ##################### .....CARGAR  DATOS A ANALIZAR....## ---------------

# definir las rutas de trabajo 

ruta_datos="~/VACIOS DE INFROMACION/GSI"
ruta_salida= "~/GBIF/VACIOS/DESCRIPTIVA"
ruta_unidad_espacial="~/VACIOS DE INFROMACION/INFO_GEO"

##subir base de  datos de registros a analisar
setwd(ruta_datos)

regcol=read.table("reg_col.txt",header=T)# base de datos original

PUNTOS=read.table("DATOS.txt",header=T) # base de datos depurada


#Puntos en formato ShapePoints (sobreposisciones)
PUNTOS2=PUNTOS
#PUNTOS2$latitud=0
#PUNTOS2$latitud=PUNTOS2$latitud1
coordinates(PUNTOS2)=~longitud+latitud


#Puntos  cordenadas unicas (funcion g_F)
PUNTOS3=as.data.frame(unique(cbind(PUNTOS$longitud,PUNTOS$latitud)))
colnames(PUNTOS3)=c("longitud","latitud")

PUNTOS4=PUNTOS3 ## Cordendas unicas en shape points
coordinates(PUNTOS4)=~longitud+latitud

# subir  unidades espaciales de analisis shp

setwd(ruta_unidad_espacial)

col=readShapePoly("COLOMBIA.shp")
mpios=readShapePoly("mun_2011_wgs84_100k.shp") 
dptos=readShapePoly("DEPTOS.shp")
PNN=readShapePoly("limspnn_100k2010.shp")
AP=readShapePoly("AreasProtegidas_RUNAP_PROJECT.shp")
CAR=readShapePoly("Car2010_Project.shp")
BIOMAS=readShapePoly("BIOMA_Dissolve.shp")
PARAMO=readShapePoly("paramo_project.shp")
mundo=readShapePoly("C:/Program Files (x86)/ArcGIS/Desktop10.1/ArcGlobeData/continent.shp")

# subir  unidades espaciales de analisis raster

ALT=raster("alt.asc")
COLOMBIA=rasterize(col,ALT)
MUNICIPIOS=raster("MUNICIPIOS_RASTER.tif")
DEPARTEMENTOS=raster("DEPARTAMENTOS.tif")
PARQUES=raster("PNN.tif")
AREA_PROTE=raster("AREA_PROT.tif")
CAR_RASTER=raster("CAR_RASTER.tif")
BIOMA_RASTER=raster("BIOMA_RAS.tif")
PARAMO_RASTER=raster("PARAMO_RAS.tif")

AREA=c(COLOMBIA,MUNICIPIOS,DEPARTEMENTOS,PARQUES,AREA_PROTE,CAR_RASTER,BIOMA_RASTER, PARAMO_RASTER) ##BIOMA_RASTER


# #####################  0. DISTRIBUCION  BASE DE DATOS ###################



setwd(ruta_salida)

pdf("MAPAS RESUMEN.pdf ")

##Distribucion base de datos original 

PORCENTAJE=nrow(PUNTOS)/nrow(regcol) *100
coordinates(regcol)=~longitud+latitud

plot(mundo,main="UBICACION DE LOS REGISTROS BASE DE DATOS INICIAL")
points(x=coordinates(regcol)[,1],y=coordinates(regcol)[,2], col=rgb(139,0,0,100,maxColorValue=255),cex=1,pch=20)

## Distribucion base de datos depurada

plot(col)
points(x=coordinates(PUNTOS2)[,1],y=coordinates(PUNTOS2)[,2], col=rgb(139,0,0,100,maxColorValue=255),cex=1,pch=20)


# ##################### 1. GENERAR TABLAS RESUMEN  ################ 
# ####Realizar sobreposiciones  ### ---------------------------------------

OV1=overlay(mpios,PUNTOS2)
cntr1=as.character(OV1$ID_MUN)
nom1=as.character(OV1$MPIOS)
cntr1.1=as.numeric(cntr1)

OV2=overlay(dptos,PUNTOS2)
cntr2=as.character(OV2$ID_DEPART)
nom2=as.character(OV2$DEPARTAMEN)
cntr2.1=as.numeric(cntr2)

OV3=overlay(PNN,PUNTOS2)
cntr3=as.double(OV3$ID)
nom3=as.character(OV3$NOM)

OV4=overlay(AP,PUNTOS2)
cntr4=as.double(OV4$IDAP)
nom4=as.character(OV4$NOMAP)


OV5=overlay(CAR,PUNTOS2)
cntr5=(OV5$ID)
nom5=as.character(OV5$IDCAR)

OV6=overlay(BIOMAS,PUNTOS2)
cntr6a=(OV6$IDBioma)
cntr6b=(OV6$IDEco)
nom6a=as.character(OV6$Bioma)

OV7=overlay(PARAMO,PUNTOS2)
cntr7=(OV7$NUM)
nom7=as.character(OV7$COMPLEJO)

# #Generar la tabla conjunta para cada  uno de los analisis ---------------

COORDENADAS=cbind(PUNTOS$longitud,PUNTOS$latitud)
colnames(COORDENADAS)=c("longitud","latitud")


TABLA=as.data.frame(cbind(COORDENADAS,cntr1,nom1,cntr2,nom2,cntr3,nom3,cntr4,nom4,cntr5,nom5,cntr6a,nom6a,cntr7,nom7 ))###valores con nombres
TABLA2=as.data.frame(cbind(COORDENADAS,cntr1.1,cntr2.1,cntr3,cntr4,cntr5,cntr6a,cntr7))##solo valores para el Moran
TABLA3=as.data.frame(cbind(nom1,nom2,nom3,nom4,nom5,nom6a,nom7))## para G and F function

colnames(TABLA)=c("LONGITUD","LATITUD","ID_MUN","MUNICIPIO","ID_DEP","DEPARTAMENTO","ID_PNN","PNN","ID_AP","AP","ID_CAR","CAR","ID_BIOMA","BIOMA","ID_PARAMO","COMPLEJO")
colnames(TABLA2)=c("LONGITUD","LATITUD","ID_MUN","ID_DEP","ID_PNN","ID_AP","ID_CAR","ID_BIOMA","ID_PARAMO")
colnames(TABLA3)=c("MUNICIPIO","DEPARTAMENTO","PNN","AP","CAR","BIOMA","COMPLEJO")

# Generar tablas resumen -------------------------------------------



## para municipio#  

res_muni=tabla_resumen(3,7,"ID_MUNICIPIO","ID_MUN",mpios,3:4)
mpios=res_muni[[1]]
res_muni3=data.frame(res_muni[[2]],(res_muni[[3]]),as.character(res_muni[[4]]),res_muni[[5]],stringsAsFactors=FALSE)
nombres=names(res_muni)
names(res_muni3)=nombres[2:5]

## para departamento#
res_dept=tabla_resumen(6,1,"ID_DEPARTAMENTO","DEPARTAMENTO",dptos,5:6)
dptos=res_dept[[1]]
res_dept3=data.frame(cbind(res_dept[[2]],as.double(res_dept[[3]]),res_dept[[4]],res_dept[[5]]),stringsAsFactors=FALSE)
nombres=names(res_dept)
names(res_dept3)=nombres[2:5]

## para PNN#
res_pnn=tabla_resumen(8,1,"ID_PNN","PNN",PNN,7:8)
PNN=res_pnn[[1]]
res_pnn3=data.frame(cbind(res_pnn[[2]],as.double(res_pnn[[3]]),res_pnn[[4]],res_pnn[[5]]),stringsAsFactors=FALSE)
names(res_pnn3)=c("PNN","CUENTA","ID_PNN","por")

## para AP#
res_AP=tabla_resumen(9,3,"ID_AP","ID_AP",AP,9:10)
AP=res_AP[[1]]
res_ap3=data.frame(cbind(res_AP[[2]],as.double(res_AP[[3]]),as.character(res_AP[[4]]),res_AP[[5]]),stringsAsFactors=FALSE)
nombres=names(res_AP)
names(res_ap3)=nombres[2:5]

## para CAR#
res_CAR=tabla_resumen(12,3,"ID_CAR","CAR",CAR,11:12)
CAR=res_CAR[[1]]
res_car3=data.frame(cbind(res_CAR[[2]],as.double(res_CAR[[3]]),res_CAR[[4]],res_CAR[[5]]),stringsAsFactors=FALSE)
names(res_car3)=c("CAR","CUENTA","ID_CAR","por")

## para BIOMAS#
res_BIO=tabla_resumen(14,2,"ID_BIO","BIOMA",BIOMAS,13:14)
BIOMAS=res_BIO[[1]]
res_bio3=data.frame(cbind(res_BIO[[2]],as.double(res_BIO[[3]]),res_BIO[[4]],res_BIO[[5]]),stringsAsFactors=FALSE)
nombres=names(res_BIO)
names(res_bio3)=nombres[2:5]


## PARA COMPLEJOS PARAMOS#
res_para=tabla_resumen(16,5,"COMPLEJO","COMPLEJO",PARAMO,15:16)
PARAMO=res_para[[1]]
res_para2=data.frame(cbind(res_para[[2]],as.double(res_para[[3]]),res_para[[4]],res_para[[5]]),stringsAsFactors=FALSE)
nombres=names(res_para)
names(res_para2)=nombres[2:5]

# PLOTTEAR LOS MAPAS RESUMEN -----------------------------------------


spplot(mpios,zcol="CUENTA",col.region=rainbow(25),main="MUNICIPIOS") ##arreglar colores y poner tabla res2 como  leyenda
spplot(dptos,zcol="CUENTA",main="DEPARTAMENTOS")
spplot(PNN,zcol="CUENTA",main="PNN")
spplot(AP,zcol="CUENTA",main="AREA PROTEGIDA")
spplot(CAR,zcol="CUENTA",main="CAR")
spplot(BIOMAS,zcol="CUENTA",main="BIOMAS")
spplot(PARAMO,zcol="CUENTA",main="COMPLEJOS_PARAMO")

dev.off()

# GUARDAR RESULTADOS TABLA RESUMEN ----------------------------------




# guardar tablas 

write.table(res_muni3,"resumen_municipio.txt",sep= "\t", row.names=F)
write.table(res_dept3,"resumen_departamento.txt",sep= "\t",row.names=F)
write.table(res_pnn3,"resumen_parques_naturales.txt",sep= "\t",row.names=F)
write.table(res_ap3,"resumen_areas_protegidas .txt",sep= "\t",row.names=F)
write.table(res_car3,"resumen_cars.txt",sep= "\t",row.names=F)
write.table(res_bio3,"resumen_biomas.txt",sep= "\t",row.names=F)
write.table(res_para2,"resumen_paramos.txt",sep= "\t",row.names=F)

# guardar  shp
writePolyShape(mpios,"MUNICIPIO_CUENTA.shp")
writePolyShape(dptos,"DEPARATAMENTO_CUENTA.shp")
writePolyShape(PNN,"PARQUES_CUENTA.shp")
writePolyShape(AP,"AP_CUENTA.shp")
writePolyShape(CAR,"CAR_CUENTA.shp")
writePolyShape(BIOMAS,"BIOMA_CUENTA.shp")
writePolyShape(PARAMO,"PARAMO_CUENTA.shp")

dev.off()



# ##################### 2. NEAREST  NEIGHBOUR G FUNTION AND  F-FUNTION  ################
NEAREST=NULL

pdf("G_F_PNN.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(2,3))

# ##para  colombia --------------------------------------------------------

t="colombia"
ventana=as(col,"owin")
AS=G_F_FUNTION(PUNTOS3,t)
A=area.owin(ventana)
n=dim(PUNTOS)[1]
SP=ppp(PUNTOS3$longitud,PUNTOS3$latitud,window=ventana)
DO=mean(nndist(SP))
DR=0.5/(sqrt(n/A))
ANN=DO/DR
SE=0.26136/(sqrt((n^2)/A))
Z=(DO-DR)/SE
p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
if (ANN<1){ patron="asociado"} else {patron="disperso"}
if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}


NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
NEAREST=rbind(NEAREST,NEAREST1)

# ###para municipio -------------------------------------------------------

NIVELES_MUN=res_muni3$ID_MUNICIPIO[res_muni3$por<=0.8] 
final_sp_mun=NULL
for (i in NIVELES_MUN){
  mask=mpios[which(mpios$ID_MUN==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_muni3$MUNICIPIO[res_muni3$ID_MUNICIPIO==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_mun=rbind(final_sp_mun,tabla_sp)
  
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("MUN",mask$DEPARTAMEN,mask$NOMBRE_ENT,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  n=dim(reg)[1]
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}

###para departamento ---------------------------------------------------
NIVELES_DEPT=res_dept3$ID_DEPARTAMENTO[res_dept3$por<=0.85] 
final_sp_dep=NULL

for (i in NIVELES_DEPT){
  mask=dptos[which(dptos$DEPARTAMEN==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_dept3$ID_DEPARTAMENTO[res_dept3$ID_DEPARTAMENTO==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_dep=rbind(final_sp_dep,tabla_sp)
  
  
  
  PUNTOS4=PUNTOS3
  coordinates(PUNTOS4)=~longitud+latitud
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("DEPARTAMENTO",mask$DEPARTAMEN,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  n=dim(reg)[1]
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}

# #### para PNN -----------------------------------------------------------
NIVELES_PNN=res_pnn3$PNN[res_pnn3$por<=0.85]  
final_sp_pnn=NULL

for (i in NIVELES_PNN){
  
  mask=PNN[which(PNN$NOM==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_pnn3$PNN[res_pnn3$PNN==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_pnn=rbind(final_sp_pnn,tabla_sp)
  
  PUNTOS4=PUNTOS3
  coordinates(PUNTOS4)=~longitud+latitud
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("PNN",mask$NOM,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  n=dim(reg)[1]
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}

# #### para AP ------------------------------------------------------------

NIVELES_AP=res_ap3$ID_AP[res_ap3$por<=0.85] ## NO FUNCIONA PARA OTUN QUIMBAYA NI PARA BARBAS BREMEN 
final_sp_ap=NULL

for (i in NIVELES_AP){
  mask=AP[which(AP$IDAP==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_ap3$AP[res_ap3$ID_AP==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_ap=rbind(final_sp_ap,tabla_sp)
  
  PUNTOS4=PUNTOS3
  coordinates(PUNTOS4)=~longitud+latitud
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("AREA_PROTEGIDA",mask$NOMAP,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  n=dim(reg)[1]
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}

# #### para CAR -----------------------------------------------------------

NIVELES_CAR=res_car3$CAR[res_car3$por<=0.85] ### bo funciona para corponor
NIVELES_CAR=NIVELES_CAR[-4]
final_sp_car=NULL

for (i in NIVELES_CAR){
  mask=CAR[which(CAR$IDCAR==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_car3$CAR[res_car3$CAR==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_car=rbind(final_sp_car,tabla_sp)
  
  PUNTOS4=PUNTOS3
  coordinates(PUNTOS4)=~longitud+latitud
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("CAR",mask$IDCAR,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  n=dim(reg)[1]
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}

# ###para PARAMOS-------------------------------------------------------

NIVELES_PARAMO=res_para3$COMPLEJO[res_para3$por<=0.85] 

final_sp_para=NULL

for (i in NIVELES_PARAMO){
  mask=PARAMO[which(PARAMO$COMPLEJO==i),]
  
  para_sp=overlay(PUNTOS2,mask)
  filas_sp=which(para_sp>=1)
  las_sp=PUNTOS[filas_sp,]
  nombre=res_para2$COMPLEJO[res_para2$ID_PARAMO==i]
  nombe_sp=unique(las_sp$nombre)
  tabla_sp=cbind(rep(nombre,length(nombe_sp)),as.character(nombe_sp))
  final_sp_para=rbind(final_sp_para,tabla_sp)
  
  PUNTOS4=PUNTOS3
  coordinates(PUNTOS4)=~longitud+latitud
  PuntosIN=overlay(PUNTOS4,mask)
  filas=which(PuntosIN>=1)
  ventana=as(mask,"owin")
  reg=PUNTOS3[filas,]
  t=paste("PARAMO",mask$COMPLEJO,sep="_")
  
  G_F_FUNTION(reg,t)
  
  ## radio de evidencia
  A=area.owin(ventana)
  n=dim(reg)[1]
  SP=ppp(reg$longitud,reg$latitud,window=ventana)
  DO=mean(nndist(SP))
  DR=0.5/(sqrt(n/A))
  ANN=DO/DR
  SE=0.26136/(sqrt((n^2)/A))
  Z=(DO-DR)/SE
  p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed.
  if (ANN<1){ patron="asociado"} else {patron="disperso"}
  if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
  NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
  NEAREST=rbind(NEAREST,NEAREST1)
}


dev.off()

# #### para BIOMAS --------------------------------------------------------


### llevaba mas de 12 horas y aun no habia podido crear la maskara , hice el dissolve
###pero queda con huecos  y me tira  error Error en owin(poly = opls) : 
####  Polygon data contain duplicated vertices, self-intersection and overlaps between polygons

# NIVELES_BIO=res_bio3$ID_BIOMA[res_bio3$por<=0.85] 
# 
# for (i in NIVELES_BIO){
#   mask=BIOMAS[which(BIOMAS$IDBioma==i),]
#   PUNTOS4=PUNTOS3
#   coordinates(PUNTOS4)=~longitud+latitud
#   PuntosIN=overlay(PUNTOS4,mask)
#   filas=which(PuntosIN>=1)
#   ventana=as(mask,"owin")
#   reg=PUNTOS3[filas,]
#   t=paste("BIOMA",mask$Bioma,sep="_")
#   
#   G_F_FUNTION(reg,t)
#   
#   ## radio de evidencia
#   A=area.owin(ventana)
#   n=dim(reg)[1]
#   SP=ppp(reg$longitud,reg$latitud,window=ventana)
#   DO=mean(nndist(SP))
#   DR=0.5/(sqrt(n/A))
#   ANN=DO/DR
#   SE=0.26136/(sqrt((n^2)/A))
#   Z=(DO-DR)/SE
#   p=(1-(pnorm(abs(Z))))*2 # The null hypothsis states that features are randomly distributed. 
#   if (ANN<1){ patron="asociado"} else {patron="disperso"}
#   if (p<=0.05){desicion="estadisticamente asociado"} else {desicion="estadisticamente aleatoreo"}
#   
#   NEAREST1=cbind(t,ANN,Z,p, patron,desicion)
#   NEAREST=rbind(NEAREST,NEAREST1)
# }

# ###GUARDAR RESULTADOS  AVARAGE NEAREST NEIGHBOR -------------------------


write.table(NEAREST,"AVARAGE NEAREST NEIGHBOR.txt", sep="\t")
write.table(final_sp_mun,"ESPECIES_MUNICIPIO.txt", sep="\t")
write.table(final_sp_dep,"ESPECIES_DEPARTAMENTO.txt", sep="\t")
write.table(final_sp_pnn,"ESPECIES_PNN.txt", sep="\t")
write.table(final_sp_ap,"ESPECIES AP.txt", sep="\t")
write.table(final_sp_car,"ESPECIES_CAR.txt", sep="\t")
write.table(final_sp_para,"ESPECIES_PARAMOS.txt", sep="\t")

# # ##################### 3.   %  DE  PIXELES CON VALOR -------------------

area_estudio=NA
resultado=NA
for (i in AREA){  
  con_punto=extract(i,PUNTOS4)  ## extraer las categorias que tienen puntos 
  categoria=unique(con_punto)
  
  
  for (j in categoria){
    
    unidad=match(i,j)
    Nmask=count(i,j) ## nuemro de pixeles en el area de estudio
    pixel=extract(unidad,PUNTOS4)  ##pixeles que tiene punto dentro del area de estudio
    Npixel=length(pixel[is.na(pixel)==F])
    Porcentaje=(Npixel/Nmask)*100
    
    nombre=paste(layerNames(i),j,sep="_")
    area_estudio=rbind(area_estudio,nombre)
    
    resultado=rbind(resultado,Porcentaje)
  }
  
  }    

Final=cbind(area_estudio,resultado)


## GUARDAR RESULTADO 

write.table(Final,"Pixel  con informacion.txt", sep="\t")



 




