rm(list=ls(all=T))
gc() ##para borrar todo lo que quede oculto en memoria

library(raster)
library("classInt")
library("maptools")
library("Matching")
library("flexmix")
library("rgdal")

################ FUNCIONES ########################

# A. FUNCION SESGO FACTORES
BIAS= function(shp,raster){
  
  ## sacar raster distancias
  raster_shp=rasterize(shp,raster)
  
  distancia=distance(raster_shp)
  ## sacar  raster reclasificacdo en 4 clases
  en_area=mask(distancia,colombia)
  var=na.omit(as.vector(en_area[,]))
  fishers=classIntervals(var,4,style="fisher")
  
  BRKS=fishers$brks
  clasificado=cut(en_area,breaks=BRKS,right = F)
  plot(clasificado,main=paste("niveles",nombre_capa,sep="_"))
  ## sacar indice de BIAS
  N=nrow(DATOS)
  valores=extract(clasificado,DATOS2)
  clases=as.data.frame(table(valores))
  sesgos=NA 
  for( i in 1:4){
    nd=clases[i,2]; if (is.na(nd)){nd=0}
    pd=count(clasificado,i)/ncell(clasificado)
    BIAS=(nd-pd*N)/(sqrt(pd*(1-pd)*N))
    sesgos=cbind(sesgos,BIAS)
  }
  sesgos=sesgos[,2:ncol(sesgos)]
  nombres=1:4
  names(sesgos)=nombres
  matris=as.matrix(sesgos)
  matris=cbind(rownames(matris),matris)
  niveles=reclassify(clasificado,matris)
  names(niveles)=nombre_capa
  plot(niveles, main=nombre_capa)
  plot(shp,add=T)
  raster::writeRaster(niveles, paste("sesgos",nombre_capa,".tif",sep=""), overwrite=TRUE) # revisar el nombre
  resultados=c(sesgos,niveles) 
  return(resultados)
}

# B. FUNCION  SESGOS AMBIENTALES

BIAS_AMBIENTAL= function(ras){   
  var= values(ras)
  var=na.omit(var)
  fishers=classIntervals(var,4,style="fisher")
  BRKS=fishers$brks
  clasificado=cut(ras,breaks=BRKS,right = F)
  
  ## sacar indice de BIAS
  #N=length(na.omit(unique(datos$celdas)))
  N=length(which(LOCALIDADES[,]>=1))
  valores=clasificado*LOCALIDADES
  #clases=as.data.frame(table(valores))
  sesgos=NA 
  for( i in 1:4){
    nd=count(valores,i)
    pd=count(clasificado,i)/length(which(clasificado[,]>=1))
    BIAS=(nd-(pd*N))/(sqrt(pd*(1-pd)*N))
    sesgos=cbind(sesgos,BIAS)
  }
  sesgos=sesgos[,2:ncol(sesgos)]
  nombres=1:4
  names(sesgos)=nombres
  matris=as.matrix(sesgos)
  matris=cbind(rownames(matris),matris)
  niveles=reclassify(clasificado,matris)
  names(niveles)=names(ras)
  plot(niveles, main=NOMBRES[k]) # revisar nombre
  raster::writeRaster(niveles, paste("ambiental",names(niveles),".tif",sep=""), overwrite=TRUE) # revisar el nombre
  resultados=c(sesgos,niveles) 
  return(resultados)
}

# C. FUNCION COMPLEMENTARIEDAD BASE DE DATOS

COMPLEMENTARIEDAD_BOOTS= function(datos, AREA){
  s=NA
  SELECCIONADAS=subset(datos,celdas==AREA)
  N=nrow(SELECCIONADAS)
  ESPECIES=unique(SELECCIONADAS$nombre)
  ESPECIES=na.omit(ESPECIES)
  Sobs=length(ESPECIES)
  for ( i in ESPECIES){
    Especie=subset(SELECCIONADAS,nombre==i)
    Pi=nrow(Especie)/N
    S_parcial=(1-Pi)^N# revisar si es  sobre el numero de especies o sobre N (N specpool )
    s=cbind(s,S_parcial)
  }   
  
  S= sum(s,na.rm=T) + Sobs
  C=Sobs/S
  names(C)=AREA
  return(C)
}


COMPLEMENTARIEDAD_JACK= function(datos, AREA){
  sj=NA
  SELECCIONADASj=subset(datos,celdas==AREA)
  ESPECIESj=unique(SELECCIONADASj$nombre)
  ESPECIESj=na.omit(ESPECIESj)
  Sobsj=length(ESPECIESj)
  n=length(unique(datos$celdas))
  cuenta=as.data.frame(table(datos$nombre))
  L=nrow(subset(cuenta, Freq==1))
  sj=Sobsj +L*((n-1)/n)
  Cj=Sobsj/sj
  names(Cj)=AREA
  return(Cj)
}


# D. FUNCION ESTANDARIZACION DE VALORES

ESTANDARIZACION=function(x){
  min_cont=min(getValues(x),na.rm=T)
  max_cont=max(getValues(x),na.rm=T)
  x=na.omit(x)
  ESTANDAR= (x-min_cont)/(max_cont - min_cont)
  return(ESTANDAR)
  
}

############### CARGAR DATOS ####################

ruta_datos="~/BIOACUSTICA/XENOCANTO_Ago30"
ruta_factores= "~/VACIOS DE INFROMACION/INFO_GEO"
ruta_ambientales="C:/Users/GIC 9/Documents/VACIOS DE INFROMACION/INFO_GEO/AMBIENTALES"

##subir base de  datos de registros a analisar
setwd(ruta_datos)
regcol=read.csv("para correr Analisis de vacios xenocanto 24.06.2013.csv",header=T)# base de datos original


DATOS=read.csv("para correr Analisis de vacios xenocanto 24.06.2013.csv",header=T) # base de datos depurada o DATOS=regcol
DATOS2=DATOS
coordinates(DATOS2)=~longitud+latitud


# subir  factores
setwd(ruta_factores)

AP=readShapePoly("AreasProtegidas_RUNAP_PROJECT.shp")
urbano=readShapePoints("casco_urbano_wgs84.shp")
rios=readShapePoly("rios_principales_wgs84.shp")
vias=readOGR(".", "vias_wgs84")
colombia=readShapePoly("COLOMBIA.shp") # corte del area de estudio
grilla=raster("grilla 10k.img")
mundo=readShapePoly("C:/Program Files (x86)/ArcGIS/Desktop10.1/ArcGlobeData/continent.shp")


# numerar celdas
nceldas=ncell(grilla)
grilla[1:nceldas]<-1:nceldas

#subir ambientales 
setwd(ruta_ambientales)
ambientales<-stack(paste("bio_",1:19,sep=""),"alt","slope")
NOMBRES=c("TEMPERATURA MEDIA ANUAL","MEDIA RANGO DIURNO","ISOTERMALIDAD","ESTACIONALIDAD TEMPERATURA","MAX T MES MAS CALIDO","MIN T MES FRIO","RANGO ANUAL T","T MEDIA  DEL CUARTO HUMEDO","T MEDIA DEL CUARTO SECO","T MEDIA DEL CUARTO CALIDO","T M DEL CUARTO FRIO","PRECIPITACIÓN ANUAL","PP MES MAS HUMEDO","PP MES MAS SECO","PP ESTACIONAL","PP CUARTO HUMEDO","PP CUARTO SECO","PP CUARTO CALIDO","PP CUARTO FRIO","ALTURA","PENDIENTE")

### arreglar capa ##
en_area=mask(grilla,colombia)
celdas=extract(en_area,DATOS2)
datos=cbind(DATOS,celdas)
localidades=rasterize(DATOS2,en_area, fun=function(x, ...) {length(unique(na.omit(x)))})
recla=as.matrix(cbind(1,localidades@data@max,1))
LOCALIDADES=reclassify(localidades,recla)



############## 1. DENSIDAD  ESTIMADA ###################

setwd(ruta_datos)

pdf("1.INDICE_densidad.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(1,1))



PORCENTAJE=nrow(DATOS)/nrow(regcol) *100

coordinates(regcol)=~longitud+latitud

plot(mundo,main="UBICACION DE LOS REGISTROS BASE DE DATOS INICIAL")
points(x=coordinates(regcol)[,1],y=coordinates(regcol)[,2], col=rgb(139,0,0,100,maxColorValue=255),cex=1,pch=20)

# Hallar la densidad de DATOS 
DATOS2$CONT=1



contador=rasterize(DATOS2,en_area,field=DATOS2$CONT, fun=sum, na.rm=T)
contador2=rasterize(DATOS2,en_area,field=DATOS2$CONT, fun=sum, na.rm=T,background=0)
contador3<-mask(contador,colombia)

plot(colombia, main="UBICACION DE LOS REGISTROS COLOMBIA")
plot(DATOS2,add=T)

plot(contador3,col= rev(topo.colors(5)),main="DENSIDAD DE  PUNTOS/ 10Km2")
plot(colombia,add=T)
density(contador3, plot=T)  ## mejorar plot  3d 
persp(contador3,xlab="X coordinates",ylab="y coordinates",zlab=
  "density",phi=35,theta=20,col="lightblue",expand=.5,ticktype="detailed")



dev.off()

raster::writeRaster(contador3,"INDICE_DENSIDAD.tif", overwrite=TRUE)
############# 2.  ANALISIS DE SESGOS ####################
pdf("2.1 INDICE_factores_sesgos.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(1,2))

DATOS_s=DATOS
DATOS_s=unique(DATOS[,5:6]) ##revisar cuales columnas son longitud y latitud
N=nrow(DATOS_s)


nombre_capa="AREA_PROTEGIDA"
sesgo_ap=BIAS(AP,grilla)

nombre_capa="RIOS"
sesgo_rios=BIAS(rios,grilla)

nombre_capa="VIAS"
sesgo_vias=BIAS(vias,grilla)

nombre_capa="CASCOS_URBANOS"
sesgo_urbano=BIAS(urbano,grilla)

dev.off()
#2.2.Sesgos ambientales 
# 2.2.1 Comparacion distribuciones varaibles ambientales kolmogorov y Kullback-Leibler Divergence

## guarda los valores en data frame
VARIABLES=ambientales[1:ncell(ambientales)]
VARIABLES=na.omit(VARIABLES) ##le quita nos na

MUESTREO=extract(ambientales,DATOS2)
MUESTREO=na.omit(MUESTREO)

pdf("2.2 INDICE_COMPARACION_VARIABLES.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(4,3))

Num_preds=dim(ambientales)[3]

RESULT=NULL

for (i in 1:Num_preds){
  VARIABLE=sample(VARIABLES[,i],nrow(MUESTREO),replace = T)
  #VARIABLE=sample(VARIABLES[,i],100)
  #MUESTREOs=sample(MUESTREO[,i],100)
  COMPARACION=ks.boot(MUESTREO[,i],VARIABLE) # ks.test( MUESTREOs,VARIABLE) #KOLMOGOROV
  VAR=NOMBRES[i]
  media_var=mean(VARIABLES[,i])
  media_DATOS=mean(MUESTREO[,i])
  des_var=sd(VARIABLES[,i])
  des_DATOS=sd(MUESTREO[,i])
  coefvar_var=des_var/media_var
  coefvar_DATOS=des_DATOS/media_DATOS
  D=COMPARACION$ks$statistic
  Pval=COMPARACION$ks.boot.pvalue
  if (Pval<=0.05){DESICION="Distribuciones diferentes"} else {DESICION="DISTRIBUCIONES IGUALES"}
  
  Nceldas=ncell(VARIABLES[,i])
  
  
  #Kullback-Leibler Divergence
  
  PI=as.data.frame(table(MUESTREO[,i])/nrow(MUESTREO))
  QTOTAL=as.data.frame(table(VARIABLES[,i])/Nceldas)
  QI=QTOTAL[match(PI[,1],QTOTAL[,1]),]
  
  
  Q=cbind(QI[,2],PI[,2])
  Q=as.matrix(Q)
  #KLdiv(Q,overlap=F, method="discrete")
  #DIVERGENCIA=KLdiv(Q,overlap=F,method="discrete",na.rm=T)[1,2]
  DENSI_VAR=density(VARIABLES[,i])
  DENSI_MUES=density(MUESTREO[,i])
  minimox=min(min(DENSI_VAR$x),min(DENSI_MUES$x))
  maxx=max(max(DENSI_VAR$x),max(DENSI_MUES$x))
  minimoy=min(min(DENSI_VAR$y),min(DENSI_MUES$y))
  maxy=max(max(DENSI_VAR$y),max(DENSI_MUES$y))
  
  
  plot(density(VARIABLES[,i]),col="black",xlab="",main=VAR,xlim=c(minimox,maxx), ylim=c(minimoy,maxy))
  lines(density(MUESTREO[,i]),col="red",main=VAR)
  
  
  RESUMEN=cbind(VAR,media_var,media_DATOS,des_var,des_DATOS,coefvar_var,coefvar_DATOS,D,Pval,DESICION)#,DIVERGENCIA)  
  
  RESULT=rbind(RESULT,RESUMEN)
  
}


## Guardar resultado

write.table(RESULT,"bondad ajuste variables.txt", sep="\t",col.names=T)

dev.off()



# 2.2.2. Capa de sesgos ambientales
pdf("2.3 INDICE_sesgos_aminbetales.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(1,2))

#sesgos_ambientales=raster(ambientales, layer=0)
sesgos_ambientales=raster(nrows= 20, ncols= 20, xmn= xmin(ambientales), xmx=xmax(ambientales), ymn=ymin(ambientales), ymx=ymax(ambientales))
res(sesgos_ambientales) <-  res(ambientales)
sesgos_ambientales[sesgos_ambientales] <- 0
resultados=NULL
i=1
Namb=dim(ambientales)[3]


for (k in 1:Namb){
  capa=ambientales[[k]]
  ses_amb=BIAS_AMBIENTAL(capa)
  sesgos_ambientales=stack(sesgos_ambientales,ses_amb[[5]])
  resultados=cbind(resultados,ses_amb[1:4])
}

resultados=as.matrix(resultados)
names(resultados)=NOMBRES

d=calc(sesgos_ambientales, sum, na.rm=T)

dev.off()

write.table(resultados,"sesgos_ambientales.txt",sep="\t",col.names=T)
raster::writeRaster(d,"INDICE_AMBIENTAL.tif", overwrite=TRUE)

############# 3.  COMPELMENTARIEDAD DE  LA BASE DE DATOS####################
pdf("3 INDICE_COMPLEMENTARIEDAD.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(1,1)) 

en_area=mask(grilla,colombia)
celdas=extract(en_area,DATOS2)
datos=cbind(DATOS,celdas)



nceldas= unique(datos$celdas)
nceldas=na.omit(nceldas)


lista<-NULL
for (k in nceldas){
  nombre= as.character(k)
  lista<-c(lista,COMPLEMENTARIEDAD_BOOTS(datos,k))
  #assign(paste(nombre,"comp",sep="_"),COMPLEMENTARIEDAD(datos,k)) 
} 
lista=as.data.frame(lista)
objeto=en_area
objeto[,]=0
objeto[nceldas]<-lista[,1]
OBJETO=mask(objeto,colombia)
plot(OBJETO, main="COMPLEMENTARIEDAD_BOOTS")
plot(colombia, add=T)
raster::writeRaster(OBJETO,"INDICE_COMPLEMENTARIEDAD_BOOTS.tif", overwrite=TRUE)



listaj<-NULL
for (t in 1:length(nceldas)){
  k=nceldas[t]
  nombre= as.character(k)
  listaj<-c(listaj,COMPLEMENTARIEDAD_JACK(datos,k))
  #assign(paste(nombre,"comp",sep="_"),COMPLEMENTARIEDAD(datos,k)) 
} 
listaj=as.data.frame(listaj)
objetoj=en_area
objetoj[,]=0
objetoj[nceldas]<-listaj[,1]
OBJETOj=mask(objetoj,colombia)
plot(OBJETOj, main="COMPLEMENTARIEDAD_JACK")
plot(colombia, add=T)
raster::writeRaster(OBJETOj,"INDICE_COMPLEMENTARIEDAD_JACK.tif", overwrite=TRUE)


dev.off()

############# 4.  GAP SELECTION INDEX ####################
pdf("4 INDICE_final.pdf ") ##comienza la grafica tipo pdf

par(mfrow=c(1,1))

# etsandarizacion de valores

DENSIDAD=ESTANDARIZACION(contador2)
AMBIENTAL=ESTANDARIZACION(d)
COMPLEMENTARIEDAD_BOOT=ESTANDARIZACION(OBJETO)
COMPLEMENTARIEDAD_JACK=ESTANDARIZACION(OBJETOj)
  
raster::writeRaster(DENSIDAD,"INDICE_DENSIDAD_EST.tif", overwrite=TRUE)
raster::writeRaster(AMBIENTAL,"INDICE_AMBIENTAL_EST.tif", overwrite=TRUE)
raster::writeRaster(COMPLEMENTARIEDAD_BOOT,"INDICE_COMPLEMENTARIEDAD_BOOT_EST.tif", overwrite=TRUE)
raster::writeRaster(COMPLEMENTARIEDAD_JACK,"INDICE_COMPLEMENTARIEDAD_JACK_EST.tif", overwrite=TRUE)

#  GSI

GSI_BOOT=(3-DENSIDAD-AMBIENTAL-COMPLEMENTARIEDAD_BOOT)/3
GSI_JACK=(3-DENSIDAD-AMBIENTAL-COMPLEMENTARIEDAD_JACK)/3

# invasora
# GSI_INVASORAS=mask((2-DENSIDAD-AMBIENTAL)/2,colombia)
# plot(GSI_INVASORAS, main="GAP SELCTION INDEX INVASORAS")

plot(colombia, main="UBICACION DE LOS REGISTROS")
plot(DATOS2,add=T)
plot(GSI_BOOT, main="GAP SELCTION INDEX_BOOT")
plot(GSI_JACK, main="GAP SELCTION INDEX_JACK")

dev.off()

raster::writeRaster(GSI_BOOT,"INDICE_GSI_BOOT.tif", overwrite=TRUE)
raster::writeRaster(GSI_JACK,"INDICE_GSI_JACK.tif", overwrite=TRUE)



# #### METRICA MEJORA VACIOS ----------------------------------------------



#subira capas anteriores
ruta_anterior="~/GBIF/VACIOS/GSI"
setwd(ruta_anterior)

densidad_anterior=raster("INDICE_DENSIDAD_EST.tif")
ambiental_anterior=raster("INDICE_AMBIENTAL_EST.tif")
complementBoot_anterior=raster("INDICE_COMPLEMENTARIEDAD_BOOT_EST.tif")
complementJack_anterior=raster("INDICE_COMPLEMENTARIEDAD_JACK_EST.tif")
GSIJack_anterior=raster("INDICE_GSI_JACK.tif")
GSIBoot_anterior=raster("INDICE_GSI_BOOT.tif")

# Mejoras y plots

par(mfrow=c(1,3))
mejora_densidad=mask(DENSIDAD-densidad_anterior,colombia)
mejora_densidad[mejora_densidad<0]=0
length(which(mejora_densidad[,]>0))/length(which(!is.na(mejora_densidad[,])))
plot(mask(densidad_anterior,colombia),main="DENSIDAD ANTERIOR");plot(mask(DENSIDAD,colombia), main="DENSIDAD ACTUAL");plot(mejora_densidad,main="MEJORA DENSIDAD",sub="mejora en el 15% del pais")
# celdas nuevas
length(which((densidad_anterior-DENSIDAD)[,]<0))
plot(mask(densidad_anterior,colombia),main="DENSIDAD ANTERIOR");plot(mask(DENSIDAD,colombia), main="DENSIDAD ACTUAL", sub="2960 nuevas localidades");plot(mejora_densidad,main="MEJORA DENSIDAD",sub="mejora en el 15% del pais")


mejora_ambiental=mask(AMBIENTAL-ambiental_anterior,colombia)
mejora_ambiental[mejora_ambiental<0]=0
plot(mask(ambiental_anterior,colombia),main="AMBIENTAL ANTERIOR");plot(mask(AMBIENTAL,colombia), main="AMBIENTAL ACTUAL");plot(mejora_ambiental,main="MEJORA AMBIENTAL")


mejora_complJack=mask(COMPLEMENTARIEDAD_JACK-complementJack_anterior,colombia)
mejora_complJack[mejora_complJack<0]=0
plot(mask(complementJack_anterior,colombia),main="JACK ANTERIOR");plot(mask(COMPLEMENTARIEDAD_JACK,colombia), main="JACK ACTUAL");plot(mejora_complJack,main="MEJORA JACK")


mejora_complBoot=mask(COMPLEMENTARIEDAD_BOOT-complementBoot_anterior,colombia)
mejora_complBoot[mejora_complBoot<0]=0
plot(mask(complementBoot_anterior,colombia),main="BOOT ANTERIOR");plot(mask(COMPLEMENTARIEDAD_BOOT,colombia), main="BOOT ACTUAL");plot(mejora_complBoot,main="MEJORA BOOT")

mejora_GSIjack=mask(GSI_JACK-GSIJack_anterior,colombia)
mejora_GSIjack[mejora_GSIjack>0]=0
#metricas
min(GSIJack_anterior[,],na.rm=T); sd(GSIJack_anterior)
min(GSI_JACK[,],na.rm=T); sd(GSI_JACK)
#cuales mejoraron 
length(which(mejora_GSIjack[,]<0))
length(which(GSI_JACK[,]<GSIJack_anterior[,]))/length(which(!is.na(mejora_GSIjack[,])))
plot(mask(GSIJack_anterior,colombia),main="GSI JACK ANTERIOR");plot(mask(GSI_JACK,colombia), main="GSI JACK ACTUAL");plot(mejora_GSIjack,main="MEJORA GSI JACK",sub="3842 celdas (28.7%) mejoraron")


mejora_GSIBOOT=mask(GSI_BOOT-GSIBoot_anterior,colombia)
#mejora_GSIBOOT[mejora_GSIBOOT<0]=0
plot(mask(GSIBoot_anterior,colombia),main="GSI BOOT ANTERIOR");plot(mask(GSI_BOOT,colombia), main="GSI BOOT ACTUAL");plot(mejora_GSIBOOT,main="MEJORA GSI BOOT")


writeRaster(mejora_densidad, "mejora_densidad.tif", format="GTiff")
writeRaster(mejora_GSIjack, "mejora_GSIjack.tif", format="GTiff")