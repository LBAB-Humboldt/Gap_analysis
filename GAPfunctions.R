# GAP functions


# C. FUNCION COMPLEMENTARIEDAD BASE DE DATOS
#lista <- c(lista, COMPLEMENTARIEDAD_BOOTS(datos, k))
#cellID  <- 28576
#sppList <- x

compBoot <- function(sppList){
  Sobs <- length(unique(sppList))
  Sexp <- Sobs + sum((1 - (table(sppList) / length(sppList))) ** length(sppList))
  return(Sexp) 
}

compJack <- function(sppList, nSamples, jackOrder = 1){
  Sobs <- length(unique(sppList))
  STable <- table(sppList)
  L <- length(which(STable == 1))
  J1 <- L * ((nSamples - 1)/nSamples)
  if (jackOrder == 1){
    Sexp <- Sobs + J1
    return(Sexp) 
  } else  if (jackOrder == 2){
    Sexp <- Sobs + ((L * (nSamples + nSamples - 3))/ nSamples) -
       ((L * (nSamples - 1) ** 2)/(nSamples * (nSamples - 1)))
    return(Sexp)
  }
}

list2Matrix <- function(inList, nRow = NULL, nCol = NULL, colNames = NULL, rowNames = NULL){
  if(!is.null(nCol)) {
    outMatrix <- matrix(unlist(inList), ncol = nCol, byrow = TRUE)
  }
  if(!is.null(nRow)) {
    outMatrix <- matrix(unlist(inList), nrow = nRow, byrow = TRUE)
  }
  if(!is.null(colNames)) {colnames(outMatrix) <- colNames}
  if(!is.null(rowNames)) {rownames(outMatrix) <- rowNames}
  outMatrix <- as.data.frame(outMatrix)
}

#simulations <- 100; nObs <- 50
compRar <- function(sppList, simulations, nObs){
  simMatrix <- matrix(0, ncol = nObs, nrow = simulations)
  for(s in 1:simulations){
    sppRaref <- sppList[sample(1:length(sppList), nObs, replace = FALSE)]
    sppSimRich <- !duplicated(sppRaref)
    acumCurve <- 1:sum(sppSimRich)
    index.j <- which(sppSimRich)
    if(max(acumCurve) != nObs){  
      filledCurve <- fillCurve(acumCurve, index.j, nObs)
      simMatrix[s, ] <- filledCurve
    } else {
      simMatrix[s, ] <- acumCurve
    }
  }
  return(data.frame(Smean = colMeans(simMatrix, na.rm = TRUE), 
                    Svar = apply(simMatrix, 2, function (x) var(x, na.rm = TRUE))
                    )
         )
}

#acum <-acumCurve; index <- index.j; lenData <- 25
fillCurve <- function(acum, index, lenData){
  fill <- rep(NA, lenData)
  fill[index] <- acum
  
  compIndex <- data.frame(cbind(index, c(index.j[-1], lenData), acum))
  compIndex <- compIndex[(compIndex[, 1] - compIndex[, 2]) < -1, ]
  if (nrow(compIndex) >0 ){
    for (c in 1:nrow(compIndex)){
      fill[(compIndex[c, 1] + 1):(compIndex[c, 2] - 1)] <- compIndex[c, 3]
    }
  }
  return(fill)
}

compJack <- function(sppList, jackOrder){
  Sobs <- length(unique(spL)) #j1.S
  sppFreq <- table(spL) #j1.1
  n <- table(table(data.j))
  #   m <- data.frame(j = names(n), n_j = n[])
  #   n <- apply(m, 2, as.integer)
  L <- length(sppFreq[sppFreq==1])
  j1.m=i1
  jack1=Sobs+L*((j1.m-1)/j1.m)

}




#sppList <- spListByCell$nombre_aceptado
#indexID <- spListByCell$celdas
#x <- sppList[indexID == 28430] # Good/28430 Bad/532
library(SPECIES)
richEst <- function(sppList, indexID){
  sEstimation <- tapply(sppList, INDEX = indexID, 
                  FUN = function(x){
                    
                    n <- table(table(x))
                    m <- data.frame(j = names(n), n_j = n[])
                    n <- apply(m, 2, as.integer)
                    
                    J <- tryCatch(jackknife(n, 1), error = function(e)  rep(NA, 5))
                    C1 <- tryCatch(chao1984(n), error = function(e)  rep(NA, 4))
                    C2 <- tryCatch(ChaoLee1992(n), error = function(e)  rep(NA, 8))
                    C3 <- tryCatch(ChaoBunge(n), error = function(e)  rep(NA, 4))
                    
                    return(c(length(x), length(unique(x)), compBoot(x),
                      unlist(J), unlist(C1), unlist(C2), unlist(C3)))
                  })
  sEstimation <- list2Matrix(sEstimation, nCol = 24,
                       colNames = c('Nobs', 'Sobs', 'Boot', 'JackknifeOrder', 'JNhat', 'JSE', 'JCI1', 'JCI2', 
                                    'ChaoNhat', 'ChaoSE', 'ChaoCI1', 'Chao1CI2', 'ChaoLNhat1', 'ChaoLNhat2',
                                    'ChaoLSE1', 'ChaoLSE2', 'ChaoLCI1', 'ChaoLCI2', 'ChaoLCI3', 'ChaoLCI4',
                                    'ChaoBNhat', 'ChaoBSE', 'ChaoBCI1', 'ChaoBCI2'),
                       rowNames = names(sEstimation)
                       )
  }


# C. FUNCION COMPLEMENTARIEDAD BASE DE DATOS
#lista <- c(lista, COMPLEMENTARIEDAD_BOOTS(datos, k))
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

normalize01 <- function(x, outDir = NULL){
  xNorm <- (x - x@data@min)/(x@data@max - x@data@min)
  if (!is.null(outDir)){
    writeRaster(xNorm, outDir, overwrite=TRUE)
  }  
  return(xNorm)
}

# A. FUNCION SESGO FACTORES
# biasLayer <- AP; rasterMask <- grilla
#biasLayer <- ambientales[[k]]; layerName <- names(ambientales)[k]; outDir <- outPlotDir; doplot <- TRUE
#biasLayer <- urbano; layerName <- "CASCOS_URBANOS"; outDir <- outPlotDir; doplot <- TRUE
#BIAS(, grilla, "CASCOS_URBANOS", outPlotDir)
BIAS <- function(biasLayer, rasterMask, layerName, outDir, doplot = TRUE){
  if (class(biasLayer) != 'RasterLayer'){
    par(mfrow = c(1, 3))
    raster_shp <- rasterize(biasLayer, rasterMask, field = 1)
    distancia <- distance(raster_shp)
    en_area <- mask(distancia, colombia)  
    if(doplot){
      plot(en_area, main = paste("Distance to\n", layerName))
    }
    var <- na.omit(en_area[])
    fishers <- classIntervals(var, 4, style = "fisher")
    clasificado <- cut(en_area, breaks = fishers$brks, right = FALSE)
  } else {
    par(mfrow = c(1, 2))
    var <- na.omit(biasLayer[])
    fishers <- classIntervals(var, 4, style = "fisher")
    clasificado <- cut(biasLayer, breaks = fishers$brks, right = FALSE)
  }
  
  if(doplot){
    plot(clasificado, main = paste("Levels\n", layerName))
  }
  
  if (class(biasLayer) == 'RasterLayer'){
    nLayer <- length(which(localidades[, ]>=1))
    valores <- clasificado * localidades  
  } else {
    valores <- extract(clasificado, DATOS2)
    nLayer <- length(which(!is.na(clasificado[])))
  }
  
  sesgos <- sapply(1:4, function(x){
    nd <- length(which(clasificado[] == x))/nLayer
    pd <- length(which(clasificado[] == x))/length(which(clasificado[] >= 1))
    BIAS <- (nd - pd * N)/(sqrt(pd * (1 - pd) * N))
  })
  
  matriz <- data.frame(ID = c(1:4), val = sesgos)
  niveles <- reclassify(clasificado, matriz, filename = paste0(outDir, "/bias_", layerName, ".tif"), overwrite = TRUE)
  names(niveles) <- layerName
  if (doplot){
    plot(niveles, main = paste("Bias\n", layerName))
    if (class(biasLayer) != 'RasterLayer'){
      plot(biasLayer, add = TRUE, border = "grey", lwd = 1)
    }
  }
  return(list(biasLayer = niveles, sesgos =  sesgos))
}

#biasLayer <- ambientales[[k]]
BIAS_AMBIENTAL <- function(biasLayer){   
  var <- na.omit(biasLayer[])
  fishers <- classIntervals(var, 4, style = "fisher")
  clasificado <- cut(biasLayer, breaks = fishers$brks, right = F)
  N <- length(which(localidades[, ]>=1))
  valores <- clasificado * localidades
  sesgos <- sapply(1:4, function(x){
    nd <- length(which(valores[] == x))
    pd <- length(which(clasificado[] == x))/length(which(clasificado[] >= 1))
    BIAS <- (nd-(pd*N))/(sqrt(pd*(1-pd)*N))
  })
  matriz <- data.frame(ID = c(1:4), val = sesgos)
  niveles <- reclassify(clasificado, matriz, filename = paste0(outDir, "/ambiental_", layerName, ".tif"), overwrite = TRUE)
  names(niveles) <- layerName
  if (doplot){
    plot(niveles, main = layerName)
    plot(shp, add = TRUE, border = "grey", lwd = 1)
  } 
  return(niveles)
}
