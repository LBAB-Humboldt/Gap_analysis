newSpRate <- function(data, layer = NULL, layerField = NULL,  simulations = 10000, outDir = FALSE){
  if (outDir != FALSE) dir.create(outDir, recursive = TRUE)
  level <- na.omit(unique(data$locality))
  summaryByLevel <- NULL
  i <- 1
  for (i in 1:length(level)){
    data.i <- data$species[which(data$locality == level[i])]
    uniqueSp <- length(unique(data.i))
    propSp <- table(data.i)/length(data.i)
    shannon <- -sum(propSp * log(propSp))
    simpson <- sum(propSp ^ 2)
    
    sim.i <- matrix(0, simulations, 2)
    raref.i <- matrix(NA, nrow = simulations, ncol = length(data.i))
    
    for (j in 1:simulations){
      data.j <- data.i[sample(1:length(data.i), length(data.i), replace = FALSE)]
      dup.j <- !duplicated(data.j)
      acum.j <- 1:sum(dup.j)
      index.j <- which(dup.j)
      lm.j <- lm(acum.j ~ index.j)
      sim.i[j, ] <- c(summary(lm.j)$coefficients[2, 1:2])
      raref.i[j, dup.j] <- acum.j
    }
    
    raref <- colMeans(raref.i, na.rm = TRUE)
    dfRaref <- data.frame(index = 1:length(raref), raref = raref)
    coefRaref <- summary(lm(raref ~ index, data = dfRaref))$coefficients[2, 1]
    coefRarefOrigin <- summary(lm(raref ~ 0 + index, data = dfRaref))$coefficients[1, 1]
    #     boxplot(data.frame(coef = sim.i[, 1]))
    #     abline(h = coefRaref)
    summaryByLevel <- rbind(summaryByLevel, c(level[i], coefRaref, coefRarefOrigin, mean(sim.i[, 1]),uniqueSp, shannon, simpson))
    
    if (outDir != FALSE){
      tryCatch(write.csv(raref, paste0(outDir,'/rarefCurve_', level[i], '.csv'), row.names = FALSE))
      colnames(sim.i) <- c('slope', 'stdError')
      tryCatch(write.csv(sim.i, paste0(outDir, '/coefTable_', level[i], '.csv'), row.names = FALSE))
    }
    cat(i, '-', length(level), ' || ')
  }
  colnames(summaryByLevel) <- c('level', 'slopeMeanCurve', 'slopeMeanCurveOrig', 'meanSlopeCurves', 'sObs', 'shannon', 'simpson')
  
  if(!is.null(layer) & !is.null(layerField)){
    layer@data$sortID <- 1:nrow(layer@data)
    layer@data <- merge(layer@data, as.data.frame(summaryByLevel), by.x = layerField, by.y = 'level', all.x = TRUE, all.y = FALSE)
    layer@data <- layer@data[order(layer@data$sortID), ]
    writeOGR(layer, dsn = outDir, layer = 'layer')
  }
  
  return(summaryByLevel)
}
