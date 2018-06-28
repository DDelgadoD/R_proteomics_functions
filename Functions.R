##########
#Funciones v0.3 rev150218
##########

###
### Funcion: NFreq(@Lecturas:dataframe @column:integer $Lecturas:dataframe)    
###
NFreq <- function(Lecturas, column){
  Lecturas<- table(unlist(Lecturas[,column]))
  Lecturas<- as.data.frame(Lecturas)
  return(Lecturas)
}

## Cálculo del coeficiente de variación
var.coef <- function(x){
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)*100
}

## Llamar a una libreria o instalar si no se tiene el paquete
installIfNot <- function (pckgName){
  if(!(require(pckgName, character.only=TRUE))){
    install.packages(pckgName)
  }
}

getNamesUniprotBasic <- function(x){
 
  installIfNot("XML")
  installIfNot("httr")
  
  protNames <- data.frame( ProtID = character(), ProtName = character(), ProtFunct = character())
  
  for (i in 1:dim(x)[1]){
    
    try(url <- paste ("http://www.uniprot.org/uniprot/", x[i,], ".xml", sep=""))
    data <- xmlParse(rawToChar(GET(url)$content))
    xml_data <- xmlToList(data)
    
    ProtName <- "Not available"
    
    try(ProtFunct<- as.character(xml_data$entry$comment$text[1]))
    try(ProtName <- as.character(xml_data$entry$gene$name$text[xml_data$entry$gene$name$.attrs=="primary"]))
    
    if (length(ProtFunct) == 0L){ProtFunct <- "not listed"}
    
    name <-data.frame( ProtID = as.character(x[i,]),ProtName=ProtName, ProtFunct=ProtFunct)
    protNames <- rbind(protNames, name)
    
    cat("     \r",round(i/(dim(x)[1])*100, 2), "%")
    
  }
  
  return(protNames)
  
}

# Funci?n que crea una tabla de correlaci?n a partir de una matriz de datos.
# @ param: x - data frame
# $ return: cor.completa - data frame
#
correlationTable<- function(x){
  aa<-x[, c(-1,-2)]
  cor.completa<-cor(t(aa), method = "pearson")
  
  # Indentificamos las columnas y las filas a?adiendo adem?s la proteina a la que pertenece el p?ptido. 
  cor.completa<- as.data.frame(cor.completa)
  colnames(cor.completa)<- x$Identificador
  cor.completa <- cbind(x$Identificador, x$Proteina, cor.completa)
  colnames(cor.completa)[1:2] <- c("Identificador", "Proteina")
  cat("     Correlacion completa...\n")
  
  #Devolvemos cor.completa como retorno de la funci?n
  return(cor.completa)
}

# Funci?n que crea una tabla con los datos de las prote?nas altamente correlacionadas (con una correlaci?n mayor a 0.9 entre 
# m?s del 50% de los p?ptidos que aparecen en cor.completa. .
# @ param: x - matriz de datos
# $ return: hasName - data frame
#
highCorrelation <- function(x,Proteinas){
  
  High.corr <- data.frame( Proteina1 = character(), Proteina2 = character(), Percent =  numeric(), Peptidos= character())
  
  cat("     Progreso altamente correlacionadas:\n")
  
  for(i in 1:(dim(Proteinas)[1]-1)){
    
    aN <- as.character(Proteinas[i,2])
    a <- as.character(Proteinas[i,1])
    
    
    Prot1 <- subset(x, Proteina == a)
    v=i+1
    
    for(s in v:dim(Proteinas)[1]){
      ProtVS <- NULL
      bN <- as.character(Proteinas[s,2])
      b <- as.character(Proteinas[s,1])
      
      Prot2 <- subset(x, Proteina == b)
      
      for (t in 1:dim(Prot2)[1]){
        
        if (is.null(ProtVS)){
          
          ProtVS <- cbind(Prot1[,1:2], Prot1[as.character(Prot2$Identificador[t])])
          
        }else{
          
          ProtVS <- cbind(ProtVS, Prot1[as.character(Prot2$Identificador[t])])
        }
      }
      
      sum = 0
      
      # Obtenemos la cantidad de peptidos con correlaci?n mayor que 0.9 en valor absoluto
      sum <- sum(abs(ProtVS[,c(-1,-2)])>0.90)
      
      
      lim <- dim(ProtVS)[1]*((dim(ProtVS)[2])-2)
      
      if (((lim*0.9) < sum) & ((dim(ProtVS)[1]>1) & ((dim(ProtVS)[2])-2>1))) {
        h <- data.frame(Proteina1 = aN, Proteina2 = bN, Percent = sum/lim*100, PeptidosTotales=paste(dim(ProtVS)[1],"X",((dim(ProtVS)[2])-2)))
        High.corr <- rbind (High.corr, h)
      }
      
      ProtVS <- NULL
      
    }
    
    cat("     \r",round(i/(dim(Proteinas)[1])*100, 2), "%")
    
  }
  cat("     \r 100% \n")
  return(High.corr)
  
  
}

RelationsNamed <- function(Prot, Correlacio){
  
  t <-  data.frame( ProtID = character(), Relations = numeric(), ProteinasRelacionades = character())
  h  <- 2
  
  for(i in 1:nrow(Prot)){
    
    a<-subset(Correlacio, Proteina1 == as.character(Prot[i,h]))[,-1]
    b<-subset(Correlacio, Proteina2 == as.character(Prot[i,h]))[,-2]
    
    colnames(a)[1]<- "Prote?na"
    colnames(b)[1]<- "Prote?na"
    
    c<-as.data.frame(rbind(a,b))
    colnames(c)[1]<-"Prote?na"
    c<-c[with(c, order(-c[2])), ]
    
    f<-''
    
    for(d in 1:nrow(c)){
      
      if(d==nrow(c)){
        f <- paste(paste(f, as.character(c[d,1])),paste("(",round(c[d, 2],2),").",sep=''),sep='')
      }else{
        f<- paste(paste(f, as.character(c[d,1])),paste("(",round(c[d, 2],2),"),",sep=''),sep='') 
      }
      
    }
    
    t <- rbind(t, data.frame(ProtID = Prot$ProtName[i], Relations= nrow(c), ProtRelacionades = f))
    
  }
  
  t<- t[t[2]>0]
  t<-t[with(t, order(-t[2])), ]
  
  return(t)
}

Datos<- function(Name, Run, Prot, carpeta){
  
  cat("Realizando c?lculos...\n")
  
  #Llamamos a las funciones para correlacionar y extraer relaciones
  corr <- correlationTable(Run)
  highCorr <- highCorrelation(corr,Prot)
  relation <- RelationsNamed(Prot,highCorr)
  
  # Guardar las correlaciones en csv para evitar el crasheo de la librer?a de excel. 
  cat("     Generando csv de correlaci?n...\n")
  write.csv(corr , paste(paste(carpeta,"/CorrelacionCSV/Correlacion ",sep=""),Name,".csv", sep=""))
  cat("     Generando xlsx de altamente correlacionados...\n")
  write.xlsx(highCorr, paste(carpeta,"/AltaCorrelacion.xlsx",sep=""), sheetName=Name, append=TRUE, row.names=FALSE)
  cat("     Generando xlsx de Relaciones con nombre de proteina...\n")
  write.xlsx(relation, paste(carpeta,"/Relaciones.xlsx",sep=""), sheetName=Name, append=TRUE, row.names=FALSE)
  
  #Creamos un nuevo objeto para el retorno
  cat("     Generando lista de datos...\n")
  listar <- list(corr, highCorr,relation)
  return(listar)
  
}


# Nombre: getNamesUniprot
# Funci?n que devuelve un data.frame con el nombre y la funci?n de una prote?na humana. Requiere connexi?n a internet.
#
# @param: x: data.frame de una sola columna donde se indique el n?mero de proteina humana que deseamos buscar
# @return: protNames: data.frame de 3 columnas donde se indica el n?mero de prote?na humana, el nombre y su principal funci?n

getNamesUniprot <- function(x){
  
  # Preparamos el data.frame para alojar los datos de salida
  protNames <- data.frame( ProtID = character(), ProtName = character(), GO_C = character(), GO_F = character(), GO_P = character(), ProtFunct=character())
  
  # Recorremos el data.frame de entrada para obtener los datos
  for (i in 1:dim(x)[1]){
    # Probamos la conexi?n
    try(url <- paste ("http://www.uniprot.org/uniprot/", x[i,], ".xml", sep=""))
    # Obtenemos los datos
    data <- xmlInternalTreeParse(url, useInternal=TRUE)
    
    
    # Los convertimos en lista para acceder a los nodos
    xml_data <- xmlToList(data)
    # Rellenamos el nombre por si no tiene correspondencia
    ProtName <- "Not available"
    #Probamos a asignar los datos
    try(ProtName <- as.character(xml_data$entry$gene$name$text[xml_data$entry$gene$name$.attrs=="primary"]))
    try(ProtFunct<- as.character(xml_data$entry$comment$text[1]))
    
    # Si el string funci?n tiene largo 0 indicamos que no hay funci?n listada
    if (length(ProtFunct) == 0L){
      ProtFunct <- "not listed"
    }
    
    Pr <- as.list(xml_data$entry)
    Pr <- Pr[names(Pr)=="dbReference"]
    Pr <- Pr[sapply(Pr, function(x) any(unlist(x) == "GO"))]
    Pr <- as.data.frame(sapply(Pr,function(x) unlist(x$property["value"])))
    
    c<-paste(gsub("C:","",Pr[grep("^C", Pr[,1]),]), collapse = ", ")
    f<-paste(gsub("F:","",Pr[grep("^F", Pr[,1]),]), collapse = ", ")
    p<-paste(gsub("P:","",Pr[grep("^P", Pr[,1]),]), collapse = ", ")
    
    # Creamos una fila 
    name <-data.frame( ProtID = as.character(x[i,]),ProtName=ProtName, GO_C=c,GO_F=f,GO_P=p,ProtFunct=ProtFunct)
    # A?adimos la fila a el data.frame que hemos creado
    protNames <- rbind(protNames, name)
    
    # Mostramos un mensaje de progreso
    cat("     \r",round(i/(dim(x)[1])*100, 2), "%   ")
    
  }
  # Devolvemos el data.frame finalizado
  return(protNames)
  
}

# Nombre: binder
# Funci?n que devuelve un data.frame con una sola columna que auna los valores de un data.frame de dos columnas
#
# @param: nX: data.frame de dos columnas 
# @return: bind: data.frame de 1 columna

binder<-function(nX){
  a <- as.data.frame(nX[,1])
  b <- as.data.frame(nX[,2])
  colnames(a)[1]<- "Prote?na"
  colnames(b)[1]<- "Prote?na"
  bind<-as.data.frame(rbind(a,b))
  
  return(bind)
}

# Aplicar filtro

filtering <- function(x, filter, where){
  
  if (where=="c"){
    # Filtrar en GO C
    f<-ProteinGO[grep(filter, ProteinGO[,3]),]
  }else if (where=="f"){
    # Filtrar en GO F
    f<-ProteinGO[grep(filter, ProteinGO[,4]),]
  } else if (where == "p"){
    # Filtrar en GO P
    f<-ProteinGO[grep(filter, ProteinGO[,5]),]
  }else{
    f <- "Use a character for 'filter' and 'where' and lowercase for 'where'. Ex: filtering(data.frame,'axon','c') "
  }
  
  return(f)
}

# Comprobar la normalidad
Normal<- function(E){
  LoQueEsNormal <- data.frame(Proteina=character(), S.W=character(), pvalor=numeric(), A.D= character(), porcentaje=numeric(), stringsAsFactors=FALSE) #Creamos una tabla vac?a
  
  for (i in 1:ncol(E)){
    # Realizamos el test
    Alpha <- shapiro.test(as.array(E[,i]))
    P <- Alpha$p.value
    ADT <- ad.test(E[,i])
    AD <- ADT$p.value
    
    if (P>=0.05){
      Wii <- data.frame(Proteina = colnames(E)[i], S.W= "Sí", pvalor=round(P, 4), A.D= "", porcentaje = AD*100)  
    }else{
      Wii <- data.frame(Proteina = colnames(E)[i], S.W = "No", pvalor=round(P, 4), A.D= "", porcentaje = AD*100) 
    }
    
    if(AD > 0.5){
      Wii$A.D = "Sí"
    }else{
      Wii$A.D = "No"
    }
    
    LoQueEsNormal <- rbind(LoQueEsNormal, Wii)
  }
  
  LoQueEsNormal
}

# Gràficos de densidad
densitedVar <- function(E){
  opar = par() 
  pareo()
  
  ytpos = (apply(E,2,min)+3*apply(E,2,max))/4
  cn = colnames(E)
  
  for(i in 1:ncol(E)){
    densityPlot(E[,i], ylab= colnames(E)[i])
    text(-1.5, ytpos[i] ,cn[i])
    
  } 
  par(opar)
}

## Adaptamos la salida gràfica para que muestre 16 gràficos
pareo <-function(){
  
  par(mfrow=c(4,4))
  par(mar=c(0.5,0.5,0.5,0.5))
  par(oma=c(1,1,1,1))
}

#Normalización por Johnson
JohnM <- function(E){
  
  for (i in 1:ncol(E)){
    x_johnson <- RE.Johnson(E[,i])
    E[,i] <- as.data.frame(x_johnson$transformed)
  }
  return(E) 
}

correlate <- function(x, met){
  
  mCor <- data.frame(Prot= character(), sign= numeric(), corr= numeric(), stringsAsFactors = FALSE )
  
  for(i in 2:ncol(x)-1){
    v <- cor.test(x[,i], x[,1], use="complete.obs", method=met)
    a <- v$estimate
    b <- v$p.value
    dat <- data.frame(Prot= colnames(x)[i], sign= round(b[1], digits=3), corr= round(a[1], digits=3),corrA= abs(round(a[1], digits=3)), stringsAsFactors = FALSE)
    mCor <-rbind(mCor, dat)
  }
  
  return(mCor)
}

pHCorr <- function(x, met){
  xF <- x
  xF$FDR = round(p.adjust(x$sign, method=met ), digits=5 )
  return(xF)
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n){
  substr(x, 1, n)
}
