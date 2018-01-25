##########
#Funciones v0.2 rev240118
##########

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

getNamesUniprot <- function(x){
  
  protNames <- data.frame( ProtID = character(), ProtName = character(), ProtFunct = character())
  
  for (i in 1:dim(x)[1]){
    
    try(url <- paste ("http://www.uniprot.org/uniprot/", x[i,], ".xml", sep=""))
    data <- xmlParse(url)
    xml_data <- xmlToList(data)
    
    ProtName <- "Not available"
    
    try(ProtFunct<- as.character(xml_data$entry$comment$text[1]))
    try(ProtName <- as.character(xml_data$entry$gene$name$text[xml_data$entry$gene$name$.attrs=="primary"]))
    
    if (length(ProtFunct) == 0L){ProtFunct <- "not listed"}
    
    name <-data.frame( ProtID = as.character(x[i,]),ProtName=ProtName, ProtFunct=ProtFunct)
    protNames <- rbind(protNames, name)
    
    cat("     \r",round(i/(dim(Proteinas)[1])*100, 2), "%")
    
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

