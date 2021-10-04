##------------------------------------------------------------------------------------------------
## fonction pour supprimer des lignes sans au moins une mesure des outcomes d'intérêt
fifi <- function(x){
  d<-0 
  len <-length(na.omit(x))
  if(len == 0){
    d <- 1
  }
  return(d)
}
#=====================================================================================
one.outcome.at.least <- function(data, nom.subject, outcomes){
  cl <- match.call()
  colnames <- colnames(data)
  data2 <- data[, c(nom.subject, outcomes)]
  p <- length(outcomes)
  if( p < 2){
    data2$x_NA <- is.na(data2[,c(outcomes)])
    data2$x_NA <- as.numeric(data2$x_NA)
  }
  if(p>=2){
    data2$x_NA <- apply(data2[,c(-which(colnames %in% nom.subject))],MARGIN = 1, FUN = fifi)
  }
  data <- data[which(data2$x_NA ==0),]
  return(data)
}
#=====================================================================================
# fonction pour discr?tiser
f_TempsDiscr <- function(Time, Delta){
  AxeT <-seq(from =min(Time), to = max(Time), by = Delta)
  Time_d <- rep(0,length(Time))
  
  for(j in 1:length(Time)){
    i <-0
    Time_d[j] <- NA
    if(!is.na(Time[j])){
      Time_d[j] <- 0
      while(Time[j] > Delta*i){
        i <- i+1
      }
      Time_d[j] <- i*Delta
    }
  }
  return(Time_d)
}

#=====================================================================================
# formatage des données brutes comme 3C pour retenir les outcomes, les covariables,
# temps discrets d'observations
# retourne un data frame
DataCreate <-function(data, outcomes, predictors = NULL, subject, Time, age0 = 0, age0.from.data = FALSE, Delta){
  
  cl <- match.call()
  colnames<-colnames(data)
  if(!(subject%in%colnames))stop("Subject should be in colnames")
  
  # v?rifier que les Time sont dans le data.frame....
  
  # veriifer que length(outcomes) == length(Time)...
  if(length(outcomes) != length(Time)){
    stop("difference between number of outcomes and number of Time variable")
  }
  ##  pré-traitement du data frame
  data <- one.outcome.at.least(data, nom.subject= subject, outcomes = outcomes)
  
  K <- length(outcomes)
  T_k <-NULL
  nameTks <- NULL
  for(k in 1:K) {
    nameTk <- paste(Time[k], k, sep = "_")
    nameTks <- c(nameTks,nameTk)
    #Discr?tisation du temps de suivi
    T_k <- cbind(T_k,assign(nameTk, f_TempsDiscr(data[,Time[k]], Delta)))
  }
  T_k <- data.frame(T_k)
  colnames(T_k) <-nameTks
  data2 <-data.frame(cbind(data[,c(subject,outcomes, predictors)],T_k))
  #calcul du T max: Temps d'observation maximal
  # Tmax <- max((data2[,nameTks]),na.rm = TRUE)
  
  # merge
  data3 <- na.omit(data2[,c(subject, outcomes[1], predictors, nameTks[1])])
  #   data3 <- data2[,c(subject, outcomes[1], predictors, nameTks[1])]
  #   data3 <- data2[,c(subject, outcomes[1], predictors, nameTks[1])]
  colname <- colnames(data3)
  colname[which(colname==nameTks[1])] <- paste(Time[1],"d", sep = "_")
  colnames(data3) <- colname
  
  if(K>1){
    for(k in 2:K){
      data_k <- na.omit(data2[,c(subject, outcomes[k], predictors, nameTks[k])])
      #             data_k <- data2[,c(subject, outcomes[k], predictors, nameTks[k])]
      # changement du nom de la colonne T_i en T pour faire le merging
      colname <- colnames(data_k)
      colname[which(colname==nameTks[k])] <- paste(Time[1],"d", sep = "_")
      colnames(data_k) <- colname
      
      data3 <- merge(data3,data_k, by=c(subject,predictors,paste(Time[1],"d", sep = "_")), all.x = TRUE, all.y = TRUE)
      data3 <- data3[order(data3[,subject], data3[,c(paste(Time[1],"d", sep = "_"))]), ]
    }
  }
  return(data3)
}