

#######function for calculating G value for the geodesic distance Test (ART)
Gvalue <- function(s1,s2,p,q){
  total <- 0 
  for(i in 1:p){
    for(j in 1:q){
      temp <-  (pi - abs(pi-abs(s1[i]-s2[j])))
      total <- total+temp}}
  return(total)}

#######function for the actual test using the  G value for the geodesic distance Test (ART)
G_test <-function(samp1,samp2){
  NR <- 999
  n1 <- length(samp1)
  n2 <- length(samp2)
  Gstat <- Gvalue(samp1,samp2,n1,n2)
  nxtrm <-1 
  n <- n1+n2
  combsample <- c(samp1,samp2)
  for(r in 1:NR){
    perm <- shuffle(1:n)
    randsamp1 <- combsample[perm[1:n1]]
    randsamp2 <- combsample[perm[(n1+1):(n1+n2)]] 
   
    Grand <- Gvalue(randsamp1,randsamp2,n1,n2)
    if(Grand >= Gstat){nxtrm<-nxtrm+1}
    
  }
  return(nxtrm/(NR+1))
  
}

