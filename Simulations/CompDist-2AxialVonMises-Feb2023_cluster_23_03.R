# Code to test which tests are more powerful in detecting differences between distributions 

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied: current 'i'", call.=FALSE)
}

message("Given argument: ", args[1])

library(circular)
library(compiler)
library(NPCirc)
library(foreach)
library(doParallel)
library(permute) #for shuffle 
options("scipen" = 10)
enableJIT(3)

#install TwoCircles package from git hub - only needed the first time 
#install.packages(c("gRbase", "circular", "devtools")) #only first time
#devtools::install_github("SMAC-Group/TwoCircles") #only first time

#library(TwoCircles) #load two Circles package 

set.seed(1)

#possible improvements run on multiple nodes: https://stackoverflow.com/questions/64984094/foreach-iterators-on-multiple-nodes 

#library(TwoCircles) #load two Circles package 
source("/gpfs/soma_fs/home/malkemper/R/Scripts/crit_val.R")
source("/gpfs/soma_fs/home/malkemper/R/Scripts/data.R")
source("/gpfs/soma_fs/home/malkemper/R/Scripts/tests.R")
source("/gpfs/soma_fs/home/malkemper/R/Scripts/Functions_Feb2023_for_publication.R")
options("scipen" = 10)
enableJIT(3)

#setting sample size for dist 1 - in two steps to reduce typing effort
nvseq=c(10,20,50,20,10,30,50)*2
nv1=c(nvseq,nvseq,nvseq) 

#setting sample size for dist 2
nvseq2=c(10,20,50,30,50,20,10)*2
nv2=c(nvseq2,nvseq2,nvseq2)

#combine sample sizes in matrix to make calling in the loop practical
n_matrix = cbind(nv1,nv2)

#preparing certain vectors of concentration parameters and direction - all 7 elements long (they need to be same length) 

kseq=c(0,0.25,0.5,1,2,4,8) #standard kappa sequence for plot 

k0=rep(0,7)
k2=rep(2,7)

mseq= c(rad(0),rad(30),rad(60),rad(90),rad(120),rad(150),rad(180))/2
m0=rep(0,7)

#this would be needed of lengths differ, at the present it doesn't matter which l is used 
lk=length(kseq) 
lm=length(mseq)
l=lk # this one is in this script used for the inner loop

#preparing the sequences in a list so they can be called in inner loop 
k_list_1=c(rep(list(kseq),7),c(rep(list(k0),7)),c(rep(list(k2),7)))
k_list_2=c(rep(list(kseq),7),rep(list(kseq),7),rep(list(k2),7))
m_list_1=c(rep(list(m0),7),rep(list(m0),7),rep(list(m0),7))
m_list_2=c(rep(list(m0),7),rep(list(m0),7),rep(list(mseq),7))


#defining the distribution type, there are always 5 combinations of sample sizes 
#therefore I prepared "chunks" of distribution and then combined them, decreasing typing effort 
type_vm=rep("vm",7) 
type_matrix_1 = c(type_vm,type_vm,type_vm)  
type_matrix_2 = c(type_vm,type_vm,type_vm)  

#again in chunks of 5: naming the file that is created (see text for write.csv)
Type_1_errorn=rep("Type_1_error",7)
Powern=rep("Power",7)

#again combining the chunks 
Name=c(Type_1_errorn,Powern,Powern)

########### set number of iterations, Sample size mean and kappa (or a sequence of any of that)
rans = 10000#iterations for  distribution

#starting the outer loop, it loops through the sample size vector, taking the respective vectors (of length 7) from the respective list
#the outer loop is stored one csv file, this means one panel in the plot. 
#for (i in 1:length(nv1)) {
for (i in args[1]:args[1]) {
message("loop i: ", i)

#grabing the correct vectros for the simulation from the list
nv1 = n_matrix[,1][i]
nv2=n_matrix[,2][i]
m1seq=m_list_1[[i]]
m2seq=m_list_2[[i]]
k1seq=k_list_1[[i]]
k2seq=k_list_2[[i]]


typeDist1 =type_matrix_1[i]
typeDist2 = type_matrix_2[i]


# make AB data.frame for the tests where the data needs to be structured in the same column
#The data frame consists of A and Bs arranged underneath each other indicating the groups (like a "treatment" column) 
Group1 <-(as.matrix( rep((as.matrix((strrep(c("A"), 1)))),nv1),dim=c(nv1, 1)))
Group2 <- (as.matrix( rep((as.matrix((strrep(c("B"), 1)))),nv2),dim=c(nv2, 1)))
GroupsAB<-as.data.frame(c(Group1[,1], Group2[,1]))

#setup parallel back-end to use many processors
cores=detectCores()
no_cores <- parallelly::availableCores(omit = 1)
message("No cores: ",no_cores)
no_workers <- min(c(no_cores,l))
message("No workers: ",no_workers)
#q() for testing
cl <- parallel::makeCluster(no_workers, setup_strategy = "sequential") #used "sequential" otherwise it crashed for me  														 
registerDoParallel(cl)

#start the inner loop. Here I use a foreach instead of for. It splits the jobs on the number of cores set up. 
#Parallelization is done with "%dopar%". The packages used in the foreach loop must be passed to the loop.
#they are not grabbed automatically from the environment.
#the results of each of the loops are combined (with rbind) in the file Power.sample
Power.sample<-foreach(e = 1:l,.combine=rbind,.packages=c("circular",
                                                "NPCirc","permute")) %dopar% {

source("/gpfs/soma_fs/home/malkemper/R/Scripts/crit_val.R")
source("/gpfs/soma_fs/home/malkemper/R/Scripts/data.R")
source("/gpfs/soma_fs/home/malkemper/R/Scripts/tests.R")

#take correct values for mean direction as well as the concentration parameter from the vector
 m1 = m1seq[e] 
  m2 = m2seq[e]                                       
  k1 = k1seq[e]
  k2 = k2seq[e]

  #First distribution generated rans times, the "units=" could be deleted. Its just here to remind me what it is, but doesn't do anything.
  
  Distribution1 <- as.data.frame(sapply(1:rans, function(x) rcircmix(nv1,model = NULL, 
                                                                     dist = c(typeDist1,typeDist1),
                                                                     param=list(p=c(0.5,0.5), 
                                                                                mu=c(m1,m1+pi), 
                                                                                con=c(k1,k1)))))
  
  ###### set up the 2nd distribution 
  Distribution2 <- as.data.frame(sapply(1:rans, function(x) rcircmix(nv2,model = NULL, 
                                                                     dist = c(typeDist2,typeDist2),
                                                                     param=list(p=c(0.5,0.5), 
                                                                                mu=c(m2,m2+pi), 
                                                                                con=c(k2,k2)))))
  
  ###### underneath we apply each of the tests used on the generated distribution 
  #using mapply (if the groups are arranged in separate column) and apply (if they are arranged underneath each other)
  
  ##watson2SampleTest
  watson1 = mapply(function(x,y) {watson.two.test(x,y)},Distribution1,Distribution2)
  watson1 <- array(as.matrix(unlist(watson1)), dim=c(5, rans))
  watson1<- watson1[1,]  # critical valueCritical Value: 0.187 
  
  #Watson U2 permutation version
  watson1.perm= sapply(1:rans,function(x) watson.two.test.perm(Distribution1[,x],Distribution2[,x]))
  
  #ART
  G.p= sapply(1:rans,function(x) G_test(Distribution1[,x],Distribution2[,x]))
 
   #"^.*?: " beginning of the string to everything before the ": ", the latter gets replaced with "" 
  Rao.p=  sapply(1:rans,function(x) 
    as.numeric(gsub("^.*?: ","",
                    capture.output(circular_test(
                      Distribution1[,x], Distribution2[,x],
                      type="mc",B = 999, test = "rao"))[6])))
  
  #####preparing DistGroups for watson wheeler
  DistGroups<-mapply(function(x,y) {c(x, y)},Distribution1,Distribution2)
  
 #Watson wheeler test
  watson.wheelerRES <- apply(DistGroups,2,FUN=function(x) watson.wheeler.test(x,GroupsAB[,1]))
  watson.wheelerP <- sapply(c(1:rans),function(x) watson.wheelerRES[[x]][["p.value"]])
  
  #Watson wheeler permutation test 
  watson.wheelerP.perm <-sapply(1:rans,function(x) watson.wheeler.test.perm(Distribution1[,x],Distribution2[,x]))
  
  #create a temp data.frame to store results
  names_col<- c('WatsonTwoSamp','WatsonTwoSamp.perm','Watson.Wheeler','Watson.Wheeler.perm',"G.test","Rao_spacing_freq" )
  
  temp <- data.frame(matrix(nrow = 1, ncol = length(names_col)))
  names(temp) <- names_col
  
  #####Calculate type I error/power in temp file 
  temp$WatsonTwoSamp<- (sum(1*(watson1>0.187)))/rans
  temp$WatsonTwoSamp.perm<- (sum(1*(watson1.perm<0.05)))/rans
  temp$Watson.Wheeler<- (sum(1*(watson.wheelerP<0.05)))/rans
  temp$Watson.Wheeler.perm<- (sum(1*(watson.wheelerP.perm<0.05)))/rans
  temp$G.test<- (sum(1*(G.p<0.05)))/rans
  temp$Rao_spacing_freq<- (sum(1*(Rao.p<0.05)))/rans
  temp # this is stored in  Power.sample
} #e
stopCluster(cl)


#the x-axis is calculated from the concentration of the second distribution, if the maximum is larger the minimum number
#otherwise it takes the mean direction vector
Power.sample$x=if(max(k2seq)>min(k2seq)) k2seq else deg(m2seq) 

Power.sample

#prepare the text of the csv. file 
text = paste("Axial_",Name[i],"_n1_",nv1,"_n2_",nv2,"_its_",rans,"_Dist1_",typeDist1,"_Dist2_",typeDist2,
             "_m1=",deg(min(m1seq)),"-",deg(max(m1seq)),"_m2_",deg(min(m2seq)),"-",deg(max(m2seq)),
             "_k1=",min(k1seq),"-",max(k1seq),"_k2_",min(k2seq),"-",max(k2seq),".csv",sep="")

#write csv file in the working directory
write.csv(Power.sample, file = text,row.names=FALSE )

}#i

