library(dplyr)
library(plyr)
?getwd()
setwd()
Col8=read.table('Downloads/Shared/EstonianPharmacoVariantsCol8.vcf',header=F,stringsAsFactors = F)
new=strsplit(Col8$V8,'\\|')
MAF_list=ldply(new,function(x){
  x[34:49]
})
colnames(MAF_list)=c('AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','ExAC_MAF','ExAC_Adj_MAF',
                     'ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF')
Exac1000G=unique(cbind(Col8[,1:5],MAF_list))
"1. Eemaldada tähed (apply)
2. Standartne formaat, kui ei siis-1 
prcomp() standard 
dimensionality reduction methods 
integrate data to arlequin 
frequency of allel info"

colnames(AC_list) = c('AC', 'AF', 'AN','DB','SET')

library(stringr)
Exac1000G["MAF_String"] <- NA
length(Exac1000G$AFR_MAF)
Exac1000G["MAF_String"] <- str_split_fixed(Exac1000G$AFR_MAF,":",2)[1:length(Exac1000G$V1)]
Exac1000G$MAF_String

Exac1000G["MAF_String"] == Exac1000G["End"]
length(Exac1000G[2])
Exac1000G$End
Exac1000G$MAF_String
Exac1000G <- tbl_df(Exac1000G)
"Sorteerime need välja mis on tühjad"
i = 0
for(i in 1:(length(Exac1000G$AFR_MAF))){
  if(Exac1000G$AFR_MAF[i] == ""){
    Exac1000G <- Exac1000G[-c(i),]
  }
}

colnames = c('AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','ExAC_MAF','ExAC_Adj_MAF',
             'ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF')

Exac1000G %>% filter(AFR_MAF != "")
"Kas kõik klapib?"
arvud = list()
i = 0
for(i in 1:(length(Exac1000G$MAF_String))){
  if(Exac1000G$V5[i] != Exac1000G$MAF_String[i]){
    arvud <- c(arvud,i)
  }
}
i = 0
j = 0
for(i in 1:(length(Exac1000G$MAF_String))){
  for(j in arvud){
      if(j == i){
        Exac1000G <- Exac1000G[-c(i),]
      }
  }
}
Exac1000G$End == Exac1000G$MAF_String
1 == 1
i = 0
for(i in 1:(length(Exac1000G$MAF_String))){
  if(grepl("&",Exac1000G$AFR_MAF[i])){
    arvud <- c(arvud,i)
  }
}
c('AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','ExAC_MAF','ExAC_Adj_MAF',
  'ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF')
i = 0
j = 0
m = 0
for(i in 1:(length(Exac1000G$MAF_String))){
  for(j in arvud){
    if(i == j){
        th <- str_split_fixed(Exac1000G$ExAC_MAF[j],"&",2)[2]
        Exac1000G$ExAC_MAF[j] <- th
    }
  }
}

"String numbers to float"
ExAC_EAS
i = 0
c = 0
d = 0
for(i in 1:(length(Exac1000G$EUR_MAF))){
  th = 0
  a = ""
  a <- str_split_fixed(Exac1000G$EUR_MAF[i],"&",2)[1]
  th <- as.numeric(str_split_fixed(a,":",2)[2])
  if(th > 0.9){
    c <- (c+1)
  }
  Exac1000G$EUR_MAF[i] <- th
}

library(xlsx) #load the package
write.csv(Exac1000G,file = "Exac1000G_Tidy.csv", row.names = FALSE)

" AEG PROOVIDA PCA-d SISSE TUUA OH BOY"

standardize <- function(x){(x-mean(x))}
scaled_classes <- apply(testing,2,function(x)(x-mean(x)))

Exac1000G$AFR_MAF <- as.numeric(Exac1000G$AFR_MAF)
Exac1000G$AMR_MAF <- as.numeric(Exac1000G$AMR_MAF)
Exac1000G$EAS_MAF <- as.numeric(Exac1000G$EAS_MAF)
Exac1000G$EUR_MAF <- as.numeric(Exac1000G$EUR_MAF)
Exac1000G$SAS_MAF <- as.numeric(Exac1000G$SAS_MAF)
Exac1000G$AA_MAF <- as.numeric(Exac1000G$AA_MAF)
Exac1000G$EA_MAF <- as.numeric(Exac1000G$EA_MAF)
Exac1000G$ExAC_MAF <- as.numeric(Exac1000G$ExAC_MAF)
Exac1000G$ExAC_Adj_MAF <- as.numeric(Exac1000G$ExAC_Adj_MAF)
Exac1000G$ExAC_AFR_MAF <- as.numeric(Exac1000G$ExAC_AFR_MAF)
Exac1000G$ExAC_AMR_MAF <- as.numeric(Exac1000G$ExAC_AMR_MAF)
Exac1000G$ExAC_EAS_MAF <- as.numeric(Exac1000G$ExAC_EAS_MAF)
Exac1000G$ExAC_FIN_MAF <- as.numeric(Exac1000G$ExAC_FIN_MAF)
Exac1000G$ExAC_NFE_MAF <- as.numeric(Exac1000G$ExAC_NFE_MAF)
Exac1000G$ExAC_OTH_MAF <- as.numeric(Exac1000G$ExAC_OTH_MAF)
Exac1000G$ExAC_SAS_MAF <- as.numeric(Exac1000G$ExAC_SAS_MAF)

testing2 <- total_data[2:4]
testing2$AFR_MAF[4]
MAF_means <- rowMeans(testing2[1:5])
testing2 <- t(testing2)
pca = prcomp(testing2)
apply(testing2,2,var)
plot(pca)
plot(pca$rotation[,3],pca$rotation[,2])

"0-d NA-deks"
"UK NL ja EST AFid kokku AF_Est, AF_NL, AF_UK
merge-da kõik kokku 
oodata lisaandmeid 
ExAC andmed eraldi tabel"





UK10KPharmacoVariants$V8 <- as.character(UK10KPharmacoVariants$V8)
NLGoPharmacoVariants$V8 <- as.character(NLGoPharmacoVariants$V8)
EstonianPharmacoVariantsCol8$V8 <- as.character(EstonianPharmacoVariantsCol8$V8)
new=strsplit(NLGoPharmacoVariants$V8,';')
NL_list=ldply(new,function(x){
  x[34:38]
})
colnames(MAF_list)=c('AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF','AA_MAF','EA_MAF','ExAC_MAF','ExAC_Adj_MAF',
                     'ExAC_AFR_MAF','ExAC_AMR_MAF','ExAC_EAS_MAF','ExAC_FIN_MAF','ExAC_NFE_MAF','ExAC_OTH_MAF','ExAC_SAS_MAF')
Exac1000G_Original=unique(cbind(Col8[,1:5],MAF_list))

colnames(UK_list) = c('AC', 'AF', 'AN','DB','SET')

NL_list <- do.call(rbind.data.frame, new)

new=strsplit(UK10KPharmacoVariants$V8,";")

UK_list <- do.call(rbind.data.frame, new)
UK_list <- cbind(UK_list,UK10KPharmacoVariants$V3)
UK_list <- subset(UK_list, select= c('c..AF.0.829807....AF.0.537688....AF.0.52949....AF.0.597593...','UK10KPharmacoVariants$V3'))
colnames(UK_list) = c('AF', 'id')

UK10KPharmacoVariants$V8 <- as.character(UK10KPharmacoVariants$V8)
EstonianPharmacoVariantsCol8$V8 <- as.character(EstonianPharmacoVariantsCol8$V8)
NLGoPharmacoVariants$V8 <- as.character(NLGoPharmacoVariants$V8)


new=strsplit(NLGoPharmacoVariants$V8,';')
NL_list <- do.call(rbind.data.frame, new)
colnames(NL_list) = c('AC', 'AF', 'AN','DB','SET')
NL_list <- cbind(NL_list,NLGoPharmacoVariants$V3)

NL_list <- subset(NL_list, select= c('AF','NLGoPharmacoVariants$V3'))

new=strsplit(EstonianPharmacoVariantsCol8$V8,';')
EST_list <- do.call(rbind.data.frame, new)
colnames(EST_list) = c('AC', 'DB', 'AN','AF','SET')
EST_list <- cbind(EST_list,EstonianPharmacoVariantsCol8$V3)

EST_list <- subset(EST_list, select= c('AF','EstonianPharmacoVariantsCol8$V3'))
colnames(NL_list) = c('NL_AF','id')
colnames(UK_list) = c('UK_AF','id')
colnames(EST_list) = c('EST_AF','id')

NL_list_AF<- NL_list$NL_AF

UK_list_AF <- UK_list$UK_AF



subset(NLGoPharmacoVariants, NLGoPharmacoVariants$V3 == '.')

UK10KPharmacoVariants<-UK10KPharmacoVariants[!(UK10KPharmacoVariants$V3 == '.'),]

EstonianPharmacoVariantsCol8$V8 <- as.character(EstonianPharmacoVariantsCol8$V8)

new=strsplit(EstonianPharmacoVariantsCol8$V8,";")
est_list <- do.call(rbind.data.frame, new)
colnames(est_list) = c('AC', 'AN','DB','AF','SET')
est_list <- do.call(rbind.data.frame, new)
EST_list_AF <- EST_list$EST_AF

colnames(total_data_interesting) <- colnames(total_data)

global_gene_data <- NULL
global_gene_data <- merge(NL_list,EST_list, by='id',type='inner')
global_gene_data <- merge(global_gene_data,UK_list, by='id',type='inner')
length(EST_list_AF) <- length(UK_list_AF)
length(NL_list_AF) <- length(EST_list_AF)
global_gene_data <- global_gene_data[!duplicated(global_gene_data), ]

colnames(UK10KPharmacoVariants) = c('V1', 'V2', 'id','AF','SET')

global_gene_data <- cbind(UK10KPharmacoVariants,NL_list_AF,UK_list_AF,EST_list_AF)

global_gene_data <- global_gene_data[complete.cases(global_gene_data), ]
global_gene_data$V4 <- NULL
global_gene_data$V5 <- NULL
global_gene_data$V6 <- NULL
global_gene_data$V7 <- NULL
global_gene_data$V8 <- NULL


colnames(global_gene_data) = c('refsnp_id','NL_AF','EST_AF','UK_AF')
library(stringr)
i = 0
global_gene_data$af4 <- 0
for(i in 1:(length(global_gene_data$EST_AF))){
  th = 0
  a = ""
  a <- str_split_fixed(global_gene_data$EST_AF[i],"=",2)[2]
  th <- as.numeric(a)
  global_gene_data$af4[i] <- th
}

i = 0
global_gene_data$af3 <- 0
for(i in 1:(length(global_gene_data$UK_AF))){
  th = 0
  a = ""
  a <- str_split_fixed(global_gene_data$UK_AF[i],"=",2)[2]
  th <- as.numeric(a)
  global_gene_data$af3[i] <- th
}
i = 0
global_gene_data$af <- 0
for(i in 1:(length(global_gene_data$NL_AF))){
  th = 0
  a = ""
  a <- str_split_fixed(global_gene_data$NL_AF[i],"=",2)[2]
  th <- as.numeric(a)
  global_gene_data$af[i] <- th
}

global_gene_data$EST_AF <- NULL
global_gene_data$UK_AF <- NULL
global_gene_data$NL_AF <- NULL

colnames(global_gene_data) = c('refsnp_id','EST_AF','UK_AF','NL_AF')

gene_data <- as.data.frame(gene_data)
global_gene_data <- merge(global_gene_data,gene_data,by='refsnp_id',type='inner')
global_gene_data  <- global_gene_data [!duplicated(global_gene_data $refsnp_id,global_gene_data$external_gene_name), ]
global_gene_data <- global_gene_data[c('refsnp_id','NL_AF','UK_AF','EST_AF','external_gene_name','VIP')]

length(global_gene_data$NL_AF) <- 4569
total_data <- merge(total_data,Exac1000G,by='refsnp_id')
install.packages("gtools")
library(gtools)

colnames(Exac1000G) = c('V1','V2','refsnp_id')
total_data$EST_AF.y <- NULL
total_data$NL_AF.y <- NULL
total_data$UK_AF.y <- NULL
total_data <- merge(global_gene_data,total_data,by='refsnp_id',type='inner')
colnames(total_data)[2] <- "EST_AF"
colnames(total_data)[3] <- "UK_AF"
colnames(total_data)[4] <- "NL_AF"

total_data_smart <- smartbind(global_gene_data,gene_data)
index <- which(duplicated(total_data))
total_data <- total_data[!duplicated(total_data), ]
global_gene_data <- as.data.frame(global_gene_data)
global_gene_data <- do.call(rbind.data.frame, global_gene_data)


total_data_clean <- total_data[!(duplicated(total_data[c("refsnp_id","NL_AF")]) | duplicated(total_data[c("refsnp_id","NL_AF")], fromLast = TRUE)), ]

write.csv(total_data, file = "total_data.csv")

"loome tabeli kus on variatsioonide sagedused"
pop_freq['UK_AF'] <- 0
UK_AF <- c(0)
NL_AF <- c(0)
EST_AF <- c(0)
pop_freq <- data.frame(UK_AF,NL_AF,EST_AF)
nimekiri <- c()
biggest_pop <- c()
pop_freq['UK_AF']<-0
pop_freq['NL_AF']<-0
pop_freq['EST_AF']<-0
for(i in 1:(length(global_gene_data$UK_AF))){
  number = max(global_gene_data$NL_AF[i],global_gene_data$UK_AF[i],global_gene_data$EST_AF[i])
  if(global_gene_data$UK_AF[i]==number ){
    nimekiri <- c(nimekiri,'UK_AF')
    arv <- pop_freq['UK_AF']
    pop_freq['UK_AF'] = (arv+1)
    biggest_pop <-c(biggest_pop,3)
  }
  else{
    if(global_gene_data$NL_AF[i]==number ){
      nimekiri <- c(nimekiri,'NL_AF')
      arv <- pop_freq['NL_AF']
      pop_freq['NL_AF'] = (arv+1)
      biggest_pop <-c(biggest_pop,2)
    }
    else{
      if(global_gene_data$EST_AF[i]==number ){
        nimekiri <- c(nimekiri,'EST_AF')
        arv <- pop_freq['UK_AF']
        pop_freq['EST_AF'] = (arv+1)
        biggest_pop <-c(biggest_pop,1)
    }
  }
  }
}

global_gene_data$biggest <- biggest_pop

total_data$V1 <- NULL
total_data$V2 <-NULL
total_data$NA.1<- NULL


barplot(prop.table(table(nimekiri)))


"VIP variatsioonide proov"
vip_pop <- c()

for(i in 1:(length(total_data$UK_AF))){
  number = max(total_data$NL_AF[i],total_data$UK_AF[i],total_data$EST_AF[i])
  if(total_data$UK_AF[i]==number & total_data$VIP[i]==1){
    vip_pop <- c(vip_pop,'UK_AF')
  }
  else{
    if(total_data$NL_AF[i]==number & total_data$VIP[i]==1){
      vip_pop <- c(vip_pop,'NL_AF')
    }
    else{
      if(total_data$EST_AF[i]==number & total_data$VIP[i]==1){
        vip_pop <- c(vip_pop,'EST_AF')
      }
    }
  }
}

barplot(prop.table(table(vip_pop)))



colnames(total_data)= c(colnames(total_data)[1:12],colnames(Exac1000G_Original)[5:21])
maf_test = total_data[14:18]
pca_maf = prcomp(t(maf_test))


plot(pca$rotation[,1],pca$rotation[,2],main="PCA rotation",xlab = "Netherlands",ylab="United Kingdom")
plot(pca$rotation[,1],pca$rotation[,3],main="PCA rotation",xlab = "Netherlands",ylab="Estonia")
plot(pca$rotation[,2],pca$rotation[,1],main="PCA rotation",xlab = "United Kingdom",ylab="Netherlands")
plot(pca$rotation[,2],pca$rotation[,3],main="PCA rotation",xlab = "United Kingdom",ylab="Estonia")
plot(pca$rotation[,3],pca$rotation[,1],main="PCA rotation",xlab = "Estonia",ylab="Netherlands")
plot(pca$rotation[,3],pca$rotation[,2],main="PCA rotation",xlab = "Estonia",ylab="United Kingdom")

plot(pca_maf$rotation[,1],pca_maf$rotation[,2])


AFR_MAF <- c(0)
AMR_MAF <- c(0)
EAS_MAF <- c(0)
EUR_MAF <- c(0)
SAS_MAF <- c(0)
ex_pop_freq <- data.frame(AFR_MAF,AMR_MAF,EAS_MAF,EUR_MAF,SAS_MAF)
nimekiri <- c()
ex_pop_freq['AFR_MAF']<-0
ex_pop_freq['AMR_MAF']<-0
ex_pop_freq['EAS_MAF']<-0
ex_pop_freq['EUR_MAF']<-0
ex_pop_freq['SAS_MAF']<-0


for(i in 1:(length(maf_test$EAS_MAF))){
  number = max(maf_test[i,])
  if(maf_test$AFR_MAF[i]==number ){
    nimekiri <- c(nimekiri,'AFR_MAF')
    arv <- ex_pop_freq['AFR_MAF']
    ex_pop_freq['AFR_MAF'] = (arv+1)
  }
  else{
    if(maf_test$AMR_MAF[i]==number ){
      nimekiri <- c(nimekiri,'AMR_MAF')
      arv <- ex_pop_freq['AMR_MAF']
      ex_pop_freq['AMR_MAF'] = (arv+1)
    }
    else{
      if(maf_test$EAS_MAF[i]==number ){
        nimekiri <- c(nimekiri,'EAS_MAF')
        arv <- ex_pop_freq['EAS_MAF']
        ex_pop_freq['EAS_MAF'] = (arv+1)
      }
      else{
        if(maf_test$EUR_MAF[i]==number ){
          nimekiri <- c(nimekiri,'EUR_MAF')
          arv <- ex_pop_freq['EUR_MAF']
          ex_pop_freq['EUR_MAF'] = (arv+1)
        }
        else{
          if(maf_test$SAS_MAF[i]==number ){
            nimekiri <- c(nimekiri,'SAS_MAF')
            arv <- ex_pop_freq['SAS_MAF']
            ex_pop_freq['SAS_MAF'] = (arv+1)
          }
          
        }
    }
  }
  }}


barplot(prop.table(table(nimekiri)))


"Uurime kas leiame variatsioone mis eestlastel unikaalne"
i = 0
d <- c()
for(i in 1:(length(total_data$EST_AF))){
  if((total_data$EST_AF[i]-0.05)>total_data$NL_AF[i]){
    if((total_data$EST_AF[i]-0.05)>total_data$UK_AF[i]){
      d <- c(d,i)
    }
  }
}
est_unique <- rbind(total_data[553,],total_data[3672,],total_data[554,],total_data[3673,],total_data[4231,])

"huvitav, vaatame sama ka inglastel ja hollandlastel"
i = 0
d <- c()
for(i in 1:(length(total_data$NL_AF))){
  if((total_data$NL_AF[i]-0.13)>total_data$EST_AF[i]){
    if((total_data$NL_AF[i]-0.13)>total_data$UK_AF[i]){
      d <- c(d,i)
    }
  }
}

nl_unique <- rbind(total_data[1329,],total_data[2384,],total_data[2384,],total_data[3426,])


i = 0
d <- c()
for(i in 1:(length(total_data$UK_AF))){
  if((total_data$UK_AF[i]-0.05)>total_data$EST_AF[i]){
    if((total_data$UK_AF[i]-0.05)>total_data$NL_AF[i]){
      d <- c(d,i)
    }
  }
}

uk_unique <- rbind(total_data[85,],total_data[332,],total_data[2763,],total_data[4389,])

global_gene_data.idless <- cbind(global_gene_data[,2],global_gene_data[,3],global_gene_data[,4])
colnames(global_gene_data.idless) <- colnames(global_gene_data)[2:4]
"Proovime clusterdada 3 rahvaste AF andmeid"

library(cluster)
global_gene_data.stand <- scale(global_gene_data[-1])
k.means.fit <- kmeans(global_gene_data[], 3) # k = 3
attributes(k.means.fit)
k.means.fit$centers
k.means.fit$cluster
k.means.fit$size

clusplot(global_gene_data, k.means.fit$cluster, main='2D representation of the Cluster solution',
         color=TRUE, shade=TRUE,
         labels=2, lines=0)

wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")}

wssplot(global_gene_data.stand, nc=6) 

"Samamoodi clusterdust teeme MAF andmetega"

maf_test.stand <- scale(maf_test[-1])
k.means.fit <- kmeans(maf_test, 5) # k = 5
attributes(k.means.fit)
k.means.fit$centers
k.means.fit$cluster
k.means.fit$size

clusplot(maf_test, k.means.fit$cluster, main='2D representation of the Cluster solution',
         color=TRUE, shade=TRUE,
         labels=2, lines=0)

"testime 0 kaotamist"

total_data_no0 <- total_data[total_data == 0] <- NA

maf_test.stand <- scale(maf_test[-1])
k.means.fit <- kmeans(maf_test, 5) # k = 5
attributes(k.means.fit)
k.means.fit$centers
k.means.fit$cluster
k.means.fit$size

total_data[is.na(total_data[])] <- 0 
clusplot(maf_test, k.means.fit$cluster, main='2D representation of the Cluster solution',
         color=TRUE, shade=TRUE,
         labels=2, lines=0)

ExAC_0 <- subset(total_data,select= AFR_MAF:ExAC_SAS_MAF)
ExAC_NA <- ExAC_0
ExAC_NA[ExAC_NA==0]<- NA


pca_exac_0 <- prcomp(na.omit(t(ExAC_0)))



plot(pca_exac_0$rotation[,2],pca_exac_0$rotation[,3],main="PCA rotation of ExAC data",xlab = "1",ylab="2")
devtools::install_github("https://github.com/vqv/ggbiplot")
summary(pca_exac_0)
summary(pca)
library(ggbiplot)

"Proovime Tõnise soovitatud dbscan() meetodit"
install.packages("fpc")
install.packages("dbscan")
library("fpc")
# Compute DBSCAN using fpc package
set.seed(123)
db <- dbscan(true_maf, eps = 0.15, MinPts = 5)
# Plot DBSCAN results
plot(db, true_maf, main = "DBSCAN", frame = FALSE)
#...Nope
install.packages("devtools")
library(devtools)
devtools::install_github("kassambara/factoextra")
install.packages('Rcpp')
install.packages('graphics')
install.packages('stats')
install.packages('methods')
#, graphics, stats, methods
library("factoextra")
fviz_cluster(db, true_maf, geom = "point")
data("iris")pl
x <- as.matrix(iris[, 1:4])
db <- dbscan(x, eps = .4, MinPts = 4)
pairs(x, col = db$cluster + 1L)
opt <- optics(x, eps = 1, minPts = 4, eps_cl = .4)

pairs(true_maf, col = db$cluster + 1L)
lof <- lof(true_maf, k = 4)
db <- dbscan(global_gene_data.idless, eps = 5, MinPts = 5)
# Plot DBSCAN results
plot(db, global_gene_data.idless, main = "DBSCAN", frame = FALSE)

################ VAHEPAUS, RAKENDAN TÕNISE GOOGLE DOC MÄRKMEID.
#### ESITEKS: HUVITAVAMATE VARIANTIDE TABEL (KUS ON EXAC ANDMED) TABEL JA SEALT UURIN EKSTREEMSEID VARIANTE
total_data_interesting<- subset(total_data, ExAC_Adj_MAF > 0)
# Leian Eesti ekstreemsemad variandid ehk need mida esinevad meil alla 0.3% ja üle 99%
est_extreme <-c(0,0)
uk_extreme <-c(0,0)
nl_extreme <-c(0,0)
extreme_est <- data.frame()
for(i in 1:(length(total_data_interesting$UK_AF))){
  if(total_data_interesting$EST_AF[i] < 0.003){
    extreme_est<- rbind(extreme_est, total_data_interesting[i,])
    est_extreme[1] <- (est_extreme[1]+1)
  }
  if(total_data_interesting$EST_AF[i] > 0.99){
    extreme_est<- rbind(extreme_est, total_data_interesting[i,])
    est_extreme[2] <- (est_extreme[2]+1)
  }
  if(total_data_interesting$UK_AF[i] < 0.003){
    uk_extreme[1] <- (uk_extreme[1]+1)
  }
  if(total_data_interesting$UK_AF[i] > 0.99){
    uk_extreme[2] <- (uk_extreme[2]+1)
  }
  if(total_data_interesting$NL_AF[i] < 0.003){
    nl_extreme[1] <- (nl_extreme[1]+1)
  }
  if(total_data_interesting$NL_AF[i] > 0.99){
    nl_extreme[2] <- (nl_extreme[2]+1)
  }
  
}

extreme_est = extreme_est[c('refsnp_id','NL_AF','UK_AF','EST_AF','external_gene_name','VIP')]

write.csv(extreme_est,file = "Extreme_EST.csv", row.names = FALSE)

## Tõnise soovitatud PCA muudatuste osa

data=replicate(15,seq(from=1,to=10,length.out=100)+rnorm(100))
colnames(data)=paste('Pop',1:15)
rownames(data)=paste('Variant',1:100)
#Suvalised VIP t2hised
VIP=sample(c(2,1),size=100,replace=T,prob=c(0.25,0.75))

pca=prcomp(data,scale=T)
summary(pca) # two PCs for cumulative proportion of >80% 
newdat<-data.frame(pca$x[,1:2])
newdat$VIP=VIP
plot(newdat$PC1,newdat$PC2,col=newdat$VIP)

pca=prcomp(t(data),scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])
plot(newdat$PC1,newdat$PC2,type='n')
text(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")

##Teeme vahepeal df kus on KÕIK variandid ja populatsioonid
all_pop_variants <- subset(total_data_interesting,select= c(NL_AF:EST_AF,AFR_MAF:ExAC_SAS_MAF))

## Ainult VIP 
colnames(global_gene_data.idless) <- c(colnames(global_gene_data[2]),colnames(global_gene_data[3]),colnames(global_gene_data[4]))
## Teeme plot kus on 3 rahva variandid ja eristame neid värvidega
pca=prcomp(global_gene_data.idless,scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:3])
newdat$id <- global_gene_data$refsnp_id
newdat$biggest = global_gene_data$biggest
newdat$VIP = global_gene_data$VIP
## Ainult VIP 
#ewdat <- subset(newdat,VIP==1)
tiff("Plot1.tiff",width=14,height = 8,units = 'in',res=1200)
plot(newdat$PC2,newdat$PC3,type='n')
text(newdat$PC2,newdat$PC3, row.names(newdat), cex=0.6, pos=4, col="red")
plot(newdat$PC2,newdat$PC3,xaxt = "n",yaxt="n",col=ifelse(newdat$biggest==3,'red'
,ifelse(newdat$biggest==2,'green','blue')),pch=ifelse(newdat$VIP==0,15,17),cex= ifelse(newdat$VIP==0,1,2),
xlab = "PC2",
ylab="PC3")
legend("bottomright",col=c("blue", "red","green"), pch=15,
       legend=c("EST_AF","UK_AF","NL"), bty="n", xjust=1,cex=1.4 ,seg.len=0.5)
legend("right",col=c("black"), pch=17,
       legend=c("VIP"),cex=1.4, bty="n", xjust=1, seg.len=0.5)
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)
dev.off()
global_gene_data.new <- merge(global_gene_data,total_data, by='refsnp_id' ,type='inner')
global_gene_data.new <- global_gene_data.new[!duplicated(global_gene_data.new$refsnp_id), ]


# PC1 = 0.58 * NL_AF + 0.58 * UK_AF + 0.58 * EST_AF
# PC2 = 0.81 * NL_AF + 0.28 * UK_AF + 0.52 * EST_AF
# PC3 = 0.14 * NL_AF + 0.76 * UK_AF + 0.63 * EST_AF


## Teema sama ka 1000g andmetega
nimekiri <- c()
biggest_pop <- c()

for(i in 1:(length(maf_test$EAS_MAF))){
  number = max(maf_test[i,])
  if(maf_test$AFR_MAF[i]==number ){
    biggest_pop <- c(biggest_pop,1)
  }
  else{
    if(maf_test$AMR_MAF[i]==number ){
      biggest_pop <- c(biggest_pop,2)
    }
    else{
      if(maf_test$EAS_MAF[i]==number ){
        biggest_pop <- c(biggest_pop,3)
      }
      else{
        if(maf_test$EUR_MAF[i]==number ){
          biggest_pop <- c(biggest_pop,4)
        }
        else{
          if(maf_test$SAS_MAF[i]==number ){
            biggest_pop <- c(biggest_pop,5)
          }
          
        }
      }
    }
  }}




pca=prcomp(maf_test,scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])

plot(newdat$PC1,newdat$PC2,type='n')
text(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")
newdat$biggest <- biggest_pop
plot(newdat$PC1,newdat$PC2,xaxt = "n",yaxt="n",col=ifelse(newdat$biggest==1,'yellow'
,ifelse(newdat$biggest==2,'red',ifelse(newdat$biggest==3,'green',ifelse(newdat$biggest==4,'blue','purple'))))
,pch=15,main="1000G population PCA",xlab = "PC1",ylab="PC2")
     

legend("topright",inset=c(-0.26,0),col=c("blue", "red","green","purple","yellow"), pch=15,cex=1.4,
       legend=c("EUR_MAF","AMR_MAF","EAS_MAF","SAS_MAF","AFR_MAF"), bty="n", border = "black",xjust=1, seg.len=0.5)

axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)

pca_eur <- subset(newdat,biggest==4)
pca_afr <- subset(newdat,biggest==1)
pca_amr <- subset(newdat,biggest==2)
pca_eas <- subset(newdat,biggest==3)
pca_sas <- subset(newdat,biggest==5)
t.test(pca_amr$PC2,pca_afr$PC2,var.equal=TRUE, paired=FALSE)

## Teema sama ka 1000g andmetega
ExAC_data <- subset(total_data,select= c(VIP,ExAC_AFR_MAF:ExAC_SAS_MAF))
ExAC_data <- subset(ExAC_data,ExAC_AMR_MAF>0)
ExAC_data_VIP <- ExAC_data$VIP
ExAC_data <- subset(ExAC_data,select= c(ExAC_AFR_MAF:ExAC_SAS_MAF))
nimekiri <- c()
biggest_pop <- c()

for(i in 1:(length(ExAC_data$ExAC_AFR_MAF))){
  number = max(ExAC_data[i,])
  if(ExAC_data$ExAC_AFR_MAF[i]==number ){
    biggest_pop <- c(biggest_pop,1)
  }
  else{
    if(ExAC_data$ExAC_AMR_MAF[i]==number ){
      biggest_pop <- c(biggest_pop,2)
    }
    else{
      if(ExAC_data$ExAC_EAS_MAF[i]==number ){
        biggest_pop <- c(biggest_pop,3)
      }
      else{
        if(ExAC_data$ExAC_FIN_MAF[i]==number ){
          biggest_pop <- c(biggest_pop,4)
        }
        else{
          if(ExAC_data$ExAC_NFE_MAF[i]==number ){
            biggest_pop <- c(biggest_pop,5)
          }
          else{
            if(ExAC_data$ExAC_OTH_MAF[i]==number ){
              biggest_pop <- c(biggest_pop,6)
            }
            else{
              if(ExAC_data$ExAC_SAS_MAF[i]==number ){
                biggest_pop <- c(biggest_pop,7)
              }
              
          
        }}}}
      }
    }
  }




pca=prcomp(subset(ExAC_data,select= c(ExAC_AFR_MAF:ExAC_SAS_MAF)),scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])

plot(newdat$PC1,newdat$PC2,type='n')
text(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")
newdat$biggest <- biggest_pop
newdat$VIP <- ExAC_data_VIP
plot(newdat$PC1,newdat$PC2,col=ifelse(newdat$biggest==1,'yellow'
                                      ,ifelse(newdat$biggest==2,'red',ifelse(newdat$biggest==3,'green',
                        ifelse(newdat$biggest==4,'blue',ifelse(newdat$biggest==5,'purple',ifelse(newdat$biggest==6,'hotpink','orange'))))))
     ,pch=ifelse(newdat$VIP==0,15,17),cex= ifelse(newdat$VIP==0,1,2),main="ExAC population PCA",xlab = "PC1",ylab="PC2")

legend("topleft",col=c("yellow","red", "green","blue","purple","hotpink","orange"), pch=15,
       legend=c("ExAC_AFR","ExAC_AMR","ExAC_EAS","ExAC_FIN","ExAC_NFE","ExAC_OTH","ExAC_SAS"), bty="n", xjust=1, seg.len=0.5)

#### Loome 2 tabelit, üks kus on kõik variandid mida on eestlastel vähem ja teine kus on kõik variandid mida
# on eestlastel rohkem

max_est <- data.frame()
total_data_nocode <- subset(total_data,select= c(refsnp_id:SAS_MAF))
for(i in 1:(length(total_data_3and1000$UK_AF))){
  number = max(total_data_3and1000$NL_AF[i],total_data_3and1000$UK_AF[i],total_data_3and1000$EST_AF[i],total_data_3and1000$AFR_MAF[i],
               total_data_3and1000$AMR_MAF[i],total_data_3and1000$EAS_MAF[i],total_data_3and1000$EUR_MAF[i],total_data_3and1000$SAS_MAF[i]
               #,total_data_3and1000$AA_MAF[i],total_data_3and1000$EA_MAF[i],total_data_3and1000$ExAC_AMR_MAF[i],total_data_3and1000$ExAC_AFR_MAF[i],total_data_3and1000$ExAC_EAS_MAF[i],
               #total_data_3and1000$ExAC_FIN_MAF[i],total_data_3and1000$ExAC_NFE_MAF[i],total_data_3and1000$ExAC_OTH_MAF[i],total_data_3and1000$ExAC_SAS_MAF[i],total_data_3and1000$ExAC_SAS_MAF[i]
  )
  if(total_data_3and1000$EST_AF[i]==number ){
    #if(total_data_3and1000$EST_AF[i]-0.015 > total_data_3and1000$NL_AF[i] & total_data_3and1000$EST_AF[i]-0.015> total_data_3and1000$UK_AF[i] & total_data_3and1000$EST_AF[i] >0.9){
      max_est <- rbind(max_est,total_data_3and1000[i,])

    #}
  }
}
max_est <- max_est[!duplicated(max_est), ]

# Miinimumid
min_est <- data.frame()
for(i in 1:(length(total_data_3and1000$UK_AF))){
  number = min(total_data_3and1000$NL_AF[i],total_data_3and1000$UK_AF[i],total_data_3and1000$EST_AF[i],total_data_3and1000$AFR_MAF[i],
               total_data_3and1000$AMR_MAF[i],total_data_3and1000$EAS_MAF[i],total_data_3and1000$EUR_MAF[i],total_data_3and1000$SAS_MAF[i]
               #,total_data_3and1000$AA_MAF[i],total_data_3and1000$EA_MAF[i],total_data_3and1000$ExAC_AMR_MAF[i],total_data_3and1000$ExAC_AFR_MAF[i],total_data_3and1000$ExAC_EAS_MAF[i],
               #total_data_3and1000$ExAC_FIN_MAF[i],total_data_3and1000$ExAC_NFE_MAF[i],total_data_3and1000$ExAC_OTH_MAF[i],total_data_3and1000$ExAC_SAS_MAF[i],total_data_3and1000$ExAC_SAS_MAF[i]
  )
  if(total_data_3and1000$EST_AF[i]==number ){
    #if(total_data_3and1000$EST_AF[i]-0.015 > total_data_3and1000$NL_AF[i] & total_data_3and1000$EST_AF[i]-0.015> total_data_3and1000$UK_AF[i] & total_data_3and1000$EST_AF[i] >0.9){
    min_est <- rbind(min_est,total_data_3and1000[i,])
    #}
  }
}
min_est <- min_est[!duplicated(min_est), ]

max_est_top <- subset(max_est,protsent > 0.12)
min_est_top <- subset(min_est,protsent > 0.13)

extreme_est <- rbind(max_est_top,min_est_top)
extreme_est = extreme_est[c('refsnp_id','NL_AF','UK_AF','EST_AF','external_gene_name','VIP')]
write.csv(extreme_est,file = "Extreme_EST.csv", row.names = FALSE)

max_est_vip <- subset(max_est,VIP > 0)
min_est_vip <- subset(min_est,VIP > 0)
extreme_est_vip <- rbind(max_est_vip,min_est_vip)
write.csv(extreme_est_vip,file = "Extreme_EST_VIP.csv", row.names = FALSE)
##### Mis on hollandlastel rohkem
max_nl <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$NL_AF[i]==number ){
    #if(total_data_new$nl_AF[i]-0.015 > total_data_new$NL_AF[i] & total_data_new$nl_AF[i]-0.015> total_data_new$UK_AF[i] & total_data_new$nl_AF[i] >0.9){
    max_nl <- rbind(max_nl,total_data_new[i,])
    #}
  }
}
max_nl <- max_nl[!duplicated(max_nl), ]

# Miinimumid
min_nl <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$NL_AF[i]==number ){
    #if(total_data_new$nl_AF[i]-0.015 > total_data_new$NL_AF[i] & total_data_new$nl_AF[i]-0.015> total_data_new$UK_AF[i] & total_data_new$nl_AF[i] >0.9){
    min_nl <- rbind(min_nl,total_data_new[i,])
    #}
  }
}
min_nl <- min_nl[!duplicated(min_nl), ]

max_est_top <
extreme_nl <- rbind(max_nl,min_nl)
extreme_nl = extreme_nl[c('refsnp_id','NL_AF','UK_AF','EST_AF','external_gene_name','VIP')]
write.csv(extreme_nl,file = "Extreme_NL.csv", row.names = FALSE)

##### Mis on inglastel rohkem
max_uk <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$UK_AF[i]==number ){
    #if(total_data_new$uk_AF[i]-0.015 > total_data_new$uk_AF[i] & total_data_new$uk_AF[i]-0.015> total_data_new$UK_AF[i] & total_data_new$uk_AF[i] >0.9){
    max_uk <- rbind(max_uk,total_data_new[i,])
    #}
  }
}
max_uk <- max_uk[!duplicated(max_uk), ]

# Miinimumid
min_uk <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$UK_AF[i]==number ){
    #if(total_data_new$uk_AF[i]-0.015 > total_data_new$uk_AF[i] & total_data_new$uk_AF[i]-0.015> total_data_new$UK_AF[i] & total_data_new$uk_AF[i] >0.9){
    min_uk <- rbind(min_uk,total_data_new[i,])
    #}
  }
}
min_uk <- min_uk[!duplicated(min_uk), ]

extreme_uk <- rbind(max_uk,min_uk)
extreme_uk = extreme_uk[c('refsnp_id','NL_AF','UK_AF','EST_AF','external_gene_name','VIP')]
write.csv(extreme_uk,file = "Extreme_UK.csv", row.names = FALSE)



pca=prcomp(t(all_pop_variants),scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])
plot(newdat$PC1,newdat$PC2,type='n')
textplot(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")


### Teeme sageduse plot kus on kõik andmed....
all_pop_variants.true <- all_pop_variants
all_pop_variants <- total_data[c("NL_AF","UK_AF","EST_AF","AFR_MAF","AMR_MAF","EAS_MAF"
                                       ,"EUR_MAF","SAS_MAF")]
## Teema sama ka 1000g andmetega
nimekiri <- c()
biggest_pop <- c()

for(i in 1:(length(all_pop_variants$EAS_MAF))){
  number = max(all_pop_variants[i,])
  if(all_pop_variants$AFR_MAF[i]==number ){
    biggest_pop <- c(biggest_pop,1)
  }
  else{
    if(all_pop_variants$AMR_MAF[i]==number ){
      biggest_pop <- c(biggest_pop,2)
    }
    else{
      if(all_pop_variants$EAS_MAF[i]==number ){
        biggest_pop <- c(biggest_pop,3)
      }
      else{
        if(all_pop_variants$EUR_MAF[i]==number ){
          biggest_pop <- c(biggest_pop,4)
        }
        else{
          if(all_pop_variants$SAS_MAF[i]==number ){
            biggest_pop <- c(biggest_pop,5)
          }
          else{
            if(all_pop_variants$EST_AF[i]==number ){
              biggest_pop <- c(biggest_pop,6)
            }
            else{
              if(all_pop_variants$NL_AF[i]==number ){
                biggest_pop <- c(biggest_pop,7)
              }
              else{
                if(all_pop_variants$UK_AF[i]==number ){
                  biggest_pop <- c(biggest_pop,8)
                }
              }
            }
          }
        }
      }  
        }
      }
    }
  

pca=prcomp(all_pop_variants,scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])

plot(newdat$PC1,newdat$PC2,type='n')
text(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")
newdat$id <- total_data_interesting$refsnp_id
newdat$biggest <- biggest_pop
newdat$VIP = total_data$VIP
plot(newdat$PC1,newdat$PC2,col=ifelse(newdat$biggest==1,'yellow'
                                      ,ifelse(newdat$biggest==2,'red',ifelse(newdat$biggest==3
                ,'green',ifelse(newdat$biggest==4,'blue',ifelse(newdat$biggest==5,'purple',
                ifelse(newdat$biggest==6,'black',ifelse(newdat$biggest==7,'orange','brown')))))))
     ,pch=ifelse(newdat$VIP==0,15,1)
,main="3 country and 1000G population PCA",xlab = "PC1",ylab="PC2")
par(xpd=TRUE)
legend("topright",inset=c(-0.1,0),col=c("blue", "red","green","purple","yellow","black","orange","brown"), pch=15,
       legend=c("EUR_MAF","AMR_MAF","EAS_MAF","SAS_MAF","AFR_MAF","EST_AF","NL_AF","UK_AF"), bty="n", xjust=1, seg.len=0.5)




### 31.12.16 muudadatused.
install.packages(c("wordcloud","tm"),repos="http://cran.r-project.org")
library(wordcloud)
library(tm)

# total_data-st eemaldan ExAC_Adj_MAF
total_data_new <- total_data
total_data_new$ExAC_MAF <- NULL
total_data_new$AA_MAF <- NULL
total_data_new$EA_MAF <- NULL
total_data_new$ExAC_Adj_MAF <- NULL

total_data_new <- subset(total_data_new, ExAC_AMR_MAF > 0)
total_data_new_numbers <- subset(total_data_new,select= c(EST_AF:NL_AF,AFR_MAF:ExAC_SAS_MAF))
par(op) ## reset
pca=prcomp(t(total_data_new_numbers),scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])
plot(newdat$PC1,newdat$PC2,type='n',width=1000,height=350,res=72)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 0, 0, 4)) 
par(xpd=TRUE)
textplot(newdat$PC1,newdat$PC2, row.names(newdat),valign="top" ,cex=1, pos=4, col="blue", byrow=FALSE)


pca_eur <- subset(newdat,biggest==4)
pca_afr <- subset(newdat,biggest==1)
pca_amr <- subset(newdat,biggest==2)
pca_eas <- subset(newdat,biggest==3)
pca_sas <- subset(newdat,biggest==5)
t.test(pca_amr$PC2,pca_sas$PC2,var.equal=TRUE, paired=FALSE)
# ainult 3 rahvast ja 1000g rahvas
other_pop <- total_data
other_pop <- subset(other_pop,select= c(EST_AF:NL_AF,AFR_MAF:SAS_MAF))

pca=prcomp(t(other_pop),scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])
plot(newdat$PC1,newdat$PC2,type='n')
textplot(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")


## Proovime teha kodeeritavatest PCA kompi


# Kirjutame sama biggest teistmoodi KÕIGILE
total_data_new <- as.data.frame(total_data_new)
nimekiri <- c()
biggest_pop <- c()

for(i in 1:(length(total_data_new_numbers$EAS_MAF))){
  number = max(total_data_new_numbers[i,])
  if(total_data_new_numbers$ExAC_NFE_MAF[i]==number ){
    biggest_pop <- c(biggest_pop,1)
  }
  else{
    if(total_data_new_numbers$ExAC_AMR_MAF[i]==number ){
      biggest_pop <- c(biggest_pop,2)
    }
    else{
      if(total_data_new_numbers$ExAC_EAS_MAF[i]==number ){
        biggest_pop <- c(biggest_pop,3)
      }
      else{
        if(total_data_new_numbers$ExAC_AFR_MAF[i]==number ){
          biggest_pop <- c(biggest_pop,4)
        }
        else{
          if(total_data_new_numbers$ExAC_FIN_MAF[i]==number ){
            biggest_pop <- c(biggest_pop,5)
          }
          else{
            if(total_data_new_numbers$ExAC_OTH_MAF[i]==number ){
              biggest_pop <- c(biggest_pop,6)
            }
            else{
              if(total_data_new_numbers$EST_AF[i]==number ){
                biggest_pop <- c(biggest_pop,7)
              }
              else{
                if(total_data_new_numbers$NL_AF[i]==number ){
                  biggest_pop <- c(biggest_pop,8)
                }
                else{
                  if(total_data_new_numbers$UK_AF[i]==number ){
                    biggest_pop <- c(biggest_pop,9)
                  }
                  else{
                    if(total_data_new_numbers$ExAC_SAS_MAF[i]==number ){
                      biggest_pop <- c(biggest_pop,10)
                    }
                    else{
                      if(total_data_new_numbers$AFR_MAF[i]==number ){
                        biggest_pop <- c(biggest_pop,11)
                      }
                      else{
                      if(total_data_new_numbers$AMR_MAF[i]==number ){
                        biggest_pop <- c(biggest_pop,12)
                      }
                      else{
                        if(total_data_new_numbers$EAS_MAF[i]==number ){
                          biggest_pop <- c(biggest_pop,13)
                        }
                        else{
                          if(total_data_new_numbers$EUR_MAF[i]==number ){
                            biggest_pop <- c(biggest_pop,14)
                          }
                          else{
                            if(total_data_new_numbers$SAS_MAF[i]==number ){
                              biggest_pop <- c(biggest_pop,15)
                            }
                    }
                }
              }
            }
          }
        }
      }  
    }
  }
}}
      }}}}
pca=prcomp(total_data_new_numbers,scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:2])
plot(newdat$PC1,newdat$PC2,type='n')
text(newdat$PC1,newdat$PC2, row.names(newdat), cex=0.6, pos=4, col="red")
#newdat$id <- total_data_interesting$refsnp_id
newdat$biggest <- biggest_pop
newdat$VIP = total_data_new$VIP
plot(newdat$PC1,newdat$PC2,yaxt='n',xaxt='n', col=ifelse(newdat$biggest==1,'lawngreen'
    ,ifelse(newdat$biggest==2,'red',ifelse(newdat$biggest==3
    ,'green',ifelse(newdat$biggest==4,'blue',ifelse(newdat$biggest==5,'purple',
    ifelse(newdat$biggest==6,'black',ifelse(newdat$biggest==7,'orange',ifelse(newdat$biggest==8,'brown',
        ifelse(newdat$biggest==9,'gray',ifelse(newdat$biggest==10,'hotpink',
       ifelse(newdat$biggest==11,'yellow',ifelse(newdat$biggest==12,'khaki',
       ifelse(newdat$biggest==13,'indianred',ifelse(newdat$biggest==14,'mediumseagreen','midnightblue'))))))))))))))
     ,pch=ifelse(newdat$VIP==0,15,17),cex= ifelse(newdat$VIP==0,1,2)
     ,xlab = "PC1",ylab="PC2")
par(mar=c(5.1,4.1,4.1,10.1), xpd=TRUE)
legend("topright",inset=c(-0.4,0),col=c("blue", "red","green","purple","lawngreen",
          "black","hotpink","orange","brown","gray","yellow","khaki","indianred","mediumseagreen","midnightblue"), pch=15,
       legend=c("ExAC_AFR_MAF","ExAC_AMR_MAF","ExAC_EAS_MAF", 
                "ExAC_FIN_MAF","ExAC_NFE_MAF","ExAC_OTH_MAF","ExAC_SAS_MAF"
                ,"EST_AF","NL_AF","UK_AF","AFR_MAF","AMR_MAF",
                "EAS_MAF","EUR_MAF","SAS_MAF"), bty="n", xjust=1,cex=1.3,x.intersp= 0.5,y.intersp= 0.6,seg.len=0.5)
#par(xpd = T, mar = par()$mar + c(0,0,0,7))
axis(1,cex.axis=1.5)
axis(2,cex.axis=1.5)

pca_fin <- subset(newdat,biggest==5)
pca_nl <- subset(newdat,biggest==8)
pca_exac_afr <- subset(newdat,biggest==11)
pca_exac_eur <- subset(newdat,biggest==1)
pca_exac_eas <- subset(newdat,biggest==3)
pca_exac_amr <- subset(newdat,biggest==2)
pca_exac_sas <- subset(newdat,biggest==10)
pca_est <- subset(newdat,biggest==7)
t.test(pca_exac_amr$PC2,pca_exac_eur$PC2,var.equal=TRUE, paired=FALSE)
### Leiame ekstreemsusi ka 1000g rahvastel

max_afr <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = max(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$AFR_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$afr_AF[i]-0.015 > total_data_nocode$afr_AF[i] & total_data_nocode$afr_AF[i]-0.015> total_data_nocode$afr_AF[i] & total_data_nocode$afr_AF[i] >0.9){
    max_afr <- rbind(max_afr,total_data_nocode[i,])
    #}
  }
}
max_afr <- max_afr[!duplicated(max_afr), ]


min_afr <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = min(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$AFR_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$afr_AF[i]-0.015 > total_data_nocode$afr_AF[i] & total_data_nocode$afr_AF[i]-0.015> total_data_nocode$afr_AF[i] & total_data_nocode$afr_AF[i] >0.9){
    min_afr <- rbind(min_afr,total_data_nocode[i,])
    #}
  }
}
min_afr <- min_afr[!duplicated(min_afr), ]


max_amr <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = max(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$AMR_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$amr_AF[i]-0.015 > total_data_nocode$amr_AF[i] & total_data_nocode$amr_AF[i]-0.015> total_data_nocode$amr_AF[i] & total_data_nocode$amr_AF[i] >0.9){
    max_amr <- rbind(max_amr,total_data_nocode[i,])
    #}
  }
}
max_amr <- max_amr[!duplicated(max_amr), ]


min_amr <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = min(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$AMR_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$amr_AF[i]-0.015 > total_data_nocode$amr_AF[i] & total_data_nocode$amr_AF[i]-0.015> total_data_nocode$amr_AF[i] & total_data_nocode$amr_AF[i] >0.9){
    min_amr <- rbind(min_amr,total_data_nocode[i,])
    #}
  }
}
min_amr <- min_amr[!duplicated(min_amr), ]


max_eas <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = max(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$EAS_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$eas_AF[i]-0.015 > total_data_nocode$eas_AF[i] & total_data_nocode$eas_AF[i]-0.015> total_data_nocode$eas_AF[i] & total_data_nocode$eas_AF[i] >0.9){
    max_eas <- rbind(max_eas,total_data_nocode[i,])
    #}
  }
}
max_eas <- max_eas[!duplicated(max_eas), ]



min_eas <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = min(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$EAS_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$eas_AF[i]-0.015 > total_data_nocode$eas_AF[i] & total_data_nocode$eas_AF[i]-0.015> total_data_nocode$eas_AF[i] & total_data_nocode$eas_AF[i] >0.9){
    min_eas <- rbind(min_eas,total_data_nocode[i,])
    #}
  }
}
min_eas <- min_eas[!duplicated(min_eas), ]



max_eur <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = max(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$EUR_MAF[i]==number && number != 0 ){
    #if(total_data_nocode$eur_AF[i]-0.015 > total_data_nocode$eur_AF[i] & total_data_nocode$eur_AF[i]-0.015> total_data_nocode$eur_AF[i] & total_data_nocode$eur_AF[i] >0.9){
    max_eur <- rbind(max_eur,total_data_nocode[i,])
    #}
  }
}
max_eur <- max_eur[!duplicated(max_eur), ]



min_eur <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = min(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$EUR_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$eur_AF[i]-0.015 > total_data_nocode$eur_AF[i] & total_data_nocode$eur_AF[i]-0.015> total_data_nocode$eur_AF[i] & total_data_nocode$eur_AF[i] >0.9){
    min_eur <- rbind(min_eur,total_data_nocode[i,])
    #}
  }
}
min_eur <- min_eur[!duplicated(min_eur), ]


max_sas <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = max(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$SAS_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$sas_AF[i]-0.015 > total_data_nocode$sas_AF[i] & total_data_nocode$sas_AF[i]-0.015> total_data_nocode$sas_AF[i] & total_data_nocode$sas_AF[i] >0.9){
    max_sas <- rbind(max_sas,total_data_nocode[i,])
    #}
  }
}
max_sas <- max_sas[!duplicated(max_sas), ]

min_sas <- data.frame()
for(i in 1:(length(total_data_nocode$UK_AF))){
  number = min(total_data_nocode$NL_AF[i],total_data_nocode$UK_AF[i],total_data_nocode$EST_AF[i],total_data_nocode$AFR_MAF[i],
               total_data_nocode$AMR_MAF[i],total_data_nocode$EAS_MAF[i],total_data_nocode$EUR_MAF[i],total_data_nocode$SAS_MAF[i]
  )
  if(total_data_nocode$SAS_MAF[i]==number  && number != 0 ){
    #if(total_data_nocode$sas_AF[i]-0.015 > total_data_nocode$sas_AF[i] & total_data_nocode$sas_AF[i]-0.015> total_data_nocode$sas_AF[i] & total_data_nocode$sas_AF[i] >0.9){
    min_sas <- rbind(min_sas,total_data_nocode[i,])
    #}
  }
}
min_sas <- min_sas[!duplicated(min_sas), ]

max_ea_maf <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$EA_MAF[i]==number  && number != 0 ){
    #if(total_data_new$ea_maf_AF[i]-0.015 > total_data_new$ea_maf_AF[i] & total_data_new$ea_maf_AF[i]-0.015> total_data_new$ea_maf_AF[i] & total_data_new$ea_maf_AF[i] >0.9){
    max_ea_maf <- rbind(max_ea_maf,total_data_new[i,])
    #}
  }
}
max_ea_maf <- max_ea_maf[!duplicated(max_ea_maf), ]


min_ea_maf <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$EA_MAF[i]==number ){
    #if(total_data_new$ea_maf_AF[i]-0.015 > total_data_new$ea_maf_AF[i] & total_data_new$ea_maf_AF[i]-0.015> total_data_new$ea_maf_AF[i] & total_data_new$ea_maf_AF[i] >0.9){
    min_ea_maf <- rbind(min_ea_maf,total_data_new[i,])
    #}
  }
}
min_ea_maf <- min_ea_maf[!duplicated(min_ea_maf), ]


max_aa_maf <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$AA_MAF[i]==number ){
    #if(total_data_new$aa_maf_AF[i]-0.015 > total_data_new$aa_maf_AF[i] & total_data_new$aa_maf_AF[i]-0.015> total_data_new$aa_maf_AF[i] & total_data_new$aa_maf_AF[i] >0.9){
    max_aa_maf <- rbind(max_aa_maf,total_data_new[i,])
    #}
  }
}
max_aa_maf <- max_aa_maf[!duplicated(max_aa_maf), ]


min_aa_maf <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$AA_MAF[i]==number ){
    #if(total_data_new$aa_maf_AF[i]-0.015 > total_data_new$aa_maf_AF[i] & total_data_new$aa_maf_AF[i]-0.015> total_data_new$aa_maf_AF[i] & total_data_new$aa_maf_AF[i] >0.9){
    min_aa_maf <- rbind(min_aa_maf,total_data_new[i,])
    #}
  }
}
min_aa_maf <- min_aa_maf[!duplicated(min_aa_maf), ]


max_exac_amr <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_AMR_MAF[i]==number ){
    #if(total_data_new$exac_amr_AF[i]-0.015 > total_data_new$exac_amr_AF[i] & total_data_new$exac_amr_AF[i]-0.015> total_data_new$exac_amr_AF[i] & total_data_new$exac_amr_AF[i] >0.9){
    max_exac_amr <- rbind(max_exac_amr,total_data_new[i,])
    #}
  }
}
max_exac_amr <- max_exac_amr[!duplicated(max_exac_amr), ]


min_exac_amr <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_AMR_MAF[i]==number ){
    #if(total_data_new$exac_amr_AF[i]-0.015 > total_data_new$exac_amr_AF[i] & total_data_new$exac_amr_AF[i]-0.015> total_data_new$exac_amr_AF[i] & total_data_new$exac_amr_AF[i] >0.9){
    min_exac_amr <- rbind(min_exac_amr,total_data_new[i,])
    #}
  }
}
min_exac_amr <- min_exac_amr[!duplicated(min_exac_amr), ]


max_exac_afr <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_AFR_MAF[i]==number ){
    #if(total_data_new$exac_afr_AF[i]-0.015 > total_data_new$exac_afr_AF[i] & total_data_new$exac_afr_AF[i]-0.015> total_data_new$exac_afr_AF[i] & total_data_new$exac_afr_AF[i] >0.9){
    max_exac_afr <- rbind(max_exac_afr,total_data_new[i,])
    #}
  }
}
max_exac_afr <- max_exac_afr[!duplicated(max_exac_afr), ]


min_exac_afr <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_AFR_MAF[i]==number ){
    #if(total_data_new$exac_afr_AF[i]-0.015 > total_data_new$exac_afr_AF[i] & total_data_new$exac_afr_AF[i]-0.015> total_data_new$exac_afr_AF[i] & total_data_new$exac_afr_AF[i] >0.9){
    min_exac_afr <- rbind(min_exac_afr,total_data_new[i,])
    #}
  }
}
min_exac_afr <- min_exac_afr[!duplicated(min_exac_afr), ]


max_exac_eas <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_EAS_MAF[i]==number ){
    #if(total_data_new$exac_eas_AF[i]-0.015 > total_data_new$exac_eas_AF[i] & total_data_new$exac_eas_AF[i]-0.015> total_data_new$exac_eas_AF[i] & total_data_new$exac_eas_AF[i] >0.9){
    max_exac_eas <- rbind(max_exac_eas,total_data_new[i,])
    #}
  }
}
max_exac_eas <- max_exac_eas[!duplicated(max_exac_eas), ]


min_exac_eas <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_EAS_MAF[i]==number ){
    #if(total_data_new$exac_eas_AF[i]-0.015 > total_data_new$exac_eas_AF[i] & total_data_new$exac_eas_AF[i]-0.015> total_data_new$exac_eas_AF[i] & total_data_new$exac_eas_AF[i] >0.9){
    min_exac_eas <- rbind(min_exac_eas,total_data_new[i,])
    #}
  }
}
min_exac_eas <- min_exac_eas[!duplicated(min_exac_eas), ]



max_exac_fin <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_FIN_MAF[i]==number ){
    #if(total_data_new$exac_fin_AF[i]-0.015 > total_data_new$exac_fin_AF[i] & total_data_new$exac_fin_AF[i]-0.015> total_data_new$exac_fin_AF[i] & total_data_new$exac_fin_AF[i] >0.9){
    max_exac_fin <- rbind(max_exac_fin,total_data_new[i,])
    #}
  }
}
max_exac_fin <- max_exac_fin[!duplicated(max_exac_fin), ]


min_exac_fin <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_FIN_MAF[i]==number ){
    #if(total_data_new$exac_fin_AF[i]-0.015 > total_data_new$exac_fin_AF[i] & total_data_new$exac_fin_AF[i]-0.015> total_data_new$exac_fin_AF[i] & total_data_new$exac_fin_AF[i] >0.9){
    min_exac_fin <- rbind(min_exac_fin,total_data_new[i,])
    #}
  }
}
min_exac_fin <- min_exac_fin[!duplicated(min_exac_fin), ]


max_exac_nfe <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_NFE_MAF[i]==number ){
    #if(total_data_new$exac_nfe_AF[i]-0.015 > total_data_new$exac_nfe_AF[i] & total_data_new$exac_nfe_AF[i]-0.015> total_data_new$exac_nfe_AF[i] & total_data_new$exac_nfe_AF[i] >0.9){
    max_exac_nfe <- rbind(max_exac_nfe,total_data_new[i,])
    #}
  }
}
max_exac_nfe <- max_exac_nfe[!duplicated(max_exac_nfe), ]


min_exac_nfe <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_NFE_MAF[i]==number ){
    #if(total_data_new$exac_nfe_AF[i]-0.015 > total_data_new$exac_nfe_AF[i] & total_data_new$exac_nfe_AF[i]-0.015> total_data_new$exac_nfe_AF[i] & total_data_new$exac_nfe_AF[i] >0.9){
    min_exac_nfe <- rbind(min_exac_nfe,total_data_new[i,])
    #}
  }
}
min_exac_nfe <- min_exac_nfe[!duplicated(min_exac_nfe), ]



max_exac_oth <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_OTH_MAF[i]==number ){
    #if(total_data_new$exac_oth_AF[i]-0.015 > total_data_new$exac_oth_AF[i] & total_data_new$exac_oth_AF[i]-0.015> total_data_new$exac_oth_AF[i] & total_data_new$exac_oth_AF[i] >0.9){
    max_exac_oth <- rbind(max_exac_oth,total_data_new[i,])
    #}
  }
}
max_exac_oth <- max_exac_oth[!duplicated(max_exac_oth), ]


min_exac_oth <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_OTH_MAF[i]==number ){
    #if(total_data_new$exac_oth_AF[i]-0.015 > total_data_new$exac_oth_AF[i] & total_data_new$exac_oth_AF[i]-0.015> total_data_new$exac_oth_AF[i] & total_data_new$exac_oth_AF[i] >0.9){
    min_exac_oth <- rbind(min_exac_oth,total_data_new[i,])
    #}
  }
}
min_exac_oth <- min_exac_oth[!duplicated(min_exac_oth), ]



max_exac_sas <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = max(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_SAS_MAF[i]==number ){
    #if(total_data_new$exac_sas_AF[i]-0.015 > total_data_new$exac_sas_AF[i] & total_data_new$exac_sas_AF[i]-0.015> total_data_new$exac_sas_AF[i] & total_data_new$exac_sas_AF[i] >0.9){
    max_exac_sas <- rbind(max_exac_sas,total_data_new[i,])
    #}
  }
}
max_exac_sas <- max_exac_sas[!duplicated(max_exac_sas), ]


min_exac_sas <- data.frame()
for(i in 1:(length(total_data_new$UK_AF))){
  number = min(total_data_new$NL_AF[i],total_data_new$UK_AF[i],total_data_new$EST_AF[i],total_data_new$AFR_MAF[i],
               total_data_new$AMR_MAF[i],total_data_new$EAS_MAF[i],total_data_new$EUR_MAF[i],total_data_new$SAS_MAF[i]
               ,total_data_new$AA_MAF[i],total_data_new$EA_MAF[i],total_data_new$ExAC_AMR_MAF[i],total_data_new$ExAC_AFR_MAF[i],total_data_new$ExAC_EAS_MAF[i],
               total_data_new$ExAC_FIN_MAF[i],total_data_new$ExAC_NFE_MAF[i],total_data_new$ExAC_OTH_MAF[i],total_data_new$ExAC_SAS_MAF[i],total_data_new$ExAC_SAS_MAF[i]
  )
  if(total_data_new$ExAC_SAS_MAF[i]==number ){
    #if(total_data_new$exac_sas_AF[i]-0.015 > total_data_new$exac_sas_AF[i] & total_data_new$exac_sas_AF[i]-0.015> total_data_new$exac_sas_AF[i] & total_data_new$exac_sas_AF[i] >0.9){
    min_exac_sas <- rbind(min_exac_sas,total_data_new[i,])
    #}
  }
}
min_exac_sas <- min_exac_sas[!duplicated(min_exac_sas), ]



## Kirjutame df millest saame histogrammi luua
total_data_new_numbers$ExAC_Adj_MAF <- NULL
extreme_pop_tabel <- rbind(colnames(total_data_new_numbers)[1:15])

colnames(extreme_pop_tabel) <- colnames(total_data_new_numbers)[1:15]

extreme_pop_tabel_nocode <- subset(extreme_pop_tabel,select= c(NL_AF:SAS_MAF))
extreme_pop_tabel_nocode <- t(extreme_pop_tabel_nocode)
extreme_pop_tabel <- t(extreme_pop_tabel)
extreme_pop_tabel <- as.data.frame(extreme_pop_tabel)
extreme_pop_tabel$AA_MAF <- NULL
extreme_pop_tabel <- cbind(extreme_pop_tabel,extreme_pop_tabel[,1])
extreme_pop_tabel_nocode <- cbind(extreme_pop_tabel_nocode,extreme_pop_tabel_nocode[,1])
colnames(extreme_pop_tabel) <- c('max','min')
colnames(extreme_pop_tabel_nocode) <- c('max','min')
extreme_pop_tabel_nocode[,1] = c(382,455,204,1363,404,1197,236,278)
extreme_pop_tabel_nocode[,2] = c(359,466,143,1684,289,1500,138,492) 
extreme_pop_tabel[,1] = c(42,36,31,145,6,69,27,38,5,6,11,66,79,75,12,3,23)
extreme_pop_tabel[,2] = c(25,28,15,186,18,126,13,37,117,118,5,31,127,37,7,9,14)

barplot(prop.table(table(t(extreme_pop_tabel))))
  plot(extreme_pop_tabel[,1],extreme_pop_tabel[,2],type='n')
extreme_pop_tabel$max <- as.numeric(extreme_pop_tabel$max)
extreme_pop_tabel$min <- as.numeric(extreme_pop_tabel$min)
extreme_pop_tabel_nocode$max <- as.numeric(extreme_pop_tabel_nocode$max)
extreme_pop_tabel_nocode$min <- as.numeric(extreme_pop_tabel_nocode$min)

extreme_pop_tabel <- as.data.frame(extreme_pop_tabel)
extreme_pop_tabel_nocode <- as.data.frame(extreme_pop_tabel_nocode)
extreme_pop_tabel_nocode[,1] <- as.numeric(extreme_pop_tabel_nocode[,1])
extreme_pop_tabel_nocode[,2] <- as.numeric(extreme_pop_tabel_nocode[,2])

hist(extreme_pop_tabel[,1])
extreme_pop_tabel.plot.bar()
extreme_pop_tabel.hist()
pop_df = pd.Dataframe()
install.packages("ggplot2")
library(ggplot2)

qplot(extreme_pop_tabel,
      geom="histogram",
      binwidth = 0.5,  
      main = "Histogram for Age", 
      xlab = "Age",  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2),
      xlim=c(20,50))

op <- par(mar = c(10,4,4,2) + 0.1)
par(las=2) # make label text perpendicular to axis
par(mar=c(5,8,4,2)) # increase y-axis margin.
revdf <- extreme_pop_tabel[rev(rownames(extreme_pop_tabel)),]
barplot(as.matrix(t(revdf)),las=2,beside=TRUE,horiz=TRUE,col=c('red','blue'))
legend("topright", c("Maksimaalne","Minimaalne"), cex=1.3, bty="n", fill=c('red','blue'))
barplot(as.matrix(t(revdf)),las=2,horiz=TRUE,beside=TRUE,col=c('red','blue'))

par(op) ## reset
## Kirjutada protsentuaalsete erinevused sisse (eestlastele)
#ja järjestada neid , ning tuua välja need mis on kõige suurema protsendiga ja VIP geenid ka 
total_data_3and1000 <- subset(total_data_new,select = c(refsnp_id:SAS_MAF))
min_est_test <- min_est
min_est_test$EST_AF <- NULL
min_est_test$refsnp_id <- NULL
min_est_test$VIP <- NULL
min_est
min_est$protsent <- min_est_test$NL_AF
for(i in 1:(length(min_est$EST_AF))){
  number = min(min_est_test[i,])
  min_est$protsent[i] <- ((number-min_est$EST_AF[i])/number)
}
min_est <- subset(min_est,protsent < 0.9)


max_est_test <- max_est
max_est_test$EST_AF <- NULL
max_est_test$refsnp_id <- NULL
max_est_test$VIP <- NULL
max_est
max_est$protsent <- max_est_test$NL_AF
for(i in 1:(length(max_est$EST_AF))){
  number = max(max_est_test[i,])
  max_est$protsent[i] <- ((number-max_est$EST_AF[i])/number)
}
max_est <- subset(max_est,protsent < 0.9)
# Teostame t-testid

# Kodeeritavad variandid kus pole VIP variandid

coded_data <- subset(total_data_new,VIP < 1)
coded_data$refsnp_id <- NULL
coded_data$VIP <- NULL

pca=prcomp(coded_data,scale=T)
summary(pca) 
newdat<-data.frame(pca$x[,1:8])

# seal kus on VIP andmed 

coded_data_vip <- subset(total_data_new,VIP > 0)
coded_data_vip$refsnp_id <- NULL
coded_data_vip$VIP <- NULL

pca_vip=prcomp(coded_data_vip,scale=T)
summary(pca_vip) 
newdat_vip<-data.frame(pca_vip$x[,1:17])

t.test(newdat$PC1,newdat_vip$PC4,var.equal=TRUE, paired=FALSE)

t.test(newdat$PC2,newdat_vip$PC2,var.equal=TRUE, paired=FALSE)

## Umbes samad arvud kuid mitte-kodeeritavatele


none_coded_data <- subset(total_data_nocode,VIP < 1)
none_coded_data <- subset(none_coded_data,select=c('refsnp_id','NL_AF','UK_AF','EST_AF',
                                                           'AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF'))
none_coded_data$refsnp_id <- NULL
none_coded_data$VIP <- NULL


pca_nocode=prcomp(none_coded_data,scale=T)

summary(pca) 
newdat_nocode<-data.frame(pca_nocode$x[,1:8])

# seal kus on VIP andmed 

none_coded_data_vip <- subset(total_data_nocode,VIP > 0)
none_coded_data_vip <- subset(none_coded_data_vip,select=c('refsnp_id','NL_AF','UK_AF','EST_AF',
                                                           'AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF'))
none_coded_data_vip$refsnp_id <- NULL

pca_nocode_vip=prcomp(none_coded_data_vip,scale=T)
summary(pca_vip) 
newdat_nocode_vip<-data.frame(pca_nocode_vip$x[,1:8])


t.test(newdat$PC1,newdat_nocode$PC1,var.equal=TRUE, paired=FALSE)

t.test(newdat_vip$PC1,newdat_nocode_vip$PC1,var.equal=TRUE, paired=FALSE)

###JÄRELDUS = ERINEVUST EI OLE KUNA KÕIK p-value ARVUD ON NULLID


### TEOSTAME Fst ARVUTUSE
# LISAN VAJALIKUD RAHVAARVUD

pop_num <- data.frame(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
colnames(pop_num) <- c('NL_AF','UK_AF','EST_AF',colnames(MAF_list))
pop_num$ExAC_MAF <- NULL
pop_num$AA_MAF <- NULL
pop_num$EA_MAF <- NULL
pop_num$ExAC_Adj_MAF <- NULL
pop_num <- as.data.frame(pop_num)
pop_num[1,] <- c(2244,2244,2244,1018,535,617,669,661,5203,5789,4327,3307,33370,8256,454)
pop_sum <- sum(pop_num[1,])


#Teeme nüüd doki 2.osa

total_data_fst <- total_data_new
total_data_fst$NL_AF_Hexp <- 1
total_data_fst$UK_AF_Hexp <- 1
total_data_fst$EST_AF_Hexp <- 1
total_data_fst$AFR_MAF_Hexp <- 1
total_data_fst$AMR_MAF_Hexp <- 1
total_data_fst$EAS_MAF_Hexp <- 1
total_data_fst$EUR_MAF_Hexp <- 1
total_data_fst$SAS_MAF_Hexp <- 1
total_data_fst$AA_MAF <- NULL
total_data_fst$EA_MAF <- NULL
total_data_fst$ExAC_AFR_MAF_Hexp <- 1
total_data_fst$ExAC_AMR_MAF_Hexp <- 1
total_data_fst$ExAC_EAS_MAF_Hexp <- 1
total_data_fst$ExAC_FIN_MAF_Hexp <- 1
total_data_fst$ExAC_NFE_MAF_Hexp <- 1
total_data_fst$ExAC_OTH_MAF_Hexp <- 1
total_data_fst$ExAC_SAS_MAF_Hexp <- 1
i = 0
for(i in 1:(length(total_data_fst$UK_AF))){
  total_data_fst$NL_AF_Hexp[i] <- (as.numeric(total_data_fst$NL_AF[i])^2 + as.numeric(1-total_data_fst$NL_AF[i])^2 )
  total_data_fst$NL_AF_Hexp[i] <- 1-total_data_fst$NL_AF_Hexp[i]
  total_data_fst$UK_AF_Hexp[i] <- ((total_data_fst$UK_AF[i])^2 + (1-total_data_fst$UK_AF[i])^2 )
  total_data_fst$UK_AF_Hexp[i] <- 1-total_data_fst$UK_AF_Hexp[i]
  total_data_fst$EST_AF_Hexp[i] <- ((total_data_fst$EST_AF[i])^2 + (1-total_data_fst$EST_AF[i])^2 )
  total_data_fst$EST_AF_Hexp[i] <- 1-total_data_fst$EST_AF_Hexp[i]
  total_data_fst$AFR_MAF_Hexp[i] <- ((total_data_fst$AFR_MAF[i])^2 + (1-total_data_fst$AFR_MAF[i])^2 )
  total_data_fst$AFR_MAF_Hexp[i] <- 1-total_data_fst$AFR_MAF_Hexp[i]
  total_data_fst$AMR_MAF_Hexp[i] <- ((total_data_fst$AMR_MAF[i])^2 + (1-total_data_fst$AMR_MAF[i])^2 )
  total_data_fst$AMR_MAF_Hexp[i] <- 1-total_data_fst$AMR_MAF_Hexp[i]
  total_data_fst$EAS_MAF_Hexp[i] <- ((total_data_fst$EAS_MAF[i])^2 + (1-total_data_fst$EAS_MAF[i])^2 )
  total_data_fst$EAS_MAF_Hexp[i] <- 1-total_data_fst$EAS_MAF_Hexp[i]
  total_data_fst$EUR_MAF_Hexp[i] <- ((total_data_fst$EUR_MAF[i])^2 + (1-total_data_fst$EUR_MAF[i])^2 )
  total_data_fst$EUR_MAF_Hexp[i] <- 1-total_data_fst$EUR_MAF_Hexp[i]
  total_data_fst$SAS_MAF_Hexp[i] <- ((total_data_fst$SAS_MAF[i])^2 + (1-total_data_fst$SAS_MAF[i])^2 )
  total_data_fst$SAS_MAF_Hexp[i] <- 1-total_data_fst$SAS_MAF_Hexp[i]
  total_data_fst$ExAC_AFR_MAF_Hexp[i] <- ((total_data_fst$ExAC_AFR_MAF[i])^2 + (1-total_data_fst$ExAC_AFR_MAF[i])^2 )
  total_data_fst$ExAC_AFR_MAF_Hexp[i] <- 1-total_data_fst$ExAC_AFR_MAF_Hexp[i]
  total_data_fst$ExAC_AMR_MAF_Hexp[i] <- ((total_data_fst$ExAC_AMR_MAF[i])^2 + (1-total_data_fst$ExAC_AMR_MAF[i])^2 )
  total_data_fst$ExAC_AMR_MAF_Hexp[i] <- 1-total_data_fst$ExAC_AMR_MAF_Hexp[i]
  total_data_fst$ExAC_EAS_MAF_Hexp[i] <- ((total_data_fst$ExAC_EAS_MAF[i])^2 + (1-total_data_fst$ExAC_EAS_MAF[i])^2 )
  total_data_fst$ExAC_EAS_MAF_Hexp[i] <- 1-total_data_fst$ExAC_EAS_MAF_Hexp[i]
  total_data_fst$ExAC_FIN_MAF_Hexp[i] <- ((total_data_fst$ExAC_FIN_MAF[i])^2 + (1-total_data_fst$ExAC_FIN_MAF[i])^2 )
  total_data_fst$ExAC_FIN_MAF_Hexp[i] <- 1-total_data_fst$ExAC_FIN_MAF_Hexp[i]
  total_data_fst$ExAC_NFE_MAF_Hexp[i] <- ((total_data_fst$ExAC_NFE_MAF[i])^2 + (1-total_data_fst$ExAC_NFE_MAF[i])^2 )
  total_data_fst$ExAC_NFE_MAF_Hexp[i] <- 1-total_data_fst$ExAC_NFE_MAF_Hexp[i]
  total_data_fst$ExAC_OTH_MAF_Hexp[i] <- ((total_data_fst$ExAC_OTH_MAF[i])^2 + (1-total_data_fst$ExAC_OTH_MAF[i])^2 )
  total_data_fst$ExAC_OTH_MAF_Hexp[i] <- 1-total_data_fst$ExAC_OTH_MAF_Hexp[i]
  total_data_fst$ExAC_SAS_MAF_Hexp[i] <- ((total_data_fst$ExAC_SAS_MAF[i])^2 + (1-total_data_fst$ExAC_SAS_MAF[i])^2 )
  total_data_fst$ExAC_SAS_MAF_Hexp[i] <- 1-total_data_fst$ExAC_SAS_MAF_Hexp[i]
  

  
}


# Arvutame kogusummad ( Hs ), doki 3.osa

i = 0
total_data_fst$summad <- 1
number <- total_data_fst[i,][17]*pop_num[1,][1]
for(i in 1:(length(total_data_fst$UK_AF))){
  number <- 0
  for(j in 1:(length(total_data_fst[1,][17:31]))){
    number <- number + as.numeric(total_data_fst[i,][16+j]) * as.numeric(pop_num[1,][j])
  }
  total_data_fst$summad[i] <- number/pop_sum
}




# Arvutame keskmise sageduse, doki 4.osa esimene pool et Ht vajab arvutamist

total_data_fst$freq_avr <- 1

for(i in 1:(length(total_data_fst$UK_AF))){
  number <- 0
  for(j in 1:(length(total_data_fst[1,][17:31]))){
    number <- number + as.numeric(total_data_fst[i,][1+j]) * as.numeric(pop_num[1,][j])
  }
  total_data_fst$freq_avr[i] <- number/pop_sum
}

i = 0
total_data_fst$Ht <- 1
for(i in 1:(length(total_data_fst$UK_AF))){
  number <- 0
  number <- as.numeric(total_data_fst$freq_avr[i])^2 + (1-as.numeric(total_data_fst$freq_avr[i]))^2 
  
  total_data_fst$Ht[i] <- 1- number
}

# Arvutame Fst ehk viimane osa

total_data_fst$Fst <- 1

for(i in 1:(length(total_data_fst$UK_AF))){
  total_data_fst$Fst[i] <-(total_data_fst$Ht[i] - total_data_fst$summad[i])/total_data_fst$Ht[i]
}

## SIIN arvutame mitte-kodeeritavate rahvaste Fst

total_data_nocode_fst <- subset(total_data_nocode, select=c('refsnp_id','NL_AF','UK_AF','EST_AF','AFR_MAF','AMR_MAF','EAS_MAF','EUR_MAF','SAS_MAF'))
total_data_nocode_fst$NL_AF_Hexp <- 1
total_data_nocode_fst$UK_AF_Hexp <- 1
total_data_nocode_fst$EST_AF_Hexp <- 1
total_data_nocode_fst$AFR_MAF_Hexp <- 1
total_data_nocode_fst$AMR_MAF_Hexp <- 1
total_data_nocode_fst$EAS_MAF_Hexp <- 1
total_data_nocode_fst$EUR_MAF_Hexp <- 1
total_data_nocode_fst$SAS_MAF_Hexp <- 1
i = 0

pop_sum <- sum(pop_num[1,][1:8] )
for(i in 1:(length(total_data_nocode_fst$UK_AF))){
  total_data_nocode_fst$NL_AF_Hexp[i] <- (as.numeric(total_data_nocode_fst$NL_AF[i])^2 + as.numeric(1-total_data_nocode_fst$NL_AF[i])^2 )
  total_data_nocode_fst$NL_AF_Hexp[i] <- 1-total_data_nocode_fst$NL_AF_Hexp[i]
  total_data_nocode_fst$UK_AF_Hexp[i] <- ((total_data_nocode_fst$UK_AF[i])^2 + (1-total_data_nocode_fst$UK_AF[i])^2 )
  total_data_nocode_fst$UK_AF_Hexp[i] <- 1-total_data_nocode_fst$UK_AF_Hexp[i]
  total_data_nocode_fst$EST_AF_Hexp[i] <- ((total_data_nocode_fst$EST_AF[i])^2 + (1-total_data_nocode_fst$EST_AF[i])^2 )
  total_data_nocode_fst$EST_AF_Hexp[i] <- 1-total_data_nocode_fst$EST_AF_Hexp[i]
  total_data_nocode_fst$AFR_MAF_Hexp[i] <- ((total_data_nocode_fst$AFR_MAF[i])^2 + (1-total_data_nocode_fst$AFR_MAF[i])^2 )
  total_data_nocode_fst$AFR_MAF_Hexp[i] <- 1-total_data_nocode_fst$AFR_MAF_Hexp[i]
  total_data_nocode_fst$AMR_MAF_Hexp[i] <- ((total_data_nocode_fst$AMR_MAF[i])^2 + (1-total_data_nocode_fst$AMR_MAF[i])^2 )
  total_data_nocode_fst$AMR_MAF_Hexp[i] <- 1-total_data_nocode_fst$AMR_MAF_Hexp[i]
  total_data_nocode_fst$EAS_MAF_Hexp[i] <- ((total_data_nocode_fst$EAS_MAF[i])^2 + (1-total_data_nocode_fst$EAS_MAF[i])^2 )
  total_data_nocode_fst$EAS_MAF_Hexp[i] <- 1-total_data_nocode_fst$EAS_MAF_Hexp[i]
  total_data_nocode_fst$EUR_MAF_Hexp[i] <- ((total_data_nocode_fst$EUR_MAF[i])^2 + (1-total_data_nocode_fst$EUR_MAF[i])^2 )
  total_data_nocode_fst$EUR_MAF_Hexp[i] <- 1-total_data_nocode_fst$EUR_MAF_Hexp[i]
  total_data_nocode_fst$SAS_MAF_Hexp[i] <- ((total_data_nocode_fst$SAS_MAF[i])^2 + (1-total_data_nocode_fst$SAS_MAF[i])^2 )
  total_data_nocode_fst$SAS_MAF_Hexp[i] <- 1-total_data_nocode_fst$SAS_MAF_Hexp[i]

  
}


# Arvutame kogusummad ( Hs ), doki 3.osa

i = 0
total_data_nocode_fst$summad <- 1
number <- total_data_nocode_fst[i,][17]*pop_num[1,][1]
for(i in 1:(length(total_data_nocode_fst$UK_AF))){
  number <- 0
  for(j in 1:(length(total_data_nocode_fst[1,][10:17]))){
    number <- number + as.numeric(total_data_nocode_fst[i,][9+j]) * as.numeric(pop_num[1,][j])
  }
  total_data_nocode_fst$summad[i] <- number/pop_sum
}




# Arvutame keskmise sageduse, doki 4.osa esimene pool et Ht vajab arvutamist

total_data_nocode_fst$freq_avr <- 1

for(i in 1:(length(total_data_nocode_fst$UK_AF))){
  number <- 0
  for(j in 1:(length(total_data_nocode_fst[1,][10:17]))){
    number <- number + as.numeric(total_data_nocode_fst[i,][1+j]) * as.numeric(pop_num[1,][j])
  }
  total_data_nocode_fst$freq_avr[i] <- number/pop_sum
}

i = 0
total_data_nocode_fst$Ht <- 1
for(i in 1:(length(total_data_nocode_fst$UK_AF))){
  number <- 0
  number <- as.numeric(total_data_nocode_fst$freq_avr[i])^2 + (1-as.numeric(total_data_nocode_fst$freq_avr[i]))^2 
  
  total_data_nocode_fst$Ht[i] <- 1-number
}

# Arvutame Fst ehk viimane osa

total_data_nocode_fst$Fst <- 1

for(i in 1:(length(total_data_nocode_fst$UK_AF))){
  total_data_nocode_fst$Fst[i] <-(total_data_nocode_fst$Ht[i] - total_data_nocode_fst$summad[i])/total_data_nocode_fst$Ht[i]
}


# ARVUTAME EESTLASTE JA TEISTE RAHVASTE VAHE

est_Fst <- total_data_fst
est_Fst$summad <- NULL
est_Fst$freq_avr <- NULL
est_Fst$Ht <- NULL
est_Fst$Fst <- NULL

est_pop_num <- pop_num
est_pop_num$EST_AF <- NULL
for(j in 1:(length(est_pop_num[1,][1:14]))){
  est_pop_num[1,][j] <- est_pop_num[1,][j]+2244
}

i = 0
est_Fst$summad <- 1


# Arvutame kogusummad ( Hs ), doki 3.osa

i = 1
est_Fst$summad <- 1
est_Fst$freq_avr <- 1
est_Fst$Ht <- 1

num_of_pop <- 3
popsum <- as.numeric(est_pop_num[,num_of_pop][1])
for(i in 1:(length(est_Fst$UK_AF))){
  number <- 0
  number <- number + as.numeric(est_Fst[i,][19])*2244 + (as.numeric(est_Fst[i,][17+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
  est_Fst$summad[i] <- number/popsum
  number <- 0
  number <- number + as.numeric(est_Fst[i,][4])*2244 + (as.numeric(est_Fst[i,][2+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
  est_Fst$freq_avr[i] <- number/popsum
  number <- 0
  number <- as.numeric(est_Fst$freq_avr[i])^2 + (1-as.numeric(est_Fst$freq_avr[i]))^2 
  est_Fst$Ht[i] <- 1- number
  est_Fst$Fst[i] <-(as.numeric(est_Fst$Ht[i]) - as.numeric(est_Fst$summad[i]))/as.numeric(est_Fst$Ht[i])
}
## Rahvas : 
colnames(est_pop_num)[num_of_pop]
## Keskmine :
sum(est_Fst$Fst)/1276

## ARVUTAME HOLLANDLASTE FST

nl_Fst <- total_data_fst
nl_Fst$summad <- NULL
nl_Fst$freq_avr <- NULL
nl_Fst$Ht <- NULL
nl_Fst$Fst <- NULL

nl_pop_num <- pop_num
nl_pop_num$NL_AF <- NULL
for(j in 1:(length(nl_pop_num[1,][1:14]))){
  nl_pop_num[1,][j] <- nl_pop_num[1,][j]+2244
}

i = 0
nl_Fst$summad <- 1


# Arvutame kogusummad ( Hs ), doki 3.osa
asi <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

suur_tabel <- data.frame(asi,asi,asi,asi,asi,asi,asi,asi,asi,asi,asi,asi,asi,asi,asi)
colnames(suur_tabel) <- c('NL_AF',colnames(nl_pop_num))
row.names(suur_tabel) <- c('NL_AF',colnames(nl_pop_num))
i = 1
nl_Fst$summad 
i = 1
nl_Fst$summad <- 1
nl_Fst$freq_avr <- 1
nl_Fst$Ht <- 1
nk_Fst$Fst <- 1
nl_pop_num[2,] <- 1 
nl_Fst$freq_avr <- 1
nl_Fst$Ht <- 1
nl_pop_num[2,] <- 1
for(g in 1:(length(uk_pop_num))){
  
num_of_pop <- g
popsum <- as.numeric(uk_pop_num[,num_of_pop][1])
for(i in 1:(length(nl_Fst$nl_AF))){
  number <- 0
  number <- number + as.numeric(nl_Fst[i,][18])*2244 + (as.numeric(nl_Fst[i,][17+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
  nl_Fst$summad[i] <- number/popsum
  number <- 0
  number <- number + as.numeric(nl_Fst[i,][3])*2244 + (as.numeric(nl_Fst[i,][2+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
  nl_Fst$freq_avr[i] <- number/popsum
  number <- 0
  number <- as.numeric(nl_Fst$freq_avr[i])^2 + (1-as.numeric(nl_Fst$freq_avr[i]))^2 
  nl_Fst$Ht[i] <- 1- number
  nl_Fst$Fst[i] <-(as.numeric(nl_Fst$Ht[i]) - as.numeric(nl_Fst$summad[i]))/as.numeric(nl_Fst$Ht[i])
}
## Rahvas : 
colnames(uk_pop_num)[num_of_pop]
## Keskmine :
suur_tabel$UK_AF[1+g] <- sum(nl_Fst$Fst, na.rm=TRUE)/1276
}
suur_tabel$UK_AF <- 1
nl_pop_num <- pop_num
nl_pop_num$NL_AF <- NULL
for(j in 1:(length(nl_pop_num[1,][1:14]))){
  nl_pop_num[1,][j] <- nl_pop_num[1,][j]+2244
}
uk_Fst <- total_data_fst
uk_Fst$summad <- 1
uk_Fst$freq_avr <- 1
uk_Fst$Ht <- 1
uk_Fst$Fst <- 1
uk_pop_num <- pop_num
uk_pop_num$UK_AF <- NULL
for(j in 1:(length(uk_pop_num[1,][1:14]))){
  uk_pop_num[1,][j] <- uk_pop_num[1,][j]+2244
}
for(g in 1:(length(uk_pop_num))){
  
  num_of_pop <- g
  popsum <- as.numeric(uk_pop_num[,num_of_pop][1])
  for(i in 1:(length(uk_Fst$uk_AF))){
    number <- 0
    number <- number + as.numeric(uk_Fst[i,][18])*2244 + (as.numeric(uk_Fst[i,][17+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
    uk_Fst$summad[i] <- number/popsum
    number <- 0
    number <- number + as.numeric(uk_Fst[i,][3])*2244 + (as.numeric(uk_Fst[i,][2+num_of_pop]) * as.numeric(pop_num[,num_of_pop+1]))
    uk_Fst$freq_avr[i] <- number/popsum
    number <- 0
    number <- as.numeric(uk_Fst$freq_avr[i])^2 + (1-as.numeric(uk_Fst$freq_avr[i]))^2 
    uk_Fst$Ht[i] <- 1- number
    uk_Fst$Fst[i] <-(as.numeric(uk_Fst$Ht[i]) - as.numeric(uk_Fst$summad[i]))/as.numeric(uk_Fst$Ht[i])
  }
  ## Rahvas : 
  colnames(uk_pop_num)[num_of_pop]
  ## Keskmine :
  suur_tabel$UK_AF[1+g] <- sum(uk_Fst$Fst, na.rm=TRUE)/1276
}

## proovime heatmapi
fst <- t(fst)
nimed <- colnames(fst)
colnames(fst) <- nimed
revdf <- fst[rev(rownames(fst)),]
proov <- data.matrix(fst)
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
color = rev(cm.colors(300))
proov_heatmap <- heatmap(fst_uus, Rowv=NA, Colv= NA, col = my_palette, scale="none",
                         margins=c(10,10),cexRow = 1.4,cexCol = 1.4)
colfunc <- colorRampPalette(c("red", "yellow", "green"))
colfunc(20)
plot(rep(1,20),col=colfunc(20),pch=15,cex=3,xaxt = "n",yaxt="n")
numbrid = seq(0.05, 1.0, by = 0.05)
axis(1, at=1:20, labels=numbrid)
xlabel("")
kolmene_tabel <- subset(kolmene_tabel, select=c(NL_AF, UK_AF,EST_AF))
kolmene_tabel <- t(kolmene_tabel)

proov_heatmap <- heatmap(kolmene_tabel, Rowv=NA, Colv=NA, col = my_palette, scale="row", margins=c(10,10))
### Pole ammu kirjutanud lol
fst_uus <- fst
fst_uus <- as.data.frame(fst_uus)
fst_uus$NL_AF <- NULL
fst_uus$UK_AF <- NULL
fst_uus[2,]
fst_uus <- t(fst_uus)
NLGoPharmacoVariants <- NLGoPharmacoVariants[!duplicated(NLGoPharmacoVariants), ]
EstonianPharmacoVariantsCol8 <- EstonianPharmacoVariantsCol8[!duplicated(EstonianPharmacoVariantsCol8), ]
UK10KPharmacoVariants<- UK10KPharmacoVariants[!duplicated(UK10KPharmacoVariants), ]
