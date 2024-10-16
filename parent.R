######Taste variable analysis
#1.
############Logistic regression for parents data#########################
setwd("~/Desktop/family result/")
install.packages("gdata")
library(gdata)
install.packages("Matrix")
library(Matrix)
install.packages("magic")
library(magic)
install.packages("geepack")
library(geepack)
#parent2 is the file where '-' values replaced by "NA" in order to be consistent.
data=read.xls("parent.xlsx",sheet=1,header=T)
dim(data)#65 rows, 124 col
attach(data)
###########
#TASTE ANALYSIS:
##########
#######################################
#Variable of interest:
#P_Sweet,P_Umami, P_Salt, P_Sour,P_Fat, PTC_Status,these are taste preference
#ST_Bitter,ST_Sweet,ST_Umami, ST_Salt,ST_Sour,ST_Fat, these are sensitivity threshold people start to feel the taste
#######################################
variable_col=(data[,7:18])
dim(variable_col)#65 rows, 12 cols
SNP_col=grep("rs",names(data))
SNPs=data[,SNP_col]
SNP_name=names(data[,SNP_col])
length(SNP_name)#106 SNPS
p_vector=c()
for (i in 1:dim(SNPs)[2]){
  Y=cbind(as.numeric(as.character(SNPs[,i])),2-as.numeric(as.character(SNPs[,i])))
  Formula=formula(paste("Y~1+Age+Sex+Ethnicity+BMI+",paste(names(variable_col),collapse = "+")))
  mod1=glm(Formula,data=data,family=binomial(link="logit"))
  p_value_ind=summary(mod1)$coefficients[-1,4]
  p_chi=1-pchisq(mod1$null.deviance-mod1$deviance, 16)
  p_value<-c(p_value_ind,p_chi)
  p_vector<-rbind(p_vector,p_value)
}
row.names(p_vector)<-SNP_name
colnames(p_vector)=c("Age","Sex","Ethnicity","BMI",names(variable_col),"Overall")
write.csv(p_vector,"/Users/nli/Desktop/family/family_pvalue_logistic.csv")
###############################################################################
############GEE for parents data########################
library(geepack)
taste_new=as.data.frame(cbind(apply(variable_col[,1:6],2,FUN=function(x){
  ifelse(x==2,1,ifelse(x==1,0,x))
  }),
  apply(variable_col[,7:12],2,FUN=function(x){
  x[x==0]=0.001
  return(log(x))
 })))
sex=data[,4]
sex=ifelse(sex==1,0,ifelse(sex==2,1,sex))
family=unique(data[,1])

lst=lapply(family,FUN=function(x){
  rows=which(data[,1]==x)
  mtx=matrix(0,length(rows),length(rows))
  diag(mtx)=1
  return(mtx)
})
blok=Reduce(adiag,lst)
cor.fixed=blok
wave=do.call(c,lapply(as.data.frame(table(data[,1]))[,2],function(x){return(1:x)}))
p_list=list()
sig_snp=list()
ct=0
ct2=0
sd_vec=c()
#################################################################
#FIRST,TEST SNP~AGE+SEX+BMI+TASTE VARIABLE USING GEE.           #
##################################################################
for (i in 1:dim(SNPs)[2]){
  for (j in 1:dim(taste_new)[2]){
     ct=ct+1
     Y=cbind(as.numeric(as.character(SNPs[,i])),2-as.numeric(as.character(SNPs[,i])))
     dat=as.data.frame(cbind(SNPs[,i],data$Family.ID,data$Age,sex,data$BMI,taste_new[,j]))
     names(dat)=c(names(SNPs)[i],"ID","Age","Sex","BMI",names(taste_new)[j])
     index=which(as.data.frame(apply(is.na(dat),1,sum))[,1]>0)
     sd_vec=rbind(sd_vec,c(sd(dat[,1],na.rm=T),sd(dat[,6],na.rm=T),names(SNPs)[i],names(taste_new)[j]))
     
     if (length(index)>0){
       zcor=fixed2Zcor(cor.fixed[-index,-index],id=factor(dat[-index,2]),waves = wave[-index])
       new_mtx=dat[-index,]
       res=Y[-index,]
     }
     else{
       zcor=fixed2Zcor(cor.fixed,id=factor(dat[,2]),waves = wave)
       new_mtx=dat
       res=Y
     }
     Formula=formula(paste("res~1+Age+Sex+BMI+",names(taste_new)[j]))
     mod1=geeglm(Formula,data=new_mtx[,-1],id=ID,family=binomial(link="logit"),corstr = "unstructured", zcor=zcor)
     mod=glm(Formula,data=new_mtx[,-1],family=binomial(link="logit"))
     p_value=as.data.frame(t(summary(mod1)$coefficients[,4]))
     names(p_value)=c("Intercept","Age","Sex","BMI",names(taste_new)[j])
     rownames(p_value)=names(SNPs)[i]
    
     if (!is.nan(p_value[,5]) & p_value[,5]<=0.05){
       ct2=ct2+1
       sig_snp[[ct2]]=p_value
     }
      p_list[[ct]]=p_value   
     }
}
write.table(sd_vec,"standard-deviation_parent.txt")
lapply(p_list, function(x) write.table(data.frame(x), 'p_list_parent.csv'  , row.names=T,col.names = TRUE,append=TRUE,sep=","))
lapply(sig_snp, function(x) write.table(data.frame(x), 'sig_p_list_parent.csv'  , row.names=T,col.names = TRUE,append=TRUE,sep=","))
##############################################################
#NEXT, TEST TASTE VARIABLE~AGE+SEX+BMI+SNP.                 #
###############################################################
for (i in 1:dim(SNPs)[2]){
  for (j in 1:dim(taste_new)[2]){
    ct=ct+1
    Y=taste_new[,j]
    dat=as.data.frame(cbind(Y,data$Family.ID,data$Age,sex,data$BMI,SNPs[,i]))
    names(dat)=c(names(taste_new)[j],"ID","Age","Sex","BMI",names(SNPs)[i])
    index=which(as.data.frame(apply(is.na(dat),1,sum))[,1]>0)
    sd_vec=rbind(sd_vec,c(sd(dat[,1],na.rm=T),sd(dat[,6],na.rm=T),names(SNPs)[i],names(taste_new)[j]))
    
    if (length(index)>0){
      zcor=fixed2Zcor(cor.fixed[-index,-index],id=factor(dat[-index,2]),waves = wave[-index])
      new_mtx=dat[-index,]
      res=Y[-index]
    }
    else{
      zcor=fixed2Zcor(cor.fixed,id=factor(dat[,2]),waves = wave)
      new_mtx=dat
      res=Y
    }
    Formula=formula(paste("res~1+Age+Sex+BMI+",names(SNPs)[i]))
    mod1=geeglm(Formula,data=new_mtx[,-1],id=ID,family=gaussian(link = "identity"),corstr = "unstructured", zcor=zcor)
    mod=glm(Formula,data=new_mtx[,-1],family=gaussian(link = "identity"))
    p_value=as.data.frame(t(summary(mod1)$coefficients[,4]))
    names(p_value)=c("Intercept","Age","Sex","BMI",names(SNPs)[i])
    rownames(p_value)=names(taste_new)[j]
    
    if (!is.nan(p_value[,5]) & p_value[,5]<=0.05){
      ct2=ct2+1
      sig_snp[[ct2]]=p_value
    }
    p_list[[ct]]=p_value   
  }
}
write.table(sd_vec,"standard-deviation_parent.txt")
lapply(p_list, function(x) write.table(data.frame(x), 'p_list_parent.csv'  , row.names=T,col.names = TRUE,append=TRUE,sep=","))
lapply(sig_snp, function(x) write.table(data.frame(x), 'sig_p_list_parent.csv'  , row.names=T,col.names = TRUE,append=TRUE,sep=","))
#family_p=read.csv("family_pvalue.csv",sep=",",header=T,row.names = T)
#family_p[,6:17]which(family_p[,6:17]<0.05)

#a=apply(family_p[,6:17],2,FUN=function(x){
#  rname=family_p[,1][which(x<=0.05)]
#  pv=x[x<=0.05]
#    return(paste(rname,pv))
# }
#)
#sum(sapply(a,length))#65
###################################################################################
#2.
#GQLS to test SNPS and tastes in Wishes List 
Wishes=read.xls("/Users/nli/Desktop/KinInbcoef/WishesList.xlsx",sheet=1)[,1:2]
dim(Wishes)#62 by 2
head(Wishes)

child=read.xls("Child.xlsx",sheet=1,header=T)
dim(child)#62 by 118
########################
#A.geno.dat
######"convert" function is used to change individual ID from "P"/"S" to all numeric values
convert=function(dat,A,B,a,b,c){
  id_vec=c()
  for (i in 1:dim(dat)[1]){
    if (substr(as.character(dat[i,2]),1,1)==A){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),a)
    }
    else if (substr(as.character(dat[i,2]),1,1)==B){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),b)
    }
    else{
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),c)
    }
  }
  id_vec=as.numeric(id_vec)
  return(id_vec)
}
#change family ID column in parent and child and wishes
id_parent=convert(data,"P","S",1,2,0)
id_child=convert(child,"A","B",3,4,5)

#data[,2]=id_parent
#child[,2]=id_child
#Wishes[,2]=id_child
###

#change family ID column in parent and child, and wishes 
familyID=function(data){
  count=1
  familyID_vec=c(count)
  for (i in 2:dim(data)[1]){
    if (data[i,1] != data[i-1,1]){
      count=count+1
      familyID_vec=c(familyID_vec,count)
    }
    else{
      familyID_vec=c(familyID_vec,count)
    }
  }
  return(familyID_vec)
}
family_child=familyID(child)
family_parent=c()
for (i in data[,1]){
  mt=match(i,child[,1])
  family_parent=c(family_parent,ifelse(is.na(family_child[mt]),i,family_child[mt]))
}
data[,1]=family_parent
child[,1]=family_child
Wishes[,1]=family_child

individualID=function(id_vec){
  ind_ID=c()
  for (i in 1:length(id_vec)){
    ID_num=as.numeric(strsplit(as.character(id_vec[i]),"")[[1]])
    num=ID_num[length(ID_num)]
    ind_ID=c(ind_ID,num)
  }
  return(ind_ID)
}
ind_parent=individualID(id_parent)
ind_child=individualID(id_child)

geno.dat=matrix(data=NA,nrow=62,ncol=112)
dim(geno.dat)#62 by 112
for (i in 1:dim(Wishes)[1]){
  index=which(Wishes[i,1]==data[,1])#refers to parent dataset
  geno.dat[i,1:2]=as.vector(Wishes[i,1:2],mode="numeric")#family ID, indi ID
  
  male_index=id_parent[index[data[index,4]==1]]
  female_index=id_parent[index[data[index,4]==2]]
  geno.dat[i,3]=ifelse(!length(male_index),0,male_index)#father ID
  geno.dat[i,4]=ifelse(!length(female_index),0,female_index)#mother ID
  
  
  geno.dat[i,5]=child[i,4]#sex
  geno.dat[i,6]=ifelse(child[i,7]==1,0,ifelse(child[i,7]==2,1,-99999))#trait is "p_sweet"
  geno.dat[i,7:112]=as.vector(child[i,13:118],mode="numeric")#106 SNPs
  geno.dat[i,7:112][is.na(geno.dat[i,7:112])]=-5
}
head(geno.dat)
colnames(geno.dat)=c("family ID","individual ID",
                     "father's ID","mother's ID","sex",
                     "trait(age)",names(child)[13:118])
write.csv(geno.dat,"geno_dat.csv",row.names = FALSE)
#############################################################################
#B.kinship coefficient:relationship between people including themselves
#input:"WishesList.xlsx"
#output: kinship coef matrix
inbred.mtx=NULL
for (i in 1:(dim(Wishes)[1])){
  inbred.mtx=rbind(inbred.mtx,c(Wishes[i,],Wishes[i,2],0))
  if (i<dim(Wishes)[1]){
    for (j in (i+1):dim(Wishes)[1]){
      if (Wishes[i,1]==Wishes[j,1]){
        inbred.mtx=rbind(inbred.mtx,c(Wishes[i,],Wishes[j,2],0.25))
      }
    }
  }
}

colnames(inbred.mtx)=c("famiy.ID","individual.ID1","individual.ID2","inbred.coef")
kin.dat=inbred.mtx
write.csv(kin.dat,"kin_dat.csv",row.names=FALSE)
##############################################################################
GQLS=function(geno.dat, kin.dat){
  F=max(geno.dat[,1])
  nf=numeric(F)
  for ( i in 1:F){
    nf[i]=sum(geno.dat[,1]==i)
  }
  kinlist=vector("list", F)
  ind=1
  for (f in 1:F){
    kinlist[[f]]=matrix(0, nf[f],nf[f])
    for ( i in 1:(nf[f])){
      for ( j in i:(nf[f])){
        kinlist[[f]][i,j]=kin.dat[ind,4]
        ind=ind+1
      }
    }
    kinlist[[f]]=2*(kinlist[[f]]+t(kinlist[[f]]))+diag(rep(1, nf[f]))
  }
  
  invkinlist=vector("list", F)
  for (f in 1:F){
    invkinlist[[f]]=solve(kinlist[[f]])
  }
  KK=ncol(geno.dat)
  muhat=rep(0,(KK-6))
  Wg.stat=rep(0,(KK-6))
  for ( i in 1:(KK-6)){
    ydatlist=vector("list", F)
    xdatlist=vector("list", F)
    klist=vector("list", F)
    invklist=vector("list", F)
    nalist=vector("list", F)
    numfam=nf
    Q1=rep(0,F)
    Q2<-rep(0,F)
    Q3<-rep(0,F)
    Q4<-rep(0,F)
    for (f in 1:F){
      xdatlist[[f]]=geno.dat[geno.dat[,1]==f, 6]
      ydatlist[[f]]=geno.dat[geno.dat[,1]==f, 6+i]
      temp=NULL
      for ( j in 1:length(ydatlist[[f]])){
        if (ydatlist[[f]][j]==-5|xdatlist[[f]][j]==-99999){
          temp=c(temp, j)
        }
      }
      if (length(temp)!=0){
        nalist[[f]]=temp
        ydatlist[[f]]=ydatlist[[f]][-nalist[[f]]]
        xdatlist[[f]]=xdatlist[[f]][-nalist[[f]]]
        klist[[f]]=kinlist[[f]][-nalist[[f]],-nalist[[f]]]
        numfam[f]<-length(ydatlist[[f]])
        if (numfam[f]!=0){
          invklist[[f]]=solve(klist[[f]])
        }
      }
      else invklist[[f]]=invkinlist[[f]]
      if (numfam[f]!=0){
        ydatlist[[f]]=ydatlist[[f]]/2
        onevec<-rep(1,numfam[f])
        Q1[f]<-t(xdatlist[[f]])%*%invklist[[f]]%*%xdatlist[[f]]
        Q2[f]<-t(xdatlist[[f]])%*%invklist[[f]]%*%onevec
        Q3[f]<-t(onevec)%*%invklist[[f]]%*%onevec
        Q4[f]<-t(onevec)%*%invklist[[f]]%*%ydatlist[[f]]
      }
    }
    muhat[i]<-sum(Q4)/sum(Q3)
    crosspro=rep(0,F)
    for (f in 1:F){
      if(numfam[f]!=0){
        onevec=rep(1,numfam[f])
        crosspro[f]<-t(ydatlist[[f]]-muhat[i]*onevec)%*%invklist[[f]]%*%xdatlist[[f]]
      }
    }
    invFmat<-2/(muhat[i]-muhat[i]^2)
    Wg.stat[i]<-(sum(crosspro)*sum(crosspro)*invFmat)/(sum(Q1)-sum(Q2)^2/sum(Q3))
  }
  pval=1-pchisq(Wg.stat,1)
  marker=seq(1:(KK-6))
  outfile=cbind(marker, muhat,Wg.stat, pval)
  return(outfile)
}
geno.dat=read.table("geno_dat.csv",header=TRUE,sep=",")
kin.dat=read.table("kin_dat.csv",header=TRUE,sep=",")
class(kin.dat)
head(kin.dat)
##################################
##GQLS function to test SNPs and five taste variables:
for (i in 7:12){
  #substitute 6th column with trait variable i
  geno.dat[,6]=ifelse(child[,i]==1,0,ifelse(child[,i]==2,1,-99999))#trait is i
  result=GQLS(geno.dat,kin.dat)
  write.csv(result,paste0(names(child)[i],".csv"),row.names=FALSE)
}
###################################
##GQLS function to test SNPs and four traits:
vec_sig=c()
for (i in 3:6){
  geno.dat[,6]=ifelse(child[,i]==1,0,ifelse(child[,i]==2,1,-99999))#trait is i
  geno.dat[,6]=ifelse(is.na(child[,i]),-99999,child[,i])
  result=GQLS(geno.dat,kin.dat)
  SNP_name_sig=names(child)[13:118][which(result[,4]<=0.05)]
  SNP_sig=result[result[,4]<=0.05,4]
  vec=cbind(names(child)[i],SNP_name_sig,round(SNP_sig,4))
  write.csv(result,paste0(names(child)[i],".csv"),row.names=FALSE)
  vec_sig=rbind(vec_sig,vec)
}
write.csv(vec_sig,"significant SNPs_four_traits-GQLS.csv",row.names=FALSE)
###########################################################################
####Diet Analysis
###########################################################################
library(gdata)
setwd("/Users/nli/Desktop/family/diet_variables/")
SNP=read.xls("SNPs.xlsx",sheet=1,header=TRUE)#gdata package
dim(SNP)#60 by 11
phy=read.xls("PHY.xlsx",sheet=1,header=TRUE)
names(phy)
dim(phy)#60 by 21
#############################################################################
#data preprocessing 
#1.normalization
par(mar=c(1,1,1,1))
layout(matrix(1:12, nrow=4, byrow=T))
for (i in 1:12){
  hist(log(phy[,10:21][,i]),main= names(phy)[10:21][i])
}
#all positive ->log-transform
phy_new=cbind(phy[,1:9],apply(phy[,10:21],2,FUN=function(x){
  x[x==0]=0.001
  return(log(x))}))

#2.block matrix
setwd("~/Desktop/family/diet_variables")
geno.fat=read.xls("geno fat.xlsx",sheet=1,header=TRUE)[,-c(4,6)]
inbred=NULL

for (i in 1:(dim(geno.fat)[1])){
  inbred=rbind(inbred,c(geno.fat[i,1:2],geno.fat[i,2],0))
  if (i<dim(geno.fat)[1]){
    for (j in (i+1):dim(geno.fat)[1]){
      if (geno.fat[i,1]==geno.fat[j,1]){
        inbred=rbind(inbred,c(geno.fat[i,1:2],geno.fat[j,2],0.25))
        
      }
    }
  }
}
colnames(inbred)=c("famiy.ID","individual.ID1","individual.ID2","inbred.coef")
inbred=matrix(unlist(inbred),nrow=length(inbred[,1]),byrow=F)
#symmatrix matrix in list:lst
n=length(geno.fat[,1])#60 ind
family=unique(geno.fat[,1])
lst=lapply(family,FUN=function(x){
  rows=which(geno.fat[,1]==x)
  mtx=matrix(NA,length(rows),length(rows))
  mtx[lower.tri(mtx,diag=T)]=inbred[which(inbred[,1]==x),4]
  mtx=forceSymmetric(mtx,"L")#it contains character, so use "as.matrix" below to convert
  return(as.matrix(diag(1,nrow=nrow(mtx),ncol=ncol(mtx))+2*mtx))
})
#block
blok=Reduce(adiag,lst)
cor.fixed=blok
#gee
library(geepack)
phenotype=phy_new[,10:21]
covariate=phy_new[,1:4]
taste=phy_new[,5:9]
taste=apply(taste,2,function(x){
  ifelse(x==1,0,ifelse(x==2,1,x))
})

attach(SNP)
attach(phy_new)
sex=covariate[,1]
sex=ifelse(sex==1,0,ifelse(sex==2,1,sex))
wave=do.call(c,lapply(as.data.frame(table(geno.fat[,1]))[,2],function(x){return(1:x)}))
###################################
#TASK1:Diet with four traits
p_vec=NULL
#for (k in 1:5){#k is the set of variables for one taste 
   for (j in 1:dim(phenotype)[2]){ #j is diet
     Y=phenotype[,j]
     name_phy=colnames(phenotype)[j]
     id=as.vector(geno.fat[,1])
     dat=as.data.frame(cbind(Y,id,sex,Age,BMI))
     index=which(as.data.frame(apply(is.na(dat),1,sum))[,1]>0)
     zcor=fixed2Zcor(cor.fixed[-index,-index],id=factor(geno.fat[-index,1]),waves = wave[-index])
     new_mtx=dat[-index,]
     names(new_mtx)=c("Y","ID","Sex","Age","BMI")
     Formula=formula(paste(c("Y~1",paste(colnames(new_mtx)[-(1:2)],collapse="+")),collapse="+"))
     mod=geeglm(Formula,family=gaussian,id=ID,data=new_mtx,corstr="unstructured", zcor=zcor)
     p_value=as.data.frame(t(summary(mod)$coefficients[,4]))
     row.names(p_value)=name_phy
     names(p_value)=rownames(summary(mod)$coefficients)
     print(c(rownames(summary(mod)$coefficients),name_phy))
     p_vec<-rbind(p_vec,p_value)
      }
write.csv(p_vec,"P-values between diet and covariates(GEE).csv")
##############################
#2.TASK2: multiple traits with one taste, one diet
snp_list=list(1,2:4,5:8,9:10,11)
phy_list=list(1:3,4,c(1,5:6),c(1,5:10),c(1,11:12))
p_list=list()
ct=0
for (k in 1:5){#k is the set of variables for one taste 
 for (j in phy_list[[k]]){ #j is diet
   for (i in snp_list[[k]]){
     ct=ct+1
     pnt=phenotype[,j]
     Y=SNP[,i]
     name_phy=colnames(phenotype)[j]
     id=as.vector(geno.fat[,1])
     dat=as.data.frame(cbind(Y,id,taste[,k],pnt))
     index=which(as.data.frame(apply(is.na(dat),1,sum))[,1]>0)
     zcor=fixed2Zcor(cor.fixed[-index,-index],id=factor(geno.fat[-index,1]),waves = wave[-index])
     new_mtx=dat[-index,]
     names(new_mtx)=c("Y","ID","taste","diet")
     res=cbind(new_mtx[,1],2-new_mtx[,1])
     Formula=formula(paste(c("res~1",paste(colnames(new_mtx)[-(1:2)],collapse="+")),collapse="+"))
     mod=geeglm(Formula,family=binomial(link="logit"),id=ID,data=new_mtx[,-1],corstr="unstructured", zcor=zcor)
     p_value=as.data.frame(t(summary(mod)$coefficients[,4]))
     row.names(p_value)=names(SNP)[i]
     names(p_value)=paste(c("Intercept",colnames(taste)[k],name_phy))
     p_list[[ct]]<-p_value
   }
 }
}
lapply(p_list, function(x) write.table( data.frame(x), 'p_list_diet.csv'  , row.names=T,col.names = TRUE,append=TRUE,sep=","))
#############################################################################################
#Combined dataset: test SNPs with taste variable:P_Sweet, P_
data=data[,-c(5,13:18)]
child=child[,-5]
#standardize BMI in data in order to be consistent with child BMI
data[,5]=(data[,5]-mean (data[,5],na.rm=T))/sd(data[,5],na.rm=T)
names(child)[5]="BMI"
names(child) %in% names(data)
length(names(child))#117
length(names(data))#117
newdata=rbind(data,child)
dim(newdata)#127 117
combine_data=newdata[order(newdata[,1]),]
variable_col=combine_data[,6:11]
taste_new=as.data.frame(cbind(apply(variable_col,2,FUN=function(x){
  ifelse(x==2,1,ifelse(x==1,0,x))
})))
sex=combine_data[,4]
sex=ifelse(sex==1,0,ifelse(sex==2,1,sex))
convert2=function(dat,A,B,C,D,a,b,c,d,e){
  id_vec=c()
  for (i in 1:dim(dat)[1]){
    if (substr(as.character(dat[i,2]),1,1)==A){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),a)
    }
    else if (substr(as.character(dat[i,2]),1,1)==B){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),b)
    }
    else if (substr(as.character(dat[i,2]),1,1)==C){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),c)
    }
    else if (substr(as.character(dat[i,2]),1,1)==D){
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),d)
    }
    else{
      id_vec[i]=paste0(substr(as.character(dat[i,2]),2,nchar(as.character(dat[i,2]))),e)
    } 
  }

  id_vec=as.numeric(id_vec)
  return(id_vec)
}
#change family ID column in parent and child 
id=convert2(combine_data,"P","S","A","B",1,2,3,4,0)
combine_data[,2]=id
inbred=NULL
for (i in 1:(dim(combine_data)[1])){
  inbred=rbind(inbred,c(combine_data[i,1:2],combine_data[i,2],0))
  if (i<dim(combine_data)[1]){
    for (j in (i+1):dim(combine_data)[1]){
      if (combine_data[i,1]==combine_data[j,1]){
        if (all(is.element(c(combine_data[i,2] %%10,combine_data[j,2] %% 10),c(1,2)))){
           inbred=rbind(inbred,c(combine_data[i,1:2],combine_data[j,2],0))
        }
        else if (all(is.element(c(combine_data[i,2] %%10,combine_data[j,2] %% 10),c(3,4)))){
          inbred=rbind(inbred,c(combine_data[i,1:2],combine_data[j,2],0.25))
        }
        else {
          inbred=rbind(inbred,c(combine_data[i,1:2],combine_data[j,2],0.5))
        }
      }
    }
  }
}
colnames(inbred)=c("famiy.ID","individual.ID1","individual.ID2","inbred.coef")
inbred=matrix(unlist(inbred),nrow=length(inbred[,1]),byrow=F)
family=unique(combine_data[,1])
lst=lapply(family,FUN=function(x){
  rows=which(combine_data[,1]==x)
  mtx=matrix(NA,length(rows),length(rows))
  mtx[lower.tri(mtx,diag=T)]=inbred[which(inbred[,1]==x),4]
  mtx=forceSymmetric(mtx,"L")#it contains character, so use "as.matrix" below to convert
  return(as.matrix(diag(1,nrow=nrow(mtx),ncol=ncol(mtx))+2*mtx))
})
blok=Reduce(adiag,lst)
cor.fixed=blok
wave=do.call(c,lapply(as.data.frame(table(combine_data[,1]))[,2],function(x){return(1:x)}))
SNP_col=grep("rs",names(combine_data))
SNPs=combine_data[,SNP_col]
SNP_name=names(combine_data[,SNP_col])
length(SNP_name)#106 SNPS
p_list=list()
sig_snp=list()
ct=0
ct2=0

for (i in 1:dim(SNPs)[2]){
  for (j in 1:dim(taste_new)[2]){
    ct=ct+1
    Y=cbind(as.numeric(as.character(SNPs[,i])),2-as.numeric(as.character(SNPs[,i])))
    dat=as.data.frame(cbind(SNPs[,i],combine_data$Family.ID,combine_data$Age,sex,taste_new[,j]))
    names(dat)=c(names(SNPs)[i],"ID","Age","Sex",names(taste_new)[j])
    index=which(as.data.frame(apply(is.na(dat),1,sum))[,1]>0)
    if (length(index)>0){
      zcor=fixed2Zcor(cor.fixed[-index,-index],id=factor(dat[-index,2]),waves = wave[-index])
      new_mtx=dat[-index,]
      res=Y[-index,]
    }
    else{
      zcor=fixed2Zcor(cor.fixed,id=factor(dat[,2]),waves = wave)
      new_mtx=dat
      res=Y
    }
    Formula=formula(paste("res~1+Age+Sex+",names(taste_new)[j]))
    mod1=geeglm(Formula,data=new_mtx[,-1],id=ID,family=binomial(link="logit"),corstr="unstructured", zcor=zcor)
    p_value=as.data.frame(t(summary(mod1)$coefficients[,4]))
    names(p_value)=c("Intercept","Age","Sex",names(taste_new)[j])
    rownames(p_value)=names(SNPs)[i]
    if (is.nan(p_value[,4])){
      print ((c(i,j,"NAN")))
    }
    if (!is.nan(p_value[,4]) & p_value[,4]==0){
      print (c(i,j,"p_value=0"))
      print(c(SNPs[i],names(taste_new)[j]))
    }
    if (!is.nan(p_value[,4]) & p_value[,4]<=0.05){
      ct2=ct2+1
      sig_snp[[ct2]]=p_value
    }
    p_list[[ct]]=p_value   
  }
}
lapply(p_list, function(x) write.table(data.frame(x), 'p_list_combine.csv', row.names=T,col.names = TRUE,append=TRUE,sep=","))
lapply(sig_snp, function(x) write.table(data.frame(x), 'sig_p_list_combine.csv', row.names=T,col.names = TRUE,append=TRUE,sep=","))
#######################################################################
#Figures
#rs713598: GG=2, GC=1, CC=0
#rs236514: GG=2, GA=1, AA=0
#1.The rs713598 and rs236514 SNPs remained significantly associated with taste outcomes in parents.  
#- The C allele of the rs173598 SNP in the TAS2R38 bitter taste receptor gene was significantly associated with PTC sensitivity. (p=0.003) (single trait)
#- The A allele of the rs236514 SNP in the KCNJ2 sour taste receptor gene was significantly associated with sour preference. (p=0.002) (single trait)
#par(mfrow=c(2,3))
#1.
par(mfrow=c(1,1))
#par(mfrow=c(3,2))
PTC=as.data.frame(cbind(data[,12],data[,names(data)=="rs713598"]))#extract that column
PTC[,2]=ifelse(PTC[,2]==0,2,ifelse(PTC[,2]==2,0,1))
colnames(PTC)=c("PTC_status","genotype")
name1=c("GG","GC","CC")
count1=table(PTC[,1],PTC[,2])
mybar1=barplot(count1,main="PTC Sensitivity in rs713598",names.arg = name1,ylim=c(0,max(count1)+3),
               xlim=c(0,13),xlab="Genotype",ylab="Frequency of PTC Taster Status",col=c("black","white"),
       legend=c("PTC Non Taster","PTC Taster"),beside=TRUE,axis.lty=1)#axis.lty=1 x asix suppresed
text(mybar1,as.vector(count)+2,as.vector(count),cex=0.9)
#2.
Sour=as.data.frame(cbind(data[,10],data[,names(data)=="rs236514"]))#extract that column
colnames(Sour)=c("P_Sour","genotype")
name2=c("AA","GA","GG")
count2=table(Sour[,1],Sour[,2])
mybar2=barplot(count2,main="Sour Preference in rs236514",names.arg = name2,ylim=c(0,max(count2)+9),
              xlim=c(0,12),xlab="Genotype",ylab="Frequency of Sour Preference",col=c("black","white"),
              beside=TRUE,axis.lty=1)#axis.lty=1 x asix suppresed
text(mybar2,as.vector(count2)+2,as.vector(count2),cex=0.9)
legend("topright",legend=c("Sour Not Preferred","Sour Preferred"),cex=0.85,fill=c("black","white"))
#3.
#Two tSNPs remained significantly associated with a taste outcome in children after applying a Bonferroni adjustment. 
#- The C allele of the rs4790522 tSNP in the TRPV1 salt taste receptor gene was associated with a significantly higher salt preference compared to the A allele in children. (p=0.001) (single trait)
#- The T allele of the rs173135 tSNP in the KCNJ2 sour taste receptor gene was associated with a significantly higher sour preference compared to the C allele in children. (p<0.001) (single trait)
#The rs9701796 SNP in the TAS1R2 sweet taste receptor gene was associated with both sweet taste preference (p=0.022) and percent energy from added sugar in the diet (p=0.047). (multiple trait)
#rs4790522: CC=2, CA=1, AA=0
Salt=as.data.frame(cbind(child[,9],child[,names(child)=="rs4790522"]))#extract that column
Salt[,2]=ifelse(Salt[,2]==0,2,ifelse(Salt[,2]==2,0,1))
#change minor allele to be major allele

colnames(Salt)=c("P_Salt","genotype")
name3=c("CC","CA","AA")
count3=table(Salt[,1],Salt[,2])
mybar3=barplot(count3,main="Salt Preference in rs4790522",names.arg = name3,ylim=c(0,max(count3)+3),
               xlim=c(0,12),xlab="Genotype",ylab="Frequency of Salt Preference",col=c("black","white"),
               beside=TRUE,axis.lty=1)#axis.lty=1 x asix suppresed
text(mybar3,as.vector(count3)+2,as.vector(count3),cex=0.9)
legend("topright",legend=c("Salt Not Preferred","Salt Preferred"),cex=0.85,fill=c("black","white"))

#4.#rs173135: CC=2, CT=1, TT=0,62 individuals from child
Sour_c=as.data.frame(cbind(child[,10],child[,names(child)=="rs173135"]))#extract that column
colnames(Sour_c)=c("P_Sour","genotype")
name4=c("TT","CT","CC")
count4=table(Sour_c[,1],Sour_c[,2])
mybar4=barplot(count4,main="Sour Preference in rs173135",names.arg = name4,ylim=c(0,max(count4)+12),
               xlim=c(0,13),xlab="Genotype",ylab="Frequency of Sour Preference",col=c("black","white"),
               beside=TRUE,axis.lty=1)#axis.lty=1 x asix suppresed
text(mybar4,as.vector(count4)+2,as.vector(count4),cex=0.9)
legend("topright",legend=c("Sour Not Preferred","Sour Preferred"),cex=0.75,fill=c("black","white"))


#5.#rs9701796: CC=2, CG=1, GG=0,60 individuals from PHY
Sweet=as.data.frame(cbind(phy_new[,8],phy_new[,15],SNP[,names(SNP)=="rs9701796"]))#extract that column
colnames(Sweet)=c("P_Sweet","Sugar_add_percent","genotype")
name5=c("GG","CG","CC")
count5=table(Sweet[,1],Sweet[,3])
mybar5=barplot(count5,main="Sweet Preference in rs9701796",names.arg = name5,ylim=c(0,max(count5)+12),
               xlim=c(0,14),xlab="Genotype",ylab="Frequency of Sweet Preference",col=c("black","white"),
               beside=TRUE,axis.lty=1)#axis.lty=1 x asix suppresed
text(mybar5,as.vector(count5)+2,as.vector(count5),cex=0.9)
legend("topright",legend=c("Sweet Not Preferred","Sweet Preferred"),cex=0.75,fill=c("black","white"))

myplot5=plot(Sweet[,3],Sweet[,2],main="Sugar Added Percentage in rs9701796",
             ylim=c(-7,2),xlab="Genotype",ylab="Log Sugar Added Percentage",xaxt="n")
axis(1, at=c(0,1,2), labels=name5) 

