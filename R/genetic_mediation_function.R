SNP_mediation <- function(b.h,se.bh,a.h,se_ah){
  mediation.h=b.h*a.h
  se_mediation.h=sqrt(se.bh^2*a.h^2+se_ah^2*b.h^2+se.bh^2*se_ah^2)
  mediation.p=pchisq((mediation.h/se_mediation.h)^2,1,lower.tail = F)
  return(as.data.frame(cbind(mediation.h,se_mediation.h,mediation.p)))
}



GWAS_mediation <- function(GWAS_M,a.h,se_ah,GWAS_Y=NULL){
  b.h=GWAS_M[,which(names(GWAS_M) %in% c("beta","BETA","b"))]
  se.bh=GWAS_M[,which(names(GWAS_M) %in% c("se","SE"))]
  mediation_table=SNP_mediation(b.h,se.bh,a.h,se_ah)
  if(!is.null(GWAS_Y)){
    b.h.Y=GWAS_Y[,which(names(GWAS_Y) %in% c("beta","BETA","b"))]
    mediation_table$mediation_proportion=mediation_table$mediation.h/b.h.Y
  }
  return(mediation_table)
}



read_BOLT <- function(BOLT_result){
  a=read.table(BOLT_result,sep="\t")
  V1 = a[apply(a,1,FUN=startsWith,prefix="Phenotype 1 variance sigma2: "),]
  V1=substr(V1,30,nchar(V1))
  V1=as.numeric(strsplit(V1,"[\\(\\)\\ ]")[[1]][1])
  V2 = a[apply(a,1,FUN=startsWith,prefix="Phenotype 2 variance sigma2: "),]
  V2=substr(V2,30,nchar(V2))
  V2=as.numeric(strsplit(V2,"[\\(\\)\\ ]")[[1]][1])
  C = a[apply(a,1,FUN=startsWith,prefix="  h2e (1,1): "),]
  C=substr(C,14,nchar(C))
  res_M=as.numeric(strsplit(C,"[\\(\\)\\ ]")[[1]][1])*V1
  C = a[apply(a,1,FUN=startsWith,prefix="  resid corr (1,2): "),]
  C=substr(C,21,nchar(C))
  COR_res_Mres_Y_est=as.numeric(as.numeric(strsplit(C,"[\\(\\)\\ ]")[[1]][1]))
  C = a[apply(a,1,FUN=startsWith,prefix="  h2e (2,2): "),]
  C=substr(C,14,nchar(C))
  res_Y=as.numeric(strsplit(C,"[\\(\\)\\ ]")[[1]][1])*V2
  COV_res_Mres_Y_est = COR_res_Mres_Y_est*sqrt(res_M*res_Y)
  alphahat_adj=COV_res_Mres_Y_est/res_M
  D= a[apply(a,1,FUN=startsWith,prefix="  entry (1,2): "),][1]
  D1 = as.numeric(strsplit(D,"[\\(\\)\\ ]")[[1]][7])
  D2 = as.numeric(strsplit(D,"[\\(\\)\\ ]")[[1]][9])
  SE = alphahat_adj/D1*D2
  return(list(alpha.h=alphahat_adj,SE_ah=SE))
}