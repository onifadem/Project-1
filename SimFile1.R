##Program for simulating files for running ROADTRIPS statistics
##Genotype datafile must be formatted such that each row corresponds to a SNP and 
#each column corresponds to an indivdual. Same format as the LTMLM
##Copying same simulation setup

numpop=2
N=2500 ##number in each population 
nSNP=5000  #number of  snps
Fst=0.01
omega=c(0.5,0.5) ##MAF
propnExtreme=0.1
nsim=100
Fst.obs=vector(length=nSNP)
pdiffs=vector(length=nSNP)
genomat=matrix(nrow=N*numpop,ncol=nSNP)
##Simulate snps for each population
for (i in 1:nSNP){
  p=runif(1,0.1,0.9) ##generating allele frequency p from uniform distribution
  alpha=p*(1-Fst)/Fst
  beta=(1-p)*(1-Fst)/Fst
  ps=rbeta(numpop,shape1=alpha,shape2=beta)
  
  for (j in 1:numpop){
    ind1=(j-1)*(N*numpop)*omega[j]+1 ##index values for individuals in population 1
    ind2=j*(N*numpop)*omega[j]  #individuals in population 2
    freqs=c(ps[j]^2,2*ps[j]*(1-ps[j]),(1-ps[j])^2)
    genomat[ind1:ind2,i]=sample(c(0,1,2),size=N*omega[j],replace=TRUE,prob=freqs)
  }
  X= genomat
}

##To simulate the putative causal locus
#set minor alelle frequencies in the two subpopulations to be p1=0.25 and p2=0.85
p1=0.25; p2 = 0.85
geno_SNP<- sample(c(0:2), (N*numpop), prob=c(p1^2, 2*p1*p2, p2^2), replace=T)
##simulate phenotype values from normal distribution.
##simulate phenotypes dependent on population since we are testing for type 1 error. 
pheno_indep <-c()
pheno1<- c(pheno_indep, rnorm((ind1-1), mean= 0.07, sd=1))
pheno2<- c(pheno_indep, rnorm((ind1:ind2), mean= -0.07, sd=1))
pheno_indep<- c(pheno1,pheno2)
##individual ID's and 
IND<- 1:(N*numpop)
combined_indep <- cbind(IND, pheno_indep, geno_SNP)
sorted_combined <- combined_indep[order(combined_indep[,2]),] ##sort to subset for EPS data
##Obtain the EPSdata 
K = propnExtreme #proportion in the extremes
Nums = nrow(sorted_combined)
keep <- c(1:(K*Nums), (Nums-(K*Nums)+1):Nums)
epsdat<- c(rep("0",K*Nums),rep("1",K*Nums))
EPS_pheno <- as.data.frame(cbind(epsdat,sorted_combined[keep,]))
dim(EPS_pheno)<- c(length(keep),4)
colnames(EPS_pheno)<- c("Bin_pheno", "Labels", "Phenotype", "Genotype")
write.table(Data_EPS, file = "EPS_pheno", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
##try running a logistic regression on this dataset
#model1<- glm(EPS_pheno[,1]~EPS_pheno[,4], family=binomial(link="logit"))  ##association between eps categories and putative locus








# write.table(Data_EPS, file = "EPS.pheno", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
# 
# #my_eps<- EPS_pheno[,-c(1:3)]
# #Data_EPS<- t(my_eps)
# #write.table(Data_EPS, file = "EPS.geno", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
# 
# EPS_geno<- EPS_pheno[,-c(1:3)]  ##genotype matrix of all SNPS
# ##rename the SNPS using SNP1 etc
# n= nrow(EPS_geno) ; prefix<- "SNP"; suffix<- seq(from=1, to=n)
# #IDs = paste(prefix,suffix,sep="")
# #names(EPS_pheno[,-c(1:3)]) <- c(paste(prefix,suffix,sep=""))
# names(EPS_geno) <- c(paste(prefix,suffix,sep=""))
# EPS_2<- t(EPS_geno)
# #randomly simulate reference and effect alleles
# Ref<- c(rep("A", nrow(EPS_2)))
# Eff<- sample(c("C","G","T"), nrow(EPS_2), replace=TRUE, prob= NULL)
# EP<- cbind(Ref,Eff, EPS_2)
# write.table(EP, file = "geno.txt", col.names = FALSE, row.names = FALSE, sep= "", quote=FALSE)
