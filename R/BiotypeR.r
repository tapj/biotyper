###################################################################
#                                                                 #
#           BiotypeR Methods                                      #
#                                                                 #
###################################################################


biotyper=function(x, ...) UseMethod("biotyper")

biotyper.data.frame=function(obs, k=3 , distance.jsd=NULL, cluster=NULL, manalysis=FALSE,  plot=FALSE, nf=3, no.unassigned=TRUE,...) {

require(ade4)
require(fpc)

if(k==1) stop("k cannot be 1")

if( is.null(distance.jsd)){
as.matrix(obs)-> obs.m
rownames(obs)-> rownames(obs.m)
dist.JSD(obs.m)-> distance.jsd
}

if (is.null(cluster) ){ cluster=cluster.biotypes(distance.jsd,k) }

obs.silhouette=mean(silhouette(cluster, distance.jsd)[,3])

obs.pca=NULL
obs.bet=NULL
obs.dpcoa=NULL

if(manalysis) {

if(no.unassigned){
dudi.pca(data.frame(t(obs[-1,])), scannf=F, nf=nf)-> obs.pca
} else { dudi.pca(data.frame(t(obs)), scannf=F, nf=nf)-> obs.pca }

bca(obs.pca, fac=as.factor(cluster), scannf=F, nf=nf)-> obs.bet
dpcoa(noise.removal(obs, percent=0.000001), distance.jsd, scannf=F, nf=nf) -> obs.dpcoa
}


biotyper.obj=list(biotypes=cluster, obs.silhouette=obs.silhouette, PCA=obs.pca, BET=obs.bet, dPCoA=obs.dpcoa)
class(biotyper.obj)="biotyper"

if (plot) { plot.biotyper(biotyper.obj) }

return(biotyper.obj)
}






plot.biotyper=function(biotyper.obj, xax=1, yax=2, potatoes=FALSE, ...) {

if (is.null(biotyper.obj$PCA)) stop ("there is no multivariate analysis inside this object\n=>restart 'biotyper' with the option 'manalysis=TRUE'")

par(mfrow=c(2,2))

#barplot
barplot(table(biotyper.obj$biotypes), xlab="biotypes", ylab="Nb samples", col=as.numeric(levels(as.factor(as.numeric(biotyper.obj$biotypes))))+1) #+1 to avoid black
box()


#PCA
plot(biotyper.obj$PCA$li[,xax], biotyper.obj$PCA$li[,yax], main="PCA", pch=16, col=as.numeric(biotyper.obj$biotypes)+1, xlab=paste("PC",xax), ylab=paste("PC",yax))
box()

if(!potatoes) {
#DPCOA
plot(biotyper.obj$dPCoA$dls[1:2], type="n", xlab=paste("PC",xax), ylab=paste("PC",yax), main="dPCoA")
s.class(biotyper.obj$dPCoA$dls, fac=factor(biotyper.obj$biotypes), add.plot=TRUE,clab=1, xax=xax, yax=yax)
points(biotyper.obj$dPCoA$dls[,xax], biotyper.obj$dPCoA$dls[,yax], col=as.numeric(biotyper.obj$biotypes)+1, cex=1, pch=16)

#BET
plot(biotyper.obj$BET$ls[,xax], biotyper.obj$BET$ls[,yax], type="n", xlab=paste("PC",xax), ylab=paste("PC",yax), main="between class")
s.class(biotyper.obj$BET$ls, fac=as.factor(biotyper.obj$biotypes), grid=F, xax=xax, yax=yax, add.plot=TRUE)
points(biotyper.obj$BET$ls[,xax], biotyper.obj$BET$ls[,yax], col=as.numeric(biotyper.obj$biotypes)+1, cex=1, pch=16)
} else {
#DPCOA
plot(biotyper.obj$dPCoA$l1[1:2], type="n", xlab=paste("PC",xax), ylab=paste("PC",yax), main="dPCoA")
s.potatoe(biotyper.obj$dPCoA$l1, fac=factor(biotyper.obj$biotypes), xax=xax, yax=yax)
s.class(biotyper.obj$dPCoA$l1, fac=factor(biotyper.obj$biotypes), add.plot=TRUE,clab=1, xax=xax, yax=yax, cell=0, cstar=0, cpoint=0)
points(biotyper.obj$dPCoA$l1[,xax], biotyper.obj$dPCoA$l1[,yax], cex=1, pch=16)

#BET
plot(biotyper.obj$BET$ls[,xax], biotyper.obj$BET$ls[,yax], type="n", xlab=paste("PC",xax), ylab=paste("PC",yax), main="between class")
s.potatoe(biotyper.obj$BET$ls, fac=as.factor(biotyper.obj$biotypes), xax=xax, yax=yax)
s.class(biotyper.obj$BET$ls, fac=as.factor(biotyper.obj$biotypes), grid=F, xax=xax, yax=yax, add.plot=TRUE, cell=0, cstar=0, cpoint=0)
points(biotyper.obj$BET$ls[,xax], biotyper.obj$BET$ls[,yax], cex=1, pch=16)

}


}



biotyper.ecosimulation=function(obs, k , ...) {

biotypes.uniform=sapply(ecosimulation$uniform, biotyper.data.frame, k=k)
print("uniform: biotyping done")

biotypes.normal=sapply(ecosimulation$normal, biotyper.data.frame, k=k)
print("normal: biotyping done")

biotypes.log.normal=sapply(ecosimulation$log.normal, biotyper.data.frame, k=k)
print("log.normal: biotyping done")

simul=list(biotypes.uniform=biotypes.uniform, biotypes.normal=biotypes.normal, biotypes.log.normal=biotypes.log.normal)
class(simul)<-"biotyper.obj.simul"
return(simul)

}


print.biotyper.obj.simul=function(biotyper.obj.simul){

sil.uniform=unlist(t(biotyper.obj.simul$biotypes.uniform)[,2])
sil.normal=unlist(t(biotyper.obj.simul$biotypes.normal)[,2])
sil.log.normal=unlist(t(biotyper.obj.simul$biotypes.log.normal)[,2])

silhouette = data.frame(sil.uniform=sil.uniform, sil.normal=sil.normal, sil.log.normal=sil.log.normal)
print(silhouette)

 }
 



summary.biotyper.obj.simul=function(biotyper.obj.simul){

sil.uniform=unlist(t(biotyper.obj.simul$biotypes.uniform)[,2])
sil.normal=unlist(t(biotyper.obj.simul$biotypes.normal)[,2])
sil.log.normal=unlist(t(biotyper.obj.simul$biotypes.log.normal)[,2])

silhouette = data.frame(sil.uniform=sil.uniform, sil.normal=sil.normal, sil.log.normal=sil.log.normal)
sum.silhouette = summary(silhouette)
print(sum.silhouette)

 }



test.biotypes=function(biotyper.obj, biotyper.obj.simul) {

if(!class(biotyper.obj)=="biotyper") {stop("this is not a biotyper object")}

if(!class(biotyper.obj.simul)=="biotyper.obj.simul") {stop("this is not a biotyper object(simulation)")}

sil=biotyper.obj$obs.silhouette

sil.uniform=unlist(t(biotyper.obj.simul$biotypes.uniform)[,2])
sil.normal=unlist(t(biotyper.obj.simul$biotypes.normal)[,2])
sil.log.normal=unlist(t(biotyper.obj.simul$biotypes.log.normal)[,2])

p.value.sil.uniform=1-(which(sort(c(sil,sil.uniform))==sil)/ length(c(sil,sil.uniform)))
p.value.sil.normal=1-(which(sort(c(sil,sil.normal))==sil)/ length(c(sil,sil.normal)))
p.value.sil.log.normal=1-(which(sort(c(sil,sil.log.normal))==sil)/ length(c(sil,sil.log.normal)))

sil.test=list(p.value.sil.uniform=p.value.sil.uniform,p.value.sil.normal=p.value.sil.normal,p.value.sil.log.normal=p.value.sil.log.normal)

silhouette.simul = list(sil.uniform=sil.uniform, sil.normal=sil.normal, sil.log.normal=sil.log.normal)
test=list(silhouette.test=sil.test, silhouette.simul=silhouette.simul, observed.sil=biotyper.obj$obs.silhouette)


class(test)<-"biotypes.test"

return(test)

}



print.biotypes.test=function(biotypes.test, ...) {


cat("Silhouette cluster validation : is species comunauty follow a continuun ?\n")
cat("here is alpha risk to say that observed value is different from expected value\n")
print(biotypes.test$silhouette.test)



}




plot.biotypes.test=function(biotypes.test, ...) {

stripchart(biotypes.test$observed.sil, pch=16, xlim=c(0,biotypes.test$observed.sil+(0.1*biotypes.test$observed.sil)), col="black", cex=2, xlab="Silhouette Coef")
stripchart(c(biotypes.test$silhouette.simul$sil.uniform), pch=16, col="blue", add=TRUE, method="jitter")
stripchart(c(biotypes.test$silhouette.simul$sil.normal), pch=16, col="red", add=TRUE, method="jitter")
stripchart(c(biotypes.test$silhouette.simul$sil.log.normal), pch=16, col="green", add=TRUE, method="jitter")
legend("topright", c("observed","uniform","normal","log.normal"), pt.cex=c(2,1,1,1), col=c("black","blue","red","green"), pch=16)

}








# if (is.null(best.species.number) ){
# best.species.number=2
# }

# rbind(
# dataframe.bet$co[order(dataframe.bet$co[,1])[1:best.species.number],],
# dataframe.bet$co[rev(order(dataframe.bet$co[,1]))[1:best.species.number],],
# dataframe.bet$co[order(dataframe.bet$co[,2])[1:best.species.number],],
# dataframe.bet$co[rev(order(dataframe.bet$co[,2]))[1:best.species.number],])-> best.species

# data=list(best.species=best.species, between=dataframe.bet, PCA=dataframe.pca, cluster=cluster)

# return(data)

# }

#########################################################
#                                                       #
#                 Cluster based on JSD distance         #
#                                                       #
#########################################################

cluster.biotypes=function(x, ...) UseMethod("cluster.biotypes")

cluster.biotypes.default=function(x,k) {

as.matrix(x)-> x.m
rownames(x)-> rownames(x.m)
dist.JSD(x.m)-> distance.JSD

cluster = as.vector(pam(as.dist(distance.JSD), k, diss=TRUE)$clustering)

return(cluster)


}

cluster.biotypes.dist=function(x,k) {


cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)

return(cluster)


}



#########################################################
#                                                       #
#                 CH index                              #
#                                                       #
#########################################################



CH.index=function(dataframe, distance.jsd=NULL, kvector=1:20, clusterSim=FALSE){

if(clusterSim) {
	require(clusterSim)
	kvector=unique(as.integer(kvector))
	nclusters=NULL

		for (k in kvector) { 
			if (k==1) {
				nclusters[k]=0
			}
		
			else {
				if(distance.jsd==NULL) stop("a distance matrix is needed")
				cluster=cluster.biotypes(distance.jsd, k)
				index.G1(t(dataframe),cluster,  d = distance.jsd, centrotypes = "medoids") -> nclusters[k]
			}

		}

	} else {
	
	nclusters=pamk(t(dataframe), criterion="ch", krange=kvector)$crit
	
	}



return(nclusters)

}






#########################################################
#                                                       #
#                 JSD Function Distance                 #
#                                                       #
#########################################################
 
dist.JSD <- function(inMatrix, pseudocount=10^(round(log10(min(as.matrix(inMatrix)[as.matrix(inMatrix) > 0])),0)-1), ...) {
	
	KLD <- function(x,y) sum(x *log(x/y))
	JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
	
	for(k in 1:matrixRowSize) { 
		for(j in 1:matrixColSize) {
			if (inMatrix[k,j] == 0) { inMatrix[k,j] <- pseudocount	}
		}
	}

	for(i in 1:matrixColSize) {
		for(j in 1:matrixColSize) { 
			JSD(as.vector(inMatrix[,i]), as.vector(inMatrix[,j]))->resultsMatrix[i,j] 
		}
	}
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
}

###################################################
#                                                 #
#  Template based communities  simulator          #
#                                                 #
###################################################


#x=matrix(rnorm(100000), nc=100)

ecosimulator=function(x,r) {

sites=names(x)

x=as.matrix(x)
pseudocount=10^(round(log10(min(as.matrix(x)[as.matrix(x) > 0])),0)-1)
x[which(x==0)]<- pseudocount
colnames(x)=sites

generator.unif=function(x){
n=dim(x)[2]
mx=apply(x,1,min)
sx=apply(x,1,max)
dx=data.frame(mx,sx)
g=function(x, c1, c2, n) {runif(n,min=x[c1], max=(x[c2]))}
simul=apply(dx,1,g,c1="mx",c2="sx", n)
rownames(simul)=sites
simul=as.data.frame(prop.table(t(simul),2))
return(simul)
}


generator.norm=function(x){
n=dim(x)[2]
mx=apply(x,1,mean)
sx=apply(x,1,sd)
dx=data.frame(mx,sx)
g=function(x, c1, c2, n) {rnorm(n,mean=x[c1], sd=(x[c2]))}
simul=apply(dx,1,g,c1="mx",c2="sx", n)
rownames(simul)=sites
simul=t(simul)
simul[which(simul < pseudocount)] <- pseudocount
simul=as.data.frame(prop.table(simul,2))
return(simul)
}


generator.lnorm=function(x){
n=dim(x)[2]

lmean=function(x) {mean(log(x[x>0]))}
lsd=function(x) {sd(log(x[x>0]))}

mx=apply(x,1,lmean)
sx=apply(x,1,lsd)
dx=data.frame(mx,sx)
g=function(x, c1, c2, n) {rlnorm(n,mean=x[c1], sd=(x[c2]))}
simul=apply(dx,1,g,c1="mx",c2="sx", n)
rownames(simul)=sites
simul=as.data.frame(prop.table(t(simul),2))
return(simul)
}

#TO DO : generator.permute=function(x) 

lx.unif=vector("list",r); for (i in 1:r) generator.unif(x)-> lx.unif[[i]]
lx.norm=vector("list",r); for (i in 1:r) generator.norm(x)-> lx.norm[[i]]
lx.lnorm=vector("list",r); for (i in 1:r) generator.lnorm(x)-> lx.lnorm[[i]]

lx=list(uniform=lx.unif, normal=lx.norm, log.normal=lx.lnorm)

class(lx)="ecosimulation"

return(lx)
}


print.ecosimulation=function(x,...){

cat(names(x),sep=", ")
cat(" = ",length(x)," lists of",length(x[[1]]), "simulated data\n")
cat("nsites =",dim(x[[1]][[1]])[2],", n[taxa or functions] =",dim(x[[1]][[1]])[1],"\n")

}




###################################################################
#                                                                 #
#                          Noise Removal                          #
#                                                                 #
###################################################################


noise.removal<-function(dataframe, percent=NULL, top=NULL, bysample=TRUE){


noise.removal.global <- function(dataframe, percent=NULL, top=NULL){
dataframe->Matrix

if (is.null(top)){

if (is.null(percent)) {
percent=100*(10^(mean(log10(rowSums(Matrix)/sum(rowSums(Matrix)))) + 1.96*sd(log10(rowSums(Matrix)/sum(rowSums(Matrix)))))) # this is percent
}
bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent  #this is percent ### noize filter
Matrix_1 <- Matrix[bigones,]
print(percent)
}

else {

bigones <- rev(order(apply(Matrix, 1, mean)))[1:top]
Matrix_1 <- Matrix[bigones,]


}


return(Matrix_1)


}


noise.removal.spec <- function(dataframe, percent=NULL, top=NULL){

dataframe->Matrix
bigones<-NULL

for(i in 1:dim(Matrix)[2]){
	
	if (is.null(top)){
		if (is.null(percent)) {
			percent<- 100*(10^(median(log10(Matrix[Matrix[,i]>0,i])) + 1.96*sd(log10(Matrix[Matrix[,i]>0,i]))))
			cat(percent,"\n")
		}
			tmp <- which(Matrix[,i] > percent/100)
			bigones<-unique(c(bigones,tmp))
	}
	else {
		tmp <- rev(order(Matrix[,i]))[1:top]
		bigones<-unique(c(bigones,tmp))
	}	

}

return(Matrix[bigones,])
print(percent)
}

if(bysample) { noise.removal.spec(dataframe, percent=percent, top=top)}

else {noise.removal.global(dataframe, percent=percent, top=top)}


}


s.potatoe<-function(dfxy, fac, xax = 1, yax = 2, col.border = seq(1:length(levels(fac)))+1, col.fill = seq(1:length(levels(fac)))+1,
shape = -0.5, open = "FALSE", ...) {
dfxy <- data.frame(dfxy)
opar <- par(mar = par("mar"))
par(mar = c(0.1, 0.1, 0.1, 0.1))
on.exit(par(opar))
x <- dfxy[, xax]
y <- dfxy[, yax]
for (f in levels(fac)) {
xx <- x[fac == f]
yy <- y[fac == f]
xc <- chull(xx, yy)
qui <- which(levels(fac) == f)
border <- col.border[qui]
col <- col.fill[qui]
xspline(xx[xc], yy[xc], shape = shape, open = open, border = border,
col = col, ...)
}
}
