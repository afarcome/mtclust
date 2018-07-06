# load libraries
require(tclust)
require(robustbase)
require(lattice)
source("m_t_clust_model.r")

# load data
data(thyroid)
data = thyroid

#- explorative anaylsis: some nice plots

#- the pairs plot colored according to the true classification labels
real_clusters = as.factor(data[,1])
X11()
pairs(data[,-1],col=as.numeric(real_clusters)+1)
parallelplot(data[,-1],groups = real_clusters)

#- the ctl curves for fixing a guess value for the trimming level and the number of clusters
curves = ctlcurves(data[,-1])
X11()
plot(curves)
points(x=c(0.05,0.05),y=c(-3500,-2140),type="l",lty="dashed")


Mtclust.mod = m_t_clust(data[,-1],K=3,alpha.fixed = 0.05, modelNames ="EVI")


Mahal_matrix = matrix(NA,nrow=length(which(Mtclust.mod$classification==0)),
                      ncol=ncol(Mtclust.mod$parameters$mean))

data_out = data[which(Mtclust.mod$classification==0),-1]

for(j in 1:ncol(Mtclust.mod$parameters$mean)){
  
  Mahal_matrix[,j] = mahalanobis(data_out,Mtclust.mod$parameters$mean[,j],Mtclust.mod$parameters$variance$sigma[,,j])
  
}

Mahalan_vec = apply(Mahal_matrix,1,min)

## Estimate the trimming level starting from 0 trimming

K_vec = c(1,2,3,4,5); p=ncol(data[,-1]); 
mod_set = matrix(list(),nrow=length(names),ncol=length(K_vec));BIC = c()

names =  c("EII","VII","EEI","VEI", "EVI","VVI", "EEE", "VEE","EVE","EEV","VVE","VEV",  "EVV","VVV")
BIC_mat = matrix(nrow=length(names),ncol=length(K_vec))

alpha_choose = 0

for(i in 1:length(names)){
  
  for(j in 1:length(K_vec)){
    
    
    mod = try(m_t_clust(data[,-1],K=K_vec[j],alpha.fixed=alpha_choose,modelNames = names[i],iter.max=50,restr.factor=50,initial.values = NULL,alpha.init = 0.10,restr.fact.init = 50))
    
    if(!inherits(mod,"try-error")){
      
      K            = K_vec[j]
      n.parameters = c(1,K,p,K+(p-1), 1+K*(p-1), p*K,(p*(p+1)/2), K+p-1+(p*(p-1)/2), 1+K*(p-1)+(p*(p-1)/2), p+(K*p*(p-1)/2), K*p+(p*(p-1))/2, K+p-1+K*p*(p-1)/2, 1+K*(p-1) + K*p*(p-1) /2,K*(p*(p+1)/2))+((K*p)+K-1)
      
      mod_set[[i,j]] = mod
      loglik         = mod_set[[i,j]]$loglikelihood[length(mod_set[[i,j]]$loglikelihood)]
      BIC_mat[i,j]   = n.parameters[i]*log(nrow(data[,-1]))-2*loglik 
    }else{print(i);print(j)}
  }
  
}

print(BIC_mat)

which(BIC_mat==min(BIC_mat,na.rm = T),arr.ind=T)

row.names(BIC_mat) = names

print(xtable(BIC_mat_inserted))
write.xlsx(BIC_mat_inserted,file="tab.xlsx")

mod_selected = mod_set[[12,3]]
mod_mine = mod_selected

## Choose the trimming level needed

dd = list()

for( i in 1:length(unique(mod_mine$classification[mod_mine$classification!=0]))){
  
  dd[[i]] =  covMcd(data[which(mod_mine$classification==i),-1],alpha=.75)
  
}

mahal = list()

for(i in 1:length(dd)){
  
  mahal[[i]] = mahalanobis(data[which(mod_mine$classification==i),-1],t(dd[[i]]$center),dd[[i]]$cov )
  
}

test = c()

for(i in 1:length(mahal)){
  
  vec = c(mahal[[i]])
  test[i] = length(which(vec>qchisq(.975,5)))
  
}

alpha_chosen = sum(test)/nrow(data)
# the obtained value look very close to the one fixed by looking at the ctl curves
print(alpha_chosen)

### Artificially contaminate data ###
rr = apply(data[,-1],2,range)
set.seed(1305)
ncont = 10
c=6
out = cbind(c(runif(ncont,rr[1,1]+rr[1,1]/c,rr[2,1]+rr[2,1]/c)),
            c(runif(ncont,rr[1,2]+rr[1,2]/c,rr[2,2]+rr[2,2]/c)),
            c(runif(ncont,rr[1,3]+rr[1,3]/c,rr[2,3]+rr[2,3]/c)),
            c(runif(ncont,rr[1,4]+rr[1,4]/c,rr[2,4]+rr[2,4]/c)),
            c(runif(ncont,rr[1,5]+rr[1,5]/c,rr[2,5]+rr[2,5]/c)) )  

out = round(out,1)            
out = data.frame(c(rep("out",ncont)),out)
colnames(out) = colnames(data)
data_cont = merge(data,out,all=TRUE)
nrow(data_cont)
pairs(data[,-1])
pairs(data_cont[,-1])

# Launch the robust version of the Mclust

Rmclust = Mclust(data_cont[,-1],G=3,
                 initialization = list(noise=c(sample(x=c(seq(1:nrow(data_cont))),size=ncont*2))))
Rmclust$loglik
table(Rmclust$classification)

# Launch the mtclust 

Rmtclust =m_t_clust(data_cont[,-1],K=3,modelNames = "EVI",alpha.fixed = .10,iter.max=50)
Rmtclust$loglikelihood
table(Rmtclust$classification)

# Launch the tclust

Rtclust = tclust(data_cont[,-1],k=3,alpha=.10)
