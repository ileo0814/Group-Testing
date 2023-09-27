library(MASS)
library(BayesLogit)
library(Matrix)
library(rgdal)
library(raster)
library(spdep)
library(invgamma)
library(usmap)
library(ggplot2)
library(coda)

###########################
sim.id = 12 #3,4,11,12
###########################

#Running on cluster
path1 = "/home/rongjie/CAPC/"
path2 = "/home/rongjie/DT_Areal/"

#if running on my computer
if(getwd()=="C:/Users/Leo/Documents"){
  path1 = 'C:/Users/Leo/OneDrive - University of South Carolina/Dr Self/CAPC/'
  path2 = "C:/Users/Leo/OneDrive - University of South Carolina/Dr Self/BIOS 890 Independent Study/results/DT_Areal/"
}

fips <- read.table(paste0(path1,'Data/County FIPS Codes.txt'),sep = "\t",header = T)
fips.sc <- fips$FIPS[fips$State=="SC"]
length(fips.sc)

#Load simulation scenarios
df.sim <- read.csv(paste0(path2,"sim.id.csv"))

#Load (centered) true spatial random effects
gamma.true <- read.table(paste0(path2,'gamma.true.txt'))
gamma.true <- gamma <- as.vector(gamma.true$V1) 

########################### Simulation Configurations ##########################

#1. Sample size
N.seq <- c(1656,3312)
############### (Choose one) #####################
  # N <- 1656
  # N <- 3312  
N <- df.sim$N[df.sim$sim.id==sim.id] #automatically choose N

#####################################

#2. Prevalence
p.seq <- c(0.05,0.09,0.15,0.5)
p.seq <- c(0.05,0.09,0.15)

############### (Choose one) ###################
  # beta <- c(-4.1,-1.5,0.2) #p=0.05
  # beta <- c(-3.2,-1.4,0.2) #p=0.09
  # beta <- c(-2.7,-1.6,0.3) #p=0.15
  # beta <- c(-1.3,-1.5,0.2) #p=0.30

p.sim <- df.sim$p[df.sim$sim.id==sim.id]
if(p.sim == 0.05){
  beta <- c(-4.1,-1.5,0.2) #p=0.05
} else if(p.sim == 0.09){
  beta <- c(-3.2,-1.4,0.2) #p=0.09
} else if(p.sim == 0.15){
  beta <- c(-2.7,-1.6,0.3) #p=0.15
} else if(p.sim == 0.30){
  beta <- c(-1.3,-1.5,0.2) #p=0.30
}


#####################################

#3. Pool size
n.size.seq <- c(6,3)
n.size.seq <- c(1)

################ (Choose one) ####################
  # n.size <- 6
  # n.size <- 3
  n.size <- df.sim$n.size[df.sim$sim.id==sim.id]

  if(n.size == 1){
    stop("This R script is for DT; use IT_Areal.R instead")
  }
    
#####################################

#4. Pool structure
  #4.1 from same county
  #4.2 randomly assigned
  G.mat.seq <- c("sequential","random") 
  
#####################################
  #G.mat.config <- "sequential"
  G.mat.config <- "random"

#####################################  

  
#5. Covariate estimation/Prevalence estimation
 param.list <- list(n.size = n.size.seq,
                   p = p.seq,
                   N = N.seq)

 wd= getwd()
 wd
 tuning.list <- expand.grid(param.list)
 tuning.list
 #write.csv(tuning.list,paste0("C:/Users/Leo/OneDrive/huangleo Leo Huang/USC Files/2022 Spring (Current)/BIOS 890 Independent Study/results/sim.list.csv"))
 
 
 # #Assign FIPS code (25 counties) so that each pool only contains test sample from individuals from same county
 n.county <- length(fips.sc) #46
 fips.sub <- fips.sc[1:n.county]
 fips.ind <- sort(rep(fips.sub,length.out=N))
 
 #Spatial random effects (from CAR model)
 #Create W matrix and D matrix
 shp <- readOGR(paste0(path1,"Data/cb_2018_us_county_500k.shp")[1]) #Need to also have .dbf and .shx in the folder w/ .shp
 
 #Subset shp with only the fips code appearing in fips.sub
 shp <- shp[shp$GEOID%in%fips.sub,]
 
 
 ########################### Data Generation ################################
 
 ###### Individual Testing (IT) ###### 
 
 #Using IT as Baseline testing
 n.sim <- 100
 n.iter <- 15000 #5000
 #n.iter <- 100
 
 #Pool testing
 #3x3x100
 #Split the individual testing data into batches
 n.pool <- N/n.size #300
 
 #Sample X from normal distribution and bernoulli distribution; N = 5000 individual
 X <- matrix(NA, nrow = N, ncol = 3)
 X[,1] <- 1
 X[,2] <- rnorm(N,0,1)
 X[,3] <- rbinom(N,1,0.5)
 
 #Specify sensitivity P(+|D) & specificity P(-|D-)
 sens <- 0.98
 spec <- 0.98
 
 #Hyperparameter for sensitivity and specificity
 a.sens <- 1 #0.095
 b.sens <- 1 #0.005
 a.spec <- 1 #0.095
 b.spec <- 1 #0.005
 
 

######################################### Spatial random effects ############################################

# n.rows <- 5
# n.cols <- 5
# n.units <- n.rows*n.cols #25
# grid.raster <- raster(ncols = n.cols,nrows = n.rows, xmn = 0, xmx = 0.1, ymn = 0, ymx = 0.1)
# shp.grid <- rasterToPolygons(grid.raster)
# shp.grid <- spTransform(shp.grid,CRS("+init=epsg:32610"))
# neighbours <- poly2nb(shp.grid)
# W <- Matrix(as(nb2mat(neighbours,style = 'B',zero.policy = T), "Matrix"),sparse = T)
# shp.grid$ID <- seq(1,n.units)
# individual.grid <- sort(rep(shp.grid$ID,length.out=N))
# D <- Diagonal(n = ncol(W),colSums(W))

neighbours <- poly2nb(shp)
W <- Matrix(as(nb2mat(neighbours,style = 'B',zero.policy = T), "Matrix"),sparse = T)
D <- Diagonal(n = ncol(W),colSums(W))

#Hyperparameters
a.tau <- 1 #0.01
b.tau <- 1 #0.01

#Prior
# tau.sq <- rinvgamma(1,shape = a.tau,scale = b.tau)
# tau.sq.true <- tau.sq

tau.sq <- 0.1
rho <- 1 #0.999
tau.sq.true <- tau.sq

# write.table(gamma,
# paste0("C:/Users/Leo/OneDrive/huangleo Leo Huang/USC Files/2022 Spring (Current)/BIOS 890 Independent Study/results/gamma.true.txt"),
# row.names = F,
# col.names = F
# )

# set.seed(123)   #save a set of centered true gamma
# gamma <- mvrnorm(1,rep(0,n.county),tau.sq*solve(D-rho*W))    #25 dim
# gamma <- gamma - mean(gamma)
# gamma.true <- gamma
# gamma.true
# range(gamma)

#Create G.mat (N by n.county) to map (N) individuals to (n.county) counties 
  
  #sequentially assigned
  G.mat.sequential <- matrix(as.numeric(sapply(unique(fips.ind),"==",fips.ind)),nrow = N,ncol=n.county) 
  G.mat.sequential <-  as(G.mat.sequential, "sparseMatrix")

  #randomly assigned
  G.mat.random <- matrix(NA,N,n.county) 
  full.set <- temp.set <- seq(1:N)
  for(i in 1:n.county){
    ind <- sample(temp.set,size = N/n.county,replace = F)
    G.mat.random[,i] <- as.numeric(full.set%in%ind)
    temp.set <- temp.set[-which(temp.set%in%ind)]
  }
  rm(ind);rm(full.set);rm(temp.set)
  G.mat.random <-  as(G.mat.random, "sparseMatrix")
  
  #For grid
  #G.mat.grid <- matrix(as.numeric(sapply(unique(individual.grid),"==",individual.grid)),nrow = N,ncol=n.county) 
  #G.mat.grid <-  as(G.mat, "sparseMatrix")

if(G.mat.config=="sequential"){
  G.mat <- G.mat.sequential
} else if(G.mat.config=="random"){
  G.mat <- G.mat.random
} 
  

    
#Compute p via logistic(beta, X)
p <- as.vector(exp(X%*%beta+G.mat%*%gamma)/(1+exp(X%*%beta+G.mat%*%gamma)))        #add spatial info
mean(p)

#Specify the true status Y.tilde
Y.tilde <- rbinom(N,1,p)
Y.tilde.true <- Y.tilde



#Generate the individual test status
Y <- rep(NA,N)
Y[Y.tilde==1] <- rbinom(length(Y.tilde[Y.tilde==1]),1,sens)
Y[Y.tilde==0] <- rbinom(length(Y.tilde[Y.tilde==0]),1,1-spec)

#Compute Bias
#Compute the sample standard deviation of the estimates (SSD)
#Compute the probability of inclusion (PI)


###### Dorfman Testing (DT) ######

#Create a list to save the indices
  #Create a list for each array to record individual's id (i.e. array id to individual id's)  
  pool.list <- list()
  for (i in 1:(N/n.size)){
    #pool.list[[i]] <- seq(i*3-2,i*3)
    pool.list[[i]] <- seq(i*n.size-(n.size-1),i*n.size)
  }
  
  #Create a list for each individual to record the pool being contributed to (i.e. individual id to pool id's) 
  ind.list <- list()
  for(i in 1:N){
    ind.list[[i]] <- which(lapply(lapply(lapply(pool.list,"%in%",i),as.numeric),sum)==1)
  }


#True status for each pool
Z.tilde <- c() #n.pool by 1
for(i in 1:n.pool){
    Z.tilde[i] <- ifelse(sum(Y.tilde[unlist(pool.list[i])])>0,1,0)
}
  
#Generate the test status (Dilution effect?)
# Z <- rep(NA, n.pool)
# Z[Y.tilde.pt==1] <- rbinom(length(Y.tilde.pt[Y.tilde.pt==1]),1,sens)
# Z[Y.tilde.pt==0] <- rbinom(length(Y.tilde.pt[Y.tilde.pt==0]),1,1-spec)

#Or equivalently
Z <- rep(NA, n.pool)
Z <- rbinom(n.pool,1,sens*Z.tilde+(1-spec)*(1-Z.tilde))
Z_0 <- Z
length(Z)
sum(Z)


if(n.size!=1){
  #Retesting (append Z with Y's in the positive pools)
  Z <- c(Z,Y[unlist(pool.list[which(Z==1)])])
  length(Z)
  sum(Z)
  
  pool.list <- append(pool.list,unlist(pool.list[which(Z==1)])) 
  
  write.table(Z_0,"Z_0.csv",row.names = F,col.names = F)  
  write.table(Z,"Z.csv",row.names = F,col.names = F)  
  
  ind.list <- list()
  for(i in 1:N){
    ind.list[[i]] <- which(lapply(lapply(lapply(pool.list,"%in%",i),as.numeric),sum)==1)
  }
  
  n.pool <- length(Z)
}

n.pool

#Gibbs sampler for beta and w
n.b <- length(beta)

#priors  
b.m <- rep(0,n.b)
b.cov <- diag(10,n.b)
b.cov.inv <- solve(b.cov)

#Initial values 
b <- matrix(NA,n.iter,n.b)
w <- matrix(NA,n.iter,N)
b[1,] <- rep(0.5,n.b)
w[1,] <- rep(0.5,n.b)
Y.tilde.mcmc <- matrix(0,n.iter,N)
Y.tilde <- rep(0,N)

sej.mcmc <- spj.mcmc <- matrix(NA,n.iter,n.pool)
sej.mcmc[1,] <- sej <- rep(sens,n.pool) #pool specific sensitivity P(+test|+true status)
spj.mcmc[1,] <- spj <- rep(spec,n.pool) #pool specific specificity P(-test|-true status)

sej <- rep(sens,n.pool)
spj <- rep(spec,n.pool)

#initial gamma
rho <- 0.999
gamma <- mvrnorm(1,rep(0,n.county),tau.sq*solve(D-rho*W))
gamma <- gamma - mean(gamma)
rho <- 1


gamma.mcmc  <- matrix(NA,n.iter,n.county)
gamma.mcmc[1,] <- gamma

tau.sq.mcmc <- rep(NA,n.iter)
tau.sq.mcmc[1] <- tau.sq

p.mcmc <- matrix(NA,n.iter,n.county)
p.mcmc[1,] <- as.vector(exp(rep(b[1,1],n.county) + gamma)/(1+exp(rep(b[1,1],n.county) + gamma)))

p.true <- as.vector(exp(rep(beta[1],n.county) + gamma.true)/(1 + exp(rep(beta[1],n.county) + gamma.true)))

n.size
N
mean(p)

run.start <- Sys.time()
for (t in 2:n.iter){
  #Step 1 Sample Y.tilde (Pool testing)  
    #p1
    p1 <- p0 <- c()
    for (i in 1:N){ 
      # j <- ind.list[[i]]
      # eta <- 1/(1+exp(-X[i,]%*%beta))
      # #p1
      # p1[i] <- eta*sej[j]^Z[j]*(1-sej[j])^(1-Z[j])
      # 
      # sij <- sum(Y.tilde[pool.list[[j]][-which(pool.list[[j]]==i)]])
      # #p0
      # p0[i] <- (1-eta)*(sej[j]^Z[j]*(1-sej[j])^(1-Z[j]))^(as.numeric(sij>0))*
      #   ((1-spj[j])^Z[j]*spj[j]^(1-Z[j]))^(as.numeric(!sij>0))
    
      ind.j <- ind.list[[i]]
      
      #eta <- 1/(1+exp(-(X[i,]%*%beta+(G.mat%*%gamma)[i])))
      eta <- 1/(1+exp(-(X[i,]%*%b[t-1,]+(G.mat%*%gamma)[i])))
      
      #p1
      p1[i] <- as.numeric(eta)*prod(sej[ind.j]^Z[ind.j]*(1-sej[ind.j])^(1-Z[ind.j]))
      #p0
      sij <- sum(Y.tilde[unlist(pool.list[ind.j][1])])-Y.tilde[i]
      if(length(pool.list[ind.j])==2){
        sij <- c(sij,0)
      }
      p0[i] <- (1-as.numeric(eta))*
        prod((sej[ind.j]^Z[ind.j]*(1-sej[ind.j])^(1-Z[ind.j]))^(as.numeric(sij>0))*
               ((1-spj[ind.j])^Z[ind.j]*spj[ind.j]^(1-Z[ind.j]))^(as.numeric(!sij>0)))
      
      Y.tilde[i] <- rbinom(1,1,prob = p1[i]/(p1[i]+p0[i]))
      }    
    
    # p.Y.tilde <- p1/(p1+p0)
    # Y.tilde <- sapply(p.Y.tilde,rbinom,n=1,size=1)
    Y.tilde.mcmc[t,] <- Y.tilde
    
  #Step 2 Compute Z.tilde based on Y.tilde
    Z.tilde <- NULL
    for(i in 1:n.pool){
      Z.tilde[i] <- ifelse(sum(Y.tilde[unlist(pool.list[i])])>0,1,0)
    }
    
  # #Step 3 Sample sens and spec (For IT only)
  #   #Sensitivity
  #    sej.mcmc[t,] <- sej
  # 
  #   #Specificity
  #    spj.mcmc[t,] <- spj
  
      
  #Step 3 Sample sens and spec
    #Sensitivity
    aem <- a.sens + sum(Z*Z.tilde)
    bem <- b.sens + sum((1-Z)*Z.tilde)
    sej <- rep(rbeta(1,aem,bem),n.pool)
    sej.mcmc[t,] <- sej

    #Specificity
    apm <- a.spec + sum((1-Z)*(1-Z.tilde))
    bpm <- b.spec + sum(Z*(1-Z.tilde))
    spj <- rep(rbeta(1,apm,bpm),n.pool)
    spj.mcmc[t,] <- spj


        
  #Step 4 Sample beta
    omg <- diag(w[t-1,])
    M <- t(X)%*%omg%*%X+b.cov.inv
    k <- Y.tilde - 1/2
    h <- k/diag(omg)
    Q <- t(t(h)%*%omg%*%X - t(gamma)%*%t(G.mat)%*%omg%*%X + t(b.m)%*%b.cov.inv) 
    b[t,] <- as.vector(mvrnorm(n = 1,chol2inv(chol(M))%*%Q,chol2inv(chol(M)))) 
    
  ##Step 5 Sample w    
    w[t,] <- apply(X = X%*%b[t,]+G.mat%*%gamma,1,FUN = rpg,num=1,h=1)
    #update omg and h 
    omg <- diag(w[t,])
    h <- k/diag(omg)
    
  #Step 6 Sample gamma
    R <- t(G.mat)%*%omg%*%G.mat + tau.sq^(-1)*(D-rho*W)
    R.inv <- solve(R)
    U <- t(t(h)%*%omg%*%G.mat - t(b[t,])%*%t(X)%*%omg%*%G.mat)  
    gamma <- mvrnorm(1,R.inv%*%U,R.inv)
    gamma <- gamma - mean(gamma)
    gamma.mcmc[t,] <- as.vector(gamma)
    
  #Step 7 Sample tau.sq
    tau.sq <- rinvgamma(1,shape = a.tau + n.county/2, rate = as.numeric(1/2*t(gamma)%*%(D-rho*W)%*%gamma+b.tau))
    tau.sq.mcmc[t] <- tau.sq
    
    
    
  #Step 8 Compute estimated prevalence for each county
    eta.est <- rep(b[t,1],n.county) + gamma
    p.est <- as.vector(exp(eta.est)/(1 + exp(eta.est)))
    p.mcmc[t,] <- p.est
    
    
    print(t)
} 

#save.image("C:/Users/Leo/OneDrive - University of South Carolina/R image/Array_testing/pool052722.Rdata")
#load("C:/Users/Leo/OneDrive - University of South Carolina/R image/Array_testing/pool052722.Rdata")

run.end <- Sys.time()
run.time <- difftime(run.end, run.start, units='mins')
run.time 


#Trace-plots (full samples)
  #Beta
  for(i in 1:length(beta)){
    plot(b[,i],type = "l", main = paste0("Beta_",i))
  }
  
  #Sensitivity
  plot(sej.mcmc[,1],type = "l")

  #Specificity
  plot(spj.mcmc[,1],type = "l")
  
  #Gamma
  for(i in 1:length(gamma)){
    plot(gamma.mcmc[,i],type = "l", main = paste0("Gamma_",i))
  }
  
  #Tau.sq
  plot(tau.sq.mcmc,main="tau.sq",type = "l") 



#Remove burn-in and compute geweke statistics
b.pred <- b[-(1:(0.5*n.iter)),]
b.geweke <- rep(NA,ncol(b.pred))
for(i in 1:ncol(b.pred)){
  b.geweke[i] <- unlist(geweke.diag(b.pred[,i]))["z.var1"]
}



Y.tilde.pred <- Y.tilde.mcmc[-(1:(0.5*n.iter)),]

sej.pred <- sej.mcmc[-(1:(0.5*n.iter)),]
sej.geweke <- unlist(apply(sej.pred,2,geweke.diag)[1])["z.var1"]


spj.pred <- spj.mcmc[-(1:(0.5*n.iter)),]
spj.geweke <- unlist(apply(spj.pred,2,geweke.diag)[1])["z.var1"]


gamma.pred <- gamma.mcmc[-(1:(0.5*n.iter)),]
gamma.geweke <- rep(NA,ncol(gamma.pred))
for(i in 1:ncol(gamma.pred)){
  gamma.geweke[i] <- unlist(geweke.diag(gamma.pred[,i]))["z.var1"]
}


tau.sq.pred <- tau.sq.mcmc[-(1:(0.5*n.iter))]
tau.sq.geweke <- unlist(geweke.diag(tau.sq.pred))["z.var1"]

p.pred <- p.mcmc[-(1:(0.5*n.iter)),]
p.geweke <- rep(NA,ncol(p.pred))
for(i in 1:ncol(p.pred)){
  p.geweke[i] <- unlist(geweke.diag(p.pred[,i]))["z.var1"]
}


#b MCMC samples 
b.hat <- apply(b.pred,2,mean)
b.upp <- apply(b.pred,2,quantile,probs = 0.975) 
b.low <- apply(b.pred,2,quantile,probs = 0.025) 
b.sd <- apply(b.pred,2,sd)


b.low
b.hat
b.upp 
beta

b.diff <- beta - b.hat
b.cover <- ifelse(beta>=b.low&beta<=b.upp,1,0)
b.out <- rbind(b.diff,b.cover,b.hat,b.sd)
b.out



#p MCMC samples
p.hat <- apply(p.pred,2,mean)
p.upp <- apply(p.pred,2,quantile,probs = 0.975) 
p.low <- apply(p.pred,2,quantile,probs = 0.025) 
p.sd <- apply(p.pred,2,sd)

p.diff <- p.true - p.hat
p.cover <- ifelse(p.true>=p.low&p.true<=p.upp,1,0)
p.out <- rbind(p.diff,p.cover,p.hat,p.sd)




#Y.tilde MCMC samples
Y.tilde.hat <- apply(Y.tilde.pred,2,mean)
Y.tilde.upp <- apply(Y.tilde.pred,2,quantile,probs = 0.975)
Y.tilde.low <- apply(Y.tilde.pred,2,quantile,probs = 0.025)

Y.tilde.diff <- Y.tilde.true - Y.tilde.hat
Y.tilde.cover <- ifelse(Y.tilde.true>=Y.tilde.low&Y.tilde.true<=Y.tilde.upp,1,0)
Y.tilde.out <- rbind(Y.tilde.diff,Y.tilde.cover,Y.tilde.true)
Y.tilde.out

plot(Y.tilde.true-colMeans(Y.tilde.pred),main="Residual plot_Y.tilde") 


#Sensitivity MCMC samples
sej.hat <- apply(sej.pred,2,mean)
sej.upp <- apply(sej.pred,2,quantile,probs=0.975)
sej.low <- apply(sej.pred,2,quantile,probs=0.025)
sej.sd <- sd(sej.pred[,1])

sej.diff <- sens - sej.hat[1]
sej.cover <- ifelse(sens>=sej.low[1]&sens<=sej.upp[1],1,0)
sej.out <- c(sej.diff,sej.cover,sej.hat[1],sej.sd)
sej.out

plot(sej.pred[,1],main="sej")

#Specificity MCMC samples
spj.hat <- apply(spj.pred,2,mean)
spj.upp <- apply(spj.pred,2,quantile,probs=0.975)
spj.low <- apply(spj.pred,2,quantile,probs=0.025)
spj.sd <- sd(spj.pred[,1])

spj.diff <- spec - spj.hat[1]
spj.cover <- ifelse(spec>=spj.low[1]&spec<=spj.upp[1],1,0)
spj.out <- c(spj.diff,spj.cover,spj.hat[1],spj.sd)
spj.out

plot(spj.pred[,1],main="spj")


#gamma MCMC samples
  gamma.hat <- apply(gamma.pred,2,mean)
  gamma.hat <- gamma.hat 
  gamma.upp <- apply(gamma.pred,2,quantile,probs=0.975) 
  gamma.low <- apply(gamma.pred,2,quantile,probs=0.025) 
  gamma.sd <- apply(gamma.pred,2,sd)
  
  gamma.diff <- gamma.true - gamma.hat
  gamma.cover <- ifelse(gamma.true>=gamma.low&gamma.true<=gamma.upp,1,0)
  gamma.out <- rbind(gamma.diff,gamma.cover,gamma.hat,gamma.sd,gamma.true)
  gamma.out

  plot(gamma.true - colMeans(gamma.pred),main="Residual plot_gamma.pred")
  
  for(i in 1:length(gamma)){
    plot(gamma.pred[,i],type = "l",main=paste0("Gamma ",i))
  }

#tau.sq MCMC samples
tau.sq.hat <- mean(tau.sq.pred)
tau.sq.upp <- quantile(tau.sq.pred,0.975)
tau.sq.low <- quantile(tau.sq.pred,0.025)
tau.sq.sd <- sd(tau.sq.pred)

tau.sq.diff <- tau.sq.true - tau.sq.hat
tau.sq.cover <- ifelse(tau.sq.true>=tau.sq.low&tau.sq.true<=tau.sq.upp,1,0)
tau.sq.out <- c(tau.sq.diff,tau.sq.cover,tau.sq.hat,tau.sq.sd)
tau.sq.out

plot(tau.sq.pred,main="tau.sq",type = "l") 


write.table(b.out,"b.out.txt",row.name=F,col.name=F)         
write.table(Y.tilde.out,"Y.tilde.out.txt",row.name=F,col.name=F)         
write.table(sej.out,"sej.out.txt",row.name=F,col.name=F)
write.table(spj.out,"spj.out.txt",row.name=F,col.name=F)
write.table(gamma.out,"gamma.out.txt",row.name=F,col.name=F) 
write.table(tau.sq.out,"tau.sq.out.txt",row.name=F,col.name = F)
write.table(p.out,"p.out.txt",row.name=F,col.name=F)  

write.table(b.geweke,"b.geweke.txt",row.name=F,col.name=F)
write.table(sej.geweke,"sej.geweke.txt",row.name=F,col.name=F)
write.table(spj.geweke,"spj.geweke.txt",row.name=F,col.name=F)
write.table(gamma.geweke,"gamma.geweke.txt",row.name=F,col.name=F) 
write.table(tau.sq.geweke,"tau.sq.geweke.txt",row.name=F,col.name = F)
write.table(p.geweke,"p.geweke.txt",row.name=F,col.name=F) 




#############################################################################################################################################################





