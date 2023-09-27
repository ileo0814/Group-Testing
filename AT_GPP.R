library(MASS)
library(BayesLogit)
library(Matrix)
library(rgdal)
library(raster)
library(spdep)
library(invgamma)
library(usmap)
library(ggplot2)
library(sf)
library(coda)
#library(mapview)
#library(ggmap)
#library(RColorBrewer)


###########################
sim.id = 116
###########################

n.sim <- 100
n.iter <- 15000 #5000


#Running on cluster
path1 = "/home/rongjie/CAPC/"
path2 = "/home/rongjie/DT_Areal/"

#if running on my computer
if(getwd()=="C:/Users/Leo/Documents"){
  path1 = 'C:/Users/Leo/OneDrive - University of South Carolina/Dr Self/CAPC/'
  path2 = "C:/Users/Leo/OneDrive - University of South Carolina/Dr Self/BIOS 890 Independent Study/results/DT_Areal/"
  n.iter = 100
}

#Load simulation scenarios
df.sim <- read.csv(paste0(path2,"sim.id.csv"))

sim.ind <- which(df.sim$sim.id==sim.id)

fips <- read.table(paste0(path1,'Data/County FIPS Codes.txt'),sep = "\t",header = T)
fips.sc <- fips$FIPS[fips$State=="SC"]
length(fips.sc)




########################### Simulation Configurations ##########################

#1. Sample size
N.seq <- c(1656,3312)

################################  
N <- df.sim$N[sim.ind] #automatically choose N
################################ 

#2. Generating surface from
gen.from <- df.sim$surface[sim.ind] 


#3. Prevalence
p.seq <- c(0.09,0.3)

p.sim <- df.sim$p[sim.ind]

if(gen.from=="trig"){
  if(p.sim==0.05){
    beta.true <- c(-2.6,-1.3,0.2) #p=0.05
  } else if(p.sim==0.09){
    beta.true <- c(-2.0,-1.3,0.7) #p=0.09
  } else if(p.sim==0.15){
    beta.true <- c(-1.5,-1.0,0.2) #p=0.15
  } else if(p.sim==0.30){
    beta.true <- c(-0.6,-0.3,0.3) #p=0.30
  }
}else if(gen.from=="GPP"){
  if(p.sim==0.05){
    beta.true <- c(-2.6,-1.3,0.2) #p=0.05
  } else if(p.sim==0.09){
    beta.true <- c(-2.0,-1.2,0.2) #p=0.09
  } else if(p.sim==0.15){
    beta.true <- c(-1.5,-1.0,0.2) #p=0.15
  } else if(p.sim==0.30){
    beta.true <- c(-0.8,-0.3,0.2) #p=0.30
  }
}


#4. Pool size
n.size.seq <- c(6,3)
#####################################

n.size <- df.sim$n.size[sim.ind]

if(n.size == 1){
  stop("This R script is for AT; use IT_GPP.R instead")
}

if(df.sim$Protocol[df.sim$sim.id==sim.id]!="AT"){
  stop("The sim.id chosen is NOT for AT")
}

#####################################

#5. Number of knots m

#####################################

m <- df.sim$m[sim.ind]

#####################################

#6. Sensitivity and specificity
#6.1 Sensitivity
sej.seq <- as.numeric(unlist(lapply(as.character(df.sim$sej),FUN=strsplit,",")[sim.ind]))
#6.2 Specificity
spj.seq <- as.numeric(unlist(lapply(as.character(df.sim$spj),FUN=strsplit,",")[sim.ind]))


#7. Pool structure
#7.1 from same county
#7.2 randomly assigned
G.mat.seq <- c("sequential","random")

#####################################
#G.mat.config <- "sequential"
G.mat.config <- "random"
#####################################

#8. Covariate estimation/Prevalence estimation
param.list <- list(n.size = n.size.seq,
                   p = p.seq,
                   N = N.seq,
                   protocol = c("IT","DT","AT"))

wd= getwd()
wd
tuning.list <- expand.grid(param.list)[,c("protocol","N","p","n.size")]


#write.csv(tuning.list,"tuning.list.csv")



########################### Data Generation ################################

###### Individual Testing (IT) ######
#Using IT as Baseline testing
n.b <- length(beta.true)

if(gen.from=="trig"){
  #(Option 1) From trig function
  X <- matrix(NA, nrow = N, ncol = n.b-1)
  X[,1] <- rbinom(N,1,0.5)
  X[,2] <- rbinom(N,1,0.5)
  b0 <- beta.true[1]
  beta <- beta.true[-1]
  
}else if(gen.from=="GPP"){
  #(Option 2) From GPP
  X <- matrix(NA, nrow = N, ncol = n.b)
  X[,1] <- 1
  X[,2] <- rbinom(N,1,0.5)
  X[,3] <- rbinom(N,1,0.5)
  beta <- beta.true
}  

n.array <- N/n.size^2 #100
n.pool <- 2*n.array*n.size


#Specify sensitivity P(+|D) & specificity P(-|D-)
#Allowing multiple assays (sensitivities & specificities)
n.assay <- df.sim$n.assay[sim.ind]

sens <- sej.seq[1:n.assay]
spec <- spj.seq[1:n.assay]


df.assay <- data.frame(id = c(1:n.assay), sens = sens, spec = spec) #sens.1 = 0.95; sens.2 = 0.98; spec.1 = 0.95; spec.2 = 0.98 spec.2 = 0.98

# for(i in 1:n.assay){  
#    assign(paste0("sens.",df.assay$id[i]),df.assay$sens[i])
#    assign(paste0("spec.",df.assay$id[i]),df.assay$spec[i])
# }


#sens <- 0.98
#spec <- 0.98


#Hyperparameter for sensitivity and specificity
a.sens <- 1 #0.095
b.sens <- 1 #0.005
a.spec <- 1 #0.095
b.spec <- 1 #0.005




######################################### Spatial random effects ############################################

# #Assign FIPS code (25 counties) so that each pool only contains test sample from individuals from same county
n.county <- length(fips.sc) #46
fips.sub <- fips.sc[1:n.county]
fips.ind <- sort(rep(fips.sub,length.out=N))

#Spatial random effects (from CAR model)
#Create W matrix and D matrix
if(!exists("shp")){
  shp <- readOGR(paste0(path1,"Data/cb_2018_us_county_500k.shp")[1]) #Need to also have .dbf and .shx in the folder w/ .shp
  #Subset shp with only the fips code appearing in fips.sub
  shp <- shp[shp$GEOID%in%fips.sub,]
}

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
a.sigma <- 1
b.sigma <- 1
a.l <- 1
b.l <- 1



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


#Sample spatial locations  
shp.st <- shp
shp.sf <- st_as_sf(shp.st)
ell <- 0.75    #0.2
ell.true <- ell
ell.sigma <- 0.1
ell.init <- 0.5

sigma.sq_SE <- 0.2
sigma.sq.true <- sigma.sq_SE
sigma.sq.init <- 0.5

loc <- st_sample(shp.sf,size = N)  
loc.coordinates <- st_coordinates(loc)
id <- sample(c(1:N),m)
knot <- loc[id]

dist.c <- st_distance(loc,knot)^2
c.mat <- sigma.sq_SE*exp(-dist.c/(2*ell)) #46 by 100
dist.C <- st_distance(knot,knot)^2
C.prime <- exp(-dist.C/(2*ell)) 
C.mat <- sigma.sq_SE*C.prime  + diag(0.00001,m)#with a nugget


#xi <- mvrnorm(1,rep(0,m),C.mat)    #m dim

# a.xi <- rnorm(1,0,1)
# b.xi <- rnorm(1,0,1)
# c.xi <- rnorm(1,0,1)
# d.xi <- rnorm(1,0,1)

#Trig function 1
# a.xi <- 0.86
# b.xi <- 0.58
# c.xi <- 0.22
# d.xi <- -0.53

#Trig function 2
a.xi <- 0.396
b.xi <- -0.526
c.xi <- -0.436
d.xi <- 1.264


prev.func <- function(a.xi,b.xi,c.xi,d.xi,lat,long){
  return(a.xi*sin(lat/b.xi) + c.xi*cos(long/d.xi))
}

if(gen.from=="trig"){
  #(Option 1) Generate data from trig function
  xi <- prev.func(a.xi,b.xi,c.xi,d.xi,loc.coordinates[,2],loc.coordinates[,1]) 
  xi.true <- xi + b0
  p <- as.vector(exp(X%*%beta + xi.true)/(1+exp(X%*%beta + xi.true)))     #from trig function
  p.true <- p
  mean(p)
}else if(gen.from=="GPP"){
  #(Option 2) Generate data from GPP
  xi <- mvrnorm(1,rep(0,m),C.mat)    #m dim
  xi.mean <- mean(xi)
  xi <- xi - xi.mean
  xi.true <- c.mat%*%solve(C.mat)%*%xi
  p <- as.vector(exp(X%*%beta + xi.true)/(1+exp(X%*%beta + xi.true)))
  p.true <- p
  mean(p)
}



#hist(as.vector(C.mat)/sigma.sq_SE)


#Make prevalence maps with grid 

# define dimensions, create grid
n.grid = 3000
cellsize <- 0.05

grid <- shp.sf %>% 
  st_make_grid(cellsize = cellsize, what = "centers", n = c(n.grid/2, n.grid/2), square = F) %>% # grid of points
  st_intersection(shp.sf)   

length(grid)

grid.sf <- st_as_sf(grid)

lat <- st_coordinates(grid.sf)[,"Y"]
long <- st_coordinates(grid.sf)[,"X"]


dist.c.grid <- st_distance(grid.sf,knot)^2
c.mat.grid <- sigma.sq.init*exp(-dist.c.grid/(2*ell)) #46 by 100


#True prevalence surface
if(gen.from=="trig"){
  #(Option 1) From trig function  
  xi.grid.true <- prev.func(a.xi,b.xi,c.xi,d.xi,lat,long)
  p.grid.true <- exp(b0 + xi.grid.true)/(1+exp(b0 + xi.grid.true))
}else if(gen.from=="GPP"){
  #(Option 2) From GPP
  xi.grid.true <- c.mat.grid%*%solve(C.mat)%*%xi
  p.grid.true <- exp(beta[1] + xi.grid.true)/(1+exp(beta[1] + xi.grid.true)) 
}


max.lat = max(coordinates(shp.st)[,2])
min.lat = min(coordinates(shp.st)[,2])

max.long = max(coordinates(shp.st)[,1])
min.long = min(coordinates(shp.st)[,1])

plot.prev.true <-  ggplot() +
  geom_point(data=grid.sf, aes(x=long,y=lat,colour = p.grid.true), shape=15,size=2) +
  scale_colour_gradient(low = "white", high="darkred",limits = c(0, max(p.grid.true))) +
  coord_sf(xlim = c(min.long-1, max.long+1), ylim = c(min.lat-1, max.lat+1), expand = FALSE) +
  geom_sf(data = shp.sf, fill = NA, color = "black", lwd = 0.3) +
  theme_bw()

plot.prev.true

mean(p.grid.true)

mean(p)


#Specify the true status Y.tilde
Y.tilde <- rbinom(N,1,p)
Y.tilde.true <- Y.tilde


#List indicating which individual tested by which assay
assay.ind.list <- sens.ind.list <- spec.ind.list <- rep(NA, n.assay)


#Random assignment of assay to individuals
assay.ind.list <- sample(1:n.assay,size = N,replace=T,prob=rep(1/n.assay,n.assay))
for(i in 1:n.assay){
  sens.ind.list[assay.ind.list==i] <- df.assay$sens[i]
  spec.ind.list[assay.ind.list==i] <- df.assay$spec[i]
}

#Generate the individual test status
Y <- rep(NA,N)
for(i in 1:n.assay){
  Y[Y.tilde==1&sens.ind.list==df.assay$sens[i]] <- rbinom(length(Y.tilde[Y.tilde==1&sens.ind.list==df.assay$sens[i]]),1,df.assay$sens[i])
  Y[Y.tilde==0&spec.ind.list==df.assay$spec[i]] <- rbinom(length(Y.tilde[Y.tilde==0&spec.ind.list==df.assay$spec[i]]),1,1-df.assay$spec[i])
}


###### Array Testing (AT) ######

#Assign unique ID for each individual
df.array <- data.frame(id=1:N ,Y=Y.tilde)

Y.array <- list()
#Convert individual true status into 3-by-3 matrices
for(i in 1:n.array){
  #Y.array[[i]] <- matrix(Y.tilde[(9*(i-1)+1):(9*i)],3,3,byrow = T)
  Y.array[[i]] <- matrix(Y.tilde[(n.size^2*(i-1)+1):(n.size^2*i)],n.size,n.size,byrow = T)
}

id.array <- list()
for(i in 1:n.array){
  #id.array[[i]] <- matrix(df.array$id[(9*(i-1)+1):(9*i)],3,3,byrow = T)
  id.array[[i]] <- matrix(df.array$id[(n.size^2*(i-1)+1):(n.size^2*i)],n.size,n.size,byrow = T)
}


#Create a list to save the indices
#Create a list for each array to record individual's id (i.e. array id to individual id's)    
array.list <- list()
byrow.i <- bycol.i <- c()

for(i in 1:n.array){
  byrow.i <- c(byrow.i,((n.size*2*(i-1)+1):(n.size*2*(i-1)+n.size))) #the indices of array by row
  bycol.i <- c(bycol.i,((n.size*2*(i-1)+n.size+1):(n.size*2*(i-1)+n.size*2))) ##the indices of array by column
}

for(i in 1:((nrow(id.array[[1]])+ncol(id.array[[1]]))*n.array)){
  if(i%in%byrow.i){
    array.list[[i]] <- id.array[[trunc((i/(n.size*2))+1)]][i%%(n.size*2),] #fill in id by row
  }
  else{
    array.list[[i]] <- id.array[[ceiling(i/(n.size*2))]][,ifelse(i%%(n.size*2)==0,n.size,i%%(n.size*2)-n.size)] #fill in id by column
  }
}



#Create a list for each individual to record the pool being contributed to (i.e. individual id to pool id's)      
ind.list <- list()
for(i in 1:N){
  ind.list[[i]] <- which(lapply(lapply(lapply(array.list,"%in%",i),as.numeric),sum)==1)
}




#List indicating which pool tested by which assay
assay.pool.list <- sens.pool.list <- spec.pool.list <- rep(NA,n.pool)
assay.pool.list <- sample(1:n.assay,size = n.pool,replace=T,prob=rep(1/n.assay,n.assay))

for(i in 1:n.assay){
  sens.pool.list[assay.pool.list==i] <- df.assay$sens[i] 
  spec.pool.list[assay.pool.list==i] <- df.assay$spec[i] 
}



#True status for each pool
Z.tilde <- c()
Z <- c()
for(i in 1:n.pool){
  Z.tilde[i] <-  ifelse(sum(Y.tilde[unlist(array.list[i])])>0,1,0)
}


#Generate the test status
Z <- rep(NA,n.pool)
for(i in 1:n.assay){
  Z[assay.pool.list==i] <- rbinom(sum(assay.pool.list==i),1,df.assay$sens[i]*Z.tilde[assay.pool.list==i]+(1-df.assay$spec[i])*(1-Z.tilde[assay.pool.list==i]))
}


Z_0 <- Z
length(Z)
sum(Z)



#Retesting (append Z with Y's in the intersection of positive pools)
ind.intersection <- as.numeric(names(table(unlist(array.list[which(Z==1)]))[table(unlist(array.list[which(Z==1)]))>1]))
Z <- c(Z,Y[ind.intersection])
Z.tilde <- c(Z.tilde,Y.tilde.true[ind.intersection])

#update list
assay.pool.list <- append(assay.pool.list,assay.ind.list[ind.intersection])


for(i in 1:n.assay){
  sens.pool.list[assay.pool.list==i] <- df.assay$sens[i] 
  spec.pool.list[assay.pool.list==i] <- df.assay$spec[i] 
}

array.list <- append(array.list,ind.intersection)

write.table(Z_0,"Z_0.csv",row.names = F,col.names = F)  
write.table(Z,"Z.csv",row.names = F,col.names = F)  

ind.list <- list()
for(i in 1:N){
  ind.list[[i]] <- which(lapply(lapply(lapply(array.list,"%in%",i),as.numeric),sum)==1)
}

n.pool <- length(array.list)
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


# for(i in 1:n.assay){
#   assign(paste0("sens.init.",df.assay$id[i]),0.90)
#   assign(paste0("spec.init",df.assay$id[i]),0.90)
# }


sens.init <- 0.90
spec.init <- 0.90

c.mat <- sigma.sq.init*exp(-dist.c/(2*ell.init)) #46 by 100
C.prime <- exp(-dist.C/(2*ell.init)) 
C.mat <- sigma.sq.init*C.prime + diag(0.00001,m) #with a nugget


sej.mcmc <- spj.mcmc <- matrix(NA,n.iter,n.pool)
sej.mcmc[1,] <- sej <- rep(sens.init,n.pool) #pool specific sensitivity P(+test|+true status)
spj.mcmc[1,] <- spj <- rep(spec.init,n.pool) #pool specific specificity P(-test|-true status)

sej <- rep(sens.init,n.pool)
spj <- rep(spec.init,n.pool)


if(gen.from=="trig"){
  xi <- mvrnorm(1,rep(0,m),C.mat)
} else if(gen.from=="GPP"){
  xi <- mvrnorm(1,rep(0,m),C.mat)
  xi <- xi - mean(xi)
}


xi.N <- c.mat%*%solve(C.mat)%*%xi
xi.mcmc  <- matrix(NA,n.iter,N)
xi.mcmc[1,] <- xi.N


#p.mcmc <- matrix(NA,n.iter,N)
#p.mcmc[1,] <-  as.vector(exp(X%*%b[1,] + xi.N)/(1+exp(X%*%b[1,] + xi.N)))

xi.grid.mcmc <- matrix(NA,n.iter, length(grid))
xi.grid.mcmc[1,] <- c.mat.grid%*%solve(C.mat)%*%xi

p.grid.mcmc <- matrix(NA,n.iter, length(grid))
p.grid.mcmc[1,] <- as.vector(exp(b[1,1] + xi.grid.mcmc[1,])/(1+exp(b[1,1] + xi.grid.mcmc[1,])))


sigma.sq.mcmc <- rep(NA,n.iter)
sigma.sq.mcmc[1] <- sigma.sq <- sigma.sq.init

ell.mcmc <- rep(NA,n.iter)
ell.mcmc[1] <- ell.init



######## Use true sej and spj for debugging only (fixing sej, spj at true values) #######

# sej = sens.pool.list
# spj = spec.pool.list

#########################################################################################    

n.size
N
mean(p)

#Y.tilde <- Y.tilde.true  #remove if update Y.tilde

ell.acc.count <- 0


run.start <- Sys.time()
  for (t in 2:n.iter){
    #Step 1 Sample Y.tilde (Array testing) 
    C.inv <- chol2inv(chol(C.mat))
    cC.inv <- c.mat%*%C.inv  
    cC.inv.xi <- (cC.inv%*%xi) 
    #p1
    p1 <- p0 <- c()
    for (i in 1:N){
      ind.j <- ind.list[[i]]
      eta <- 1/(1+exp(-(X[i,]%*%b[t-1,]+(cC.inv.xi)[i])))
      #p1
      p1[i] <- as.numeric(eta)*prod(sej[ind.j]^Z[ind.j]*(1-sej[ind.j])^(1-Z[ind.j]))
      #p0
      sij <- colSums(sapply(array.list[ind.j][1:2],function(i,a) return(a[i]),Y.tilde))-Y.tilde[i]
      if(length(array.list[ind.j])==3){
        sij <- c(sij,0)
      }
      p0[i] <- (1-as.numeric(eta))*
        prod((sej[ind.j]^Z[ind.j]*(1-sej[ind.j])^(1-Z[ind.j]))^(as.numeric(sij>0))*
               ((1-spj[ind.j])^Z[ind.j]*spj[ind.j]^(1-Z[ind.j]))^(as.numeric(!sij>0)))
      
      Y.tilde[i] <- rbinom(1,1,prob = p1[i]/(p1[i]+p0[i]))
    }
    
    Y.tilde.mcmc[t,] <- Y.tilde
  
    
    #Step 2 Compute Z.tilde based on Y.tilde
    Z.tilde <- c()
    for(i in 1:n.pool){
      Z.tilde[i] <-  ifelse(sum(Y.tilde[unlist(array.list[i])])>0,1,0)
    }
    
    
    #Step 3 Sample sens and spec
    #Sensitivity & Specificity
    for(i in 1:n.assay){
      aem <- a.sens + sum(Z[assay.pool.list==i]*Z.tilde[assay.pool.list==i])
      bem <- b.sens + sum((1-Z[assay.pool.list==i])*Z.tilde[assay.pool.list==i])
      
      apm <- a.spec + sum((1-Z[assay.pool.list==i])*(1-Z.tilde[assay.pool.list==i]))
      bpm <- b.spec + sum(Z[assay.pool.list==i]*(1-Z.tilde[assay.pool.list==i]))
      
      sej[assay.pool.list==i] <- rep(rbeta(1,aem,bem),sum(assay.pool.list==i))
      spj[assay.pool.list==i] <- rep(rbeta(1,apm,bpm),sum(assay.pool.list==i))
    }
    
    sej.mcmc[t,] <- sej
    spj.mcmc[t,] <- spj
    
    
    #Step 4 Sample beta
    omg <- diag(w[t-1,])
    M <- t(X)%*%omg%*%X+b.cov.inv
    k <- Y.tilde - 1/2
    h <- k/diag(omg)
    Q <- t(X)%*%omg%*%(h-cC.inv.xi) + b.cov.inv%*%b.m
    M.inv <- chol2inv(chol(M))   
    b[t,] <- as.vector(mvrnorm(n = 1,M.inv%*%Q,M.inv))
    #b[t,] <- beta #remove if update b
    
    
    ##Step 5 Sample w
    w[t,] <- apply(X = X%*%b[t,] + cC.inv%*%xi,1,FUN = rpg,num=1,h=1)
    #update omg and h
    omg <- diag(w[t,])
    k <- Y.tilde - 1/2
    h <- k/diag(omg)
    
    
    #Step 6 Sample xi
    omg.cC.inv <- omg%*%cC.inv
    V <- t(cC.inv)%*%omg.cC.inv + C.inv
    d <- t(t(h-X%*%b[t,])%*%omg.cC.inv)
    #d <- as.vector((t(h)-t(b[t,])%*%t(X))%*%omg%*%cC.inv) #equivalent to above
    #V.inv <- solve(V)
    V.inv <- chol2inv(chol(V))
    
    
    if(gen.from=="trig"){
      xi <- mvrnorm(n = 1,mu= V.inv%*%d,Sigma = V.inv)
    } else if(gen.from=="GPP"){
      xi <- mvrnorm(n = 1,mu= V.inv%*%d,Sigma = V.inv)
      xi <- xi - mean(xi)
    } 
    
    cC.inv.xi <- cC.inv%*%xi
    xi.mcmc[t,] <- cC.inv.xi
    
    c.mat.grid <- sigma.sq*exp(-dist.c.grid/(2*ell)) #46 by 100
    # dist.C <- st_distance(knot,knot)^2
    # C.prime <- exp(-dist.C/(2*ell)) 
    # C.mat <- sigma.sq_SE*C.prime  + diag(0.00001,m)#with a nugget
    
    
    xi.grid <- c.mat.grid%*%C.inv%*%xi
    xi.grid.mcmc[t,] <- xi.grid
    
    
    if(gen.from=="trig"){
      eta.grid <- xi.grid
      exp.eta.grid <- exp(eta.grid)
    }else if(gen.from=="GPP"){
      eta.grid <- b[t,1] + xi.grid
      exp.eta.grid <- exp(eta.grid)
    }  
    
    #p.mcmc[t,] <- as.vector(exp(X%*%b[t,] + cC.inv.xi)/(1+exp(X%*%b[t,] + cC.inv.xi)))
    p.grid.mcmc[t,] <- as.vector(exp.eta.grid/(1+exp.eta.grid))
    
    
    #Step 7 Sample sigma.sq
    #C.prime.inv <- solve(C.prime)
    C.prime.inv <- chol2inv(chol(C.prime + diag(0.00001,m)))
    sigma.sq <- rinvgamma(1,
                          shape = a.sigma + m/2,
                          rate = b.sigma + 1/2*t(xi)%*%C.prime.inv%*%xi)
    sigma.sq.mcmc[t] <- sigma.sq
    
    
    
    #Step 8 Sample ell (Metropolis-Hastings Step)
    u <- runif(1,0,1)
    ell.log <- rnorm(1,log(ell),ell.sigma)
    ell.new <- exp(ell.log)
    
    C.mat.current <- sigma.sq*(exp(-dist.C/(2*ell.mcmc[t-1]))) + diag(0.00001,m)
    C.mat.new <- sigma.sq*(exp(-dist.C/(2*ell.new))) + diag(0.00001,m)
    
    eta <- X%*%b[t,]+ xi.mcmc[t,]
    
    
    log.p.log.ell.current <- -1/2*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) - 1/2*determinant(C.mat.current,logarithm=T)$modulus[1]-
      1/2*t(xi)%*%chol2inv(chol(C.mat.current))%*%xi -
      a.l*log(ell.mcmc[t-1]) - b.l/ell.mcmc[t-1] 
    log.p.log.ell.new <- -1/2*(t(h)%*%omg%*%h - 2*t(h)%*%omg%*%eta + t(eta)%*%omg%*%eta) -1/2*determinant(C.mat.new,logarithm=T)$modulus[1]-1/2*t(xi)%*%chol2inv(chol(C.mat.new))%*%xi - a.l*ell.log - b.l/ell.new
    
    acc <- min(1,exp(log.p.log.ell.new - log.p.log.ell.current))
    if(u < acc){
      ell.mcmc[t] <- ell <- ell.new
      C.mat.current <- C.mat.new
      c.mat <- sigma.sq*exp(-dist.c/(2*ell.new)) 
      C.prime <- exp(-dist.C/(2*ell.new)) 
      C.mat <- sigma.sq*C.prime + diag(0.00001,m)  #with a nugget
      
      ell.acc.count <- ell.acc.count + 1
      
    } else {
      ell.mcmc[t] <- ell <- ell.mcmc[t-1]
      c.mat <- sigma.sq*exp(-dist.c/(2*ell.mcmc[t-1])) 
      C.prime <- exp(-dist.C/(2*ell.mcmc[t-1])) 
      C.mat <- sigma.sq*C.prime + diag(0.00001,m) #with a nugget
      
    }
    
    #print(t)
    print(paste0(t,", ell.acc.count = ", ell.acc.count))
  } 
    
    
run.end <- Sys.time()
run.time <- difftime(run.end, run.start, units='mins')
run.time 

    
    

#Trace-plots (full samples)
#Beta
for(i in 1:length(beta)){
  plot(b[,i],type = "l", main = paste0("Beta_",i))
}

#Sensitivity

assay.pool.list

sem.mcmc <- matrix(NA,n.iter,n.assay)
spm.mcmc <- matrix(NA,n.iter,n.assay)   
    
    
for(i in 1:n.assay){
  sem.mcmc[,i] <- sej.mcmc[,which(assay.pool.list==i)[1]] 
  spm.mcmc[,i] <- spj.mcmc[,which(assay.pool.list==i)[1]] 
  plot(sem.mcmc[,i],type = "l",main = paste0("sem_",i))
  plot(spm.mcmc[,i],type = "l",main = paste0("spm_",i))
}


#xi
for(i in 1:length(xi)){
  plot(xi.mcmc[,i],type = "l", main = paste0("xi_",i))
}

#ell
plot(ell.mcmc,main="ell.sq",type = "l") 

#sigma.sq
plot(sigma.sq.mcmc,main="sigma.sq.sq",type = "l") 


#Remove burn-in
#Remove burn-in and compute geweke statistics
b.pred <- b[-(1:(0.5*n.iter)),]
b.geweke <- rep(NA,ncol(b.pred))
for(i in 1:ncol(b.pred)){
  b.geweke[i] <- unlist(geweke.diag(b.pred[,i]))["z.var1"]
}



Y.tilde.pred <- Y.tilde.mcmc[-(1:(0.5*n.iter)),]

sem.pred <- as.matrix(sem.mcmc[-(1:(0.5*n.iter)),])
spm.pred <- as.matrix(spm.mcmc[-(1:(0.5*n.iter)),])

sem.geweke <- rep(NA,n.assay)
spm.geweke <- rep(NA,n.assay)
for(i in 1:n.assay){
  sem.geweke[i] <- unlist(geweke.diag(sem.pred[,i]))["z.var1"]
  spm.geweke[i] <- unlist(geweke.diag(spm.pred[,i]))["z.var1"]
}



xi.pred <- xi.mcmc[-(1:(0.5*n.iter)),]
xi.geweke <- rep(NA,ncol(xi.pred))
for(i in 1:ncol(xi.pred)){
  xi.geweke[i] <- unlist(geweke.diag(xi.pred[,i]))["z.var1"]
}


xi.grid.pred <- xi.grid.mcmc[-(1:(0.5*n.iter)),]
# xi.grid.geweke <- rep(NA,ncol(xi.grid.pred))
# for(i in 1:ncol(xi.grid.pred)){
#   xi.grid.geweke[i] <- unlist(geweke.diag(xi.grid.pred[,i]))["z.var1"]
# }


# #p
# p.pred <- p.mcmc[-(1:(0.5*n.iter)),]
# p.geweke <- rep(NA,N)
# for(i in 1:N){
#   p.geweke[i] <- unlist(geweke.diag(p.pred[,i]))["z.var1"]
# }


#p.grid
p.grid.pred <- p.grid.mcmc[-(1:(0.5*n.iter)),]
# p.grid.geweke <- rep(NA,ncol(p.grid.pred))
# for(i in 1:ncol(xi.grid.pred)){
#   p.grid.geweke[i] <- unlist(geweke.diag(p.grid.pred[,i]))["z.var1"]
# }


ell.pred <- ell.mcmc[-(1:(0.5*n.iter))]
ell.geweke <- unlist(geweke.diag(ell.pred))["z.var1"]

sigma.sq.pred <- sigma.sq.mcmc[-(1:(0.5*n.iter))]
sigma.sq.geweke <- unlist(geweke.diag(sigma.sq.pred))["z.var1"]  

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
sem.hat <- apply(sem.pred,2,mean)
sem.upp <- apply(sem.pred,2,quantile,probs=0.975)
sem.low <- apply(sem.pred,2,quantile,probs=0.025)
sem.sd <- apply(sem.pred,2,sd)

sem.diff <- df.assay$sens - sem.hat
#sem.cover <- ifelse(df.assay$sens[unique(assay.pool.list)]>=sem.low&df.assay$sens[unique(assay.pool.list)]<=sem.upp,1,0)
sem.cover <- ifelse(df.assay$sens>=sem.low&df.assay$sens<=sem.upp,1,0)
sem.out <- rbind(sem.diff,sem.cover,sem.hat,sem.sd)
#sem.out


#Specificity MCMC samples
spm.hat <- apply(spm.pred,2,mean)
spm.upp <- apply(spm.pred,2,quantile,probs=0.975)
spm.low <- apply(spm.pred,2,quantile,probs=0.025)
spm.sd <- apply(spm.pred,2,sd)

spm.diff <- df.assay$spec - spm.hat
#spm.cover <- ifelse(df.assay$spec[unique(assay.pool.list)]>=spm.low&df.assay$spec[unique(assay.pool.list)]<=spm.upp,1,0)
spm.cover <- ifelse(df.assay$spec>=spm.low&df.assay$spec<=spm.upp,1,0)
spm.out <- c(spm.diff,spm.cover,spm.hat,spm.sd)
#spm.out



#ell MCMC samples
ell.hat <- mean(ell.pred)
ell.upp <- quantile(ell.pred,0.975)
ell.low <- quantile(ell.pred,0.025)
ell.sd <- sd(ell.pred)

ell.diff <- ell.true - ell.hat
ell.cover <- ifelse(ell.true>=ell.low&ell.true<=ell.upp,1,0)
ell.out <- c(ell.diff,ell.cover,ell.hat,ell.sd)
ell.out


#sigma.sq MCMC samples
sigma.sq.hat <- mean(sigma.sq.pred)
sigma.sq.upp <- quantile(sigma.sq.pred,0.975)
sigma.sq.low <- quantile(sigma.sq.pred,0.025)
sigma.sq.sd <- sd(sigma.sq.pred)

sigma.sq.diff <- sigma.sq.true - sigma.sq.hat
sigma.sq.cover <- ifelse(sigma.sq.true>=sigma.sq.low&sigma.sq.true<=sigma.sq.upp,1,0)
sigma.sq.out <- c(sigma.sq.diff,sigma.sq.cover,sigma.sq.hat,sigma.sq.sd)
sigma.sq.out



#xi MCMC samples
dist.c <- st_distance(loc,knot)^2
c.mat <- sigma.sq.hat*exp(-dist.c/(2*ell.hat)) #46 by 100
dist.C <- st_distance(knot,knot)^2
C.prime <- exp(-dist.C/(2*ell.hat)) 
C.mat <- sigma.sq.hat*C.prime  + diag(0.00001,m)#with a nugget


xi.hat <- as.vector(apply(xi.pred,2,mean))
xi.upp <- as.vector(apply(xi.pred,2,quantile,probs=0.975))
xi.low <- as.vector(apply(xi.pred,2,quantile,probs=0.025))
xi.sd <- as.vector(apply(xi.pred,2,sd))

xi.diff <- as.vector(xi.true - xi.hat)
xi.cover <- as.vector(ifelse(xi.true>=xi.low&xi.true<=xi.upp,1,0))
xi.out <- rbind(xi.diff,xi.cover,xi.hat, xi.sd)
xi.out

plot(xi.true - colMeans(xi.pred),main="Residual plot_xi.pred")

xi.grid.hat <- as.vector(apply(xi.grid.pred,2,mean))
xi.grid.upp <- as.vector(apply(xi.grid.pred,2,quantile,probs=0.975))
xi.grid.low <- as.vector(apply(xi.grid.pred,2,quantile,probs=0.025))
xi.grid.sd <- as.vector(apply(xi.grid.pred,2,sd))

xi.grid.diff <- as.vector(xi.grid.true - xi.grid.hat)
xi.grid.cover <- as.vector(ifelse(xi.grid.true>=xi.grid.low&xi.grid.true<=xi.grid.upp,1,0))
xi.grid.out <- rbind(xi.grid.diff,xi.grid.cover,xi.grid.hat, xi.grid.sd)
xi.grid.out


# #p
# p.hat <-  as.vector(apply(p.pred,2,mean))
# p.upp <- as.vector(apply(p.pred,2,quantile,probs=0.975))
# p.low <- as.vector(apply(p.pred,2,quantile,probs=0.025))
# p.sd <- as.vector(apply(p.pred,2,sd))
# 
# p.diff <- as.vector(p.true - p.hat)
# p.cover <- as.vector(ifelse(p.true>=p.low&p.true<=p.upp,1,0))
# p.out <- rbind(p.diff,p.cover,p.hat, p.sd)
# p.out




#p.grid
p.grid.hat <-  as.vector(apply(p.grid.pred,2,mean))
p.grid.upp <- as.vector(apply(p.grid.pred,2,quantile,probs=0.975))
p.grid.low <- as.vector(apply(p.grid.pred,2,quantile,probs=0.025))
p.grid.sd <- as.vector(apply(p.grid.pred,2,sd))

p.grid.diff <- as.vector(p.grid.true - p.grid.hat)
p.grid.cover <- as.vector(ifelse(p.grid.true>=p.grid.low&p.grid.true<=p.grid.upp,1,0))
p.grid.out <- rbind(p.grid.diff,p.grid.cover,p.grid.hat, p.grid.sd)
p.grid.out  




write.table(b.out,"b.out.txt",row.name=F,col.name=F)         
write.table(Y.tilde.out,"Y.tilde.out.txt",row.name=F,col.name=F)         
write.table(sem.out,"sem.out.txt",row.name=F,col.name=F)
write.table(spm.out,"spm.out.txt",row.name=F,col.name=F)
write.table(xi.out,"xi.out.txt",row.name=F,col.name=F) 
write.table(xi.grid.out,"xi.grid.out.txt",row.name=F,col.name=F) 
write.table(p.grid.out,"p.grid.out.txt",row.name=F,col.name=F) 
write.table(ell.out,"ell.out.txt",row.name=F,col.name = F)
write.table(sigma.sq.out,"sigma.sq.out.txt",row.name=F,col.name = F)

write.table(lat,"lat.txt",row.name=F,col.name=F)
write.table(long,"long.txt",row.name=F,col.name=F)
write.table(assay.pool.list,"assay.pool.list.txt",row.name=F,col.name=F)

write.table(b.geweke,"b.geweke.txt",row.name=F,col.name=F)
write.table(sem.geweke,"sem.geweke.txt",row.name=F,col.name=F)
write.table(spm.geweke,"spm.geweke.txt",row.name=F,col.name=F)
write.table(xi.geweke,"xi.geweke.txt",row.name=F,col.name=F) 
#write.table(xi.grid.geweke,"xi.grid.geweke.txt",row.name=F,col.name=F) 
#write.table(p.grid.geweke,"p.grid.geweke.txt",row.name=F,col.name=F) 
write.table(ell.geweke,"ell.geweke.txt",row.name=F,col.name = F)
write.table(sigma.sq.geweke,"sigma.sq.geweke.txt",row.name=F,col.name=F) 
# write.table(p.out,"p.out.txt",row.name=F,col.name=F) 
# write.table(p.geweke,"p.geweke.txt",row.name=F,col.name=F) 







