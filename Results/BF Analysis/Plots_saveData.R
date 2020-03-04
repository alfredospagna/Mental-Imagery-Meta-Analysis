############################################
##  Code to produce plots as in the paper ##
############################################
rm(list = ls())
# Set seed 
set.seed(50)
# Set working directory 
setwd("C:/Users/jlu6/Desktop/code_vmi_50000_24_rest") #AS_change this
require(oro.nifti)
norm_vec <- function(x) sqrt(sum(x^2))  

#########################
##  Read in the data   ##
#########################
data <- read.table('foci.txt', header = F, sep = "\t")  #move to folder with foci
data <- data[, c(2:4, 1)]
names(data) <- c('X','Y', 'Z', 'study')          # AS_add header Coordinates and study number
attach(data)
data[1:10, ]

N <- length(unique(data[, 4]))
ni <- tabulate(data$study)    		  # Number of foci per study
exclude <- which(ni == 0)           # Contrasts to exclude because no foci are left
ni <- ni[which((tabulate(data$study) != 0) == TRUE)]
ID <- rep(1:N, ni)
d <- dim(data)[2] - 1  				      # Number of dimensions
data <- cbind(data[, c(1:d)], ID)

##################
##  Create axes ##
##################
mask_full = readNIfTI('brainmask.nii')  # Finer grid used for plotting: 2 x 2 x 2 mm brain mask (voxels 91 x 109 x 91)
# attributes(mask_full)
# mask_centre <- c(90, -126, -72)
Ix <- seq(1, 182, by = 2)        
Iy <- seq(1, 218, by = 2)      
Iz <- seq(1, 182, by = 2)       
A <- (Iy[2] - Iy[1])*(Ix[2] - Ix[1])*(Iz[2] - Iz[1])       # Volume of each voxel

# Axes to reproduce plots similar to those in the paper (make it a squared figure)
Ixplotting <- seq(-109, 107, by = 2)
Iyplotting <- seq(-127, 89, by = 2)
Izplotting <- seq(-91, 125, by = 2)

Ixresc <- seq(-91, 89, by = 2)
Iyresc <- seq(-127, 89, by = 2)
Izresc <- seq(-73, 107, by = 2)



########################################################
##  Define radial cubic B-splines or Gaussian kernels ##
########################################################
# ND: mask is 91 x 109 x 91, but mask is null above z = 82 
# Discard slices above z = 79 provided there are no foci falling outside the mask				
Grid <- c()  	
grid_slice <- as.matrix(expand.grid(Ix,Iy))
for (n in 1:91){	# 1:91
  map = mask_full[,,n]
  map = map > 0.5
  msk = as.vector(map)
  Grid <- rbind(Grid, cbind(grid_slice[msk, ], Iz[n]))
}
dim(Grid)
rm(grid_slice)

##################################################
##   Restrict computation to points in the mask ##
##    to avoid bias from off-mask estimation    ##
##################################################
xx <- seq(40, 145, length.out = 8)
yy <- seq(40, 180, length.out = 8)
zz <- seq(38, 90, length.out = 7)
kernels <- expand.grid(xx, yy, zz)
knots <- kernels

# At every odd numbered z slice, shift the x coordinate to re-create a chess-like kernel grid
for (n in 1:floor(length(zz)/2)){
  kernels[which(kernels[, 3] == zz[n*2]), 1] <- kernels[which(kernels[, 3] == zz[2*n]), 1] + 10
}

knots <- kernels[- which(kernels[, 1] > 145), ] 

nu <- 1/(2*256)
d <- 3
B.pred <- matrix(0, dim(Grid)[1], dim(knots)[1]) # Observed matrix of basis functions
for (i in 1:dim(B.pred)[1]){
  obs.knot <- matrix(as.numeric(Grid[i, 1:d]), dim(knots)[1], d, byrow = TRUE) - knots
  B.pred[i, ] <- as.matrix(exp(-nu * apply(obs.knot, 1, norm_vec)^2))
}
B.pred <- cbind(1, B.pred)    				 # Dense basis for integral evaluation

for (i in 1:dim(B.pred)[2]){
  ind = which(B.pred[, i] < 1e-35)
  if (length(ind) > 0) B.pred[ind, i] = 0
}
p = dim(B.pred)[2]


#######################################
##      Define global constants      ##
#######################################																													 
nrun = 50000;                       
burn = 25000;                       
thin = 25;                       
every = 100;                      
start = 250;		
sp = (nrun - burn)/thin;	
maxk = 50


######################################
##      Reading covariates in       ##
######################################
Z = read.table('covariates.txt');
r = dim(Z)[2];          							    # Number of covariates per study

######################################
##      Reading study-type in       ##
######################################
outcome <- read.table('studytype.txt'); Y = outcome
								 

train <- as.matrix(as.numeric(read.table(file = "train.txt", header = FALSE, sep = ",")))
test <- setdiff(1:N, train)


##################################################
##    Code for 3D plotting of the intensities   ##
##  Need to extract the basis evaluated at the  ##
##               matching z-slice               ##
##################################################
library(fields)
theta <- as.matrix(read.table(file = "thetacbma_post.txt", header = FALSE, sep = " "))

library(oro.nifti)

####################################################################
# Visual mental imagery
####################################################################
type1 <- which(Y == 1);     #AS_extract all the values from study type = 1 (VM)
type1_cbma <- matrix(0, sp, p)
for (i in 1:length(type1)){
  j = type1[i]
  ind <- c()
  for (l in 1:p){
    k = l - 1
    ind[l] <- k*N + j
  }
  type1_cbma <- type1_cbma + theta[, ind] 
}
type1_cbma <- type1_cbma/length(type1)
intensity1_cbma <- apply((B.pred %*% t(type1_cbma)), 1, mean)

img_mat_3d <- array(NA,dim=c(91,109,91))

for (n in 1:91){
  ind = which(Grid[, 3] == Iz[n]) 
  if (length(ind) > 0) { # This checks the mask is non-null 
					 
    # This is the z coordinate in mms of what I'm plotting (the center of the voxel)
    map = mask_full[,,n]
    map = map > 0.5
    msk = as.vector(map)
    mask.pos <- which(map == TRUE, arr.ind = TRUE) 
    
    int1 <- intensity1_cbma[ind]
    image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))  #AS_log intensity (not posterior prob)
      
    for (f in 1:dim(mask.pos)[1]){
      pos <- mask.pos[f, ]
      image_matrix[pos[1], pos[2]] <- (int1[f]) 
    } 
    
    img_mat_3d[,,n] <- image_matrix
  } 	   
}

writeNIfTI(img_mat_3d,"img_50000_25_vmi_rest")  #AS_Change

h0 = exp(img_mat_3d)                      #AS_exponential of log of posterior internsity of h0

h1 = mean(exp(img_mat_3d),na.rm = TRUE)   #AS_exponential of log of posterior internsity of h1

bf = h0/h1
writeNIfTI(bf,"bf_50000_25_vmi_rest")         #AS_Change

#Back to Matlab and transform back from 91 x 109 x 91 into MNI space


## Choose some axial slices ##
set = c( 20, 24, 28, 32, 38, 45)    #AS_change here for more slices

for (i in 1:length(set)){
  n = set[i]
  # This is the z coordinate in mms of what I'm plotting (the center of the voxel)
  map = mask_full[,,n]
  map = map > 0.5
  msk = as.vector(map)
  mask.pos <- which(map == TRUE, arr.ind = TRUE) 
  int1 <- intensity1_cbma[which(Grid[, 3] == Iz[n])]
  image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
  for (f in 1:dim(mask.pos)[1]){
    pos <- mask.pos[f, ]
    image_matrix[pos[1], pos[2]] <- (int1[f]) 
  }
  image_matrix1 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
  scale1 <- range((image_matrix1), na.rm = TRUE)
  
  out_cbma <- c(min(scale1), max(scale1))
  
  pdf(file = paste("Z = ", Izresc[n], ".pdf"), width = 14, height = 4.5)
  m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.20, 0.90), c(0, 1, 0, 0.20))
  split.screen(m)
  screen(1)
  par(mar = c(0,1,1,1), oma = c(0, 0, 1, 0))
  title(main = paste("Post mean log intensities Z = ", Izresc[n]))  #AS_title of plot
  
  split.screen(c(1,5), screen=2)-> ind2
  j = 1
  screen(ind2[j])
  par(pty = "s", mar = c(0, 0, 1, 0), oma = c( 0, 0, 1, 0 ) )
  # par(mfrow = c(2,3), pty = "s", mar = c(2,4,2,3))
  image(c(Ixplotting), c(Iyplotting), matrix((image_matrix1), length(Ixplotting), length(Ixplotting), byrow = FALSE), xlab = ' ', ylab = ' ',
        zlim = out_cbma, main = "vmi", col = (gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2)), axes = F)
 
  screen(3)
  par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
  image.plot(zlim = out_cbma, horizontal = T, legend.only = TRUE, 
             legend.width = .6, col=gray.colors(64, start = 0.9, end = 0.1, gamma = 2.2))
  close.screen( all=TRUE)
  
  dev.off()
}

# # # # # # # # # # # # # # # # # 
# # -- Examine learnt bases -- # 
# # # # # # # # # # # # # # # #
factor <- as.matrix(read.table(file = "Factor.txt", header = FALSE, sep = " ", skip = burn/every))
Lambdaout <- as.matrix(read.table(file = "Lambda.txt", header = FALSE, sep = " "))

j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out = out + B.pred[ind, ] %*% lambda_i
}
out <- out / sp

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out2 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out2 = out2 + B.pred[ind, ] %*% lambda_i
}
out2 <- out2 / sp

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out3 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out3 = out3 + B.pred[ind, ] %*% lambda_i
}
out3 <- out3 / sp

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out4 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out4 = out4 + B.pred[ind, ] %*% lambda_i
}
out4 <- out4 / sp

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out5 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out5 = out5 + B.pred[ind, ] %*% lambda_i
}
out5 <- out5 / sp

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
out6 <- matrix(0, length(ind), maxk)
for(i in 1:sp){
  lambda_i = matrix(Lambdaout[i, ], nrow = p, ncol = maxk, byrow = F)
  out6 = out6 + B.pred[ind, ] %*% lambda_i
}
out6 <- out6 / sp

pdf(file = paste("learnt1.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[1])))
# close.screen(all = TRUE)

split.screen(c(1,6), screen = 2)-> ind2

l = 1
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out1range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

g = 1
screen(ind2[g])
j = 1
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 2
screen(ind2[g])
j = 2
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 3
screen(ind2[g])
j = 3
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 4
screen(ind2[g])
j = 4
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 5
screen(ind2[g])
j = 5
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

g = 6
screen(ind2[g])
j = 6
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out1range, main = paste("Z=",Izresc[set[j]]))

screen(3) 
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out1range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))
close.screen(all = TRUE)
dev.off()


pdf(file = paste("learnt2.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[2])))
# close.screen(all = TRUE)

split.screen(c(1,6), screen = 2)-> ind2
l = 2
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3) 
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))
close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt3.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.90, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[3])))

split.screen(c(1,6), screen = 2)-> ind2
l = 3
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt4.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.80, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[4])))

split.screen(c(1,6), screen = 2)-> ind2
l = 4
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()

pdf(file = paste("learnt5.pdf"), width = 14, height = 4)
m <- rbind(c(0, 1, 0.80, 1), c(0, 1, 0.2, 0.90), c(0, 1, 0.15, 0.25))
split.screen(m)
screen(1)
par(mar = c(0,1,2,1), oma = c(0, 0, 1, 0))
title(main = expression(paste(tilde(phi)[5])))

split.screen(c(1,6), screen = 2)-> ind2
l = 5
j = 1
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out[j, l]
}
image_matrix1 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale1 = range(image_matrix1, na.rm = T)

j = 2
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE)  
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out2[j, l]
}
image_matrix2 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale2 = range(image_matrix2, na.rm = T)

j = 3
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out3[j, l]
}
image_matrix3 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale3 = range(image_matrix3, na.rm = T)

j = 4
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out4[j, l]
}
image_matrix4 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale4 = range(image_matrix4, na.rm = T)

j = 5
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out5[j, l]
}
image_matrix5 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale5 = range(image_matrix5, na.rm = T)

j = 6
ind = which(Grid[, 3] == Iz[set[j]])
map = mask_full[,,set[j]]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE) 
maxk <- 50
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))
for (j in 1:dim(mask.pos)[1]){
  pos <- mask.pos[j, ]
  image_matrix[pos[1], pos[2]] <- out6[j, l]
}
image_matrix6 <- rbind(matrix(NA, 9, length(Iy)), image_matrix, matrix(NA, 9, length(Iy)))
scale6 = range(image_matrix6, na.rm = T)

out2range = c(min(min(scale1), min(scale2), min(scale3), min(scale4), min(scale5), min(scale6)),
              max(max(scale1), max(scale2), max(scale3), max(scale4), max(scale5), max(scale6)))

l = 1
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix1, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 2
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix2, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 3
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix3, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 4
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix4, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 5
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix5, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

l = 6
screen(ind2[l])
par(pty = "s", mar = c(0, 0, 0, 0), oma = c( 0, 0, 2, 0 ))
image(c(Ixplotting), c(Iyplotting), image_matrix6, xlab = ' ', ylab = ' ', axes = F,
      col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)), zlim = out2range, main = paste("Z=",Izresc[set[l]]))

screen(3)
par(mfrow = c(1,1),mar = c(2,0,0,0), oma = c(0, 0, 0, 0))
image.plot(zlim = out2range, horizontal = T, legend.only = TRUE, 
           legend.width = .6, col=gray.colors(64, start = 0, end = 1, gamma = 2.2))

close.screen(all = TRUE)
dev.off()


#####################################
##    Study-specific intensities   ## The following lines convert data to/from different spaces 
#####################################
# -- # Convert coordinates to voxel space # -- #
origin = c(90,-126,-72)
data$X = round((origin[1]-data$X)/2)
data$Y = round((data$Y-origin[2])/2)
data$Z = round((data$Z-origin[3])/2)

# foci that fall outside [1,91]x[1,109]x[1,91] are removed
keep = (data$X>=1)&(data$X<=91)&(data$Y>=1)&(data$Y<=109)&(data$Z>=1)&(data$Z<=91)
# Restrict data to the amygdalae
# keep = (data$X>=1)&(data$X<=91)&(data$Y>=1)&(data$Y<=109)&(data$Z>=20)&(data$Z<=45)
data = data[keep,]

# Eliminate foci that fall outside of the mask
mask = readNIfTI('brainmask.nii')
mask = (mask>0)*1
for (i in dim(data)[1]:1) {
  if (mask[ data[i,1] , data[i,2] , data[i,3] ]==0){
    data = data[-i,]
  }
}

# Convert data to mm space 
xaxis <- seq(1, 182, by = 2)
yaxis <- seq(1, 218, by = 2)
zaxis <- seq(1, 182, by = 2)
data$X = xaxis[data$X]
data$Y = yaxis[data$Y]
data$Z = zaxis[data$Z]

# Rescale data: NOTE: might have to change z-coordinate if plotting something other than axial slices
resdata <- matrix(NA, sum(ni), 3)
for (i in 1:sum(ni)){
  ind <- c(which(Ix == data[i, 1]), which(Iy == data[i, 2]), which(Iz == data[i, 3]))
  resdata[i, ] <- c(Ixplotting[ind[1]+9], Iyresc[ind[2]], Izresc[ind[3]])
}
data <- cbind(resdata, ID)
data <- as.data.frame(data)
names(data) <- c('x','y', 'z', 'study')          # Coordinates of the foci and study number
attach(data)

# -- # Choose a study # -- #
i <- sample(1:N, 1)
th <- matrix(0, sp, p)
ind <- c()
for (l in 1:p){
  k = l - 1
  ind[l] <- k*N + i
}
th <- theta[, ind] 

# -- # Choose a slice: this is the z coordinate in mms of what I'm plotting (the center of the voxel) # -- #
n <- 40
grid_slice <- as.matrix(expand.grid(Ix,Iy))
map = mask_full[,,n]
map = map > 0.5
msk = as.vector(map)
mask.pos <- which(map == TRUE, arr.ind = TRUE)  
mask_plot <- grid_slice[msk, ]
coords <- cbind(mask_plot, Iz[n])

B.plot <- matrix(0, dim(mask_plot)[1], dim(knots)[1]) # Observed matrix of basis functions
for (j in 1:dim(B.plot)[1]){
  obs.knot <- matrix(as.numeric(coords[j, 1:d]), dim(knots)[1], d, byrow = TRUE) - knots
  B.plot[j, ] <- as.matrix(exp(-nu * apply(obs.knot, 1, norm_vec)^2)) 
}
B.plot <- cbind(1, B.plot)
intensity <- apply(exp(B.plot %*% t(th)), 1, mean)
image_matrix <- matrix(NA, nrow = length(Ix), ncol = length(Iy))

for (f in 1:dim(mask.pos)[1]){
  pos <- mask.pos[f, ]
  image_matrix[pos[1], pos[2]] <- intensity[f] 
}
image_matrix2 <- rbind(matrix(NA, nrow = 9, ncol = length(Iy)), image_matrix, matrix(NA, nrow = 9, ncol = length(Iy)))
out <- range(image_matrix2, na.rm = TRUE)

pp <- data[ID == i, 1:d]; pp <- pp[order(pp$z), ]
image.plot(c(Ixplotting), c(Iyplotting), image_matrix2, nlevel = 64, xlab = 'x', ylab = 'y', zlim = out,
           legend.shrink = 1, legend.width = 1, col = (gray.colors(64, start = 0, end = 1, gamma = 2.2)))
if (sum(which(pp[, 3] == Izresc[n])) > 0)  points(pp[which(pp[, 3] == Izresc[n]), 1:2], pch = 16)
title(main = paste("Est int study", i, ", slice z = ", Izresc[n]+1, "mm"))

