# FMMcsVS
Bayesian Finite Mixure Regression Models with Cluster-Specific Variable Selections

## Download and install the package:
devtools::install_github("zwang400/FMMcsVS")

## load the package
library(FMMcsVS)

## generate simulation data
### independent covariates with default parameter values
sim_data_1 <- data_gen_func()

### for D=6, M=3, bloacked-wise correlation among X1-X2, X3-X4, X5-X6, correation of 0.2,0.3,0.4
sim_data_2 <- data_gen_func(rho=c(0.2, 0.3, 0.4))

### correlated covariates with pre-specified covariance matrix
sim_data_3 <- data_gen_func(cor_mtx = matrix(0.8, 6, 6))

### generate data for split model, independent covariates (W)
sim_data_split_1 <- data_gen_split()

### data for split model, pair-wise correlation of 0.5 among W
sim_data_split_2 <- data_gen_split(rho=0.5)

## run MCMC simulations

### FDMM with fixed \gamma=1, with VS
sim_fdmm_1 <- simulation_func(sim_data_1$X, sim_data_1$y, gamma_hyperprior=F, gamma_fixed=1)

### FBMM without VS, b_bel with hyperprior
sim_fbmm_1 <- simulation_func(sim_data_1$X, sim_data_1$y, prior="Bessel")

### Dynamic FDMM with VS
sim_dyn_fdmm_1 <- simulation_func(sim_data_1$X, sim_data_1$y, prior="Dynamic_FDMM")

### FUMM with S~Unif(0.2, 1)
sim_fumm_1 <- simulation_func(sim_data_1$X, sim_data_1$y, prior="Uniform", a_unif=0.2)

### split model with FBMM with VS
sim_split_fbmm_1 <- simulation_split(sim_data_split_1$W, sim_data_split_1$Z, sim_data_split_1$y, prior="Bessel")

### RPMS (DPM) with VS
sim_rpms_1 <- simulation_func_rpms(sim_data_1$X, sim_data_1$y)

## posterior inference
calculate ARI, MSE, VS errors

### FBMM 
post_inf_fbmm_1 <- post_inf(sim_res=sim_fbmm_1, data=sim_data_1)

### RPMS
post_rpms_1 <- post_inf_rpms(sim_res=sim_rpms_1, data=sim_data_1)

## real data analysis
For country, coffee, flea data, datasets been loaded automatically after loading the package, no need to use data(...)

### for coffee data
X <- coffee[, -c(1,2,10)]
y <- coffee[,10]
variety.coffee <- coffee[,1]

original data after standardization

X_scl <- scale(X)
y_std <- y/sd(y)

#delete two covariates
X_select <- X[, -c(5,11)]
X_sel_scl <- scale(X_select)

#run FBMM with all covariates/selected covaeiates
set.seed(020301)
coffee.fbmm <- simulation_func(X_sel_scl, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5) 
coffee.fbmm.nvs <- simulation_func(X_sel_scl, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5, SS=F) 

coffee.fbmm.all <- simulation_func(X_scl, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5) 
coffee.fbmm.nvs.all <- simulation_func(X_scl, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5, SS=F) 


#run split model with FBMM w. VS
W_scl_1 <- scale(X[, c(5,11)])
Z_scl_1 <- scale(X[, -c(5,11)])
coffee.fbmm.split_1 <- simulation_split(W_scl_1, Z_scl_1, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5) 
coffee.fbmm.split.nvs_1 <- simulation_split(W_scl_1, Z_scl_1, y_std, prior = "Bessel", M_init = 12, Lambda = 3, a_alpha = 5, SS=F) 

#calclulate ppd
size.test <- 5 #size of test data
scl = 50001:1e5 #after burn-in
N = 1e5
thin = 50
yy <- seq(-10, 10, 0.05) #list of y_new values to evaluate the ppd at
#testing data index
test.indx <- sample(1:length(y), size.test, replace = F)

X_train <- scale(X_select[-test.indx, ])
center.x.train <- attr(X_train, "scaled:center")
scale.x.train <- attr(X_train, "scaled:scale")

y_train <- y[-test.indx]/sd(y[-test.indx])
center.y.train <- mean(y[-test.indx])
scale.y.train <- sd(y[-test.indx])

#apply the same scaling of training data to test data of X
X_test_sel <- X_select[test.indx, ]
X_test <- as.matrix(sweep(sweep(X_test_sel, 2, center.x.train), 2, scale.x.train, FUN = "/"))
y_test <- y[test.indx]

### run simulation with training data
coffee.fbmm.train <- simulation_func(X_train, y_train, prior = "Bessel", M_init = 12, Lambda = 7, a_alpha = 5)
coffee.fbmm.train.nss <- simulation_func(X_train, y_train, prior = "Bessel", M_init = 12, Lambda = 7, a_alpha = 5, SS = F)
coffee.fdmm.train <- simulation_func(X_train, y_train, M_init = 12, Lambda = 7, a_alpha = 5)
coffee.fdmm.train.nss <- simulation_func(X_train, y_train, M_init = 12, Lambda = 7, a_alpha = 5, SS = F)
coffee.dpm.train <- simulation_func_rpms(X_train, y_train, k_init = 12, a_alpha = 5)
coffee.dpm.train.nss <- simulation_func_rpms(X_train, y_train, k_init = 12, a_alpha = 5, SS = F)



## calculate ppd
### for M3
sim.res <- coffee.fbmm.train
coffee.y.pred <- matrix(0, ncol = length(yy), nrow = size.test)
for(i in 1:size.test){
    for(j in 1:length(yy)){
        coffee.y.pred[i, j] <- ppd_sim(sim.res, X_test[i,], yy[j] , scl, thin)
    }
}
yy.scl <- yy*scale.y.train + center.y.train
coffee.y.pred.trans <- coffee.y.pred / scale.y.train

#point estimate of y_new and mse value
est.y.test <- apply(coffee.y.pred.trans.scl, 1, function(s) sum(yy.scl * s))
mse.3 <- mean((est.y.test - y_test)^2)





