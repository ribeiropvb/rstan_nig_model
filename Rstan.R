
library(tidyverse)
library(magrittr)
library(MASS)
library(GGally)
library(patchwork)

library(rstan)
library(bayesplot)
#library(coda)
#library(shinystan)
#library(rstanarm)
#library(loo)
#library(tidybayes)



setwd('C:/Users/Pedro/Desktop/Pedro/UFG/8 Periodo/Modelos Lineares/Regressao Bayesiana')
theme_clean <- theme_light() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

# Simulação dos dados ####

set.seed(2022)

# Dados

X <- mvrnorm(
  n = 1000
  , mu = c(3, 6)
  , Sigma = matrix(c(1.5,0,0,3.5),2,2)
) %>% as.data.frame() %>% 
  set_colnames(c('x1','x2'))
X <- model.matrix(~x1+x2,X)

x_1 <- X[,2]
x_2 <- X[,3]

# Priori

sigma <- 1/rgamma(n = 1, shape = 3, rate = 2)
e <- rnorm(n = 1000, mean = 0, sd = sqrt(sigma))
mu_0 <- c(0,0,0)
V <- matrix(sigma*c(12,0,0,0,10,0,0,0,9), ncol = 3, nrow = 3, byrow = T)

betas <- mvrnorm(n = 1, mu = mu_0, Sigma = V)
# The regression model
y <- betas[1] + (betas[2] * x_1) + (betas[3] * x_2) + e
da <- data.frame(y, x_1, x_2)

da %>%
  ggpairs(
    upper = list(continuous = "density", combo = "box_no_facet"),
    lower = list(continuous = "points", combo = "dot_no_facet"))+
  theme_clean

model1 <- 'data {
  int<lower=0> N; 		// Número observações
  int<lower=0> K;		// Número de colunas
  real y[N]; 			// Variável Resposta
  matrix[N,K] X; 		// Matriz modelo
//  matrix[N2,K] new_X; 	// Matriz valores preditos
//  int N2; 				// Tamanho da matriz x_new
}
parameters {
  vector[K] beta; //the regression parameters
  real sigma; //the standard deviation
}
transformed parameters {
  vector[N] linpred;
  linpred = X*beta;
}
model {
	sigma ~ inv_gamma(0.25, 0.5);
	for(i in 1:K){
		beta[i] ~ normal(0,2.5);
	}
	y ~ normal(linpred, sigma);
}'

fit_Chain1 <- stan(
  model_code = model1
  #file = 'my_model.stan'
  , data = list(
    N = dim(X)[1]
    , K = dim(X)[2]
    , X = X
    , y = y
  ), seed = 2203
  , chains = 1
)

model2 <- 'data {
  int<lower=0> N; 		// Número observações
  int<lower=0> K;		// Número de colunas
  real y[N]; 			// Variável Resposta
  matrix[N,K] X; 		// Matriz modelo
//  matrix[N2,K] new_X; 	// Matriz valores preditos
//  int N2; 				// Tamanho da matriz x_new
}
parameters {
  vector[K] beta; //the regression parameters
  real sigma; //the standard deviation
}
transformed parameters {
  vector[N] linpred;
  linpred = X*beta;
}
model {
	sigma ~ inv_gamma(0.05, 0.15);
	for(i in 1:K){
		beta[i] ~ normal(0,2.5);
	}
	y ~ normal(linpred, sigma);
}'

fit_Chain2 <- stan(
  model_code = model2
  #file = 'my_model.stan'
  , data = list(
    N = dim(X)[1]
    , K = dim(X)[2]
    , X = X
    , y = y
  ), seed = 2203
  , chains = 1
)

#save(fit_Chain1, file = 'fit_Chain1')
#save(fit_Chain2, file = 'fit_Chain2')

## Traceplot - Chain 1
p1 <- mcmc_trace(fit_Chain1, pars = c("sigma")) +
  labs (
    x ='Número de Interações'
    , y = expression(paste("Valor de ", sigma^2))
    , title = expression(
      paste(
        'Traceplot - ', sigma^2
          ,': Valores Iniciais: '
          , sigma^2
          , '~ inv_gamma(0.25, 0.5)'
      )
    )
  )+
  theme(plot.title = element_text(hjust = 0.5))

## Traceplot
p2 <- mcmc_trace(fit_Chain2, pars = c("sigma")) +
  labs (
    x ='Número de Interações'
    , y = expression(paste("Valor de ", sigma^2))
    , title = expression(
      paste(
        "Traceplot - ", sigma^2
          ,': Valores Iniciais: '
          , sigma^2
          , '~ inv_gamma(0.05, 0.15)'
      )
    )
  )+
  theme(plot.title = element_text(hjust = 0.5))

p1/p2

params_chain1 <- extract(fit_Chain1)
params_chain2 <- extract(fit_Chain2)


bind_rows(
  data.frame(
    chain = 'Cadeia 1'
    , value = params_chain1$sigma
  )
  ,data.frame(
    chain = 'Cadeia 2'
    , value = params_chain2$sigma
  )
) %>% 
  as_tibble() %>% 
  ggplot(aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#FF6666")+
  facet_wrap(~chain)+
  theme_light()

par <- cbind(params_chain1$beta,params_chain1$sigma) %>% 
  magrittr::set_colnames(
    c('beta_0','beta_1','beta_2','sigma')
  ) %>% 
  as.data.frame()

mcmc_trace(fit_Chain2, pars = c("beta[1]", "beta[2]", "beta[3]", "sigma")) +
  labs (
    x ='Número de Interações'
    , y = expression(paste("Valor de ", sigma^2))
  )+
  theme(plot.title = element_text(hjust = 0.5))

bind_cols(
  tibble(params_chain1$beta[,1]),
  tibble(params_chain1$beta[,2]),
  tibble(params_chain1$beta[,3]),
  tibble(params_chain1$sigma)
) %>% 
  magrittr::set_colnames(
    c('Beta_0','Beta_1','Beta_2','Sigma')
  ) %>% 
  mcmc_hist()

# Intervalo de Credibilidade ####

library(posterior)
marg_dist<- as_draws_df(fit_Chain1)

sum_fit <- summarise_draws(
  marg_dist[,1:4], "mean", "median", "sd", "mad"
  ,  ~quantile(.x, probs = c(0.005,0.025,0.995,0.975))
  , default_convergence_measures()
)

sum_fit

# Cobertura ####

model_cobertura <- 'data {
  int<lower=0> N; 		// Número observações
  int<lower=0> K;		// Número de colunas
  real y[N]; 			// Variável Resposta
  matrix[N,K] X; 		// Matriz modelo
//  matrix[N2,K] new_X; 	// Matriz valores preditos
//  int N2; 				// Tamanho da matriz x_new
}
parameters {
  vector[K] beta; //the regression parameters
  real sigma; //the standard deviation
}
transformed parameters {
  vector[N] linpred;
  linpred = X*beta;
}
model {
	sigma ~ inv_gamma(0.25, 0.5);
	for(i in 1:K){
		beta[i] ~ normal(0,2.5);
	}
	y ~ normal(linpred, sigma);
}'

fit_cobertura <- stan(
  model_code = model_cobertura
  #file = 'my_model.stan'
  , data = list(
    N = dim(X)[1]
    , K = dim(X)[2]
    , X = X
    , y = y
  ), seed = 2203
)


params <- extract(fit_cobertura)
n <- length(X[,1])
lm_model <- lm(y~X[,2]+X[,3])

tibble(
  `Parâmetro` = c(
    'Beta_0','Beta_1'
    , 'Beta_2','sigma^2'
  ),
  `Parâmetros Originais` = c(
    betas[1],betas[2]
    ,betas[3], sigma
  ),
  `Modelo Bayesiano` = c(
    mean(params$beta[,1])
    , mean(params$beta[,2])
    , mean(params$beta[,3])
    , mean(params$sigma)
  ),
  `Modelo Linear` = c(
    lm_model$coefficients[1]
    , lm_model$coefficients[2]
    , lm_model$coefficients[3]
    , sum(lm_model$residuals^2)/(n-3)
  )
)

marg_cobertura<- as_draws_df(fit_cobertura)

sum_cobertura <- summarise_draws(
  marg_cobertura[,1:4], "mean", "median", "sd", "mad"
  ,  ~quantile(.x, probs = c(0.005,0.025,0.975,0.995))
  , default_convergence_measures()
)

params <- extract(fit_cobertura)

bind_cols(
  tibble(
    beta_0 = params$beta[,1]
    , beta_1 = params$beta[,2]
    , beta_2 = params$beta[,3]
    , sigma = params$sigma
  ) %>% mutate(
    beta_0_cob_5 = between(beta_0, sum_cobertura$`2.5%`[1], sum_cobertura$`97.5%`[1]) %>% 
      as.numeric(),
    beta_1_cob_5 = between(beta_1, sum_cobertura$`2.5%`[2], sum_cobertura$`97.5%`[2]) %>% 
      as.numeric(),
    beta_2_cob_5 = between(beta_2, sum_cobertura$`2.5%`[3], sum_cobertura$`97.5%`[3]) %>% 
      as.numeric(),
    sigma_cob_5 = between(sigma, sum_cobertura$`2.5%`[4], sum_cobertura$`97.5%`[4]) %>% 
      as.numeric()
  ) %>% summarise(
    beta_0_cob_5 = sum(beta_0_cob_5)/4000,
    beta_1_cob_5 = sum(beta_1_cob_5)/4000,
    beta_2_cob_5 = sum(beta_2_cob_5)/4000,
    sigma_cob_5 = sum(sigma_cob_5)/4000
  ) %>% reshape2::melt(),
  tibble(
    beta_0 = params$beta[,1]
    , beta_1 = params$beta[,2]
    , beta_2 = params$beta[,3]
    , sigma = params$sigma
  ) %>% mutate(
    beta_0_cob_1 = between(beta_0, sum_cobertura$`0.5%`[1], sum_cobertura$`99.5%`[1]) %>% 
      as.numeric(),
    beta_1_cob_1 = between(beta_1, sum_cobertura$`0.5%`[2], sum_cobertura$`99.5%`[2]) %>% 
      as.numeric(),
    beta_2_cob_1 = between(beta_2, sum_cobertura$`0.5%`[3], sum_cobertura$`99.5%`[3]) %>% 
      as.numeric(),
    sigma_cob_1 = between(sigma, sum_cobertura$`0.5%`[4], sum_cobertura$`99.5%`[4]) %>% 
      as.numeric()
  ) %>% summarise(
    beta_0_cob_1 = sum(beta_0_cob_1)/4000,
    beta_1_cob_1 = sum(beta_1_cob_1)/4000,
    beta_2_cob_1 = sum(beta_2_cob_1)/4000,
    sigma_cob_1 = sum(sigma_cob_1)/4000
  ) %>% reshape2::melt() %>% 
    dplyr::select(-variable)
) %>% magrittr::set_colnames(
  c('Variável', 'Cobertura 95%', 'Cobertura 99%')
) %>% mutate(
  `Variável` = c('Beta 0', 'Beta 1','Beta 2','Sigma^2')
)




# %>% mutate(
#   cob_beta_0_5 = ifelse(
#     beta_0 > sum_cobertura$`2.5%`[1] & beta_0 < sum_cobertura$`97.5%`[1],
#     1, 0
#   ),
#   cob_beta_1_5 = ifelse(
#     beta_1 > sum_cobertura$`2.5%`[2] & beta_1 < sum_cobertura$`97.5%`[2],
#     1, 0
#   ),
#   cob_beta_2_5 = ifelse(
#     beta_2 > sum_cobertura$`2.5%`[3] & beta_1 < sum_cobertura$`97.5%`[3],
#     1, 0
#   )
# )


# Intervalo de Estimação ####

library(DepthProc)

model_est <- 'data {
  int<lower=0> N; 		// Número observações
  int<lower=0> K;		// Número de colunas
  real y[N]; 			// Variável Resposta
  matrix[N,K] X; 		// Matriz modelo
//  matrix[N2,K] new_X; 	// Matriz valores preditos
//  int N2; 				// Tamanho da matriz x_new
}
parameters {
  vector[K] beta; //the regression parameters
  real sigma; //the standard deviation
}
transformed parameters {
  vector[N] linpred;
  linpred = X*beta;
}
model {
	sigma ~ inv_gamma(0.25, 0.5);
	for(i in 1:K){
		beta[i] ~ normal(0,2.5);
	}
	y ~ normal(linpred, sigma);
}'

fit_est <- stan(
  model_code = model_est
  #file = 'my_model.stan'
  , data = list(
    N = dim(X)[1]
    , K = dim(X)[2]
    , X = X
    , y = y
  ), seed = 2203
#  , chains = 1
)

params_est <- extract(fit_est)
dim(X)
dim(params_est$beta)

y_est <- matrix(NA, nrow = dim(X)[1], ncol = 3)
for (i in 1:dim(X)[1]) {
  y_est[i,1] <- sum(X[i,1]*params_est$beta[i,1])
  y_est[i,2] <- sum(X[i,2]*params_est$beta[i,2])
  y_est[i,3] <- sum(X[i,3]*params_est$beta[i,3])
}
functional_boxplot <- fncBoxPlot(y_est, bands = c(0, 0.95, 0.99))
functional_boxplot$data

functional_boxplot + theme_clean







