#' Bootstrap Weights Methods for Survey Data
#'
#'The function `boot.weights.stsrs` applies one of the following bootstrap weights methods on complete (full response) survey data selected under either SRSWOR or STSRSWOR: Rao, Wu and Yue (1992), Bertail and Combris (1997), Chipperfield and Preston (2007) and Beaumont and Patak (2012)
#'
#' @param data A vector, matrix or data frame. If it is a matrix or data frame then the column of study variable has to be named `study.variable`.  If the sampling design is STSRSWOR, a column identifying strata named `stratum` has to be included.
#' @param population.size A vector of stratum population sizes
#' @param R The number of bootstrap replicates
#' @param parameter One of the following population parameters can be applied: `"total"` (population total), `"mean"` (population mean), `"quartile.25"` (population 1st quartile), `"quartile.50"` or `"median"` (population median) or `"quartile.75"` (population 3rd quartile). If the parameter of interest is the population mean or total, the HT-estimator is applied. If the parameter of interest is a population quartile, the estimator in Sarndal, Swensson, and Wretman (1992, Chapter 5) is applied. The default is the population total.
#' @param bootstrap.method One of the following bootstrap methods can be applied: `"Rao.Wu.Yue"` (Rao, Wu and Yue, 1992),`"Bertail.Combris"` (Bertail and Combris, 1997), `"Chipperfield.Preston"` (Chipperfield and Preston, 2007) or `"Beaumont.Patak"` (Beaumont and Patak, 2012). The default is `"Rao.Wu.Yue"`.
#' @param boot.sample.size A vector of bootstrap sample sizes within strata only required for the method of Rao, Wu and Yue (1992). The length of this vector has to be the same as the number of strata. The default is NULL. If the method of Rao, Wu and Yue (1992) is applied and `boot.sample.size` is not specified, the bootstrap sample size will be `nh-1` within each stratum, where `nh` is the original sample size within stratum `h`.
#' @param distribution.adjust The default is NULL. A distribution should be specified for the method of Bertail and Combris (1997) and Beaumont and Patak (2012) to generate the bootstrap weight adjustments if epsilon is NULL. One of the following distribution can be used: `"Normal"`, `"Lognormal"`, `"Exponential"` or `"Uniform"`.
#' @param epsilon The default is NULL. If either Bertail and Combris (1997) or Beaumont and Patak (2012) is applied and `distribution.adjust` is NULL, a value must be given to epsilon so that Eq(5) in Beaumont and Patak (2012) can be used to generate the bootstrap weight adjustments.
#'
#' @return
#'
#' `boot.statistic` A vector of bootstrap statistics
#'
#' `boot.var` The bootstrap variance estimator of the estimator of parameter of interest.
#'
#' `boot.mean` The average of the bootstrap estimator of the parameter of interest.
#'
#' `boot.sample` A list of results for each iteration. That includes a column of original sample values, a column of bootstrap weight adjustments, a column of bootstrap weights and a column of stratum identifier.
#'
#' @export
#'
#' @import MASS
#'
#' @importFrom stats aggregate rbinom rexp rnorm runif var
#'
#' @references
#'
#' Beaumont, J.-F. and Patak, Z. (2012). On the generalized bootstrap for sample surveys with special attention to Poisson sampling. International Statistical Review 80 (1), 127–148.
#'
#' Bertail, P. and Combris, P. (1997). Bootstrap généralisé d’un sondage. Annales d’économie et de statistique 46, 49–83.
#'
#' Chipperfield, J. and Preston, J. (2007). Efficient bootstrap for business surveys. Survey Methodology 33 (2), 167–172.
#'
#' Rao, J. N. K., Wu, C. F. J. and Yue, K. (1992). Some recent work on resampling methods for complex surveys. Survey Methodology 18 (2), 209–217.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). Model-Assisted Survey Sampling. NewYork: Springer.
#'
#'
#' @examples
#'
#' R<- 20
#'
#' data(data_samp_srs)
#' population_size<- 6000
#' # The sampling fraction is about 30%.
#' # data_samp_srs is a sample taken from data_pop available in the package.
#'
#' boot.RWY<- boot.weights.stsrs(data_samp_srs, population_size, R)
#' boot.RWY$boot.var
#'
#' boot.CP<- boot.weights.stsrs(data_samp_srs, population_size, R,
#'            bootstrap.method="Chipperfield.Preston")
#' boot.CP$boot.var
#'
#' boot.BP.med<- boot.weights.stsrs(data_samp_srs, population_size, R,
#'                parameter="median", bootstrap.method="Beaumont.Patak",
#'                distribution.adjust="Exponential")
#' boot.BP.med$boot.var
#' boot.BP.med$boot.sample[[5]]
#'
#'
#' data(data_samp_stsrs)
#' population_size_st<- c(4500, 6300, 3500, 2000, 1500)
#' # The overall sampling fraction is about 30%.
#' # data_samp_stsrs is a sample taken from data_pop_st available in the package.
#'
#' boot.RWY.st<- boot.weights.stsrs(data_samp_stsrs, population_size_st, R)
#' boot.RWY.st$boot.var
#' boot.RWY.st$boot.statistic
#'
#'
#'
boot.weights.stsrs<- function(data, population.size, R, parameter="total", bootstrap.method="Rao.Wu.Yue", boot.sample.size=NULL, distribution.adjust=NULL, epsilon=NULL){

  Total.HT<- function(y, w, q)  sum(y*w)
  Quantiles<- function(y, w, q){
    n<- nrow(y)
    y.w<- cbind(y,w/sum(w))
    y.w.order<- y.w[order(y),]
    i<-1
    Sum.Weight<- 0
    while(Sum.Weight < q){
      Sum.Weight<- Sum.Weight + y.w.order[i,2]
      i<- i+1
    }
    if(Sum.Weight == q) (y.w.order[i-1, 1]+y.w.order[i, 1])/2
    else y.w.order[i-1, 1]
  }

  if(is.null(data) || is.null(population.size) || is.null(R) || is.null(parameter) || is.null(bootstrap.method))
    stop("The following arguments have to be nonnull: data, population.size, R, parameter, bootstrap.method.")

  if(!is.null(data) & !is.null(population.size) & !is.null(R) & !is.null(parameter) & !is.null(bootstrap.method)){
    q<- NULL
    if(parameter=="total" || parameter=="mean") statistic<- Total.HT
    if(parameter=="quartile.25"){
      statistic<- Quantiles
      q<- 0.25
    }
    if(parameter=="quartile.50" || parameter=="median"){
      statistic<- Quantiles
      q<- 0.5
    }
    if(parameter=="quartile.75"){
      statistic<- Quantiles
      q<- 0.75
    }

    if(is.vector(data)) {
      data<- as.data.frame(data)
      names(data)<- "study.variable"
    }
    if(is.matrix(data)) data<- as.data.frame(data)

    if(!is.null(parameter) & parameter!="total" & parameter!="mean" & parameter!="median" & parameter!="quartile.25" & parameter!="quartile.50" & parameter!="quartile.75")
      stop("The parameter of interest has to be: total, mean, median (or quartile.50), quartile.25 or quartile.75. The default is the population total.")

    if(!is.null(bootstrap.method) & bootstrap.method!="Rao.Wu.Yue" & bootstrap.method!="Chipperfield.Preston" & bootstrap.method!="Beaumont.Patak" & bootstrap.method!="Bertail.Combris")
      stop("One of the following bootstrap methods should be specified: 'Rao.Wu.Yue' (Rao, Wu and Yue, 1992), 'Chipperfield.Preston' (Chipperfield and Preston, 2007), Beaumont.Patak (Beaumont and Patak, 2012) and 'Bertail.Combris' (Bertail and Combris, 1997)'.")

    if(is.null(data$study.variable))
      stop("The data set must contain a culumn with the values of the variable of interest called 'study.variable' and a culumn that identifies the strata called 'stratum' if the survey design is STSRSWOR.")

    if(is.null(data$stratum)) data$stratum<- 1

    data<- data[with(data, order(stratum)),  ]

    strata<- unique(data$stratum)
    strata.size<- table(data$stratum)
    L<- length(strata)

    if(length(population.size)!=L)
      stop("The length of the arguments population.size and the number of strata must be the same. Each element of population.size is the number of units in each subpopulation (or stratum).")

    if(parameter=="total" || parameter=="mean" || parameter=="median" || parameter=="quartile.25" || parameter=="quartile.50" || parameter=="quartile.75"){
      if(length(population.size)==L){

        if(bootstrap.method=="Rao.Wu.Yue" || bootstrap.method=="Chipperfield.Preston"){
          boot.sample<- list()
          boot.statistic<- NULL

          for (b in 1:R){

            boot.sample.adjustment<- NULL

            for(h in 1:L){
              n<- strata.size[h]
              N<- population.size[h]
              f<- n/N

              if(bootstrap.method=="Rao.Wu.Yue"){
                if(is.null(boot.sample.size)) n.bootstrap<- n-1
                else n.bootstrap<- boot.sample.size[h]
                sample.index<- sample(1:n, n.bootstrap, replace = TRUE)
                m.bootstrap<- tabulate(sample.index, nbins=n)
                bootstrap.adjustment<- 1+sqrt((1-f)*n.bootstrap/(n-1)) * (n*m.bootstrap/n.bootstrap-1)
              }

              if(bootstrap.method=="Chipperfield.Preston"){
                n.bootstrap<- floor(n/2)
                sample.index<- sample(1:n, n.bootstrap, replace = FALSE)
                m.bootstrap<- tabulate(sample.index, nbins=n)
                bootstrap.adjustment<- 1+sqrt( (1-f)*n.bootstrap/(n-n.bootstrap)) * (n*m.bootstrap/n.bootstrap-1)
              }

              boot.sample.adjustment<- rbind(boot.sample.adjustment, cbind(data$study.variable[data$stratum==strata[h]], bootstrap.adjustment, bootstrap.adjustment*rep(N/n,n), rep(strata[h], n)) )
            }

            boot.statistic[b]<- statistic(as.numeric(boot.sample.adjustment[,1]), as.numeric(boot.sample.adjustment[,3]), q)
            if(parameter=="mean") boot.statistic[b]<- boot.statistic[b]/sum(population.size)

            colnames(boot.sample.adjustment)<- c("study.variable", "weight.adjustment", "bootstrap.weight", "stratum")
            boot.sample[[b]]<- as.data.frame(boot.sample.adjustment)
          }

          boot.mean<- mean(boot.statistic)
          boot.var<- var(boot.statistic)

          return<- list(
            boot.statistic=boot.statistic,
            boot.var=boot.var,
            boot.mean=boot.mean,
            boot.sample=boot.sample,
            parameter=parameter,
            number.bootstrap.replicates=R
          )
        }

        if(bootstrap.method=="Beaumont.Patak" || bootstrap.method=="Bertail.Combris"){

          if(is.null(distribution.adjust) && is.null(epsilon))
            stop("distribution.adjust and epsilon are NULL by default. However, if either the method of Bertail and Combris (1997) or that of Beaumont and Patak (2012) is applied, either a positive value must be given to epsilon or one of the following distribution should be chosen in order to generate the bootstrap weight adjustments: Normal, Uniform, Exponential or Lognormal.")

          if(!is.null(distribution.adjust) && !is.null(epsilon))
            stop("Either distribution.adjust or epsilon must be NULL.")

          if(!is.null(distribution.adjust) && is.null(epsilon)){
            if(distribution.adjust!="Normal" && distribution.adjust!="Uniform" && distribution.adjust!="Exponential" && distribution.adjust!="Lognormal")
              stop("if either the method of Bertail and Combris (1997) or that of Beaumont and Patak (2012) is applied and epsilon is NULL, one of the following distribution should be chosen in order to generate the bootstrap weight adjustments: Normal, Uniform, Exponential or Lognormal.")
            else{
              boot.sample<- list()
              boot.statistic<- NULL

              for (b in 1:R){

                boot.sample.adjustment<- NULL

                for(h in 1:L){
                  n<- strata.size[h]
                  N<- population.size[h]
                  f<- n/N

                  if(distribution.adjust=="Normal"){
                    Sigma<- matrix(-(1-f)/(n-1), nrow=n, ncol=n)
                    diag(Sigma)<- 1-f
                    bootstrap.adjustment<- mvrnorm(1, rep(1,  n), Sigma)
                  }
                  if(distribution.adjust=="Uniform"){
                    interval<- 1+c(-1, 1) * sqrt(3*(1-n/N))
                    bootstrap.adjustment<- runif(n, interval[1], interval[2])
                  }
                  if(distribution.adjust=="Exponential"){
                    ExpV<- rexp(n, 1)
                    bootstrap.adjustment<- 1+(ExpV-1)*sqrt(1-n/N)
                  }
                  if(distribution.adjust=="Lognormal"){
                    lognormV<- rnorm(n, (-0.5)*log(2-n/N), sqrt(log(2-n/N)))
                    bootstrap.adjustment<- exp(lognormV)
                  }
                  tau<- 1.1
                  while(sum(bootstrap.adjustment<0)>0){
                    bootstrap.adjustment[bootstrap.adjustment<0]<- (bootstrap.adjustment[bootstrap.adjustment<0]+tau-1)/tau
                    tau<- tau+0.1
                  }
                  boot.sample.adjustment<- rbind(boot.sample.adjustment, cbind(data$study.variable[data$stratum==strata[h]], bootstrap.adjustment, bootstrap.adjustment*rep(N/n,n), rep(strata[h], n)) )
                }

                boot.statistic[b]<- statistic(as.numeric(boot.sample.adjustment[,1]), as.numeric(boot.sample.adjustment[,3]), q)
                if(parameter=="mean") boot.statistic[b]<- boot.statistic[b]/sum(population.size)

                colnames(boot.sample.adjustment)<- c("study.variable", "weight.adjustment", "bootstrap.weight", "stratum")
                boot.sample[[b]]<- as.data.frame(boot.sample.adjustment)
              }

              boot.mean<- mean(boot.statistic)
              boot.var<- var(boot.statistic)

              return<- list(
                boot.statistic=boot.statistic,
                boot.var=boot.var,
                boot.mean=boot.mean,
                boot.sample=boot.sample,
                parameter=parameter,
                number.bootstrap.replicates=R
              )
            }
          }

          if(is.null(distribution.adjust) && !is.null(epsilon)){

            non.positive.semidefinite<- 0
            for(h in 1:L){
              n<- strata.size[h]
              N<- population.size[h]
              f<- n/N
              Sigma<- matrix((-1)*(1-f)/(n-1), nrow=n, ncol=n)
              diag(Sigma)<- 1-f
              Matrix.Eigen <- eigen(Sigma)$vectors
              Matrix.Diag <- diag(eigen(Sigma)$values)
              if((sum(Matrix.Diag<0)>0)) non.positive.semidefinite<- non.positive.semidefinite+1
            }

            if(non.positive.semidefinite==0){

              boot.sample<- list()
              boot.statistic<- NULL

              for (b in 1:R){

                boot.sample.adjustment<- NULL

                for(h in 1:L){
                  n<- strata.size[h]
                  N<- population.size[h]
                  f<- n/N
                  Sigma<- matrix((-1)*(1-f)/(n-1), nrow=n, ncol=n)
                  diag(Sigma)<- 1-f
                  Matrix.Eigen <- eigen(Sigma)$vectors
                  Matrix.Diag <- diag(eigen(Sigma)$values)
                  Sigma.sqrt<- Matrix.Eigen %*% sqrt(Matrix.Diag) %*% solve(Matrix.Eigen)
                  Ind<- rbinom(n, 1, 1/(1+epsilon^2))
                  a.tilde<- (-1)*epsilon*Ind + 1/epsilon*(1-Ind)
                  bootstrap.adjustment<- 1 + Sigma.sqrt %*% a.tilde

                  if(sum(bootstrap.adjustment<0)>0){
                    tau<- 1.1
                    while(sum(bootstrap.adjustment<0)>0){
                      bootstrap.adjustment[bootstrap.adjustment<0]<- (bootstrap.adjustment[bootstrap.adjustment<0]+tau-1)/tau
                      tau<- tau+0.1
                    }
                  }
                  boot.sample.adjustment<- rbind(boot.sample.adjustment, cbind(data$study.variable[data$stratum==strata[h]], bootstrap.adjustment, bootstrap.adjustment*rep(N/n,n), rep(strata[h], n)) )
                }

                boot.statistic[b]<- statistic(as.numeric(boot.sample.adjustment[,1]), as.numeric(boot.sample.adjustment[,3]), q)
                if(parameter=="mean") boot.statistic[b]<- boot.statistic[b]/sum(population.size)

                colnames(boot.sample.adjustment)<- c("study.variable", "weight.adjustment", "bootstrap.weight", "stratum")
                boot.sample[[b]]<- as.data.frame(boot.sample.adjustment)
              }

              boot.mean<- mean(boot.statistic)
              boot.var<- var(boot.statistic)

              return<- list(
                boot.statistic=boot.statistic,
                boot.var=boot.var,
                boot.mean=boot.mean,
                boot.sample=boot.sample,
                parameter=parameter,
                number.bootstrap.replicates=R
              )
            }

            if(non.positive.semidefinite>0)
              stop("The diagonal matrix of eigenvalues contains negative values, therefore, another distribution should be used to generate the bootstrap weight adjustments.")
          }
        }
      }
    }
  }
  return
}
