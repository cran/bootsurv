#' Direct Bootstrap Methods for Survey Data
#'
#'
#' The function `direct.boot.stsrs` applies one of the following bootstrap methods on complete (full response) survey data selected under either SRSWOR or STSRSWOR: Efron (1979), McCarthy and Snowden (1985), Rao and Wu (1988) and Sitter (1992, JASA).
#'
#' @param data A vector, matrix or data frame. If it is a matrix or data frame then the column of study variable has to be named `study.variable`.  If the sampling design is STSRSWOR, a column identifying strata named `stratum` has to be included.
#' @param population.size A vector of stratum population sizes
#' @param R The number of bootstrap replicates
#' @param parameter One of the following population parameters can be applied: `"total"` (population total), `"mean"` (population mean), `"quartile.25"` (population 1st quartile), `"quartile.50"` or `"median"` (population median) or `"quartile.75"` (population 3rd quartile). If the parameter of interest is the population mean or total, the HT-estimator is applied. If the parameter of interest is a population quartile, the estimator in Sarndal, Swensson, and Wretman (1992, Chapter 5) is applied. The default is the population total.
#' @param bootstrap.method One of the following bootstrap methods can be applied: `"Efron"` (Efron, 1979), `"McCarthy.Snowden"` (McCarthy and Snowden, 1985), `"Rao.Wu"` (Rao and Wu, 1988) or `"Sitter.BMM"` (Sitter, 1992). The default is `"Rao.Wu"`.
#' @param boot.sample.size If the method of Rao and Wu (1988) is applied, a vector of bootstrap sample sizes for each stratum may be specified. The length of this vector must match the number of strata. By default, if 'boot.sample.size' is not specified, the bootstrap sample size within each stratum will be 'nh-3', where 'nh' is the original sample size in stratum 'h'.
#'
#' @return
#'
#' `boot.statistic` A vector of bootstrap statistics
#'
#' `boot.var` The bootstrap variance estimator of the estimator of the parameter of interest
#'
#' `boot.mean` The average of the bootstrap estimator of the parameter of interest
#'
#' `boot.sample` For each iteration, a list of results is generated, including three columns: bootstrap values (which may be rescaled values if resampling is done on a rescaled version of the original sample), selected indices in each stratum, and a stratum identifier column.
#'
#' @export
#'
#' @importFrom stats aggregate rbinom rexp rnorm runif var
#'
#' @references
#' Efron, B. (1979). Bootstrap methods: another look at the jackknife. The Annals of Statistics 7 (1), 1–26.
#'
#' McCarthy, P. J. and C. B. Snowden (1985). The bootstrap and finite population sampling. Vital and Health Statistics, Series 2, No. 95. DHHS Publication No. (PHS) 85–1369. Public Health Service. Washington. U.S. Government Printing Office.
#'
#' Rao, J. N. K. and C. F. J. Wu (1988). Resampling inference with complex survey data. Journal of the American Statistical Association 83 (401), 231–241.
#'
#' Särndal, C.-E., Swensson, B. & Wretman, J. (1992). Model-Assisted Survey Sampling. NewYork: Springer.
#'
#' Sitter, R. R. (1992). A resampling procedure for complex survey data. Journal of the American Statistical Association 87 (419), 755–765.
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
#' boot.RW<- direct.boot.stsrs(data_samp_srs, population_size, R)
#' boot.RW$boot.var
#'
#' boot.Efron<- direct.boot.stsrs(data_samp_srs, population_size, R,
#'               parameter="total", bootstrap.method="Efron")
#' boot.Efron$boot.var
#'
#' boot.RW.med<- direct.boot.stsrs(data_samp_srs, population_size, R,
#'                parameter="median")
#' boot.RW.med$boot.var
#'
#' data(data_samp_stsrs)
#' population_size_st<- c(4500, 6300, 3500, 2000, 1500)
#' # The overall sampling fraction is about 30%.
#' # data_samp_stsrs is a sample taken from data_pop_st available in the package.
#'
#' boot.RW.st<- direct.boot.stsrs(data_samp_stsrs, population_size_st, R,
#'               parameter="total", bootstrap.method="Rao.Wu")
#' boot.RW.st$boot.statistic
#'
#'
#'
direct.boot.stsrs<- function(data, population.size, R, parameter="total", bootstrap.method="Rao.Wu", boot.sample.size=NULL) {

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

    if(!is.null(bootstrap.method) & bootstrap.method!="Rao.Wu" & bootstrap.method!="Efron" & bootstrap.method!="Sitter.BMM" & bootstrap.method!="McCarthy.Snowden")
      stop("One of the following bootstrap methods should be specified: Efron, Rao.Wu, Sitter.BMM or McCarthy.Snowden. The default is Rao.Wu.")

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
      if(bootstrap.method=="Rao.Wu" || bootstrap.method=="Efron" || bootstrap.method=="Sitter.BMM" || bootstrap.method=="McCarthy.Snowden"){
        if(length(population.size)==L){

          boot.sample<- list()
          boot.statistic<- NULL

          for (b in 1:R){

            boot.sample.index<- NULL
            boot.values<- NULL

            for(h in 1:L){
              n<- strata.size[h]
              N<- population.size[h]
              f<- n/N

              if(bootstrap.method=="Efron"){
                Boot.Rescaling.factor<- 1
                n.subboot<- 1
                k.boot<- n
              }

              if(bootstrap.method=="McCarthy.Snowden"){
                Boot.Rescaling.factor<- 1
                n.subboot<- 1
                k.boot<- floor((n-1)/(1-f) + 0.5)
              }

              if(bootstrap.method=="Rao.Wu"){
                if(is.null(boot.sample.size)) k.boot<- n-3
                else k.boot<- boot.sample.size[h]
                Boot.Rescaling.factor<- sqrt((1-f)*k.boot/(n-1))
                n.subboot<- 1
              }

              if(bootstrap.method=="Sitter.BMM"){
                n.subboot<- floor(n/(2-f))
                Boot.Rescaling.factor<- 1
                f.subboot<- n.subboot/n
                k<- (n*(1-f.subboot))/(n.subboot*(1-f))
                v<- (1/floor(k)-1/k)/(1/floor(k)-1/ceiling(k))
                Ind<- rbinom(1,1,v)
                k.boot<- floor(k) + Ind
              }

              n.bootstrap<- n.subboot*k.boot

              mean.data<- mean(data$study.variable[data$stratum==strata[h]])
              Rescaled.data<- mean.data + Boot.Rescaling.factor * (data$study.variable - mean.data)

              sample.index<- NULL
              for(i in 1:k.boot)  sample.index<- c(sample.index, sample(1:n, n.subboot, replace = FALSE))
              boot.sample.index<- rbind(boot.sample.index, cbind(Rescaled.data[sample.index], sample.index, rep(strata[h], n.bootstrap)) )
              boot.values<- rbind(boot.values, cbind(Rescaled.data[sample.index], rep(population.size[h]/n.bootstrap, n.bootstrap)) )
            }

            boot.statistic[b]<- statistic(boot.values[ ,1], boot.values[ ,2], q)
            if(parameter=="mean") boot.statistic[b]<- statistic(boot.values[ ,1], boot.values[ ,2], q)/sum(population.size)
            colnames(boot.sample.index)<- c("bootstrap.sample", "ID", "stratum")
            boot.sample[[b]]<- as.data.frame(boot.sample.index)
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
    }
  }
  return
}

