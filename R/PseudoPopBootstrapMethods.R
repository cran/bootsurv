#' Pseudo-population Bootstrap Methods for Survey Data
#'
#' The function `pseudopop.boot.stsrs` applies one of the following pseudo-population bootstrap methods on complete (full response) survey data selected under either SRSWOR or STSRSWOR: Bickel and Freedman (1984), Chao and Lo (1985), Sitter (1992, CJS), Booth, Butler and Hall (1994) and Chao and Lo (1994).
#'
#' @param data A vector, matrix or data frame. If it is a matrix or data frame then the column of study variable has to be named `study.variable`.  If the sampling design is STSRSWOR, a column identifying strata named `stratum` has to be included.
#' @param population.size A vector of stratum population sizes
#' @param R.pop The number of bootstrap replicates to create bootstrap pseudo-populations
#' @param R.samp The number of bootstrap replicates to draw bootstrap samples from each bootstrap pseudo-population
#' @param parameter One of the following population parameters can be applied: `"total"` (population total), `"mean"` (population mean), `"quartile.25"` (population 1st quartile), `"quartile.50"` or `"median"` (population median) or `"quartile.75"` (population 3rd quartile). If the parameter of interest is the population mean or total, the HT-estimator is applied. If the parameter of interest is a population quartile, the estimator in Sarndal, Swensson, and Wretman (1992, Chapter 5) is applied. The default is the population total.
#' @param bootstrap.method One of the following bootstrap methods can be applied: `"Bickel.Freedman"` (Bickel and Freedman, 1984),`"Chao.Lo.1985"` (Chao and Lo, 1985), `"Sitter.BWO"` (Sitter, 1992), `"Booth.Butler.Hall"` (Booth, Butler and Hall, 1994) or `"Chao.Lo.1994"` (Chao and Lo, 1994). The default is `"Booth.Butler.Hall"`.
#'
#' @return
#'
#' `boot.statistic` A vector of bootstrap statistics
#'
#' `boot.parameter` A vector of bootstrap parameters computed on bootstrap pseudo-populations
#'
#' `boot.var` The bootstrap variance estimator of the estimator of parameter of interest
#'
#' `boot.mean` The average of the bootstrap estimator of the parameter of interest
#'
#' `boot.sample` A list of size `R.pop`. Each list contains a list of results on each generated bootstrap pseudo-population. This includes three columns: bootstrap values, selected indices in each stratum, and a stratum identifier column.
#'
#' @export
#'
#' @importFrom stats aggregate rbinom rexp rnorm runif var
#'
#' @references
#' Bickel, P. J. and Freedman, D. A. (1984). Asymptotic normality and the bootstrap in stratified sampling. The Annals of Statistics 12, 470–82.
#'
#' Booth, J. G., Butler, R. W. and Hall, P. (1994). Bootstrap methods for finite populations. Journal of the American Statistical Association 89 (428), 1282–1289.
#'
#' Chao, M. T. and Lo, S.-H. (1985). A bootstrap method for finite population. Sankhya: The Indian Journal of Statistics, Series A 47, 399–405.
#'
#' Chao, M. T. and Lo, S.-H. (1994). Maximum likelihood summary and the bootstrap method in structured finite populations. Statistica Sinica 4 (2), 389–406.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). Model-Assisted Survey Sampling. NewYork: Springer.
#'
#' Sitter, R. R. (1992). Comparing three bootstrap methods for survey data. The Canadian Journal of Statistics 20 (2), 135–154.
#'
#' @examples
#'
#' R.pop<- 5
#' R.samp<- 10
#'
#' data(data_samp_srs)
#' population_size<- 6000
#' # The sampling fraction is about 30%.
#' # data_samp_srs is a sample taken from data_pop available in the package.
#'
#' boot.Booth<- pseudopop.boot.stsrs(data_samp_srs, population_size, R.pop, R.samp)
#' boot.Booth$boot.var
#'
#' boot.BF<- pseudopop.boot.stsrs(data_samp_srs, population_size, R.pop, R.samp,
#'            bootstrap.method="Bickel.Freedman")
#' boot.BF$boot.var
#'
#' boot.Sitter.med<- pseudopop.boot.stsrs(data_samp_srs, population_size, R.pop,
#'                    R.samp, parameter="median", bootstrap.method="Sitter.BWO")
#' boot.Sitter.med$boot.var
#' boot.Sitter.med$boot.sample[[2]][[5]]
#'
#' data(data_samp_stsrs)
#' population_size_st<- c(4500, 6300, 3500, 2000, 1500)
#' # The overall sampling fraction is about 30%.
#' # data_samp_stsrs is a sample taken from data_pop_st available in the package.
#'
#' boot.Booth.st<- pseudopop.boot.stsrs(data_samp_stsrs, population_size_st, R.pop, R.samp)
#' boot.Booth.st$boot.statistic
#'
#'
#'
pseudopop.boot.stsrs<- function(data, population.size, R.pop, R.samp, parameter="total", bootstrap.method="Booth.Butler.Hall") {

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
  G.CL.1985<- function(t, n) (1-n/t)*t*(n-1)/((t-1)*n)
  stratified.srs<- function(Nh, nh, design){
    L<- length(Nh)
    sample.index<- NULL
    for(h in 1:L){
      if(design=="SRSWR") sample.index<- rbind(sample.index, cbind(sample(1:Nh[h], nh[h], replace = TRUE), rep(h,nh[h])) )
      if(design=="SRSWOR") sample.index<- rbind(sample.index, cbind(sample(1:Nh[h], nh[h], replace = FALSE), rep(h,nh[h])) )
    }
    sample.index
  }

  if(is.null(data) || is.null(population.size) || is.null(R.pop) || is.null(R.samp) || is.null(parameter) || is.null(bootstrap.method))
    stop("The following arguments have to be nonnull: data, population.size, R.pop, R.samp, parameter, bootstrap.method.")

  if(!is.null(data) & !is.null(population.size) & !is.null(R.pop) & !is.null(R.samp) & !is.null(parameter) & !is.null(bootstrap.method)){
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

    if(!is.null(bootstrap.method) & bootstrap.method!="Booth.Butler.Hall" & bootstrap.method!="Bickel.Freedman" & bootstrap.method!="Chao.Lo.1985" & bootstrap.method!="Chao.Lo.1994" & bootstrap.method!="Sitter.BWO")
      stop("One of the following bootstrap methods should be specified: Bickel.Freedman (Bickel and Freedman, 1984), Chao.Lo.1985 (Chao and Lo, 1985), Chao.Lo.1994 (Chao and Lo, 1994), Sitter.BWO (Sitter, 1992, CJS) or Booth.Butler.Hall (Booth, Butler and Hall, 1994). The default is Booth.Butler.Hall.")

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
      if(bootstrap.method=="Bickel.Freedman" || bootstrap.method=="Booth.Butler.Hall" || bootstrap.method=="Sitter.BWO" || bootstrap.method=="Chao.Lo.1994" || bootstrap.method=="Chao.Lo.1985"){
        if(length(population.size)==L){

          boot.sample<- list()
          boot.statistic<- matrix(0, R.pop, R.samp)
          boot.parameter<- NULL

          for(a in 1:R.pop){

            PseudoPop.Bootstrap<- NULL
            n.bootstrap<- NULL
            PseudoPop.size<- NULL

            for(h in 1:L){
              n<- strata.size[h]
              N<- population.size[h]
              f<- n/N

              if(bootstrap.method == "Booth.Butler.Hall" || bootstrap.method == "Chao.Lo.1994" || bootstrap.method == "Chao.Lo.1985" || bootstrap.method == "Bickel.Freedman"){
                n.bootstrap<- c(n.bootstrap, n)
                k.boot<- floor(N/n)
              }

              if(bootstrap.method == "Chao.Lo.1985") Prob.PseudoPop<- (G.CL.1985(N, n)-G.CL.1985(n*(k.boot+1), n))/(G.CL.1985(n*k.boot, n)-G.CL.1985(n*(k.boot+1), n))
              if(bootstrap.method == "Bickel.Freedman") Prob.PseudoPop<- (1-(N-n*k.boot)/n)*(1-(N-n*k.boot)/(N-1))
              if(bootstrap.method == "Sitter.BWO") {
                k.boot<- floor(N*(n-(1-f))/n^2)
                a1<- (n*k.boot+1)/(n*(n-1)*(n*k.boot+n-1))
                a2<- (k.boot-1)/(n*(n*k.boot-1))
                Prob.PseudoPop<- 1-(((1-f)/(n*(n-1)))-a2)/(a1-a2)
              }

              PseudoPop.Fixed<- rep(1:n, k.boot)

              PseudoPop.Random<- NULL
              if(bootstrap.method == "Booth.Butler.Hall") PseudoPop.Random<- sample(1:n, N-n*k.boot, replace = FALSE)
              if(bootstrap.method == "Chao.Lo.1994") PseudoPop.Random<- sample(1:n, N-n*k.boot, replace = TRUE)
              if(bootstrap.method == "Chao.Lo.1985" || bootstrap.method == "Bickel.Freedman" || bootstrap.method == "Sitter.BWO"){
                Ind.PseudoPop.Random<- rbinom(1, 1, Prob.PseudoPop)
                if(Ind.PseudoPop.Random == 1) PseudoPop.Random<- NULL
                if(Ind.PseudoPop.Random == 0) PseudoPop.Random<- 1:n
              }

              if(bootstrap.method == "Sitter.BWO") n.bootstrap<- c(n.bootstrap, n - (1 - Ind.PseudoPop.Random))

              PseudoPop.index.h<- c(PseudoPop.Fixed, PseudoPop.Random)
              PseudoPop.size.h<- length(PseudoPop.index.h)
              PseudoPop.size<- c(PseudoPop.size, PseudoPop.size.h)
              PseudoPop.Bootstrap.h<- data.frame(data$study.variable[data$stratum==strata[h]][PseudoPop.index.h], PseudoPop.index.h, rep(strata[h], PseudoPop.size.h))
              PseudoPop.Bootstrap<- rbind(PseudoPop.Bootstrap, PseudoPop.Bootstrap.h)
            }

            boot.parameter[a]<- statistic(PseudoPop.Bootstrap[,1], rep(1, sum(PseudoPop.size)), q)
            if(parameter=="mean") boot.parameter[a]<- boot.parameter[a]/PseudoPop.size

            boot.sample.pp<- NULL
            for(b in 1:R.samp){

              boot.sample.index<- NULL
              boot.values<- NULL

              sample.index<- stratified.srs(PseudoPop.size, n.bootstrap, "SRSWOR")
              for(h in 1:L){
                boot.sample.h<- PseudoPop.Bootstrap[PseudoPop.Bootstrap[, 3]==strata[h], ][sample.index[,1][sample.index[,2]==h], ]
                boot.sample.index<- rbind(boot.sample.index, boot.sample.h)
              }
              boot.values<- rbind(boot.values, cbind(boot.sample.index[,1], rep(PseudoPop.size/n.bootstrap, n.bootstrap)) )

              boot.statistic[a, b]<- statistic(boot.values[ ,1], boot.values[ ,2], q)
              if(parameter=="mean") boot.statistic[a, b]<- boot.statistic[a, b]/PseudoPop.size

              colnames(boot.sample.index)<- c("bootstrap.sample", "ID", "stratum")
              boot.sample.pp[[b]]<- as.data.frame(boot.sample.index)
            }

            boot.sample[[a]]<- boot.sample.pp
          }

          boot.mean<- mean(boot.statistic)
          boot.var<- mean(apply(boot.statistic, 1, var))

          return<- list(
            boot.statistic=boot.statistic,
            boot.parameter=boot.parameter,
            boot.var=boot.var,
            boot.mean=boot.mean,
            boot.sample=boot.sample,
            parameter=parameter,
            number.bootstrap.pop.replicates=R.pop,
            number.bootstrap.samp.replicates=R.samp
          )
        }
      }
    }
  }

  return
}

