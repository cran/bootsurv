#' Bootstrap methods for two-stage sampling designs
#'
#' The function `boot.twostage` applies one of the following bootstrap methods on complete (full response) survey data selected under stratified two-stage cluster sampling SRSWOR/SRSWOR: Rao and Wu (1988), Rao, Wu and Yue (1992), the modified version of Sitter (1992, CJS) (see Chen, Haziza and Mashreghi, 2022), Funaoka, Saigo, Sitter and Toida (2006), Chauvet (2007) or Preston (2009).
#' This function also applies the method of Rao, Wu and Yue (1992) on complete survey data selected under stratified two-stage cluster sampling IPPSWOR/SRSWOR or the method of Chauvet (2007) on complete survey data selected under stratified two-stage cluster sampling CPS/SRSWOR.
#'
#' @param data A vector, matrix or data frame. The column of study variable has to be a numeric column named `study.variable` and a column identifying clusters named `cluster` has to be included. If the population is stratified, a column identifying strata named `stratum` has to be included. If an IPPS design is applied on the first stage a column of first stage inclusion probability named `Pi1` has to be included.
#' @param no.cluster A vector of the number of clusters within strata.
#' @param cluster.size The number of elements within the selected clusters within each stratum. The length of this vector must be the same as the number of all selected clusters in all strata.
#' @param R The number of bootstrap replicates. For the Chauvet (2007) method, `R` is a vector with two values: `(R.pop, R.samp)` representing the number of pseudo-populations and the number of bootstrap samples drawn from each pseudo-population.
#' @param parameter One of the following population parameters can be applied: `"total"` (population total), `"mean"` (population mean), `"quartile.25"` (population 1st quartile), `"quartile.50"` or `"median"` (population median) or `"quartile.75"` (population 3rd quartile). If the parameter of interest is the population mean or total, the HT-estimator is applied. If the parameter of interest is a population quartile, the estimator in Sarndal, Swensson, and Wretman (1992, Chapter 5) is applied. The default is the population total.
#' @param survey.design It can be either `"IPPS"` only if the method of Rao, Wu and Yue (1992) is applied or `"CPS"` only if the method of Chauvet (2007) is applied or `"SRSWOR"`. The default is `"SRSWOR"`.
#' @param bootstrap.method One of the following bootstrap methods can be applied in the case of statratified SRS/SRS: `"Rao.Wu"` (Rao and Wu, 1988), `"Rao.Wu.Yue"` (Rao, Wu and Yue, 1992), `"Modified.Sitter"` (the modified version of Sitter 1992 discussed in Chen, Haziza and Mashreghi, 2022), `"Funaoka.etal"` (Funaoka, Saigo, Sitter and Toida, 2006), `"Chauvet"` (Chauvet, 2007) or `"Preston"` (Preston, 2009).
#' @param population.size A vector of stratum population sizes.
#' @param boot.sample.size A vector of bootstrap sample sizes within strata. The bootstrap sample size is required only for the method of Rao, Wu and Yue (1988). If it is not specified, the bootstrap sample size will be `nh-1` within each stratum, where `nh` is the original sample size within stratum `h`.
#'
#' @return
#'
#' `boot.statistic` A vector of bootstrap statistics of size R.
#'
#' `boot.var` The bootstrap variance estimator of the estimator of parameter of interest.
#'
#' `boot.mean` The average of the bootstrap estimator of the parameter of interest.
#'
#' `boot.sample` A list of results for each iteration. That includes a column of original sample values, a column of cluster identifier and a column of stratum identifier. More information is availble depending on the bootstrap method.
#'
#' @export
#'
#' @importFrom stats aggregate rbinom rexp rnorm runif var
#'
#' @references
#'
#' Chauvet, G. (2007). Méthodes de bootstrap en population finie. PhD thesis, École Nationale de Statistique et Analyse de l’Information, Bruz, France.
#'
#' Chen, S., Haziza, D. and Mashreghi, Z., (2022). A Comparison of Existing Bootstrap Algorithms for Multi-Stage Sampling Designs. Stats, 5(2), pp.521-537.
#'
#' Funaoka, F., Saigo, H., Sitter, R.R., Toida, T. (2006). Bernoulli bootstrap for stratified multistage sampling. Survey Methodology, 32, 151–156.
#'
#' Rao, J.N.K., Wu, C.F.J. (1998). Resampling inference with complex survey data. Journal of the American Statistical Association, 83, 231–241.
#'
#' Rao, J.N.K., Wu, C.F.J., Yue, K. (1992). Some recent work on resampling methods for complex surveys. Survey Methodology, 18, 209–217.
#'
#' Särndal, C.-E., Swensson, B. and Wretman, J. (1992). Model-Assisted Survey Sampling. NewYork: Springer.
#'
#' Sitter, R.R. (1992). Comparing three bootstrap methods for survey data. The Canadian Journal of Statistics, 20, 135–154.
#'
#' Preston, J. (2009). Rescaled bootstrap for stratified multistage sampling. Survey Methodology, 35, 227–234.
#'
#'
#' @examples
#'
#' R<- 20
#'
#' data(data_samp_clust)
#' data(data_pop_clust)
#' no_cluster<- 200
#' cluster_size<- table(data_pop_clust$cluster)[unique(data_samp_clust$cluster)]
#'
#' # The first stage sampling fraction is about 20% and the overall second stage sampling is about 15%.
#' # data_samp_clust is a sample taken from data_pop_clust available in the package.
#'
#' boot.RWY<- boot.twostage(data_samp_clust, no_cluster, cluster_size, R)
#' boot.RWY$boot.var
#'
#' boot.Pr<- boot.twostage(data_samp_clust, no_cluster, cluster_size, R, bootstrap.method="Preston")
#' boot.Pr$boot.var
#'
#' boot.RWY.med<- boot.twostage(data_samp_clust, no_cluster, cluster_size, R, parameter="median")
#' boot.RWY.med$boot.var
#' boot.RWY.med$boot.sample[[5]]
#'
#' boot.Ch<- boot.twostage(data_samp_clust, no_cluster, cluster_size, R=c(5, 10),
#'            bootstrap.method="Chauvet")
#' boot.Ch$boot.mean
#'
#' data(data_samp_stclust)
#' data(data_pop_stclust)
#' # The first stage sampling fraction is about 20% and the overall second stage sampling is about 15%.
#' # data_samp_stclust is a sample taken from data_pop_stclust available in the package.
#'
#' no_cluster_stclust<- c(100, 125, 65)
#' cluster_size_pop_st<- aggregate(data_pop_stclust$cluster,
#'  by=list(data_pop_stclust$stratum), table)[[2]]
#' L<- length(unique(data_samp_stclust$stratum))
#' cluster_size_st<- NULL
#' for(h in 1:L) cluster_size_st<- c(cluster_size_st,
#'  cluster_size_pop_st[[h]][unique(data_samp_stclust$cluster[data_samp_stclust$stratum==h])])
#'
#' boot.RWY.st<- boot.twostage(data_samp_stclust, no_cluster_stclust, cluster_size_st, R)
#' boot.RWY.st$boot.statistic
#'
#'
#'
boot.twostage<- function(data, no.cluster, cluster.size, R, parameter="total", bootstrap.method="Rao.Wu.Yue", survey.design="SRSWOR", population.size=NULL, boot.sample.size=NULL){

  cps.poisson<- function(N, U.pi, n=NULL){
    epsilon<- runif(N)
    s.index<-1*(epsilon< U.pi)
    if(n!="NULL"){
      while(sum(s.index)!=n){
        epsilon<- runif(N)
        s.index<- 1*(epsilon < U.pi)
      }
    }
    return(s.index)
  }
  stratified.srs<- function(Nh, nh, design){
    L<- length(Nh)
    sample.index<- NULL
    for(h in 1:L){
      if(design=="SRSWR") sample.index<- rbind(sample.index, cbind(sort(sample(1:Nh[h], nh[h], replace = TRUE)), rep(h,nh[h])) )
      if(design=="SRSWOR") sample.index<- rbind(sample.index, cbind(sort(sample(1:Nh[h], nh[h], replace = FALSE)), rep(h,nh[h])) )
    }
    sample.index
  }

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

  if(is.null(data) || is.null(no.cluster) || is.null(cluster.size) || is.null(R) || is.null(parameter) || is.null(bootstrap.method) || is.null(survey.design))
    stop("The following arguments have to be nonnull: data, no.cluster, cluster.size, R, parameter, bootstrap.method, survey.design.")

  if(!is.null(R) & bootstrap.method=="Chauvet" & length(R)!=2)
    stop("For the method of Chauvet (2007), R is a vector with two values: (R.pop, R.samp) representing the number of pseudo-populations and the number of bootstrap samples drawn from each pseudo-population.")

  if(!is.null(data) & !is.null(no.cluster) & !is.null(cluster.size) & !is.null(R) & !is.null(parameter) & !is.null(bootstrap.method) & !is.null(survey.design)){
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

    if(!is.null(bootstrap.method) & bootstrap.method!="Rao.Wu" & bootstrap.method!="Rao.Wu.Yue" & bootstrap.method!="Modified.Sitter" & bootstrap.method!="Funaoka.etal" & bootstrap.method!="Chauvet" & bootstrap.method!="Preston")
      stop("One of the following bootstrap methods should be specified: 'Rao.Wu', 'Rao.Wu.Yue', 'Modified.Sitter', 'Funaoka.etal', 'Chauvet' or 'Preston'. The default is 'Rao.Wu.Yue'")

    if(!is.null(survey.design) & survey.design!="SRSWOR" & survey.design!="IPPS" & survey.design!="CPS")
      stop("This function works for two-stage cluster sampling with SRSWOR in both stages. It can also applies the method of Rao, Wu and Yue (1992) if an IPPS design and SRSWOR design are respectively applied on the first and second stages and the method of Chauvet (2007) if conditional Poisson sampling (CPS) is applied on the first stage and SRSWOR is applied on the second stage.")

    if(is.null(data$study.variable))
      stop("The data set must contain a culumn with the values of the variable of interest called 'study.variable', a culumn that identifies the clusters called 'cluster' and a culumn that identifies the strata called 'stratum' if the survey design is STSRSWOR.")

    if(is.null(data$stratum)) data$stratum<- 1
    if(is.null(data$cluster)) data$cluster<- 1

    if(bootstrap.method=="Rao.Wu"){
      if(is.null(population.size))
        stop("In the bootstrap method of Rao and Wu (1988), the total population size is required. The length of the arguments 'no.cluster' and 'population.size' must be the same. Each element of 'no.cluster' and 'population.size' show the number of clusters and the number of units in each subpopulation (or stratum), respectively.")
      if(!is.null(population.size) & length(no.cluster)!=length(population.size))
        stop("In the bootstrap method of Rao and Wu (1988), the length of the arguments 'no.cluster' and 'population.size' must be the same. Each element of 'no.cluster' and 'population.size' show the number of clusters and the number of units in each subpopulation (or stratum), respectively.")
    }

    if(is.null(population.size) & parameter=="mean")
      stop("To compute the bootstrap statistics for the case of the population mean, the population size is required.")
    else{
    data<- data[with(data, order(stratum, cluster)),  ]

    strata<- unique(data$stratum)
    strata.size<- table(data$stratum)
    L<- length(strata)

    if(!is.numeric(cluster.size) || length(cluster.size)!=sum(table(data$stratum,data$cluster)!=0))
      stop("'cluster.size' has to be a numeric vector. Each element of this vector identifies the number of elements within each selected cluster. The length of this vector must be the same as the total number of selected clusters presented in 'data'.")


    if(!is.null(parameter) || parameter=="total" || parameter=="mean" || parameter=="median" || parameter=="quartile.25" || parameter=="quartile.50" || parameter=="quartile.75"){
      if(is.numeric(cluster.size) & length(cluster.size)==sum(table(data$stratum,data$cluster)!=0)){

        M<- list()
        data$fpc2<- rep(0, dim(data)[1])
        data$ID<- rep(0, dim(data)[1])
        fpc1<- NULL
        if(!is.null(population.size)) population.mean.size<- population.size/no.cluster

        j<- 1
        pseudo.pop.PPBSitter<- NULL
        for(h in 1:L){
          cluster.h <- data$cluster[data$stratum==strata[h]]
          sampled.cluster.h<- unique(cluster.h)
          Nh<- no.cluster[h]

          no.sampled.cluster.h<- length(sampled.cluster.h)
          M[[h]]<- cluster.size[j:(j+no.sampled.cluster.h-1)]
          Mh<- as.vector(M[[h]])
          j<- j+no.sampled.cluster.h

          mh<- table(cluster.h)
          nh<- length(sampled.cluster.h)

          fpc1<- c(fpc1, rep(nh/Nh, sum(mh)))

          for(i in 1:nh){
            data$ID[data$stratum==strata[h] & data$cluster==sampled.cluster.h[i]]<- 1:mh[i]
            data$fpc2[data$stratum==strata[h] & data$cluster==sampled.cluster.h[i]]<- rep(mh[i]/Mh[i], mh[i])
          }

          if(bootstrap.method=="Modified.Sitter"){
            data.h<- as.data.frame(data[data$stratum==strata[h], ])
            k1.PPBSitter.Mod<- round(Nh/nh)
            k2.PPBSitter.Mod<- round(Mh/mh)

            nh.PPBSitter.Mod<- round(Nh/(1+(1-nh/Nh)*(k1.PPBSitter.Mod*nh-1)/(nh*k1.PPBSitter.Mod*(nh-1)/Nh)))
            mh.PPBSitter.Mod<- round(Mh/(1+((nh/Nh)*(Mh/mh-1)*nh.PPBSitter.Mod*(k2.PPBSitter.Mod*mh-1)/(nh*k2.PPBSitter.Mod*(mh-1))) )  )

            ## 1st stage pseudo-population:
            pseudo.ssu.PPBSitter<- NULL
            for(j in 1:nh){
              r=0
              repeat{
                r=r+1
                pseudo.ssu.PPBSitter<- rbind(pseudo.ssu.PPBSitter, data.frame(data.h[data.h$cluster==sampled.cluster.h[j], ]) )
                if(r==k2.PPBSitter.Mod[j]){break}
              }
            }

            ## 2nd stage pseudo-population:
            pseudo.psu.PPBSitter<- NULL
            new.cluster<- 0
            for(j in 1:nh){
              t<- 0
              repeat{
                new.cluster=new.cluster+1
                t<- t+1
                pseudo.psu.PPBSitter<- rbind(pseudo.psu.PPBSitter, data.frame(pseudo.ssu.PPBSitter[pseudo.ssu.PPBSitter$cluster==sampled.cluster.h[j], ], new.cluster) )
                if(t==k1.PPBSitter.Mod){break}
              }
            }

            ## Final pseudo-population
            fpc1.Sitter<- rep(nh.PPBSitter.Mod/(nh*k1.PPBSitter.Mod), sum(mh*k1.PPBSitter.Mod*k2.PPBSitter.Mod))
            fpc2.Sitter<- rep(mh.PPBSitter.Mod/(mh*k2.PPBSitter.Mod), mh*k1.PPBSitter.Mod*k2.PPBSitter.Mod)
            pseudo.psu.PPBSitter$fpc1.Sitter<- fpc1.Sitter
            pseudo.psu.PPBSitter$fpc2.Sitter<- fpc2.Sitter
            pseudo.psu.PPBSitter$weight.Sitter<- 1/(fpc1.Sitter*fpc2.Sitter)

            pseudo.pop.PPBSitter<- rbind(pseudo.pop.PPBSitter, pseudo.psu.PPBSitter)
          }
        }

        data$fpc1<- fpc1
        data$Pi2<- data$fpc2

        if(bootstrap.method=="Rao.Wu" & !is.null(population.size) & length(no.cluster)==length(population.size)){
          if(survey.design=="SRSWOR"){
            data$Pi1<- data$fpc1
            data$weight<- 1/(data$Pi1*data$Pi2)

            boot.sample<- list()
            boot.statistic<- NULL

            for (b in 1:R){

              boot.sample.index<- NULL
              boot.values<- NULL

              for(h in 1:L){

                Nh<- no.cluster[h]
                Mh<- as.vector(M[[h]])
                data.h<- as.data.frame(data[data$stratum==strata[h], ])
                cluster.h<- data.h$cluster
                sampled.cluster.h<- unique(cluster.h)
                mh<- table(cluster.h)
                nh<- length(sampled.cluster.h)

                population.size.h<- population.size[h]

                ## 1st stage resampling:
                sample.psu.index.RW<- sort(sample(1:nh, nh, replace=TRUE))
                sample.psu.RW<- sampled.cluster.h[sample.psu.index.RW]
                sample1.RW<- NULL
                new.cluster<- 0
                for(i in 1:nh){
                  new.cluster<- new.cluster+1
                  sample1.RW<- rbind(sample1.RW, data.frame(data.h[cluster.h==sample.psu.RW[i], ], new.cluster))
                }
                mh.RW<- mh[sample.psu.index.RW]
                Mh.RW<- Mh[sample.psu.index.RW]

                ## 2nd stage resampling:
                sample2.index<- stratified.srs(mh.RW, mh.RW, "SRSWR")
                sample2.RW<- NULL
                for(i in 1:nh) sample2.RW<- rbind(sample2.RW, data.frame(sample1.RW[sample1.RW$new.cluster==i, ][sample2.index[, 1][sample2.index[, 2]==i], ]))

                ## Rescaling the selected y-values

                lambda1<- sqrt(nh*(1-nh/Nh)/(nh-1))
                lambda2i<- sqrt(sample2.RW$fpc1*(1-sample2.RW$fpc2)*rep(mh.RW, mh.RW)/(rep(mh.RW, mh.RW)-1))

                #sample2.RW$study.variable<- as.numeric(sample2.RW$study.variable)
                emean<- sum(data.h$study.variable*data$weight)/population.size.h
                etotal.cluster.RW<- Mh.RW*aggregate(sample2.RW$study.variable, list(sample2.RW$new.cluster), FUN=mean)[,2]

                rescaled.y.RW<- emean + lambda1*(rep(etotal.cluster.RW, mh.RW)/population.mean.size[h]-emean) + lambda2i*(sample2.RW$study.variable*rep(Mh.RW, mh.RW)- rep(etotal.cluster.RW, mh.RW))/population.mean.size[h]

                boot.sample.index<- rbind(boot.sample.index, data.frame(rescaled.y.RW, sample2.RW$ID, sample2.RW$cluster, sample2.RW$stratum ))
                boot.values<- rbind(boot.values, cbind(rescaled.y.RW, rep(population.size.h/(nh*mh.RW), mh.RW)) )
              }

              boot.statistic[b]<- statistic(as.numeric(boot.values[ ,1]), as.numeric(boot.values[ ,2]), q)
              if(parameter=="mean" & !is.null(population.size)) boot.statistic[b]<- boot.statistic[b]/sum(population.size)
              colnames(boot.sample.index)<- c("bootrstap.values", "ID", "cluster", "stratum")
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
          else stop("This function can apply one of the following methods only if SRSWOR is applied on the first stage: Rao and Wu (1988), Rao, Wu and Yue (1992), the modified version of Sitter (1992, CJS) (see Chen, Haziza and Mashreghi, 2022), Funaoka, Saigo, Sitter and Toida (2006), Chauvet (2007) and Preston (2009).")
        }
        if(bootstrap.method=="Modified.Sitter" || bootstrap.method=="Funaoka.etal"|| bootstrap.method=="Preston"){
          if(survey.design=="SRSWOR"){
            data$Pi1<- data$fpc1
            data$weight<- 1/(data$Pi1*data$Pi2)

            boot.sample<- list()
            boot.statistic<- NULL

            for (b in 1:R){

              boot.sample.index<- NULL
              boot.values<- NULL

              for(h in 1:L){

                Nh<- no.cluster[h]
                Mh<- as.vector(M[[h]])
                data.h<- as.data.frame(data[data$stratum==strata[h], ])
                cluster.h<- data.h$cluster
                sampled.cluster.h<- unique(cluster.h)
                mh<- table(cluster.h)
                nh<- length(sampled.cluster.h)

                if(bootstrap.method=="Modified.Sitter"){
                  k1.PPBSitter.Mod<- round(Nh/nh)
                  k2.PPBSitter.Mod<- round(Mh/mh)

                  nh.PPBSitter.Mod<- round(Nh/(1+(1-nh/Nh)*(k1.PPBSitter.Mod*nh-1)/(nh*k1.PPBSitter.Mod*(nh-1)/Nh)))
                  mh.PPBSitter.Mod<- round(Mh/(1+((nh/Nh)*(Mh/mh-1)*nh.PPBSitter.Mod*(k2.PPBSitter.Mod*mh-1)/(nh*k2.PPBSitter.Mod*(mh-1))) )  )

                  pseudo.pop.PPBSitter.h<- pseudo.pop.PPBSitter[pseudo.pop.PPBSitter$stratum==strata[h], ]

                  ## 1st stage resampling:
                  sample.psu.PBBSitter<- sort(sample(1:(k1.PPBSitter.Mod*nh), nh.PPBSitter.Mod, replace = FALSE))
                  sample1.PPBSitter<- pseudo.pop.PPBSitter.h[pseudo.pop.PPBSitter.h$new.cluster %in% sample.psu.PBBSitter,]

                  ## 2nd stage resampling:
                  Mh.PPBSitter<- table(as.numeric(pseudo.pop.PPBSitter.h$new.cluster))[sample.psu.PBBSitter]
                  mh.cluster.Mod<- rep(mh.PPBSitter.Mod, each=k1.PPBSitter.Mod)[sample.psu.PBBSitter]

                  sample2.index<- stratified.srs(Mh.PPBSitter, mh.cluster.Mod, "SRSWOR")
                  sample2.PPBSitter<- NULL
                  for(i in 1:nh.PPBSitter.Mod) sample2.PPBSitter<- rbind(sample2.PPBSitter, data.frame(sample1.PPBSitter[sample1.PPBSitter$new.cluster==unique(sample1.PPBSitter$new.cluster)[i], ][sample2.index[, 1][sample2.index[, 2]==i], ]) )

                  boot.sample.index<- rbind(boot.sample.index, cbind(sample2.PPBSitter$study.variable, sample2.PPBSitter$ID, sample2.PPBSitter$cluster, sample2.PPBSitter$stratum ))
                  boot.values<- rbind(boot.values, cbind(sample2.PPBSitter$study.variable, sample2.PPBSitter$weight.Sitter) )
                }

                if(bootstrap.method=="Funaoka.etal"){

                  nh.BernFunaoka<- nh-1
                  Prob1.BernFunaoka<- 1-(1-nh/Nh)/(2*(1-1/nh))

                  mh.BernFunaoka<- mh-1
                  Prob2.BernFunaoka<- 1-data.h$fpc1*(1-data.h$fpc2)/(2*(1-Prob1.BernFunaoka)^{-1}*(1-rep(mh^{-1}, mh)) )

                  ##Taking a sample of size n-1 WR from psu's
                  sampled.psu.BernFunaoka<- sort(sample(sampled.cluster.h, nh.BernFunaoka, replace=TRUE))

                  ##The psu that will be kept with Prob1 and replaced with 1-Prob1:
                  Bern.clusters<- rbinom(nh, 1, prob=Prob1.BernFunaoka)
                  Kept.clusters<- sampled.cluster.h[Bern.clusters==1]
                  Replaced.clusters<- sampled.cluster.h[Bern.clusters==0]
                  No.Replaced.clusters<- sum(Bern.clusters==0)

                  ##Replaced clusters
                  Replaced.cluster.BernFunaoka<- sort(sample(sampled.psu.BernFunaoka, No.Replaced.clusters, replace=TRUE))
                  sample2.BernFunaoka<- NULL
                  new.cluster<- 0
                  if(No.Replaced.clusters>0){
                    for(i in 1:No.Replaced.clusters){
                      new.cluster<- new.cluster+1
                      sample2.BernFunaoka<- rbind(sample2.BernFunaoka, data.frame(data.h[data.h$cluster==Replaced.cluster.BernFunaoka[i], ], new.cluster))
                    }
                  }

                  ##Kept clusters - Go to Step II
                  No.Kept.clusters<- nh-No.Replaced.clusters

                  if(No.Kept.clusters>0){
                    for(i in 1:No.Kept.clusters){
                      new.cluster<- new.cluster+1
                      sampled.ssu.BernFunaoka<- sample(mh[names(mh)==Kept.clusters[i]], mh.BernFunaoka[names(mh.BernFunaoka)==Kept.clusters[i]], replace = TRUE)
                      Bern.ssu<- rbinom(mh[names(mh)==Kept.clusters[i]], 1, prob=Prob2.BernFunaoka[data.h$cluster==Kept.clusters[i]])
                      No.Replaced.ssu<- sum(Bern.ssu==0)
                      sample.ssu.BernFunaoka<- data.h[data.h$cluster==Kept.clusters[i], ]
                      if(No.Replaced.ssu>0)  sample.ssu.BernFunaoka[Bern.ssu==0, ]<- sample.ssu.BernFunaoka[sample(sampled.ssu.BernFunaoka, No.Replaced.ssu, replace = TRUE), ]
                      sample2.BernFunaoka<- rbind(sample2.BernFunaoka, cbind(data.frame(sample.ssu.BernFunaoka), new.cluster))
                    }
                  }

                  boot.sample.index<- rbind(boot.sample.index, cbind(sample2.BernFunaoka$study.variable, sample2.BernFunaoka$ID, sample2.BernFunaoka$cluster, sample2.BernFunaoka$stratum ))
                  boot.values<- rbind(boot.values, cbind(sample2.BernFunaoka$study.variable, sample2.BernFunaoka$weight) )
                }

                if(bootstrap.method=="Preston"){
                  nh.Preston<- floor(nh/2)
                  mh.Preston<- floor(mh/2)
                  lambda1<- sqrt(nh.Preston*(1-nh/Nh)/(nh-nh.Preston))
                  lambda2i<- sqrt(mh.Preston*(nh/Nh)*(1-mh/Mh)/(mh-mh.Preston))

                  sample.psu.index.Preston<- sort(sample(1:nh, nh.Preston, replace = FALSE))
                  sample.psu.ind.Preston<- rep(0, nh)
                  sample.psu.ind.Preston[sample.psu.index.Preston]<- 1
                  data.h$psu.indicator.Preston<- rep(sample.psu.ind.Preston, mh)

                  data.h$w1.Preston<- (1+lambda1*(nh* data.h$psu.indicator.Preston/nh.Preston-1))/data.h$Pi1
                  data.h$weight.Preston<- data.h$weight

                  data.h$ssu.indicator.Preston<- rep(0, sum(mh))
                  for(i in 1:nh){
                    sample.ssu.index.Preston<- sort(sample(1:mh[i], mh.Preston[i], replace = FALSE))
                    sample.ssu.ind.Preston<- rep(0, mh[i])
                    sample.ssu.ind.Preston[sample.ssu.index.Preston]<- 1
                    data.h$ssu.indicator.Preston[data.h$cluster==sampled.cluster.h[i]]<- sample.ssu.ind.Preston
                    data.h$weight.Preston[data.h$cluster==sampled.cluster.h[i]]<- (1 + lambda1*(nh*sample.psu.ind.Preston[i]/nh.Preston-1) + lambda2i[i] * sqrt(nh/nh.Preston) * sample.psu.ind.Preston[i] * (mh[i]*sample.ssu.ind.Preston/mh.Preston[i]-1) ) * data.h$weight[data.h$cluster==sampled.cluster.h[i]]#  /data.h$w1.Preston[data.h$cluster==sampled.cluster.h[i]]
                  }

                  boot.sample.index<- rbind(boot.sample.index, cbind(data.h$stratum, data.h$cluster, data.h$ID, data.h$psu.indicator.Preston, data.h$ssu.indicator.Preston))
                  boot.values<- rbind(boot.values, cbind(data.h$study.variable, data.h$weight.Preston) )
                }
              }

              boot.statistic[b]<- statistic(as.numeric(boot.values[ ,1]), as.numeric(boot.values[ ,2]), q)
              if(parameter=="mean" & !is.null(population.size)) boot.statistic[b]<- boot.statistic[b]/sum(population.size)

              if(bootstrap.method=="Preston") colnames(boot.sample.index)<- c("stratum", "cluster", "ID", "indicator.selected.psu", "indicator.selected.ssu")
              if(bootstrap.method=="Modified.Sitter" || bootstrap.method=="Funaoka.etal") colnames(boot.sample.index)<- c("bootrstap.values", "ID", "cluster", "stratum")
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
          else stop("This function can apply one of the following methods only if SRSWOR is applied on the first stage: Rao and Wu (1988), Rao, Wu and Yue (1992), the modified version of Sitter (1992, CJS) (see Chen, Haziza and Mashreghi, 2022), Funaoka, Saigo, Sitter and Toida (2006), Chauvet (2007) and Preston (2009).")
        }
        if(bootstrap.method=="Rao.Wu.Yue"){
          if(survey.design=="SRSWOR" || survey.design=="IPPS"){
            if(is.null(data$Pi1) || survey.design=="SRSWOR") data$Pi1<- data$fpc1
            data$weight<- 1/(data$Pi1*data$Pi2)

            boot.sample<- list()
            boot.statistic<- NULL

            for (b in 1:R){

              boot.sample.index<- NULL
              boot.values<- NULL

              for(h in 1:L){

                Nh<- no.cluster[h]
                Mh<- as.vector(M[[h]])
                data.h<- data[data$stratum==strata[h], ]
                cluster.h<- data.h$cluster
                sampled.cluster.h<- unique(cluster.h)
                mh<- table(cluster.h)
                nh<- length(sampled.cluster.h)

                if(is.null(boot.sample.size)) nh.RWY<- nh-1
                if(!is.null(boot.sample.size)) nh.RWY<- boot.sample.size[h]
                sample.psu.index.RWY<- sample(1:nh, nh.RWY, replace = TRUE)
                No.Rep.RWY<- tabulate(sample.psu.index.RWY, nh)
                data.h$No.Repetition.RWY<- rep(No.Rep.RWY, mh)

                RescFac.RWY<- 1-sqrt(nh.RWY/(nh-1)) + sqrt(nh.RWY/(nh-1))*nh*No.Rep.RWY/nh.RWY
                data.h$weight.RWY<- rep(RescFac.RWY, mh)*data.h$weight

                boot.sample.index<- rbind(boot.sample.index, data.frame(data.h$study.variable, data.h$ID, data.h$cluster, data.h$stratum, data.h$weight.RWY, data.h$No.Repetition.RWY))
                boot.values<- rbind(boot.values, cbind(data.h$study.variable, data.h$weight.RWY) )
              }

              boot.statistic[b]<- statistic(as.numeric(boot.values[ ,1]), as.numeric(boot.values[ ,2]), q)
              if(parameter=="mean" & !is.null(population.size)) boot.statistic[b]<- boot.statistic[b]/sum(population.size)

              colnames(boot.sample.index)<- c("study.variable", "ID", "cluster", "stratum", "bootstrap.weights", "No.selected.psu")
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
          else stop("This function can apply the method of Rao, Wu, Yue (1992) if either SRSWOR or IPPS is applied on the first stage.")
        }
        if(bootstrap.method=="Chauvet"){
          if(survey.design=="SRSWOR" || survey.design=="CPS"){
            if(is.null(data$Pi1) || survey.design=="SRSWOR") data$Pi1<- data$fpc1
            data$weight<- 1/(data$Pi1*data$Pi2)

            boot.sample<- list()
            boot.statistic<- matrix(0, R[1], R[2])

            for (a in 1:R[1]){


              N.PPBChauvet<- NULL
              n<- NULL
              pseudo.pop.PPBChauvet<- NULL
              U.pi.Chauvet<- list()
              m.Chauvet<- list()
              cluster.size.PPBChauvet<- list()

              for(h in 1:L){

                Nh<- no.cluster[h]
                Mh<- as.vector(M[[h]])
                data.h<- data[data$stratum==strata[h], ]
                cluster.h<- data.h$cluster
                sampled.cluster.h<- unique(cluster.h)
                mh<- table(cluster.h)
                nh<- length(sampled.cluster.h)
                n<- c(n, nh)

                ### Creating the pseudo-population
                ## Repeating the elements of the selected ssu
                pseudo.ssu.PPBChauvet<- data.h[rep(1:sum(mh), round(1/data.h$Pi2)), ]
                rownames(pseudo.ssu.PPBChauvet)<- 1:dim(pseudo.ssu.PPBChauvet)[1]

                ## Repeating the selected clusters
                Pi1.Chauvet<- aggregate(pseudo.ssu.PPBChauvet$Pi1, by=list(pseudo.ssu.PPBChauvet$cluster), FUN=unique)[,2]
                rept.psu.add<- rep(0, nh)
                if(survey.design=="SRSWOR")  rept.psu.add[sample(1:nh, Nh-nh*floor(Nh/nh), replace = FALSE)]<- 1
                if(survey.design=="CPS") rept.psu.add<- cps.poisson(nh, 1/Pi1.Chauvet-floor(1/Pi1.Chauvet), NULL)
                rept.psu<- floor(1/Pi1.Chauvet) + rept.psu.add

                pseudo.psu.PPBChauvet<- NULL
                new.cluster<- 0
                for(i in 1:nh){
                  r=0
                  repeat{
                    r=r+1
                    new.cluster<- new.cluster+1
                    pseudo.psu.PPBChauvet<- rbind(pseudo.psu.PPBChauvet, cbind(data.frame(pseudo.ssu.PPBChauvet[pseudo.ssu.PPBChauvet$cluster==sampled.cluster.h[i], ]), new.cluster))
                    if(r==rept.psu[i]){break}
                  }
                }

                pseudo.pop.PPBChauvet.h<- table(as.numeric(pseudo.psu.PPBChauvet$new.cluster))
                cluster.size.PPBChauvet[[h]]<- pseudo.pop.PPBChauvet.h
                mh.Chauvet<- rep(mh, rept.psu)
                m.Chauvet[[h]]<- as.vector(mh.Chauvet)
                fpc2.Chauvet<- mh.Chauvet/pseudo.pop.PPBChauvet.h

                if(survey.design=="CPS"){
                  Uh.pi.Chauvet<- nh*pseudo.pop.PPBChauvet.h/sum(pseudo.pop.PPBChauvet.h)
                  while(sum(Uh.pi.Chauvet>1) >= 1){
                    for(j in 1:Nh.PPBChauvet){
                      if(Uh.pi.Chauvet[j]>1) Uh.pi.Chauvet[j]<-1
                      if(Uh.pi.Chauvet[j]<1) Uh.pi.Chauvet[j]<-(nh-sum(Uh.pi.Chauvet>=1))*pseudo.pop.PPBChauvet.h[j]/sum(pseudo.pop.PPBChauvet.h[Uh.pi.Chauvet<1])
                    }
                  }
                  U.pi.Chauvet[[h]]<- as.vector(Uh.pi.Chauvet)
                  pseudo.psu.PPBChauvet$Pi1<- rep(Uh.pi.Chauvet, pseudo.pop.PPBChauvet.h)
                  pseudo.psu.PPBChauvet$weight<- 1/(pseudo.psu.PPBChauvet$Pi1*pseudo.psu.PPBChauvet$Pi2)
                }

                pseudo.pop.PPBChauvet<- rbind(pseudo.pop.PPBChauvet, pseudo.psu.PPBChauvet)#, stratum=rep(strata[h], PseudoPop.size.h)) )

                N.PPBChauvet<- c(N.PPBChauvet, sum(rept.psu))
              }

                boot.sample.pp<- NULL
                sample.PPBChauvet<- NULL

                for (b in 1:R[2]){

                  boot.sample.index<- NULL
                  boot.values<- NULL

                  for(h in 1:L){

                    data.h<- data[data$stratum==strata[h], ]

                    n<- c(n, nh)
                    nh<- n[h]
                    Nh.PPBChauvet<- N.PPBChauvet[h]
                    mh.Chauvet<- m.Chauvet[[h]]
                    cluster.size.PPBChauvet.h<- cluster.size.PPBChauvet[[h]]

                    pseudo.pop.PPBChauvet.h<- pseudo.pop.PPBChauvet[pseudo.pop.PPBChauvet$stratum==strata[h], ]

                    #Inclusion Probability pik Proportional to M.Chauvet if the design is IPPS ###
                    if(survey.design=="CPS"){
                      Uh.pi.Chauvet<- U.pi.Chauvet[[h]]
                      sI1.PPBChauvet<- sort((1:Nh.PPBChauvet)[cps.poisson(Nh.PPBChauvet, Uh.pi.Chauvet, nh)==1])
                      # (1:Nh.PPBChauvet)[UPrandomsystematic(Uh.pi.Chauvet)==1]
                      #if(survey.design=="Poisson") sI1.PPBChauvet<- (1:Nh)[cps.poisson(Nh.PPBChauvet, Uh.pi.Chauvet, NULL)==1]
                      sample1.PPBChauvet<- pseudo.pop.PPBChauvet.h[pseudo.pop.PPBChauvet.h$new.cluster %in% sI1.PPBChauvet, ]
                      sampled.cluster.Chauvet<- unique(sample1.PPBChauvet$new.cluster)
                    }

                    if(survey.design=="SRSWOR"){
                      sI1.PPBChauvet<- sort(sample(1:Nh.PPBChauvet, nh, replace = FALSE))
                      sample1.PPBChauvet<- pseudo.pop.PPBChauvet.h[pseudo.pop.PPBChauvet.h$new.cluster %in% sI1.PPBChauvet, ]
                      sampled.cluster.Chauvet<- unique(sample1.PPBChauvet$new.cluster)
                    }

                    ## draw 2nd stage sample:
                    mh.Chauvet.sample1<- mh.Chauvet[sI1.PPBChauvet]
                    sample2.index<- stratified.srs(cluster.size.PPBChauvet.h[sI1.PPBChauvet], mh.Chauvet.sample1, design = "SRSWOR")
                    sample2.PPBChauvet<- NULL
                    for(i in 1:nh) sample2.PPBChauvet<- rbind(sample2.PPBChauvet, data.frame(sample1.PPBChauvet[sample1.PPBChauvet$new.cluster==sI1.PPBChauvet[i], ][sample2.index[, 1][sample2.index[, 2]==i], ]) )

                    ## Replace the selected sample with the original sample with pi=n/N:
                    u=0
                    for(i in 1:nh){
                      I.Sboot<-rbinom(1, 1, as.numeric(aggregate(sample2.PPBChauvet$Pi1, by=list(sample2.PPBChauvet$new.cluster), FUN=unique)[i,2]))
                      if(I.Sboot==0) sample2.PPBChauvet[(u+1):(u+mh.Chauvet.sample1[i]), -dim(sample2.PPBChauvet)[2]]<- data.h[data.h$cluster==unique(sample2.PPBChauvet$cluster[(u+1):(u+mh.Chauvet.sample1[i])]), ]
                      u=u+mh.Chauvet.sample1[i]
                    }

                    boot.sample.index<- rbind(boot.sample.index, cbind(sample2.PPBChauvet$study.variable, sample2.PPBChauvet$ID, sample2.PPBChauvet$cluster, sample2.PPBChauvet$stratum, sample2.PPBChauvet$weight ))
                    boot.values<- rbind(boot.values, cbind(sample2.PPBChauvet$study.variable, sample2.PPBChauvet$weight) )
                  }

                  boot.statistic[a, b]<- statistic(as.numeric(boot.values[ ,1]), as.numeric(boot.values[ ,2]), q)
                  if(parameter=="mean" & !is.null(population.size)) boot.statistic[a, b]<- boot.statistic[a, b]/sum(population.size)

                  colnames(boot.sample.index)<- c("bootrstap.values", "ID", "cluster", "stratum", "bootstrap.survey.weights")
                  boot.sample.pp[[b]]<- as.data.frame(boot.sample.index)
                }

                boot.sample[[a]]<- boot.sample.pp

                }

            boot.mean<- mean(boot.statistic)
            boot.var<- mean(apply(boot.statistic, 1, var))

            return<- list(
              boot.statistic=boot.statistic,
              boot.var=boot.var,
              boot.mean=boot.mean,
              boot.sample=boot.sample,
              parameter=parameter,
              number.bootstrap.replicates=R
            )

          }
          else stop("This function can apply the method of Cheuvet (2007) if either SRSWOR, conditional poisson sampling is applied on the first stage.")
        }
      }
    }

  }
}
  return
}
