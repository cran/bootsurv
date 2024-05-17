\name{data_pop, data_pop_st, data_pop_clust, data_pop_stclust, data_samp_srs, data_samp_stsrs, data_samp_clust, data_samp_stclust}
\alias{data_pop}
\alias{data_pop_st}
\alias{data_pop_clust}
\alias{data_pop_stclust}
\alias{data_samp_srs}
\alias{data_samp_stsrs}
\alias{data_samp_clust}
\alias{data_samp_stclust}
\title{Populations and samples gerenated in the \code{bootsurv} package}
\description{
  This package contains multiple datasets described below.
}
\section{Datasets}{
  \describe{
    \item{\code{data_pop}}{This is a population of size 6,000. This data set contains a column of generated study variable, labeled as \code{study.variable}.}

    \item{\code{data_pop_st}}{This dataset represents a population of size 17,800, divided into 5 strata. It includes a column for the generated study variable, labeled as \code{study.variable}, and a column identifying the strata, labeled as \code{stratum}. The subpopulation sizes within each stratum are as follows: 4,500, 6,300, 3,500, 2,000, and 1,500, respectively.}

    \item{\code{data_pop_clust}}{This dataset represents a population consisting of 10,048 units distributed across 200 clusters. The number of units within each cluster was generated using a Poisson distribution with a mean of 50. It includes columns for the generated study variable, labeled as \code{study.variable}, and cluster identification, denoted as \code{cluster}.}

    \item{\code{data_pop_stclust}}{This dataset represents a population with 14,511 units distributed across three strata, consisting of 100, 125, and 65 clusters, respectively. The number of units within each cluster was generated using a Poisson distribution with a mean of 50. It includes columns of the generated study variable, labeled as \code{study.variable}, stratum identification, labeled as \code{stratum}, and cluster identification within each stratum, labeled as \code{cluster}.}

    \item{\code{data_samp_srs}}{This dataset comprises a sample of size 1,850, obtained through simple random sampling without replacement from the \code{data_pop} dataset.}

    \item{\code{data_samp_stsrs}}{This dataset represents a sample of size 5,350 obtained through stratified simple random sampling without replacement from the stratified population \code{data_pop_st}. The sample consists of subsample sizes of 1,350, 1,900, 1,050, 600, and 450.}

    \item{\code{data_samp_clust}}{This sample was drawn using a two-stage cluster sampling method, with simple random sampling without replacement applied at each stage. The sample is drawn from the \code{data_pop_clust} dataset. In the first stage, approximately 20\% of clusters were selected. Subsequently, within each selected cluster, approximately 15\% of units were sampled.}

    \item{\code{data_samp_stclust}}{A stratified two-stage cluster sampling method is applied to draw this sample from the \code{data_pop_stclust} dataset. In each stratum, simple random sampling without replacement is applied at each stage. The first stage sampling fraction is approximately 20\%, and the overall second stage sampling is approximately 15\%.}
  }
}
