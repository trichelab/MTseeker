% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/newFilterMT.R
\name{newFilterMT}
\alias{newFilterMT}
\title{Filter SummarizedExperiment, DataFrame, MVRanges, or MVRangesList using user specified constraints}
\usage{
newFilterMT(
  vars,
  minTotalDepth = 20,
  minAltDepth = 2,
  minVAF = 0.9,
  minCovg = 20,
  fpFilter = TRUE,
  NuMT = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{vars}{variants, can be MVRanges\link{List} or DataFrame/SummarizedExperiment with colData()$\code{mtCovg}}

\item{minTotalDepth}{minimum number of total reads (ref + alt reads) (20)}

\item{minAltDepth}{minimum number of reads supporting a mutation needed (2)}

\item{minVAF}{minimum VAF to be considered (0.90)}

\item{minCovg}{minimum total depth (20, cf. Griffin, Genetics in Medicine 2014)}

\item{fpFilter}{apply Triska's homopolymer false positive filter? Only applicable for rCRS (FALSE)}

\item{NuMT}{apply the 0.03 VAF NuMT filter from Ju (GR 2015)? (FALSE)}

\item{verbose}{prints the filters applied (FALSE)}
}
\value{
\preformatted{       a filtered SE, data.frame, MVRanges, or MVRangesList 
}
}
\description{
Griffin et al. (Genetics in Medicine 2014) recommends 20x coverage for mtDNA
sequencing to have comparable error rates to Sanger sequencing.  By default,
that is the cutoff applied here to ensure halfway decent variant annotation.
}
\details{
Triska (Cancer Res, in revision) suggests a small number of masked regions
where homopolymers can be a problem; these are avoided if fpFilter

The NuMT filtration step (Ju, in eLife 2014, suggests a variant allele cutoff
of 0.03 to avoid false positive calls from nuclear-mitochondrial translocated
or 'NuMT' fragments) is also a useful tool to cut down on nonsensical calls,
although it may be important to use caution as low heteroplasmy can also
resolve into apparent near-homoplasmy at the single-cell level, at least in
our (TJT & co) experience.

As a consequence of the Wild West nature for published methods of high-
throughput mitochondrial sequence variant analysis at the time of writing
(2018), the default for this function is to filter on coverage only, as
the user is expected to determine what additional filters to apply. We
could envision changing these defaults down the road as standards congeal.
}
\examples{
filterMT(data.frame(sample="foo", mtCovg=1000))

}
