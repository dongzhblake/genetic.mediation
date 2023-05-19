\name{GWAS_mediation}
\alias{GWAS_mediation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GWAS_mediation
}
\description{
GWAS_mediation conducts mediation analysis across the genome by directly taking GWAS summary statistics
}
\usage{
GWAS_mediation(GWAS_M, a.h, se_ah, GWAS_Y = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GWAS_M}{
GWAS data for the mediator
}
  \item{a.h}{
causal effect estimate of mediator on the outcome
}
  \item{se_ah}{
SE of the causal effect
}
  \item{GWAS_Y}{
GWAS data for the outcome (not necessary, for calculating mediated proportion only)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{mediation.h}{
  mediation effect estimate of variant(s)}
  \item{se_mediation.h}{
  SE of mediation effect estimate of variant(s)}
  \item{mediation.p}{
  P value of mediation effect estimate of variant(s)}
    \item{mediation_proportion}{
  Estimated mediated proportion of variant(s)}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
GWAS_M=read.table("example_folder/sumstats_M_test",header=T)
GWAS_Y=read.table("example_folder/sumstats_Y_test",header=T)
GWAS_mediation(GWAS_M,1,0.5)
GWAS_mediation(GWAS_M,1,0.5,GWAS_Y)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
