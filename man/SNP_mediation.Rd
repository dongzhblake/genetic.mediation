\name{SNP_mediation}
\alias{SNP_mediation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
SNP_mediation
}
\description{
SNP_mediation conducts mediation analysis at variant level.
}
\usage{
SNP_mediation(b.h,se.bh,a.h,se_ah)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.h}{
BETA estimate of the variant on the mediator (can be a vector with the same length of se.bh)
}
  \item{se.bh}{
SE of the BETA estimate (can be a vector with the same length of b.h)
}
  \item{a.bh}{
causal effect estimate of mediator on the outcome
}
  \item{se.ah}{
SE of the causal effect
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
SNP_mediation(0.2,0.05,1,0.5)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
