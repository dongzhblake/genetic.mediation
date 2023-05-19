\name{read_BOLT}
\alias{read_BOLT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read_BOLT
}
\description{
read_BOLT extracts the causal effect estimate and its standard error from BOTL-REML output log file.
}
\usage{
read_BOLT(BOLT_result)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BOLT_result}{
The file name of BOLE-REML log file
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{alphahat_adj}{Extracted causal effect of the mediator to the outcome}
\item{SE}{Standard error of the estimated causal effect}
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
read_BOLT("example_folder/output.log")
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
