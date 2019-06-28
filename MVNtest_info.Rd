\name{MVNTest}
\alias{MVNTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cramer's Method Test of Multivariate Normality
}
\description{
Performs a test for multivariate normality based on Cramer's thereom,
that if all linear projections of the p-dimensions in a dataset X are univariate normal then X is multivariate normal.
A high number of projections is tested. A q-value is returned after controlling for the false discovery rate.
}
\usage{
MVNTest(X, numproj)
}

\arguments{
  \item{X}{
A dataset.
}
\item{numproj}
{
The number of desired projections.
}
}
\details{
The features in the dataset must be in the columns. The default number of projections is set to 100000.
}
\value{
The function returns a q-value for the test af
}
\references{
}
\author{
Ranjan Maitra
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
load(iris)
setosa = iris[c(1:50),c(1:4)]
MVNtest(as.matrix(setosa))
MVNtest(as.matrix(setosa), 1000)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
