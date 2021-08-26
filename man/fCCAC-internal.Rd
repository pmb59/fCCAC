\name{fCCAC}
\title{fCCAC internal functions}
\description{Internal undocumentation functions}

\alias{multiplot}

\keyword{internal}

\examples{
  library(ggplot2)
  # This example uses the ChickWeight dataset, which comes with ggplot2
  # http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
  p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) + geom_line() + ggtitle("Growth curve for individual chicks")
  multiplot(p1, p1, cols=1)
}
