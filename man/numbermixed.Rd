%                               -*- Mode: Rd -*- 
% numbermixed.Rd --- 
% Author          : Claus Dethlefsen
% Created On      : Sat May 25 23:03:38 2002
% Last Modified By: Claus Dethlefsen
% Last Modified On: Mon Jan 05 14:37:29 2009
% Update Count    : 15
% Status          : Unknown, Use with caution!
% 
\encoding{latin1}
\name{numbermixed}
\alias{numbermixed}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{The number of possible networks}

\description{Calculates the number of different directed acyclic graphs
  for a set of discrete and continuous nodes. }

\usage{
numbermixed(nd,nc)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{nd}{an integer, which gives the number of discrete nodes.}
  \item{nc}{an integer, which gives the number of continuous nodes.}
}

\details{
  No arrows are allowed from continuous nodes to discrete nodes. Cycles
  are not allowed. The
  number of networks is given by Bottcher (2003), using the result in
  Robinson (1977).

  When nd+nc>15, the procedure is quite slow.
}

\value{
  A numeric containing the number of directed acyclic graphs with the given node
  configuration. 
}

\references{
  Bottcher, S.G. (2003). Learning Conditional Gaussian
  Networks. Aalborg University, 2003.

  Robinson, R.W. (1977). Counting unlabeled acyclic digraphs, Lecture
  Notes in Mathematics, 622: Combinatorial Mathematics.

}


\author{
  Susanne Gammelgaard Bottcher, \cr
  Claus Dethlefsen \email{rpackage.deal@gmail.com}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\examples{
numbermixed(2,2)
\dontrun{numbermixed(5,10)}
}

\keyword{models}
