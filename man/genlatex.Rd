%                               -*- Mode: Rd -*- 
% genlatex.Rd --- 
% Author          : Claus Dethlefsen
% Created On      : Sat May 25 23:02:07 2002
% Last Modified By: Claus Dethlefsen
% Last Modified On: Thu Dec 04 13:39:54 2008
% Update Count    : 20
% Status          : Unknown, Use with caution!
% 

\name{genlatex}
\alias{genlatex}
\alias{genpicfile}
\encoding{latin1}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{From a network family, generate LaTeX output}

\description{The networks in a network family is arranged as
  pictex-graphs in a LaTeX-table. 
}

\usage{
genlatex(nwl,outdir="pic/",prefix="scoretable",picdir="",picpre="pic",
         ncol=5,nrow=7,width=12/ncol,vadjust=-1.8)
genpicfile (nwl,outdir="pic/",prefix="pic",w=1.6,h=1.6,bigscale=3) 
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{nwl}{object of class \code{networkfamily} containing a list of
    objects of class \code{network}.}
  \item{outdir}{character string, the directory for storing output.}
  \item{prefix}{character string, the filename (without extension) of the LaTeX file. The
    filenames of the picfiles begin with the given prefix.} 
  \item{picdir}{character string, the directory where pic-files are stored.}
  \item{picpre}{character string, prefix for pic-files.}
  \item{ncol}{integer, the number of columns in LaTeX table.}
  \item{nrow}{integer, the number of rows in LaTeX table.}
  \item{width}{numeric, the width of each cell in the LaTeX table.}
  \item{vadjust}{numeric, the vertical adjustment in LaTeX table.}
  \item{w}{numeric, the width of pictex objects}
  \item{h}{numeric, the height of pictex objects}
  \item{bigscale}{numeric, the scaling of the best network, which is output in 'nice.tex'}
  
}

\value{Files:

  \item{\{outdir\}\{picpre\}xx.tex}{one pictex file for each network in the
  network family, indexed by xx.}

  \item{\{outdir\}\{prefix\}.tex}{LaTeX file with table including all pictex files.}

  \item{\{outdir\}\{picpre\}nice.tex}{pictex file with the best network.}
}

\seealso{\code{\link{networkfamily}}, \code{\link{pictex}}}

% \references{
%   Further information about \bold{deal} can be found at:\cr
%   \url{http://www.math.aau.dk/~dethlef/novo/deal}.
% }

\author{
  Susanne Gammelgaard Bottcher, \cr
  Claus Dethlefsen \email{rpackage.deal@gmail.com}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\examples{
data(rats)
allrats <- getnetwork(networkfamily(rats,network(rats)))
allrats <- nwfsort(allrats)

\dontrun{dir.create("c:/temp")}
\dontrun{genpicfile(allrats,outdir="c:/temp/pic/")}
\dontrun{genlatex(allrats,outdir="c:/temp/pic/",picdir="c:/temp/pic/")}

## LATEX FILE:
#\documentclass{article}
#\usepackage{array,pictex}
#\begin{document}
#\input{scoretable}
#\input{picnice}
#\end{document}

#data(ksl)
#ksl.nw  <- network(ksl)
#ksl.prior <- jointprior(ksl.nw,64)
#mybanlist <- matrix(c(5,5,6,6,7,7,9,
#                    8,9,8,9,8,9,8),ncol=2)
#banlist(ksl.nw) <- mybanlist                  
#ksl.nw <- getnetwork(learn(ksl.nw,ksl,ksl.prior))
#ksl.search <- autosearch(ksl.nw,ksl,ksl.prior,
#                      trace=TRUE)
#ksl.searchlist <- makenw(ksl.search$table,ksl.search$nw)
#ksl.searchlist <- nwfsort(ksl.searchlist)
\dontrun{genpicfile(ksl.searchlist)}
\dontrun{genlatex(ksl.searchlist)}

}

\keyword{iplot}
