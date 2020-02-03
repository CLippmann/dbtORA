\name{GOroot2TermPaths}
\alias{GOroot2TermPaths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to get all paths from the gene ontology root of the DAG to one target term.
}
\description{
Function to get all paths from the gene ontology root of the DAG to one target term.
}
\usage{
GOroot2TermPaths(TargetTerm, AdjMatrix, GOtermNr, GOroot = 8150)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TargetTerm}{
Numeric;
GO term number of the term that the paths should be found for.}
  \item{AdjMatrix}{
Numeric;
Adjacency matrix of (only) all ancestors of TargetTerm and TargetTerm itself.
AdjMatrix[i,j] == 1 iff there exists an edge in 
GO-DAG from node i to node j (i is parent of j).}
  \item{GOtermNr}{
Numeric;
Vector of GO term numbers corresponding to rows and columns of AdjMatrix.
Has to contain TargetTerm.}
  \item{GOroot}{
Numeric; Default = 8150
The GO term number of the root of the ontology in which the target term 
and all GOtermNrs are. One of: 
8150 for biological process, 
3674 for molecular function or 
5575 for cellular component.}
}

\value{
\item{Taxonomy}{Numeric;
A list of indices such that GOtermNr[Taxonomy[[i]]] is the i-th path from 
\code{GOroot} to \code{TargetTerm}.}
}

\author{
CL
}

