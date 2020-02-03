\name{termpathsHeadlines}
\alias{termpathsHeadlines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates headlines for each path from root to GO detail terms.
}
\description{
Calculates headlines for each path from root to GO detail terms. Headlines
represent most important nodes in these paths.
}
\usage{
termpathsHeadlines(AdjMatrix, GOtermNr, Importance, OntologyNr = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AdjMatrix}{
Numeric;
Adjacency matrix of GO terms. AdjMatrix[i,j] == 1 if there exists an edge in 
GO-DAG from node i to node j (i is parent of j).}
  \item{GOtermNr}{
Numeric;
Vector of all GO term numbers that are in the considered DAG. 
Corresponding to rows and columns of AdjMatrix.}
  \item{Importance}{
Numeric; Values between 0 and 100.
Vector of importance values (= rowwise minimum of (partial shannon) information 
value and certainty) for GO terms.}
  \item{OntologyNr}{
Numeric; Default = 1
To select the ontology in which the GOtermNrs are. One of: 
1 for biological process, 
2 for molecular function or 
4 for cellular component.}
}
\details{
PLEASE NOTE: 
All given GO terms have to be in the specified ontology and all 
ancestors from GO terms to ontology root (including the root-term itself) must 
be included in the given vector of GO terms.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{Headlines}{Numeric;
Vector of GO term numbers that are headlines, the GO term numbers with maximum
importance value for each path from detail to GOroot.}
\item{AllTaxonomies}{List of paths from GOroot to all details given by AdjMatrix 
in form of GO term numbers vectors.
(E.g. AllTaxonomies[[1]]==c(8150, 44699, 44763, 22402, 51231, 22) is a path from 
BP-root "GO:0008150" to detail "GO:0000022"via nodes 44699, 44763, 22402 and 51231
in BP gene ontology.)}
\item{MaxImportanceInd}{Vector of indices which indicate the position of the highest value of 
importance in the paths in AllTaxonomies.
(E.g. MaxImportanceInd==c(2,4,5) then AllTaxonomies[[1]][2] would be the GO term 
with the highest importance value (i.e. headline) in the first path, 
AllTaxonomies[[2]][4] the one in the second path and AllTaxonomies[[3]][5] the one
in the third path.)}
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{GOroot2TermPaths}}.
}


