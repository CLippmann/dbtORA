\name{termDescription}
\alias{termDescription}
\title{
Yields description of inputtet GO-term.
}
\description{
Yields description of inputtet GO-term ID.
}
\usage{
termDescription(GOtermId)
}
\arguments{
  \item{GOtermId}{
Vector\code{[1:n]} of GO-term IDs, e.g. \code{"GO:0008150"}.
}
}
\details{
Requires package \code{GO.db}.
}
\value{
\item{GOTermDescription}{Vector\code{[1:n]} of strings that denote the GO terms,
e.g. \code{"biological process"}.}
}

\author{
ALU, MT, CL
}
