\name{updateORAdatabase}
\alias{updateORAdatabase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to update the database used for ORA. 
}
\description{
With this function user can update database used for ORA. Determining the point in
time this should be done himself. (E.g., to keep the database untill project is done.)
The new files will be saved in system.file('extdata', package = 'ORA'). Updating the
data is only necessary if user updates R to newer version and packages GO.db and 
org.Hs.eg.db are updated. Timestamp of data provided in the package ORA: 31 January 2018.
Function has then to be called every time after restart of R session.
}
\usage{
updateORAdatabase()
}

\author{
CL
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{dbtORA}}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
