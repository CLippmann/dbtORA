sortby <- function(mat,column,desc=FALSE){
  if(column>dim(mat)[2]|column<=0){
    return("The argument column must be 0<column<=columnnumber.");
  }
  if(class(mat)=="matrix"){
    # order returns a permutation(index) and rearranges the columns.
    return(mat[order(mat[,column],decreasing=desc),])
  }
  else{
    print("The first argument must be a matrix.")
  }
}#end sortby
