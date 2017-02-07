empty.df <-
function(header){
        df<-data.frame(matrix(matrix(rep(1,length(header)),1),1))
  		colnames(df)<-header
        return(df[NULL,])
}


