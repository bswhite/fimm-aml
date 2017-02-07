############################################################################################
#Developer name: Bhagwan Yadav
#Developed at Institute for Molecular Medicine Finland (FIMM)
############################################################################################

DSS <-
function(rdata,y=10,DSS.type=3,concn.scale=NULL,log=TRUE){
#rdata should be in the format containing IC50, SLOPE, MAX,MIN.Concentration,MAX.Concentration

if(is.null(concn.scale)& log){
stop("please set e.g. concn.scale=1e-9; 1e-9 for nano molar concentrations, 1e-6 for micro molar concentrations and so on.")
}

if(y<=0){
stop("Percent inhibition must be greater than Zero. Please, see manual for details.")
}

DSS <- empty.df(list('DSS'))

for (row in 1:nrow(rdata)){
  temp<-as.matrix(rdata[row,])
  if(isTRUE(ncol(rdata)==5)){
  ic50=as.numeric(as.character(str_replace(temp[1],',','.'))) #IC50
  b= as.numeric(as.character(str_replace(temp[2],',','.')))#slope;
  d=0 # min response
  a=as.numeric(as.character(str_replace(temp[3],',','.')))#max response
  Min.Conc_raw<- as.numeric(as.character(str_replace(temp[4],',','.')))#Min Concentration tested
  Max.Conc<- as.numeric(as.character(str_replace(temp[5],',','.')))#Max Concentration tested
  }
if(isTRUE(ncol(rdata)==6)){
  ic50=as.numeric(as.character(str_replace(temp[1],',','.'))) #IC50
  b= as.numeric(as.character(str_replace(temp[2],',','.')))#slope;
  d=as.numeric(as.character(str_replace(temp[3],',','.')))#min response
  a=as.numeric(as.character(str_replace(temp[4],',','.')))#max response
  Min.Conc_raw<- as.numeric(as.character(str_replace(temp[5],',','.')))#Min Concentration tested
  Max.Conc<- as.numeric(as.character(str_replace(temp[6],',','.')))#Max Concentration tested
  a=a-d
  d=d-d
}
  if(log){
  	Min.Conc<- log10(Min.Conc_raw*concn.scale) #
  	x2<-log10(Max.Conc*concn.scale)
  }else{
  	Min.Conc<- Min.Conc_raw #
  	x2<-Max.Conc
  }
    
  if(is.na(ic50)||is.na(b)||is.na(a)){
    DSS[row]<-NA
  }
  else if(isTRUE(ic50>=Max.Conc)){
    DSS[row]<-0
  }
  else if(b<0){
    DSS[row]<-0
  }
  else{
  	if(log){
 		c<-log10(ic50*concn.scale)
 	}else{
 		c<-ic50}
    if(a>y){
      if(y!=0){
        x1<-(c - ((log(a-y)-log(y-d))/(b*log(10))))

        if(isTRUE(x1 < Min.Conc)){x1<-Min.Conc}
  		else if(isTRUE(x1 > x2)){x1<-x2}
        }
      else {x1<-Min.Conc}

    # This is a logistic function used in Dotmatics.com
    # y = d+(a-d)/(1+10^(b*(c-x)))
    #inverse function
    # x = c - ((log(a-y)-log(d-y))/(b*log(10)))
      
    int_y=(((((a-d)*log(1+10^(b*(c-x2))))/(b*log(10)))+a*x2)-((((a-d)*log(1+10^(b*(c-x1))))/(b*log(10)))+a*x1)) - (y*(x2-x1))
    if(int_y<0){int_y<-0}
    total_area<-(x2-Min.Conc)*(100-y)
    
    if(DSS.type==1){
    	norm_area<-((int_y/total_area)*100)#DSS1
    }
    
    if(DSS.type==2){
    	norm_area<-((int_y/total_area)*100)/log10(a)#DSS2
    }
    
	if(DSS.type==3){
    	norm_area<-((int_y/total_area)*100)*(log10(100)/log10(a))*((x2-x1)/(x2-Min.Conc)) #DSS3
    }
    
      DSS[row]<-round(norm_area,digits=8)

    } else {DSS[row]<-0} 
   } 
  }
return (data.frame(DSS))
}

