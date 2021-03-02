########################################################################
## Building predictive Random Forest models for biodiversity recovery ##
########################################################################

## This script was used to build a predictive model for the recovery of 
## invertebrate abundance in second-growth forests. The same analysis was 
## repeated for invertebrate richness and for abundance and richness of 
## vertebrates and plants 

## Load biodiversity data
datafull<-read.table("database.full.txt", h=T, stringsAsFactors=FALSE, fileEncoding="latin1")

## Select Invertebrates data
data.invert<-subset(datafull, datafull$Invertebrates==1)

## Invertebrates data from landscapes restored through natural regeneration in tropical & subtropical regions
data.invert.trop<-subset(data.invert, data.invert$Latitude<35 & data.invert$Latitude>-35 & data.invert$Restoration_activity=="passive")

## Load Ecological & socioeconomic factors measured at five distinct buffer sizes
buffer<-readRDS("Dados.buffer.RData")

## Remove extreme values (Outlyers)
invert.so<-subset(data.invert.trop, data.invert.trop$RR<3 & data.invert.trop$RR>-3)

## Select abundance data
invert.abund<-subset(invert.so, invert.so$Ecological_metric=="abundance")

## Calculating landscape deviation (here called RR_var) in biodiversity recovery
data.invert<-invert.abund
data.invert$N<-1
data.invert2<-data.invert2[order(data.invert2$Site),]
N<-as.data.frame(tapply(data.invert2$N, data.invert2$Site, sum))
colnames(N)<-"N"
data.invert2$RR<-data.invert2$RR*data.invert2$RR
RR<-as.data.frame(tapply(data.invert2$RR, data.invert2$Site, sum))
RR_var<-RR/(N-1)
colnames(RR_var)<-"RR_var"
RR_var$Site<-rownames(RR_var)
data.invert3<-merge(data.invert2, RR_var, by="Site")

obj.invert<-as.data.frame(data.invert3[,c(1,2,3,4,8)])
colnames(obj.invert)
RRvar.buffer<-merge(obj.invert, buffer, by="Site")
RRvar.buffer<-RRvar.buffer[order(RRvar.buffer$Site),]
RRvar.buffer[1:20,]
obj.invert2<-RRvar.buffer[ !duplicated(RRvar.buffer$Site), ]
class(obj.invert2$Longitude)
obj.invert2$Longitude<-as.numeric(obj.invert2$Longitude)

obj.invert3 <-obj.invert2[!is.infinite(rowSums(obj.invert2)),]
obj.invert3<-na.omit(obj.invert3)

response.invert<-obj.invert3[1:5]
predictor<-obj.invert3[,c(6:269)]

## List of factors to be passed to random forest
var.list = list('AnualTemp' = grep('AnualTemp', colnames(predictor)),
                'BAI1317' = grep('BAI1317', colnames(predictor)),
                'CEC' = grep('CEC', colnames(predictor)),
                'Cropland' = grep('Cropland', colnames(predictor)),
                'DriestQuartes' = grep('DriestQuartes', colnames(predictor)),
                'ForestUntl10' = grep('ForestUntl10', colnames(predictor)),
                'GrossDef' = grep('GrossDef', colnames(predictor)),
                # 'HFPrint' = grep('HFPrint', colnames(predictor)),
                'IDH03' = grep('IDH03', colnames(predictor)),
                'IDH9010' = grep('IDH9010', colnames(predictor)),
                'NBR1317' = grep('NBR1317', colnames(predictor)),
                'NPP0010' = grep('NPP0010', colnames(predictor)),
                'NPP03' = grep('NPP03', colnames(predictor)),
                'OpCostIISMean' = grep('OpCostIISMean', colnames(predictor)),
                'PastAgric' = grep('PastAgric', colnames(predictor)),
                'PercUrbArea09' = grep('PercUrbArea09', colnames(predictor)),
                'PrecSeanslty' = grep('PrecSeanslty', colnames(predictor)),
                'RuralPop' = grep('RuralPop', colnames(predictor)),
                'RuralPvty' = grep('RuralPvty', colnames(predictor)),
                'SustainablePA' = grep('SustainablePA', colnames(predictor)),
                'TotalRoadDensity' = grep('TotalRoadDensity', colnames(predictor)),
                'WDeficit03' = grep('WDeficit03', colnames(predictor)),
                'WDeficit8910' = grep('WDeficit8910', colnames(predictor)),
                'YrlPrec' = grep('YrlPrec', colnames(predictor)),
                'bldfie' = grep('bldfie', colnames(predictor)),
                'cdd' = grep('cdd', colnames(predictor)),
                'cropIIS' = grep('cropIIS', colnames(predictor)),
                'croppastureIIS' = grep('croppastureIIS', colnames(predictor)),
                'deltaGHSL9015' = grep('deltaGHSL9015', colnames(predictor)),
                'deltaGPWCount0010' = grep('deltaGPWCount0010', colnames(predictor)),
                'deltaGPWDensity0010' = grep('deltaGPWDensity0010', colnames(predictor)),
                'elevation' = grep('elevation', colnames(predictor)),
                'f2003' = grep('f2003', colnames(predictor)),
                'fty' = grep('fty', colnames(predictor)),
                'govrnce' = grep('govrnce', colnames(predictor)),
                'pastoIIS' = grep('pastoIIS', colnames(predictor)),
                'phihox' = grep('phihox', colnames(predictor)),
                'sa' = grep('sa', colnames(predictor)),
                'slope' = grep('slope', colnames(predictor)),
                'sndppt' = grep('sndppt', colnames(predictor)),
                'strictlyPA' = grep('strictlyPA', colnames(predictor)),
                'treecover' = grep('treecover', colnames(predictor)),
                'urb' = grep('urb', colnames(predictor)),
                'wfire' = grep('wfire', colnames(predictor)))

var.list

## Function that receives a list of factors at multiple landscape sizes 
## (i.e. multiple buffer sizes) and returns the strongest scale of effect 
## of each factor using Random Forest algorithm

library(randomForest)

engine = function(iter.vars){
  set.seed(201)
  res = randomForest(response.invert[,5]~., data=iter.vars)
  varImp = as.data.frame(varImpPlot(res))
  node.purity = as.data.frame(varImp$IncNodePurity)
  node.purity[,2] = rownames(varImp)
  colnames(node.purity) = c("NodePurity", "Varnames")
  node.purity = node.purity[order(node.purity$NodePurity, decreasing=T),]
  buffer.size = node.purity$Varnames[1]
  ColNumber = grep(buffer.size, colnames(predictor))
  out.list = list(res=res, buffer.size=buffer.size, ColNumber=ColNumber)
  return(out.list)
  
}

iter.res = as.list(1:length(var.list))
names(iter.res) = names(var.list)

## Loop to apply engine to each member of var.list
for (iter.names in names(var.list)){
  iter.cols = var.list[[iter.names]]
  iter.vars = predictor[,iter.cols]
  #print('start')
  #print(iter.vars)
  iter.res[[iter.names]] = engine(iter.vars)
  iter.buffer = iter.res$buffer.size
  
}

names(iter.res)

buffer.sizes = data.frame(buffer.sizes = sapply(names(var.list), function(x){iter.res[[x]]$buffer.size}))
col.numbers = data.frame(col.numbers = sapply(names(var.list), function(x){iter.res[[x]]$ColNumber}))

buffer.sizes
col.num<-col.numbers$col.numbers
var.sel.buff<-predictor[,col.num]

## Factor selection for prediction with VSURF varying the parameter mtry 
## (number of factors tested for in each node of a tree)

library(VSURF)

## mytr=default (= 1/3 of the total number of factors)
set.seed(201)
vsurf.invert.buff6<-VSURF(var.sel.buff, response.invert[,3], ntree=5000)
vsurf.invert.buff6$varselect.pred

## mytr=2
set.seed(201)
vsurf.invert.buff<-VSURF(var.sel.buff, response.invert[,3], mtry=2, ntree=5000)
vsurf.invert.buff$varselect.pred

## mytr=4
set.seed(201)
vsurf.invert.buff2<-VSURF(var.sel.buff, response.invert[,3], mtry=4, ntree=5000) 
vsurf.invert.buff2$varselect.pred

## mytr=6
set.seed(201)
vsurf.invert.buff3<-VSURF(var.sel.buff, response.invert[,3], mtry=6, ntree=5000) 
vsurf.invert.buff3$varselect.pred

## mytr=8
set.seed(201)
vsurf.invert.buff4<-VSURF(var.sel.buff, response.invert[,3], mtry=8, ntree=5000)
vsurf.invert.buff4$varselect.pred

## mytr=10
set.seed(201)
vsurf.invert.buff5<-VSURF(var.sel.buff, response.invert[,3], mtry=10, ntree=5000)
vsurf.invert.buff5$varselect.pred


## Final Random Forest model with selected variables
set.seed(201)
RF.LandVar_invertebrates<-randomForest(response.invert[,5]~., data=var.sel.buff[,c(2,31)], ntree=2000, importance=T)
print(RF.LandVar_invertebrates)
getTree(RF.LandVar_invertebrates, labelVar=T)
varImp<-varImpPlot(RF.LandVar_invertebrates)

Dados.RF.Invert.abundance<-as.data.frame(cbind(response.invert[,c(1,3,4,5)], var.sel.buff[,c(2,31)]))
write.csv(Dados.RF.Invert.abundance, file="Dados.RF.Invert.abundance.csv")
Dados.RF.Invert.abundance<-read.csv("Dados.RF.Invert.abundance.csv")

## Calculate mean restoration age in each landscape
age<-as.data.frame(cbind(data.invert.trop$Site, data.invert.trop$Time_restored))
colnames(age)<-c("Site","Age")
mean_age<-as.data.frame(tapply(age$Age, age$Site, mean))
mean_age$Site<-rownames(mean_age)
colnames(mean_age)[1]<-"Mean_age"

age_abundance_invert<-merge(mean_age, Dados.RF.Invert.abundance, by="Site")
hist(age_abundance_invert$Mean_age)
mean(age_abundance_invert$Mean_age, na.rm = T)
median(age_abundance_invert$Mean_age, na.rm = T)
range(age_abundance_invert$Mean_age, na.rm = T)

#################################################################################
