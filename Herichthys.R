#######################################
### cleannig geographic data points####
#######################################

library(spThin)
coor<-read.csv("D:/H_carpintis_geo_points.csv")
thinned_dataset_full<-thin( loc.data = coor, lat.col = "lat", long.col = "lon", spec.col = "esp", thin.par = 2, reps = 10, locs.thinned.list.return = TRUE, write.files = TRUE, max.files = 1, out.dir = "D:/Models/H_carpintis_geo_points.csv", out.base = "H_carpintis", write.log.file = TRUE,log.file = "H_carpintis.csv")

#############################################################################
### Variance inflation factor, used to reduce the number climatic variables##
#############################################################################

library(usdm)
inputDir<-"D:/Climatic_variables/SET current 1"
files <- list.files(inputDir, pattern = '.asc$', full.names = TRUE)
r <- stack(files)
v1<-vifcor(r, th=0.9)
v1
v2<-vifstep(r,th=10)
v2

######################
#####    ENM     #####
######################

library(kuenm)
library(ENMeval)
library(raster)
library(sp)
library(maptools)
library(rgeos)
library(dismo)
library(rJava)
library(jsonlite)
library(grDevices)
library(devtools)
library(rmaxent)
library(sdmpredictors)

system.file('java', package="dismo")
maxent()
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
if (file.exists(jar) & require(rJava)) {
  fnames<- list.files(path="D:/HERICHTHYS/H car/Current/Set 1",
                      pattern = '.asc$', full.names = TRUE)
  predictors <- stack(fnames)
  crs(predictors) <- "+proj=utm +zone=1 +datum=WGS84"
  ocurrencias<- list.files(path="D:/Models/H_carpintis_geo_points",pattern = "*.csv", full.names= TRUE)
  occ <- read.table(ocurrencias, header=TRUE, sep=',')[,-1]
  fold <- kfold(occ, k=10)
  occtest <- occ[fold == 1, ]
  occtrain <- occ[fold != 1, ]
  me<- maxent(predictors, occtrain)
  plot(me)
  r <- predict(me, predictors)
  plot(r)
  points(occ)
  bg<- data.frame(randomPoints(predictors, 10000))
  pvtest <- data.frame(extract(predictors, occtest))
  avtest <- data.frame(extract(predictors, bg))
  e2 <- evaluate(me, p=pvtest, a=avtest)
  testp <- predict(me, pvtest)
  head(testp)
  testa <- predict(me, avtest)
  e3 <- evaluate(p=testp, a=testa)
  e3
  threshold(e3)
  plot(e3, 'ROC')
}
tiff("D:/HERICHTHYS/H car/Current/Results/Variables_contribution.tiff")
plot(me)
dev.off()
names(bg) <- names(occ)
tune.args <- list(fc = c("L","LQ","LQH","H"), rm = 1:5)
eval2 <- ENMevaluate(occ, predictors, bg,tune.args = tune.args, method='jackknife', RMvalues=c(1,2, 5, 10, 15, 20),parallel = TRUE,clamp = TRUE, algorithm='maxent.jar')
eval2
eval2@results
eval2@results[which(eval2@results$delta.AICc==0),]
aic.opt <- eval2@models[[which(eval2@results$delta.AICc==0)]]
aic.opt
aic.opt@results
#eval.variable.importance(eval2)
var.importance(aic.opt)
eval2@predictions
write.csv(eval2@results, file= "D:/HERICHTHYS/H car/Current/Results/EVAL_results.csv")
write.csv(aic.opt@results, file= "D:/HERICHTHYS/H car/Current/Results/AIC_results_bestmodel.csv")
write.csv(aic.opt@results, file= "D:/HERICHTHYS/H car/Current/Results/aic.opt@results.csv")
writeRaster((eval2@predictions[[which(eval2@results$delta.AICc==0)]]),filename="D:/HERICHTHYS/H car/Current/Results/H car_raw.asc")

######################
##### Partial ROC#####
######################

model <- raster("D:/HERICHTHYS/H car/Current/Results/H car_raw.asc")
thres <- 5
rand_perc <- 50
iterac <- 500
p_roc <- kuenm_proc(occ.test = occtest, model = model, threshold = thres,rand.percent = rand_perc, iterations = iterac)
p_roc$pROC_summary
p_roc$pROC_results
write.csv(p_roc$pROC_results, file= "D:/HERICHTHYS/H car/Current/Results/pROC_results.csv")


