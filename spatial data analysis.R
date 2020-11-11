library(geoR)
library(geospt)
library(gstat)
library(sgeostat)
library(nortest)  ##for anderson test
library(sp)    ##for spplot

##########################
##### spatial plot ######
##########################
Densidade<-as.geodata(wrc,data.col=3)
W<-data.frame(Densidade)
x<-W$CoordX
y<-W$CoordY
z<-W$data
plot(y~x,typ="n")
text(y~x,labels=z,cex=0.6)
##########################
#### variogram cloud #####
##########################
vc=variog(Densidade,bin.cloud=TRUE,estimator.type="classic")
plot(vc, bin.cloud=TRUE)
##########################
##normal bodan-plot###
##########################
summary(Densidade)
plot(Densidade)
data<-as.data.frame(Densidade)
data
coord<-Densidade$coords
coord
z<-matrix(Densidade$data,nc=1)
Densidadee<-cbind(coord,z)
dimnames(Densidadee)<-list(paste(),paste(c("x","y","z")))
##########################
###### H-parakonesh ######
##########################
Densidade <- as.geodata(wrc, data.col=3)
as.data.frame(Densidade)
summary(Densidade)
plot(Densidade)
data<-as.data.frame(Densidade)
data
coord<-Densidade$coords
coord
z<-matrix(Densidade$data,nc=1)
Densidade<-cbind(coord,z)
dimnames(Densidade)<-list(paste(),paste(c("x","y","z")))
##########################
###### H-parakonesh ######
##########################
ca20<-as.data.frame(Densidade)
coordinates(ca20) = ~x+y
hscat(z~1,ca20, c(5,10,15,20,25,30,35),col=6,pch=16)
##########################
####### plot data ########
##########################
Plot<-function(data){
plot(data[,c(1,2)],main="plot x~y")      
stem(data[,3],scale=2)		
X11()
par(mfrow=c(2,2))
Data<-as.geodata(data)  
plot(Data)  
}
Plot(Densidade)
##########################
###### test normal #######
##########################
normality<-function(data){
shapiro<-shapiro.test(data[,3])         	
anderson<-ad.test(data[,3]) 
qqnorm(data[,3])
qqline(data[,3])
return(list(shapiro,anderson))
}
normality(data)

######### remove trend ###########
##########################
####### regression #####
##########################
# regression
plot(data$CoordY,data$data, xlab="Z",ylab="Y",bg="yellow")
lmST<-lm(data~CoordY,data=data)
abline(lmST, col = "red",lty=1,lwd=2)
error1<-c()
for(i in 1:250){
error1[i]<-predict(lm(data~CoordX,data=data[-i,]),new=data[i,])-data[i,3]
}
cv_lm<-mean(error1^2)
summary(lmST)
##########################
######### spline ########
##########################
smooth<-smooth.spline(data[,2],data[,3])
lines(smooth,col=4,lty=2,lwd=2)
legend("topright",c("LM","Spline"),col=c("red","blue"),lty=c(1,2),lwd=c(2,2))
error2<-c()
for(i in 1:250){
error2[i]<-predict(smooth.spline(data[-i,1],data[-i,3]),data[i,1])$y-data[i,3]
}
cv_Spline<-mean(error2^2)
cv_lm-cv_Spline
cv_Spline
cv_lm
predictZ<-predict(smooth,data[,2])$y
newZ<-data[,3]-predictZ
data[,3]<-newZ
NewData<-data
Plot(NewData)
normality(NewData)
##########################
######### hining ########
##########################
H<-function(data){
q1<-quantile(data[,3],probs=.25)
q3<-quantile(data[,3],probs=.75)
fu<-q3+1.5*(q3-q1)
fl<-q1-1.5*(q3-q1)
for(i in data[,3]){
if(i>fu||i<fl)print(i)                      
}
}
H(Densidade)
##########################
### mean&median plot ###
##########################
MM<-function(Data){
data<-as.matrix(Data)
x<-data[!duplicated(data[,1]),1]
meanx<-c()
medianx<-c()
Q<-c()->sigmahat
u1<<-c()
for(i in 1:length(x)){
y1<-data[data[,1]==x[i],]
meanx[i]<-apply(y1,2,mean)[3]
medianx[i]<-apply(y1,2,median)[3]
Q[i]<-apply(y1,2,IQR)[3]
sigmahat[i]<-Q[i]/(2*0.6745)
u1[i]<<-(sqrt(length(y1[,3]))*(meanx[i]-medianx[i]))/(0.7555*sigmahat[i])
if(u1[i]>3)print(u1[i])
}
y<-data[!duplicated(data[,2]),2]
meany<-c()
mediany<-c()
Q<-c()->sigmahat
u2<<-c()               
for(i in 1:length(y)){
z1<-data[data[,2]==y[i],]
meany[i]<-apply(z1,2,mean)[3]
mediany[i]<-apply(z1,2,median)[3]
Q[i]<-apply(z1,2,IQR)[3]
sigmahat[i]<-Q[i]/(2*0.6745)
u2[i]<<-(sqrt(length(z1[,3]))*(meany[i]-mediany[i]))/(0.7555*sigmahat[i])
if(u2[i]>3)print(u2[i])
}
par(mfrow=c(1,2))            
plot(x,meanx,type="p",pch=10,main="column summaries",col=6,xlab="x-coordinate",ylab="mean and median",ylim=c(min(meanx,medianx),max(meanx,medianx)))
lines(x,medianx,type="p",pch=16,col=4)
plot(y,meany,type="p",pch=10,main="row summaries",col=6,xlab="y-coordinate",ylab="mean and median")
lines(y,mediany,type="p",pch=16,col=4)
}
MM(Densidade)
u2
u1
##########################
###### hamsan gardi ######
##########################
ca20<-as.data.frame(NewData)
coordinates(ca20) = ~CoordX+CoordY
variog<-variogram(data~CoordX+CoordY,ca20,alpha=c(0,45,90,135))
plot(variog,type='l',col="darkviolet")
plot(variogram(z~x+y,ca20,cloud=TRUE),col=4)

##########################
####### variofit #########
##########################
vario.sp<-variog(as.geodata(NewData),estimator.type= "modulus")
plot(vario.sp,type="p",col=2,pch=16,main="modulus  semivariogram ",xlab="||h||",ylab=expression(gamma(h)))
## Gaussian
x<-eyefit(vario.sp)
ols <- variofit(vario.sp, ini=x, fix.nug=FALSE, wei="equal",cov.model="gaussian")
rE1<-summary(ols)$sum.of.squares
wls <- variofit(vario.sp, ini=x, fix.nug=FALSE,cov.model="gaussian")
rE2<-summary(wls)$sum.of.squares
RE<-cbind(rE1,rE2)
RE
## exponential
x<-eyefit(vario.sp)
olsM <- variofit(vario.sp, ini=x, fix.nug=FALSE, wei="equal",cov.model="exponential")
rM1<-summary(olsM)$sum.of.squares
wlsM <- variofit(vario.sp, ini=x, fix.nug=FALSE,cov.model="exponential")
rM2<-summary(wlsM)$sum.of.squares
RM<-cbind(rM1,rM2)
RM
##########################
####### kriging ##########
##########################
ca20<-as.data.frame(NewData)
x<-ca20$CoordX
y<-ca20$CoordY
z<-ca20$data
c=cbind(x,y,z)
c=data.frame(c)
cooord=cbind(x,y)
colnames(c)=c("x","y","z")
coordinates(c) = ~x+y
cooord=as.data.frame(cooord)
cooord=expand.grid(seq(4,18,0.5),seq(3,13,0.5))
colnames(cooord)=c("x","y") 
gridded(cooord) = ~x+y
m <-vgm(0.01,"Gau",83.26,0.003)
m
pishgooee<-krige(z~x+y,c,cooord,model=m)
spplot(pishgooee["var1.pred"], main = "Universal kriging predictions")
spplot(pishgooee["var1.var"],  main = "Universal kriging variance")
pis<- krige(z~1,c, cooord, model = m)
spplot(pis["var1.pred"], main = "ordinary kriging predictions")
spplot(pis["var1.var"],  main = "ordinary kriging variance")
pis<-krige(z~1,c, cooord, model = m,beta= 0.0441576)
spplot(pis["var1.pred"], main = "simple kriging predictions")
spplot(pis["var1.var"],  main = "simple kriging variance")

##############
####cross validation
summary(pishgooee)
res <- as.data.frame(pishgooee)$residual
sqrt(mean(res^2))
mean(res)
mean(res^2/as.data.frame(pishgooee))$var1.var)
mean(NewData$data)
var(NewData$data)



##############################################################
########## krig bayes  ##########

krige.bayes(NewData, coords = cooord, data = NewData$data,
            locations = "no", borders, m, prior)

prior.control(beta.prior = c("flat", "normal", "fixed"),
              beta = NULL, beta.var.std = NULL,
              sigmasq.prior = c("reciprocal", "uniform",
                                "sc.inv.chisq", "fixed"),
              sigmasq = NULL, df.sigmasq = NULL,
              phi.prior = c("uniform", "exponential","fixed",
                            "squared.reciprocal", "reciprocal"),
              phi = NULL, phi.discrete = NULL,
              tausq.rel.prior = c("fixed", "uniform", "reciprocal"),
              tausq.rel, tausq.rel.discrete = NULL)

krigBz=krige.bayes(NewData, coords = cooord, data = NewData$data,
            model = model.control(cov.model = "gaussian",
            psill =  0.0058,range = 1, nugget = 0.9159),
            prior = prior.control(beta.prior = "normal",sigmasq.prior ="sc.inv.chisq")
plot(krigBz, type="h", tausq.rel = FALSE)


#############################
####### correlation #########
#############################
cor(wrc$Pr81,NewData$data)
cor(wrc$Pr100,NewData$data)
cor(wrc$Pr3060,NewData$data)
cor(wrc$Pr15300,NewData$data)
cor(wrc$Pr5,NewData$data)
cor(wrc$Pr10,NewData$data)
cor(wrc$Pr60,NewData$data)
cor(wrc$Pr100,NewData$data)
cor(wrc$Pr306,NewData$data)
cor(wrc$Pr816,NewData$data)


#########################
##### cokriging ########
#########################
kv.ok <- krige(data ~ 1,c, coords, m)
kv.ck <- predict.gstat(g, coords)
SpatialPointsDataFrame(coordinates(kv.ok), data=as.data.frame(kv.diff))
