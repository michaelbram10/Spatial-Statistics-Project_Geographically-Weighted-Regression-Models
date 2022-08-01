##### TUGAS 3: GWR (Geographically Weighted Regression) #####

### Import Packages ###
library(GWmodel)
library(car)
library(lmtest)
library(maptools)
library(RColorBrewer)
library(rgdal)
library(ggmap)
library(tidyr)
library(tmap)
library(rgeos)
library(foreign)
library(sp)
library(lattice)
library(classInt)
library(e1071)
library(shapefiles)

### Analisis Regresi Global ###
### Import Data ###
library(readxl)
Data <- read_excel("ipm_jabar.xlsx", 
                   sheet = "Data GWR digeser")
View(Data)
head(Data)

### Describe Data ###
str(Data)
dim(Data)
names(Data)

### Define Variabel ###
Y = Data$IPM
X1 = Data$`Pengeluaran Per Kapita (Ribu Rupiah/Orang/Tahun)`
X2 = Data$`Harapan Lama Sekolah`
X3 = Data$`Persentase Penduduk Miskin`

### Model Regresi Global ###
mod <- lm(Y~X1+X2+X3, data=Data)
summary(mod)
names(mod)
#er.ols
resid <- data.frame(ID = c(0:26), mod$residuals)
resid
data.frame(resid$mod.residuals)

### Skor AIC ###
AIC(mod)

##UJI ASUMSI RESIDUAL
#Menampilkan plot
par(mfrow=c(2,2)) 
plot(mod)

#1. UJI ASUMSI NORMALITAS (Shapiro Wilk Test)
#Catatan : plot ambil dari Normal Q-Q

#Menampilkan  Normal Q-Q Plot
plot(mod, 2,
     main = "Uji Asumsi Normalitas",
     col = "red", pch = 16, cex = 1.2)

#Cara 1: Uji Kolmogorov-Smirnov (pake ini aja)
mod.er=residuals(mod)
ks.test(rnorm(35, mean=0),mod.er) 

#Menampilkan Normal Q-Q Plot
qqnorm(mod.er, ylab="Residuals", xlab="Normal Scores")
qqline(mod.er)

#Cara 2: Uji Shapiro Wilk
shapiro.test(residuals(mod))

#2. UJI ASUMSI HETEROSKESDATISITAS (Breusch Pagan Test)
#Catatan : plot ambil dari Residual vs Leverage atau Scale Location

#Menampilkan plot
plot(mod,5,
     main = "Uji Asumsi Homogenitas",
     col = "red", pch = 16, cex = 1.2)

#Cara 1: Uji Breusch Pagan (pake ini aja)
library("quantmod")
library("lmtest")
bptest(mod)
spreadLevelPlot(mod)

#Cara 2:
library("quantmod")
library("lmtest")
bptest(mod, varformula = NULL, studentize = TRUE)

#3. UJI ASUMSI MULTIKOLONIERITAS
vif(mod)

#4. UJI AUTOKORELASI (Durbin Watson Test)
#Catatan : plot ambil dari Residuals vs Fitted

plot(mod, 4,
     main = "Uji Asumsi AutokORELASI",
     col = "red", pch = 16, cex = 1.2)

#Cara 1: Uji Durbin Watson (pake ini aja)
dwtest(Y~X1+X2+X3, data=Data)

#Cara 2:
durbinWatsonTest(mod)

#Cara 3:
dwtest(mod,alternative="two.sided")

#ga kepake
x <- readOGR("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 3 (GWR)/PETA_JAWA_BARAT/petajawabaratbeneran.shp")
x <- x[x$PROVINSI %in% c("JAWA BARAT"),] #ga kepake

### Mendapatkan matriks jarak antar lokasi observasi ###
Data.spdf <- SpatialPointsDataFrame(Data[,3:4],Data) #3 dan 4 adalah kolom Lintang Bujur
head(Data.spdf)
m.jarak <- gw.dist(dp.locat=coordinates(Data.spdf))
m.jarak

### Metode Pseudo-stepwise ###
#Untuk model GWR dengan fungsi Kernel Gaussian
DeVar<-"Y"
InDeVars<-c("X1", "X2", "X3")
model.sel <- model.selection.gwr(DeVar, InDeVars, data=Data.spdf, kernel="gaussian", adaptive=FALSE, bw=1, approach="CV", dMat=m.jarak)
model.sel
sorted.models<-model.sort.gwr(model.sel,numVars=length(InDeVars), ruler.vector=model.sel[[2]][,2])
sorted.models
model.list <- sorted.models[[1]]
model.list
model.view.gwr(DeVar, InDeVars, model.list=model.list)
plot(sorted.models[[2]][,2],col="black", pch=20, lty=5, 
     main = "Alternative view of GWR model selection procedure", 
     ylab = "CV", 
     xlab = "Model number", type = "b")

### Untuk model GWR dengan fungsi Kernel Bisquare ###
DeVar <- "Y"
InDeVars<-c("X1", "X2", "X3")
model.sel<-model.selection.gwr(DeVar, InDeVars, data=Data.spdf, kernel="gaussian", adaptive=FALSE, bw=1, approach="CV", dMat=m.jarak)
model.sel
sorted.models<-model.sort.gwr(model.sel,numVars=length(InDeVars), ruler.vector=model.sel[[2]][,2])
sorted.models
model.list<-sorted.models[[1]]
model.list
model.view.gwr(DeVar, InDeVars, model.list=model.list)
plot(sorted.models[[2]][,2],col="black", pch=20, lty=5, 
     main = "Alternative view of GWR model selection procedure", 
     ylab = "CV", 
     xlab = "Model number", type = "b")

### Menentukan bandwidth optimum menggunakan CV (Cross Validation) ###
#Untuk model GWR dengan fungsi Kernel Gaussian
bw.mod.gwr.gauss<-bw.gwr(Y~X1+X2+X3 , data=Data.spdf, 
                         approach="CV", 
                         kernel="gaussian", 
                         adaptive=F)
bw.mod.gwr.gauss

#Untuk model GWR dengan fungsi Kernel Bisquare
bw.mod.gwr.bisq<-bw.gwr(Y~X1+X2+X3, data=Data.spdf, 
                       approach="CV", 
                       kernel="bisquare", 
                       adaptive=F)
bw.mod.gwr.bisq

### Melakukan estimasi parameter model GWR (Pemilihan Model Terbaik) ###
#Untuk model GWR dengan fungsi Kernel Gaussian
gwr.gaus<-gwr.basic(Y~X1+X2+X3, data=Data.spdf, 
                    bw=bw.mod.gwr.gauss, 
                    kernel="gaussian", 
                    adaptive=F)
print(gwr.gaus)
summary(gwr.gaus)

### Untuk model GWR dengan fungsi Kernel Bisquare ###
gwr.bisq<-gwr.basic(Y~X1+X2+X3, data=Data.spdf, 
                    bw=bw.mod.gwr.bisq, 
                    kernel="bisquare", 
                    adaptive=F)
print(gwr.bisq)
summary(gwr.bisq)

### Data untuk Peta2an ###
#Catatan: yang dipake gwr.gaus karena memiliki AIC terkecil dan R^2 terbesar (model terbaik)
df.int <- data.frame(ID = c(0:26), gwr.gaus$SDF$Intercept)
df.x1 <- data.frame(ID = c(0:26), gwr.gaus$SDF$X1)
df.x2 <- data.frame(ID = c(0:26), gwr.gaus$SDF$X2)
df.x3 <- data.frame(ID = c(0:26), gwr.gaus$SDF$X3)
df.res <- data.frame(ID = c(0:26), gwr.gaus$SDF$residual)
df.studres <- data.frame(ID = c(0:26), gwr.gaus$SDF$Stud_residual)

### Menggabungkan Data Peta2an ###
library(tidyverse)
df.list <- list(df.int, df.x1, df.x2, df.x3, df.res, df.studres)
dats <- data.frame(df.list %>% reduce(full_join), by = 'ID')
head(dats)

resid <- data.frame(ID = c(0:26), mod$residuals)
dats2 <- merge(dats, resid, by = 'ID')
View(dats2)

#Jadiin Data Txt
write.table(dats2, file = "datapeta.txt")

### Melakukan uji signifikansi parameter GWR ###
gwr.t.adjust(gwr.gaus) #pake ini karena model terbaik
gwr.t.adjust(gwr.bisq) 

sign <- gwr.t.adjust(gwr.bisq)
df.beta <- data.frame(sign$results$t)
df.beta

################################################################################
##### Peta #####
### Import Packages ###
library(maptools)
library(RColorBrewer)
library(sp)
library(foreign)
library(lattice)
library(rgdal)
library(classInt)
library(class)
library(e1071)
library(shapefiles)

### Import Shapefile ###
Jabar <- readOGR("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 3 (GWR)/petajawabaratminuswaduk/petajawabaratminuswaduk.shp")
plot(Jabar, main="Peta Jawa Barat")
head(Jabar)

### Import Data Peta ###
library(readxl)
datapeta <- read_excel("datapeta.xlsx", sheet = "datapeta digeser")
View(datapeta)
head(datapeta)

### Menampilkan Plot Peta Jawa Barat ###
plot(Jabar, density=16, col="grey", axes=T, cex.axis=.75)
title(main="Peta Jawa Barat", sub="Created with Rstudio",font.sub=2)
title(xlab="Longitude", ylab='Latitude',cex.lab=.75,line=2.25)
text(coordinates(Jabar), labels=Jabar$ADM2_EN, cex=.5)
#ADM2_EN adalah nama kabupaten/kota Jawa Barat di Shapefile

### Merubah variabel pada data Shapefile menjadi variabel numerik ###
Jabar@data$ID <- as.numeric(row.names(Jabar@data))
Jabar@data$ID
Jabar@data$row <- as.numeric(row.names(Jabar@data))
Jabar@data$row 

str(Jabar@data)
length(Jabar@data$ADM2_EN) #ADM2_EN adalah nama kabupaten/kota Jawa Barat di Shapefile
length(Jabar@data$ID)

### Menggabungkan Shapefile dengan Data di Excel ###
temp <- merge(Jabar@data, Data, by="ID", all.Jabar=T, sort=F)
temp
Jabar@data <- temp[order(temp$row),]
Jabar@data

### Menggabungkan Shapefile dan Data di Excel dengan Data Peta ###
temp1 <- merge(Jabar@data, datapeta, by="ID", all.Jabar=T, sort=F)
temp1
Jabar@data <- temp1[order(temp1$row),]
Jabar@data

### Layout ###
Data <- read_excel("ipm_jabar.xlsx")
Data.spdf <- SpatialPointsDataFrame(Jabar@data[,19:18], Jabar@data) #19 dan 18 adalah kolom Lintang Bujur
head(Data.spdf)

kabkot <- list('sp.pointLabel', Data.spdf, label=Jabar@data$ADM2_EN, cex=0.5, col='black')
titik <- list('sp.points', Data.spdf, pch=16, cex=.8, col='blue')

#Plot Peta Residual Regresi Global
spplot(Jabar, 'er_ols', col.regions = brewer.pal(9,"Oranges"), 
       main="Plot Peta Residual Regresi Global", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=8)   #er_ols dari data Jabar@data

#Plot Peta Residual Model GWR dengan Fungsi Pembobot Kernel Gaussian
spplot(Jabar, 'gwr_res', col.regions = brewer.pal(9,"PuRd"), 
       main="Plot Peta Residual Model GWR dengan Fungsi Pembobot Kernel Gaussian", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=8)

#Plot Peta IPM
spplot(Jabar, "y", col.regions=brewer.pal(9,"YlOrRd"), cuts=4,  
       main="IPM Jawa Barat 2021", 
       scales = list(draw = TRUE), sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude") 
       #at=seq(0, 100, 5))

#Plot Region (Plot Peubah Signifikan)
spplot(Jabar, "gwr_t", col.regions=brewer.pal(5,"PuBu"),
       main="Peta Sebaran Peubah Signifikan", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=3)

#Plot variabel X1
spplot(Jabar, "gwr_x1", col.regions=brewer.pal(9,"Reds"), 
       main="Sebaran Koefisien X1", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=8)
       #at=seq(-.006, .006, .001))

#Plot variabel X2
spplot(Jabar, "gwr_x2",col.regions=brewer.pal(9,"Greens"), 
       main="Sebaran Koefisien X2", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=8)

#Plot variabel X3
spplot(Jabar, "gwr_x3",col.regions=brewer.pal(9,"Purples"), 
       main="Sebaran Koefisien X3", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       cuts=8) 
       #at=seq(-.005, .005, .001))