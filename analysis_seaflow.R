library(flowCore)
library(splancs)
library(plotrix)
library(caroline)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

#####################
## 1.DOWNLOAD DATA ##
#####################

dat://cdfef982ea4032592e454c1a39b0a3855738b309d7e78ef8b2d0152adc5ffd02


############################
## 2. BATCH FILES SEAFLOW ##
############################
library(popcycle)
path.to.data <- "~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data"
setwd(path.to.data)

file.list <- dir(".", pattern = "00-00$", recursive=T, full.names=T); inst <- 740
file.list <- dir(".", pattern = "07-00$", recursive=T, full.names=T); inst <- 751 # small
file.list <- dir(".", pattern = "08-00$", recursive=T, full.names=T); inst <- 751 # large


summary.table <- NULL
draw.gate <- TRUE

for (file in file.list) {
    print(paste("processing file:",file))
#file <- file.list[20]

### read EVT
evt <- readSeaflow(file, transform=T)
evt <- evt[evt$fsc_small > 1 | evt$D1 > 1 | evt$D2 > 1, ]
evt <- subset(evt, D2 < D1 + 2500 & D1 < D2 + 2500)

evt$pop <- 0

### GATING
#1. BEADS
x <- subset(evt, fsc_small > 2 & pe > 5)
if(draw.gate) plot.cytogram(x, "fsc_small","pe", main="NOISE & BEADS & SYN")
print("Gating Beads")
if(draw.gate) poly.beads <- getpoly(quiet=TRUE)
beads <- subset(x,inout(x[,c("fsc_small","pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
evt[row.names(beads),'pop'] <- "beads"


#2. cultures
x <- subset(evt, pop==0 & fsc_small > 2)
plot.cytogram(x, "fsc_small","chl_small", main="Pico")
points(beads$fsc_small, beads$chl_small, cex=0.3, col=2, pch=16)
print("gating PICO")
if(draw.gate) poly.pico <- getpoly(quiet=TRUE)
pico <- subset(x,inout(x[,c("fsc_small","chl_small")],poly=poly.pico, bound=TRUE, quiet=TRUE))
evt[row.names(pico),'pop'] <- "picoeuk"

### SAVE PLOT
png(paste0(file,".png"),width=9, height=12, unit='in', res=100)
par(mfrow=c(2,2))
plot.vct.cytogram(evt, "fsc_small","chl_small")
plot.vct.cytogram(evt, "fsc_small","pe")
plot.vct.cytogram(evt, "chl_small","pe")
#plot.vct.cytogram(evt, "fsc_small","fsc_perp")
dev.off()

### SUMMARY
stat.table <- NULL
for(i in unique(evt$pop)){
#print(i)
  if(i == 0) next
  p <- subset(evt, pop == i)
  n <- nrow(p)
  if(n ==0) {
      fsc <- 0
      chl <- 0
          }else{
      fsc <- round(median(p$fsc_small))
      chl <- round(median(p$chl_small))
      pe <- round(median(p$pe))
      }
    var <- cbind(i,n,fsc,chl,pe)
    stat.table <- rbind(stat.table, var)
}

table <- data.frame(cbind(stat.table, file=basename(file)))
summary.table <- rbind(summary.table, table)

}

write.csv(summary.table,file=paste0(inst, "-summary2.csv", sep=""), row.names=FALSE)




#############################
#### 3. FSC NORMALIZATION ###
#############################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

### CHOOSE SeaFlow serial number
inst <- 740
inst <- 751

summary <- read.csv(file=paste0(inst,"-summary.csv", sep=""))
beads <- subset(summary, i == 'beads')
cultures <- subset(summary, i == 'picoeuk')

print(paste("Does rows of beads and cultures match?",unique(cultures$file == beads$file)))
cultures$norm.fsc <- round(cultures$fsc/beads$fsc,2)
cultures$norm.chl <- round(cultures$fsc/beads$chl,2)


write.csv(cultures,file=paste0(inst,"-cultures.csv", sep=""), row.names=FALSE)



###########################################
### 4. MERGE FCM data with POC/PON data ###
###########################################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

### CHOOSE SeaFlow serial number
inst <- 740
inst <- 751


poc <- read.csv("Qc-cultures.csv")
cultures <- read.csv(paste0(inst,"-cultures.csv"))

if(inst == 740){ cultures$Sample.ID <- c(rep("TW 3365",2),"-" ,rep("NAV",2), "-" ,rep("TAPS 1135",2), rep("-",2),rep("TAPS 3367",2), "-" ,rep("PT 632",2), "-" ,rep("MICRO",2),"-" ,
                                          rep("MED4",2),rep("-",3), rep("AS9601",2),"-" , rep("1314",2),"-" ,rep("MED4",2),"-" ,rep("NAT12A",2),"-" ,rep("WH8102",2), "-" ,rep("7803",2),"-")
          }
if(inst == 751){ cultures$Sample.ID <- c(rep("MED4",2), "-", rep("AS9601",2),"-" , rep("1314",2),"-" ,rep("NAT12A",2),"-" ,rep("WH8102",2), "-" ,rep("7803",2),"-",
                                          rep("TW 3365",2),"-" ,rep("NAV",2), "-" ,rep("TAPS 1135",2), rep("-",2),rep("TAPS 3367",2), "-" ,rep("PT 632",2), rep("-",3),rep("MICRO",2),"-")
          }

cultures.sd <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=sd)
cultures <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=mean)
cultures$Sample.ID <- cultures$Group.1
cultures$norm.fsc.sd <- cultures.sd$norm.fsc
cultures$norm.chl.sd <- cultures.sd$norm.chl


### MERGE POC with Cell Abundance
merge <- merge(poc, cultures, by='Sample.ID')

write.csv(merge[,c("Sample.ID","norm.fsc","norm.fsc.sd","norm.chl","norm.chl.sd","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file=paste0(inst,"-Qc-cultures.csv"), row.names=FALSE)


############################
### 5. LINEAR REGRESSION ###
############################
library(lmodel2)
library(scales)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

### CHOOSE SeaFlow serial number
inst <- 740
inst <- 751

merge <- read.csv(paste0(inst,"-Qc-cultures.csv"))


merge2 <- subset(merge, Sample.ID !="PT 632")

# linear regression type II
reg <- lmodel2(pgC.cell ~ norm.fsc, data=log(merge2[,c("pgC.cell","norm.fsc")],10))


png(paste0(inst,"-Qc-scatter.png"),width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n',cex=2,bg=alpha(.rainbow.cols(nrow(merge2)),0.5), pch=21, ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste0("SeaFlow #",inst))
with(merge2, arrows(norm.fsc, pgC.cell - merge2$pgC.cell.sd, norm.fsc, pgC.cell + merge2$pgC.cell.sd,  code = 3, length=0))
with(merge2, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0))
axis(2, at=c(0.1,1,10,100,1000), labels=c(0.1,1,10,100,1000))
par(new=T)
plot(log10(merge2$norm.fsc), log10(merge2$pgC.cell), yaxt='n',xaxt='n',xlab=NA, ylab=NA,pch=NA, bty='n')
abline(b=reg$regression.results$Slope[1], a=reg$regression.results$Intercept[1], col=2,lwd=2)
abline(b=reg$confidence.intervals[1,4], a=reg$confidence.intervals[1,2], col='grey',lwd=2)
abline(b=reg$confidence.intervals[1,5], a=reg$confidence.intervals[1,3], col='grey',lwd=2)
#text(log(merge$norm.fsc,10), log(merge$pgC.cell,10), labels=merge$Sample.ID)
legend("topleft", legend=bquote(paste("Qc=",.(round(10^reg$regression.results$Intercept[1],3)),"(scatter"^{.(round(reg$regression.results$Slope[1],3))},")")), bty='n',cex=2)
legend("bottomright", legend=merge2$Sample.ID, pt.bg=alpha(.rainbow.cols(nrow(merge2)),0.5),pch=21,bty='n')

dev.off()


##################
### 6. SUMMARY ###
##################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

merge1 <- read.csv("740-Qc-cultures.csv")
merge2 <- read.csv("751-Qc-cultures.csv")
merge1 <- subset(merge1, Sample.ID !="PT 632")
merge2 <- subset(merge2, Sample.ID !="PT 632")
  reg <- lmodel2(pgC.cell ~ norm.fsc, data=log(merge1[,c("pgC.cell","norm.fsc")],10))
df1 <- data.frame(inst=740, expo=reg$regression.results$Slope[1],expo_97.5=reg$confidence.intervals[1,5],expo_2.5=reg$confidence.intervals[1,4],coeff=10^reg$regression.results$Intercept[1], coeff_97.5=10^reg$confidence.intervals[1,3],coeff_2.5=10^reg$confidence.intervals[1,2])
  reg <- lmodel2(pgC.cell ~ norm.fsc, data=log(merge2[,c("pgC.cell","norm.fsc")],10))
df2 <- data.frame(inst=751, expo=reg$regression.results$Slope[1],expo_97.5=reg$confidence.intervals[1,5],expo_2.5=reg$confidence.intervals[1,4],coeff=10^reg$regression.results$Intercept[1], coeff_97.5=10^reg$confidence.intervals[1,3],coeff_2.5=10^reg$confidence.intervals[1,2])

df <- rbind(df1,df2)

write.csv(df,file=paste0("seaflow_qc-calibration.csv"), row.names=FALSE)
