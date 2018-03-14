library(flowCore)
library(splancs)
library(plotrix)
library(caroline)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))


plot.cytogram <- function (evtopp, para.x = "FSC.small.stuff", para.y = "X692.40.small.stuff", ...)
{
    cols <- colorRampPalette(c("blue4", "royalblue4", "deepskyblue3",
        "seagreen3", "yellow", "orangered2", "darkred"))
    par(pty = "s")
    plot(evtopp[, c(para.x, para.y)], pch = 16, cex = 0.3,col = densCols(log10(evtopp[, c(para.x, para.y)]),
          colramp = cols), xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy',...)
          }

plot.vct.cytogram <- function (opp, para.x = "fsc_small", para.y = "chl_small", ...)
          {
        plot(opp[, c(para.x, para.y)], pch = 16, cex = 0.3, col = as.numeric(as.factor(opp$pop)),
                    xlim = c(1, 10^4), ylim = c(1, 10^4), log='xy', ...)
              legend("topleft", legend = (unique(opp$pop)), col = unique(as.numeric(as.factor(opp$pop))),
                  pch = 16, pt.cex = 0.6, bty = "n")
                    }

#####################
## 1.DOWNLOAD DATA ##
#####################

dat://cdfef982ea4032592e454c1a39b0a3855738b309d7e78ef8b2d0152adc5ffd02



###########################
## 2. BATCH FILES INFLUX ##
###########################
path.to.data <- "~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data"
setwd(path.to.data)

file.list <- dir(".", pattern = ".fcs$", recursive=T, full.names=T)


summary.table <- NULL
draw.gate <- TRUE

for (file in file.list[1:32]) {
    print(paste("processing file:",file))
#file <- file.list[2]

### read FCS
fcs <- read.FCS(file, transformation=T, emptyValue=F)
opp <- tab2df(exprs(fcs))
opp$pop <- 0


### GATING
#1. NOISE & BEADS
x <- subset(opp, pop==0)
if(draw.gate) plot.cytogram(x, "X692.40.small.stuff", "X580.30", main="NOISE & BEADS & SYN")
print("Gating Beads")
if(draw.gate) poly.beads <- getpoly(quiet=TRUE)
beads <- subset(x,inout(x[,c("X692.40.small.stuff","X580.30")],poly=poly.beads, bound=TRUE, quiet=TRUE))
opp[row.names(beads),'pop'] <- "beads"

#2. cultures
x <- subset(opp, pop==0)
if(draw.gate) {
    plot.cytogram(x, "FSC.small.stuff", "X692.40.small.stuff", main="PRO & PicoEuk")
        points(beads$FSC.small.stuff, beads$X692.40.small.stuff, col="grey", pch=16, cex=0.4)
        }

print("gating PICO")
if(draw.gate) poly.pico <- getpoly(quiet=TRUE)
pico <- subset(x,inout(x[,c("FSC.small.stuff","X692.40.small.stuff")],poly=poly.pico, bound=TRUE, quiet=TRUE))
opp[row.names(pico),'pop'] <- "picoeuk"

### SAVE PLOT
png(paste0(file,".png"),width=9, height=12, unit='in', res=100)
par(mfrow=c(2,2))
plot.vct.cytogram(opp, "FSC.small.stuff","X692.40.small.stuff")
plot.vct.cytogram(opp, "FSC.small.stuff","X580.30")
plot.vct.cytogram(opp, "X692.40.small.stuff","X580.30")
plot.vct.cytogram(opp, "SSC","X692.40.small.stuff")
dev.off()

### SUMMARY
stat.table <- NULL
for(i in unique(opp$pop)){
#print(i)
  if(i == 0) next
  p <- subset(opp, pop == i)
  n <- nrow(p)
  if(n ==0) {
      fsc <- 0
      chl <- 0
          }else{
      fsc <- round(median(p$FSC.small.stuff))
      chl <- round(median(p$X692.40.small.stuff))
      ssc <- round(median(p$SSC))
      pe <- round(median(p$X580.30))
      }
    var <- cbind(i,n,fsc,chl,ssc,pe)
    stat.table <- rbind(stat.table, var)
}

table <- data.frame(cbind(stat.table, file=basename(file)))
summary.table <- rbind(summary.table, table)

}

write.csv(summary.table,file=paste("influx-summary.csv", sep=""), row.names=FALSE)










#########################
### 3. MERGE with LOG ###
#########################
path.to.data <- "~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data/preparation"
setwd(path.to.data)
summary.table <- read.csv(file=paste("summary.csv", sep=""))

### PREPARATION
### FEB 10
#volume <- c(rep(50,4), rep(50,4), rep(50,4), rep(50,4),rep(150,4),rep(50,4),rep(100,8),rep(100,4),rep(50,4))
#name <- c("TAPS-1335", "E. hux","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Minutocellus","Navicula","TW")
### FEB 15
volume <- c(rep(50,4), rep(100,4), rep(100,4), rep(50,4),rep(50,4),rep(100,4),rep(100,4),rep(100,4),rep(50,4),rep(50,4))
name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Minutocellus","Navicula")
### FEB 17
volume <- c(rep(100,2), rep(150,2), rep(100,2), rep(100,2),rep(100,2),rep(200,2),rep(100,2),rep(100,2),rep(100,2))
name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Navicula")
### FEB 21
volume <- c(rep(50,2), rep(150,2), rep(50,2), rep(50,2),rep(50,2),rep(150,2),rep(100,2),rep(100,2),rep(50,2))
name <- c("TAPS-1335", "E. hux","TW","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus","Navicula")
### FEB 23
volume <- c(rep(50,2), rep(100,2), rep(50,2), rep(50,2),rep(100,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2),rep(50,2))
name <- c("TAPS-1335", "E. hux","TAPS-3367","Phaeo", "Licmophora","Micromonas","Minutocellus-D","Minutocellus","Navicula-D","Navicula","TW")
### FEB 23 - POC
volume <- c(rep(300.45,2), rep(300,2), rep(200.18,2), rep(200.1,2),rep(300.28,2),rep(302,2),rep(200.57,2),rep(199.97,2),rep(300.25,2),rep(300.24,2))
name <- rep(c("Navicula","PT-632","TAPS-1335","TAPS-3367","TW-3365"),each=2)
### FEB 24
volume <- rep(100,6)
name <- c("E. hux","Licmophora","Micromonas")
### APRIL 25
volume <- rep(100,12)
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")
### APRIL 27
volume <- rep(50,12)
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")
### MAY 2
volume <- c(rep(100,8), 140,140, 100,100)
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")
### MAY 3
volume <- c(75,75,rep(50,10))
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")
### MAY 4
volume <- rep(50,12)
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")
### MAY 5
volume <- rep(50,12)
name <- c("Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102")






### HARVEST
path.to.data  <- "~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data/harvest"

setwd(path.to.data)

summary.table <- read.csv(file=paste("influx-summary.csv", sep=""))

volume <- c(rep(500.3,2), rep(500.5,2),rep(496,2),rep(500.12,2),rep(300.3,2),rep(300.22,2),rep(300.45,2), rep(300,2), rep(150.46,2), rep(200.18,2),rep(150.2,2), rep(200.1,2),rep(300.28,2),rep(302,2),rep(200.57,2),rep(199.97,2),rep(300.25,2),
            rep(300.24,2),rep(99.9,2),rep(100.13,2),rep(99.87,2),rep(99.84,2),rep(100.05,2),rep(99.98,2),rep(100.2,2),rep(100.2,2),rep(99.93,2),rep(99.9,2),rep(100.07,2),rep(100.07,2))
name <- rep(c("EHux","Licmo","Micro","Navicula","PT-632","PT-632","TAPS-1335","TAPS-3367","TW-3365","Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102"),each=4)




### CALCULATE CELL ABUNDANCE
summary.table$volume.uL <- volume
summary.table$abundance_cells.mL <- 1000 * summary.table$n / summary.table$volume.uL # cells / mL

# WARNING: CYANOBAACTERIA cultures were diluted 100X before counting.
summary.table$abundance_cells.mL[37:60] <- summary.table$abundance_cells.mL[37:60] * 100

write.csv(summary.table,file=paste("influx-stats.csv", sep=""), row.names=FALSE)









#############################
#### 4. FSC NORMALIZATION ###
#############################

path <- paste("~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data/harvest/influx-stats.csv")
summary <- read.csv(file=path)
beads <- subset(summary, i == 'beads')
cultures <- subset(summary, i == 'picoeuk')

print(paste("Does rows of beads and cultures match?",unique(cultures$file == beads$file)))
cultures$norm.fsc <- round(cultures$fsc/beads$fsc,2)
cultures$norm.chl <- round(cultures$fsc/beads$chl,2)


write.csv(cultures[,c("file","n","volume.uL","abundance_cells.mL","norm.fsc", "norm.chl")],file=paste("~/Documents/DATA/Codes/fsc-poc-calibration/influx-cultures.csv", sep=""), row.names=FALSE)


###########################################
### 5. MERGE FCM data with POC/PON data ###
###########################################
poc <- read.csv("~/Documents/DATA/Codes/fsc-poc-calibration/poc-data.csv")
cultures <- read.csv("~/Documents/DATA/Codes/fsc-poc-calibration/influx-cultures.csv")

#cultures$Sample.ID <- c(rep("EHUX",2), rep("LICMO",2), rep("MICRO",2),rep("NAV", 2), rep("PT 632",4), rep("TAPS 1135",2), rep("TAPS 3367",2), rep("TW 3365",2), rep("1314",2), rep("7803",2), rep("AS9601",2), rep("MED4",2),rep("NAT12A",2),rep("WH8102",2))
#write.csv(cultures,file=paste("~/Documents/DATA/Codes/fsc-poc-calibration/influx-cultures.csv", sep=""), row.names=FALSE)

poc.sd <- aggregate(poc, by=list(poc$Sample.ID), FUN=sd)
poc <- aggregate(poc, by=list(poc$Sample.ID), FUN=mean)
poc$Sample.ID <- poc$Group.1
poc$Date <- poc$Group.2
poc$C..ug.ml.sd <- poc.sd$C..ug.ml.
poc$N..ug.ml.sd <- poc.sd$N..ug.ml.

cultures.sd <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=sd)
cultures <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=mean)
cultures$Sample.ID <- cultures$Group.1
cultures$abundance_cells.mL.sd <- cultures.sd$abundance_cells.mL
cultures$norm.fsc.sd <- cultures.sd$norm.fsc
cultures$norm.chl.sd <- cultures.sd$norm.chl

### MERGE POC with Cell Abundance
merge <- merge(poc, cultures, by='Sample.ID')

### Calculate Quotas
merge$pgC.cell <- 10^6*(merge$C..ug.ml.)/(merge$abundance_cells.mL) # pgC.cell-1
merge$pgN.cell <- 10^6*(merge$N..ug.ml.)/(merge$abundance_cells.mL) # pgN.cell-1
merge$pgC.cell.sd <- sqrt((merge$C..ug.ml.sd/merge$C..ug.ml.)^2 + (merge$abundance_cells.mL.sd/merge$abundance_cells.mL)^2)
merge$pgN.cell.sd <- sqrt((merge$N..ug.ml.sd/merge$N..ug.ml.)^2 + (merge$abundance_cells.mL.sd/merge$abundance_cells.mL)^2)

#write.csv(merge[,c("Sample.ID","norm.fsc","norm.fsc.sd","norm.chl","norm.chl.sd","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file="~/Documents/DATA/Codes/fsc-poc-calibration/Influx-Qc-cultures.csv", row.names=FALSE)

#### CELL QUOTAS reference
#write.csv(merge[,c("Sample.ID","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file="~/Documents/DATA/Codes/fsc-poc-calibration/Qc-cultures.csv", row.names=FALSE)


############################
### 6. LINEAR REGRESSION ###
############################
library(lmodel2)

merge <- read.csv(paste0("~/Documents/DATA/Codes/fsc-poc-calibration/Influx-Qc-cultures.csv"))

### WARNING !!! ###
# 1. We excluded all samples collected on 2/23/17 due to abnormal low values of POC (hypothesis: filters were porus and didn't retain all cells)
# 2. We excluded data from Micromonas due to the very low biomass collected on filters (barely above blank).
merge2 <- subset(merge, Sample.ID !="PT 632" & Sample.ID !="EHUX" & Sample.ID !="LICMO")


# linear regression
reg <- lmodel2(pgC.cell ~ norm.fsc, data=log(merge2[,c("pgC.cell","norm.fsc")],10))
#reg2 <- lmodel2(pgC.cell ~ norm.fsc, data=log(merge[,c("pgC.cell","norm.fsc")],10))



png("~/Documents/DATA/Codes/fsc-poc-calibration/Influx-Qc-scatter.png",width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(merge$norm.fsc,merge$pgC.cell, log='xy', yaxt='n',cex=2,bg='grey', pch=21,ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main='Influx')
points(merge2$norm.fsc,merge2$pgC.cell, cex=2,bg='red3', pch=21)
with(merge, arrows(norm.fsc, pgC.cell - merge$pgC.cell.sd, norm.fsc, pgC.cell + merge$pgC.cell.sd,  code = 3, length=0))
with(merge, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0))
axis(2, at=c(0.1,1,10,100,1000), labels=c(0.1,1,10,100,1000))
par(new=T)
plot(log10(merge$norm.fsc), log10(merge$pgC.cell), yaxt='n',xaxt='n',xlab=NA, ylab=NA,pch=NA, bty='n')
abline(b=reg$regression.results$Slope[1], a=reg$regression.results$Intercept[1], col=2,lwd=2)
abline(b=reg$confidence.intervals[1,4], a=reg$confidence.intervals[1,2], col='grey',lwd=2)
abline(b=reg$confidence.intervals[1,5], a=reg$confidence.intervals[1,3], col='grey',lwd=2)
#text(log(merge$norm.fsc,10), log(merge$pgC.cell,10), labels=merge$Sample.ID)
legend("topleft", legend=bquote(paste("Qc=",.(round(reg$regression.results$Intercept[1],3)),"(scatter"^{.(round(reg$regression.results$Slope[1],3))},")")), bty='n',cex=2)

dev.off()
