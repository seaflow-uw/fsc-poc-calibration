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

for (file in file.list) {
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

write.csv(summary.table,file=paste("influx_raw_summary.csv", sep=""), row.names=FALSE)










#############################
#### 3. FSC NORMALIZATION ###
#############################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
path.to.data <- "~/Documents/DATA/Codes/fsc-poc-calibration/fsc-poc-calibration-data"
setwd(path.to.data)
summary.table <- read.csv(file=paste("influx_raw_summary.csv", sep=""))

volume <- c(rep(500.3,2), rep(500.5,2),rep(496,2),rep(500.12,2),rep(300.3,2),rep(300.22,2),rep(300.45,2), rep(300,2), rep(150.46,2), rep(200.18,2),rep(150.2,2), rep(200.1,2),rep(300.28,2),rep(302,2),rep(200.57,2),rep(199.97,2),rep(300.25,2),
            rep(300.24,2),rep(99.9,2),rep(100.13,2),rep(99.87,2),rep(99.84,2),rep(100.05,2),rep(99.98,2),rep(100.2,2),rep(100.2,2),rep(99.93,2),rep(99.9,2),rep(100.07,2),rep(100.07,2))
name <- rep(c("EHux","Licmo","Micro","Navicula","PT-632","PT-632","TAPS-1335","TAPS-3367","TW-3365","Pro 1314", "Syn 7803", "Pro AS9601", "Pro Med4","Pro Nat12A", "Syn WH8102"),each=4)

### CALCULATE CELL ABUNDANCE
summary.table$volume.uL <- volume
summary.table$abundance_cells.mL <- 1000 * summary.table$n / summary.table$volume.uL # cells / mL

# WARNING: Small size phytoplankton cultures were dilutd before counting (see Influx-notebook/harvest1.jpg, harvest2.jpg Influx-notebook/harvest3.jpg).

summary.table$abundance_cells.mL[c(37:40,45:56)] <- summary.table$abundance_cells.mL[c(37:40,45:56)] * 50 # Prochlorococcus cultures were diluted 50X
summary.table$abundance_cells.mL[c(41:44,57:60)] <- summary.table$abundance_cells.mL[c(41:44,57:60)] * 100 # Synechococcus cultures were diluted 100X
summary.table$abundance_cells.mL[9:12] <- summary.table$abundance_cells.mL[9:12] * 5 #Micromonas pusilla was dilutd 5 times


beads <- subset(summary.table, i == 'beads')
cultures <- subset(summary.table, i == 'picoeuk')
print(paste("Does rows of 'beads' and rows of 'cultures' match?",unique(cultures$file == beads$file)))
cultures$norm.fsc <- round(cultures$fsc/beads$fsc,2)
cultures$norm.chl <- round(cultures$fsc/beads$chl,2)


write.csv(cultures[,c("file","n","volume.uL","abundance_cells.mL","norm.fsc", "norm.chl")],file=paste(path.to.git.repository,"/influx-cultures.csv", sep=""), row.names=FALSE)


###########################################
### 4. MERGE FCM data with POC/PON data ###
###########################################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)
poc <- read.csv("poc-data.csv")
cultures <- read.csv("influx-cultures.csv")
cultures$Sample.ID <- c(rep("EHUX",2), rep("LICMO",2), rep("Micromonas pusilla",2),rep("Navicula transitans", 2), rep("Phaeodactylum tricornutum",4), rep("Thalassiosira pseudonana (1135)",2),
                        rep("Thalassiosira pseudonana (3367)",2), rep("TW 3365",2), rep("Prochlorococcus (1314)",2), rep("Synechococcus (7803)",2),
                        rep("Prochlorococcus (AS9601)",2), rep("Prochlorococcus (MED4)",2),rep("Prochlorococcus (NAT12A)",2),rep("Synechococcus (WH8102)",2))

poc.sd <- aggregate(poc, by=list(poc$Sample.ID), FUN=sd)
poc.mean <- aggregate(poc, by=list(poc$Sample.ID), FUN=mean)
poc.mean$Sample.ID <- poc.mean$Group.1
poc.mean$C..ug.ml.sd <- poc.sd$C..ug.ml.
poc.mean$N..ug.ml.sd <- poc.sd$N..ug.ml.

cultures.sd <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=sd)
cultures.mean <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=mean)
cultures.mean$Sample.ID <- cultures.mean$Group.1
cultures.mean$abundance_cells.mL.sd <- cultures.sd$abundance_cells.mL
cultures.mean$norm.fsc.sd <- cultures.sd$norm.fsc
cultures.mean$norm.chl.sd <- cultures.sd$norm.chl

### MERGE POC with Cell Abundance
merge <- merge(poc.mean, cultures.mean, by='Sample.ID')

### Calculate Quotas
merge$pgC.cell <- 10^6*(merge$C..ug.ml.)/(merge$abundance_cells.mL) # pgC.cell-1
merge$pgN.cell <- 10^6*(merge$N..ug.ml.)/(merge$abundance_cells.mL) # pgN.cell-1
merge$pgC.cell.sd <- merge$pgC.cell * sqrt((merge$C..ug.ml.sd/merge$C..ug.ml.)^2 + (merge$abundance_cells.mL.sd/merge$abundance_cells.mL)^2)
merge$pgN.cell.sd <- merge$pgN.cell * sqrt((merge$N..ug.ml.sd/merge$N..ug.ml.)^2 + (merge$abundance_cells.mL.sd/merge$abundance_cells.mL)^2)

### CELL QUOTAS vs NORM FSC for INFLUX
write.csv(merge[,c("Sample.ID","norm.fsc","norm.fsc.sd","norm.chl","norm.chl.sd","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file="Influx-Qc-cultures.csv", row.names=FALSE)


#### CELL QUOTAS reference
write.csv(merge[,c("Sample.ID","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file="Qc-cultures.csv", row.names=FALSE)


############################
### 5. LINEAR REGRESSION ###
############################
library(scales)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)
merge <- read.csv("Influx-Qc-cultures.csv")
mie <- read.csv("INFLUXcalibrated-mie.csv")

### WARNING !!! ###
# 1. For regression purpose, we excluded elongated cell type, specifically Phaedactylum tricornutum since forward scatter is sensitive to the cell width (not cell length), underestimating the true cell size
merge2 <- subset(merge, Sample.ID !="Phaeodactylum tricornutum")
merge2 <- merge2[order(merge2$norm.fsc),]

# linear regression
reg <- lm(pgC.cell ~ poly(norm.fsc,1,raw=T) , data=log(merge2[,c("pgC.cell","norm.fsc")],10))
summary(reg)

png("Influx-Qc-scatter.png",width=12, height=12, unit='in', res=100)

plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=NA,xlim=c(0.002,20), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)")
lines(mie$scatter, mie[,paste0('Qc_mid')], col='red3', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_upr')], col='grey', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_lwr')], col='grey', lwd=2)
with(merge2, arrows(norm.fsc, pgC.cell - pgC.cell.sd, norm.fsc, pgC.cell + pgC.cell.sd,  code = 3, length=0, col='grey', lwd=2))
with(merge2, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0,col='grey',lwd=2))
points(merge2$norm.fsc,merge2$pgC.cell,bg=alpha(.rainbow.cols(nrow(merge2)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(as.vector(merge2$Sample.ID),"Mie-based model (n = 1.38 +/- 0.3)"), pch=c(rep(21,nrow(merge2)),NA, NA), lwd=c(rep(NA,nrow(merge2)),2, NA), bty='n',
          pt.bg=alpha(.rainbow.cols(nrow(merge2)),0.5), col=c(rep(1,nrow(merge2)),'red3'), text.font=c(rep(3,nrow(merge2)),1))


dev.off()
