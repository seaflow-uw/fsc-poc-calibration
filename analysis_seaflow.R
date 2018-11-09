library(flowCore)
library(splancs)
library(plotrix)
library(caroline)
library(popcycle)
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

file.list <- dir("740_Small_cyano", pattern = "00-00$", recursive=T, full.names=T); inst <- 740 # small
file.list <- dir("740_Large_picoeuks", pattern = "00-00$", recursive=T, full.names=T); inst <- 740 # large
file.list <- dir(".", pattern = "07-00$", recursive=T, full.names=T); inst <- 751 # small
file.list <- dir(".", pattern = "08-00$", recursive=T, full.names=T); inst <- 751 # large

# gating
summary.table <- NULL
draw.gate <- TRUE

    for (file in file.list) {
        print(paste("processing file:",file))
    #file <- file.list[11]

    ### read EVT and get OPP
    evt <- readSeaflow(file, transform=F)
    evt <- subset(evt, fsc_small > 2)
    evt <- evt[round(seq(1,nrow(evt), length.out=100000)),]
    ip <- inflection.point(evt)
    filter.params <- create.filter.params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2, width=10000)
    # plot.filter.cytogram(evt, filter.params[2,])

    evt <- readSeaflow(file, transform=T)
    opp <- filter.notch(evt, filter.params[2,])$opp
    opp$pop <- 0

    ### GATING
    #1. BEADS
    x <- subset(opp, pop==0)
    if(draw.gate) plot.cytogram(x, "fsc_small","pe", main="NOISE & BEADS & SYN")
    print("Gating Beads")
    if(draw.gate) poly.beads <- getpoly(quiet=TRUE)
    beads <- subset(x,inout(x[,c("fsc_small","pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))
    opp[row.names(beads),'pop'] <- "beads"

    x <- subset(opp, pop==0)
    plot.cytogram(x, "fsc_small", "chl_small", main="Pico")
    points(beads$fsc_small, beads$chl_small, cex=0.3, col=2, pch=16)
    print("gating PICO")
    if(draw.gate) poly.pico <- getpoly(quiet=TRUE)
    pico <- subset(x,inout(x[,c("fsc_small","chl_small")],poly=poly.pico, bound=TRUE, quiet=TRUE))
    opp[row.names(pico),'pop'] <- "picoeuk"

    ### SAVE PLOT
    png(paste0(file,".png"),width=9, height=12, unit='in', res=100)
    par(mfrow=c(2,2))
    plot.vct.cytogram(opp, "fsc_small","chl_small")
    plot.vct.cytogram(opp, "fsc_small","pe")
    plot.vct.cytogram(opp, "chl_small","pe")
    #plot.vct.cytogram(evt, "fsc_small","fsc_perp")
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

write.csv(summary.table,file=paste0(unique(dirname(file.list)),"/",inst,"_raw_summary.csv", sep=""), row.names=FALSE)





#############################
#### 3. FSC NORMALIZATION ###
#############################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"


for(inst in c(740,751)){
    print(inst)
    list <- dir(".", pattern = paste0(inst,"_raw_summary.csv"), recursive=T, full.names=T)

    DF <- NULL
    for(file in list){
        df <- read.csv(file)
        beads <- subset(df, i == 'beads')
        cultures <- subset(df, i == 'picoeuk')

        print(paste("Does rows of beads and cultures match?",unique(cultures$file == beads$file)))
        cultures$norm.fsc <- round(cultures$fsc/mean(beads$fsc),2)
        cultures$norm.chl <- round(cultures$fsc/mean(beads$chl),2)
        DF <- rbind(DF, cultures)
    }
    time <- as.POSIXct(DF$file, format = "%FT%H-%M-%S", tz = "GMT")
    DF[order(time),]
    write.csv(DF,file=paste0(path.to.git.repository,"/",inst,"-cultures.csv", sep=""), row.names=FALSE)

}

###########################################
### 4. MERGE FCM data with POC/PON data ###
###########################################
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)
poc <- read.csv("Qc-cultures.csv")

for(inst in c(740,751)){

    #inst <- 740
    cultures <- read.csv(paste0(inst,"-cultures.csv"))

    if(inst == 740){ cultures$Sample.ID <- c(rep("TW 3365",2),rep("NAV",2), rep("TAPS 1135",2), rep("TAPS 3367",2), rep("PT 632",2),rep("MICRO",2),
                                              rep("AS9601",2), rep("1314",2),rep("MED4",2),rep("NAT12A",2),rep("WH8102",2),rep("7803",2))
              }
    if(inst == 751){ cultures$Sample.ID <- c(rep("TW 3365",2),rep("NAV",2), rep("TAPS 1135",2), rep("TAPS 3367",2),rep("PT 632",2), rep("MICRO",2),
                                              rep("MED4",2), rep("AS9601",2),rep("1314",2),rep("NAT12A",2), rep("WH8102",2), rep("7803",2))
              }

    cultures.sd <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=sd)
    cultures.mean <- aggregate(cultures, by=list(cultures$Sample.ID), FUN=mean)
    cultures.mean$Sample.ID <- cultures.mean$Group.1
    cultures.mean$norm.fsc.sd <- cultures.sd$norm.fsc
    cultures.mean$norm.chl.sd <- cultures.sd$norm.chl


    ### MERGE POC with Cell Abundance
    merge <- merge(poc, cultures.mean, by='Sample.ID')

    write.csv(merge[,c("Sample.ID","norm.fsc","norm.fsc.sd","norm.chl","norm.chl.sd","abundance_cells.mL","abundance_cells.mL.sd","pgC.cell","pgN.cell","pgC.cell.sd","pgN.cell.sd")],file=paste0(inst,"-Qc-cultures.csv"), row.names=FALSE)

}



############################
### 5. LINEAR REGRESSION ###
############################
library(scales)
library(viridis)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

### CHOOSE SeaFlow serial number
inst <- 740
inst <- 751


merge <- read.csv(paste0(inst,"-Qc-cultures.csv"))
merge2 <- subset(merge, Sample.ID !="PT 632" )#& Sample.ID !="TAPS 3367" & Sample.ID !="TAPS 1135" & Sample.ID !="NAV")
merge2 <- merge2[order(merge2$norm.fsc),]

mie <- read.csv("calibrated-mie.csv")

png(paste0(inst,"-Qc-scatter.png"),width=12, height=12, unit='in', res=100)

par(mfrow=c(1,1), pty='s',cex=1.4)
plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=NA,xlim=c(0.002,10), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste("#",inst))
with(merge2, arrows(norm.fsc, pgC.cell - pgC.cell.sd, norm.fsc, pgC.cell + pgC.cell.sd,  code = 3, length=0, col='grey', lwd=2))
with(merge2, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0,col='grey',lwd=2))
lines(mie$scatter, mie[,paste0("Qc_",inst,"_mid")], col='red3', lwd=2)
lines(mie$scatter, mie[,paste0("Qc_",inst,"_upr")], col='grey', lwd=2)
lines(mie$scatter, mie[,paste0("Qc_",inst,"_lwr")], col='grey', lwd=2)
points(merge2$norm.fsc,merge2$pgC.cell,bg=alpha(viridis(nrow(merge2)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(as.vector(merge2$Sample.ID),"Mie-based model (index refraction = 1.031 +/- 0.014)"), pch=c(rep(21,nrow(merge2)),NA), lwd=c(rep(NA,nrow(merge2)),2), bty='n',
          pt.bg=alpha(viridis(nrow(merge2)),0.5), col=c(rep(1,nrow(merge2)),'red3'), text.font=c(rep(3,nrow(merge2)),1))

dev.off()
