
library(readODS)
library(corrplot)
library(RColorBrewer)
library(pals)
library(agricolae)
library('plot.matrix')
palette_people<-brewer.pal(3, "Dark2")

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

data<-read_ods("1HRV-MIRJAM_LM_0.2.ods")

#plotting all correlations
corrmatrix<-cor(data[,c(-2, -14, -16,-17,-18,-22,-24, -25, -26)])
corrplot.mixed(corr = corrmatrix, lower = "ellipse", upper = "number", tl.col = "black")
corrplot(corr = corrmatrix, type = "upper", tl.pos = "td", diag = FALSE)



# calculating p value of correlations
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(corrmatrix)


png("Correlations.png", width = 2500, height = 2500, res=300)
corrplot(corr = corrmatrix, type = "upper", tl.pos = "td", diag = FALSE, p.mat = p.mat, sig.level = 0.05, insig = "blank", col=coolwarm(20), tl.col="black")
dev.off()


png("Correlations_values.png", width = 2500, height = 2500, res=300)
corrplot(corrmatrix, method="color", col=coolwarm(20),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", #Text label color and rotation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE, number.cex=0.65
)
dev.off()




# MET
png("MET.png", width = 4000, height = 1500, res=300)
par(mfrow=c(1,3))
plot(data$MET, data$KOFFEIN, col=palette_people[data$person], pch=data$person+14, ylab="MET", xlab="Koffein")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$KOFFEIN ~data$MET)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$MET, data$Ansträngning, col=palette_people[data$person], pch=data$person+14, ylab="MET", xlab="Ansträngning")
legend("bottomleft", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Ansträngning ~data$MET)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topright", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$MET, data$Återhämtning, col=palette_people[data$person], pch=data$person+14, ylab="MET", xlab="Återhämtning")
legend("bottomleft", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Återhämtning ~data$MET)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topright", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

dev.off()



#SBT and DBT
png("SBT_DBT.png", width = 3500, height = 1700, res=300)
par(mfrow=c(1,2))
plot(data$SBT, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="SBT", xlab="Avkoppling")
legend("bottomleft", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~data$SBT)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topright", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$DBT, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="DBT", xlab="Avkoppling")
legend("bottomleft", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~data$DBT)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topright", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))


dev.off()





#Styrka+kardio
png("kardio_Styrka.png", width = 4500, height = 2500, res=300)
palette_styrka<-brewer.pal(4,"Accent")
par(mfrow=c(2,3))
boxplot(data$SDNN ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="SDNN", ylim=c(0,75), col=palette_styrka)
boxplot(data$RMSSD ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="RMSSD", ylim=c(0,75), col=palette_styrka)
boxplot(data$`Total Power` ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="Total Power", ylim=c(0,11), col=palette_styrka)
boxplot(data$`LF Power` ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="LF Power", ylim=c(0,10), col=palette_styrka)
boxplot(data$`HF Power` ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="HF Power", ylim=c(0,10), col=palette_styrka)
boxplot(data$`LF/HF ratio` ~ as.factor(data$`Kardio+Styrka`), xlab="Kardio+Styrka", ylab="LF/HF ratio", ylim=c(0,3), col=palette_styrka)
dev.off()


#ANOVA on Styrka+kardio
aov_data<-data
colnames(aov_data)[26]<-"KardStyr"

groups<-mat.or.vec(4,13)
rownames(groups)<-levels(aov_data$KardStyr)
colnames(groups)<-c("Kardio+Styrka group",
                    "SDNN", "SDNN_groups",
                    "RMSSD", "RMSSD_groups",
                    "Tot_pwr", "Tot_pwr_groups",
                    "LF_pwr", "LF_pwr_groups",
                    "HF_pwr", "SDNN_groups",
                    "LF_HF_ratio", "LF_HF_ratio_groups")
groups_significance<-mat.or.vec(2,6)
rownames(groups_significance)<-c("P value", "F value")
colnames(groups_significance)<-c("SDNN",
                    "RMSSD", 
                    "Tot_pwr", 
                    "LF_pwr", 
                    "HF_pwr", 
                    "LF_HF_ratio")



aov_data$KardStyr<-as.factor(aov_data$KardStyr)

amod1<-aov(SDNN ~ KardStyr, data=aov_data)
HSD1<-HSD.test(amod1, "KardStyr", group=T)
amod2<-aov(RMSSD ~ KardStyr, data=aov_data)
HSD2<-HSD.test(amod2, "KardStyr", group=T)
amod3<-aov(`Total Power` ~ KardStyr, data=aov_data)
HSD3<-HSD.test(amod3, "KardStyr", group=T)
amod4<-aov(`LF Power`  ~ KardStyr, data=aov_data)
HSD4<-HSD.test(amod4, "KardStyr", group=T)
amod5<-aov(`HF Power`  ~ KardStyr, data=aov_data)
HSD5<-HSD.test(amod5, "KardStyr", group=T)
amod6<-aov(`LF/HF ratio` ~ KardStyr, data=aov_data)
HSD6<-HSD.test(amod6, "KardStyr", group=T)

groups[,1]<-rownames(HSD1$groups)[order(rownames(HSD1$groups))]
groups[,2:3]<-as.matrix(HSD1$groups[order(rownames(HSD1$groups)),])
groups[,4:5]<-as.matrix(HSD2$groups[order(rownames(HSD2$groups)),])
groups[,6:7]<-as.matrix(HSD3$groups[order(rownames(HSD3$groups)),])
groups[,8:9]<-as.matrix(HSD4$groups[order(rownames(HSD4$groups)),])
groups[,10:11]<-as.matrix(HSD5$groups[order(rownames(HSD5$groups)),])
groups[,12:13]<-as.matrix(HSD6$groups[order(rownames(HSD6$groups)),])

groups_significance[1,1]<-summary(amod1)[[1]][1,5] #P-value
groups_significance[2,1]<-summary(amod1)[[1]][1,4] #F-value

groups_significance[1,2]<-summary(amod2)[[1]][1,5] #P-value
groups_significance[2,2]<-summary(amod2)[[1]][1,4] #F-value

groups_significance[1,3]<-summary(amod3)[[1]][1,5] #P-value
groups_significance[2,3]<-summary(amod3)[[1]][1,4] #F-value

groups_significance[1,4]<-summary(amod4)[[1]][1,5] #P-value
groups_significance[2,4]<-summary(amod4)[[1]][1,4] #F-value

groups_significance[1,5]<-summary(amod5)[[1]][1,5] #P-value
groups_significance[2,5]<-summary(amod5)[[1]][1,4] #F-value

groups_significance[1,6]<-summary(amod6)[[1]][1,5] #P-value
groups_significance[2,6]<-summary(amod6)[[1]][1,4] #F-value

write_ods(as.data.frame(groups), "Kardio_Styrka_ANOVA_groups.ods", sheet="ANOVA_groups")
write_ods(as.data.frame(groups_significance), "Kardio_Styrka_ANOVA_groups.ods", sheet="ANOVA_significance", append=T, row_names=T)





#Avkoppling
png("Avkoppling.png", width = 4000, height = 2500, res=300)
par(mfrow=c(2,3))
plot(data$SDNN, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="SDNN")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~data$SDNN)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$RMSSD, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="RMSSD")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~ data$RMSSD)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`Total Power`, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="Total Power")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~ data$`Total Power` )
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF Power`, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="LF Power")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~ data$`LF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`HF Power`, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="HF Power")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~ data$`HF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF/HF ratio`, data$`Avkoppling(Tid)`, col=palette_people[data$person], pch=data$person+14, ylab="Avkoppling", xlab="LF/HF Ratio")
legend("bottomright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$`Avkoppling(Tid)` ~ data$`LF/HF ratio`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

dev.off()


#Träning
png("Träning.png", width = 4000, height = 2500, res=300)
par(mfrow=c(2,3))
plot(data$SDNN, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="SDNN")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$SDNN)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$RMSSD, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="RMSSD")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$RMSSD)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`Total Power`, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="Total Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$`Total Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF Power`, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="LF Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$`LF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`HF Power`, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="HF Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$`HF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF/HF ratio`, data$Träning, col=palette_people[data$person], pch=data$person+14, ylab="Träning", xlab="LF/HF ratio")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Träning ~data$`LF/HF ratio`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

dev.off()


#Stressad
png("Stressad.png", width = 4000, height = 2500, res=300)
par(mfrow=c(2,3))
plot(data$SDNN, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="SDNN")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$SDNN)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$RMSSD, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="RMSSD")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$RMSSD)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`Total Power`, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="Total Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$`Total Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF Power`, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="LF Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$`LF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`HF Power`, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="HF Power")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$`HF Power`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))

plot(data$`LF/HF ratio`, data$Stressad, col=palette_people[data$person], pch=data$person+14, ylab="Stressad", xlab="LF/HF ratio")
legend("topright", c("Person 1", "Person 2", "Person 3"), col=palette_people, pch=seq(1:3)+14, bty="n")
regr<-lm(data$Stressad ~data$`LF/HF ratio`)
abline(regr, lty=2)
if(summary(regr)$coefficients[,4][2]<0.05){textcol="darkgreen"} else if (summary(regr)$coefficients[,4][2]<0.1){textcol="orange"} else {textcol="firebrick"}
legend("topleft", c(paste("R squared=", round(summary(regr)$r.squared,3)),
                    paste("P=", round(summary(regr)$coefficients[,4][2],3))), bty="n", text.col=c("black", textcol))


dev.off()


