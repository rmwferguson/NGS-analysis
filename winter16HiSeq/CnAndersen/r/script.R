myOtus <- read.delim("otu_table.txt", header = T, row.names = 1, skip = 1) #load in OTU table

head(myOtus)
library(vegan) 

taxa <- data.frame(otuID = rownames(myOtus), taxonomy = myOtus$taxonomy) # make an object for taxa information
myOtus <- myOtus[, -ncol(myOtus)] # remove taxa from OTU table
myOtus <- as.data.frame(t(myOtus)) #transpose OTU table

#normalise seq depth to min
otuNorm <- rrarefy(myOtus, min(rowSums(myOtus))) 

# OTU richness per sample
otuRichness <- specnumber(otuNorm) 

print(otuRichness)
write.table(otuRichness, "otuRichness.txt", sep="\t")

# shannon div

otuShannon <- diversity(otuNorm, "shannon") 
write.table(otuShannon, "otuShannon.txt", sep="\t")

# simpsons

otuSimpsons <- diversity(otuNorm, "simpson") 
write.table(otuSimpsons, "otuSimpsons.txt", sep="\t")

# even

OtuEven <- otuShannon/log(specnumber((otuNorm)))

#plot

otuRichPlot <- ggplot(data = envData, aes(Stage, otuRichness, fill = Stage)) +
  geom_boxplot()+
  labs(y = "OTU richness (S)", x = "SParticle size (um)")+
  #coord_cartesian(ylim = c(1,6))+
  theme(legend.position="none")+
  theme_bw()+
  #ggtitle("Gram negative") +
  theme(plot.title=element_text(face="italic"),legend.position="none", axis.title = element_text(size=7))


otuRichPlot


ShanPlot <- ggplot(data = envData, aes(Stage, otuShannon, fill = Stage)) +
  geom_boxplot()+
  labs(y = "Shannon Winer Index (H)", x = "Particle size (um)")+
  #coord_cartesian(ylim = c(1,6))+
  theme(legend.position="none")+
  theme_bw()+
  #ggtitle("Gram negative") +
  theme(plot.title=element_text(face="italic"),legend.position="none", axis.title = element_text(size=7))


ShanPlot


SimpPlot <- ggplot(data = envData, aes(Stage, otuSimpsons, fill = Stage)) +
  geom_boxplot()+
  labs(y = "Simpson's Index (1-D)", x = "Particle size (um)")+
  #coord_cartesian(ylim = c(1,6))+
  theme(legend.position="none")+
  theme_bw()+
  #ggtitle("Gram negative") +
  theme(plot.title=element_text(face="italic"),legend.position="none", axis.title = element_text(size=7))

theme(legend.position="none")

SimpPlot

EvenPlot  <- ggplot(data = envData, aes(Stage, OtuEven, fill = Stage)) +
  geom_boxplot()+
  labs(y = "Pielou's evenness (J)", x = "Particle size (um)")+
  theme(legend.position="none")+
  theme_bw()+
  #ggtitle("Gram negative") +
  theme(plot.title=element_text(face="italic"),legend.position="none", axis.title = element_text(size=7))

EvenPlot


multiplot(otuRichPlot, ShanPlot, SimpPlot, EvenPlot, cols=2)

tiff("div.tiff", height = 5, width = 7, units = "in", res=300)
#postscript("Figure1.tiff", width = 672, height = 672)
multiplot(otuRichPlot, ShanPlot, SimpPlot, EvenPlot, cols=2)
dev.off()


# alpha div test statistics

envData <- read.csv("env.csv", header = T) 

modRichness <- lm(otuRichness ~ envData$stage) 
modSimpson <- lm(otuSimpsons ~ envData$stage) 
modShannon <- lm(otuShannon ~ envData$stage) 


summary(modRichness)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)      1201.99      48.56   24.75  < 2e-16 ***
#  envData$stage  -153.91      12.47  -12.34 4.11e-14 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 127.8 on 34 degrees of freedom
#Multiple R-squared:  0.8176,	Adjusted R-squared:  0.8122 
#F-statistic: 152.4 on 1 and 34 DF,  p-value: 4.106e-14

summary(modSimpson)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    1.004277   0.005185 193.674  < 2e-16 ***
#  envData$stage -0.010282   0.001331  -7.722 5.57e-09 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.01364 on 34 degrees of freedom
#Multiple R-squared:  0.6369,	Adjusted R-squared:  0.6262 
#F-statistic: 59.63 on 1 and 34 DF,  p-value: 5.566e-09

summary(modShannon)

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)    6.18168    0.13493   45.81  < 2e-16 ***
#  envData$stage -0.42111    0.03465  -12.15 6.32e-14 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.355 on 34 degrees of freedom
#Multiple R-squared:  0.8129,	Adjusted R-squared:  0.8074 
#F-statistic: 147.7 on 1 and 34 DF,  p-value: 6.316e-14

# beta div

# make dist matrix

commSim <- vegdist(otuNorm, "jaccard") 

par(mfrow = c(1, 2)) # make heat map and cluster
heatmap(as.matrix(commSim)) 
plot(hclust(commSim)) 
par(mfrow = c(1, 1)) # make heat map and cluster

# permanover to test cluters

permResult <- adonis(commSim ~ envData$set + envData$stage, permutations = 1000) 
permResult

#Permutation: free
#Number of permutations: 1000

#Terms added sequentially (first to last)

#Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
#envData$set    2    0.9281 0.46404  2.6949 0.10714 0.002997 ** 
#  envData$stage  1    2.2244 2.22437 12.9179 0.25678 0.000999 ***
#  Residuals     32    5.5102 0.17219         0.63609             
#Total         35    8.6626                 1.00000             
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#NMds of data

nmdsResult <- metaMDS(otuNorm, "jaccard") 

plot(nmdsResult, type = "n") 
points(nmdsResult$points, display = "sites") 

nmdsResult


# make better nmdsPlot with sets as shapes and colour as stage


with(envData, levels(stage))
scl <- 6
colvec <- c("coral1", "darkgoldenrod4","chartreuse4","cyan4","cornflowerblue","hotpink")
treatvec <- c(15,16,17)

plot(nmdsResult, type = "n", scaling = scl)

with(envData, points(nmdsResult, display = "sites", col = colvec[stage], scaling = scl, pch = 19))

#cant get legend to work?
legend('topright', col=colvec[stage], legend=levels(envData$site), pch = 16, cex = 0.7)                  

tiff("nMDSCol.tiff", height = 7, width = 7, units = "in", res=600)

dev.off()

envData$set
