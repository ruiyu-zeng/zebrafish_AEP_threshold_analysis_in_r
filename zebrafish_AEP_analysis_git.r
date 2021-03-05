# load the necessary libraries for visualization and computation 
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(ggpubr)

zebrafish_r <- read_xlsx("zebrafish_raw_11102020.xlsx")
str(zebrafish_r)
typeof(zebrafish_r)
zebrafish_r$freq <- as.factor(zebrafish_r$freq)
zebrafish_r$age_group <- as.factor(zebrafish_r$age_group)
# generate contingency table
fish_contigency <- table(zebrafish_r$age_group,zebrafish_r$frequency)

#visualize data
# update on 2/13/2020, make the tuning curve with x on a log scale 
ggplot(data = zebrafish_r, aes(x=frequency, y = threshold, color = age_group))+geom_line()
# make boxplot
ggboxplot(zebrafish_r,x="frequency",y="threshold",color="age_group")
ggboxplot(zebrafish_r,x="age_group",y="threshold")

# make lineplot (code below produces plot I need to print)
tuning_curve <- ggline(zebrafish_r,numeric.x.axis = TRUE, x="freq",y="threshold",color="age",palette= c("black","chocolate3","slateblue"),size=1,add=c("mean_se"),position="dodge")
tuning_curve <- tuning_curve + geom_jitter()
ggpar(tuning_curve,xscale = "log10",ticks = TRUE,tickslab = TRUE, xticks.by = 100)
postscript("AEP_tuning_curve.eps",width=20,height =17,horizontal = FALSE, onefile = FALSE)
pdf("AEP_tuning.pdf",width = 20, height = 17)
tuning_curve
dev.off()



# Fit a linear mixed effects model using lme4
library(lme4)
# remove the attributes in age group -- which are previously established contrasts 
attr(zebrafish_r$age_group,"contrasts") <- NULL

lme_mod1 <- lmer(threshold ~  age + freq + (1+freq|ID), data= zebrafish_r,REML= FALSE)
lme_mod2 <- lmer(threshold ~  freq + (1+freq|ID),data=zebrafish_r, REML = FALSE)
anova(lme_mod1,lme_mod2)
lme_mod3 <- lmer(threshold ~ age_group + (1|ID), data= zebrafish_r,REML= FALSE)
lme_mod4 <- lmer(threshold ~ frequency + age_group + (1|ID), data= zebrafish_r,REML= FALSE)
anova(lme_mod3,lme_mod4)

#Use the nparLD packages to do non-parametric testing
library(nparLD)
boxplot(threshold ~ age_group * frequency,data=zebrafish_r, names = FALSE, col = c("grey",2,3), lwd=2)
axis(1,at =2, labels = "115Hz",font =1, cex =2)
axis(1, at =5, labels = "800Hz",font =1, cex=2)
axis(1, at =8, labels = "1850Hz",font =1, cex=2)
axis(1, at =11, labels = "3700Hz",font =1, cex=2)
axis(1, at =14, labels = "4500Hz",font =1, cex=2)

ex.f1f1_SPL <- f1.ld.f1(zebrafish_r$threshold,time=zebrafish_r$freq,group = zebrafish_r$age,subject = zebrafish_r$ID, time.name = "frequency",group.name = "age_group",description = TRUE)

