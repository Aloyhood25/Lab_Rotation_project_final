#age effect on sperm motility in Drosophila melanogaster
#Oliver Otti
#10.05.22

#BEFORE WE START SOME THINGS TO REMEMBER
# always comment code
# indent and spacing, use control+i to tell R to indent your code
# keep yourself to 80 characters per line, if possible: split comments and 
# function calls if necessary
# never save the environment!

# Loading packages and data -----------------------------------------------

#clear R's brain
rm(list=ls())


#load additional packages to have all the necessary functions for graphs and
# available
library(tidyverse)#this package contains all we need, is this not nice
library(wesanderson)

#load data into R
my_data<-read.csv("01 Raw data/sperm motility SF_sperm age interaction(1).csv",
                  header=T)

#verify if your data loaded and if all variables are correctly assigned
str(my_data)

# Sperm motility plots ---------------------------------------------------

#add "min" to the measurement timepoint to make it look nices in the figure
my_data$time.factor<-paste(my_data$time_recording,"min")

#order the levels
my_data$time.factor <- ordered(my_data$time.factor, levels = c("0 min", 
                                                               "2 min", "4 min","6 min"))

#Estimated marginal means plot ------------------------
#analysis of slope differences between age groups and treatments
library(lme4)
library(lmerTest)
library(emmeans)

# Assuming your model is something like: Motility ~ Time * Treatment +
#(1 | Individual)

m1<-lmer(mean_motility_3frames~time.factor+treatment_combination+
           (1|flyID)+time.factor:treatment_combination,
         data=my_data)

library(DHARMa)
sim_m1<-simulateResiduals(m1)
plot(sim_m1)

summary(m1)
library(car)
Anova(m1)
corrs<-vcov(m1)

write.csv(Anova(m1),file=
            "02 clean data/glmm sperm motility SF_sperm age interaction.csv")

emm <- emmeans(m1, specs = ~ time.factor*treatment_combination)
emm
emm_plot <- as.data.frame(emm)

#pairwise comparisons within each time of measurement
simp <- pairs(emm, simple = "treatment_combination")
simp

#was does this test do?
joint_tests(emm)
joint_tests(emm, by="time.factor")

levels(emm_plot$time.factor) <- c("0", "2", "4","6")

pdf("04 R plots/estimated marginal means sperm motility age and SF.pdf")
# Plot the EMMs
levels(my_data$time.factor) <- c("0", "2", "4","6")

ggplot()+
  geom_line(data=emm_plot, aes(x = time.factor, y = emmean, 
                               group = treatment_combination, 
                               color = treatment_combination),linewidth=1, 
            position=position_dodge(0.4))+
  geom_point(data=emm_plot, aes(x = time.factor, y = emmean,
                                group = treatment_combination, 
                                color =  treatment_combination),size=3, 
             position=position_dodge(0.4))+
  geom_errorbar(data=emm_plot, aes(x = time.factor, y = emmean, 
                                   group = treatment_combination, 
                                   color = treatment_combination, 
                                   ymin  =  emmean-SE,
                                   ymax  =   emmean+SE), width =  0.1, 
                linewidth  =  0.5, position=position_dodge(0.4))+
  geom_point(data=my_data, aes(x = time.factor, y = mean_motility_3frames, 
                               group = treatment_combination, 
                               color = treatment_combination),alpha=0.1,
             position=position_dodge(0.4))+
  labs(x = "Time in minutes", 
       y = expression(paste("Sperm motility (temporal noise ",sigma,")")),
       color="Treatment")+
  scale_color_manual(values=c("darkred","#E69F00","#0072B2", "#56B4E9"))+
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))+
  theme_classic()+
  theme(text=element_text(size=20),
        axis.text.x = element_text(size=16),legend.title=element_text(size=16),
        legend.text=element_text(size=13),legend.position="bottom")

dev.off()