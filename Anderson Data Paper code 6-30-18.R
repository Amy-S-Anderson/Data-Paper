### Amy Anderson Data Paper code, 6/30/18 ###
#This script is for the analyzing the association of reproductive state, hookworm infection, hemoglobin, and peripheral WBCs.



#clear workspace memory
rm(list=ls())

#set working directory
setwd("C:\\Users\\Amy\\Desktop\\Anemia.Project\\Tsimane.data")
#check that R is looking in the folder you told it to
getwd()

####Import giant data file ----
ad<-read.csv("C:\\Users\\Amy\\Desktop\\Anemia.Project\\Tsimane.data\\Biochem-Flow-Medical-Anthropometry-Reproduction-GPS-ADB.16.05.24-19Dec15PIDs.csv",header=TRUE,stringsAsFactors=FALSE)
str(ad)

#Collapse 'trimester' variable to create FemaleReproductiveState variable
ad$FRS <- ad$Trimester
preg <- c("1", "2", "3")
ad.orig <- ad
ad <- ad.orig
ad[ad$Trimester %in% preg, "FRS"] <- "Pregnant"
typeof(ad$Trimester)


#create "anemic" variable with demographic-specific WHO cutoffs
{ad$anemic<-NA
  ad<-ad[!is.na(ad$Age),]
  ad<-ad[!is.na(ad$male),]
  #non-pregnant women
  ad$anemic[ad$Age >= 15 & ad$male==0 & (ad$FRS!="Pregnant" | is.na(ad$FRS))] <- ((ad$Hb <= 11.9) + (ad$Hb <= 10.9) + (ad$Hb < 8))[ad$Age >= 15 & ad$male==0 & (ad$FRS!="Pregnant" | is.na(ad$FRS))]
  #pregnant women
  ad$anemic[which(ad$Age >= 15 & ad$male==0 & (ad$FRS=="Pregnant"))] <- ((ad$Hb <= 10.9) + (ad$Hb <= 9.9) + (ad$Hb < 7))[which(ad$Age >= 15 & ad$male==0 & ad$FRS=="Pregnant")]
  #men
  ad$anemic[ad$Age >= 15 & ad$male==1] <- ((ad$Hb <= 12.9) + (ad$Hb <= 10.9) + (ad$Hb < 8))[ad$Age >= 15 & ad$male==1]}

#create binary yes/no anemic variable
ad$anemic.binary <- ifelse(ad$anemic >0, 1,0)
#create binary yes/no tested for hookworm variable
ad$Tested.Hookworm <- ifelse(is.na(ad$Present.Hookworm) == TRUE,0,1)


#create a LymphCnt variable (lymphocytes per nanoliter of blood) to match EosinCnt and NeutroCnt variables
ad$LymphCnt <- (ad$WBC * ad$Lymphocitos)/100
ad$Neutrophils <- ad$NeutroCnt*100/ad$WBC



#####clean up data set ----
#delete obvious (physiologically impossible) errors in bloodwork
ad$WBC[ad$WBC>40]<-NA
ad$HTO[ad$HTO<20]<-NA
ad$Age[ad$Age<0]<-NA
#create MCHC variable from hemoglobin and hematocrit; create micro/normo/macrocytic categorical variable
ad$MCHC<-(ad$Hb*100)/ad$HTO
ad$MCHCCat<-0+(ad$MCHC<=32)+(ad$MCHC>=36)*2
#rename sex categorical variable from binary to male/female
ad$male[ad$male==1]<-"Male"
ad$male[ad$male==0]<-"Female"
#48.55% female; 51.45% male in a sample of 40,416 clinic visits
ad$comunid <- as.factor(ad$comunid)
#replace NA's under NumPartos for men with 'Male'
ad$FRS[ad$male=="Male"]<-"Male"
ad$Present.Hookworm <- ifelse(ad$Present.Hookworm>0, 1,0) #for some reason a few individuals have hookworm values of 0.5
library(dplyr)
#remove several hundred irrelevant medical diagnosis rows
ad <- ad %>%
 dplyr::select(-starts_with("ADBCat"), -starts_with("CCSL"), -starts_with("Med."), -starts_with("orina"), -starts_with("Med."), -starts_with("Treat."), -starts_with("Ninez"), -starts_with("Nivel"))





####SUBSETTING DATA ----




#Exclusion criteria

#data frame of all individuals 15-45 yrs with necessary variables
ad1 <- ad[-which(rowSums( is.na(ad[,(c("Hb", "pid", "Age", "male", "BMI.TZ", "Date", "FRS","Present.Hookworm", "comunid", "anemic.binary", "EosinCnt"))]) ) > 0), ] %>%
  filter(!duplicated(.)) %>%
  filter(Age>14.99 & Age <=45)




#How many women/cases in the THLHP database are between 15 and 45 yrs old?
all.repro.women <- ad %>%
  filter(male=="Female") %>%
  filter(Age>14.99 & Age <=45) %>%
  filter(!duplicated(.)) %>%
  filter(!is.na(FRS)) %>%
  #calculate bmi for observations that include weight and height but not bmi
  mutate(bmi = ifelse(!is.na(Weight) & !is.na(Height),(Weight/(Height/100)/(Height/100)), BMI))

#calculate effect of trimester on bmi, controlling for age and PID
library(lme4)
tc <- lmer(bmi ~ Trimester + Age + (1|pid), data = all.repro.women)
summary(tc)

#subtract trimester effect from each bmi in dataset. tc='trimester-corrected'. #use parameters for trimester, then calculate trimester-corrected bmi - individual bmi, minus trimester effect.
all.repro.women$tc.bmi <- ifelse(all.repro.women$Trimester=="Cycling", all.repro.women$bmi- 0.383172, ifelse(all.repro.women$Trimester=="2", all.repro.women$bmi- 0.564493, ifelse(all.repro.women$Trimester=="3", all.repro.women$bmi - 1.545059, ifelse(all.repro.women$Trimester=="Lactating", all.repro.women$bmi - 0.625295, ifelse(all.repro.women$Trimester=="1", all.repro.women$bmi, NA)))))
#use repeat measures on women to create average standardized bmi  for each woman who has at least 1 data point for bmi
all.repro.women <- all.repro.women %>%
  group_by(pid) %>% mutate(tcm.bmi = mean(tc.bmi, na.rm=T), tcm2.bmi = ifelse(!is.na(tc.bmi),tc.bmi,tcm.bmi)) ##tcm2.bmi has  bmi for non-missing observations and mean trimester-corrected bmi per person for missing cases of women who have at least one measure
#149 women are missing tc.bmi, but all 1016 have tcm.bmi2


#exclude women without complete cases for variables of interest
var.limited.women<- all.repro.women %>%
  filter(!is.na(pid)) %>%
  filter(!is.na(Date)) %>%
filter(!is.na(Hb)) %>%
  filter(!is.na(Present.Hookworm)) %>%
  filter(!is.na(EosinCnt)) %>%
  filter(!is.na(comunid)) %>%
  filter(FRS!= "Menopause") %>%
  filter(FRS!="Pre-Menst")
  
 
length(unique(var.limited.women$FRS=="Pre-Menst"))

#could not get filter command to return the correct number of cases for nulliparous women or women whose last birth was >10 yrs ago, so resorted to this bulky method instead.
var.limited.women <- var.limited.women[-which(var.limited.women$Age>25 & var.limited.women$NumPartos==0),]

#Refactor comunid
var.limited.women$comunid <-factor(var.limited.women$comunid)
nrow(var.limited.women) 


#rename to shorten variable name
women <- var.limited.women



#format data collection dates
women$Date2<-as.POSIXct(women$Date,format="%m/%d/%Y")
women$Year<-as.POSIXlt(women$Date2)$year+1900
#Refactor Year to get rid of years with zero observations (which were stopping models from converging)
women$Year <- as.factor(women$Year)
#create 'years since last birth' variable
women$DateLastBirth2 <- as.POSIXct(women$DateLastBirth, format="%m/%d/%Y")
women$yrs.since.last.birth <- as.numeric(difftime(women$Date2, women$DateLastBirth2, units = "weeks")/52)
women$days.since.last.birth <- as.numeric(difftime(women$Date2, women$DateLastBirth2, units = "days"))
women$DateLastBirth2 <- as.numeric(women$DateLastBirth2)
women <- women[-which(women$yrs.since.last.birth>10),] 
#format date of last menses as a date rather than a character
women$FUM <- as.POSIXct(women$FUM, format="%m/%d/%Y")
#Center and Scale (standardize) Age to make models run more smoothly
women$Age.TZ <- scale(women$Age)
#create 'time since last menses' variable
#ad$weeks.since.menses <- difftime(ad$Date2, ad$FUM, units = "weeks")
women$Trimester = factor(women$Trimester,levels =c("Cycling","1", "2","3","Lactating"))
#create days sequence across for graphing purposes
women$graphdays <- ifelse(women$FRS=="Cycling", -30, ifelse(women$FRS=="Pregnant", women$PregDays, women$days.since.last.birth))

#final exclusion criteria - trimester-corrected bmi
women <- women %>% 
  filter(!is.na(tcm2.bmi))


#log erythrocyte sed rate to normalize distribution
women$ESR <- log(women$VSG)
summary(is.na(women$ESR)) #5 women do not have ESR
women$ESR[women$ESR == 0] <- NA


#code for double checking # of women and # of cases against the results of the filter command  
  length(unique(women$pid)) #614 women; 1016 obs
length(unique(var.limited.women$pid[which(var.limited.women$FRS=="Pre-Menst")]))
length(unique(women$pid[which(is.na(women$tcm2.bmi))]))
table(is.na(var.limited.women$Hb))
table(is.na(all.repro.women$Hb))




#age-and-village-matched men
men <- ad1 %>% #already age-limited, and with Hb, BMI, hookworm, and EosinCnt data
  filter(male=="Male") 
men <- men[which(men$comunid %in% var.limited.women$comunid ),]
#refactor comunid
men$comunid <- factor(men$comunid)
#format data collection dates
men$Date2<-as.POSIXct(men$Date,format="%m/%d/%Y")
men$Year<-as.POSIXlt(men$Date2)$year+1900
#refactor year
men$Year <- factor(men$Year)
men$tcm2.bmi <- men$BMI
men$FUM <- as.POSIXct(men$FUM, format="%m/%d/%Y")
men$Age.TZ <- scale(men$Age)
#number of cases; number of men
length(unique(men$pid)) #778 observations; 489 men

#Remove 'Male Juvenile' category from male Trimester
table(men$Trimester)
men$Trimester <- "Male"

#reorder levels in the Trimester variable so the algorithm doesn't compare other reproductive states to first trimester as baseline
women$Trimester <- as.factor(women$Trimester)
women$Trimester = factor(women$Trimester,levels =c("Cycling","1", "2","3","Lactating"))


#combined sample of women and men
men.and.women <- bind_rows(women, men)
men.and.women$Trimester <- as.factor(men.and.women$Trimester)
men.and.women$Trimester = factor(men.and.women$Trimester,levels =c("Cycling","1", "2","3","Lactating", "Male"))
men.and.women$Present.Hookworm <- ifelse(as.numeric(men.and.women$Present.Hookworm)>0, 1, 0)
men.and.women$Present.Hookworm <- as.factor(men.and.women$Present.Hookworm)

#All sample subsets
str(women)
str(men)
str(men.and.women)


#dataframe with Present.hookworm NAs allowed for bias testing
ad2 <- ad[-which(rowSums( is.na(ad[,(c("Hb", "pid", "Age", "male", "BMI.TZ", "Date2", "Year", "FRS","Tested.Hookworm", "comunid", "anemic.binary"))]) ) > 0), ] %>%
  filter(!duplicated(.)) %>%
  filter(Age>14.99 & Age <=45) %>%
  filter(male=="Female")
 
nrow(ad2) # n = 4,220; 

hist(ad[which(ad$FRS=="Male"),]$BMI.TZ)
summary(ad[which(ad$FRS=="Male"),]$BMI.TZ)
nrow(ad2[which(ad2$FRS=="Male"),])






#DESCRIPTIVES ----
nrow(women) #1016 observations
length(unique(women$pid)) #614 women
hist(women$Age)
summary(women$Age) #median = 36.19, mean = 33.39
table(women$FRS)
#Cycling Lactating  Pregnant 
#560       142       163 


table(women[!duplicated(women$pid),]$FRS)
table(women$Trimester) #number of cases in each trimester
# 1   2   3
#70  78  69
women$comunid <- as.factor(women$comunid)

description.table <- men.and.women %>% 
  group_by(Trimester) %>%
  summarise(num_obs= n(), 
            num_women = length(unique(pid)), mean(Age), sd(Age), mean(Hb), sd(Hb), mean(Present.Hookworm), sd(Present.Hookworm), mean(BMI), sd(BMI), median(EosinCnt), quantile(EosinCnt,  probs=0.05, na.rm=T), quantile(EosinCnt,  probs=0.95, na.rm=T), median(NeutroCnt), quantile(NeutroCnt,  probs=0.05, na.rm=T), quantile(NeutroCnt,  probs=0.95, na.rm=T), median(LymphCnt), quantile(LymphCnt, probs=0.05, na.rm=T), quantile(LymphCnt, probs=0.95, na.rm=T), median(VSG, na.rm=T), quantile(VSG, probs=0.05, na.rm=T), quantile(VSG, probs=0.95, na.rm=T)) %>%
  t(.)
description.table
write.csv(data.frame(summary(description.table)$fixed[,c(1,3,4)]), file="description_table.csv")

preg.description.table <- women %>%
  filter(FRS=="Pregnant") %>%
  summarise(num_obs= n(), 
            num_women = length(unique(pid)), mean(Age), sd(Age), mean(Hb), sd(Hb), mean(Present.Hookworm), sd(Present.Hookworm), mean(BMI), sd(BMI), median(EosinCnt), quantile(EosinCnt,  probs=0.05, na.rm=T), quantile(EosinCnt,  probs=0.95, na.rm=T), median(NeutroCnt), quantile(NeutroCnt,  probs=0.05, na.rm=T), quantile(NeutroCnt,  probs=0.95, na.rm=T), median(LymphCnt), quantile(LymphCnt, probs=0.05, na.rm=T), quantile(LymphCnt, probs=0.95, na.rm=T), median(VSG, na.rm=T), quantile(VSG, probs=0.05, na.rm=T), quantile(VSG, probs=0.95, na.rm=T)) %>%
  t(.)
preg.description.table



cyc <- women %>%
  filter(FRS=="Cycling", IPI<=12) %>%
  summarise(count=n()) #106 cycling women gave birth in the last year
nrow(cyc)

cyc2 <-women %>%
  filter(FRS=="Cycling", IPI>12, IPI<=24) %>%
summarise(count=n()) #124 cycling women gave birth 1-2 years ago

lact <- women %>%
  filter(FRS=="Lactating", IPI<=12) %>%
  summarise(count=n())#95 lactating women gave birth in the last year

lact2 <- women %>%
  filter(FRS=="Lactating", IPI>12, IPI<=24) %>%
  summarise(count=n())#28 lactating women gave birth 1-2 years ago





#Look at crp as a potential indicator of inflammatory response
women.crp <- (women[which(women$crp!="NA"),]) %>% #n = 243
filter(ESR !="NA") #n = 240
women.crp$crp <- log(women.crp$crp)
plot(women.crp$ESR ~ women.crp$crp) #esr and crp are highly correlated, but we have ESR on far more women

#explore gastro symptoms as potential indicator of high worm load
library(MASS)
chisq.test(women$Present.Hookworm, women$gastro)#after this analysis I excluded medical diagnosis columns from the dataframe to make it a bit more manageable to work with, so this code no longer runs
#X-squared = 0.52759, df = 1, p-value = 0.4676
#no apparent relationship between hookworm infection and GI symptoms

#point biserial to test whether worst cases of anemia are correlated with GI symptoms
cor.test(women$Hb, women$gastro)
t = -0.44128, df = 839, p-value = 0.6591
alternative hypothesis: true correlation is not equal to 0
95 percent confidence interval:
  -0.08275015  0.05242375
sample estimates:
  cor 
-0.0152328 
#no relationship there either



  
#BIAS TESTS ----
#Do I need to redo these with brms? Does it matter?
library(lmerTest)
library(MCMCglmm)
#create binary tested for hookworm yes/no variable
#Do people tested for hookworm vary systematically in their Hb measurements from those not tested for hookworm?
 
bias.check.Hb <- MCMCglmm(Hb ~ Tested.Hookworm + Age + FRS, random=~comunid + Year + pid, data = ad2) 
summary(bias.check.Hb) #no


#Do hookworm-tested individuals have different WBC counts?
bias.check.WBC <- MCMCglmm(WBC ~ Tested.Hookworm + Age + FRS, random=~comunid + Year + pid, data = ad2)
summary(bias.check.WBC) #no

#Does Hookworm prevalence vary by year?
library(ggplot2)
library(tidyr)
hookworm.year <- aggregate(women$Present.Hookworm, by=list(women$Year), FUN=mean)
#rename columns so Year column can be used to merge with sample size dataframe
colnames(hookworm.year) <- c("Year", "HwPrevalence")
#Dataframe of # of ppl sampled each year
pid.by.year <- colSums(xtabs(~pid + Year, data = women))
#use tidyr to change rows to columns 
pid.by.year <- gather(pid.by.year, Year, pid)
#merge year and hookworm prevalence with year and sample size
hw.yr <- merge(hookworm.year, pid.by.year, by="Year")









###MODELS ----

library(brms)
#run models twice - does including averaged bmi data change the results?

#Is there a difference in effect of hookworm on Hb for men vs cycling women?
men.and.cycling <- men.and.women %>%
  filter(FRS!= "Pregnant") %>%
  filter(FRS!= "Lactating")
mod1.brms <- brm(Hb ~ Present.Hookworm + male + Age + tcm2.bmi + male:Present.Hookworm +  (1|comunid) + (1|Year) + (1|pid), data = men.and.cycling, chains = 2, control = list(max_treedepth = 12, adapt_delta=0.9), iter=4000)
summary(mod1.brms)
write.csv(data.frame(summary(mod1.brms)$fixed[,c(1,3,4)]), file="mod1.brms.csv")
save(hb1.brms, file = "mod1.brms.rda")


#MODEL 1
#For all adults,is hookworm infection a predictor of hemoglobin level?
hb1.brms <- brm(Hb ~ Present.Hookworm + male + Age + tcm2.bmi + male:Present.Hookworm +  (1|comunid) + (1|Year) + (1|pid), data = men.and.women, chains = 2, control = list(max_treedepth = 12, adapt_delta=0.9), iter=4000)
summary(hb1.brms)
write.csv(data.frame(summary(hb1.brms)$fixed[,c(1,3,4)]), file="hb1.brms.csv")
save(hb1.brms, file = "hb1.brms.rda")


#MODEL 2
#For women, does the effect of hookworm on Hb vary by reproductive state (cycling, pregnant, lactating)?
hb2.brms <- brm(Hb ~ Present.Hookworm + FRS + Age + tcm2.bmi + FRS:Present.Hookworm + (1|comunid) + (1|pid) + (1|Year), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)

#run again without added individuals from calculating mean bmi's
hb2.brms.check <- brm(Hb ~ Present.Hookworm + FRS + Age + tc.bmi + FRS:Present.Hookworm + (1|comunid) + (1|pid) + (1|Year), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(hb2.brms.check)

summary(hb2.brms)
write.csv(data.frame(summary(hb2.brms)$fixed[,c(1,3,4)]), file="hb2.brms.csv")
save(hb2.brms, file = "hb2.brms.rda")

#max tree depth prob doesn't need to be over 12. 


#MODEL 3
#For pregnant women, does the effect of hookworm on Hb vary by Trimester?
hb3.brms <- brm(Hb ~ Present.Hookworm + Trimester + Age + tcm2.bmi + Trimester:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women, chains=2, control= list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(hb3.brms)
write.csv(data.frame(summary(hb3.brms)$fixed[,c(1,3,4)]), file="hb3.brms.csv")
save(hb3.brms, file = "hb3.brms.rda")
#####THIS WILL SAVE YOU - marginal values function in the brms package
women$fit.hb3 <- fit.hb3[,1]
  
  fit.hb3 <- fitted(hb3.brms)         
df <- data.frame(Present.Hookworm=c(0,1,0,1,0,1,0,1,0,1), Trimester=c("Cycling","Cycling", "1","1","2","2","3","3","Lactating","Lactating"), Age=mean(women$Age), tcm2.bmi=mean(women$tcm2.bmi))
fitmeans.hb3 <- as.data.frame(fitted(hb3.brms, newdata = df, re_formula=NA))
hbplotdat <- as.data.frame(cbind(c(Trimester=c("Cycling","Cycling", "1","1","2","2","3","3","Lactating","Lactating")), Present.Hookworm=c(0,1,0,1,0,1,0,1,0,1), fitmeans.hb3))
colnames(hbplotdat) <- c("Trimester", "Hookworm", "Hb", "error", "Lower95", "Upper95")

#run again without added individuals from calculating mean bmi's
hb3.brms.check <- brm(Hb ~ Present.Hookworm + Trimester + Age + tc.bmi + Trimester:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women, chains=2, control= list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(hb3.brms.check)


#MODEL 4
#Does prevalence of hookworm infection vary by reproductive state?
 hw1.brms <- brm(Present.Hookworm ~ Trimester + Age + tcm2.bmi +  (1|comunid) + (1|Year) + (1|pid), family = binomial("logit"), data = women, chains=2, control= list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
 s.hw1 <- summary(hw1.brms)
 write.csv(data.frame(exp(summary(hw1.brms)$fixed[,c(1,3,4)])), file="hw1.brms.csv")
 save(hw1.brms, file = "hw1.brms.rda")
 
 #run again without added individuals from calculating mean bmi's
 hw1.brms.check <- brm(Present.Hookworm ~ Trimester + Age + tc.bmi +  (1|comunid) + (1|Year) + (1|pid), family = binomial("logit"), data = women, chains=2, control= list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
 summary(hw1.brms.check)#Doesn't seem to change model outcome. Still s
 
 #isolate cases that were added with tcm.bmi2 variable to look at prevalence of hookworm infection
 additions <- women %>%
   filter(is.na(tc.bmi))
 table(additions$Trimester, additions$Present.Hookworm)
 w.recbmi <- women %>%
   filter(!is.na(tc.bmi))
 
 exp(summary(hw1.brms)$fixed[,c(1,3,4)])
 hw1fe <- summary(hw1.brms)$fixed[,c(1,3,4)]
 
 wbmi <- mean(women$tcm2.bmi)#23.78779
 wage <- mean(women$Age)#33.39469
 
 #log-odds of hookworm infection for a cycling woman of average age and bmi
 cycling.odds <- exp(hw1fe[1,1] + wbmi*hw1fe[6,1] + wage*hw1fe[7,1])
 cycling.odds/(1+cycling.odds) # predicted probability that 28.220% of cycling women have hookworm. Should be closer to 50%
 
 
 table(women$Trimester, women$Present.Hookworm)
 prop.test(table(additions$Trimester, additions$Present.Hookworm))
 prop.test(table(w.recbmi$Trimester, w.recbmi$Present.Hookworm)) #Proportions of hookworm infection by trimester seem to folllow slightly different patterns in these two groups of women, but neither follow a pattern of increasing hookworm prevalence later in pregnancy.
 
 
#Odds ratio for pregnancy as a predictor of hookworm - exp(coef)
exp(0.32) #1.377 <-- 0.32 = beta for FRSPregnant, hw1.brms
#upper 95% CI
exp(0.32 + 1.96 * 0.19)
exp(0.7) 
#lower 95% CI
exp(0.32 - 1.96 *0.19)
exp(-0.06)
exp(cbind("Odds ratio" = coef(hwprev), confint.default(hwprev, level = 0.95))) ##try this


#MODEL 5
#Does the significance of hookworm infection as a predictor of anemia diagnosis vary by reproductive state?  
hw2.brms <- brm(anemic.binary ~ Present.Hookworm + FRS + Age + tcm2.bmi + FRS:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), family = binomial("logit"), data = women, chains = 2, control= list(max_treedepth = 10, adapt_delta=0.9),iter=4000)

summary(hw2.brms) 
write.csv(data.frame(summary(hw2.brms)$fixed[,c(1,3,4)]), file="hw2.brms.csv")
save(hw2.brms, file = "hw2.brms.rda")





exp(0.435)
exp(-0.807)
exp(2.005)


#MODEL 6
#Does Hookworm predict Eosinophil count for all groups?
#for both men and women
eosin1.brms <- brm(EosinCnt ~ Present.Hookworm + Age + tcm2.bmi + FRS + FRS:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = men.and.women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(eosin1.brms)
write.csv(data.frame(summary(eosin1.brms)$fixed[,c(1,3,4)]), file="eosin1.brms.csv")
save(eosin1.brms, file = "eosin1.brms.rda")


#MODEL 7
#and for just women, does the association between hookworm and eosinophils vary by reproductive state?
eosin2.brms <- brm(EosinCnt ~ Present.Hookworm + Age + tcm2.bmi + FRS + FRS:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(eosin2.brms)
write.csv(data.frame(summary(eosin2.brms)$fixed[,c(1,3,4)]), file="eosin2.brms.csv")
save(eosin2.brms, file = "eosin2.brms.rda")


#MODEL 8
#Are there differences by trimester?
eosin3.brms <- brm(EosinCnt ~ Present.Hookworm + Age + tcm2.bmi + Trimester + Trimester:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(eosin3.brms)
write.csv(data.frame(summary(eosin3.brms)$fixed[,c(1,3,4)]), file="eosin3.brms.csv")
save(eosin3.brms, file = "eosin3.brms.rda")

#run again without 'imputed' bmi to check results
eosin3.brms.check <- brm(EosinCnt ~ Present.Hookworm + Age + tc.bmi + Trimester + Trimester:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(eosin3.brms.check)


#checking association of hookworm with other WBC subtypes
#Lymphocytes
lympho.brms <- brm(LymphCnt ~ Present.Hookworm + Age + tcm2.bmi + FRS + (1|comunid) + (1|Year) + (1|pid), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(lympho.brms)
write.csv(data.frame(summary(lympho.brms)$fixed[,c(1,3,4)]), file="lympho.brms.csv")
save(lympho.brms, file = "lympho.brms.rda")
#Neutrophils
NeutroCnt.brms <- brm(NeutroCnt ~ Present.Hookworm + Age + tcm2.bmi  + FRS + (1|comunid) + (1|Year) + (1|pid), data = women, chains = 2, cores = 3, control = list(max_treedepth = 10, adapt_delta=0.9),iter=4000)
summary(lympho.brms)
write.csv(data.frame(summary(NeutroCnt.brms)$fixed[,c(1,3,4)]), file="neutro.brms.csv")
save(NeutroCnt.brms, file = "NeutroCnt.brms.rda")


#Does ESR respond to hookworm, and is this moderated by reproductive state/trimester?
#Look at relationship between hookworm, ESR, and FRS
#1016 women have ESR data; 7 women don't.
#MODEL 9
women$ESR<- log(women$VSG)
esr1.brms <- brm(ESR ~ FRS + Age + tcm2.bmi + Present.Hookworm + FRS:Present.Hookworm + (1|pid) + (1|comunid) + (1|Year), data = women)
summary(esr1.brms)
write.csv(data.frame(summary(esr1.brms)$fixed[,c(1,3,4)]), file="esr1.brms.csv")
save(esr1.brms, file = "esr1.brms.rda")


#MODEL 10
#look at changes in ESR response to hookworm across trimesters
esr2.brms <- brm(ESR ~ Trimester + Age + tcm2.bmi + Present.Hookworm + Trimester:Present.Hookworm + (1|pid) + (1|comunid) + (1|Year), data = women)
summary(esr2.brms)
write.csv(data.frame(summary(esr2.brms)$fixed[,c(1,3,4)]), file="esr2.brms.csv")
save(esr2.brms, file = "esr2.brms.rda")



########REPORTING MODEL RESULTS######
load("Hb1.brms.rda")
load("Hb2.brms.rda")
load("Hb3.brms.rda")
load("eosin1.brms.rda")
load("eosin2.brms.rda")
load("eosin3.brms.rda")
load("Hw1.brms.rda")
load("Hw2.brms.rda")
load("lympho.brms.rda")
load("NeutroCnt.brms.rda")
load("esr1.brms.rda")
load("esr2.brms.rda")
     

###calculate posterior samples
library(brms)
#hemoglobin - model 3
post<-posterior_samples(hb3.brms)
sum(post$`b_Present.Hookworm:TrimesterLactating`>0)/nrow(post) #0.8695 % of posterior sample >0
sum(post$`b_Present.Hookworm:Trimester1`<0)/nrow(post) #0.86675% of posterior sample <0
sum(post$`b_Present.Hookworm:Trimester2`>0)/nrow(post) 
sum(post$`b_Present.Hookworm:Trimester3`>0)/nrow(post) 

#eosinophils - model 8
post1<-posterior_samples(eosin3.brms)
sum(post1$`b_Present.Hookworm:TrimesterLactating`<0)/nrow(post1) #78.325% of posterior sample < 0
sum(post1$`b_Present.Hookworm:Trimester2`<0)/nrow(post1) #78.675
########check all interaction effects for model 8



#hookworm infection - model 4
#CI for cycling women OR of infection
post4 <- posterior_samples(hw1.brms)
sum(post4$`b_Trimester1`<0)/nrow(post4) # 47%
sum(post4$`b_Trimester2`<0)/nrow(post4) #14.4%
sum(post4$`b_Trimester3`<0)/nrow(post4) #29.98%
sum(post4$`b_TrimesterLactating`<0)/nrow(post4) #7.10%
#plot the predicted probability

fitted(hw1.brms, )
exp(-0.051)
exp(0.691)
#CI for lactation women OR of infection
exp(-0.109)
exp(0.688)


#MODEL 5 - LOGISTIC REGRESSION FOR ANEMIA 
#calculate predicted probability of anemia for a woman of average age and bmi if she is cycling, pregnant, or lactating
hw2fe <- summary(hw2.brms)$fixed[,c(1,3,4)]
cycling.nohw <- exp(hw2fe[1,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1])
cycling.nohw/(cycling.nohw+1) #8.95% chance of anemia

cycling.hw <- exp(hw2fe[1,1] + hw2fe[2,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1])
cycling.hw/(cycling.hw+1) #13.18% chance of anemia if cycling and hookworm

preg.nohw <- exp(hw2fe[1,1] + hw2fe[4,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1])
preg.nohw/(preg.nohw+1) #4.25% chance of anemia if pregnant, no hookworm

preg.hw <- exp(hw2fe[1,1] + hw2fe[2,1] + hw2fe[4,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1] + hw2fe[8,1])
preg.hw/(preg.hw+1) #8.35% chance of anemia if pregnant and hookworm

lact.nohw <- exp(hw2fe[1,1] + hw2fe[3,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1])
lact.nohw/(lact.nohw+1)#21.12% chance of anemia if amenorrhea

lact.hw <- exp(hw2fe[1,1] + hw2fe[2,1] + hw2fe[3,1] + mean(women$Age)*hw2fe[5,1] + mean(women$tcm2.bmi)*hw2fe[6,1] + hw2fe[7,1])
lact.hw/(lact.hw+1) #24.55% chance of anemia if amenorrhea and hookworm


#model 2 - interaction effect for pregnancy/hookworm
post3 <- posterior_samples(hb2.brms)
sum(post3$`b_Present.Hookworm:FRSPregnant`<0)/nrow(post3)#78.6% of posterior distribution below zero

#model 6 - eosinophils and hookworm, men and women
post5 <- posterior_samples(eosin1.brms)
sum(post5$`b_Present.Hookworm1:FRSMale`<0)/nrow(post5)#46.35% of posterior distribution below zero for males compared to cycling females

#model 11 - ESR, FRS
post6 <- posterior_samples(esr1.brms)
sum(post6$`b_Present.Hookworm`>0)/nrow(post6) #96.75% of posterior is above zero

#model 12 - ESR, Trimester
post7 <- posterior_samples(esr2.brms)
sum(post7$`b_Trimester1:Present.Hookworm`>0)/nrow(post7)#0.9435



hist(post$"b_Present.Hookworm:TrimesterLactating")
#Calculate marginal values
e3 <- summary(eosin3.brms)

marg.eosin.c <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) #1.768 = eosin for average cycling woman without hookworm
marg.eosin.chw <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) + e3$fixed["Present.Hookworm","Estimate"] #1.965 = eosin for average cycling woman with hookworm
marg.eosin.l <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) + e3$fixed["TrimesterLactating","Estimate"] #2.079 = eosin for average lactating woman without hookworm
marg.eosin.lhw <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) + e3$fixed["Present.Hookworm","Estimate"] + e3$fixed["TrimesterLactating","Estimate"] + e3$fixed["Present.Hookworm:TrimesterLactating","Estimate"]#2.142 = eosin for average lactating woman with hookworm
marg.eosin.2 <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) + e3$fixed["Trimester2","Estimate"] #1.390 = eosin for average 2nd trimester woman without hookworm
marg.eosin.2hw <- e3$fixed["Intercept","Estimate"] + e3$fixed["Age","Estimate"]*mean(women$Age) + e3$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi) + e3$fixed["Present.Hookworm","Estimate"] + e3$fixed["Trimester2","Estimate"] + e3$fixed["Present.Hookworm:Trimester2","Estimate"]#1.391 = eosin for average 2nd trimester woman with hookworm
lactdif <- marg.eosin.lhw - marg.eosin.l #0.063
cycdif <- marg.eosin.chw - marg.eosin.c #0.197
dif2 <- marg.eosin.2hw - marg.eosin.2 #0.0014

hb3.brms <- brm(Hb ~ Present.Hookworm + Trimester + Age + tcm2.bmi + Trimester:Present.Hookworm + (1|comunid) + (1|Year) + (1|pid), data = women
                

###calculate marginal Hemoglobin for men and women (using model hb1.brms) for graphing                
                #basic formula = observed values - observed.variable.values*b_fixed.variable + b_fixed.variable*mean(variable)
men.and.women$marg.Hb <- men.and.women$Hb  - men.and.women$Age*summary(hb1.brms)$fixed["Age","Estimate"] + summary(hb1.brms)$fixed["Age","Estimate"]*mean(men.and.women$Age) - men.and.women$tcm2.bmi*summary(hb1.brms)$fixed["tcm2.bmi","Estimate"] + summary(hb1.brms)$fixed["tcm2.bmi","Estimate"]*mean(men.and.women$tcm2.bmi)

grouped <- men.and.women %>%
    group_by(Present.Hookworm) %>%
  summarize(mean(marg.Hb))
  
#marginal ESR
women$marg.ESR <- women$ESR  - women$Age*summary(esr2.brms)$fixed["Age","Estimate"] + summary(esr2.brms)$fixed["Age","Estimate"]*mean(women$Age) - women$tcm2.bmi*summary(esr2.brms)$fixed["tcm2.bmi","Estimate"] + summary(esr2.brms)$fixed["tcm2.bmi","Estimate"]*mean(women$tcm2.bmi)

#marginal eosinophils
men.and.women$marg.Eosin <- men.and.women$EosinCnt  - men.and.women$Age*summary(eosin1.brms)$fixed["Age","Estimate"] + summary(eosin1.brms)$fixed["Age","Estimate"]*mean(men.and.women$Age) - men.and.women$tcm2.bmi*summary(eosin1.brms)$fixed["tcm2.bmi","Estimate"] + summary(eosin1.brms)$fixed["tcm2.bmi","Estimate"]*mean(men.and.women$tcm2.bmi)

mean(men.and.women$marg.Hb)  
  

#Calculating model-predicted summary values
#men without hookworm
#SAMPLE code from Carmen's model example
post<-posterior_samples(hb3.brms)

form.int<-function(x){
  paste0(format(median(x),digits=2,nsmall=2)," (",format(quantile(x,c(0.025)),digits=2,nsmall=2),"-",format(quantile(x,c(0.975)),digits=2,nsmall=2),")")
}

getcomparisons<-function(post){
  out<-data.frame(Tcyc=NA,Ncyc=NA,CycDiff=NA,Tpreg=NA,Npreg=NA,PregDiff=NA,TDpreg=NA,NDPreg=NA,DPregDiff=NA)
  
  #compare cycling
  estcycT<-post$b_Intercept
  estcycN<-post$b_Intercept+post$b_populationNHANES
  out$Tcyc<-formint(estcycT)
  out$Ncyc<-formint(estcycN)
  out$CycDiff<-formint(post$b_populationNHANES)
  
  #compare pregnant
  estpreT<-post$b_Intercept+post$b_repstatusPregnant
  estpreN<-post$b_Intercept+post$b_repstatusPregnant+post$b_populationNHANES+post$`b_repstatusPregnant:populationNHANES`
  out$Tpreg<-formint(estpreT)
  out$Npreg<-formint(estpreN)
  
  pregdiff<-estpreN-estpreT
  out$PregDiff<-formint(pregdiff)


library(lsmeans)
library(tidyr)
means.trim <- lsmeans(hb2, specs = "Trimester", data = men.and.women)
means.trim
means.hook <- lsmeans(hb2, specs = "Present.Hookworm", data = men.and.women)
means.hook
means.int <- lsmeans(hb2, specs = c("Trimester", "Present.Hookworm"), data = men.and.women)
means.int <- summary(means.int) 
  means.int <-as.data.frame(means.int[,1:4])


#*find example with linear model specifying contrasts in lsmeans*
#a priori contrasts: 

  





########    GRAPHS     #########

  library(ggplot2)
  library("wesanderson")
#Fig. 1 Association between hookworm and hemoglobin varies by Trimester/R.S.
#****how to do this with marginal means?
interaction.plot(x.factor     = means.int$Present.Hookworm,
                 trace.factor = means.int$Trimester, 
                 response     = means.int$lsmean, 
                 fun = mean,
                 type="b",
                 col=c("black","blue","darkgreen"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")


#change graphdays. if >500, = 500
#Fig. 2  Covariance among hemoglobin, eosinophils, and ESR across reproductive states
library(mgcv)
mHb<-gam(Hb~s(graphdays, by=as.factor(Present.Hookworm)) + Age + tcm2.bmi, data=women)


#predict: median aged woman.
median(women$Age)
mean(women$tcm2.bmi)
Days<-seq(-100,500,0.1)
pHbneg<-predict(mHb, newdata = data.frame(graphdays=Days,Age=36.19, tcm2.bmi=23.79, Present.Hookworm=0),se.fit=TRUE)
pHbpos<-predict(mHb, newdata = data.frame(graphdays=Days,Age=36.19, tcm2.bmi=23.79, Present.Hookworm=1),se.fit=TRUE)


plot.data<-data.frame(Days,pHbneg)
plot.data2<-data.frame(Days,pHbpos)

plot(plot.data$Days, plot.data$Estimate, col="blue")
lines(Days, plot.data$X97.5.ile, lty = 'dashed', col = 'blue')
lines(Days, plot.data$X2.5.ile, lty = 'dashed', col = 'blue')
par(new=TRUE)
plot(plot.data2$Days, plot.data2$Estimate, col="red")
lines(Days, pHbpos$fit+pHbpos$se.fit, lty = 'dashed', col = 'red')
lines(Days, pHbpos$fit-pHbpos$se.fit, lty = 'dashed', col = 'red')

#Alternative fig. 1 - boxplot of hb levels by hw status, across trimesters
men.and.women$Trimester = factor(men.and.women$Trimester,levels =c( "Male","Cycling","1", "2","3","Lactating")) #reorder levels of Trimester factor for plotting purposes

tiff(file="Fig.1 Hb-Hookworm boxplot by trimester.tiff", compression = "lzw")
ggplot(data = men.and.women,  
       aes(x = Trimester, y = marg.Hb, fill = Present.Hookworm)) +
  ylim(c(8, NA)) + 
  geom_jitter(aes(color = Present.Hookworm, alpha = 0.5)) +
  scale_color_manual(values = c("steelblue1", "red2")) + #jitter color palette
  geom_boxplot(outlier.shape=NA, alpha=0.7) + #outlier.shape=NA gets rid of outliers on boxplot, b/c they are doubled in the jitter
    scale_fill_manual(values = c("steelblue1", "red2")) #boxplot color palette
dev.off()


#Fig. 2? Eosinophils by trimester
tiff(file="Fig.2 Eosin-Hookworm boxplot by trimester.tiff", compression = "lzw")
ggplot(data = men.and.women,  
       aes(x = Trimester, y = marg.Eosin, fill = Present.Hookworm)) +
  ylim(c(0, 6)) + 
  geom_jitter(aes(color = Present.Hookworm, alpha = 0.5)) +
  scale_color_manual(values = c("steelblue1", "red2")) + #jitter color palette
  geom_boxplot(outlier.shape=NA, alpha=0.7) + #outlier.shape=NA gets rid of outliers on boxplot, b/c they are doubled in the jitter
  scale_fill_manual(values = c("steelblue1", "red2")) #boxplot color palette
dev.off()


#Fig. 3? ESR boxplot by reproductive state and hookworm

women$expmarg.ESR <- exp(women$marg.ESR)

tiff(file="Fig.3 ESR-Hookworm boxplot trimester.tiff", compression = "lzw")
ggplot(data = women,  
       aes(x = Trimester, y = expmarg.ESR, fill = factor(Present.Hookworm))) +
    geom_jitter(aes(color = factor(Present.Hookworm), alpha = 0.5), position=position_jitterdodge()) + 
  scale_color_manual(values = c("steelblue1", "red2")) + #jitter color palette
  geom_boxplot(outlier.shape=NA, alpha=0.5, fill = factor(Present.Hookworm)) + #outlier.shape=NA gets rid of outliers on boxplot, b/c they are doubled in the jitter
  scale_fill_manual(values = c("steelblue1", "red2"), 
                    name = "Hookworm Infection",
                    labels = c("Absent", "Present")) #boxplot color palette

dev.off()

#attempt at jitter of fitted values with forest plot of parameter estimates and CIs for Hb (model 3)
ggplot(data = women,  
       aes(x = Trimester, y = fit.hb3, fill = factor(Present.Hookworm))) +
  geom_jitter(aes(color = factor(Present.Hookworm), alpha = 0.5), position=position_jitterdodge()) +
  scale_color_manual(values = c("steelblue1", "red2")) +
  stat_summary(fun.data=mean, fun.args = list(mult=1), 
               geom="pointrange", color="red")
  
  
  
  geom_pointrange(data=hbplotdat, mapping=aes(x=Trimester, y=Hb, ymin=Lower95, ymax=Upper95, group=factor(Hookworm)), position=position_dodge(width=0.5),size=1, color="blue", fill="white") +
  geom_col(position=dodge) +
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95), position="dodge")
fit.hb3 <- data.frame(fit.hb3)
ggplot(data = women,
       aes(y = fit.hb3, group=Trimester)) +
  geom_jitter()
  
  
  geom_pointrange()
  

ylim(c(0, 6)) + 
  position=position_dodge(width=0.5)
ymin=2.5%ile, ymax=97.5%ile



#calculate points shown in plot above (mean Hb in each group)
FRS.hw.hb <- women %>%
  group_by(FRS, Present.Hookworm) %>%
  summarize(hb.m = mean(Hb)) #hookworm prevalence by FRS

##******How to plot the differences between hook-pos and hook-neg means for each level in Trimester?

#try it in ggplot - a line graph
#limit dataframe to relevant variables and calculate sem
Hb2 <- men.and.women[,c("marg.Hb", "Present.Hookworm", "Trimester", "pid")] %>%
    group_by(Trimester, Present.Hookworm) %>%
  summarise(hb_groups = mean(marg.Hb),
            hb_sem = (sd(marg.Hb)/sqrt(length(marg.Hb)))) #This is raw data. Use Marginal means. 
Hb2 #Tibble with Hb means and sem by trimester and hookworm

Hb2$Present.Hookworm <- as.factor(Hb2$Present.Hookworm)

tiff(file="Fig.1 Hb-Hookworm interaction x trimester.tiff", compression = "lzw")
ggplot(data = Hb2,
  aes(x = Present.Hookworm, y = hb_groups, color = Trimester)) +
  geom_line(aes(group = Trimester),size=1.25) +
  geom_point() +
  geom_errorbar(aes(x = Present.Hookworm, ymin = hb_groups - hb_sem, ymax = hb_groups + hb_sem), width = .03) +
  labs(title = "Hookworm-Hemoglobin association by Reproductive State",
       subtitle = "Error bars indicate standard error of the mean",
       x = "Hookworm status",
       y = "Hemoglobin (g/dL)") +
    theme(axis.line = element_line(colour = "black"), text=element_text(size=16), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), title=element_text(size=12)) 
dev.off()



#eosin interaction plot
ggplot(data = men.and.women,  
       aes(x = MalePregNot, y = EosinCnt, color = Present.Hookworm)) +
  geom_boxplot() +
  geom_point() +
  geom_errorbar(aes(x = Present.Hookworm, ymin = hb_groups - hb_sem, ymax = hb_groups + hb_sem), width = .03)

factor(plot.data$Trimester, levels=c("Cycling","1", "2", "3", "Lactating")) %>%

#try it again in ggplot, but this time as a bar graph
  ggplot(data=hb2) +
  aes(x=Trimester, y=hb_groups, fill=as.factor(Present.Hookworm)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=hb_groups - hb_sem, ymax=hb_groups + hb_sem),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9)) +
  coord_cartesian(ylim=c(11,13.25)) +
  labs(x = "Trimester/Reproductive State",
       y = "Mean Hemoglobin (g/dL)") +
  scale_fill_discrete(name="Hookworm", labels=c("Absent","Present")) +
  theme_bw()

  
#Alternative Fig. 2 - boxplot of marginal effects, hw/no hw across Trimesters








  #blank white background, no grid lines
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 




  

  
  
  
  






#plot interaction effect: FRSXHW = Hb
interaction.plot(x.factor     = women[which(women$EosinCnt!="NA"),]$Present.Hookworm,
                 trace.factor = women[which(women$EosinCnt!="NA"),]$FRS, 
                 response     = women[which(women$EosinCnt!="NA"),]$EosinCnt, 
                 fun = mean,
                 type="b",
                 col=c("black","blue","darkgreen"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")

#calculate points shown in plot above (mean eosincnt in each group)
FRS.hw.eosin <- women %>%
  filter(EosinCnt!="NA") %>%
  group_by(FRS, Present.Hookworm) %>%
  summarize(eosin.m = mean(EosinCnt)) #hookworm prevalence by FRS


#Forest plot
library(dotwhisker)
library(broom)
library(dplyr)





interaction.plot(x.factor     = women[which(women$EosinCnt!="NA"),]$Present.Hookworm,
                 trace.factor = women[which(women$EosinCnt!="NA"),]$Trimester, 
                 response     = women[which(women$EosinCnt!="NA"),]$EosinCnt, 
                 fun = mean,
                 type="b",
                 col=c("black","blue","darkgreen"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")




interaction.plot(x.factor     = women[which(women$ESR!="NA"),]$Present.Hookworm,
                 trace.factor = women[which(women$ESR!="NA"),]$Trimester, 
                 response     = women[which(women$ESR!="NA"),]$ESR, 
                 fun = mean,
                 type="b",
                 col=c("black","blue","darkgreen"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")
  
  #51 women in 1st trimester; 57 women in2nd trimester; 51 women in 3rd trimester




####PLOTS---######
library(tidyr)
library(dplyr)
library(ggplot2)

#Hb loss by trimester, bar graph with error bars
hb.trim.barplot <- ggplot(women, aes(x=Trimester, y=Hb, fill=as.factor(Present.Hookworm))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=EosinCnt-se, ymax=EosinCnt+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

ggplot(data = women, aes(x = Trimester, y = Hb, fill = Present.Hookworm)) + 
         geom_col(position = "dodge")
       
#That should give you adjacent bars for each trimester, allowing us to easily see the difference between hookworm statuses.
       
#Being picky: I would relabel your hookworm status as present / absent. If you follow the above code then it'll be:
         + 
         scale_fill_manual(labels = c("Absent", "Present"))
       
#look at ESR by trimester
esr.trim <- women %>%
  group_by(Present.Hookworm, Trimester) %>%
  filter(ESR!="NA") %>%
  filter(Trimester=="1"|"2"|"3")
plot(esr.trim$ESR ~ esr.trim$Trimester)
esr.boxplot <- ggplot(esr.trim, aes(x=Trimester, y=ESR, fill=as.factor(Present.Hookworm))) +
  geom_boxplot()

eosin.trim.boxplot <- ggplot(women, aes(x=Trimester, y=EosinCnt, fill=as.factor(Present.Hookworm))) +
  geom_boxplot()

hb.trim.boxplot <- ggplot(women, aes(x=Trimester, y=Hb, fill=as.factor(Present.Hookworm))) +
  geom_boxplot()

#Look at CRP by trimester
crp.trim.boxplot <- ggplot(women.crp, aes(x=Trimester, y=crp, fill=as.factor(Present.Hookworm))) +
  geom_boxplot()

#bar plot of hookworm in pregnant women
barplot(women[which(women$FRS=="Pregnant"),]$Present.Hookworm ~ women[which(women$FRS=="Pregnant"),]$Trimester)
#bar plot with error bars
eosin.trim.barplot <- ggplot(women, aes(x=Trimester, y=EosinCnt, fill=as.factor(Present.Hookworm))) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=EosinCnt-se, ymax=EosinCnt+se),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

summarySE(women, measurevar="EosinCnt", groupvars=c("Present.Hookworm","Trimester"))


