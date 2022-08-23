#--------------------------------------#
# GDS Project
# 27/07/2022
#--------------------------------------#


#- Get the packages ----------------------------------------------------------#
source("./Scripts/packages.R") # This script will install/load in packages that you need for this project
library(ggplot2)
library(readxl)
#- Get the data --------------------------------------------------------------#
### Load in the data from the file in the /Input/ folder
dat_raw <- read_sav("./Input/GDS_V1_Dataset_20220708_(N=6,800 with 799 cases from 23 studies).sav")
population <- read_excel('./Input/GDS_SummaryRows_WithVerification_2022-08-08.xlsx')
population <- subset(population,population$DIAGNOSTIC_INTERVIEW =='SCID')
population <- population[,c(4,7)]
population <- population[order(population$`AUTHOR, YEAR`),]
#- Calculate the GDS scores over the rows ------------------------------------#
### This will calculate the GDS-15 scores for each patient
##### There are also two columns (GDS30_TOTAL, GDS15_TOTAL) that give the GDS totals BUT
#####   only for studies where they reported ONLY the total and not the individual item scores
dat <- dat_raw %>%
  mutate(
    # Get the sum for the individual questions
    GDS15_sum = GDS_Q1+GDS_Q2+GDS_Q3+GDS_Q4+GDS_Q7+GDS_Q8+GDS_Q9+GDS_Q10+GDS_Q12+GDS_Q14+GDS_Q15+GDS_Q17+GDS_Q21+GDS_Q22+GDS_Q23,
    # If the sums are missing, then use the GDS15_TOTAL variable 
    GDS15 = if_else(is.na(GDS15_sum), GDS15_TOTAL, GDS15_sum))

#- Analysis ------------------------------------------------------------------#
### Lu, here you can try now to calculate the values which we talked about last week. 
### They are also listed in the protocol in the data analysis section. 
### Use the GDS15 variable which I created above for all analyses.

#Select the studies that used SCID as diagnostic criterion
dat <- subset(dat,dat$DEP_CRITERION == 1)
#Remove the participants that did not have GDS-15 scores.
dat <- dat %>% drop_na(c('GDS15'))
write.csv(dat,'dat.csv',row.names = FALSE)


#Select the columns that are needed for main analysis
table<- data.frame(dat$STUDY_AUTHOR_YEAR,dat$COUNTRY,dat$GDS15,dat$MDD)
#Change the column names
colnames(table) <- c('Author','Country','GDS15','MDD')

# Replace the number with the name of the country that it represents
table$Country <- gsub('16','Nigeria',table$Country)
table$Country <- gsub('15','Turkey',table$Country)
table$Country <- gsub('14','Korea',table$Country)
table$Country <- gsub('13','Greece',table$Country)
table$Country <- gsub('12','Iran',table$Country)
table$Country <- gsub('11','Spain',table$Country)
table$Country <- gsub('10','Italy',table$Country)
table$Country <- gsub('9','India',table$Country)
table$Country <- gsub('8','Brazil',table$Country)
table$Country <- gsub('7','Australia',table$Country)
table$Country <- gsub('6','China',table$Country)
table$Country <- gsub('5','Netherlands',table$Country)
table$Country <- gsub('4','Germany',table$Country)
table$Country <- gsub('3','UK',table$Country)
table$Country <- gsub('2','USA',table$Country)
table$Country <- gsub('1','Canada',table$Country)


#Calculate the total number of participants in each study.
total <- table %>% group_by(Author)%>% tally()

#Select the columns for age and sex analysis in table 1.
table_age_sex<- data.frame(dat$STUDY_AUTHOR_YEAR,dat$AGE_CONT,dat$SEX,dat$GDS15)
colnames(table_age_sex) <- c('Author','Age','Sex','GDS15')
table_age_sex<- table_age_sex %>% drop_na()

#Calculate the mean age of the participants in each study.
mean_age <- table_age_sex %>% group_by(Author)%>% summarise(Mean=mean(Age))

#Calculate the percentage of females in each study.
female_percentage <- table_age_sex %>% group_by(Author)%>% summarise(female_percentage=sum(Sex<1),total=sum(Sex<2))
female_percentage$female_percentage<- female_percentage$female_percentage/female_percentage$total
female_percentage$female_percentage <- female_percentage$female_percentage*100


#Calculate the percentage of participants who scored above cutoff 5.
cutoff_5<- table %>% group_by(Author)%>% summarise(cutoff_5=sum(GDS15>=5),cutoff_5_percentage=sum(GDS15>=5))
cutoff_5[, 3] <- data.frame(lapply(cutoff_5[, 3], function(x) x/total[, 2]))
cutoff_5[, 3]<- cutoff_5[, 3]*100

#Calculate the percentage of participants that were diagnosed with MDD.
MDD_percentage <- table %>% group_by(Author,Country)%>% summarise(diagnosed_n=sum(MDD),diagnosed_percentage=sum(MDD))
MDD_percentage[, 4]  <- data.frame(lapply(MDD_percentage[, 4], function(x) x/total[, 2]))
colnames(MDD_percentage)[1] <- 'Author, year'
MDD_percentage$diagnosed_percentage <- MDD_percentage$diagnosed_percentage*100

#Calculate the difference of these percentages.
difference_5 <- data.frame(cutoff_5$cutoff_5_percentage - MDD_percentage$diagnosed_percentage)
difference_5 <- cbind(MDD_percentage[,1],difference_5)
colnames(difference_5) <- c('Author, year','difference_5')

#Calculate the ratio of these percentages.
ratio_5 <- data.frame(cutoff_5$cutoff_5_percentage / MDD_percentage$diagnosed_percentage)
ratio_5 <- cbind(MDD_percentage[,1],ratio_5)
colnames(ratio_5) <- c('Author, year','ratio_5')

#Calculate the pooled prevalence at each cutoff.
pool_cutoff<- table %>% summarise(cutoff_0=sum(GDS15>=0),cutoff_1=sum(GDS15>=1),cutoff_2=sum(GDS15>=2),cutoff_3=sum(GDS15>=3),cutoff_4=sum(GDS15>=4),cutoff_5=sum(GDS15>=5),cutoff_6=sum(GDS15>=6),
                                  cutoff_7=sum(GDS15>=7),cutoff_8=sum(GDS15>=8),cutoff_9=sum(GDS15>=9),cutoff_10=sum(GDS15>=10),cutoff_11=sum(GDS15>=11),cutoff_12=sum(GDS15>=12),cutoff_13=sum(GDS15>=13),
                                  cutoff_14=sum(GDS15>=14),cutoff_15=sum(GDS15>=15))
pool_cutoff <- pool_cutoff/sum(total$n)
#Calculate the pooled prevalence of MDD.
pool_MDD <- data.frame(sum(table$MDD)/sum(total$n))
colnames(pool_MDD)[1] <- 'prevalence of MDD'

#calculate the pooled difference. We can see that 'prevalence-match cutoff' is 8.
pool_difference <- data.frame(pool_cutoff-pool_MDD$`prevalence of MDD`)

#Calculate the percentage of participants who scored above cutoff 8.
cutoff_8<- table %>% group_by(Author)%>% summarise(cutoff_8=sum(GDS15>=8),cutoff_8_percentage=sum(GDS15>=8))
cutoff_8[, 3] <- data.frame(lapply(cutoff_8[, 3], function(x) x/total[, 2]))
cutoff_8[, 3]<- cutoff_8[, 3]*100

#Calculate the difference (cutoff 8)
difference_8 <- data.frame(cutoff_8$cutoff_8_percentage - MDD_percentage$diagnosed_percentage)
difference_8 <- cbind(MDD_percentage[,1],difference_8)
colnames(difference_8) <- c('Author, year','difference_8')

#Calculate the ratio (cutoff 8)
ratio_8 <- data.frame(cutoff_8$cutoff_8_percentage / MDD_percentage$diagnosed_percentage)
ratio_8 <- cbind(MDD_percentage[,1],ratio_8)
colnames(ratio_8) <- c('Author, year','ratio_8')



#Create a data.frame to hold the study-level data
table1<-as.data.frame(matrix(NA, nrow=14,ncol=13))
colnames(table1)<-c('Author, year','Country', 'Population','N total','N (%) major depression','Mean age','% female','N (%) GDS-15 ≥ 5',
'% difference: GDS-15 ≥ 5 – major depression','Ratio: GDS-15 ≥ 5 / major depression','N (%) GDS-15 ≥ 8','% difference: GDS-15 ≥ 8 – major depression',
'Ratio: GDS-15 ≥ 8 / major depression')

#Fill in the data.frame
table1$`Author, year` <- cutoff_percentage$`Author, year`
table1$Country <- MDD_percentage$Country
table1$Population <- population$SAMPLE_DESCRIPTION
table1$`N total` <-paste(total$n)
table1$`N (%) major depression` <- paste(MDD_percentage$diagnosed_n,' (',format(round(MDD_percentage$diagnosed_percentage,digits = 1),nsmall=1),'%)',sep = '')
table1$`Mean age` <- as.numeric(paste(format(round(mean_age$Mean,digits=1),nsmall=1)))
table1$`% female` <- paste(round(female_percentage$female_percentage,digits=1),'%',sep = '')
table1$`N (%) GDS-15 ≥ 5` <- paste(cutoff_5$cutoff_5,' (',format(round(cutoff_5$cutoff_5_percentage,digits = 1),nsmall=1),'%)',sep = '')
table1$`% difference: GDS-15 ≥ 5 – major depression` <- paste(format(round(difference_5$difference_5,digits = 1),nsmall=1),'%',sep = '')
table1$`Ratio: GDS-15 ≥ 5 / major depression` <- paste(format(round(ratio_5$ratio_5,digits = 1),nsmall=1))
table1$`N (%) GDS-15 ≥ 8` <- paste(cutoff_8$cutoff_8,' (',format(round(cutoff_8$cutoff_8_percentage,digits = 1),nsmall=1),'%)',sep = '')
table1$`% difference: GDS-15 ≥ 8 – major depression` <- paste(format(round(difference_8$difference_8,digits = 1),nsmall=1),'%',sep = '')
table1$`Ratio: GDS-15 ≥ 8 / major depression` <- paste(format(round(ratio_8$ratio_8,digits = 1),nsmall=1))
table1 <- format(table1,nsmall=1)

#Export the table 1 as a csv file
write.csv(table1,'table1.csv',row.names = FALSE)

#Graph 1: study-level differences by sample size
plot1 <- data.frame(cbind(difference_8$`Author, year`,difference_8$difference_8,total$n))
colnames(plot1) <- c('Author,year','difference','sample_size')
plot1$difference <- as.numeric(plot1$difference)
plot1$sample_size <- as.numeric(plot1$sample_size)
plot1$difference<- round(plot1$difference,digits = 1)

ggplot(plot1, aes(x=difference, y=sample_size)) +
  geom_point(size=1.5)+theme_bw()+labs(x='Difference GDS-15 cutoff 8 or more - SCID',y='Sample size')+
  scale_x_continuous(breaks = seq(-20,20,5))+scale_y_continuous(breaks = seq(0, 2100, 200))
ggsave("study-level differences by sample size.png", units="in", width=5, height=4, dpi=300)

#Graph2: Prevalence estimates and 95% CI based on each GDS-15 cutoff threshold from ≥ 0 to ≥ 15.
plot2 <- data.frame(t(pool_cutoff))
colnames(plot2)[1] <- 'Percentage'
plot2$Percentage <- round(plot2$Percentage*100,digits=1)
plot2$cutoff <- 0:15

ggplot(plot2, aes(x=cutoff, y=Percentage)) +
  geom_bar(stat='identity',fill='grey')+theme_bw()+labs(x='Cutoff',y='Proportion (%) at or above the cutoff')+scale_x_continuous(breaks = seq(0,15,1))+geom_hline(yintercept=14.2)+ geom_text(aes(0,14.2,label = 'SCID prevalence = 14.2%', vjust = -1, hjust = -0.003))
ggsave("Pooled prevalence estimates above each cutoff.png", units="in", width=5, height=4, dpi=300)


