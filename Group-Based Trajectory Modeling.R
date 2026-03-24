library(gbmt)
library(lcmm)
library(tidyr)
library(dplyr)
library(ggplot2)
#数据清洗######
data <- read.csv('~/Documents/r/sepsis/sepsis.csv',header = T)
data <- data %>%
  filter(transplant == 0,
    longterm_steroid == 0,
    anti_cancer_therapy == 0,
    hematology_disease == 0,
    aids == 0,
    metastatic_solid_tumor == 0,
    malignant_cancer == 0,
    rheumatic_disease == 0
  )  %>%
  select(id=stay_id,
         gender = gender,
         age = admission_age,
         los_hospital = los_hospital,
         los_icu = los_icu,
         sofa = sofa_score,
         gcs = gcs_min,
         apsiii = apsiii,
         lods = lods,
         oasis = oasis,
         sapsii = sapsii,
         sirs = sirs,
         meld = meld,
         aki = aki_stage,
         myocardial_infarct = myocardial_infarct,
         congestive_heart_failure = congestive_heart_failure,
         peripheral_vascular_disease = peripheral_vascular_disease,
         cerebrovascular_disease = cerebrovascular_disease,
         dementia = dementia,
         chronic_pulmonary_disease = chronic_pulmonary_disease,
         peptic_ulcer_disease = peptic_ulcer_disease,
         mild_liver_disease = mild_liver_disease,
         severe_liver_disease = severe_liver_disease,
         diabetes_without_cc = diabetes_without_cc,
         diabetes_with_cc = diabetes_with_cc,
         paraplegia = paraplegia,
         charlson_comorbidity_index = charlson_comorbidity_index,
         hospital_outcome = hospital_outcome,
         icu_outcome = icu_outcome,
         survival_time_28d = survival_time_28d,
         lactate = lactate_mean,
         lymphocyte_1 = lymphocyte_1d_min,
         lymphocyte_2 = lymphocyte_2d_min,
         lymphocyte_3 = lymphocyte_3d_min,
         lymphocyte_4 = lymphocyte_4d_min,
         lymphocyte_5 = lymphocyte_5d_min,
         lymphocyte_6 = lymphocyte_6d_min,
         lymphocyte_7 = lymphocyte_7d_min
  )
data <- data %>%
  filter(rowSums(!is.na(select(., lymphocyte_1:lymphocyte_7))) >= 2) %>%
  filter(!apply(select(., lymphocyte_1:lymphocyte_7), 1, function(x) any(x[!is.na(x)] > 10))) %>%
  mutate(outcome = ifelse(survival_time_28d >= 28, 0, 1))

#分组分段线图1-7Day######
mean_outcome_1 <- data %>%
  filter(outcome == 1) %>%
  summarise(across(starts_with("lymphocyte_"), mean, na.rm = TRUE))
mean_outcome_0 <- data %>%
  filter(outcome == 0) %>%
  summarise(across(starts_with("lymphocyte_"), mean, na.rm = TRUE))
mean_outcome_1_long <- mean_outcome_1 %>%
  pivot_longer(cols = everything(), names_to = "Day", values_to = "Mean") %>%
  mutate(Group = "Outcome 1")
mean_outcome_0_long <- mean_outcome_0 %>%
  pivot_longer(cols = everything(), names_to = "Day", values_to = "Mean") %>%
  mutate(Group = "Outcome 0")
mean_data_long <- bind_rows(mean_outcome_1_long, mean_outcome_0_long)
ggplot(mean_data_long, aes(x = Day, y = Mean, color = Group, group = Group)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Lymphocyte Counts Over 7 Days",
       x = "Day",
       y = "Average Lymphocyte Count") +
  theme_minimal() +
  scale_x_discrete(labels = c("lymphocyte_1" = "Day 1", "lymphocyte_2" = "Day 2", "lymphocyte_3" = "Day 3",
                              "lymphocyte_4" = "Day 4", "lymphocyte_5" = "Day 5", "lymphocyte_6" = "Day 6",
                              "lymphocyte_7" = "Day 7"))
#GBTM#####
set.seed(123)
data_long <- data %>%
  pivot_longer(
    cols = starts_with("lymphocyte_"),
    names_to = "time_point",
    names_prefix = "lymphocyte_",
    values_to = "lymphocyte"
  ) %>%
  mutate(
    time_point = as.numeric(time_point)
  )
df<-data_long
df<-as.data.frame(df)
varNames <- c("lymphocyte")
#分别建立1组、2组、3组、4组、5组的3阶模型，scaling=0表示Y不做任何转换。
m1_3 <- gbmt(x.names=varNames, unit="id", time="time_point", d=3, ng=1, data=df,scaling=0)
m2_3 <- gbmt(x.names=varNames, unit="id", time="time_point", d=3, ng=2, data=df,scaling=0)
m3_3 <- gbmt(x.names=varNames, unit="id", time="time_point", d=3, ng=3, data=df,scaling=0)
m4_3 <- gbmt(x.names=varNames, unit="id", time="time_point", d=3, ng=4, data=df,scaling=0)
m5_3 <- gbmt(x.names=varNames, unit="id", time="time_point", d=3, ng=5, data=df,scaling=0)
#模型评价#####
ICs<-rbind(m1_3$ic,m2_3$ic, m3_3$ic,m4_3$ic,m5_3$ic)
ICs<-as.data.frame(ICs)
#分别显示三组的average posterior probability.
data1 <- m1_3$appa
data2 <- m2_3$appa
data3 <- m3_3$appa
data4 <- m4_3$appa
data5 <- m5_3$appa
max_length <- max(length(data1), length(data2), length(data3), length(data4), length(data5))
traj_1 <- c(data1, rep(NA, max_length - length(data1)))
traj_2 <- c(data2, rep(NA, max_length - length(data2)))
traj_3 <- c(data3, rep(NA, max_length - length(data3)))
traj_4 <- c(data4, rep(NA, max_length - length(data4)))
traj_5 <- c(data5, rep(NA, max_length - length(data5)))
avePP <- data.frame(traj_1,traj_2, traj_3, traj_4,traj_5 )
avePP_t <- t(avePP)
colnames(avePP_t) <- c("G1", "G2", "G3","G4","G5")
avePP_t<-as.data.frame(avePP_t)
avePP_t$model <- row.names(avePP_t)
row.names(avePP_t) <- NULL
avePP_t<-avePP_t[,c(6,1,2,3,4,5)]
##把信息准则和平均后验概率合并
table<-cbind(avePP_t,ICs)
print(table)
##绘制轨迹图######
##绘制轨迹图——二轨迹
par(bty = "n")
plot(m2_3,bands=F,xlab="time",ylim=c(0.5,2),ylab = "Lymphocyte counts",titles = "Lymphocyte trajectories",transparency=75,
     add.grid = F,add.legend = FALSE,equal.scale = T,xaxt = "n",col=c("#CCCCFF", "#CCCCCC"))
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("topright", legend = c("trajectory 1", "trajectory 2"), 
       col = c("#CCCCFF", "#CCCCCC"), lty = 1, lwd = 2, bty = "n")
##绘制轨迹图——三轨迹
par(bty = "n")
plot(m3_3,bands=T,xlab="time",ylim=c(0.5,2),ylab = "Lymphocyte counts",titles = "Lymphocyte trajectories",transparency=75,
     add.grid = F,add.legend = FALSE,equal.scale = T,xaxt = "n",col=c("#CCCCFF", "#CCCCCC","#CCCC99"))#n.ahead = 3
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("top", legend = c("trajectory 1: High level-Rapid drop-Fluctuation", 
                         "trajectory 2: Stable",
                         "trajectory 3: Low level-Slow increase"), 
       col = c("#CCCCFF", "#CCCCCC","#CCCC99"), lty = 1, lwd = 2, bty = "n")
##绘制轨迹图——四轨迹
par(bty = "n")
plot(m4_3,bands=F,xlab="time",ylim=c(0.5,2),ylab = "Lymphocyte counts",titles = "Lymphocyte trajectories",transparency=75,
     add.grid = F,add.legend = FALSE,equal.scale = T,xaxt = "n",col=c("#CCCCFF", "#CCCCCC","#CCCC99","purple"))
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("topright",legend = c("trajectory 1", "trajectory 2","trajectory 3","trajectory 4"), 
       col = c("#CCCCFF", "#CCCCCC","#CCCC99","purple"), lty = 1, lwd = 2, bty = "n")
##绘制轨迹图——五轨迹
par(bty = "n")
plot(m5_3,bands=F,ylim=c(0.5,2),xlab="time",ylab = "Lymphocyte counts",titles = "Lymphocyte trajectories",transparency=75,
     add.grid = F,add.legend = FALSE,equal.scale = T,xaxt = "n",col=c("#CCCCFF", "#CCCCCC","#CCCC99","purple",'red'))
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("topright", legend = c("trajectory 1", "trajectory 2","trajectory 3","trajectory 4",'trajectory 5'), 
       col = c("#CCCCFF", "#CCCCCC","#CCCC99","purple",'red'), lty = 1, lwd = 2, bty = "n")#inset = c(0, -0.05),y.intersp =0.8
#根据模型的assign参数，对data的每个id赋予相应的轨迹分类，进行后续的回归分析
# head(m3_3$assign)
# 30000484 30004018 30005366 30007565 30015288 30020307 
# 1        1        2        2        2        3 

#HLME######
set.seed(123)
m1 = hlme(fixed =lymphocyte~poly(time_point,degree = 3,raw = TRUE),
          random=~poly(time_point,degree = 3,raw = TRUE),
          subject = 'id',
          ng=1,
          data=df,nproc=5)
m2 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(fixed =lymphocyte~poly(time_point,degree = 3,raw = TRUE),
                      random=~poly(time_point,degree = 3,raw = TRUE),
                      mixture =~poly(time_point,degree = 3,raw = TRUE), 
                      subject = 'id',
                      ng=2,
                      data=df,nproc=5))
m3 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(fixed =lymphocyte~poly(time_point,degree = 3,raw = TRUE),
                      random=~poly(time_point,degree = 3,raw = TRUE),
                      mixture =~poly(time_point,degree = 3,raw = TRUE), 
                      subject = 'id',
                      ng=3,
                      data=df,nproc=5))
m4 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(fixed =lymphocyte~poly(time_point,degree = 3,raw = TRUE),
                      random=~poly(time_point,degree = 3,raw = TRUE),
                      mixture =~poly(time_point,degree = 3,raw = TRUE), 
                      subject = 'id',
                      ng=4,
                      data=df,nproc=5))
m5 <- gridsearch(rep = 5, maxiter = 10, minit = m1,
                 hlme(fixed =lymphocyte~poly(time_point,degree = 3,raw = TRUE),
                      random=~poly(time_point,degree = 3,raw = TRUE),
                      mixture =~poly(time_point,degree = 3,raw = TRUE), 
                      subject = 'id',
                      ng=5,
                      data=df,nproc=5))
save(m1,m2,m3,m4,m5,file = 'null_steroid_7d_3D_hlme_models.RData')
summarytable(m1,m2,m3,m4,m5,
             which = c("G", "loglik", "conv", "npm", "AIC", 
                       "BIC", "SABIC", "entropy", "%class"))
plot(m3, which = "fit", var.time = "time_point", 
     bty = "l", ylab = "Lymphocyte",
     xlab = "Time points", lwd = 2,
     marg = T,legend=NULL,shades = T,xaxt = "n")
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("topright", legend = c("trajectory 1", "trajectory 2",'trajectory 3'), 
       col = 1:3, lty = 1, lwd = 2, bty = "n")
#根据模型的prrob数据框内的class参数，对data的每个id赋予相应的轨迹分类，进行后续的回归分析
# head(m3$pprob)
# id class         prob1     prob2        prob3
# 1 30000484     2 9.491845e-121 0.9999995 5.478746e-07
# 2 30005366     2 4.404715e-156 1.0000000 8.539077e-11
# 3 30007565     2 7.978929e-161 1.0000000 1.318323e-13
# 4 30020307     2 1.023600e-144 0.9971219 2.878100e-03
# 5 30021727     2 1.754165e-140 1.0000000 1.718468e-10
# 6 30030798     2 8.514781e-123 1.0000000 9.618152e-11

#JLCM#####
library(lcmm)
set.seed(123)
jlcm1 <- Jointlcmm(
  fixed = lymphocyte ~ poly(time_point,degree = 3,raw = TRUE),
  random = ~ poly(time_point,degree = 3,raw = TRUE),
  subject = "id",
  ng = 1,
  data = df,
  survival = Surv(survival_time_28d, outcome) ~ sofa+sirs+charlson_comorbidity_index,
  nproc = 8
)
summary(jlcm1)

jlcm2 <- gridsearch(rep = 5, maxiter = 10, minit = jlcm1,Jointlcmm(
  fixed = lymphocyte ~ poly(time_point,degree = 3,raw = TRUE),
  random = ~ poly(time_point,degree = 3,raw = TRUE),
  mixture =~ poly(time_point,degree = 3,raw = TRUE),
  subject = "id",
  ng = 2,
  data = df,
  survival = Surv(survival_time_28d, outcome) ~ sofa+sirs+charlson_comorbidity_index,
  nproc = 8
))

jlcm3 <- gridsearch(rep = 5, maxiter = 10, minit = jlcm1,Jointlcmm(
  fixed = lymphocyte ~ poly(time_point,degree = 3,raw = TRUE),
  random = ~ poly(time_point,degree = 3,raw = TRUE),
  mixture =~ poly(time_point,degree = 3,raw = TRUE),
  subject = "id",
  ng = 3,
  data = df,
  survival = Surv(survival_time_28d, outcome) ~ sofa+sirs+charlson_comorbidity_index,
  nproc = 8
))

jlcm4 <- gridsearch(rep = 5, maxiter = 10, minit = jlcm1,Jointlcmm(
  fixed = lymphocyte ~ poly(time_point,degree = 3,raw = TRUE),
  random = ~ poly(time_point,degree = 3,raw = TRUE),
  mixture =~ poly(time_point,degree = 3,raw = TRUE),
  subject = "id",
  ng = 4,
  data = df,
  survival = Surv(survival_time_28d, outcome) ~ sofa+sirs+charlson_comorbidity_index,
  nproc = 8
))

jlcm5 <- gridsearch(rep = 5, maxiter = 10, minit = jlcm1,Jointlcmm(
  fixed = lymphocyte ~ poly(time_point,degree = 3,raw = TRUE),
  random = ~ poly(time_point,degree = 3,raw = TRUE),
  mixture =~ poly(time_point,degree = 3,raw = TRUE),
  subject = "id",
  ng = 5,
  data = df,
  survival = Surv(survival_time_28d, outcome) ~ sofa+sirs+charlson_comorbidity_index,
  nproc = 8
))
save(jlcm1,jlcm2,jlcm3,jlcm4,jlcm5,file = 'null_steroid_7d_3D_jlcm_min_models.RData')
summarytable(jlcm1,jlcm2,jlcm3,jlcm4,jlcm5, which = c("G", "loglik", "conv", "npm", "AIC", 
                                                      "BIC", "SABIC", "entropy", "%class"))

plot(jlcm2, which = "fit", var.time = "time_point", 
     bty = "l", ylab = "Lymphocyte",
     xlab = "Time points", lwd = 2,
     marg = T,legend=NULL,shades = T,xaxt = "n")
axis(side = 1, at = c(1,2,3,4,5,6,7), labels = c("Day1", "Day2","Day3","Day4","Day5",'Day6','Day7'), las = 1)
legend("topright", legend = c("trajectory 1", "trajectory 2",'trajectory 3'), 
       col = 1:3, lty = 1, lwd = 2, bty = "n")

#根据模型的prrob数据框内的class参数，对data的每个id赋予相应的轨迹分类，进行后续的回归分析
# head(jlcm3$pprob)
# id class    probYT1      probYT2      probYT3
# 1 30000484     2 0.00394447 0.9958241388 0.0002313909
# 2 30004018     3 0.00222400 0.0001336513 0.9976423483
# 3 30005366     1 0.51066730 0.2626743983 0.2266583011
# 4 30007565     1 0.50037215 0.3534362786 0.1461915697
# 5 30015288     1 0.48406005 0.3964752991 0.1194646538
# 6 30020307     1 0.56693478 0.0308405533 0.4022246659

