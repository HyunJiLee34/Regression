## load csv file
df_young <-read.csv( "/Users/hyunji/Downloads/finalfinal_dfyoung(2).csv",header=T)
df_young
##install
#install.packages("corrplot")
#install.packages("qpcR")
## library
library(corrplot)
library(leaps)
library(lmtest)
library(car)
library(qpcR)
options(max.print=1000000)
## df X,Y지정##

Y <- df_young[,1]
X1 <- df_young[,2] #심혈관질환
X2 <- df_young[,3]# 정신질환
X3 <- df_young[,4]# 주관적 체형인식
X4 <- df_young[,5] #기타지역
X5 <- df_young[,6] #광역시
X6 <- df_young[,7] #여성
X7 <- df_young[,8] #가구소득

full_model <- lm(운동량~심혈관질환+정신질환+주관적.체형인식+기타지역+광역시+female+가구소득,data=df_young)
full_model <- lm(Y~X1+X2+X3+X4+X5+X6+X7)
summary(full_model)

## df spearman 상관관계##
par(mfrow=c(1,1))
X <- cbind(Y,X1,X2,X3,X4,X5,X6,X7)
R <- cor(X,method="spearman");R
mcor<-round(R,3)
corrplot(mcor, method="shade", shade.col=NA, tl.col="black", tl.srt=45, addCoef.col="black")



#최적 모형선택 수정된 결정계수 기준 #
#dfall_reg <- df[,c("운동량","심혈관질환","정신질환","주관적.체형인식","기타지역","광역시","female","가구소득")]
#regfit <- regsubsets(x=운동량~.,data=df_young,method="exhaustive",nbest=7)
#result_regfit <- summary(regfit)
#result_regfit$adjr2 # adj_R2

# plot of adj_R2
#plot(result_regfit$adjr2,pch=19,cex=2,ylab="adj_R2",xlab="model",type="b")
#result_regfit
# final model 부분 f 검정
#final_model_R2 <- lm(운동량 ~ 정신질환+ 기타지역+ 광역시 + female,data=df)
#summary(final_model_R2) 

# comparison with null model using partial-F test
#null_model_R2 <- lm(운동량 ~ 1,data=df)
#anova(final_model_R2,null_model_R2)
#부분f검정에 대한 pvalue가나옴
#귀무가설은 정신질환 기타지역 광역시 female 이 모두 0이다
#pvalue가 작기 때문에기각 따라서 귀무가설은틀렷음  


## forward 운동량 ~ female + 광역시 + 기타지역
null_model <- lm(운동량~1,data=df_young)
fit_forward <- step(null_model, scope = ~ 심혈관질환+정신질환+주관적.체형인식+기타지역+광역시+female+가구소득,direction="forward",test="F")
# backward 운동량 ~ female + 광역시 + 기타지역
fit_backward <-step(full_model,direction="backward",test="F")
# stepwise 운동량 ~ female + 광역시 + 기타지역
fit_stepwise <-step(null_model,scope = ~ 심혈관질환+정신질환+주관적.체형인식+기타지역+광역시+female+가구소득,direction="both",test="F")

#다중공선성#
## df spearman 상관관계##
X <- cbind(Y,X2,X4,X5,X6)
R <- cor(X,method="spearman");R
mcor<-round(R,3)
corrplot(mcor, method="shade", shade.col=NA, tl.col="black", tl.srt=45, addCoef.col="black")

reduced_model <- lm(운동량~정신질환+기타지역+광역시+female,data=df_young)
reduced_model <- lm(Y~X2+X4+X5+X6)
summary(reduced_model)

par(mfrow=c(2,2))
plot(reduced_model)

###잔차###
#hat matrix
X <- as.matrix(cbind(1,df_young[,-1]));X
X <-X[,-2];X
X <-X[,-3];X
X <-X[,-6];X
H <- X%*%solve(t(X)%*%X)%*%t(X)
Y_hat <- predict(reduced_model) 
res <- df_young$운동량-Y_hat;res
#standardized residuals
s <- sqrt(sum(res^2)/(dim(df_young)[1]-5)) # Estimate of sigma 
std_res <- res/(s*sqrt(1-diag(H)));std_res


# residual plot
par(mfrow=c(1,2))
plot(std_res,pch=19,cex=1,ylab="std residual",xlab="Index",ylim=c(-2.5,2.5)) 
abline(h=-2,col="red",lty=2)
abline(h=2,col="red",lty=2)

###지렛점###
leverage <- diag(H);leverage
2*(dim(X)[2])/dim(X)[1]

# leverage plot
plot(leverage,pch=19,cex=1,ylab="leverage value",xlab="Index") 
abline(h=0.01709402,col="blue",lty=2)

###영향력측도###
influence.measures(reduced_model)
cooks.distance(reduced_model) # cooks distance 
dfbetas(reduced_model) # DFBETAS 
dffits(reduced_model) # DFFITS 
covratio(reduced_model) # COVRATIO
# standard of measures
p <- 5
n <- 585
cook_standard <- 3.67/(n-p);cook_standard
DFBETAS_standard <- 2/sqrt(n);DFBETAS_standard 
DFFITS_standard <- 2*sqrt(p/n);DFFITS_standard 
COVRATIO_standard <- 3*p/n; COVRATIO_standard

# 1) Cooks distance
par(mfrow=c(1,1))
plot(cooks.distance(reduced_model),pch=19,ylab="cooks distance") 
abline(h=cook_standard,col="red",lty=2)
dfbetas(full_model)

# 2) DFBETAS
par(mfrow=c(2,2))
plot(dfbetas(reduced_model)[,2],pch=19,ylab="DFBETAS_1") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)
dfbetas(reduced_model)[,2]

plot(dfbetas(reduced_model)[,3],pch=19,ylab="DFBETAS_2") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)

plot(dfbetas(reduced_model)[,4],pch=19,ylab="DFBETAS_3") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)

plot(dfbetas(reduced_model)[,5],pch=19,ylab="DFBETAS_3") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)


# 3) DFFITS
par(mfrow=c(1,1)) 
plot(dffits(reduced_model),pch=19,ylab="DFFITS") 
abline(h=DFFITS_standard,col="red",lty=2) 
abline(h=-DFFITS_standard,col="red",lty=2)

# 4) COVRATIOS
plot(abs(covratio(reduced_model)-1),pch=19,ylab="|COVRATIO-1|")
abline(h=COVRATIO_standard,col="red",lty=2)

df_new<-df_young[-c(487,272,167,344,381,517,169,426,202,511,525,232,538,340,270,284,215,108,385,523),c(1,3,5,6,7)]
df_new


reduced_data_model <- lm(운동량~정신질환+기타지역+광역시+female, data = df_new)
summary(reduced_data_model)
par(mfrow=c(2,2))
plot(reduced_data_model)



###After scaling 잔차###
#hat matrix
X <- as.matrix(cbind(1,df_new[,-1]));X

H <- X%*%solve(t(X)%*%X)%*%t(X)
Y_hat <- predict(reduced_data_model) 
res <- df_new$운동량-Y_hat;res
#standardized residuals
s <- sqrt(sum(res^2)/(dim(df_new)[1]-5)) # Estimate of sigma 
std_res <- res/(s*sqrt(1-diag(H)));std_res


# residual plot
par(mfrow=c(1,2))
plot(std_res,pch=19,cex=1,ylab="std residual",xlab="Index",ylim=c(-2.5,2.5)) 
abline(h=-2,col="red",lty=2)
abline(h=2,col="red",lty=2)

###지렛점###
leverage <- diag(H);leverage
2*(dim(X)[2])/dim(X)[1]

# leverage plot
plot(leverage,pch=19,cex=1,ylab="leverage value",xlab="Index") 
abline(h=0.01769912,col="blue",lty=2)

###영향력측도###
influence.measures(reduced_data_model)
cooks.distance(reduced_data_model) # cooks distance 
dfbetas(reduced_data_model) # DFBETAS 
dffits(reduced_data_model) # DFFITS 
covratio(reduced_data_model) # COVRATIO
# standard of measures
p <- 5
n <- 565
cook_standard <- 3.67/(n-p);cook_standard
DFBETAS_standard <- 2/sqrt(n);DFBETAS_standard 
DFFITS_standard <- 2*sqrt(p/n);DFFITS_standard 
COVRATIO_standard <- 3*p/n; COVRATIO_standard

# 1) Cooks distance
par(mfrow=c(1,1))
plot(cooks.distance(reduced_data_model),pch=19,ylab="cooks distance") 
abline(h=cook_standard,col="red",lty=2)
dfbetas(full_model)

# 2) DFBETAS
par(mfrow=c(2,2))
plot(dfbetas(reduced_data_model)[,2],pch=19,ylab="DFBETAS_1") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)
dfbetas(reduced_data_model)[,2]

plot(dfbetas(reduced_data_model)[,3],pch=19,ylab="DFBETAS_2") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)

plot(dfbetas(reduced_data_model)[,4],pch=19,ylab="DFBETAS_3") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)

plot(dfbetas(reduced_data_model)[,5],pch=19,ylab="DFBETAS_3") 
abline(h=DFBETAS_standard,col="red",lty=2) 
abline(h=-DFBETAS_standard,col="red",lty=2)


#모형의 확인#
set.seed(1234)
rn <- sample(x=c(1:565),size=565,replace=F);rn
df_new$rn <- rn
head(df_new)
# 0.7 0.3
train_dat <- df_new[df_new$rn>396,]
test_dat <- df_new[df_new$rn<=396,] 
# predicted error
train_model <- lm(운동량~정신질환+기타지역+광역시+female ,data=train_dat)
summary(train_model)

predict_value <- predict(train_model,newdata=test_dat[,c('정신질환','기타지역','광역시','female')])
predict_error <- sum((test_dat$운동량-predict_value)^2)
predict_error
#[1] 1881.638
  
#press
final_model <- lm(운동량 ~ 정신질환+ 기타지역+ 광역시 + female,data=df_new)
PRESS_Final <- PRESS(final_model)
PRESS_Final$stat #2583.109


# SST, SSE, SSR
SST <- sum((df_new$운동량-mean(df_new$운동량))^2) 
SSE <- sum(resid(final_model)^2) #[1] 2545
SSR <- SST-SSE
#SSE와 press가 거의 비슷함 -> 예측의 정확도가 좋음 
# R2 vs R2_predic 
1-(SSE/SST)  # R2 = 0.04525281

1-(PRESS_Final$stat/SST) # R2_predic = 0.03100814
# press를 이용하여 구한 r2값ㄷ과 예측값을 구하지 않ㄴ은 r2도 비슷함 
# 설명력이 급격히 떨어지진 x
# 모형이 예측모형의 확인 면에서 좋음
dwtest(reduced_data_model)





