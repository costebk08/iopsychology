######################################
###PredBias Predictive Bias Analysis##
######Authored by Brian Costello######
#Modified by Jeff Cucina for use with#
######local empiricla data files######
######################################

####Load packages
###If any error messages appear associated with 
###these package, then you You need to install them

library(dplyr)
library(readr)
library(haven)
library(tidyr)
library(stringr)
library(openxlsx)
library(readxl)
library(tidyverse)
library("broom")
library("MASS")
library("lm.beta")
library("tidyverse")
library("car")
library(e1071) 
library(effsize)
library(DescTools)
library(rootSolve)

###Define PredBias Function

PredBias <-function(Combined, x, y, output_file_name, 
                    scatterplot_file_name, predictor_histogram_file_name, 
                    criterion_histogram_file_name) {

  Combined$x <- Combined[x]
  Combined$y <- Combined[y]
  Combined$x<-unlist(Combined$x)
  Combined$y<-unlist(Combined$y)
  print(Combined)
  
    print(Combined$x)
    print(Combined$y)

    
Combined<-mutate(Combined,case=1)

reference_dat<-Combined[which(Combined$ref_foc=="a_reference"),]

focal_dat<-Combined[which(Combined$ref_foc=="b_focal"),]
  
####Getting Gulliksen-Wilks info for entire sample####
####Grand means and N
##GMX
GMX<-mean(Combined$x)
##GMY
GMY<-mean(Combined$y)
##N
n_total<-nrow(Combined)

####Getting GW stats
gw_test<-
  ##Using the "combined" dataset
  Combined%>%
  ##Calculations by group
  group_by(ref_foc)%>%
  ##Getting initial summary stats
  summarise(
    #test means
    mean_x=mean(x),
    #criterion means
    mean_y=mean(y),
    #sub-group Ns
    n=sum(case),
    #test SEs
    ex=sum((x-mean_x)^2),
    #criterion SEs
    ey=sum((y-mean_y)^2),
    #covariance Ses
    exy=sum((x-mean_x)*(y-mean_y))
  )%>%
  ##Adding C7
  mutate(
    C7=ey-((exy^2)/ex)
  )%>%
  ##Step 3 metrics
  ungroup()%>%
  mutate(
    C8=sum(C7)/sum(n),
    C9=sum(n)*log(C8),
    C10=C7/n,
    C11=n*log(C10),
    GWA=C9-sum(C11),
    E=sum(exy),
    AB=mean_y-GMY,
    CD=mean_x-GMX,
    SDOT=sum(ey)-((E^2)/sum(ex)),
    GWB=sum(n)*(log(SDOT/sum(n))-log(sum(C7)/sum(n))),
    SA=sum(ey)+sum(n*(AB^2)),
    SB=sum(exy)+sum(n*AB*CD),
    SC=sum(ex)+sum(n*(CD^2)),
    SDASH=SA-((SB^2)/SC),
    GWC=n_total*(log(SDASH/n_total)-log(SDOT/n_total))
  )%>%
  ##3 chi square tests and p values
  summarise(
    GWA=mean(GWA),
    GWA_p=pchisq(GWA,1,lower.tail=FALSE),
    GWB=mean(GWB),
    GWB_p=pchisq(GWB,1,lower.tail=FALSE),
    GWC=mean(GWC),
    GWC_p=pchisq(GWC,1,lower.tail=FALSE),
  )




####Cleary Test####
####Reg 1
reg1<-lm(y~x,data=Combined)
reg1_summary<-summary(reg1)
step1_r2<-reg1_summary$r.squared
step1_beta<-lm.beta(reg1)
##output metrics
step1_r<-sqrt(step1_r2)



####Reg 2
reg2<-lm(y~x+ref_foc,data=Combined)
reg2_summary<-summary(reg2)
step2_r2<-reg2_summary$r.squared
step2_r<-sqrt(step2_r2)
step2_beta<-lm.beta(reg2)
step2_step1_comp<-anova(reg1,reg2)["2",]
##output metrics
step2_delta_r<-step2_r-step1_r
step2_delta_r2<-step2_r2-step1_r2
step2_b_intercept<-step2_beta$coefficients["(Intercept)"]
step2_b_main<-step2_beta$coefficients["ref_focb_focal"]
step2_beta_main<-step2_beta$standardized.coefficients["ref_focb_focal"]
#step2_beta_intercept<-step2_beta$standardized.coefficients["(Intercept)"]
step2_delta_r_f<-step2_step1_comp$`F`
step2_delta_r_p<-step2_step1_comp$`Pr(>F)`


####Reg 3
reg3<-lm(y~x*ref_foc,data=Combined)
reg3_summary<-summary(reg3)
step3_r2<-reg3_summary$r.squared
step3_r<-sqrt(step3_r2)
step3_beta<-lm.beta(reg3)
step3_step2_comp<-anova(reg2,reg3)["2",]
##output metrics
step3_delta_r<-step3_r-step2_r
step3_delta_r2<-step3_r2-step2_r2
step3_b_interaction<-step3_beta$coefficients["x:ref_focb_focal"]
step3_beta_interaction<-step3_beta$standardized.coefficients["x:ref_focb_focal"]
step3_delta_r_f<-step3_step2_comp$`F`
step3_delta_r_p<-step3_step2_comp$`Pr(>F)`


####Effect sizes
## Cohen
cohen_f2<-step3_delta_r2/(1-(step1_r2+step2_delta_r2+step3_delta_r2))
##Liu & Yuan
liu_yuan<-step3_delta_r2/(step1_r2+step3_delta_r2)


####putting results together
cleary_test<-
  tibble(
    step1_r=step1_r,
    step2_delta_r=step2_delta_r,
    step2_delta_r2=step2_delta_r2,
    step2_b_intercept=step2_b_intercept,
    step2_b_main=step2_b_main,
    step2_beta_main=step2_beta_main,
    #step2_beta_intercept=step2_beta_intercept,
    step2_delta_r_f=step2_delta_r_f,
    step2_delta_r_p=step2_delta_r_p,
    step3_delta_r=step3_delta_r,
    step3_delta_r2=step3_delta_r2,
    step3_b_interaction=step3_b_interaction,
    step3_beta_interaction=step3_beta_interaction,
    step3_delta_r_f=step3_delta_r_f,
    step3_delta_r_p=step3_delta_r_p,
    cohen_f2=cohen_f2,
    liu_yuan=liu_yuan
  )




####Separate Regressions for Each Group####
####reference
ref_reg<-lm(y~x,data=reference_dat)
ref_reg_beta<-lm.beta(ref_reg)

ref_constant<-ref_reg_beta$coefficients["(Intercept)"]
ref_b<-ref_reg_beta$coefficients["x"]
ref_beta<-ref_reg_beta$standardized.coefficients["x"]
ref_p<-tidy(ref_reg)%>%filter(term=="x")%>%pull(p.value)


####focal
focal_reg<-lm(y~x,data=focal_dat)
focal_reg_beta<-lm.beta(focal_reg)

focal_constant<-focal_reg_beta$coefficients["(Intercept)"]
focal_b<-focal_reg_beta$coefficients["x"]
focal_beta<-focal_reg_beta$standardized.coefficients["x"]
focal_p<-tidy(focal_reg)%>%filter(term=="x")%>%pull(p.value)


####binding together results
sep_reg_tests<-
  tibble(
    ref_constant=ref_constant,
    ref_b=ref_b,
    ref_beta=ref_beta,
    ref_p=ref_p,
    
    focal_constant=focal_constant,
    focal_b=focal_b,
    focal_beta=focal_beta,
    focal_p=focal_p
  )



####SD(Y) Difference Test####
levene_test_mean<-
  tidy(
    leveneTest(y~as.factor(ref_foc),Combined, center=mean)
  )%>%
  filter(!is.na(statistic))%>%
  transmute(
    levene_mean_f=statistic,
    levene_mean_p=p.value
  )

levene_test_median<-
  tidy(
    leveneTest(y~as.factor(ref_foc),Combined)
  )%>%
  filter(!is.na(statistic))%>%
  transmute(
    levene_median_f=statistic,
    levene_median_p=p.value
  )


n<-sum(Combined$case)

####Educational Measurement Differential Validity/Prediction Approach Output (residual metrics)####
residual_metrics<-
  tibble(
    ref_foc=Combined$ref_foc,
    residuals=step1_beta$residuals,
    residuals_z=scale(residuals)[1:(n)]
  )%>%
  mutate(
    ref_foc=ifelse(ref_foc=="a_reference","reference","focal")
  )%>%
  group_by(ref_foc)%>%
  summarise(
    residuals_m=mean(residuals),
    residuals_sd=sd(residuals),
    residuals_z_m=mean(residuals_z),
    residuals_z_sd=sd(residuals_z)
  )%>%
  gather(
    "key",
    "value",
    -ref_foc
  )%>%
  transmute(
    metric=paste(ref_foc,key,sep="_"),
    value=value
  )%>%
  spread(
    metric,
    value
  )

##Standard Error of Z-test comparing two independent correlations
n_reference<-nrow(reference_dat)
n_focal<-nrow(focal_dat)
se_z_corr<-sqrt( (1/(n_focal-3)) + (1/(n_reference-3)) )
se_z_corr_results<-tibble(
  se_z_corr=se_z_corr
)

####Fisher's Z test
##Fisher z-test for ref and focal
fisher_z_reference<-DescTools::FisherZ(ref_beta)
fisher_z_focal<-DescTools::FisherZ(focal_beta)
##z-test comparing correlations
fisher_z_corr_compare_z<-(fisher_z_reference-fisher_z_focal)/se_z_corr
fisher_z_corr_compare_p<-2*(1-pnorm(abs(fisher_z_corr_compare_z)))
fisher_z_corr_compare_q<-fisher_z_reference-fisher_z_focal
##assemble results
fisher_z_test<-
  tibble(
    fisher_z_reference=fisher_z_reference,
    fisher_z_focal=fisher_z_focal,
    fisher_z_corr_compare_z=fisher_z_corr_compare_z,
    fisher_z_corr_compare_p=fisher_z_corr_compare_p,
    fisher_z_corr_compare_q=fisher_z_corr_compare_q
  )




####Getting Descriptives for separate samples####



MeanXreference<-mean(reference_dat$x)
MeanYreference<-mean(reference_dat$y)
MeanXfocal<-mean(focal_dat$x)
MeanYfocal<-mean(focal_dat$y)

MinXreference<-min(reference_dat$x)
MinYreference<-min(reference_dat$y)
MinXfocal<-min(focal_dat$x)
MinYfocal<-min(focal_dat$y)

MaxXreference<-max(reference_dat$x)
MaxYreference<-max(reference_dat$y)
MaxXfocal<-max(focal_dat$x)
MaxYfocal<-max(focal_dat$y)

SDXreference<-sd(reference_dat$x)
SDYreference<-sd(reference_dat$y)
SDXfocal<-sd(focal_dat$x)
SDYfocal<-sd(focal_dat$y)

skewnessXreference<-skewness(reference_dat$x)
skewnessYreference<-skewness(reference_dat$y)
skewnessXfocal<-skewness(focal_dat$x)
skewnessYfocal<-skewness(focal_dat$y)

kurtosisXreference<-kurtosis(reference_dat$x)
kurtosisYreference<-kurtosis(reference_dat$y)
kurtosisXfocal<-kurtosis(focal_dat$x)
kurtosisYfocal<-kurtosis(focal_dat$y)

####t-tests
t_test_results<-
  ##predictor
  tidy(t.test(x~ref_foc,data=Combined))%>%
  transmute(
    t_test_x=statistic,
    t_test_x_p=p.value
  )%>%
  bind_cols(
    ##criterion
    tidy(t.test(y~ref_foc,data=Combined))%>%
      transmute(
        t_test_y=statistic,
        t_test_y_p=p.value
      )
  )


####effect size calculations for predictor and criterion
##predictor

cohen_d_x<-cohen.d(as.numeric(unlist(reference_dat[x])), as.numeric(unlist(focal_dat[x])))

hedges_g_x<-cohen.d(as.numeric(unlist(reference_dat[x])), as.numeric(unlist(focal_dat[x])), hedges.correction=TRUE)

cohen_d_y<-cohen.d(as.numeric(unlist(reference_dat[y])), as.numeric(unlist(focal_dat[y])))

hedges_g_y<-cohen.d(as.numeric(unlist(reference_dat[y])), as.numeric(unlist(focal_dat[y])), hedges.correction=TRUE)


##assemble results
cohen_d_x_results<-
  tibble(
    cohen_d_x_estimate=cohen_d_x$estimate,
    cohen_d_x_within_group_sd=cohen_d_x$sd,
        )

cohen_d_y_results<-
  tibble(
    cohen_d_y_estimate=cohen_d_y$estimate,
    cohen_d_y_within_group_sd=cohen_d_y$sd,
  )

hedges_g_x_results<-
  tibble(
    hedges_g_x_estimate=hedges_g_x$estimate
      )

hedges_g_y_results<-
  tibble(
    hedges_g_y_estimate=hedges_g_y$estimate
  )

####tests for Homoscedasticity
homoscedasticity_tests<-
  tidy(lmtest::bptest(reg1))%>%
  transmute(
    breusch_pagan_lm1=statistic,
    breusch_pagan_p_lm1=p.value
  )%>%
  bind_cols(
    tidy(lmtest::bptest(reg2))%>%
      transmute(
        breusch_pagan_lm2=statistic,
        breusch_pagan_p_lm2=p.value
      )
  )%>%
  bind_cols(
    tidy(lmtest::bptest(reg3))%>%
      transmute(
        breusch_pagan_lm3=statistic,
        breusch_pagan_p_lm3=p.value
      )
  )




####binding together results

referenceValidity<-ref_beta

focalValidity<-focal_beta

sep_descriptives<-
  tibble(
    n_reference,
    MeanXreference,
    SDXreference,
    MinXreference,
    MaxXreference,
    skewnessXreference,
    kurtosisXreference,
    MeanYreference,
    SDYreference,
    MinYreference,
    MaxYreference,
    skewnessYreference,
    kurtosisYreference,
    referenceValidity,
    n_focal,
    MeanXfocal,
    SDXfocal,
    MinXfocal,
    MaxXfocal,
    skewnessXfocal,
    kurtosisXfocal,
    MeanYfocal,
    SDYfocal,
    MinYfocal,
    MaxYfocal,
    skewnessYfocal,
    kurtosisYfocal,
    focalValidity,
      )

####getting final data for output####
final_dat<-
  bind_cols(x,
    y,        
    sep_descriptives,
    t_test_results,
    cohen_d_x_results,
    hedges_g_x_results,
    cohen_d_y_results,
    hedges_g_y_results,
    gw_test,
    cleary_test,
    sep_reg_tests,
    fisher_z_test,
    se_z_corr_results,
    residual_metrics,
    levene_test_mean,
    levene_test_median,
    homoscedasticity_tests
  )

print(final_dat)

##write output file *NOTE* you can change this to whatever file name you want to use
#output_file_name

write_csv(final_dat, output_file_name)


####Generating scatterplot####


####Get stats
##reference
ref_r<-round((ref_reg_beta$standardized.coefficients["x"]),digits=2)
ref_x<-round(ref_reg$coefficients["x"],digits=2)
ref_intercept<-round(ref_reg$coefficients["(Intercept)"],digits=2)

##focal
focal_r<-round((focal_reg_beta$standardized.coefficients["x"]),digits=2)
focal_x<-round(focal_reg$coefficients["x"],digits=2)
focal_intercept<-round(focal_reg$coefficients["(Intercept)"],digits=2)

##open a png
png(
  filename=paste0(
    getwd(),"/",
    scatterplot_file_name
  ), 
  width=480,
  height=480,
  pointsize=12,
  bg="white",
  res=NA
)
##make graph
scatterplot<-
  PredBiasData%>%
  ggplot(aes(x=test,y=criterion,color=`ref_foc`))+
  geom_point()+ 
  geom_smooth(method=lm,se=FALSE,fullrange=TRUE)+
  labs(caption=paste0("Ref: R = ",ref_r,"; y = ",ref_x,"*x + ",ref_intercept,
                                      "; ", "Focal: R = ",focal_r,"; y = ",focal_x,"*x + ",focal_intercept)
  )
  
print(scatterplot + theme_bw())
##close png
dev.off()

####Histograms####

##data cleaning
combined_dat1<-Combined%>%
  mutate(
    ref_foc=case_when(
      ref_foc=="a_reference"~"Reference",
      ref_foc=="b_focal"~"Focal"
    )
  )%>%
  rename(`Focal/Reference`=ref_foc)


####predictor histogram
##open a png
png(
  filename=paste0(
    getwd(),"/",
    predictor_histogram_file_name
  ), 
  width=480,
  height=480,
  pointsize=12,
  bg="white",
  res=NA
)

hist_x<-combined_dat1%>%
  ggplot(aes(x=x, fill=`Focal/Reference`, color=`Focal/Reference`)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(
    title="Predictor Histogram",
    x="Predictor Values",
    y="Count"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
print(hist_x)

dev.off()


####criterion histogram
##open a png
png(
  filename=paste0(
    getwd(),"/",
    criterion_histogram_file_name
  ), 
  width=480,
  height=480,
  pointsize=12,
  bg="white",
  res=NA
)

hist_y<-combined_dat1%>%
  ggplot(aes(x=y, fill=`Focal/Reference`, color=`Focal/Reference`)) +
  geom_histogram(position="identity", alpha=0.5) +
  labs(
    title="Criterion Histogram",
    x="Criterion Values",
    y="Count"
  )+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))
print(hist_y)
dev.off()


}