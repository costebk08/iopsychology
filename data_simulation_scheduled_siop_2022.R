#######################################################################
##########DATA SIMULATION WITH DIFFERENTIAL PREDICTION ANALYSIS########
#######################################################################
options(java.parameters = "Inf")


####packages
library(dplyr)
library(broom)
library(readr)
library(MASS)
library(lm.beta)
library(tidyverse)
#library(car)
#library(linpk)
library(tictoc)
library(parallel)
#library(ggplot2)
####Setting up arguments to feed function####
####To control randomization
#set.seed(888)


####custom user arguments####
####non-varying
##out of function arguments 
output_location<-"C:/Users/bcostell/Desktop/Project-Tickets/Monte Carlo Data Simulation/data_output"
do_you_want_to_use_parallel_processing<-"No"
do_you_want_to_output_sim_conditions<-"No"

##in function arguments
reference_m_x<-0
reference_sd_x<-1
reference_sd_y<-1
focal_m_x<-0
focal_sd_x<-1
iteration<-1
do_you_want_to_output_simulated_data<-"Yes"
do_you_want_to_output_scatterplots<-"Yes"
do_you_want_to_output_histograms<-"Yes"
do_you_want_to_output_residual_plots<-"Yes"
x_var_name<-"x"
y_var_name<-"y"


####varying
##validity combos
#reference_validity<-seq(0,0.5,by=0.1)
#focal_validity<-seq(0,0.5,by=0.1)
reference_validity<-seq(0.1,0.2,by=0.1)
focal_validity<-seq(0.1,0.2,by=0.1)

##focal sd y
#focal_sd_y<-seq(0.5,1.5,by=0.1)
focal_sd_y<-seq(0.9,1.1,by=0.1)

##Intercept differences
#intercept_diffs<-seq(-0.3,0.3,by=0.1)
intercept_diffs<-seq(0.2,0.3,by=0.1)

##D
#D<-c(1.0,0.7)
D<-c(0.7)

##Sample size
#n<-c(500,1000,2000)
n<-c(5000)

##sample proportion conditions
#population_percentage<-seq(0.05,0.95,by=0.05)
population_percentage<-seq(0.5,0.5,by=0.05)


####Creating condition matrix
conditions_iterations<-
  expand_grid(
    
    ##custom non-varying arguments
    reference_m_x=reference_m_x,
    reference_sd_x=reference_sd_x,
    reference_sd_y=reference_sd_y,
    x_var_name=x_var_name,
    y_var_name=y_var_name,
    
    focal_m_x=focal_m_x,
    focal_sd_x=focal_sd_x,
    do_you_want_to_output_simulated_data=do_you_want_to_output_simulated_data,
    do_you_want_to_output_scatterplots=do_you_want_to_output_scatterplots,
    do_you_want_to_output_histograms=do_you_want_to_output_histograms,
    do_you_want_to_output_residual_plots=do_you_want_to_output_residual_plots,
    
    ##varying user-specified arguments
    reference_validity=reference_validity,
    focal_validity=focal_validity,
    focal_sd_y=focal_sd_y,
    intercept_diffs=intercept_diffs,
    D=D,
    n=n,
    population_percentage=population_percentage,
    
    ##number of iterations
    iteration=1:iteration
  )%>%
  
  ##get each unique condition
  group_by_at(vars(iteration))%>%
  mutate(condition=1:n())%>%
  ungroup()%>%
  mutate(z=1:n())




####function to simulate data####
####create output location
date_cleaned<-
  gsub(
    "\\-",
    "_",
    Sys.Date()
  )

##create folders dynamically based on date
folder_length<-str_split(output_location, pattern = "/")%>%
  unlist%>%length()

create_output_folders<-
  function(x, folder_length){
    if(
      !dir.exists(path=paste0(x,"/",date_cleaned,"_1"))
    ){
      
      dir.create(paste0(x,"/",date_cleaned,"_1"))
      
    }else{
      
        
      dir.create(
        paste0(
          x,"/",date_cleaned,"_",
          max(as.numeric(
            sapply(
              str_split(
                sapply(
                  str_split(
                    grep(
                      date_cleaned,
                      list.dirs(path=paste(x,sep="/")),
                      value=TRUE
                    ),
                    pattern="/"
                  ),
                  "[",folder_length + 1
                ),
                pattern="_"
              ),
              "[",4
            )
          ))+1
        )
      )
    }
  }

lapply(output_location, create_output_folders, folder_length)

##set date of interest dynamically
date_of_interest<-
  sapply(
    str_split(
      grep(
        date_cleaned,
        list.dirs(
          path=paste(
            output_location,sep="/")
        ),
        value=TRUE
      ),
      pattern="/"
    ),
    "[", folder_length + 1
  )%>%max()

##get file output folder
file_output_folder<-
  paste0(output_location, "/", date_of_interest)

##output sim conditions
if(do_you_want_to_output_sim_conditions=="Yes"){
  readr::write_csv(
    conditions_iterations,
    paste0(
      file_output_folder,
      "/sim_conditions_",
      date_of_interest,
      ".csv"
    )
  )
}




####function for simulation and analysis####
####for showing example
#z<-24
#dat<-conditions_iterations


####start the function
monte_carlo_data_simulation<-
  function(z, dat, date_of_interest, output_location, file_output_folder){
    
    ####varying arguments
    conditions1<-dat[z,]
    
    condition<-conditions1$condition
    iteration<-conditions1$iteration
    reference_validity<-conditions1$reference_validity
    focal_validity<-conditions1$focal_validity
    focal_sd_y<-conditions1$focal_sd_y
    intercept_diffs<-conditions1$intercept_diffs
    D<-conditions1$D
    #include_in_analysis<-conditions1$include_in_analysis
    n<-conditions1$n
    population_percentage<-conditions1$population_percentage
    focal_n<-round(population_percentage*n,digits = 0)
    reference_n<-n-focal_n
    
    ####fixed arguments
    do_you_want_to_output_simulated_data<-conditions1$do_you_want_to_output_simulated_data
    do_you_want_to_output_scatterplots<-conditions1$do_you_want_to_output_scatterplots
    do_you_want_to_output_histograms<-conditions1$do_you_want_to_output_histograms
    do_you_want_to_output_residual_plots<-conditions1$do_you_want_to_output_residual_plots
    reference_m_x<-conditions1$reference_m_x
    reference_sd_x<-conditions1$reference_sd_x
    reference_sd_y<-conditions1$reference_sd_y
    
    focal_m_x<-conditions1$focal_m_x
    focal_sd_x<-conditions1$focal_sd_x
    x_var_name<-conditions1$x_var_name
    y_var_name<-conditions1$y_var_name
    
    
    ####Getting reference and focal data
    combined_dat<-
      ##reference
      tibble(
        ref_foc="a_reference",
        e=rnorm(reference_n,0,1),
        x=rnorm(reference_n,reference_m_x,reference_sd_x)
      )%>%
      mutate(
        #x=x-0,
        y=reference_sd_y*(0+reference_validity*x+sqrt(1-reference_validity^2)*e)
      )%>%
      bind_rows(
        
        ##getting focal dat
        tibble(
          ref_foc="b_focal",
          e=rnorm(focal_n,0,1),
          x=rnorm(focal_n,focal_m_x,focal_sd_x)
        )%>%
          mutate(
            x=x-D,
            y=focal_sd_y*(intercept_diffs+focal_validity*x+sqrt(1-focal_validity^2)*e)
          )
      )%>%
      mutate(
        condition=condition,
        iteration=iteration,
        dummy_column=1
      )%>%
      
      ##dynamically convert predictor names to x and y
      #x
      rename_at(vars(all_of(x_var_name)),list(~gsub(x_var_name,"x",.))) %>% 
      #y
      rename_at(vars(all_of(y_var_name)),list(~gsub(y_var_name,"y",.)))
  
    
  
      
    #####observed sample stats####
    manipulation_check<-
      combined_dat%>%
      group_by(ref_foc)%>%
      summarise(
        x_m=mean(x),
        x_sd=sd(x),
        x_min=min(x),
        x_max=max(x),
        x_skewness=moments::skewness(x),
        x_kurtosis=moments::kurtosis(x),
        y_m=mean(y),
        y_sd=sd(y),
        y_min=min(y),
        y_max=max(y),
        y_skewness=moments::skewness(y),
        y_kurtosis=moments::kurtosis(y),
        x_y_cor=cor(x,y),
        .groups = "drop"
      )%>%
      gather(
        "metric",
        "value",
        -ref_foc
      )%>%
      transmute(
        ref_foc_metric=paste(ref_foc,metric,sep="_"),
        value=value,
        ref_foc_metric=
          paste0(
            "observed_",
            str_remove_all(ref_foc_metric,"a_|b_")
          )
      )%>%
      spread(
        ref_foc_metric,
        value
      )%>%
      mutate(
        observed_pooled_sd_x=
          sqrt((observed_reference_x_sd^2+observed_focal_x_sd^2)/2),
        observed_d=
          (observed_reference_x_m-observed_focal_x_m)/observed_pooled_sd_x
      )
    
    
    ####t-tests
    t_test_results<-
      ##predictor
      tidy(t.test(x~ref_foc,data=combined_dat))%>%
      transmute(
        t_test_x=statistic,
        t_test_x_p=p.value
      )%>%
      bind_cols(
        ##criterion
        tidy(t.test(y~ref_foc,data=combined_dat))%>%
          transmute(
            t_test_y=statistic,
            t_test_y_p=p.value
          )
        )
    
    
    ####effect size calculations for predictor and criterion
    ##predictor
    sd_within_group_x<-
      sqrt(
        (
          (focal_n-1)*(manipulation_check$observed_reference_x_sd)*
            (manipulation_check$observed_reference_x_sd)+(reference_n-1)*
            (manipulation_check$observed_reference_x_sd)*(manipulation_check$observed_reference_x_sd)
          )/((focal_n-1)+(reference_n-1))
        )
    cohens_d_x<-
      (manipulation_check$observed_reference_x_m-manipulation_check$observed_focal_x_m)/sd_within_group_x
    hedges_g_x<-cohens_d_x*(1-3/(4*(reference_n+focal_n)-9))
    
    ##criterion
    sd_within_group_y<-
      sqrt(
        (
          (focal_n-1)*(manipulation_check$observed_reference_y_sd)*
            (manipulation_check$observed_reference_y_sd)+(reference_n-1)*
            (manipulation_check$observed_reference_y_sd)*(manipulation_check$observed_reference_y_sd)
        )/((focal_n-1)+(reference_n-1))
      )
    cohens_d_y<-
      (manipulation_check$observed_reference_y_m-manipulation_check$observed_focal_y_m)/sd_within_group_y
    hedges_g_y<-cohens_d_y*(1-3/(4*(reference_n+focal_n)-9))
    
    ##assemble together
    effect_size_calcs<-
      tibble(
        sd_within_group_x=sd_within_group_x,
        cohens_d_x=cohens_d_x,
        hedges_g_x=hedges_g_x,
        sd_within_group_y=sd_within_group_y,
        cohens_d_y=cohens_d_y,
        hedges_g_y=hedges_g_y
      )
    
    
    ####Getting Gulliksen-Wilks info for entire sample####
    ####Grand means and N
    ##GMX
    GMX<-mean(combined_dat$x)
    ##GMY
    GMY<-mean(combined_dat$y)
    ##N
    n_total<-nrow(combined_dat)
    
    ####Getting GW stats
    gw_test<-
      ##Using the "combined" dataset
      combined_dat%>%
      ##Calculations by group
      group_by(ref_foc)%>%
      ##Getting initial summary stats
      summarise(
        #test means
        mean_x=mean(x),
        #criterion means
        mean_y=mean(y),
        #sub-group Ns
        n=sum(dummy_column),
        #test SEs
        ex=sum((x-mean_x)^2),
        #criterion SEs
        ey=sum((y-mean_y)^2),
        #covariance Ses
        exy=sum((x-mean_x)*(y-mean_y)),
        .groups = "drop"
      )%>%
      ##Adding C7
      mutate(
        C7=ey-((exy^2)/ex)
      )%>%
      ##Step 3 metrics
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
        .groups = "drop"
      )
    
    
    
    ####Cleary Test####
    ####Reg 1
    reg1<-lm(y~x,data=combined_dat)
    reg1_summary<-summary(reg1)
    step1_r2<-reg1_summary$r.squared
    step1_beta<-lm.beta(reg1)
    ##output metrics
    step1_r<-sqrt(step1_r2)
    
    
    ####Reg 2
    reg2<-lm(y~x+ref_foc,data=combined_dat)
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
    reg3<-lm(y~x*ref_foc,data=combined_dat)
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
    #car::ncvTest(reg1)
    
    ####Additional effect sizes
    ##Possible Cohen
    cohen_f2<-step3_delta_r2/(1-(step1_r2+step2_delta_r2+step3_delta_r2))
    ##Liu & Yuan
    liu_yuan<-step3_delta_r2/(step1_r2+step3_delta_r2)
    ##Standard Error of Z-test comparing two independent correlations
    se_z_corr<-sqrt( (1/(focal_n-3)) + (1/(reference_n-3)) )
    
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
        liu_yuan=liu_yuan,
        se_z_corr=se_z_corr
      )
    
    
    ####Fisher's Z test
    ##Fisher z-test for ref and focal
    fisher_z_reference<-DescTools::FisherZ(manipulation_check$observed_reference_x_y_cor)
    fisher_z_focal<-DescTools::FisherZ(manipulation_check$observed_focal_x_y_cor)
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
    
    ####Separate Regressions for Each Group####
    ####reference
    ref_reg<-lm(y~x,data=subset(combined_dat,ref_foc=="a_reference"))
    ref_reg_beta<-lm.beta(ref_reg)
    
    ref_constant<-ref_reg_beta$coefficients["(Intercept)"]
    ref_b<-ref_reg_beta$coefficients["x"]
    ref_beta<-ref_reg_beta$standardized.coefficients["x"]
    ref_p<-tidy(ref_reg)%>%filter(term=="x")%>%pull(p.value)
    
    
    ####focal
    focal_reg<-lm(y~x,data=subset(combined_dat,ref_foc=="b_focal"))
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
    ##mean
    levene_test_mean<-
      tidy(
        car::leveneTest(y~as.factor(ref_foc),combined_dat,center=mean)
      )%>%
      filter(!is.na(statistic))%>%
      transmute(
        levene_mean_f=statistic,
        levene_mean_p=p.value
      )
    
    ##median
    levene_test_median<-
      tidy(
        car::leveneTest(y~as.factor(ref_foc),combined_dat,center=median)
      )%>%
      filter(!is.na(statistic))%>%
      transmute(
        levene_median_f=statistic,
        levene_median_p=p.value
      )
    
    
    ####Educational Measurement Differential Validity/Prediction Approach Output (residual metrics)####
    residual_metrics<-
      tibble(
        ref_foc=gsub("a_|b_","",combined_dat$ref_foc),
        residuals=step1_beta$residuals,
        residuals_z=scale(residuals)[1:(n)]
      )%>%
      group_by(ref_foc)%>%
      summarise(
        residuals_m=mean(residuals),
        residuals_sd=sd(residuals),
        residuals_z_m=mean(residuals_z),
        residuals_z_sd=sd(residuals_z),
        .groups = "drop"
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
    
    
    ####getting final data for output####
    final_dat<-
      tibble(
        x_var_name=x_var_name,
        y_var_name=y_var_name,
        condition=condition,
        iteration=iteration,
        #include_in_analysis=include_in_analysis,
        n=n,
        population_percentage=population_percentage,
        focal_n=focal_n,
        reference_n=reference_n,
        reference_sd_y=reference_sd_y,
        focal_sd_y=focal_sd_y,
        reference_validity=reference_validity,
        focal_validity=focal_validity,
        intercept_diffs=intercept_diffs,
        d=D
      )%>%
      bind_cols(
        manipulation_check,
        t_test_results,
        effect_size_calcs,
        gw_test,
        cleary_test,
        homoscedasticity_tests,
        sep_reg_tests,
        residual_metrics,
        levene_test_mean,
        levene_test_median,
        fisher_z_test
      )
    print(z)
    
    
    
    ####Outputting simulated data
    if(do_you_want_to_output_simulated_data=="Yes"){
      
      ##create output location 
      sim_dat_path<-paste0(file_output_folder,"/","sim_dat/")
      
      ##create output folder if needed
      if(!dir.exists(sim_dat_path)){
        dir.create(sim_dat_path)
      }
      
      ##output data
      readr::write_csv(
        combined_dat,
        paste0(
          sim_dat_path,
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".csv"
        )
      )
    }
    
    
    ####Generating scatterplot####
    if(do_you_want_to_output_scatterplots=="Yes"){
    
      ##create output location 
      scatterplot_path<-paste0(file_output_folder,"/","scatterplots/")
      
      ##create output folder if needed
      if(!dir.exists(scatterplot_path)){
        dir.create(scatterplot_path)
        }
    
    ####Get stats
    ##reference
    ref_r2<-round((manipulation_check$observed_reference_x_y_cor)^2,digits=2)
    ref_x<-round(ref_reg$coefficients["x"],digits=2)
    ref_intercept<-round(ref_reg$coefficients["(Intercept)"],digits=2)
    
    ##focal
    focal_r2<-round((manipulation_check$observed_focal_x_y_cor)^2,digits=2)
    focal_x<-round(focal_reg$coefficients["x"],digits=2)
    focal_intercept<-round(focal_reg$coefficients["(Intercept)"],digits=2)
    
    ##data cleaning
    combined_dat1<-combined_dat%>%
      mutate(
        ref_foc=case_when(
          ref_foc=="a_reference"~"Reference",
          ref_foc=="b_focal"~"Focal"
        )
      )%>%
      rename(`Focal/Reference`=ref_foc)
    
    ##open a png
    png(
      filename=paste0(
        scatterplot_path,
        "c",condition,
        "_i_",iteration,
        "_rv_",reference_validity,
        "_fv_",focal_validity,
        "_fsdy_",focal_sd_y,
        "_intdiffs_",intercept_diffs,
        "_D_",D,
        ".png"
        ), 
      width=480,
      height=480,
      pointsize=12,
      bg="white",
      res=NA
      )
    ##make graph
    scatterplot<-
      combined_dat1%>%
      ggplot(aes(x=x,y=y,color=`Focal/Reference`))+
      geom_point(alpha=0.2)+ 
      geom_smooth(method=lm,se=FALSE,fullrange=TRUE)+
      annotate(
        geom="text",
        x=0,
        y=4,
        label=paste0("Ref: R^2 = ",ref_r2,"; y = ",ref_x,"x + ",ref_intercept)
      )+
      annotate(
        geom="text",
        x=0,
        y=3.5,
        label=paste0("Focal: R^2 = ",focal_r2,"; y = ",focal_x,"x + ",focal_intercept)
      )
    print(scatterplot)
    ##close png
    dev.off()
    
    }
    
    
    ####Histograms####
    if(do_you_want_to_output_histograms=="Yes"){
      
      ##create output location 
      histogram_path<-paste0(file_output_folder,"/","histograms/")
      
      ##create output folder if needed
      if(!dir.exists(histogram_path)){
        dir.create(histogram_path)
      }
      
      ##data cleaning
      combined_dat1<-combined_dat%>%
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
          histogram_path,
          "x_hist_",
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".png"
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
          histogram_path,
          "y_hist_",
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".png"
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
    
    
    ####Residual plots####
    if(do_you_want_to_output_residual_plots=="Yes"){
      
      ##create output location 
      residual_plot_path<-paste0(file_output_folder,"/","residual_plots/")
      
      ##create output folder if needed
      if(!dir.exists(residual_plot_path)){
        dir.create(residual_plot_path)
      }
      
      ##Model 1
      png(
        filename=paste0(
          residual_plot_path,
          "model1_resid_",
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".png"
        ), 
        width=480,
        height=480,
        pointsize=12,
        bg="white",
        res=NA
      )
      reg1_resid<-augment(reg1)
      reg1_resid_plot<-
        ggplot(reg1_resid, aes(x, y)) +
        geom_point() +
        stat_smooth(method = lm, se = FALSE) +
        geom_segment(aes(xend = x, yend = .fitted), color = "red", size = 0.3)
      print(reg1_resid_plot)
      dev.off()
    
      ##Model 2
      png(
        filename=paste0(
          residual_plot_path,
          "model2_resid_",
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".png"
        ), 
        width=480,
        height=480,
        pointsize=12,
        bg="white",
        res=NA
      )
      reg2_resid<-augment(reg2)
      reg2_resid_plot<-
        ggplot(reg2_resid, aes(x, y)) +
        geom_point() +
        stat_smooth(method = lm, se = FALSE) +
        geom_segment(aes(xend = x, yend = .fitted), color = "red", size = 0.3)
      print(reg2_resid_plot)
      dev.off()
      
      ##Model 3
      png(
        filename=paste0(
          residual_plot_path,
          "model3_resid_",
          "c",condition,
          "_i_",iteration,
          "_rv_",reference_validity,
          "_fv_",focal_validity,
          "_fsdy_",focal_sd_y,
          "_intdiffs_",intercept_diffs,
          "_D_",D,
          ".png"
        ), 
        width=480,
        height=480,
        pointsize=12,
        bg="white",
        res=NA
      )
      reg3_resid<-augment(reg3)
      reg3_resid_plot<-
        ggplot(reg3_resid, aes(x, y)) +
        geom_point() +
        stat_smooth(method = lm, se = FALSE) +
        geom_segment(aes(xend = x, yend = .fitted), color = "red", size = 0.3)
      print(reg3_resid_plot)
      dev.off()
    }
    
    
    return(final_dat)
    
    

    
    ####end function
  }




####applying overall sim
##get conditions and iterations
sim_conditions<-conditions_iterations$z
#sim_conditions<-1:1000

####Apply function####
if(do_you_want_to_use_parallel_processing=="Yes"){
  
  ##parallelized
  {tic()
    monte_carlo_output<-
      bind_rows(
        mclapply(
          sim_conditions,
          monte_carlo_data_simulation,
          mc.cores = detectCores()-1,
          conditions_iterations,
          date_of_interest,
          output_location,
          file_output_folder
        )
      )
    toc()}
}else if(do_you_want_to_use_parallel_processing=="No"){
  
  ##normal
  {tic()
    monte_carlo_output<-
      bind_rows(
        lapply(
          sim_conditions,
          monte_carlo_data_simulation,
          conditions_iterations,
          date_of_interest,
          output_location,
          file_output_folder
        )
      )
    toc()}
  }

##output results
readr::write_csv(
  monte_carlo_output,
  paste0(
    file_output_folder,
    "/mc_sim_results_",
    date_of_interest,
    ".csv"
  )
)

