#!/usr/bin/env Rscript
# Population-weighted exposure sensitivity; the original M5 is not overwritten.
suppressPackageStartupMessages({library(data.table);library(splines);library(sandwich)})
options(encoding='UTF-8')
args0=grep('^--file=',commandArgs(FALSE),value=TRUE)
root=dirname(normalizePath(sub('^--file=','',args0),winslash='/'))
out=file.path(root,'outputs','population_exposure');dir.create(out,showWarnings=FALSE,recursive=TRUE)
heat=c('absolute_cumulative_heat','relative_cumulative_heat','mixed_cumulative_heat',
 'daytime_cumulative_heat','nighttime_cumulative_heat','compound_day_night_cumulative_heat',
 'annual_mean_humidex','annual_mean_tmean')
covs=c('mean_value','SEV_agestd_air_hap','SEV_agestd_air_pm','SEV_agestd_diet_fruit',
 'SEV_agestd_diet_veg','SEV_agestd_metab_bmi','mean_LDL')
cols=c('measure_id','metric_id','age_id','sex_id','location_id','location','year','val',
 'lower','upper','population','superregion','region_wb','fid',covs,heat)
panel=fread(file.path(root,'data','input_panel.csv'),select=cols,encoding='UTF-8',showProgress=FALSE)
dat=panel[metric_id==1L&age_id==22L&sex_id==3L&measure_id%in%c(1L,2L,5L,6L)]
stopifnot(dat[,anyDuplicated(paste(location_id,year,measure_id))]==0L)
dat=dat[is.finite(val)&val>=0&is.finite(lower)&is.finite(upper)&lower>=0&lower<=val&upper>=val&
 is.finite(population)&population>0&!is.na(superregion)&trimws(superregion)!=''&
 !is.na(region_wb)&trimws(region_wb)!='']
valid_rows=dat[,Reduce(`&`,lapply(.SD,is.finite)),.SDcols=c(covs,heat)]
dat=dat[valid_rows]
stopifnot(nrow(dat)==24552L,uniqueN(dat$location_id)==198L)
cat('Main eligible sample:', nrow(dat), 'rows;', uniqueN(dat$location_id), 'countries\n')
sc=fread(file.path(root,'outputs','exposure_scaling.csv'))[heat_variable=='compound_day_night_cumulative_heat']
d=dat[measure_id==1L]
d[,heat_iqr:=(compound_day_night_cumulative_heat-sc$center_median)/sc$iqr]
f=as.formula(paste('val ~ heat_iqr +',paste(covs,collapse=' + '),'+ ns(year,df=3):factor(superregion) + offset(log(population))'))
m=glm(f,family=quasipoisson(link='log'),data=d,control=glm.control(maxit=100,epsilon=1e-8))
r=fread(file.path(root,'outputs','did_results_full.csv'))[model_id=='M5'&measure_id==1L&heat_variable=='compound_day_night_cumulative_heat']
stopifnot(abs(coef(m)['heat_iqr']-r$beta)<1e-8)
cat('Original mortality CDNCH coefficient reproduced:',coef(m)['heat_iqr'],'\n')
