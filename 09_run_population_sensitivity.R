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
ex=fread(file.path(out,'country_year_exposures.csv'),encoding='UTF-8')
stopifnot(ex[,anyDuplicated(paste(fid,year))]==0L)
joined=merge(dat,ex,by=c('fid','year'),all.x=TRUE,sort=FALSE)
stopifnot(nrow(joined)==nrow(dat))
popvars=paste0(heat,'_population_weighted')
keys=unique(joined[year>=2000,c('location_id','year','sampling_method','population_weight_sum',popvars),with=FALSE])
keys[,valid_population_exposure:=Reduce(`&`,lapply(.SD,is.finite))&population_weight_sum>0,.SDcols=popvars]
keys[is.na(valid_population_exposure),valid_population_exposure:=FALSE]
fwrite(keys[,.(location_id,year,sampling_method,population_weight_sum,valid_population_exposure)],file.path(out,'sample_eligibility.csv'),bom=TRUE)
valid=keys[valid_population_exposure==TRUE,.(location_id,year)]
matched=merge(joined[year>=2000],valid,by=c('location_id','year'))
stopifnot(matched[,.N,by=.(location_id,year)][,all(N==4L)])
scale=fread(file.path(root,'outputs','exposure_scaling.csv'),encoding='UTF-8')
ref=fread(file.path(root,'outputs','did_results_full.csv'),encoding='UTF-8')[model_id=='M5']
labels=c('1'='Mortality','2'='DALYs','5'='Prevalence','6'='Incidence')
specifications=c('main_area_1990_2020','matched_area_2000_2020','population_2000_2020')
fitted=list();diagnostics=list();idx=0L
for (s in specifications) {
 for (m in c(1L,2L,5L,6L)) {
  d=if(s=='main_area_1990_2020') copy(joined[measure_id==m]) else copy(matched[measure_id==m])
  for(v in heat) {
   value=if(s=='population_2000_2020') paste0(v,'_population_weighted') else v
   sc=scale[heat_variable==v]
   d[,heat_iqr:=(get(value)-sc$center_median)/sc$iqr]
   formula=as.formula(paste('val ~ heat_iqr +',paste(covs,collapse=' + '),
    '+ ns(year, df=3):factor(superregion) + offset(log(population))'))
   fit=glm(formula,family=quasipoisson(link='log'),data=d,
           control=glm.control(maxit=100,epsilon=1e-8))
   stopifnot(isTRUE(fit$converged))
   V=vcovCL(fit,cluster=d$location_id,type='HC1',cadjust=TRUE)
   b=unname(coef(fit)['heat_iqr']);se=unname(sqrt(diag(V))['heat_iqr'])
   nc=uniqueN(d$location_id);crit=qt(.975,nc-1L)
   idx=idx+1L
   fitted[[idx]]=data.table(specification=s,measure_id=m,outcome=unname(labels[as.character(m)]),
    heat_variable=v,increment=sc$iqr,increment_unit=sc$unit,scaling='Full main-sample area-weighted IQR',
    beta=b,std_error_clustered=se,ci_low=b-crit*se,ci_high=b+crit*se,
    rate_ratio=exp(b),rr_ci_low=exp(b-crit*se),rr_ci_high=exp(b+crit*se),
    p_value=2*pt(abs(b/se),nc-1L,lower.tail=FALSE),n_observations=nobs(fit),
    n_countries=nc,first_year=min(d$year),last_year=max(d$year),
    own_exposure_iqr=IQR(d[[value]]),n_nearest_country_years=sum(d$sampling_method=='representative_point_nearest'),
    formula=deparse(formula,width.cutoff=500))
   diagnostics[[idx]]=data.table(specification=s,measure_id=m,heat_variable=v,
    converged=fit$converged,rank=fit$rank,columns=ncol(model.matrix(fit)),
    dispersion=summary(fit)$dispersion,iterations=fit$iter,
    n_observations=nobs(fit),n_countries=nc,n_regions=uniqueN(d$superregion))
  }
 }
 message('Completed ',s)
}
result=rbindlist(fitted)
result[,p_adjust_bh:=p.adjust(p_value,method='BH'),by=specification]
result[,`:=`(beta_per_own_iqr=beta*own_exposure_iqr/increment,
 rr_per_own_iqr=exp(beta*own_exposure_iqr/increment),
 rr_low_per_own_iqr=exp(ci_low*own_exposure_iqr/increment),
 rr_high_per_own_iqr=exp(ci_high*own_exposure_iqr/increment))]
check=merge(result[specification=='main_area_1990_2020'],ref,
 by=c('heat_variable','measure_id'),suffixes=c('_recomputed','_stored'))
stopifnot(nrow(check)==32L,max(abs(check$beta_recomputed-check$beta_stored))<1e-8,
 max(abs(check$std_error_clustered_recomputed-check$std_error_clustered_stored))<1e-8)
fwrite(check[,.(heat_variable,measure_id,beta_difference=beta_recomputed-beta_stored,
 se_difference=std_error_clustered_recomputed-std_error_clustered_stored)],file.path(out,'main_model_reproduction_qc.csv'),bom=TRUE)
fwrite(result,file.path(out,'sensitivity_results.csv'),bom=TRUE)
fwrite(rbindlist(diagnostics),file.path(out,'model_diagnostics.csv'),bom=TRUE)
area=result[specification=='matched_area_2000_2020']
pop=result[specification=='population_2000_2020']
comparison=merge(area,pop,by=c('heat_variable','measure_id','outcome'),suffixes=c('_area','_population'))
comparison[,`:=`(beta_difference=beta_population-beta_area,
 percent_change_absolute_beta=100*(abs(beta_population)-abs(beta_area))/abs(beta_area),
 direction_reversal=sign(beta_population)!=sign(beta_area))]
fwrite(comparison,file.path(out,'matched_comparison.csv'),bom=TRUE)
exposure_stats=rbindlist(lapply(heat,function(v){
 p=unique(matched[,c('location_id','year',v,paste0(v,'_population_weighted')),with=FALSE]);
 a=p[[v]];w=p[[paste0(v,'_population_weighted')]]
 data.table(heat_variable=v,n_country_years=nrow(p),n_countries=uniqueN(p$location_id),
 area_median=median(a),population_median=median(w),area_iqr=IQR(a),population_iqr=IQR(w),
 pearson=cor(a,w),spearman=cor(a,w,method='spearman'),mean_difference=mean(w-a),
 median_difference=median(w-a),rmse=sqrt(mean((w-a)^2)))
}))
fwrite(exposure_stats,file.path(out,'exposure_comparison.csv'),bom=TRUE)
cat('Countries:',uniqueN(matched$location_id),'country-years:',nrow(matched)/4,'\n')
print(comparison[heat_variable=='compound_day_night_cumulative_heat',
 .(outcome,rate_ratio_area,rr_ci_low_area,rr_ci_high_area,rate_ratio_population,
 rr_ci_low_population,rr_ci_high_population,p_value_area,p_value_population)])
