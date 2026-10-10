#!/usr/bin/env Rscript
# Nested country-intercept sensitivity, preserving original primary outputs.
suppressPackageStartupMessages({library(data.table); library(splines); library(sandwich)})
args0 <- grep('^--file=', commandArgs(FALSE), value=TRUE)
root <- dirname(normalizePath(sub('^--file=', '', args0), winslash='/'))
out <- file.path(root, 'outputs', 'm5_country_intercept')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sc <- fread(file.path(root, 'outputs', 'exposure_scaling.csv'))
heat <- sc$heat_variable
covs <- c('mean_value','SEV_agestd_air_hap','SEV_agestd_air_pm',
 'SEV_agestd_diet_fruit','SEV_agestd_diet_veg','SEV_agestd_metab_bmi','mean_LDL')
cols <- c('measure_id','metric_id','age_id','sex_id','location_id','year','val',
 'lower','upper','population','superregion','region_wb',covs,heat)
input <- Sys.getenv('DID_INPUT_FILE', file.path(root,'data','input_panel.csv'))
d <- fread(input, select=cols, showProgress=FALSE)
d <- d[metric_id==1L & age_id==22L & sex_id==3L & measure_id %in% c(1L,2L,5L,6L)]
stopifnot(d[,anyDuplicated(paste(location_id,year,measure_id))]==0L)
d <- d[is.finite(val)&val>=0&is.finite(lower)&is.finite(upper)&lower>=0&lower<=val&upper>=val&
 is.finite(population)&population>0&!is.na(superregion)&trimws(superregion)!=''&
 !is.na(region_wb)&trimws(region_wb)!='']
d <- d[Reduce(`&`,lapply(d[,c(covs,heat),with=FALSE],is.finite))]
stopifnot(nrow(d)==24552L, uniqueN(d$location_id)==198L)
adj <- c(M3a='factor(region_wb)',M3b='factor(superregion)',
 M5='ns(year, df=3):factor(superregion)',
 'M5*'='ns(year, df=3):factor(superregion) + factor(location_id)')
rhs <- paste(paste(covs,collapse=' + '),adj,sep=' + ')
names(rhs) <- names(adj)
definitions <- data.table(model_id=names(rhs),formula=paste('val ~ heat_iqr +',rhs,'+ offset(log(population))'))
fwrite(definitions,file.path(out,'model_definitions.csv'),bom=TRUE)
res <- list(); diag <- list(); variation <- list(); k <- 0L
for(m in c(1L,2L,5L,6L)) {
 dm <- copy(d[measure_id==m])
 for(v in heat) {
  s <- sc[heat_variable==v]
  dm[,heat_iqr:=(get(v)-s$center_median)/s$iqr]
  within <- dm[,get(v)-mean(get(v)),by=location_id]$V1
  for(id in names(rhs)) {
   k <- k+1L
   f <- as.formula(definitions[model_id==id,formula])
   fit <- glm(f,data=dm,family=quasipoisson('log'),control=glm.control(maxit=100,epsilon=1e-8))
   stopifnot(isTRUE(fit$converged),nobs(fit)==6138L)
   V <- vcovCL(fit,cluster=dm$location_id,type='HC1',cadjust=TRUE)
   b <- unname(coef(fit)['heat_iqr']); se <- unname(sqrt(diag(V))['heat_iqr'])
   stopifnot(is.finite(b),is.finite(se),se>0)
   crit <- qt(.975,197)
   res[[k]] <- data.table(measure_id=m,heat_variable=v,heat_label=s$heat_label,model_id=id,
    beta=b,std_error_clustered=se,ci_low=b-crit*se,ci_high=b+crit*se,
    p_value=2*pt(abs(b/se),197,lower.tail=FALSE),rate_ratio=exp(b),
    n_observations=nobs(fit),n_countries=198L,increment=s$iqr)
   diag[[k]] <- data.table(measure_id=m,heat_variable=v,model_id=id,converged=fit$converged,
    rank=fit$rank,design_columns=ncol(model.matrix(fit)),iterations=fit$iter,dispersion=summary(fit)$dispersion)
   # Unweighted projection is a design diagnostic, not GLM Fisher information.
   projection <- lm(as.formula(paste(v,'~',rhs[[id]])),data=dm)
   variation[[k]] <- data.table(measure_id=m,heat_variable=v,model_id=id,
    sd_heat_raw=sd(dm[[v]]),sd_heat_within_country=sd(within),
    sd_heat_residual=sd(residuals(projection)),sd_heat_residual_iqr=sd(residuals(projection))/s$iqr)
  }
 }
 message('Completed outcome ',m)
}
r <- rbindlist(res); r[,p_adjust_bh:=p.adjust(p_value,'BH'),by=model_id]
vr <- rbindlist(variation)
r <- merge(r,vr,by=c('measure_id','heat_variable','model_id'))
ref <- fread(file.path(root,'outputs','did_results_full.csv'))[model_id %in% c('M3a','M3b','M5')]
qc <- merge(r[model_id!='M5*'],ref,by=c('measure_id','heat_variable','model_id'),suffixes=c('_new','_stored'))
qc <- qc[,.(measure_id,heat_variable,model_id,beta_difference=beta_new-beta_stored,se_difference=std_error_clustered_new-std_error_clustered_stored)]
stopifnot(nrow(qc)==96L,max(abs(qc$beta_difference))<1e-8,max(abs(qc$se_difference))<1e-8)
comparisons <- rbindlist(lapply(c('M5','M3a','M3b'),function(base){
 x <- merge(r[model_id=='M5*'],r[model_id==base],by=c('measure_id','heat_variable'),suffixes=c('_star','_base'))
 x[,baseline:=base]
 x[,`:=`(delta_beta=beta_star-beta_base,delta_se=std_error_clustered_star-std_error_clustered_base,
  se_ratio=std_error_clustered_star/std_error_clustered_base,
  absolute_beta_change_pct=100*(abs(beta_star)-abs(beta_base))/abs(beta_base),
  direction_reversal=sign(beta_star)!=sign(beta_base),
  residual_sd_change_pct=100*(sd_heat_residual_star/sd_heat_residual_base-1))]
 x
}))
fwrite(r,file.path(out,'model_results.csv'),bom=TRUE)
fwrite(comparisons,file.path(out,'model_comparisons.csv'),bom=TRUE)
fwrite(vr,file.path(out,'exposure_variation.csv'),bom=TRUE)
fwrite(rbindlist(diag),file.path(out,'model_diagnostics.csv'),bom=TRUE)
fwrite(qc,file.path(out,'original_model_reproduction_qc.csv'),bom=TRUE)
fmt <- function(x) formatC(x,digits=5,format='f')
tab <- function(x) c(paste0('| ',paste(names(x),collapse=' | '),' |'),
 paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
 apply(x,1,function(z) paste0('| ',paste(z,collapse=' | '),' |')))
labels <- c('1'='死亡','2'='DALYs','5'='患病','6'='发病')
report <- c('# M5* 国家截距敏感性分析', '',
 'M5* = 原 M5 + factor(location_id)。保留7项协变量、超级区域特异3 df自然样条时间趋势、准Poisson/log及log(population) offset；国家截距是固定效应。',
 '每模型每结局均为1990–2020年、198国、6,138国家年。使用相同原始暴露IQR、国家聚类HC1标准误及197自由度t推断。M3a为世行区域截距，M3b为GBD超级区域截距。', '',
 'Δβ = β(M5*) − β(基准)，ΔSE = SE(M5*) − SE(基准)。|β|变化百分比用于区分关联减弱与带符号系数下降；方向反转单独标记。基准β接近零时百分比不稳定。',
 '原始SD(Heat)在同一样本中不随模型改变。国家内SD为国家均值去除后的样本SD；残差SD为Heat对该模型全部协变量及固定效应/趋势进行不加权OLS投影后的残差样本SD。该指标衡量设计中的剩余暴露变异，不等同于准Poisson模型的加权信息量。人口offset不是待估调整项，不进入此投影。', '',
 '## 比较表')
for(base in c('M5','M3a','M3b')) {
 x <- comparisons[baseline==base]
 t <- x[,.(`结局`=unname(labels[as.character(measure_id)]),`指标`=heat_label_star,
  `基准β`=fmt(beta_base),`M5* β`=fmt(beta_star),`Δβ`=fmt(delta_beta),
  `基准SE`=fmt(std_error_clustered_base),`M5* SE`=fmt(std_error_clustered_star),`ΔSE`=fmt(delta_se),
  `SE比`=fmt(se_ratio),`|β|变化%`=fmt(absolute_beta_change_pct),`反转`=ifelse(direction_reversal,'是','否'),
  `原始SD`=fmt(sd_heat_raw_star),`基准残差SD`=fmt(sd_heat_residual_base),`M5*残差SD`=fmt(sd_heat_residual_star))]
 report <- c(report,'',paste('### M5* 与',base),'',tab(t))
}
report <- c(report,'','## 解释边界','',
 '优先比较M5*与M5，因为两者只差国家截距。与M3a/M3b比较还改变了时间趋势结构，不能把差异全部归于国家截距。',
 'β变化小、SE增加且残差暴露SD下降，与信息损失/估计精度降低相容；这不是机制的证明。|β|明显减弱而SE变化有限，与国家稳定差异解释部分原关联相容，应谨慎解释。若β和SE均明显改变或方向反转，则属于混合结果，不应强行二选一。',
 '未预设“小变化/大幅下降”的阈值，不将连续差异事后转成确定性诊断。标准误还受模型均值、离散程度、参数数目及聚类校正影响。此分析不能单独确证混杂消除或因果效应。',
 '',sprintf('原模型复现：96项M3a/M3b/M5结果，最大|Δβ|=%.3g，最大|ΔSE|=%.3g。',max(abs(qc$beta_difference)),max(abs(qc$se_difference))))
summary_text <- c('', '## 本次结果解释', '',
 sprintf('M5*相对M5：32项估计中%d项|β|减小，%d项方向反转，%d项SE增加。',
  comparisons[baseline=='M5',sum(absolute_beta_change_pct<0)],
  comparisons[baseline=='M5',sum(direction_reversal)],
  comparisons[baseline=='M5',sum(se_ratio>1)]),
 '核心复合昼夜热指标：死亡、DALYs、患病、发病的|β|分别减少96.0%、97.1%、80.3%、93.4%；SE分别减少82.5%、83.5%、83.0%、84.0%。残差SD从21.39675降至19.01813（下降11.1%），原始SD均为28.32898。',
 '这不是“β基本不变且SE增加”的信息损失模式。核心指标的β和SE同时大幅下降，也不是严格的“SE变化小”模式；结果更符合国家稳定差异可能解释部分原关联，但不能据此确定差异来源或证明因果混杂。SE下降提示国家截距吸收结局异质性后精度也发生了改变，应同时报告全部指标的混合结果。')
writeLines(c(report,summary_text),file.path(out,'M5star_comparison_report.md'),useBytes=TRUE)
print(comparisons[baseline=='M5' & heat_variable=='compound_day_night_cumulative_heat',
 .(measure_id,beta_base,beta_star,delta_beta,std_error_clustered_base,std_error_clustered_star,se_ratio,absolute_beta_change_pct,sd_heat_residual_base,sd_heat_residual_star)])

