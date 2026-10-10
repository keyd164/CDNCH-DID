project <- 'D:/DID'
audit <- file.path(project, 'editorial/numeric_validation')
dir.create(audit, recursive=TRUE, showWarnings=FALSE)
script <- readLines(file.path(project,'01_run_did.R'), encoding='UTF-8')
# Reuse the original sample selection and scaling without executing output writes.
init <- paste(script[1:205], collapse='\n')
init <- sub('project_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE))', 'project_dir <- "D:/DID"', init, fixed=TRUE)
eval(parse(text=init))
expressions <- parse(text=script)
needed <- c('cluster_vcov','scale_heat','model_rhs','simulate_attributable_burden_ci','spline_basis_contrast','fit_compound_spline','standard_nonlinear_rr_definition','standard_nonlinear_rr_source','attributable_uncertainty_draws','attributable_uncertainty_seed_base')
for (e in expressions) if (is.call(e) && identical(e[[1]],as.name('<-')) && is.symbol(e[[2]]) && as.character(e[[2]]) %in% needed) eval(e)
checks <- list(); fitted <- list()
for (id in expected_measures) {
 message('Re-fitting outcome ',id)
 dat <- analysis_data[measure_id==id]
 fit <- fit_compound_spline(dat,'M5',NA_real_)
 fitted[[as.character(id)]] <- fit
 scaled <- scale_heat(dat,'compound_day_night_cumulative_heat')
 basis <- ns(scaled$heat_iqr,df=4); bn <- paste0('heat_ns',1:4)
 scaled[,(bn):=as.data.table(basis)]
 fm <- as.formula(paste0('val ~ ',model_rhs(paste(bn,collapse=' + '),'M5'),' + offset(log(population))'))
 mod <- glm(fm,family=quasipoisson(link='log'),data=scaled,control=glm.control(maxit=100L,epsilon=1e-8))
 mm <- model.matrix(mod)
 # Independent full-model prediction contrast at identical covariate values.
 ref <- copy(scaled); ref[,(bn):=as.data.table(matrix(rep(as.numeric(predict(basis,0)),each=nrow(ref)),ncol=4))]
 independent_rr <- exp(predict(mod,newdata=scaled,type='link')-predict(mod,newdata=ref,type='link'))
 independent_af <- pmax(0,(independent_rr-1)/independent_rr)
 checks[[as.character(id)]] <- data.table(measure_id=id, n=nrow(dat),countries=uniqueN(dat$location_id),regions=uniqueN(dat$superregion),median=median(dat$compound_day_night_cumulative_heat),model_columns=ncol(mm),rank=mod$rank,converged=mod$converged,independent_rr_max_error=max(abs(independent_rr-fit$attributable$relative_risk)),independent_burden=sum(dat$val*independent_af),rr_le_one_count=sum(independent_rr<=1),formula=paste(deparse(fm),collapse=' '),separate_region_intercept_columns=sum(grepl('^factor\\(superregion\\)',colnames(mm))))
 if (id==1) writeLines(colnames(mm),file.path(audit,'design_matrix_columns.txt'))
}
new <- rbindlist(lapply(fitted,`[[`,'attributable_model_uncertainty'))
old <- fread(file.path(project,'outputs/attributable_burden_summary_m5.csv'))
cols <- c('cumulative_attributable_burden','cumulative_attributable_burden_model_ci_low','cumulative_attributable_burden_model_ci_high','cumulative_attributable_fraction','cumulative_attributable_fraction_model_ci_low','cumulative_attributable_fraction_model_ci_high')
comparison <- merge(new[,c('measure_id',cols),with=FALSE],old[,c('measure_id',cols),with=FALSE],by='measure_id',suffixes=c('_rerun','_existing'))
for (v in cols) comparison[,(paste0(v,'_difference')):=get(paste0(v,'_rerun'))-get(paste0(v,'_existing'))]
fwrite(new,file.path(audit,'recomputed_summary.csv'));fwrite(comparison,file.path(audit,'existing_output_comparison.csv'));fwrite(rbindlist(checks),file.path(audit,'model_and_prediction_checks.csv'))
rows <- rbindlist(lapply(fitted,`[[`,'attributable'))
fwrite(rows,file.path(audit,'recomputed_country_year.csv'))
fwrite(rbindlist(lapply(fitted,`[[`,'curve')),file.path(audit,'recomputed_curves.csv'))
writeLines(capture.output(sessionInfo()),file.path(audit,'R_session.txt'))
print(new[,c('measure_id',cols),with=FALSE]);print(rbindlist(checks))
