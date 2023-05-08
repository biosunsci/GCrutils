# library("survival")
# library("survminer")

get_threshold = function(data,z_score=3){
    mean_d = mean(data)
    std_d = sd(data)
    list(mean_d - z_score*std_d, mean_d + z_score*std_d)
}
#' Title
#'
#' @param os_table
#' @param cover
#'
#' @return
#' @export
#'
#' @examples
os_left_align = function(os_table,cover=12){
    if (is.numeric(os_table$survivalEvent)){
        mask1 = os_table$survivalEvent == 1
    }else if(is.logical(os_table$survivalEvent)){
        mask1 = os_table$survivalEvent == TRUE
    }else
        stop('wrong dtype of survivalEvent, must be int or bool')
    os_table = os_table[(mask1)|(os_table$survivalMonth >= cover),]
    os_table
}
#' Title
#'
#' @param os_table
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
os_right_align = function(os_table,cutoff=NULL){
    if (is.null(cutoff)){
        cutoff = get_threshold(os_table.survivalMonth)[[2]]
    }
    cutoff = round(cutoff)
    mask = os_table$survivalMonth > cutoff
    os_table[mask,'survivalMonth'] = cutoff
    os_table[mask,'survivalEvent'] = 0
    os_table$survivalEvent = as.integer(os_table$survivalEvent)
    os_table
}
#' Title
#'
#' @param os_table
#' @param left
#' @param right
#'
#' @return
#' @export
#'
#' @examples
os_align = function(os_table,left=12,right=60){
    #     os_table = os_table[['patient_id','survivalEvent','survivalMonth']]
    if (os_table$survivalEvent %>% is.logical)
        os_table$survivalEvent = os_table$survivalEvent %>% as.numeric
    os_table = os_left_align(os_table,cover=left)
    os_table = os_right_align(os_table,cutoff=right)
    os_table
}


#' Title
#'
#' @param data
#' @param covariates
#' @param key
#' @param test.method  test.method in c('logrank','wald','likelihood'), deault wald
#' @param tags
#' @param auto_anno
#' @param join
#'
#' @return
#' @export
#'
#' @examples
library(ggfortify)
yplot_survival = function(data, covariates
                          ,key=NULL
                          ,test.method='wald'
                          ,tags=c('month','status')
                          ,auto_anno=TRUE
                          ,join='_'
){

    if (! is.null(key)){
        templ = stringr::str_flatten(c('survival::Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'survival::Surv(survivalMonth,survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(x) {
                                res = list()
                                formu=as.formula(paste(templ, x))
                                res[['formula']] = formu
                                res[['feature']] = x
                                res[['model']] = survival::survfit(formu, data = data,)
                                diff = survival::coxph(formu, data = data) %>% summary
                                if (test.method=='logrank'){
                                    p.value = diff$sctest[['pvalue']]
                                }else if (test.method=='wald'){
                                    p.value = diff$waldtest[['pvalue']]
                                }else if (test.method=='likelihood'){
                                    p.value = diff$logtest[['pvalue']]
                                }
                                res[['pval']] = signif(p.value,4)
                                res[['test.method']] = test.method
                                res
                            },simplify = FALSE)
    plots = lapply(univ_formulas, function(li){
        model = li$model
        # fn = feature name
        fn = li$feature
        pval = li$pval
        if (pval < 0.05){
            pval.color = 'red'
        }else{
            pval.color = 'black'
        }
        # cal annotation y depend on group number of feature fn
        grps = data[[fn]] %>% unique
        n = length(grps)

        gg = autoplot(model) +
            labs(color=fn,fill=fn) +
            scale_y_continuous(limits = c(0,1),labels = scales::percent)
        gg = gg +  annotate(
            geom = "text", x=0,y = 0.042*n,
            color=pval.color,
            label = paste('p =',pval,'\n'),
            hjust = "inward"
        )
        if (auto_anno==TRUE){
            df = univ_formulas[[fn]]$model %>% summary %>% .$table
            text_label = paste('[',df[,'records'],']',df %>% row.names,' median survival: ',df[,'median'])
            gg = gg + annotate(
                geom = "text", x=0,y= 0.02*n,
                label = text_label %>% stringr::str_flatten('\n'),
                hjust = "inward"
                #                     ,vjust = 'inward'
            )
        }

        gg
    })
    names(plots) = names(univ_formulas)
    list(data=univ_formulas,plots=plots)
}


get_threshold = function(data,z_score=3){
    mean_d = mean(data)
    std_d = sd(data)
    list(mean_d - z_score*std_d, mean_d + z_score*std_d)
}
#' Title
#'
#' @param os_table
#' @param cover
#'
#' @return
#' @export
#'
#' @examples
os_left_align = function(os_table,cover=12){
    if (is.numeric(os_table$survivalEvent)){
        mask1 = os_table$survivalEvent == 1
    }else if(is.logical(os_table$survivalEvent)){
        mask1 = os_table$survivalEvent == TRUE
    }else
        stop('wrong dtype of survivalEvent, must be int or bool')
    os_table = os_table[(mask1)|(os_table$survivalMonth >= cover),]
    os_table
}
#' Title
#'
#' @param os_table
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
os_right_align = function(os_table,cutoff=NULL){
    if (is.null(cutoff)){
        cutoff = get_threshold(os_table.survivalMonth)[[2]]
    }
    cutoff = round(cutoff)
    mask = os_table$survivalMonth > cutoff
    os_table[mask,'survivalMonth'] = cutoff
    os_table[mask,'survivalEvent'] = 0
    os_table$survivalEvent = as.integer(os_table$survivalEvent)
    os_table
}
#' Title
#'
#' @param os_table
#' @param left
#' @param right
#'
#' @return
#' @export
#'
#' @examples
os_align = function(os_table,left=12,right=60){
    #     os_table = os_table[['patient_id','survivalEvent','survivalMonth']]
    if (os_table$survivalEvent %>% is.logical)
        os_table$survivalEvent = os_table$survivalEvent %>% as.numeric
    os_table = os_left_align(os_table,cover=left)
    os_table = os_right_align(os_table,cutoff=right)
    os_table
}


#' Title
#'
#' @param data
#' @param covariates
#' @param key
#' @param test.method  test.method in c('logrank','wald','likelihood'), deault wald
#' @param tags
#' @param auto_anno
#' @param join
#'
#' @return
#' @export
#'
#' @examples
library(ggfortify)
yplot_survival = function(data, covariates
                          ,key=NULL
                          ,test.method='wald'
                          ,tags=c('month','status')
                          ,auto_anno=TRUE
                          ,join='_'
){

    if (! is.null(key)){
        templ = stringr::str_flatten(c('survival::Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'survival::Surv(survivalMonth,survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(x) {
                                res = list()
                                formu=as.formula(paste(templ, x))
                                res[['formula']] = formu
                                res[['feature']] = x
                                res[['model']] = survival::survfit(formu, data = data,)
                                diff = survival::coxph(formu, data = data) %>% summary
                                if (test.method=='logrank'){
                                    p.value = diff$sctest[['pvalue']]
                                }else if (test.method=='wald'){
                                    p.value = diff$waldtest[['pvalue']]
                                }else if (test.method=='likelihood'){
                                    p.value = diff$logtest[['pvalue']]
                                }
                                res[['pval']] = signif(p.value,4)
                                res[['test.method']] = test.method
                                res
                            },simplify = FALSE)
    plots = lapply(univ_formulas, function(li){
        model = li$model
        # fn = feature name
        fn = li$feature
        pval = li$pval
        if (pval < 0.05){
            pval.color = 'red'
        }else{
            pval.color = 'black'
        }
        # cal annotation y depend on group number of feature fn
        grps = data[[fn]] %>% unique
        n = length(grps)

        gg = autoplot(model) +
            labs(color=fn,fill=fn) +
            scale_y_continuous(limits = c(0,1),labels = scales::percent)
        gg = gg +  annotate(
            geom = "text", x=0,y = 0.042*n,
            color=pval.color,
            label = paste('p =',pval,'\n'),
            hjust = "inward"
        )
        if (auto_anno==TRUE){
            df = univ_formulas[[fn]]$model %>% summary %>% .$table
            text_label = paste('[',df[,'records'],']',df %>% row.names,' median survival: ',df[,'median'])
            gg = gg + annotate(
                geom = "text", x=0,y= 0.02*n,
                label = text_label %>% stringr::str_flatten('\n'),
                hjust = "inward"
                #                     ,vjust = 'inward'
            )
        }

        gg
    })
    names(plots) = names(univ_formulas)
    list(data=univ_formulas,plots=plots)
}



ydo_single_factor_coxph = function(data, covariates, key=NULL, join='_', tags=c('month','status')){
    "MAKE SURE that data values consist of int, which 0 mean no mut, others mean mut
 except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)
 mut_cases = column$values != 0
 wt_cases = total - mut_cases
"

    if (! is.null(key)){
        templ = stringr::str_flatten(c('Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'Surv(survivalMonth,survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(x) {
                                res = list()
                                formu=as.formula(paste(templ, x))
                                res[['formula']] = formu
                                res[['feature']] = x
                                res[['model']] = coxph(formu, data = data,)
                                res
                            },simplify = FALSE)

    univ_results <- lapply(univ_formulas,
                           function(li){
                               feature = li$feature
                               col = data %>% dplyr::pull(!!feature)
                               mut_cases = sum(col!=0,na.rm = TRUE)
                               wt_cases = sum(col==0,na.rm = TRUE)
                               all_cases = length(col)
                               if (all_cases > mut_cases + wt_cases){
                                   # col contains NA
                                   note = paste('contains',all_cases - mut_cases - wt_cases,'NA values')
                               }else{
                                   note = ""
                               }
                               ratio = signif(mut_cases/(mut_cases+wt_cases),digits = 3) * 100
                               x = li$model
                               x <- summary(x)
                               #获取p值
                               p.value<-signif(x$wald["pvalue"], digits=4)
                               p.sig = ymark_psig(p.value)
                               #获取HR
                               HR <-signif(x$coef[2], digits=2);
                               #获取95%置信区间
                               HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                               HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                               HR_95_CI_lower = HR.confint.lower
                               HR_95_CI_upper = HR.confint.upper
                               #feature = r$conf.int %>% row.names

                               res<-c(mut_cases, wt_cases, ratio, p.value
                                      , p.sig
                                      , HR, HR_95_CI_lower, HR_95_CI_upper,note)
                               names(res)<-c('mut_cases','wt_cases','mut_ratio%', "p_value"
                                             , 'p_sig'
                                             , "HR", 'HR_95_CI_lower', 'HR_95_CI_upper','note')
                               return(res)
                           })
    #转换成数据框，并转置
    res = t(as.data.frame(univ_results, check.names = FALSE))
    as.data.frame(res) %>% arrange(p_value)
}


ydo_single_factor_logrank = function(data, covariates, key=NULL, join='_', tags=c('month','status'), dev=FALSE){
    "MAKE SURE that data values consist of int, which 0 mean no mut, others mean mut
 except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)
 mut_cases = column$values != 0
 wt_cases = total - mut_cases
"
    if (! is.null(key)){
        templ = stringr::str_flatten(c('Surv(',key,join,tags[[1]],',',key,join,tags[[2]],')~'))
    }else{
        templ = 'Surv(survivalMonth, survivalEvent)~'
    }
    all_cases = data %>% nrow
    univ_formulas <- sapply(covariates,
                            function(feature) {
                                res = list()
                                formu=as.formula(paste(templ, feature))
                                res[['formula']] = formu
                                res[['feature']] = feature
                                res
                            },simplify = FALSE)
    if (dev ==TRUE){
        return(univ_formulas)
    }
    univ_results <- lapply(univ_formulas,
                           function(li){
                               feature = li$feature
                               col = data %>% dplyr::pull(!!feature)
                               mut_cases = sum(col!=0,na.rm = TRUE)
                               wt_cases = sum(col==0,na.rm = TRUE)
                               all_cases = length(col)
                               if (all_cases > mut_cases + wt_cases){
                                   # col contains NA
                                   note = paste('contains',all_cases - mut_cases - wt_cases,'NA values')
                               }else{
                                   note = ""
                               }
                               if (mut_cases==0 || wt_cases==0){
                                   res<-c(mut_cases, wt_cases, NA, NA
                                          , NA
                                          ,'A(1)/B(0)', NA,NA,NA)
                               }else{
                                   ratio = signif(mut_cases/(mut_cases+wt_cases),digits = 3) * 100
                                   model = survdiff(li[['formula']], data = data)
                                   # model is like
                                   #     survdiff(formula = formu, data = data)
                                   #
                                   #                 N Observed Expected (O-E)^2/E (O-E)^2/V
                                   # Radiotherapy=0 11        4    0.798      12.8      16.1
                                   # Radiotherapy=1 36        0    3.202       3.2      16.1
                                   #
                                   # Chisq= 16.1  on 1 degrees of freedom, p= 6e-05

                                   #获取p值
                                   p.value = 1 - pchisq(model$chisq, length(model$n) - 1)
                                   p.value = signif(p.value,4)
                                   p.sig = ymark_psig(p.value)
                                   #获取HR = A/B = (O(A)/E(B))/(O(B)/E(B))
                                   # here A is feature=1, B is feature=0
                                   HR = (model$obs[2]/model$exp[2])/(model$obs[1]/model$exp[1])
                                   #获取95%置信区间
                                   low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/model$exp[2]+1/model$exp[1]))
                                   up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/model$exp[2]+1/model$exp[1]))
                                   HR_95_CI_lower = signif(low95,2)
                                   HR_95_CI_upper = signif(up95,2)
                                   #feature = r$conf.int %>% row.names

                                   res<-c(mut_cases, wt_cases, ratio, p.value
                                          , p.sig
                                          ,'A(1)/B(0)', HR, HR_95_CI_lower, HR_95_CI_upper,note)
                               }

                               names(res)<-c('mut_cases','wt_cases','mut_ratio%', "p_value"
                                             , 'p_sig'
                                             ,'HR_formular', "HR", 'HR_95_CI_lower', 'HR_95_CI_upper',"note")
                               return(res)
                           })
    #转换成数据框，并转置
    res = t(as.data.frame(univ_results, check.names = FALSE))
    as.data.frame(res) %>% arrange(p_value)
}


ydo_single_factor_coxph2_helper = function(GRPS = NULL){
    " if GRPS == NULL: generate helper message
if GRPS is character: assign a global var ITER and returns a vector COLNAME
"
    if (is.null(GRPS)){
        print(" -- demo -- ")
        print("global var ITER is needed to run ysingle_factor_coxph2")
        print("GRPS = c('Good','Poor','Others');    ITER = combn(GRPS %>% sort,2,simplify = FALSE)")
        print('COLNAME is like:')
        print("COLNAME = c(paste0('cases_',GRPS,'_d'),'N_d',paste0('ratio%_',GRPS,'_d'),'p_value_d', 'p_sig', paste0('pval_',NAMES,'_d'),'Compares','HR.mean_d', paste0('HR_',NAMES,'_d'),paste0('CI_',NAMES),'note')")
    }else{
        ITER = combn(GRPS %>% sort,2,simplify = FALSE)
        NAMES = c()
        for (i in ITER){
            NAMES = c(NAMES,i %>% str_flatten(collapse = '_vs_'))
        }
        COLNAME = c(
            paste0('cases_',GRPS,'_d'),
            'N_d',
            paste0('ratio%_',GRPS,'_d'),
            "p_value_d", "p_sig",
            paste0('pval_',NAMES,'_d'),
            'Compares',
            "HR.mean_d",
            paste0('HR_',NAMES,'_d'),
            paste0('CI_',NAMES),
            "note"
        )
        assign('ITER',ITER)
        COLNAME
    }
}
ydo_single_factor_coxph2 = function (data, covariates, key = NULL, test.method='coxph', join = "_", tags = c("month", "status")) {
    "MAKE SURE that data values consist of int,
    which 0 mean no mut, others mean mut\n
    except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)\n
    mut_cases = column$values != 0\n wt_cases = total - mut_cases\n
    test.method in c('coxph','logrank','likelihood')
        coxph: wald-test
        logrank: logrank-test
"
    `+` = .Primitive("+")
    if (is.null(ITER)){
        ITER = data[,covariates]
    }
    if (!is.null(key)) {
        templ = stringr::str_flatten(c("Surv(", key, join, tags[[1]],
                                       ",", key, join, tags[[2]], ")~"))
    }
    else {
        templ = "Surv(survivalMonth,survivalEvent)~"
    }

    univ_formulas <- sapply(covariates, function(x) {
        res = list()
        formu = as.formula(paste(templ, x))
        res[["formula"]] = formu
        res[["feature"]] = x
        res
    }, simplify = FALSE)
    univ_results <- lapply(univ_formulas, function(li) {
        feature = li$feature
        col = data %>% dplyr::pull(!!feature)
        formu = li$formu

        grps = col %>% unique %>% .[!is.na(.)]

        ps = c()
        cmps = c()
        CIs = c()
        HRs = c()

        N = length(grps)

        nums = c(
            Good=sum(col == 'Good', na.rm = TRUE),
            Poor=sum(col == 'Poor', na.rm = TRUE),
            Others=sum(col=='Others', na.rm = TRUE)
        )
        all_cases = length(col)

        if (all_cases > sum(nums)) {
            note = paste("contains", all_cases - sum(nums), "NA values")
        } else {
            note = ""
        }
        ratios = signif(nums/all_cases,digits = 3) * 100

        for (x in ITER){
            a = x[[1]]
            b = x[[2]]
            if (nums[[a]]==0 | nums[[b]]==0){
                p.value = NA
                HR = NA
                HR_95_CI_lower = NA
                HR_95_CI_upper = NA
            }else{
                d = data[col %in% c(a,b),]
                model = coxph(formu, data = d )
                x <- summary(model)
                if (test.method == "logrank") {
                    p.value = x$sctest[["pvalue"]]
                }
                else if (test.method == "coxph" | test.method=='wald') {
                    p.value = x$waldtest[["pvalue"]]
                }
                else if (test.method == "likelihood") {
                    p.value = x$logtest[["pvalue"]]
                }
                p.value <- signif(p.value, digits = 4)
                HR <- signif(x$coef[2], digits = 2)
                HR_95_CI_lower <- signif(x$conf.int[, "lower .95"], 2)
                HR_95_CI_upper <- signif(x$conf.int[, "upper .95"], 2)
            }
            ps = ps %>% append(p.value)
            cmps = cmps %>% append(paste0(a,'_vs_',b))
            CIs = CIs %>% append(paste0('[',HR_95_CI_lower,' - ',HR_95_CI_upper,']'))
            HRs = HRs %>% append(HR)
        }

        p.value = mean(ps,na.rm = TRUE)
        # all_pval = ps %>% str_flatten(collapse = ', ')
        p.sig = ymark_psig(p.value)
        cmp.str = cmps %>% str_flatten(collapse = ', ')
        # CI = CIs %>% str_flatten(collapse = ', ')
        HR.mean = mean(HRs,na.rm = TRUE)
        # HR = HRs %>% str_flatten(collapse = ', ')


        res <- c(nums, N,ratios,
                 p.value, p.sig , ps, cmp.str,
                 HR.mean, HRs, CIs, note)
        names(res) = NULL
        #         names(res) = COLNAME
        return(res)
    })
    #     as.data.frame(univ_results, check.names = FALSE)
    return(univ_results)
    res = t(as.data.frame(univ_results, check.names = FALSE))  %>% as.data.frame
    colnames(res) = COLNAME
    res
}



ydo_single_factor_coxph3 = function (data, covariates, key = NULL, test.method = "coxph",
                                     join = "_", tags = c("month", "status"))
{
    "MAKE SURE that data values consist of int, \n    which 0 mean no mut, others mean mut\n \n    except for columns c(patient_id,survivalEvent,survivalMonth,Clin_classification)\n \n    mut_cases = column$values != 0\n wt_cases = total - mut_cases\n\n    test.method in c('coxph','logrank','likelihood')\n        coxph: wald-test\n        logrank: logrank-test\n"
    `+` = .Primitive("+")


    if (!is.null(key)) {
        templ = stringr::str_flatten(c("Surv(", key, join, tags[[1]],
                                       ",", key, join, tags[[2]], ")~"))
    }
    else {
        templ = "Surv(survivalMonth,survivalEvent)~"
    }
    univ_formulas <- sapply(covariates, function(x) {
        res = list()
        formu = as.formula(paste(templ, x))
        res[["formula"]] = formu
        res[["feature"]] = x
        res
    }, simplify = FALSE)
    univ_results <- lapply(univ_formulas, function(li) {
        feature = li$feature
        col = data %>% dplyr::pull(!!feature)
        formu = li$formu
        grps = col %>% unique %>% .[!is.na(.)]
        ps = c()
        cmps = c()
        CIs = c()
        HRs = c()
        N = length(grps)
        nums = c()
        for (g in grps){
            nums = c(nums,sum(col == g, na.rm = TRUE))
        }
        # 使用grps_names,为了绕开下述 R语言设计缺陷:
        # 用0作为names的时候,不能读取
        # x = list('a')
        # names(x) = 0 # OK
        # x[[0]] # raise Error: attempt to select less than one element in get1index <real>
        grps_names = grps %>% make.names
        names(nums) = grps_names

        all_cases = length(col)
        if (all_cases > sum(nums)) {
            note = paste("contains", all_cases - sum(nums), "NA values")
        }
        else {
            note = ""
        }
        # R语言设计缺陷:
        # 用0作为names的时候,不能读取
        # x = list('a')
        # names(x) = 0 # OK
        # x[[0]] # raise Error: attempt to select less than one element in get1index <real>
        ratios = signif(nums/all_cases, digits = 3) * 100
        ITER = combn(grps_names,2,simplify = FALSE)
        # print(paste(feature,'-',nums))
        for (x in ITER) {
            x = x %>% sort
            a = x[[1]]
            b = x[[2]]
            if (nums[[a]] == 0 | nums[[b]] == 0) {
                p.value = NA
                HR = NA
                HR_95_CI_lower = NA
                HR_95_CI_upper = NA
            }
            else {
                a_ori = grps[match(a, grps_names)]
                b_ori = grps[match(b, grps_names)]
                d = data[col %in% c(a_ori,b_ori), ]
                model = coxph(formu, data = d)
                x <- summary(model)
                if (test.method == "logrank") {
                    p.value = x$sctest[["pvalue"]]
                }
                else if (test.method == "coxph" | test.method ==
                         "wald") {
                    p.value = x$waldtest[["pvalue"]]
                }
                else if (test.method == "likelihood") {
                    p.value = x$logtest[["pvalue"]]
                }
                p.value <- signif(p.value, digits = 4)
                HR <- signif(x$coef[2], digits = 2)
                HR_95_CI_lower <- signif(x$conf.int[, "lower .95"],
                                         2)
                HR_95_CI_upper <- signif(x$conf.int[, "upper .95"],
                                         2)
            }
            ps = ps %>% append(p.value)
            cmps = cmps %>% append(paste0(b_ori, "_vs_", a_ori))
            CIs = CIs %>% append(paste0("[", HR_95_CI_lower,
                                        " - ", HR_95_CI_upper, "]"))
            HRs = HRs %>% append(HR)
        }
        p.value = mean(ps, na.rm = TRUE)
        p.sig = ymark_psig(p.value)
        cmp.str = cmps %>% str_flatten(collapse = ", ")
        HR.mean = mean(HRs, na.rm = TRUE)
        res <- c(cmp.str, nums, N, ratios,test.method, p.value, p.sig, ps,
                 HR.mean, HRs, CIs, note)
        names(res) = NULL
        return(res)
    })
    # return(univ_results)
    res = t(as.data.frame(univ_results, check.names = FALSE)) %>%
        as.data.frame
    colnames(res) = c('Cmp_Groups','LeftGroup_cases','Right_Group_cases','SubGroupNumbers'
                      ,'LeftGroup_ratio%', 'Right_Group_ratio%','stat_method',"p_value"
                      , 'p_sig','p_val_list'
                      , "HR.mean", 'HR_list', 'HR_95_CI_list',"note")
    res
}


