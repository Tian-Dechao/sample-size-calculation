rm(list=ls())
library(pwr)
library(MASS)
library(ggplot2)
library(reshape2)
library(lme4)
library(nlme)
suppressMessages(library(docopt))
suppressMessages(source("/hive/dechaot/Coevolve/src/main/r/myfunctions.R"))

script.name = scriptname()

"
scriptname
Calculate sample size and do power analysis. 

Usage:
    scriptname [chi | t] [options]
    scriptname methylo (frideman | lmm) <N> <t> <delta> <r>  [options]
    scriptname methyloFig <i> <o> [options]

Options:
    --help <> help files
    --xlim <> help files
" -> opt

opt = gsub('scriptname', script.name, opt)
opts = docopt(opt)

if(opts$chi){
    w = 0.8
    #df=14
    #df=2
    df=6
    alpha=0.05
    power=0.95
    N = pwr.chisq.test(w = w, df = df, sig.level = alpha, power = power)$N
    print(N)
    print(floor(N/2  + 1) * 2)
    alpha = alpha/100
    power = pwr.chisq.test(w = w, N = floor(N/2 +1) * 2, df = df, sig.level = alpha)$power
    print(power)
    alpha = alpha/10
    power = pwr.chisq.test(w = w, N = floor(N/2 +1) * 2, df = df, sig.level = alpha)$power
    print(power)
}

correlation_matrix = function(n, decay=0.9){
    c = matrix(0, ncol=n, nrow=n)
    for(i in 1:n){
        for(j in i:n){
            if(i==j){ c[i, j]=1 }
            else{
                d = j-i
                c[i,j] = decay^d
                c[j,i] = decay^d
            }
        }
    }
    return(c)
}

generate_data = function(N, t, delta){
    m = cumsum(c(0, rep(delta, t-1)))
    sigma = correlation_matrix(n=t, decay=0.95) 
    # generate the data
    subj = as.factor(rep(1:N, t))
    Timepoint = rep(1:t, rep(N, t))
    Y = mvrnorm(n=N, mu=m, Sigma=sigma)
    # as.vector bycol bydefault
    Y.vect = as.vector(Y)
    Data = data.frame(Y=Y.vect, subj=subj, time=Timepoint, time2 = as.factor(Timepoint), stringsAsFactors = F)
    return(Data)
    # prepare the model and return pvalue
    #pval = friedman.test(Y~time|subj, data=Data)$p.value
    #return(pval)
}

frideman_test_1 = function(df){
    pval = friedman.test(Y~time|subj, data=df)$p.value
    return(pval)
}
lmm_test_1 = function(df){
    # reference 
    # A very basic tutorial for performing linear mixed effects analyses
    #  pval:  https://stats.stackexchange.com/questions/22988/how-to-obtain-the-p-value-check-significance-of-an-effect-in-a-lme4-mixed-mode
    #lmm.fit = lmer(Y~time + (1|subj), data=df)
    #print(fixef(lmm.fit))
    #print(summary(lmm.fit))
    m1 = lme(Y~time2, random=~1|subj, data=df)
    m1.coef = anova(m1)
    pval = m1.coef[2, 'p-value']
}

if(opts$methylo){
    "
    Repeat-measure with 3-4 time points
    Data from the same patients are correlated 
    differentiating CTC samples
    collected at 3-4 time points during therapy for patients
    Each time point have unique average signals
    must specify the variances (or standard deviations) and the correlations among the repeated measurements
    correla-tions among repeated measurements decline exponen- tially with time
    ANOVA could be help
    Global methylation levels
    "
    N = as.numeric(opts$N)
    t = as.numeric(opts$t)
    delta = as.numeric(opts$delta)
    r = as.numeric(opts$r)
    alpha = c(0.01, 10^(-3), 10^(-4), 10^(-5))
    # generate the data 
    if(opts$frideman){
        pvals = sapply(1:r, function(z){
            data = generate_data(N=N, t=t, delta=delta)
            pvals = frideman_test_1(df=data)
            return(pvals)
        })
    } else if(opts$lmm){
        pvals = sapply(1:r, function(z){
            data = generate_data(N=N, t=t, delta=delta)
            pvals = lmm_test_1(df=data)
            return(pvals)
        })
    } else {
        print('Check the test method')
        q(save='')
    }
    # choose the test, Frideman vs linear mixed effect model
    #pvals = sapply(1:r, function(z) frideman_test_1(N=N, t=t, delta=delta))
    pow = c()
    for(i in 1:length(alpha)){
        pow[i] = sum(pvals <= alpha[i]) / r
    }
    res = c(delta, pow)
    write(res, file='', sep='\t')
}
nl = c('alpha_1'=1e-2,
          'alpha_2'=1e-3,
         'alpha_3'=1e-4,
        'alpha_4'=1e-5)
rename_features = function(x, nl){
    y = x
    names(y) = x
    y[names(nl)] = nl
    return(y)

}
if(opts$methyloFig){
    x = read.table(opts$i, sep='\t')
    # different alpha corrresponding to multiple 
    colnames(x) = c('delta', paste('alpha_', 1:4, sep=''))
    print(x[x[,1] == 0.3, ])
    x.long = melt(x, id.vars=c('delta'))
    if(!is.null(opts$xlim)){
        xlim2 = as.numeric(opts$xlim)
        ind = x.long[, 'delta'] <= xlim2
        x.long = x.long[ind, ]
    }
    #delta.sub = seq(0.1, 0.5, by=0.01)
    #x.long.sub = x.long[x.long[, 'delta'] %in% delta.sub, ]
    pdf(opts$o, width=3, height=3)
    #plot(x.long.sub[, c('delta', 'value')])
    fig = ggplot(data=x.long, aes(x=delta, y=value, colour=variable, group=variable)) +
      geom_line() +
      xlab('Mean difference in methylation\nbetween 2 time points') + ylab('Power') +
      scale_color_discrete(name="Alpha", labels=rename_features('', nl)) +
      theme(legend.position = c(1,0), legend.justification = c(1, 0))

    print(fig)
    dev.off()
}
