library(gamlss)

# runs a zero-inflated beta distribution GAM using gamlss
# returns the model
run_mtxDE <- function(formula, data){
    mod <- gamlss::gamlss(as.formula(formula),
                        data=data,
                        family=gamlss.dist::BEZI)
    return(mod)
}
