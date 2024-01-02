library(gamlss)

# runs a zero-inflated beta distribution GAM using gamlss
# returns the model
run_model <- function(data){
    mod<-gamlss(outcome~1, sigma.formula=~1, nu.formula=~1,
                data=data,
                family=BEZI)
    return(mod)
}
