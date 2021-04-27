#Code adapted from https://github.com/kdyson/R_Scripts/blob/master/AICc_PERMANOVA.R by Mark. A. F. Gillingham.

# Function to calculate AICc for PERMANOVA. Requires input from adonis {vegan}

AICc.PERMANOVA <- function(adonis.model) {
    
    # check to see if object is an adonis model...
    
    if (!(adonis.model$aov.tab[1,1] >= 1))
        stop("object not output of adonis {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model
    # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
    
    RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "SumsOfSqs"]
    MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
    
    k <- ncol(adonis.model$model.matrix)# + 1 # add one for error variance
    
    nn <- nrow(adonis.model$model.matrix)
    
    # AIC : 2*k + n*ln(RSS)
    # AICc: AIC + [2k(k+1)]/(n-k-1)

    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    

    loglik<--(nn+(nn*log(2*pi))+(nn*log(RSS/nn)))/2
    AIC<-nn+(nn*log(2*pi))+(nn*log(RSS/nn))+(2*(k+1))
    BIC<-nn+(nn*log(2*pi))+(nn*log(RSS/nn))+(log(nn)*(k+1))
    AICc<-AIC+(2*k*(k + 1))/(nn - k - 1)


    output <- list("AIC" = AIC, "AICc" = AICc, "k" = k, "loglik"=loglik,"BIC"=BIC)
    
    return(output)   

}


AICc.NULL <- function(adonis.model) {
    
    # check to see if object is an adonis model...
    
    if (!(adonis.model$aov.tab[1,1] >= 1))
        stop("object not output of adonis {vegan} ")
    
    # Ok, now extract appropriate terms from the adonis model
    # Calculating AICc using residual sum of squares (RSS) since I don't think that adonis returns something I can use as a liklihood function...
    
    RSS <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Total", "SumsOfSqs"]
    #MSE <- adonis.model$aov.tab[rownames(adonis.model$aov.tab) == "Residuals", "MeanSqs"]
    
    k <- 1# + 1 # add one for error variance
    
    nn <- nrow(adonis.model$model.matrix)
    
    # AIC : 2*k + n*ln(RSS)
    # AICc: AIC + [2k(k+1)]/(n-k-1)

    # based on https://en.wikipedia.org/wiki/Akaike_information_criterion;
    # https://www.researchgate.net/post/What_is_the_AIC_formula;
    # http://avesbiodiv.mncn.csic.es/estadistica/ejemploaic.pdf
    
    # AIC.g is generalized version of AIC = 2k + n [Ln( 2(pi) RSS/n ) + 1]
    # AIC.pi = k + n [Ln( 2(pi) RSS/(n-k) ) +1],
    
    
    loglik<--(nn+(nn*log(2*pi))+(nn*log(RSS/nn)))/2
    AIC<-nn+(nn*log(2*pi))+(nn*log(RSS/nn))+(2*(k+1))
    BIC<-nn+(nn*log(2*pi))+(nn*log(RSS/nn))+(log(nn)*(k+1))
    AICc<-AIC+ (2*k*(k + 1))/(nn - k - 1)

    output <- list("AIC" = AIC, "AICc" = AICc, "k" = k, "loglik"=loglik,"BIC"=BIC)
    
    return(output)   

}


