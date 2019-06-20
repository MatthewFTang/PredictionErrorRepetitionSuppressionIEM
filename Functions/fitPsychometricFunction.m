

function  logLikelihood= fitPsychometricFunction(pGuess,results)
    
y = Weibull(pGuess,results.intensity)

logLikelihood = sum(results.response.*log(y) + (1-results.response).*log(1-y))


%%
