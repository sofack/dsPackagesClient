#' @title Fit a cox proportional hazard Model (coxph) with pooling via Study Level Meta-Analysis (SLMA)
#' @description Fits a cox proportional hazard Model (coxph) on data from single or multiple sources
#' with pooled co-analysis across studies being based on SLMA (Study Level Meta Analysis).
#' @details \code{ds.coxSLMA} specifies the structure of a cox proportional hazard Model (coxph)
#' to be fitted separately on each study or data source. Calls the serverside functions
#' coxSLMADS1 (aggregate) and coxSLMADS2 (aggregate).ds.coxSLMA sends
#' a command to every data source to fit the model required but each separate source
#' simply fits that model to completion (ie undertakes all iterations until
#' the model converges) and the estimates (regression coefficients) and their Variance covariance
#' matrices from each source are sent back to the client and are then pooled using SLMA
#' across studies using the mixmeta function (from the mixmeta package) using fixed
#' optimisation method.But once the estimates and Variance-covariances are on the clientside, the user
#' can alternatively choose to use another meta-analysis package in any way he/she wishes,
#' to pool the coefficients across studies.
#' In \code{formula} Most shortcut notation for formulas allowed under R's standard \code{coxph()}
#' function is also allowed by \code{ds.coxSLMA}.
#'
#' coxph can be fitted very simply using a formula such as:
#'
#' \deqn{y~a+b+c+d}
#'
#' which simply means fit a coxph with \code{y} as the outcome variable (Survival object
#' calculated separately using ds.Surv and stored in the server) and
#' \code{a}, \code{b}, \code{c} and \code{d} as covariates.
#' By default all such models also include an intercept (regression constant) term.
#' The \code{dataName} argument avoids you having to specify the name of the
#' data frame in front of each covariate in the formula.
#' For example, if the data frame is called \code{DataFrame} you
#' avoid having to write: \eqn{DataFrame$y~DataFrame$a+DataFrame$b+DataFrame$c+DataFrame$d}
#'
#' The \code{checks} argument verifies that the variables in the model are all defined (exist)
#' on the server-site at every study
#' and that they have the correct characteristics required to fit the model.
#' It is suggested to make \code{checks} argument TRUE only if an unexplained
#' problem in the model fit is encountered because the running process takes several minutes.
#'
#' Server functions called: \code{coxSLMADS1}, and \code{coxSLMADS2}.
#' @param formula an object of class formula describing
#' the model to be fitted. For more information see
#' \strong{Details}.
#' @param weights a character string specifying the name of a variable containing
#' prior regression weights for the fitting process. \code{ds.coxSLMA} does not allow a weights vector to be
#' written directly into the coxph formula.
#' @param mixmeta If TRUE and numstudies > 1, the regression coefficient and Variance-Covariance matrix are pooled across
#' studies using fixed-effects meta-analysis framework.
#' @param dataName a character string specifying the name of an (optional) data frame
#' that contains all of the variables in the coxph formula.
#' @param checks logical. If TRUE \code{ds.coxSLMA} checks the structural integrity
#' of the model. Default FALSE. For more information see \strong{Details}.
#' @param maxit a numeric scalar denoting the maximum number of iterations that
#' are permitted before \code{ds.coxSLMA} declares that the model has failed to converge.
#' For more information see \strong{Details}.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return The serverside aggregate functions \code{coxSLMADS1} and \code{coxSLMADS2} return
#' output to the clientside.
#' This is precisely the same as the coxph object that is usually created by a call to coxph() in native R and it
#' contains all the same elements (see help for coxph in native R). Because it is a serverside
#' object, no disclosure blocks apply. However, such disclosure blocks do apply to the information
#' passed to the clientside. In consequence, rather than containing all the components of a
#' standard coxph object in native R, the components of the coxph object that are returned by
#' \code{ds.coxSLMA} include: a mixture of non-disclosive elements of the coxph object
#' reported separately by study included in a list object called \code{output.summary}; and
#' a series of other list objects that represent inferences aggregated across studies.
#' @return the study specific items include:
#' @return \code{coefficients}: a matrix with 5 columns:
#'    \itemize{
#'    \item{First}{: the names of all of the regression parameters (coefficients) in the model}
#'    \item{second}{: the estimated values}
#'    \item{third}{: the exponentials of the estimated values}
#'    \item{fourth}{: corresponding standard errors of the estimated values}
#'    \item{fifth}{: the ratio of estimate/standard error}
#'    \item{sixth}{: the p-value treating that as a standardised normal deviate}
#' }
#' @return \code{formula}: model formula, see description of formula as an input parameter (above).
#' @return \code{CoefMatrix}: the matrix of parameter estimates.
#' @return \code{vcovmatrix}: the variance-covariance matrix of parameter estimates.
#' @return \code{weights}: the name of the vector (if any) holding regression weights.
#' @return \code{Nmissing}: the number of missing observations in the given study.
#' @return \code{Nvalid}: the number of valid (non-missing) observations in the given study.
#' @return \code{Ntotal}: the total number of observations in the given study
#'                        (\code{Nvalid} + \code{Nmissing}).
#' @return \code{data}: equivalent to input parameter \code{dataName} (above).
#' @return \code{call}:  summary of key elements of the call to fit the model.
#' @return \code{na.action}:  chosen method of dealing with missing values. This is
#' usually, \code{na.action = na.omit} - see help in native R.
#' @return \code{iter}: the number of iterations required to achieve convergence
#' of the coxph model in each separate study.
#' @return Once the study-specific output has been returned, \code{ds.coxSLMA}
#' returns a series of lists relating to the aggregated inferences across studies.
#' These include the following:
#' @return \code{num.valid.studies}: the number of studies with valid output
#' included in the combined analysis
#' @author Sofack, Ghislain.(Based on ds.glmSLMA by Paul Burton for DataSHIELD Development Team)

#' @examples
#' \dontrun{
#'
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'
#'   # Example 1: Fitting coxph for survival analysis
#'   # For this analysis we need to load survival data from the server
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING1", driver = "OpalDriver")
#'   builder$append(server = "study2",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'
#'  # Log onto the remote Opal training servers
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D")
#'
#'
#'  # The 'cens' variable should be an interger/ numeric
#'
#' ds.asInteger(x.name = "D$cens",
#'              newobj = "CENS",
#'              datasources = connections)
#'
#' # Create the serverside survival object
#
#' ds.Surv(time = "D$survtime",
#'         event = "D$cens",
#'         newobj = "Survobj"
#'         datasources = connections)
#'
#'
#' ds.coxSLMA(formula = Survobj ~ noise.56 + pm10.16 + bmi.26 + age.60 ,
#'            data = "D",
#'            weights = NULL,
#'            checks = FALSE,
#'            maxit = 20,
#'            datasources = connections)
#'
#'  # Clear the Datashield R sessions and logout
#'
#' datashield.logout(connections)
#'
#' }
#' @export


ds.coxSLMA<-function(formula=NULL, weights=NULL,dataName=NULL, checks=FALSE, maxit=30,
                     combine.with.mixmeta = TRUE, datasources=NULL) {

  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  # verify that 'formula' was set
  if(is.null(formula)){
    stop(" Please provide a valid regression formula!", call.=FALSE)
  }

  # check if user gave weights directly in formula, if so the argument 'weights'
  # to provide name of offset or weights variable
  if(sum(as.numeric(grepl('weights', formula, ignore.case=TRUE)))>0)
  {
    cat("\n\n WARNING: you may have specified a regression weights")
    cat("\n as part of the model formula.")
    cat("\n you must specify  weights separately from the formula")
    cat("\n using the weights argument.\n\n")
  }

  formula <- stats::as.formula(formula)

  # if the argument 'dataName' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(dataName))){
    defined <- dsBaseClient:::isDefined(datasources, dataName)
  }

  #MOVE ITERATION COUNT BEFORE ASSIGNMENT OF beta.vect.next
  #Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration

  iteration.count<-0

  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)

  #ARBITRARY LENGTH FOR START BETAs AT THIS STAGE BUT IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
  beta.vect.next <- rep(0,5)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")

  #Identify the correct dimension for start Betas via calling first component of coxSLMADS

  cally1 <- call('coxSLMADS1', formula, weights, dataName)

  study.summary.0 <- DSI::datashield.aggregate(datasources, cally1)


  at.least.one.study.data.error<-0
  at.least.one.study.valid<-0


  num.par.coxph<-NULL

  coef.names<-NULL

  for(hh in 1:numstudies) {
    if(study.summary.0[[hh]]$errorMessage!="No errors")
    {
      at.least.one.study.data.error<-1
    }else{
      at.least.one.study.valid<-1
      num.par.coxph<-study.summary.0[[hh]][[1]][[2]]
      coef.names<-study.summary.0[[hh]][[2]]
    }
  }

  y.invalid<-NULL
  Xpar.invalid<-NULL
  w.invalid<-NULL
  coxph.saturation.invalid<-NULL
  errorMessage<-NULL

  for(ss in 1:numstudies)
  {
    y.invalid<-c(y.invalid,study.summary.0[[ss]][[3]])
    Xpar.invalid<-rbind(Xpar.invalid,study.summary.0[[ss]][[4]])
    w.invalid<-c(w.invalid,study.summary.0[[ss]][[5]])
    coxph.saturation.invalid <-c(coxph.saturation.invalid,study.summary.0[[ss]][[6]])
    errorMessage<-c(errorMessage,study.summary.0[[ss]][[7]])
  }


  y.invalid<-as.matrix(y.invalid)
  sum.y.invalid<-sum(y.invalid)
  dimnames(y.invalid)<-list(names(datasources),"Y VECTOR")

  Xpar.invalid<-as.matrix(Xpar.invalid)
  sum.Xpar.invalid<-sum(Xpar.invalid)
  dimnames(Xpar.invalid)<-list(names(datasources),coef.names)

  w.invalid<-as.matrix(w.invalid)
  sum.w.invalid<-sum(w.invalid)
  dimnames(w.invalid)<-list(names(datasources),"WEIGHT VECTOR")


  coxph.saturation.invalid<-as.matrix(coxph.saturation.invalid)
  sum.coxph.saturation.invalid<-sum(coxph.saturation.invalid)
  dimnames(coxph.saturation.invalid)<-list(names(datasources),"MODEL OVERPARAMETERIZED")

  errorMessage<-as.matrix(errorMessage)
  dimnames(errorMessage)<-list(names(datasources),"ERROR MESSAGES")


  output.blocked.information.1<-"EVERY STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:"
  output.blocked.information.2<-"Any values of 1 in the following tables denote potential disclosure risks."
  output.blocked.information.3<-"Please use the argument <datasources> to include only valid studies."
  output.blocked.information.4<-"Errors by study are as follows:"


  #CASE 1 - NO STUDIES VALID
  if(!at.least.one.study.valid)
  {
    message("\n\nEVERY STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:\n",
            "Any values of 1 in the following tables denote potential disclosure risks.\n",
            "Errors by study are as follows:\n")

    return(list(
      output.blocked.information.1,
      output.blocked.information.2,
      output.blocked.information.4,
      y.vector.error=y.invalid,
      X.matrix.error=Xpar.invalid,
      weight.vector.error=w.invalid,
      coxph.overparameterized=coxph.saturation.invalid,
      errorMessage=errorMessage
    ))
  }


  #CASE 2 - AT LEAST ONE STUDY VALID AND AT LEAST ONE INVALID
  if(at.least.one.study.data.error)
  {
    message("\n\nAT LEAST ONE STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:\n",
            "Any values of 1 in the following tables denote potential disclosure risks.\n",
            "No analytic results are returned for potentially disclosive studies and\n",
            "pooled co-estimates across studies are based only on the valid studies.\n",
            "You may also choose to exclude invalid studies from\n",
            "the whole analysis using the <datasources> argument.\n",
            "Errors by study are as follows:\n")

  }

  beta.vect.next <- rep(0,num.par.coxph)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")


  #Provide arbitrary starting value for deviance to enable subsequent calculation of the
  #change in deviance between iterations
  dev.old<-9.99e+99

  #Convergence state needs to be monitored.
  converge.state<-FALSE

  #Define a convergence criterion. This value of epsilon corresponds to that used
  #by default for GLMs in R (see section S3 for details)
  epsilon<-1.0e-08

  f<-NULL

  #Now calling the second component of coxSLMADS to generate score vectors and information matrices

  cally2 <- call('coxSLMADS2', formula, weights, dataName)

  study.summary <- DSI::datashield.aggregate(datasources, cally2)


  #NOW ONLY WORKING WITH SITUATIONS WITH AT LEAST ONE VALID STUDY

  numstudies<-length(datasources)


  #IF combine.with.mixmeta == TRUE, FIRST CHECK THAT THE MODELS IN EACH STUDY MATCH
  #IF THERE ARE DIFFERENT NUMBERS OF PARAMETERS THE ANALYST WILL
  #HAVE TO USE THE RETURNED MATRICES FOR coefs AND vcovs TO DETERMINE WHETHER
  #COMBINATION ACROSS STUDIES IS POSSIBLE AND IF SO, WHICH PARAMETERS GO WITH WHICH
  #ALSO DETERMINE WHICH STUDIES HAVE VALID DATA

  if (combine.with.mixmeta == TRUE & numstudies > 1){

    study.include.in.analysis <- NULL
    study.with.errors<-NULL
    all.studies.valid<-1
    no.studies.valid<-1

    #MAKE SURE THAT IF SOME STUDIES HAVE MORE PARAMETERS IN THE
    #FITTED coxph (eg BECAUSE OF ALIASING) THE FINAL RETURN MATRICES
    #HAVE ENOUGH ROWS TO FIT THE MAXIMUM LENGTH

    numcoefficients.max<-0

    for(g in numstudies){
      if(length(study.summary[[g]]$coefficients[,1])>numcoefficients.max){
        numcoefficients.max<-length(study.summary[[g]]$coefficients[,1])
      }
    }

    numcoefficients<-numcoefficients.max

    coefmatrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)

    vcovmatrix<- NULL


    for(k in 1:numstudies){
      coefmatrix[,k]<-study.summary[[k]]$coefficients[,1]

      vcovmatrix[k]<-list(study.summary[[k]]$vcov)
    }

    # transpose the coefmatrix
    coefmatrix = t(coefmatrix)

    #Annotate output matrices with study indicators

    study.names.list<-NULL
    coef.study.names.list<-NULL
    vcov.study.names.list<-NULL

    for(v in 1:numstudies){

      study.names.list<-c(study.names.list,paste0("study",as.character(v)))
      coef.study.names.list<-c(coef.study.names.list,paste0("coef study ",as.character(v)))
      vcov.study.names.list<-c(vcov.study.names.list,paste0("vcov study ",as.character(v)))
    }

    colnames(coefmatrix) = dimnames(study.summary[[1]]$coefficients)[[1]]
    rownames(coefmatrix) =  coef.study.names.list

    output.summary.text<-paste0("list(")

    for(u in 1:numstudies){
      output.summary.text<-paste0(output.summary.text,"study",as.character(u),"=study.summary[[",as.character(u),"]],"," ")
    }

    output.summary.text.save<-output.summary.text
    output.summary.text<-paste0(output.summary.text,"input.coef.matrix.for.SLMA=as.matrix(coefmatrix),input.vcov.matrix.for.SLMA=as.list(vcovmatrix))")


    output.summary<-eval(parse(text=output.summary.text))


    ##########END OF ANNOTATION CODE ################


    ## MULTIVARIATE METAANALYSE

    coef.matrix.for.SLMA<-as.matrix(coefmatrix)
    vcov.matrix.for.SLMA<-as.list(vcovmatrix)

    #SELECT VALID COLUMNS ONLY (THERE WILL ALWAYS BE AT LEAST ONE)

    usecols<-NULL

    for(ut in 1:(dim(coef.matrix.for.SLMA)[2]))
    {
      if(!is.na(coef.matrix.for.SLMA[1,ut])&&!is.null(coef.matrix.for.SLMA[1,ut]))
      {
        usecols<-c(usecols,ut)
      }
    }

    coefmatrix.valid<-coef.matrix.for.SLMA[,usecols]

    usecols1<-NULL

    for(ut1 in 1:length(vcov.matrix.for.SLMA))
    {
      if(!is.na(vcov.matrix.for.SLMA[ut])&&!is.null(vcov.matrix.for.SLMA[ut]))
      {
        usecols1<-c(usecols1,ut1)
      }
    }
    vcovmatrix.valid<-vcov.matrix.for.SLMA[usecols1]


    #Check for matched parameters

    num.valid.studies<-as.numeric(dim(as.matrix(coefmatrix.valid))[1])

    coefficient.vectors.match<-TRUE

    if(!coefficient.vectors.match){
      cat("\n\nModels in different sources vary in structure\nplease match coefficients for meta-analysis individually\n")
      cat("nYou can use the DataSHIELD generated estimates and vcov matrix as the basis for a meta-analysis\nbut carry out the final pooling step independently of DataSHIELD using whatever meta-analysis package you wish\n\n")
      return(list(output.summary=output.summary))
    }


    ## Mixmeta analysis

    #mix <- mixmeta::mixmeta(formula, S = NULL,  method= NULL)

    mix.fixed <- mixmeta::mixmeta(formula = coefmatrix.valid, S = vcovmatrix.valid, method = "fixed")
    # mix.ml <- mixmeta::mixmeta(formula = coefmatrix.valid, S = vcovmatrix.valid, method = "ml")
    # mix.reml <- mixmeta::mixmeta(formula = coefmatrix.valid, S = vcovmatrix.valid, method = "reml")
    # mix.mm <- mixmeta::mixmeta(formula = coefmatrix.valid, S = vcovmatrix.valid, method = "mm")
    # mix.vc <- mixmeta::mixmeta(formula = coefmatrix.valid, S = vcovmatrix.valid, method = "vc")

    return(list(output.summary=output.summary,num.valid.studies = num.valid.studies,
                coefmatrix.valid=coefmatrix.valid, vcovmatrix.valid = vcovmatrix.valid,
                Fixed = summary(mix.fixed)))
  }

  else {
    return(study.summary)
  }

}
# ds.coxSLMA
