#' @title Client side function calling \code{SurvDS}
#' @description This function is similar to R function \code{Surv}.
#' @details This function calculates the survival object
#' usually used as a response variable in the cox proportional hazard model formula..
#' @param time a numeric vector indicating the start or follow up time.
#' @param time2 a numeric vector indicating the ending time of the interval for interval censored or counting process data only.
#' @param event a numeric vector providing the name of the status indicator(event).Usually binary
#' e.g 0=alive, 1=dead.For multiple endpoint data the event variable will be a factor, whose
#' first level is treated as censoring. Although unusual, the event indicator can be omitted,
#' in which case all subjects are assumed to have an event.
#' @param type character string specifying the type of censoring.Possible
#' values are "right", "left", "counting", "interval", "interval2" or "mstate"..
#' @param newobj a character string specifying the name of the object to which the survival object
#' on the serverside in each study is to be written.
#' If no <newobj> argument is specified, the output
#' object defaults to "Survobj".
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return the results of the survival object stored on the server.
#' @author Sofack, Ghislain. (based on corTestDS by Demetris Avraam, for DataSHIELD Development Team)
#' @export
#'

ds.Surv <- function( time = NULL,
                     time2= NULL,
                     event = NULL,
                     type = NULL,
                     newobj= NULL,
                     datasources=NULL) {

  # if no opal login details are provided look for 'opal' objects in the environment

  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }

  if(is.null(time)){
    stop("Please provide a valid survival start time!", call.=FALSE)
  }

  if(is.null(event)){
    stop("Please provide a valid event parameter!", call.=FALSE)
  }

  # the input variable might be given as column table (i.e. D$object)
  # or just as a vector not attached to a table (i.e. object)
  # we have to make sure the function deals with each case

  objects <- c(time , event)
  xnames <- dsBaseClient:::extract(objects)
  varnames <- xnames$elements
  obj2lookfor <- xnames$holders

  # check if the input object(s) is(are) defined in all the studies
  for(i in 1:length(varnames)){
    if(is.na(obj2lookfor[i])){
      defined <- dsBaseClient:::isDefined(datasources, varnames[i])
    }else{
      defined <- dsBaseClient:::isDefined(datasources, obj2lookfor[i])
    }
  }

  # call the internal function that checks the input object(s) is(are) of the same class in all studies.
  for(i in 1:length(objects)){
    typ <- dsBaseClient:::checkClass(datasources, objects[i])
  }

  # create a name by default if user did not provide a name for the new variable
  if(is.null(newobj)){
    newobj <- "Survobj"
  }

  # call the server side function
  calltext <- call("SurvDS", time, time2, event, type)

  DSI::datashield.assign(datasources, newobj, calltext)



  #############################################################################################################

  #DataSHIELD CLIENTSIDE MODULE: CHECK KEY DATA OBJECTS SUCCESSFULLY CREATED                                  #

  #SET APPROPRIATE PARAMETERS FOR THIS PARTICULAR FUNCTION                                                 	#
  test.obj.name<-newobj


  # CALL SEVERSIDE FUNCTION                                                                                	#
  calltext <- call("testObjExistsDS", test.obj.name)


  object.info<-DSI::datashield.aggregate(datasources, calltext)


  # CHECK IN EACH SOURCE WHETHER OBJECT NAME EXISTS
  # AND WHETHER OBJECT PHYSICALLY EXISTS WITH A NON-NULL CLASS

  num.datasources<-length(object.info)


  obj.name.exists.in.all.sources<-TRUE
  obj.non.null.in.all.sources<-TRUE


  for(j in 1:num.datasources){
    if(!object.info[[j]]$test.obj.exists){
      obj.name.exists.in.all.sources<-FALSE
    }
    if(is.null(object.info[[j]]$test.obj.class) || object.info[[j]]$test.obj.class=="ABSENT"){														 	#
      obj.non.null.in.all.sources<-FALSE
    }
  }

  if(obj.name.exists.in.all.sources && obj.non.null.in.all.sources){

    return.message<-
      paste0("A data object <", test.obj.name, "> has been created in all specified data sources")

  }else{

    return.message.1<-
      paste0("Error: A valid data object <", test.obj.name, "> does NOT exist in ALL specified data sources")

    return.message.2<-
      paste0("It is either ABSENT and/or has no valid content/class,see return.info above")

    return.message.3<-
      paste0("Please use ds.ls() to identify where missing")


    return.message<-list(return.message.1,return.message.2,return.message.3)

  }

  calltext <- call("messageDS", test.obj.name)
  studyside.message<-DSI::datashield.aggregate(datasources, calltext)

  no.errors<-TRUE
  for(nd in 1:num.datasources){
    if(studyside.message[[nd]]!="ALL OK: there are no studysideMessage(s) on this datasource"){
      no.errors<-FALSE
    }
  }


  if(no.errors){
    validity.check<-paste0("<",test.obj.name, "> appears valid in all sources")
    return(list(is.object.created=return.message,validity.check=validity.check))
  }

  if(!no.errors){
    validity.check<-paste0("<",test.obj.name,"> invalid in at least one source. See studyside.messages:")
    return(list(is.object.created=return.message,validity.check=validity.check,
                studyside.messages=studyside.message))
  }

  ###END OF CHECK OBJECT CREATED CORECTLY MODULE
}
#ASSIGN FUNCTION
# ds.Surv
