#' rFSA: Feasible Solution Algorithm (FSA) for Linear Models
#' @description  A function using a Feasible Solution Algorithm to find a set of feasible solutions for a linear model of a specific form that could include mth-order interactions (Note that these solutions are optimal in the sense that no one swap to any of the variables will increase the criterion function.)
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. See help(lm) for details.
#' @param data a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. 
#' @param fixvar a variable to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.
#' @param quad to include quadratic terms or not.
#' @param m order of terms to potentially include. If interactions is set to TRUE then m is the order of interactions to be considered. Defaults to 2. For Subset selection (interaction=F), m is the size of the subset to examine. Default is 2. 
#' @param numrs number of random starts to perform.
#' @param cores number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details. If function detects a Windows user it will automatically set cores=1. 
#' @param interactions T or F for whether to include interactions in model. Defaults to FALSE. 
#' @param criterion which criterion function to either maximize or minimize. For linear models one can use: r.squared, adj.r.squared, cv5.lmFSA (5 Fold Cross Validation error), cv10.lmFSA (10 Fold Cross Validation error), apress (Allen's Press Statistic), int.p.val (Interaction P-value), AIC, BIC.
#' @param minmax whether to minimize or maximize the criterion function
#' @param ... arguments to be passed to the lm function
#' @details PLEASE NOTE: make sure categorical variables are factors or characters otherwise answers will not reflect the variable being treated as a continuous variable.
#' @return returns a list of solutions and table of unique solutions.
#' $solutions is a matrix of fixed terms, start position, feasible solution, criterion function value (p-value of interaction), and number of swaps to solution.
#' $table is a matrix of the unique feasible solutions and how many times they occured out of the number of random starts chosen. It also returns any warning messages with these solutions in the last column.
#' $efficiency is text comparing how many models you ran during your FSA search compared to how many you would have done with exhaustive search. Note: The FSA algorithm takes additional time to run on top of the model checks that were done during the algorithm. This additional time is approximately 15% more time than if you had just ran the model checks. 
#' @importFrom parallel mclapply
#' @importFrom graphics par plot
#' @importFrom methods show
#' @importFrom stats AIC anova as.formula cov fitted formula glm influence lm predict resid var
#' @importFrom utils capture.output tail
#' @export
#' @examples
#' #use mtcars package see help(mtcars)
#' data(mtcars)
#' colnames(mtcars)
#' fit<-lmFSA(formula="mpg~hp+wt",data=mtcars,fixvar="hp",
#'                quad=FALSE,m=2,numrs=10,cores=1)
#' print(fit) #print formulas of fitted models
#' summary(fit) #review
lmFSA = function(formula,data,fixvar = NULL,quad = FALSE,m = 2,numrs = 1,
                 cores = 1,interactions = TRUE,criterion = r.squared,minmax = "max",...) {
  if(identical(criterion,bdist)){return(show("Sorry the criterion function you listed cannot be used with lmFSA."))}
  formula <- as.formula(formula)
  fit <- lm(formula,data = data,...)
  yname <- all.vars(formula)
  if (!all(c(yname,fixvar) %in% colnames(data))) {
    return(
      show(
        "Sorry, one of the variables you specified in your formula or fixvar is not a name for a column in the data you specified. Please try again."
      )
    )
  }
  
  originalnames <- colnames(data)
  data <- data.frame(data)
  lhsvar <- yname[1]
  
  if (.Platform$OS.type == "unix") {
  } else {
    cores = 1
  }
  
  
  ypos <- which(colnames(data) == lhsvar)
  startvar <- NULL
  xdata <- data[,-ypos]
  ydata <- data[,ypos]
  newdata <- data.frame(cbind(ydata,xdata))
  fixpos <- which(colnames(xdata) %in% fixvar)
  if (length(fixpos) == 0) {
    fixpos = NULL
  }
  
  history <- matrix(rep(NA,numrs * (2 * m + 3)),ncol = ((2 * m + 3)))
  history[,1:m] <- rstart(m = m,nvars = (dim(newdata)[2] - 1),numrs = numrs)
  curpos <- which(colnames(xdata) %in% startvar[-1])
  if (length(curpos) != 0) {
    history <- rbind(c(curpos,rep(NA,length(curpos) + 2)),history)
  }
  
  fsa <- function(i,history,...) {
    cur <- history[i,1:m]
    last <- rep(NA,m)
    numswap <- 0
    memswap <- NULL
    if (minmax == "max") {
      last.criterion <- (-Inf)
    }
    if (minmax == "min") {
      last.criterion <- (Inf)
    }
    checks <- 0
    while (!identical(cur,last) && !identical(c(cur[2],cur[1]),last)) {
      last <- cur
      if (numswap == 0) {
        moves <- swaps(cur = cur,n = dim(xdata)[2],quad = quad)
      }
      if (numswap > 0) {
        moves <-
          nextswap(
            curpos = cur,n = dim(xdata)[2],quad = quad,prevpos = memswap
          )$nswaps
      }
      if (dim(moves)[2] == 0) {
        moves <- t(t(last))
      }
      if (interactions == T) {
        form <-
          function(j)
            formula(paste0(
              colnames(newdata)[1],"~",paste0(fixvar,sep = "+"),paste(colnames(xdata)[moves[,j]],collapse = "*")
            ),sep = "")
      }
      if (interactions == F) {
        form <-
          function(j)
            formula(paste0(
              colnames(newdata)[1],"~",paste0(fixvar,sep = "+"),paste(colnames(xdata)[moves[,j]],collapse = "+")
            ),sep = "")
      }
      tmp <-
        parallel::mclapply(
          X = 1:dim(moves)[2],FUN = function(k)
            criterion(lm(form(k),data = newdata,...)),mc.cores = cores
        )
      checks <- checks + dim(moves)[2]
      if (minmax == "max") {
        cur <- moves[,which.max.na(unlist(tmp))[1]]
        cur.criterion <- unlist(tmp[which.max.na(unlist(tmp))[1]])
        if (last.criterion > cur.criterion) {
          cur <- last.pos
          cur.criterion <- last.criterion
        }
      }
      if (minmax == "min") {
        cur <- moves[,which.min.na(unlist(tmp))[1]]
        cur.criterion <- unlist(tmp[which.min.na(unlist(tmp))[1]])
        if (last.criterion < cur.criterion) {
          cur <- last.pos
          cur.criterion <- last.criterion
        }
      }
      numswap <- numswap + 1
      last1 <- last
      last.criterion <- cur.criterion
      last.pos <- cur
      memswap <- unique(c(memswap,last1))
    }
    history[i,(1 + m):(2 * m)] <- cur
    history[i,(dim(history)[2] - 2)] <- cur.criterion
    history[i,(dim(history)[2] - 1)] <- numswap - 1
    history[i,(dim(history)[2])] <- checks
    return(history[i,])
  }
  solutions <-
    matrix(unlist(lapply(
      1:numrs,FUN = function(i)
        fsa(i,history)
    )),ncol = dim(history)[2],byrow = T)
  solutions[,1:(2 * m)] <-
    matrix(colnames(newdata)[c(solutions[,1:(2 * m)] + 1)],ncol = (2 * m))
  solutions <- data.frame(solutions)
  colnames(solutions) <-
    c(
      paste("start",1:m,sep = "."),paste("best",1:m,sep = "."),"criterion","swaps","checks"
    )
  solutions$criterion <-
    as.numeric(levels(solutions$criterion))[solutions$criterion]
  solutions$swaps <-
    as.numeric(levels(solutions$swaps))[solutions$swaps]
  solutions$checks <-
    as.numeric(levels(solutions$checks))[solutions$checks]
  if (length(fixvar) != 0) {
    solutions <-
      data.frame(fixvar = matrix(
        rep(x = fixvar,dim(solutions)[1]),nrow = dim(solutions)[1],byrow = T
      ),solutions)
  }
  solutions <- solutions
  a <- solutions[,(length(fixvar) + m + 1):(length(fixvar) + m + 1 + m)]
  b <- unique(t(apply(a,sort,MARGIN = 1)),MARGIN = 1)
  a <- t(apply(a,sort,MARGIN = 1))
  c <- cbind(b,0)
  for (i in 1:dim(b)[1]) {
    for (j in 1:dim(a)[1]) {
      c[i,(m + 2)] <-
        sum(as.numeric(c[i,(m + 2)]) + as.numeric(identical(a[j,],b[i,])))
    }
  }
  tableres <- data.frame(cbind(c),stringsAsFactors = F)
  colnames(tableres)[(dim(tableres)[2])] <- "times"
  colnames(tableres)[2:(dim(tableres)[2] - 1)] <-
    paste("Var",1:m,sep = "")
  colnames(tableres)[1] <- "criterion"

  call <- mget(names(formals()),sys.frame(sys.nframe()))
  ls <-
    list(
      originalfit = fit,call = call,solutions = solutions,table = tableres,efficiency =
        paste(
          "You did:",sum(solutions$checks)," model checks compared to ",choose(n = dim(xdata)[2],k = m)," checks you would have done with exahstive search."
        )
    )
  class(ls) <- "FSA"
  invisible(ls)
  return(ls)
}
