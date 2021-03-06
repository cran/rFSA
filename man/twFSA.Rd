% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twFSA.R
\name{twFSA}
\alias{twFSA}
\title{twFSA}
\usage{
twFSA(
  formula,
  data,
  fitfunc = lm,
  fixvar = NULL,
  quad = FALSE,
  cores = 1,
  criterion = AIC,
  minmax = "min",
  checkfeas = NULL,
  var4int = NULL,
  min.nonmissing = 1,
  ...
)
}
\arguments{
\item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{a data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{fitfunc}{the method that should be used to fit the model. For Example: lm, glm, or other methods that rely on formula, data, and other inputs.}

\item{fixvar}{variable(s) to fix in the model. Usually a covariate that should always be included (Example: Age, Sex). Will still consider it with interactions. Default is NULL.}

\item{quad}{Include quadratic terms or not. Logical.}

\item{cores}{number of cores to use while running. Note: Windows can only use 1 core. See mclapply for details. If function detects a Windows user it will automatically set cores=1.}

\item{criterion}{which criterion function to either maximize or minimize. For linear models one can use: r.squared, adj.r.squared, cv5.lmFSA (5 Fold Cross Validation error), cv10.lmFSA (10 Fold Cross Validation error), apress (Allen's Press Statistic), int.p.val (Interaction P-value), AIC, BIC.}

\item{minmax}{whether to minimize or maximize the criterion function}

\item{checkfeas}{vector of variables that could be a feasible solution. These variables will be used as the last random start.}

\item{var4int}{specification of which variables to check for marginal feasiblilty. Default is NULL}

\item{min.nonmissing}{the combination of predictors will be ignored unless this many of observations are not missing}

\item{...}{other arguments passed to fitfunc.}
}
\value{
matrix of results
}
\description{
A function for termwise feasiblity
}
