\name{pasad_test}
\alias{pasad_test}
\title{
A test function of Process-Aware Stealthy Attack Detection(PASAD)
}
\description{
Singular value decomposition of log covariance matrix (Trajectory matrix). This is a test phase of pasad
}
\usage{
pasad_test(obj, test_idx, newdata, r = 3, calib = 1, thres = NULL, movn = 10, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A pasad model object result from \code{pasad_train}.}
  \item{test_idx}{A test index for pasad. For example, \code{test_idx = c(1:4801)}.}
  \item{newdata}{A test data for pasad.}
  \item{r}{A cardinal number of eigen value. Generally, r is same by the training phase.}
  \item{calib}{A threshold multiplier. since the threshold calculated during the training phase is a rather small, so the calib can be multiplied to the threshold. (default is 1, should be biggner than 0).}
  \item{thres}{If user already know the threshold, user can input the threshold.}
  \item{movn}{A moving average window value to see the trend of the score.}
  \item{plot}{Whether to draw a anomaly score or not. (default : TRUE).}
}

\value{
An object of class \code{pasad_test}.

  \item{total_scores}{Anomaly score results from the pasad.}
  \item{tr_score}{A pasad score of the training index.}
  \item{te_score}{A pasad score of the test index.}
  \item{extraced}{A noise-reduced signal.}
  \item{threshold}{A (calculated) threshold.}
  \item{outidx}{An index of a anomaly score greater than the threshold.}
  
}
\references{
Wissam Aoudi, Mikel Iturbe, and Magnus Almgren (2018) <DOI:10.1145/3243734.3243781>. \emph{In Proceedings of the 2018 ACM SIGSAC Conference on Computer and Communications Security}, pp. 817-831.
\cr
\url{https://github.com/mikeliturbe/pasad}
\url{https://github.com/rahulrajpl/PyPASAD}
}

\author{
Donghwan Kim\cr
\email{donhkim9714@korea.ac.kr}
\email{dhkim2@bistel.com}
}



\seealso{\code{\link{pasad_train}}}
\examples{
# data input
fpath = system.file("extdata", "sa.csv", package="pasadr")
sa = read.csv(fpath)

# check data
sig = sa$V5
plot(sig)

# training using pasad and check the scree plot
train_idx = c(1:500)  
obj = pasad_train(x = sig, 
                  train_idx = train_idx,
                  r = 1, 
                  ws = length(train_idx)/2,
                  scree_plot = TRUE)
                  
# test whole data and check the anomaly score of test data
pred = pasad_test(obj = obj, 
                  test_idx = 1:4801,
                  newdata = sig, 
                  r = 1, 
                  calib = 1,  
                  plot = TRUE)
                  
# check the structure of test results
str(pred)

}
\keyword{pasad}
\keyword{trajectory matrix}
