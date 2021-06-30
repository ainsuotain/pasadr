\name{pasad_train}
\alias{pasad_train}
\title{
A training function of Process-Aware Stealthy Attack Detection(PASAD)
}
\description{
Singular value decomposition of log covariance matrix (Trajectory matrix). This is a training phase of pasad
}
\usage{
pasad_train(x, train_idx, r = 3, ws, scree_plot = FALSE)
}

\arguments{
  \item{x}{A signal data for inspectation.}
  \item{train_idx}{A training index fot pasad. For example, \code{train_idx = c(1:500)}.}
  \item{r}{A cardinal number of eigen value. Generally r is smaller than 3. (default : 3).}
  \item{ws}{A length of lag for creating covariance matrix. (default is a half of training length).}
  \item{scree_plot}{Whether to draw a scree_plot or not. (default : TRUE).}
}

\value{
An object of class \code{pasad_train}.
  \item{N}{A length of signal data.}
  \item{L}{A length of lag for creating covariance matrix.}
  \item{U}{The r leading eigenvectors.}
  \item{X}{A trajectory matrix.}
  \item{x}{An original signal.}
  \item{ws}{A length of lag for creating covariance matrix.}
  \item{train_idx}{A training index fot pasad.}
  \item{x_train}{A data used for training.}
  \item{singulars}{A transpose of singular matrix}
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
\seealso{\code{\link{pasad_test}}}
\examples{
# data input
fpath = system.file("extdata", "sa.csv", package="pasadr")
sa = read.csv(fpath)

## NOT RUN:
## data(package = "pasadr")

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
                  
# check the pasad model objects
str(obj)
                  
}
\keyword{pasad}
\keyword{trajectory matrix}
