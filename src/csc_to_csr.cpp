#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
Rcpp::List csc_tocsr(
    const int n_row,
    const int n_col, 
    const std::vector<int> &Ap, 
    const std::vector<int> &Ai, 
    const std::vector<double> &Ax
) {
  // source code in-parts taken from 
  // https://github.com/scipy/scipy/blob/3b36a574dc657d1ca116f6e230be694f3de31afc/scipy/sparse/sparsetools/csr.h#L379-L424
  const int nnz = Ap[n_col];
  std::vector<int> Bp(n_row+1), Bj(nnz);
  std::vector<double> Bx(nnz);
  for (int n = 0; n < nnz; n++){            
    Bp[Ai[n]]++;
  }
  for(int row = 0, cumsum = 0; row < n_row; row++){     
    int temp  = Bp[row];
    Bp[row] = cumsum;
    cumsum += temp;
  }
  Bp[n_row] = nnz; 
  for(int col = 0; col < n_col; col++){
    for(int jj = Ap[col]; jj < Ap[col+1]; jj++){
      int row  = Ai[jj];
      int dest = Bp[row];
      Bj[dest] = col;
      Bx[dest] = Ax[jj];
      Bp[row]++;
    }
  } 
  for(int row = 0, last = 0; row <= n_row; row++){
    int temp  = Bp[row];
    Bp[row] = last;
    last    = temp;
  }
  
  Rcpp::List L = Rcpp::List::create(Rcpp::Named("p") = Bp , Rcpp::Named("i") = Bj, Rcpp::Named("x") = Bx);
  return L;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
library(SeuratObject)
data("pbmc_small")

m <- GetAssayData(pbmc_small, "counts")
m2.list <- csc_tocsr(m@Dim[[1]], m@Dim[[2]], m@p, m@i, m@x)
m2 <- new(
  Class = 'dgRMatrix',
  p = as.integer(m2.list$p),
  j = as.integer(m2.list$i),
  x = m2.list$x,
  Dim = slot(object = m, name = 'Dim')
)

mr <- as(object = as.matrix(x = m), Class = 'dgRMatrix')

all.equal(target = slot(object = m2, name = 'p'), current = slot(object = mr, name = 'p'))
all.equal(target = slot(object = m2, name = 'j'), current = slot(object = mr, name = 'j'))
all.equal(target = slot(object = m2, name = 'x'), current = slot(object = mr, name = 'x'))
all.equal(target = slot(object = m2, name = 'Dim'), current = slot(object = mr, name = 'Dim'))
*/