#include <Rcpp.h>
using namespace Rcpp;

//' Euclidean climatch scores
//'
//' Vector of the climatch scores within the recipient region
//'
//' @param recipient A data.frame of climatic variables for the recipient region.
//' @param source A data.frame of climatic variables for the source region.
//' @param globvar A vector of the global variance of each climate variable, in the same order as the columns of source and recipient region data.frames.
//'
//' @return A vector of climatch scores corresponding to each grid cell within recipient region, i.e., each row in the recipient data.frame.
//'
//' @usage climatch_vec(recipient, source, globvar)
//'
//' @export
//' @references Crombie, J., Brown, L., Lizzio, J., & Hood, G. (2008). "Climatch user manual"
//'
//' @examples
//' i <- as.data.frame(matrix(runif(n=180, min=1, max=20), nrow=60)) # Fake source climate data
//' j <- as.data.frame(matrix(runif(n=300, min=10, max=40), nrow=100)) # Fake recipient data
//' variance <- c(600, 800, 450) # Fake global variance
//'
//' climatch_vec(recipient = j, source = i, globvar = variance)
//'
// [[Rcpp::export]]
NumericVector climatch_vec(DataFrame recipient, DataFrame source, NumericVector globvar){
  // Convert dataframes to matrices
  NumericMatrix jarea = internal::convert_using_rfunction(recipient, "as.matrix");
  NumericMatrix iarea = internal::convert_using_rfunction(source, "as.matrix");

  double var = globvar.size(); // Number of variables being matched
  int mn = source.nrow();      // Number of rows in source region
  int jn = recipient.nrow();   // Number of rows in recipient region

  NumericVector jmax(jn);

  // The climatch algorithm
  for(int j=0; j<jn; j++){
    double imin = 0;
    NumericVector biosum(mn);
    NumericVector ja = jarea(j,_); // Select the row of the recipient region

    // For each recipient row run the match to every source row
    for(int i=0; i<mn; i++){
      double biovar = 0.0;
      NumericVector ia = iarea(i,_); // Select the row of the source region

      // Computes the match between the given pair of source and recipient rows
      for(int m=0; m<var; m++){
        biovar = (pow((ia[m] - ja[m]),2)/globvar[m]) + biovar; // Run the match for each variable
      }

      biosum[i] = sqrt(1/var*biovar); // Summarize score across all variables
    }

    imin = min(biosum); // select closest match
    jmax[j] = (1-imin)*10; // Convert match to 0-10 score
  }

 NumericVector max = ifelse(jmax<0, 0, jmax);

  return max;
}
