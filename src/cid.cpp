// CIDnetworks C++ code snippets.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix makeEdgeList (int nn) {
  int rows = nn*(nn-1)/2;
  IntegerMatrix out(rows,2);
  int ii,jj,cc;
  cc=0;
  for (ii=1; ii<nn; ++ii) for (jj=ii+1; jj<nn+1; ++jj) {out(cc,0)=ii; out(cc,1)=jj; ++cc;}
  return(out);
}

// [[Rcpp::export]]
IntegerMatrix makeEdgeListSelfies (int nn) {
  int rows = nn*(nn+1)/2;
  IntegerMatrix out(rows,2);
  int ii,jj,cc;
  cc=0;
  for (ii=1; ii<nn+1; ++ii) for (jj=ii; jj<nn+1; ++jj) {out(cc,0)=ii; out(cc,1)=jj; ++cc;}
  return(out);
}

// [[Rcpp::export]]
NumericMatrix symBlock (NumericVector entries) {
	int dd = (sqrt(1+8*entries.size())-1)/2;
	NumericMatrix output (dd,dd);
	IntegerMatrix els = makeEdgeListSelfies (dd);
	int ii;
	for (ii=0; ii<els.rows(); ++ii) {
		output(els(ii,0)-1, els(ii,1)-1) = entries[ii];
		output(els(ii,1)-1, els(ii,0)-1) = entries[ii];
	}
	return(output);
}

// [[Rcpp::export]]
NumericVector eldc (NumericMatrix latentSpacePos, IntegerMatrix edgelist) {
  int rr, cc;
  NumericVector out(edgelist.rows());
  for (rr=0; rr<edgelist.rows(); ++rr) {
    out[rr] = 0;
    for (cc=0; cc < latentSpacePos.cols(); ++cc) {
      out[rr] += pow(latentSpacePos(edgelist(rr,0)-1,cc) - latentSpacePos(edgelist(rr,1)-1,cc), 2);
    }
    out[rr] = sqrt(out[rr]); 
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector cosineClosenessC (NumericMatrix latentSpacePos, IntegerMatrix edgelist) {
  int rr, cc;
  NumericVector out(edgelist.rows());
  for (rr=0; rr<edgelist.rows(); ++rr) {
    out[rr] = 0;
    for (cc=0; cc < latentSpacePos.cols(); ++cc) {
      out[rr] += latentSpacePos(edgelist(rr,0)-1,cc)*latentSpacePos(edgelist(rr,1)-1,cc);
    }
    //out[rr] = sqrt(out[rr]); 
  }
  return(out);
}



// [[Rcpp::export]]
IntegerMatrix coincidence (IntegerMatrix edgelist, int nodes) {

	int ii,jj,kk;
	IntegerMatrix out(nodes, nodes);
	for (jj=0; jj<nodes; ++jj) 	for (kk=0; kk<nodes; ++kk) out(jj,kk) = 0;

	for (ii=0; ii<edgelist.rows(); ++ii) {
		out(edgelist(ii,0)-1, edgelist(ii,1)-1)++;   
		out(edgelist(ii,1)-1, edgelist(ii,0)-1)++;   
		out(edgelist(ii,0)-1, edgelist(ii,0)-1)++;   
		out(edgelist(ii,1)-1, edgelist(ii,1)-1)++;   
	}

	return(out);

}


// [[Rcpp::export]]
NumericVector xty (IntegerMatrix edgelist, NumericVector outcome, int nodes) {

	int ii,jj;
	NumericVector out(nodes, 0.0);
//	for (jj=0; jj<nodes; ++jj) out[jj] = 0;

	for (ii=0; ii<edgelist.rows(); ++ii) {
		out[edgelist(ii,0)-1] += outcome[ii];   
		out[edgelist(ii,1)-1] += outcome[ii];   
	}

	return(out);

}




