#include <string>
#include <vector>
#include <algorithm>

using namespace std;



// group alpha and lambda together so that sorting lambda sorts alpha
typedef struct s_alphalambda {
	double alpha;
	double lambda;
} alphalambda;

// struct for keeping track of the beta matrix better
typedef struct s_betarow {
	double alpha;
	double lambda;
	int row;
	string col;
} betarow;

// for sorting
bool lambdaLessThan (alphalambda left, alphalambda right)  {
	return (left.lambda < right.lambda);
}


void read_cp(double *alpha, double *lam,  int *row, char **col, double gamma, int length, int k, double *beta){
	vector<betarow> matrix;

	for (int i = 0; i < length; i++)  {
		betarow b;
		b.alpha = alpha[i];
		b.lambda = lam[i];
		b.row = row[i];
		b.col = string (col[i]);

		matrix.push_back (b);
	}

	int current_item=0;
	int betaCounter = 0;
	string current_col="";
	vector< vector<alphalambda> > accumulated;
	for (int i = 0; i < k; i++)
		accumulated.push_back (vector<alphalambda> ());

	for (int i = 0; i < length; i++){
		if(matrix[i].col != current_col){
			if(current_col != ""){
				// for each of the k lists, sort and calculate
				for (int j = 0; j < k; j++)  {
					sort (accumulated[j].begin (), accumulated[j].end (), lambdaLessThan);
					int upperBoundIndex = 1;
					while (gamma > accumulated[j][upperBoundIndex].lambda)  {
						upperBoundIndex++;
						if (upperBoundIndex >= accumulated[j].size ())
							break;
					}

					double x0 = accumulated[j][upperBoundIndex-1].lambda;
					double x1 = (upperBoundIndex < accumulated[j].size () ? accumulated[j][upperBoundIndex].lambda : 0.0);
					double y0 = accumulated[j][upperBoundIndex-1].alpha;
					double y1 = (upperBoundIndex < accumulated[j].size () ? accumulated[j][upperBoundIndex].alpha : 0.0);

					double ygamma = y0;
					if (upperBoundIndex < accumulated[j].size ())
						ygamma = y0 + (y1 - y0) * (gamma - x0) / (x1 - x0);

					beta[betaCounter] = ygamma;
					betaCounter++;
				}
			}

			// then clear
			for(int j=0; j<k; j++)
				accumulated[j].clear ();
		}
		//for each row
		//
		// the row is row - 1 so that a value of 1 goes to index 0, etc.
		int myrow=matrix[i].row - 1;
		alphalambda al;
		al.alpha = matrix[i].alpha;
		al.lambda = matrix[i].lambda;
		accumulated[myrow].push_back (al);

		current_col = matrix[i].col;
		current_item++;
	}

	// sort and calculate final batch
	for (int j = 0; j < k; j++)  {
		sort (accumulated[j].begin (), accumulated[j].end (), lambdaLessThan);
		int upperBoundIndex = 1;
		while (gamma > accumulated[j][upperBoundIndex].lambda)  {
			upperBoundIndex++;
			if (upperBoundIndex >= accumulated[j].size ())
				break;
		}

		double x0 = accumulated[j][upperBoundIndex-1].lambda;
		double x1 = (upperBoundIndex < accumulated[j].size () ? accumulated[j][upperBoundIndex].lambda : 0.0);
		double y0 = accumulated[j][upperBoundIndex-1].alpha;
		double y1 = (upperBoundIndex < accumulated[j].size () ? accumulated[j][upperBoundIndex].alpha : 0.0);

		double ygamma = y0;
		if (upperBoundIndex < accumulated[j].size ())
			ygamma = y0 + (y1 - y0) * (gamma - x0) / (x1 - x0);

		beta[betaCounter] = ygamma;
		betaCounter++;
	}
}
//
//
// R wrappers
extern "C" {
  void read_cp_R(double *alpha, double *lam,  int *row, char **col, double *gamma, int *length, int *k, double *beta) {
       read_cp(alpha, lam, row, col, *gamma, *length, *k, beta);
      }
  }
