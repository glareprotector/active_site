   #include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <time.h>
#include <sys/time.h>
#include "globals.h"
#include <fstream>

using namespace std;



void Error (const char *error_msg);

int _ASSERT_FAILED (char *filename, int line_number, const char *error_msg);

#ifdef NDEBUG
#define Assert(test,error_msg)
#else
#define Assert(test,error_msg) (test ? 0 : _ASSERT_FAILED(__FILE__, __LINE__, error_msg))
#endif

bool toggle_error = false;

int _ASSERT_FAILED (char *filename, int line_number, const char *error_msg){
  if (toggle_error) return 0;
  cerr << "Assertion failed in file \"" << filename << "\", line " << line_number << ": " << error_msg << endl;
  abort();
  return 0;
}

void Error (const char *error_msg){
  cerr << "ERROR: " << error_msg << endl;
  toggle_error = true;
  //exit (1);
}

double rt = 1e-8;

class Minimizer {

 protected:

  void register_pys(PyObject* pMaker, PyObject* pParams, bool recalculate){
    pMaker_cur = pMaker;
    pParams_cur = pParams;
    recalculate_cur = recalculate;
  }
  
  void unregister_pys(){
    pMaker_cur = NULL;
    pParams_cur = NULL;
    recalculate_cur = -1;
  }
  
  PyObject* get_pMaker(){
    assert(pMaker_cur != NULL);
    return pMaker_cur;
  }

  PyObject* get_pParams(){
    assert(pParams_cur != NULL);
    return pParams_cur;
  }


  bool get_recalculate(){
    assert(recalculate_cur != -1);
    return recalculate_cur;
  }

  PyObject* pMaker_cur;
  PyObject* pParams_cur;
  int recalculate_cur;

 private:

  bool use_preconditioner;
  int N;
  int M;			   
  double term_ratio;
  double alpha;
  double beta;
  double gamma1;
  double gamma2;
  double min_improvement_ratio;
  int max_line_search_evaluations;

  virtual void Report (const vector<double> &theta, int iteration, double objective, double step_length, int which_infer) = 0;
  virtual void Report (const string &s) = 0;
  //virtual void ComputeGradient (vector<double> &g, const vector<double> &x, bool bCalculateGate) = 0;
  virtual void ComputeGradient (vector<double> &g, const vector<double> &x, int which_obj, int which_reg) = 0;
  virtual void ComputeHessianDiagonal (vector<double> &h, const vector<double> &x){ h = x; }
  virtual double ComputeFunction  (const vector<double> &x, int which_obj, int which_reg) = 0;

  double LineSearch (const vector<double> &x, const vector<double> &d, const vector<double> &g, double &f_curr, int which_obj, int which_reg);

public:
  Minimizer(bool use_preconditioner,
	    const int    M                           = 20,       /* number of previous gradients to remember                          */
	    const double term_ratio                  = 0.000001,  /* required ratio of gradient norm to parameter norm for termination */
	    const double alpha                       = 0.000001,    /* minimum improvement ratio for sufficient decrease                 */
	    const double beta                        = 0.5,      /* default step size                                                 */
	    const double gamma1                      = 0.01,     /* maximum step size                                                 */
	    const double gamma2                      = 0.8,      /* minimum step size                                                 */
	    const double min_improvement_ratio       = 1.000001,     /* minimum improvement ratio after sufficient decrease               */
	    const int    max_line_search_evaluations = 10);      /* maximum number of line search function evaluations                */
  
  virtual ~Minimizer(){}
  
  arbi_array<num1d> LBFGS (PyObject* pMaker, PyObject* pParams, bool recalculate, vector<double> &x0,                      /* initial guess of solution    */
			   const int max_iterations, int which_obj, int which_reg, int which_infer);         /* maximum number of iterations */
  void ApproximateGradient (vector<double> &g, const vector<double> &x, const double EPSILON, int which_obj, int which_reg);
};

/* Standard linear algebra */
double DotProduct (const vector<double> &x, const vector<double> &y);
double Norm (const vector<double> &x);
const vector<double> operator/(const vector<double> &x, const vector<double> &y);
const vector<double> operator+(const vector<double> &x, const vector<double> &y);
const vector<double> operator-(const vector<double> &x, const vector<double> &y);
const vector<double> &operator+= (vector<double> &x, double c);
const vector<double> operator*(const vector<double> &x, double c);

/* Constructor */
Minimizer::Minimizer(bool use_preconditioner, const int M, const double term_ratio,
		     const double alpha, const double beta, const double gamma1, const double gamma2, 
			 const double min_improvement_ratio, const int max_line_search_evaluations) :
	use_preconditioner (use_preconditioner), M (M), term_ratio (term_ratio),
	alpha (alpha), beta (beta), gamma1 (gamma1), gamma2 (gamma2), 
	min_improvement_ratio (min_improvement_ratio), 
	max_line_search_evaluations (max_line_search_evaluations) 
{
}
  
/* Modified cubic backtracking line search */
double Minimizer::LineSearch (const vector<double> &x, const vector<double> &d, 
			      const vector<double> &g, double &f_curr, int which_obj, int which_reg)
{
  int num_evaluations = 0;
  const double dot_prod = DotProduct (d, g);
  bool increasing_step = false;
  bool sufficient_decrease = false;
  double best_t = 0;
  double best_f = f_curr;

  /* First, try a full Newton step. */
  double t_new1 = 1;
  double f_new1 = ComputeFunction (x + d * t_new1, which_obj, which_reg);
  if (f_new1 < best_f){ best_f = f_new1; best_t = t_new1; }
  if (f_new1 <= f_curr + alpha * t_new1 * dot_prod) sufficient_decrease = true;

  /* If a sufficient decrease is found, then we'll allow the multiplier to get larger.
   * Otherwise, we'll force the multiplier to get gradually smaller. */
  if (sufficient_decrease) increasing_step = true;
    
  /* Now perform quadratic interpolation of:
   *  
   *    f_curr + dot_prod * t + ((f_new1 - f_curr) / t_new1^2 - dot_prod / t_new1) * t^2
   *
   * Note that this function is equal to f_curr at t = 0 and f_new1 at t = t_new1.  Note
   * also that the quadratic fit works only if the coefficient of t^2 is positive; this
   * will always be the case when a sufficient decrease has not been found since
   *
   *    f_new1 > f_curr + alpha * t_new1 * dot_prod 
   *
   * implies
   *
   *    f_new1 > f_curr + t_new1 * dot_prod
   */
  double t_new2 = t_new1;
  double f_new2 = f_new1;
  
  t_new1 = -dot_prod / (2 * ((f_new2 - f_curr) / t_new2 - dot_prod) / t_new2);
  
  /* Check to make sure the minimization of the quadratic was valid.  If not, try scaling the
   * t value instead. */
  if (f_new2 <= f_curr + dot_prod * t_new2){
    if (increasing_step) t_new1 = t_new2 / beta;
    else t_new1 = t_new2 * beta;  /* This case is not really necessary, as explained above. */
  }

  /* If we're doing a decreasing step, clip the prediction to a restricted range. */
  if (!increasing_step) t_new1 = max (gamma1 * t_new2, min (gamma2 * t_new2, t_new1));
  
  /* Compute the new function value, check for sufficient decrease, and check for termination. */
  f_new1 = ComputeFunction (x + d * t_new1, which_obj, which_reg);
  if (f_new1 < best_f){ best_f = f_new1; best_t = t_new1; }
  if (f_new1 <= f_curr + alpha * t_new1 * dot_prod) sufficient_decrease = true;
  if (sufficient_decrease && f_new1 >= f_new2 * min_improvement_ratio){ f_curr = best_f; return best_t; }

  while (true){

    /* Now perform cubic interpolation of
     *
     *    f_curr + dot_prod * t + b * t^2 + a * t^3
     */
    
    double a = 1 / (t_new1 - t_new2) * 
      ((f_new1 - f_curr - dot_prod * t_new1) / (t_new1 * t_new1) -
       (f_new2 - f_curr - dot_prod * t_new2) / (t_new2 * t_new2));
    
    double b = 1 / (t_new1 - t_new2) * 
      (-(f_new1 - f_curr - dot_prod * t_new1) * t_new2 / (t_new1 * t_new1) +
       (f_new2 - f_curr - dot_prod * t_new2) * t_new1 / (t_new2 * t_new2));

    t_new2 = t_new1;
    f_new2 = f_new1;
    
    t_new1 = (-b + sqrt(b*b - 3*a*dot_prod)) / (3*a);

    /* Check to make sure the minimization of the cubic was valid.  If not, try scaling the
     * t value instead. */
    if (b*b - 3*a*dot_prod <= 0){
      if (increasing_step) t_new1 = t_new2 / beta;
      else t_new1 = t_new2 * beta;
    }

    /* If we're doing a decreasing step, clip the prediction to a restricted range. */
    if (!increasing_step) t_new1 = max (gamma1 * t_new2, min (gamma2 * t_new2, t_new1));

    /* Compute the new function value, check for sufficient decrease, and check for termination. */
    f_new1 = ComputeFunction (x + d * t_new1, which_obj, which_reg);
    if (f_new1 < best_f){ best_f = f_new1; best_t = t_new1; }
    if (f_new1 <= f_curr + alpha * t_new1 * dot_prod) sufficient_decrease = true;
    if (sufficient_decrease && f_new1 * min_improvement_ratio >= f_new2){ f_curr = best_f; return best_t; }
    
    if (++num_evaluations >= max_line_search_evaluations){ f_curr = best_f; return best_t; }
  }
}

/* Finite-difference--based gradient */
//void Minimizer::ApproximateGradient (vector<double> &g, 
void Minimizer::ApproximateGradient (vector<double> &og, 
				     const vector<double> &x, 
				     const double EPSILON, int which_obj, int which_reg){

  /*
  ofstream myfile;
  myfile.open("theta_asdf");
  for(int i = 0; i < x.size(); i++){
    myfile<<x[i]<<',';
    }*/



  double base = ComputeFunction (x, which_obj, which_reg);
  vector<double> x_copy = x;
  for(int i = 0; i < x_copy.size(); i++){
    //cout<<x[i]<<" ";

  }
  for (int i = 0; i < 103; i++){
    x_copy[i] += EPSILON;
    //g[i] = (ComputeFunction (x_copy) - base) / EPSILON;
    double f = ComputeFunction (x_copy, which_obj, which_reg);
	double g = (f - base) / EPSILON;
	x_copy[i] = x[i];
	//og[i]=g;
	if (proc_id==0)
		cerr << i << ": " << og[i] <<"::" << g << " " << (og[i]-g)/g <<"; f=" << f <<"; base=" << base <<endl;
  }
  //exit(1);
}

/* LBFGS routine */
arbi_array<num1d> Minimizer::LBFGS (PyObject* pMaker, PyObject* pParams, bool recalculate, vector<double> &x0, const int max_iterations, int which_obj, int which_reg, int which_infer){


  register_pys(pMaker, pParams, recalculate);



	time_t start, end;
	time(&start);

	/* Initialization */
	ostringstream oss;
	int N = x0.size();
	vector<vector<double> > x (2, vector<double>(N, 0.0));   /* iterates                    */
	vector<vector<double> > g (M, vector<double>(N, 0.0));   /* gradients                   */
	vector<vector<double> > y (M, vector<double>(N, 0.0));   /* y[k] = g[k+1] - g[k]        */
	vector<vector<double> > s (M, vector<double>(N, 0.0));   /* s[k] = x[k+1] - x[k]        */
	vector<double> d (N, 0.0);                               /* d[k] = -H[k] g[k]           */
	vector<double> rho (M, 0.0);                             /* rho[k] = 1 / (y[k]^T s[k])  */
	vector<double> a (M, 0.0);
	vector<double> b (M, 0.0);
	vector<double> h (N, 0.0);                               /* hessian diagonal            */

	double f_prev = 0;
	double f_curr = ComputeFunction (x0, which_obj, which_reg);

if (proc_id==0)
{
	time(&end);
	double dif = difftime (end,start);
	cerr << "ComputeFunction (" << f_curr << ") " << dif << " seconds." << endl;
	time(&start);
}

	int k = 0, iterations = 0;
	x[0] = x0;

	int num_consec_small_steps = 0;
	bool progress_made = true;
	//bool bCalculateGate=true;
	Assert (x0.size() != 0, "Empty initial vector specified.");

	while (true){

		/* STEP ONE: Compute new gradient vector */
		//ComputeGradient (g[k%M], x[k%2], bCalculateGate);
	  ComputeGradient (g[k%M], x[k%2], which_obj, which_reg);

if (proc_id==0)
{
	time(&end);
	double dif = difftime (end,start);
	cerr << "iteration-" << iterations << ": ComputeGradient " << dif << " seconds. " << k << endl;
	time(&start);
}

	
if (iterations%1000  == 100)
//if(true)
{
  

	for (int i = 7; i < 8; i++)
	  //ApproximateGradient (aa[i], x[k%2], pow(10.0,(double)-i));
	  ApproximateGradient (g[k%M], x[k%2], pow(10.0,(double)-i), which_obj, which_reg);

	//for (int i = 0; i < 200; i++){//g[k%M].size(); i++){
	//  cerr << g[k%M][i];
	//  for (int j = 4; j < 5; j++) cerr << "   " << aa[j][i];
	//  cerr << endl;
	//}
	//cerr << ComputeFunction (x0) << endl;

	//if (iterations==3)exit(0);
}
	

	if (use_preconditioner) ComputeHessianDiagonal (h, x[k%2]);

		if (k > 0) y[(k-1+M)%M] = g[k%M] - g[(k-1+M)%M];
		if (k > 0) s[(k-1+M)%M] = x[k%2] - x[(k-1+2)%2];

		if (k > 0) rho[(k-1+M)%M] = 1.0 / DotProduct (y[(k-1+M)%M], s[(k-1+M)%M]);
		d = g[k%M];
		if ((k > 0 && DotProduct (y[(k-1+M)%M], s[(k-1+M)%M]) <= 0) || k >= 1000){
		  /* Delete the old gradient info */
		  g[0] = g[k%M];
		  x[0] = x[k%2];
		  k = 0;
		}
		//else 
		{
			/* Use Nocedal's recursion to compute H_k g_k */
			for (int j = k-1; j >= max(0,k-M); j--){
				a[j%M] = DotProduct (s[j%M], d) * rho[j%M];
				d = d - y[j%M] * a[j%M];
			}
		  
			/* Apply preconditioner (inverse Hessian diagonal) */

			if (use_preconditioner){
				d = d / h;
			} else {
				d = d * (1.0/Norm(g[k%M]));
			}
		  
			/* Continue using recursion formula */
			for (int j = max(0,k-M); j <= k-1; j++){
				b[j%M] = DotProduct (y[j%M], d) * rho[j%M];
				d = d + s[j%M] * (a[j%M] - b[j%M]);
			}
		}
    		//if (proc_id==0){cerr << "&& d:";for (int i=0; i<15; i++) cerr << d[i]<<","; cerr << endl;}
		d = d * -1.0;
		f_prev = f_curr;
		double step = LineSearch (x[k%2], d, g[k%M], f_curr, which_obj, which_reg);
		x[(k+1)%2] = x[k%2] + d * step;

		iterations++;
		k++;

if (proc_id==0)
{
	time(&end);
	double dif = difftime (end,start);
	cerr << "LineSearch " << dif << " seconds.(" << f_curr << ")" << endl;
	time(&start);
}

                Report (x[k%2], iterations, f_curr, step, which_infer);
		//bCalculateGate = false;
		if (iterations >= max_iterations){ 
		  oss.str("");
		  oss << "Termination: maximum number of iterations reached"; 
		  Report (oss.str());
		  break; 
		}

		if (f_curr == 0){
		  oss.str("");
		  oss << "Termination: Zero reached.";
		  Report (oss.str());
		  break;
		}

		if (fabs(f_prev - f_curr) / f_prev < rt) num_consec_small_steps++;
		else {
			num_consec_small_steps = 0;
			progress_made = true;
		}
		progress_made = true;  

		if (num_consec_small_steps == 9){
			if (progress_made){
				progress_made = false;
				num_consec_small_steps = 0;
				oss.str("");
				oss << "Restart: Too many consecutive small steps";
				Report (oss.str());
				g[0] = g[k%M];
				x[0] = x[k%2];
				rt = max(rt/2,1e-8);
				k = 0;
			}
			else {
				oss.str("");
				oss << "Termination: Too many consecutive small steps";
				Report (oss.str());
				break;
			}
		}
if (proc_id==0)
{
	time(&end);
	double dif = difftime (end,start);
	cerr << "STEP SIX: Check termination conditions " << dif << " seconds." << endl;
	time(&start);
}

	}

	x0 = x[k%2];



// copy x into an arbi_array
int x0size = x0.size();
arbi_array<num1d> ans(x0size);
for(int i = 0; i < x0.size(); i++){
  ans(i) = x0[i];
 }

unregister_pys();

return ans;

}

/* Standard linear algebra */
double DotProduct (const vector<double> &x, const vector<double> &y){
  Assert (x.size() == y.size(), "Vector size mismatch.");
  double ret = 0;
  for (int i = 0; i < (int) x.size(); i++) ret += x[i] * y[i];
  return ret;
}

double Norm (const vector<double> &x){
  return sqrt(DotProduct (x,x));
}

const vector<double> operator/(const vector<double> &x, const vector<double> &y){
  vector<double> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] /= y[i];
  return ret;
}

const vector<double> operator+(const vector<double> &x, const vector<double> &y){
  vector<double> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] += y[i];
  return ret;
}

const vector<double> operator-(const vector<double> &x, const vector<double> &y){
  vector<double> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] -= y[i];
  return ret;
}

const vector<double> &operator+= (vector<double> &x, double c){
  for (int i = 0; i < (int) x.size(); i++) x[i] += c;
  return x;
}

const vector<double> operator*(const vector<double> &x, double c){
  vector<double> ret (x);
  for (int i = 0; i < (int) ret.size(); i++) ret[i] *= c;
  return ret;
}
