#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

double difr(double r, double z, double(*U)(double,double), double h1)
{
    return (U(r+h1,z) - U(r-h1, z))/(2*h1);
}

double difz(double r, double z, double(*U)(double,double), double h2)
{
    return (U(r,z+h2) - U(r, z-h2))/(2*h2);
}

vector <vector < vector <double> > > & lambda1(vector <vector < vector <double> > > const & old_vector, double h1)
{
    vector < vector < vector <double> > > * new_vector =
            new vector < vector < vector <double> > > (old_vector);

    for(unsigned i = 0; i < old_vector.size(); ++i)
        for(unsigned j = 1; j < old_vector[i].size() - 1; ++j)
            for(unsigned k = 1; k < old_vector[i][j].size() - 1; ++k)
                (*new_vector)[i][j][k] = (old_vector[i][j + 1][k] - old_vector[i][j - 1][k])/(h1 * 2);
    for(unsigned j = 0; j < old_vector[i].size(); ++j) {
        (*new_vector)[0][j][0] = 0;
        (*new_vector)[1][j][0] = p[j][0];
        (*new_vector)[2][j][0] = 0;
        (*new_vector)[0][j][old_vector[i][j].size() - 1] = 0;
        (*new_vector)[1][j][old_vector[i][j].size() - 1] = 0;
        (*new_vector)[2][j][old_vector[i][j].size() - 1] = 0;
    }
    for(unsigned j = 0; j < old_vector[i][0].size(); ++j) {
        (*new_vector)[0][0][j] = 0;
        (*new_vector)[1][0][j] = 0;
        (*new_vector)[2][0][j] = 0;
        (*new_vector)[0][old_vector[i].size() - 1][j] = 0;
        (*new_vector)[1][old_vector[i].size() - 1][j] = 0;
        (*new_vector)[2][old_vector[i].size() - 1][j] = 0;
    }

    return *new_vector;
}

vector <vector < vector <double> > > & lambda2(vector <vector < vector <double> > > const & old_vector, double h2)
{
    vector < vector < vector <double> > > * new_vector =
            new vector < vector < vector <double> > > (old_vector);

    for(unsigned i = 0; i < old_vector.size(); ++i)
        for(unsigned j = 1; j < old_vector[i].size() - 1; ++j)
            for(unsigned k = 1; k < old_vector[i][j].size() - 1; ++k)
                (*new_vector)[i][j][k] = (old_vector[i][j][k + 1] - old_vector[i][j][k - 1])/(h2 * 2);

    for(unsigned i = 0; i < old_vector.size(); ++i) {
        for(unsigned j = 0; j < old_vector[i].size(); ++j) {
            (*new_vector)[i][j][0] = 0;
            (*new_vector)[i][j][old_vector[i][j].size() - 1] = 0;
        }
        for(unsigned j = 0; j < old_vector[i][0].size(); ++j) {
            (*new_vector)[i][0][j] = 0;
            (*new_vector)[i][old_vector[i].size() - 1][j] = 0;
        }
    }

    return *new_vector;
}

double fun_vr(double r, double z)
{
    if (r >= 0 && z >= 0 && r < 1 && z < 1)
        return ;
    if (z < 0)
        return ;
    if (z > 1)
        return ;
    if (r < 0)
        return fun_vr(-r,z);
    if (r > 1)
        return fun_vr(1-r,z);
}

double fun_vz(double r, double z)
{
    if (r >= 0 && z >= 0 && r < 1 && z < 1)
        return 1-r;
    if (z < 0)
        return ;
    if (z > 1)
        return ;
    if (r < 0)
        return ;
    if (r > 1)
        return fun_vz(1-r,z);
}

int main(int argc, char *argv[])
{
    unsigned N1 = 5, N2 = 7; N3 =3;
    double h1   = 1.0/N1, h2 = 1.0/N2;
    double L    = 1;
    double H    = 1;
    double z11  = H/L;
    double z22  = H;
    double Re   = 20;

    vector <vector <double> > p     (N1, vector <double> (N2));
    vector <vector <double> > vr    (N1, vector <double> (N2));
    vector <vector <double> > vz    (N1, vector <double> (N2));
    vector <vector <double> > tau11 (N1, vector <double> (N2));
    vector <vector <double> > tau12 (N1, vector <double> (N2));
    vector <vector <double> > tau22 (N1, vector <double> (N2));
    vector <double> nu(N1);
    vector <vector <vector <double> > > W1  (3, vector <vector <double> > (N1, vector <double> (N2, 0)));
    vector <vector <vector <vector <double> > > > U (3, vector<vector<vector<double> > > (N1, vector <vector <double> > (N2, vector <double> (N3, 0))));
    vector <vector <vector <double> > > W2  (3, vector <vector <double> > (N1, vector <double> (N2, 0)));
    vector <vector <vector <double> > > W   (3, vector <vector <double> > (N1, vector <double> (N2, 0)));
    vector <vector <vector <double> > > L1W1(3, vector <vector <double> > (N1, vector <double> (N2, 0)));
    vector <vector <vector <double> > > L2W2(3, vector <vector <double> > (N1, vector <double> (N2, 0)));
    for (unsigned i = 0; i < vr.size(); ++i)
        for (unsigned j = 0; j < vr[i].size(); ++j){
            vr[i][j]    = fun_vr(i*h1, j*h2);
            vz[i][j]    = fun_vz(i*h1, j*h2);
            p[i][j]     = 1;
        }

    for (unsigned i = 0; i < nu.size(); ++i)
        nu[i] = 1/Re;

    for (unsigned i = 0; i < tau11.size(); ++i)
        for(unsigned j = 0; j < tau11[i].size(); ++j)
        {
            tau11[i][j] = 2*nu[i]*z11*difr(i*h1, j*h2, fun_vr, h1);
            tau12[i][j] = nu[i]*(z11*difr(i*h1, j*h2, fun_vz, h1) + z22*difz(i*h1, j*h2, fun_vr, h2));
            tau22[i][j] = 2*nu[i]*z22*difz(i*h1, j*h2, fun_vz, h2);
            //cout<< tau11[i][j]<< "     "<< tau12[i][j]<< "     "<< tau22[i][j]<< endl;
        }

    for (unsigned i = 0; i < W1[0].size(); ++i)
        for (unsigned j = 0; j < W1[1][i].size(); ++j){
            W1[0][i][j] = vr[i][j];
            W1[1][i][j] = (vr[i][j] * vr[i][j] - tau11[i][j] + p[i][j]);
            W1[2][i][j] = (vr[i][j] * vz[i][j] - tau12[i][j]);

            W2[0][i][j] = vz[i][j];
            W2[1][i][j] = (vr[i][j] * vz[i][j] - tau12[i][j]);
            W2[2][i][j] = (vz[i][j] * vz[i][j] - tau22[i][j] + p[i][j]);

            cout<< W1[0][i][j] <<"     "<< W1[1][i][j] <<"     "<< W1[2][i][j] << endl;
            cout<< W2[0][i][j] <<"     "<< W2[1][i][j] <<"     "<< W2[2][i][j] << endl;
        }

    L1W1 = lambda1(W1, h1);
    L2W2 = lambda2(W2, h2);

    for (unsigned i = 0; i < W1[0].size(); ++i){
        cout<<endl;
        for (unsigned j = 0; j < W1[0][i].size(); ++j){
            L1W1[0][i][j] = z11 * L1W1[0][i][j];
            L1W1[1][i][j] = z11 * L1W1[1][i][j];
            L1W1[2][i][j] = z11 * L1W1[2][i][j];

            L2W2[0][i][j] = z22 * L2W2[0][i][j];
            L2W2[1][i][j] = z22 * L2W2[1][i][j];
            L2W2[2][i][j] = z22 * L2W2[2][i][j];

            W[0][i][j]    = L1W1[0][i][j] + L2W2[0][i][j];
            W[1][i][j]    = L1W1[1][i][j] + L2W2[1][i][j];
            W[2][i][j]    = L1W1[2][i][j] + L2W2[2][i][j];

        }
    }
}
