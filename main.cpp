#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Mat : public vector<double>
{
  public:
    const int nrows, ncols;
    Mat(int nrows, int ncols)
    : nrows(nrows), ncols(ncols)
    , vector<double>(nrows * ncols)
    { }
    Mat(int nrows, int ncols, double initval)
    : nrows(nrows), ncols(ncols)
    , vector<double>(nrows * ncols, initval)
    { }
    double & operator()(int i, int j)
    {
        return operator[](i * ncols + j);
    }
    double operator()(int i, int j) const
    {
        return operator[](i * ncols + j);
    }
};

// ici je vais refaire le projectmymy file mais avec un signe moins devant le calcul
int main()
{
    const int N = 42, M = 500;
    // PLus on augmente Cv plus le transfèrt thermique est difficile
    double Cv = 45;
    const double Tfin = 0.5 * M / (N * N);
    const double dx = 1.0 / N;
    const double dy = 1.0 / N;
    const double dt = Tfin / M;
    vector<double> x(N+1), y(N+1);

    //cout << Tfin << " " << M << " " << dt << endl;
    for (int i=0; i<N+1; i++)
    {
        x[i] = i*dx;
        y[i] = i*dx;
    }

    // Check de : condition de stabilité
    if (dt > 0.5 * dx * dx)
    {
        cerr << "bad dt" << endl;
        cerr << "(ignoring)" << endl;
    }

    // Deuxième partie : comment calculer T(xi,yj) avec i€[0,N] et j€[0,N] T(x,y)

    // Initialisation :

    Mat T (N+1, N+1, 0.0);
    Mat T0(N+1, N+1, 0.0);
    //vector<vector<double> > T (N+1, vector<double>(N+1, 0.0));
    //vector<vector<double> > T0(N+1, vector<double>(N+1, 0.0));

    for (int i=N/3; i<2*N/3; i++)
    for (int j=N/3; j<2*N/3; j++)
    {
        // condition initiale centrale (zone centrale plus chaude mais elle peut varier)
        T(i, j) = 3;
    }

    // Calcul de T
    // j'ai enlevé la partie nulle pour alléger l'équation :
    // (-A)*(T[i+1][j]-T[i-1][j]) /dx-B*( T[i][j+1]-T[i][j-1] )/dy-
    float A = 0;
    float B = 0;
    float C = 2;
    //A=dk/dx B=dk/dy C=K(x,y) //On fait la supposition que K homogène et constant
    ofstream mymy("temperature.gnu");
    double Td = 250;

    for (int j=0; j<N+1; j++)
    {
        T(0, j) = 0.0;
        T(j, 0) = 0.0;
        T(N, j) = 0.0;
        T(j, N) = 0.0;
    }

    for (int l=0; l<500; l++)
    {
        //for (int i=0; i<N+1; i++)
        for (int i=1; i<N; i++)
        {
            for (int j=1; j<N; j++)
            {
                if (T(i, j) > Td)
                {
                    Cv = 24;
                } // pour sortir de dlg il faut que >i 215
                else
                {
                    Cv = T(i, j);
                } // avant cetait 24

                const float hihi =
                    -C * (T(i+1, j  )
                        - T(i  , j  ) * 4
                        + T(i-1, j  )
                        + T(i  , j+1)
                        + T(i  , j-1));

                T(i, j) += -dt*hihi / (Cv*dx*dx);

                mymy << l*dt << " "
                     << x[i] << " "
                     << y[j] << " "
                     << T(i, j) << endl;
            }
            mymy << endl ;
        }
        mymy << endl ;
    }

    //set term post
    //set out '|lpr'
    //replot
    //q

    return 0;

}

