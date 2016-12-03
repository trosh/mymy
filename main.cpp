#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

//ici je vais refaire le projectmymy file mais avec un signe moins devant le calcul
int main()
{

    int N=42, M=500;
    int i,j,k;
    double Cv=45;//PLus on augmente Cv plus le transfèrt thermique est difficile
    double Tfin=0.5*M*pow(1./N,2) ;
    double dx=1.0/N, dy=1/N ;
    double x[N+1];
    double y[N+1];
    float temps[M];//sert à rien
    double dt=Tfin/M;


    //cout << Tfin << " " << M << " " << dt << endl;
    for (i = 0 ; i < N+1 ; i++)
    {
        x[i]=i*dx;
        y[i]=i*dx;
    }

    //for (i = 0 ; i< M ; i++)
    //    temps[i]=i/M;


    //Check de : condition de stabilité
    if (dt > 0.5 * pow(dx,2))
    {
        printf("bad dt");
    }

    //Deuxième partie : comment calculer T(xi,yj) avec i€[0,N] et j€[0,N] T(x,y)

    //Initialisation :

    double T[N+1][N+1], T0[N+1][N+1];

    for (i=0 ; i<N+1 ; i++)
    {
        for (j=0 ; j<N+1 ; j++)
        {
          T[i][j]=0; //condition initiale commune
        }
    }

    for (i=N/3 ; i<2*N/3 ; i++)
    {
        for (j=N/3 ; j<2*N/3 ; j++)
        {
            T[i][j]=3; //condition initiale centrale (zone centrale plus chaude mais elle peut varier)
        }
    }


    //Calcul de T:
    float A=0, B=0; //j'ai enlevé la partie nulle pour alléger l'équation :  (-A)*(T[i+1][j]-T[i-1][j]) /dx-B*( T[i][j+1]-T[i][j-1] )/dy-
    float C=2;
    int l;
    float inter, blabla, hihi;
    //A=dk/dx B=dk/dy C=K(x,y)//On fait la supposition que K homogène et constant
    l=0;//Compteur
    ofstream mymy("temperature.gnu");
    double Td=250;
    double *pointeurT=NULL;

    for (j=0 ; j<N+1 ; j++)
    {
        T[0][j]=0;
        T[j][0]=0;
        T[N][j]=0;
        T[j][N]=0;
    }

    for  (l=0  ; l<10000; l++)
    {
        for (i=1 ; i<N ; i++)//avant, c'était i=0 et i<N+1
        {
            for (j=1 ; j<N ; j++)
            {
                pointeurT=&T[i][j];
                if (T[i][j]>Td) {Cv=24;}//pour sortir de dlg il faut que >215
                if (T[i][j]<=Td) {Cv=(T[i][j]);}//avant cetait 24
                hihi=-C*(T[i+1][j]-4*T[i][j]+T[i-1][j]+T[i][j+1]+T[i][j-1]);
                T[i][j]= -(dt/Cv)*( hihi/pow(dx,2) ) +T[i][j] ;

                mymy << l*dt << " " << x[i] << " " << y[j] << " " << T[i][j] << endl;
            }
            mymy  << endl ;
        }
        mymy  << endl ;
    }


    //set term post
    //set out '|lpr'
    //replot
    //q

}

