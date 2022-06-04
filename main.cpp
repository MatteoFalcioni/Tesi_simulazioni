#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 500;  //Kuramoto parameters   
double K = 50;                         

double k = 1;   //Cucker-Smale parameters
double sigma = 1; 
double beta = 8.0 ; 
double R = 5.0;

double n_c = 7.0; //Parisi parameter  (# of topological neighbours)

double t = 0.0;    //time related parameters
size_t nSteps = 500;
double dt = 0.1;   //time step
double T = 5*dt;

double L = 30.0;  //box dimension
double dx = 0.1;  //spatial step

//**********************************************************************************************************************************************
//**********************//
int model_type = 2;     ////////  choose model: 0 = global interaction. 1 = Cucker-Smale, metric interaction. 2 = Parisi, topological.   ////////  
//**********************//
//**********************************************************************************************************************************************


int main(){

    std::vector<double> x(N);   // Initial condition, vector of N elements (N ODEs)    

    std::vector<double> Int(N);  //interaction terms 

    int n = int(N);

    std::vector<double> Phases = Phases_generator(N);
    std::vector<double> Xpos = Positions_generator(L, n);
    std::vector<double> Ypos = Positions_generator(L, n);

    std::vector<double> Xplus(n);  //contiene la copia delle lucciole tra 0 ed R e tra L-R ed L
    std::vector<double> Yplus(n);

    double num = 1e5; 

    for (int i=0; i<n; i++){    //setting xi and yi to a big value instead of 0 in order to check if some rij are wrong
        Xplus[i] = num;
        Yplus[i] = num;
    }

    std::fstream positions;   //printing positions
    positions.open("Positions.txt", std::ios::out);
    for (int i=0; i<n; i++){   
        positions << Xpos[i] << '\t' << Ypos[i] << '\n';
    }
    positions.close();


    double Adj[n][n];   //Adjacency matrix
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            Adj[i][j] = 0;
        }
    }

    if (model_type == 0) {      //global interaction
        for (int i=0; i<n; ++i) {
            for (int j=0; j<n; ++j) {

                if ( i!=j ) {
                    Adj[i][j] = 1;
                }
                else if ( i == j ) {
                    Adj[i][j] = 0;
                }

            }
        }
    }

        /*
            if (model_type == 1){   //Cucker-Smale

            double avg_neighbours=0;
            for (int i=0; i<n; i++){    //filling Adjacency matrix 

                double xi = Xpos[i];
                double yi = Ypos[i];

                std::vector<double> Ri(n);

                for (int j=0; j<n; j++){

                    double xj = Xpos[j];
                    double yj = Ypos[j];

                    if( j<i ) {
                        Adj[i][j] = CS_entries(xi, xj, yi, yj, R, k, sigma, beta); 
                        Adj[j][i] = CS_entries(xi, xj, yi, yj, R, k, sigma, beta); 
                    }


                    //voglio controllare in media quante vicine vede ognuna una volta fissato R, perché almeno deve vederne n_c per comparare i due modelli, poi non voglio che ne veda una marea sennò troppo easy
                    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                    double mod_r = std::sqrt(r_ij);                 
                    Ri[j] = mod_r;
                    if ( j==i ) { Ri[j] = 0; }

                }

                double i_neighbours = 0; 
                double dr = 0.001;
                double r0 = dr*2;
                double Rmax = R;    

                for (double r=0; r<Rmax; r += r0) {    //da togliere nel for col movimento
                    for (int j=0; j<n; j++) { 
                        if ( j != i ) {
                            if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                i_neighbours += 1;
                            }
                            
                        }
                    }                
                }

                avg_neighbours += (i_neighbours / N);
            }
            std::cout << " # of average neighbours in metric interaction was " << avg_neighbours << '\n';

        } 

        else if (model_type == 2){  //Parisi
            for (int i=0; i<n; i++){
                double xi = Xpos[i];
                double yi = Ypos[i];

                std::vector<double> Ri(n);   //vector of relative distances from i

                for (int j=0; j<n; j++){
                    double xj = Xpos[j];
                    double yj = Ypos[j];

                    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                    double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                    Ri[j] = mod_r;
                    if ( j==i ) { Ri[j] = 0; }
                }
                std::cout<< "relative distances from " <<i<< " which is in ( " << Xpos[i] << " , " << Ypos[i] << " ) are: " << '\n';
                for (int j=0; j<n; j++){
                    std::cout << Ri[j] << '\t';
                }
                std::cout <<'\n';

                int i_neighbours = 0;  //conta i vicini di i
                double dr = 0.001;
                double r0 = dr*2;

                for (double r=0; r<R; r += r0) {   //NB: dr != r0 needed
                    if ( i_neighbours < n_c ){                     
                        for (int j=0; j<n; j++) { 
                            if ( j != i ) {

                                if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                    //std::cout << "for " <<i<< " " <<j<< " was seen as a close neighbour as r = " << r <<'\n';
                                    i_neighbours += 1;
                                    Adj[i][j] = 1/n_c;
                                    //if ( i_neighbours == n_c ) { std::cout<< "# of neighbours reached for " <<i<< '\n'; }
                            
                                }
                            
                            }
                        }
                    }                
                }
            }
          
        } */

/*
    double sum=0;
    for (int j=0; j < n; j++) {         //normalizing to have sum_j A_ij = 1 (Markov matrix) NB: when the fireflies will be moving, normalize just for t=0
        for (int i=0; i < n; i++) {
            sum += Adj[i][j];
        }
        for (int i=0; i < n; i++) {
            if( Adj[i][j] != 0 ) { Adj[i][j] /= sum; }
            if( isnan( Adj[i][j] ) ) { Adj[i][j] = 0; }
        }
        sum = 0;
    } 


    if ( n <= 20 ) { 
        for (int i=0; i<n; i++){        //printing Adjacency matrix
            for (int j=0; j<n; j++){
                std::cout<< Adj[i][j] << '\t';
            }
            std::cout << '\n';
        } 
    }

*/

    for (int i=0; i<N; i++){   //setting initial values for x[i]
        x[i] = Phases[i];
        Int[i] = 0;
    }

    std::fstream fout;
    fout.open("solutions.txt", std::ios::out); 
    std::fstream sync;
    sync.open("Synchronization.txt", std::ios::out); 
    std::fstream s0;
    s0.open("Initial_conditions.txt", std::ios::out);

    int counter=0;

    for ( int ii=0; ii<nSteps; ++ii ){  

        //printing m-l/N : m = # of fireflies in 0, l = # of fireflies in 1 
        double m = 0;
        double l = 0;
        for (int i=0; i<N; i++){
            if ( normalizer(x[i]) == 0 ) {
                m += 1;
            }
            if ( normalizer(x[i] ) == 1 ) {
                l += 1;
            }
        }
        if (t==0) {
            s0 << m << '\t' << l << '\t' << '\n'; 
        }
        sync << t <<'\t'<< (m-l)/N << '\n';

        fout << t << '\t';      //print solution at time t    
        for (int i=0; i<N; i++) {
            fout << x[i] << '\t';
        }
        fout << '\n';      


        //****************************interazione a t = T, t = T-dt***************************************//

        if ( ( t!=0 ) && ( t == 0 +(T-dt) || dmod(t , T , 10) == 0 || dmod(t-(T-dt) , T , 10) == 0 ) ) {     
            std::cout << "interaction at t = " <<t << '\n';

            for(int i=0;i<n;i++){
                if ( Xpos[i] > 0 && Xpos[i] < R ) {
                    double xright = Xpos[i] + L;
                    Xplus[i] = xright;
                }
                if ( Xpos[i] > L-R && Xpos[i] < L ) {
                    double xleft = Xpos[i] - L;
                    Xplus[i] = xleft;
                }
            }
            for(int i=0;i<n;i++){
                if ( Ypos[i] > 0 && Ypos[i] < R ) {
                    double yup = Ypos[i] + L;
                    Yplus[i] = yup;
                }
                if ( Ypos[i] > L-R && Ypos[i] < L ) {
                    double ydown = Ypos[i] - L;
                    Yplus[i] = ydown;
                }
            }            
            //ora se ho una lucciola che sta tra L-R ed L, ad esempio, devo andare a guardare in Xpos E in Xplus
            
            if (model_type == 1){   //Cucker-Smale

                for (int i=0; i<n; i++){    //filling Adjacency matrix

                    double xi = Xpos[i];
                    double yi = Ypos[i];

                    for (int j=0; j<n; j++){        //tutte devono controllare Xpos, Ypos

                        double xj = Xpos[j];
                        double yj = Ypos[j];
              
                        if( j<i ) {

                            double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  //modulo quadro di ri-rj
                            double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                            if ( mod_r <= R ) {
                                Adj[i][j] = ( k / ( std::pow( ((sigma*sigma) + r_ij) , beta ) ) );
                                Adj[j][i] = Adj[i][j];
                            } else if ( mod_r > R ) {
                                Adj[i][j] = 0;
                                Adj[j][i] = 0;
                            } 
                            
                        }
                    }
                    //in più, se stanno vicino alle pareti controllano anche le altre strisce
                    if ( (xi > L-R && xi < L) || (xi > 0 && xi < R) ) {
                        //std::cout << "xi, i= " <<i<< " ,was in " <<xi<< " so it should check on the other side " <<'\n';
                        for (int j=0; j<n; j++){      
                            //std::cout << "Xplus[" <<j<< "] = " <<Xplus[j] << '\n'; 
                            if (Xplus[j] != num ) {                                
                                double xj = Xplus[j];
                                double yj = Ypos[j];    //y è la stessa
              
                                if( j!=i ) {

                                    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  
                                    double mod_r = std::sqrt(r_ij); 

                                    if ( mod_r <= R ) {
                                        //std::cout << "xj was in " <<xj<<  " ;rij was " << mod_r << " while R was " <<R<< " and it was considered a neighbour of xi = " << xi <<'\n';
                                        Adj[i][j] = ( k / ( std::pow( ((sigma*sigma) + r_ij) , beta ) ) );
                                    } else if ( mod_r > R ) {
                                        //std::cout << "xj was in " <<xj<<  " ;rij was " << mod_r << " while R was " <<R<< " ,so it was not considered a neighbour of xi = " << xi <<'\n';
                                        Adj[i][j] = 0;
                                    } 
                                }
                            }
                        }

                    }
                    
                    if ( (yi > L-R && yi < L) || (yi > 0 && yi < R) ) {
                        for (int j=0; j<n; j++){      

                            if (Yplus[j] != num ) {
                                double xj = Xpos[j];        //x è la stessa
                                double yj = Yplus[j];    
              
                                if( j!=i ) {

                                    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  
                                    double mod_r = std::sqrt(r_ij); 

                                    if ( mod_r <= R ) {
                                        Adj[i][j] = ( k / ( std::pow( ((sigma*sigma) + r_ij) , beta ) ) );
                                    } else if ( mod_r > R ) {
                                        Adj[i][j] = 0;
                                    } 
                            
                                }
                            }
                        }

                    }                   


                }

            } 

            else if (model_type == 2){  //Parisi
                for (int i=0; i<n; i++){
                    double xi = Xpos[i];
                    double yi = Ypos[i];

                    std::vector<double> Ri(n);   //vector of relative distances from i

                    for (int j=0; j<n; j++){
                        double xj = Xpos[j];
                        double yj = Ypos[j];

                        double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                        double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                        Ri[j] = mod_r;
                        if ( j==i ) { Ri[j] = 0; }
                    }                    
                    
                    int i_neighbours = 0;  //conta i vicini di i
                    double dr = 0.001;
                    double r0 = dr*2;

                    for (double r=0; r<R; r += r0) {   //NB: dr != r0 needed
                        if ( i_neighbours < n_c ){                     
                            for (int j=0; j<n; j++) { 
                                if ( j != i ) {

                                    if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                        i_neighbours += 1;
                                        Adj[i][j] = 1/n_c;                            
                                    }                            
                                }
                            }
                        }                
                    }

                    //qui dovrò fare un Rplus_i se xi o yi stanno vicino ai bordi
                    //ATTENZIONE NON E' COSI' SEMPLICE: è un'interazione topologica. se prima gli faccio controllare i vicini in Xpos e poi quelli in Xplus quelli in xplus non li conta, raggiunge i_neighbour in Xpos...
                    if ( (xi > L-R && xi < L) || (xi > 0 && xi < R) ) {
                        
                        std::vector<double> RXplus_i(n);   //vector of relative distances from i

                        for (int j=0; j<n; j++){
                            if (Xplus[j] != num ) {

                                double xj = Xplus[j];
                                double yj = Ypos[j];    //y è la stessa

                                double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                                double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                                RXplus_i[j] = mod_r;
                                if ( j==i ) { RXplus_i[j] = 0; }

                            }

                        }

                        int i_neighbours = 0;  //conta i vicini di i                   

                        for (double r=0; r<R; r += r0) {   //NB: dr != r0 needed
                            if ( i_neighbours < n_c ){                     
                                for (int j=0; j<n; j++) { 
                                    if ( j != i ) {

                                        if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                            i_neighbours += 1;
                                            Adj[i][j] = 1/n_c;                            
                                        }                            
                                    }
                                }
                            }
                        }                
                        

                    }
                                             
                    

                }
            }
            /*
            if ( n <= 20 ) { 
                for (int i=0; i<n; i++){        //printing Adjacency matrix
                    for (int j=0; j<n; j++){
                        std::cout<< Adj[i][j] << '\t';
                    }
                    std::cout << '\n';
                } 
            }
            */          

            for (int i=0; i<n; ++i) {

                for (int j=0; j<n; ++j){

                    Int[i] += (1/N)* Adj[i][j] * ( tanh(x[j]-x[i])  )  ;   //saving interaction terms

                    /*        
                    if (t>2 && t<3){
                        if (Adj[i][j] != 0) {         
                            std::cout <<"Adj["<<i<<"]["<<j<<"]"<< " was: " << Adj[i][j] <<'\n';
                            std::cout <<"the term added to Int["<<i<<"]"<< " was " << (1/N) * Adj[i][j] * ( tanh(x[j]-x[i])  ) <<'\n';
                        }
                    } */ 
                } 
                //if (t>2 && t<3) { std::cout << "interaction term for " <<i<< " at time " <<t<< " was " << Int[i] << '\n'; }   
            }
            
        } 

        for (int i=0; i<n; i++) {
            if( Int[i] < -1e-8 ) {  
                counter += 1;
            }    
        }
        
        if ( counter > 0 ) { std::cout<< "at time " <<t<< "there were " <<counter<< " negative interaction terms on a total of " <<n<< '\n'; }
        if ( counter >= n-1 ) { 
            std::cout << "******ERROR****** : every interaction term was negative; so every firefly stayed in her state and synchronization was impossible to achieve" <<'\n'; 
        }
        counter = 0;  
        
        t += dt;    //adjourn current time   

        for (int i=0; i<n; ++i){             
            if ( Int[i] > -1e-8 ) {     //termini di interazione quando si sincronizza dell'ordine di -1e-5 -> -1e-7 trascurabile 

                x[i] += 0.2;

            }

            Int[i] = 0;
        }

        for (int i=0; i<n; ++i ){
            if ( x[i] > 1.9999 ) { x[i] = 0; }
        }

        move(Xpos, L, dx);
        move(Ypos, L, dx);

    }
    s0.close();
    fout.close();  
    sync.close();

}