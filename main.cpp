#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include <algorithm> 
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 500;  //Kuramoto parameter                           

double k = 1;   //Cucker-Smale parameters
double sigma = 1; 
double beta = 7.0 ; 
double R = 5.0;

double n_c = 1.0001; //Parisi parameter  (# of topological neighbours)

size_t nSteps = 500;
double dt = 0.1;   //time step
double T = 5*dt;

double L = 24.5;  //box dimension
double dx = 0.1;  //spatial step

double nSim = 500;  //# di simulazioni da eseguire

//**********************************************************************************************************************************************
//**********************//
int model_type = 1;     ////////  choose model: (0) = global interaction. (1) = Cucker-Smale, metric interaction. (2) = Parisi, topological.   ////////  
//**********************//
//**********************************************************************************************************************************************


int main(){

    double nc_max = 8.0;
    double beta_max = 13.0;
    std::cout << "starting loop" <<'\n';

    while (beta < beta_max) {
        std::cout << "beta = " << beta <<'\n';

    std::vector<int> t_average;

    std::fstream fout;
    std::fstream sync;
    std::fstream positions;  

    double probability = 0;   //probabilità che il sistema si sincronizzi sulle nSim simulazioni
    double N_sync = 0;   //numero di volte in cui il sistema si sincronizza sulle nSim simulazioni
    
    for (int s = 0; s<nSim; s++) { 

        double t = 0.0;

        fout.open("solutions.txt", std::ios::out); 
        sync.open("Synchronization.txt", std::ios::out);
        positions.open("Positions.txt", std::ios::out); 

        std::vector<double> x(N);   // Initial condition, vector of N elements (N ODEs)    

        std::vector<double> Int(N);  //interaction terms 

        int n = int(N);

        std::vector<double> Phases = Phases_generator(N);
        std::vector<double> Xpos = Positions_generator(L, n); //posizioni delle lucciole
        std::vector<double> Ypos = Positions_generator(L, n);

        std::vector<double> Xcheck(n);    //prolungamenti della scatola per avere pareti periodiche
        std::vector<double> Ycheck(n);

        std::vector<double> System(nSteps);  //vettore su cui scrivo m-l/N, ogni elemento è uno step temporale

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

        */

        for (int i=0; i<N; i++){   //setting initial values for x[i]
            x[i] = Phases[i];
            Int[i] = 0;
        }

        int counter=0;

        for ( int ii=0; ii<nSteps; ++ii ){  

            //printing m-l/N : m = # of fireflies in 0, l = # of fireflies in 1 
            double m = 0;
            double l = 0;
            for (int i=0; i<N; i++){
                if ( normalizer(x[i]) == 0 ) {
                    m += 1;
                }
                if ( normalizer(x[i]) == 1 ) {
                    l += 1;
                }
            }
            double state = (m-l)/N;

            sync << t <<'\t'<< state << '\n';
            System[ii] = state;

            fout << t << '\t';      //print solution at time t    
            for (int i=0; i<N; i++) {
                fout << x[i] << '\t';
            }
            fout << '\n';      

            double avg_neighbours = 0;


            //****************************interazione a t = T, t = T-dt***************************************//

            if ( ( t!=0 ) && ( t == 0 +(T-dt) || dmod(t , T , 10) == 0 || dmod(t-(T-dt) , T , 10) == 0 ) ) {     
                //std::cout << "interaction at t = " <<t << '\n';

                for (int i=0; i<n; i++){

                    for (int k=0; k<n; k++) {
                        Xcheck[k] = Xpos[k];
                        Ycheck[k] = Ypos[k];
                    } 

                    double xi = Xpos[i];
                    double yi = Ypos[i];

                    if ( xi > L-R && xi < L ) {
                        for (int j=0; j<n; j++) {
                            if ( Xcheck[j] > 0 && Xcheck[j] < R ) {
                                Xcheck[j] += L;
                            }
                        }
                    }
                    else if ( xi > 0 && xi < R ) {
                        for (int j=0; j<n; j++) {
                            if ( Xcheck[j] > L-R && Xcheck[j] < L ) {
                                Xcheck[j] -= L;
                            }
                        }
                    }
                    if ( yi > L-R && yi < L ){
                        for (int j=0; j<n; j++) {
                            if ( Ycheck[j] > 0 && Ycheck[j] < R ) {
                                Ycheck[j] += L;
                            }
                        }                        
                    }
                    else if ( yi > 0 && yi < R ) {
                        for (int j=0; j<n; ++j){
                            if ( Ycheck[j] > L-R && Ycheck[j] < L ) {
                                Ycheck[j] -= L;
                            }
                        }
                    }
            
                    if (model_type == 1){   //Cucker-Smale
                    
                        for (int j=0; j<n; j++){        //tutte devono controllare Xpos, Ypos

                            double xj = Xcheck[j];
                            double yj = Ycheck[j];
              
                            if( j!=i ) {

                                double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  //modulo quadro di ri-rj
                                double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                                if ( mod_r <= R ) {
                                    //std::cout << "xi in (" << xi <<" , "<< yi<< ") was considered a neighbour of xj in (" <<xj<< " , " <<yj<< ") because their distance was " <<mod_r <<'\n';
                                    Adj[i][j] = ( k / ( std::pow( ((sigma*sigma) + r_ij) , beta ) ) );
                                } else if ( mod_r > R ) {
                                    //std::cout << "xi in (" << xi <<" , "<< yi<< ") was NOT considered a neighbour of xj in (" <<xj<< " , " <<yj<< ") in fact their distance was " <<mod_r <<'\n';
                                    Adj[i][j] = 0;
                                } 
                            
                            }
                        }
                 
                    } 
            
                    else if (model_type == 2){  //Parisi

                        std::vector<double> Ri(n);   //vector of relative distances from i

                        for (int j=0; j<n; j++){
                            double xj = Xcheck[j];
                            double yj = Ycheck[j];

                            double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi); 
                            double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

                            Ri[j] = mod_r;
                            if ( j==i ) { Ri[j] = 0; }
                        }                    
                    
                        int i_neighbours = 0;  //conta i vicini di i
                        double dr = 0.01;
                        double r0 = dr*2;
                        int nc = trunc(n_c + 0.1);

                        for (double r=0; r<R; r += r0) {   //NB: dr != r0 needed
                            if ( i_neighbours < nc ){                     
                                for (int j=0; j<n; j++) { 
                                    if ( j != i ) {

                                        if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                            //std::cout << "xj in (" << Xcheck[j] <<" , "<< Ycheck[j]<< ") was considered a neighbour of xi in (" << xi <<" , "<< yi<< ") because their distance was " << Ri[j] <<'\n';
                                            i_neighbours += 1;
                                            Adj[i][j] = 1/n_c;                            
                                        }                            
                                    }
                                }
                            }                
                        }

                        /*
                        double i_n = 0;
                        double Rmax = R;    


                        for (double r=0; r<Rmax; r += r0) {    //da togliere nel for col movimento
                        for (int j=0; j<n; j++) { 
                            if ( j != i ) {
                                if ( r < Ri[j]+dr && r > Ri[j]-dr ) {
                                    i_n += 1;
                                }
                            
                            }
                        }                
                        }

                        avg_neighbours += (i_n / N);   */                 
                
                    }
                
                }
                    //std::cout << " # of average neighbours was " << avg_neighbours << '\n';

                    /*
                    if ( n <= 20 ) { 
                    std::cout << "printing adjacency matrix at time t = " <<t<< '\n\n';
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
        
                //if ( counter > 0 ) { std::cout<< "at time " <<t<< "there were " <<counter<< " negative interaction terms on a total of " <<n<< '\n'; }
                if ( counter >= n-1 ) { 
                    std::cout << "******ERROR****** : every interaction term was negative; so every firefly stayed in her state and synchronization was impossible to achieve" <<'\n'; 
                }
                counter = 0;  
        
                t += dt;    //adjourn current time   

                for (int i=0; i<n; ++i){             
                    if ( Int[i] > -1e-8 ) {     //termini di interazione quando si sincronizzano dell'ordine di -1e-5 -> -1e-7 trascurabile 

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

        fout.close();  
        sync.close();

        std::fstream data;
        data.open("data.txt", std::ios::app); 

        std::vector<int> times;

        for (int i=0; i<nSteps; i++) {      //sistema sincronizzato quando ho 5 step con +1 e 5 con -1 --> salvo su isto i tempi (i) in cui accade
            if ( i < nSteps - 6 ) {
                if (System[i] == 1 && System[i+1] == 1 && System[i+2] == 1 && System[i+3] == 1 && System[i+4] == 1 && System[i+5] == -1) {
                    times.push_back(i);      
                } 
                else if ( System[i] == -1 && System[i+1] == -1 && System[i+2] == -1 && System[i+3] == -1 && System[i+4] == -1 && System[i+5] == 1 ) {
                    times.push_back(i);
                }
            }
        }

        if (times.size() != 0) {
            data << times[0] << '\n';        //è il primo valore di times quello importante (prima volta in cui si sincronizza il sistema)
            t_average.push_back(times[0]);
            N_sync += 1;
            probability += 1/nSim;
        } 
        else {
            data << nSteps << '\n';
        }
        data.close();

    }

    int nc = trunc(n_c);

    std::cout << nSim << " simulations performed" <<'\n';
    if (model_type == 1) { std::cout << " For beta = " << beta; }
    if (model_type == 2) { std::cout << " For n_c = " << nc; }

    std::cout << " the probability of the system sychronizing turned out to be P = " << probability <<'\n';

    //media temporale per ricavare <t> in cui il sistema si sincronizza e la deviazione standard (errore su <t>)
    double _dt_ = 0;
    double Sum = 0;
    for(int i = 0; i<t_average.size(); i++) {
        _dt_ += t_average[i] / N_sync ;
    }
    for(int i = 0; i<t_average.size(); i++) {
        Sum += (t_average[i] - _dt_) * (t_average[i] - _dt_) ;
    }

    double sigma_t = std::sqrt( Sum/N_sync );

    if (model_type == 1) {
        std::cout << "beta = " <<beta<< "; <dt> = ( " << _dt_ << " +/- " << sigma_t <<" ) " <<'\n';
    }
    else if (model_type == 2) {
        std::cout << "n_c = " <<nc<< "; <dt> = ( " << _dt_ << " +/- " << sigma_t <<" ) " <<'\n';
    }
    

    beta += 1.0;

    }

    std::cout << "loop ended" <<'\n';

}