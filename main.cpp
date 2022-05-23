#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <cmath>
#include <fstream>
#include <numeric>
#include "sync.h"

using namespace boost::numeric::odeint;
double N = 20;  //Kuramoto parameters   
double K = 50;                         

double k = 0.1; //Cucker-Smale parameters
double sigma = 1; 
double beta = 1/7;
double R = 10;

double n_c = 5; //Parisi parameter  (# of topological neighbours)

double t = 0.0;    //time related parameters
size_t nSteps = 200;
double dt = 0.1;
double T = 5*dt;
double maxdiff = 0.001;
int Tr = 0;

double L = 30;  //box dimension

//**********************//
int Model_type = 1;     ////////  choose model: 0 = Cucker-Smale, metric interaction. 1 = Parisi, topological  ////////  
//**********************//

typedef std::vector<double> state_type;

struct MCU{             //functor to be passed to do_step for the integration
    void operator() (const state_type &x, state_type &dxdt, const double t) {
        
        for(int i=0; i<N; ++i){    
            dxdt[i] = 1/T; 
        }

    }
};

int main(){

    state_type x(N);   // Initial condition, vector of N elements (N ODEs)
    state_type x1(N);      

    state_type x_t(N);  //needed to save the solution before it's overwritten
    state_type Int(N);  //interaction terms 

    int n = int(N);

    state_type Phases = Phases_generator(N);
    state_type Xpos = Positions_generator(L, n);
    state_type Ypos = Positions_generator(L, n); 

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

    if (Model_type == 0){   //Cucker-Smale
        for (int i=0; i<n; i++){    //filling Adjacency matrix 

            double xi = Xpos[i];
            double yi = Ypos[i];

            for (int j=0; j<n; j++){

                double xj = Xpos[j];
                double yj = Ypos[j];
              
                if( j!=i ) {
                    Adj[i][j] = CS_entries(xi, xj, yi, yj, R, k, sigma, beta); 
                } else if ( j==i ) {
                    Adj[i][j] = 0;
                }
        
            }
        }
    } 

    else if (Model_type == 1){  //Parisi
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
            /*std::cout<< "relative distances from " <<i<< " which is in ( " << Xpos[i] << " , " << Ypos[i] << " ) are: " << '\n';
            for (int j=0; j<n; j++){
                std::cout << Ri[j] << '\t';
            }
            std::cout <<'\n';*/

            int i_neighbours = 0;  //conta i vicini di i
            double dr = 0.001;
            double r0 = dr*2;
            double Rmax = L * 1.5;    //"limit" for r; topological interaction has no metric limit, it's only needed in the for loop. L*sqrt(2) would be enough (the box is L*L), 1.5 chosen for certainty

            for (double r=0; r<Rmax; r += r0) {   //NB: dr != r0 needed
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
    
    if ( n <= 20 ) { 
        for (int i=0; i<n; i++){        //printing Adjacency matrix
            for (int j=0; j<n; j++){
                std::cout<< Adj[i][j] << '\t';
            }
            std::cout << '\n';
        } 
    }

    for (int i=0; i<N; i++){   //setting initial values for x[i]
        x[i] = Phases[i];
    }

    std::fstream fout;
    fout.open("solutions.txt", std::ios::out); 
    std::fstream sync;
    sync.open("Synchronization.txt", std::ios::out); 
    std::fstream check;
    check.open("Initial_conditions.txt", std::ios::out);

    //create stepper:
    runge_kutta4<state_type> rk4; 

    int counter=0;
    std::vector<double> x_new = Phases_generator(n);

    for ( int ii=0; ii<nSteps; ++ii ){  //Integration loop

        if (t==0) {             //print initial conditions
            fout << t << '\t';
            for (int i=0; i<N; i++) {
                fout << x[i] << '\t';
            }
            fout << '\n';   
        }

        if ( t >= 0 ) {             //printing m-l/N : m = # of fireflies in 0, l = # of fireflies in 1 
            double m = 0;
            double l = 0;

            for (int i=0; i<N; i++){
                if ( normalizer(x[i] ) == 0 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore m+=1 " <<'\n';
                    m += 1;
                }
                if ( normalizer(x[i] ) == 1 ) {
                    //std::cout << "x["<<i<<"] at time t = " <<t<< " was " <<x[i] <<" .Its normalized value was " <<normalizer(x[i]) << " and therefore l+=1 " <<'\n';
                    l += 1;
                }
            }

            if (t==0) {
                check << m << '\t' << l << '\t' << '\n'; 
            }
            
            sync << t <<'\t'<< (m-l)/N << '\n';

        }
        
        t += dt;    //adjourn current time 
        
        //rk4.do_step(MCU(), x , t, dt);     //perform one integration step. Solution x is overwritten  
        for(int i=0; i<n; ++i){     //without integration
            x[i] += 0.2;
        }   

        for (int i=0; i<n; ++i){            
            if (Int[i] < 0) {
                x[i] -= 0.2;
            }

            Int[i] = 0;
        }   

        /*for (int i=0; i<N; ++i){
            //std::cout << "Int["<<i<<"] = " <<Int[i] << '\n';

            if( Int[i] < 0 ) {         //if the state is incoherent i-th element will stay there longer (aka the state will still be that of x_t)
                x[i] = x_t[i];   //******************+QUESTO E' SBAGLIATO*****************
                //std::cout << i << " was reset in its t-1 state at time " << t << " as Int[i] was " << Int[i] << '\n';
            }

            Int[i] = 0;   //reset Interaction for the new step
            

        }*/


        fout << t << '\t';      //print solution at time t    
        for (int i=0; i<N; i++) {
            fout << x[i] << '\t';
        }
        fout << '\n';   

        //*****************interazione a t = T, t = T-dt***************************************

        if ( ( t!=0 ) && ( t == 0.4 || dmod(t , T , 10) == 0 || dmod(t-(4*dt) , T , 10) == 0 ) ) {     
            std::cout << "interaction at t = " <<t << '\n';

            for (int i=0; i<N; ++i) {
                
                //x_t[i] = x[i]; //saving states in x_t[i]   //ERROREE NON FARLO

                //std::cout << "evaluating interaction term for " <<i<< '\n';

                for (int j=0; j<N; ++j){

                    //if (i != j) {                             //probably not needed as Chi(i,i) = +1 but Adj[i][i] = 0 

                        //if( Adj[i][j] != 0 ) {

                            Int[i] += /*(1/N)**/ ( Adj[i][j] * Chi(x[i] , x[j], maxdiff) ) ;   //saving interaction terms
                            /*if (t>130 && t<140){
                                
                                    std::cout <<"Adj["<<i<<"]["<<j<<"]"<< " was: " << Adj[i][j] <<'\n';
                                    std::cout <<"Chi(i,j) was " << Chi(x[i] , x[j], maxdiff) << " as i was " <<x[i]<< " and j was " <<x[j] <<'\n';
                                    std::cout <<"the term added to Int["<<i<<"]"<< " was " << Adj[i][j] * Chi(x[i] , x[j], maxdiff ) <<'\n';
                            } */
                        //}
                    //}   
                }
                //std::cout << "interaction term for " <<i<< " at time " <<t<< " was " << Int[i] << '\n';
            }
        } 

        for (int i=0; i<n; i++) {
            if( Int[i] < 0 ) {  
                counter += 1;
            }      
        }
        /*std::cout<< "at time " <<t<< "there were " <<counter<< " negative interaction terms on a total of " <<n<< '\n';
        if ( counter >= n-1 ) { 
            std::cout << "******ERROR****** : every interaction term was negative; so every firefly stayed in her state and synchronization was impossible to achieve" <<'\n'; 
            std::cout << "trying reset" <<'\n';
            state_type x_new = Phases_generator(n);
            for (int i = 0; i<n; ++i){
                x[i] = x_new[i];          ///////////**************************ECCOLO PORCOD** ERA STA MERDATAAAAAAAAAAA RESETTARE FACENDO x[i] = x_t[i] sminchia tutto *********************************
            }
        }
        counter = 0;    */     



        //*********************interazione a t random*************************  si può fare in maniera molto più semplice e più pulita con dmod
        /*
        //std::cout<< "t = " <<t<<'\n';
        double s = trunc(t + 0.001);
        double t1 = t - s;                 //get t1 beetween 0 and 1 (t=198.7 -> t1 = 198.7 - 198 = 0.7)
        //std::cout<<" t1 = " <<t1<<'\n';
        if ( t1 < 0 ) { t1 = -t1; }
        int tf = (t1+0.001)*10;
        //std::cout << "tf = " <<tf<<'\n';
        

        if ( t1 < 0.001 ) {         //which is every second (every 10dt), for example 196.0 -> 0.0, 197.0 -> 0.0 ...

            std::random_device rd;  
            std::mt19937 seed(rd()); 
            std::uniform_int_distribution<int> t_dist(1, 9);
            Tr += t_dist(seed);

        }
        //std::cout<< " Tr = " <<Tr<<'\n';
        //std::cout<< " Tf-T = " <<tf-Tr<<'\n';
        
        if ( (Tr != 0) && (tf>=Tr) && ( tf - Tr < 0.001 ) ) {
            std::cout << "interaction at t = " <<t << '\n';

            for (int i=0; i<n; ++i) {

                //jumpStart[i] = 0;
                
                x_t[i] = x[i]; //saving states in x_t[i]

                for (int j=0; j<N; ++j){

                    if (i != j) {                             //needed for Chi(i,i) = +1, not needed for Adj[i][i] = 0 

                        if( Adj[i][j] != 0 ) {

                            Int[i] += (1/N) *   ( Adj[i][j] * Chi(x[i] , x[j], maxdiff) ) ;   //saving interaction terms 

                        }
                    }   
                }

            }
            Tr=0;

            for (int i=0; i<n; i++) {
                if( Int[i] < 0 ) {  
                    counter += 1;
                }      
            }
            //std::cout<< "at time " <<t<< "there were " <<counter<< " negative interaction terms on a total of " <<n<< '\n';
            if ( counter >= n-1 ) { 
                //std::cout << "******ERROR****** : every interaction term was negative; so every firefly would have stayed in her state and synchronization would have been impossible to achieve" <<'\n'; 
                
                for (int k=0; k<n; ++k) {   //resettando alcuni termini di interazione a +1 non si sistema nulla; all'interazione successiva vengono di nuovo valutati come negativi 
                    if (k % 10 == 0) {
                        Int[k] = +1;
                        //std::cout<< "Int[" <<k<< "] was set to +1 to fix error" <<'\n'; 
                        jumpStart[k] = 1;    //the k-th term of Int[i] will not be reset to 0 until next interaction
                    }
                }
                std::cout << "some interaction terms were set to +1 in order to make synchronization achievable again" <<'\n';
                

            }
            counter = 0;                   

        } */

    }
    check.close();
    fout.close();  
    sync.close();

}


//SEMPRE LO STESSO PROBLEMA: IL TERMINE DI INTERAZIONE RISULTA NEGATIVO PER TUTTE: QUESTO LE PORTA A FERMARSI CONTEMPORANEAMENTE NELLO STATO IN CUI SONO. COSI' NON SI SINCRONIZZANO MAI...
//COSA LO CAUSA??