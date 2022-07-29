#include "sync.h"

#include <vector>
#include <random>
#include <cmath>

int dmod(double t_, double T_, double res) {

    if ( T_-0.01 >t_ ) { return -1; }    //all the 0.01 are needed to account for double imprecision

    else {
        double t1 = t_+0.01;
        double T1 = T_+0.01;
        int t = trunc(t1*res);
        int T = trunc(T1*res);
        return t % T;
    }

}

int normalizer(double Theta) {
    int theta = trunc(Theta + 0.001);       //truncs theta in the proper range (2.7 --> 2)
    if (theta % 2 == 0) { theta = 0; }      //odd range -> off, even range -> on
    else if (theta % 2 != 0) { theta = 1; }
    return theta;
}

std::vector<double> Phases_generator(int N){
    std::vector<double> Phases(N);

    std::random_device rd;  
    std::mt19937 seed(rd()); 
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for (int i=0; i<N; i++){
        double rand = distribution(seed);
        if (rand <= 0.5){ Phases[i] = 0; } 
        else if (rand > 0.5) { Phases[i] = 1; }
    }
    return Phases;
}

std::vector<double> Positions_generator(double L, int N){
    std::vector<double> Positions(N);

    std::random_device rd;  
    std::mt19937 seed(rd()); 

    std::uniform_real_distribution<double> pos_dist(0.0, L);

    for (int i=0; i<N; i++){
        double Xrand = pos_dist(seed);
        Positions[i] = Xrand;
    }

    return Positions;    
}

double CS_entries(double xi, double xj, double yi, double yj, double R, double k, double sigma, double beta){

    double r_ij = (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi);  //modulo quadro di ri-rj
    double mod_r = std::sqrt(r_ij);  //modulo di ri-rj

    if ( mod_r <= R ) {
        return ( k / ( std::pow( ((sigma*sigma) + r_ij) , beta ) ) );
    } else if ( mod_r > R ) {
        return 0;
    }

}

void move (std::vector<double>& pos, double L, double dx){

    int n = pos.size();

    std::random_device rd;  
    std::mt19937 seed(rd()); 
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for (int i=0; i<n; i++){
        double rand = distribution(seed);
        double way = distribution(seed);
        
        if (rand <= 0.5){       //50% di probabilità che si muova o meno
            if (way < 0.5) {    //50% di probabilità che vada nella direzione positiva o negativa
                pos[i] += dx;
            }
            else if (way > 0.5){
                pos[i] -= dx;
            }
        } 

        if ( pos[i] > L ) { pos[i] -= L; }      //pareti periodiche
        if ( pos[i] < 0 ) { pos[i] += L; }
        
    }

}


















////////////matrici di adiacenza indipendenti dal tempo (senza move)
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




//////come avevo fatto il controllo delle strisce per il metrico:
/*


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

                    } */       





/*struct MCU{             //functor to be passed to do_step for the integration
    void operator() (const state_type &x, state_type &dxdt, const double t) {
        
        for(int i=0; i<N; ++i){    
            dxdt[i] = 1/T; 
        }

    }
};*/


    //create stepper:
    //runge_kutta4<state_type> rk4; 

/*
                if (t>95 && t<97){
                    std::cout << "at time t: " <<t<< " negative interaction term for i = " <<i<< '\n';
                    std::cout << " infact its state was" << x[i] << " while its neighbours were ";
                    for (int j=0; j<n; ++j){
                        if (Adj[i][j] != 0) {
                            std::cout << x[j] <<'\t';
                        }
                    }
                    std::cout <<'\n';
                }

*/





//////////////definizione di Chi che fa sempre o +1 o -1 a seconda che siano sincronizzate o meno
/*bool sameRange(double Theta_i, double Theta_j) {  

    int theta_i = normalizer(Theta_i);
    int theta_j = normalizer(Theta_j); 

    int Phase_diff = theta_i - theta_j;

    if (Phase_diff == 0) { return true; } //coherent phases (0,0) or (1,1)
    if (Phase_diff != 0) { return false; } //incoherent phases (1,0) or (0,1)
}

double Chi(double theta_i, double theta_j, double maxdiff) {  //to be synchronized the fireflies need to be in the same range AND differ less than maxdiff

    bool samerange = sameRange(theta_i, theta_j);

    double xi = trunc(theta_i + 0.001);
    double xj = trunc(theta_j + 0.001);
    theta_i -= xi;
    theta_j -= xj;

    double phase_diff = theta_i - theta_j;
    if ( phase_diff < 0 ) { phase_diff = -phase_diff; }

    if ( phase_diff < maxdiff + 0.001 && samerange ) {
        return 0;
    } 
    if ( phase_diff > maxdiff + 0.001 || !samerange ) {
        return -1;
    } 
} */



    //(*)
        //rk4.do_step(MCU(), x , t, dt);     //perform one integration step. Solution x is overwritten  
        /*for (int i=0; i<N; ++i){
            //std::cout << "Int["<<i<<"] = " <<Int[i] << '\n';

            if( Int[i] < 0 ) {         //if the state is incoherent i-th element will stay there longer (aka the state will still be that of x_t)
                x[i] = x_t[i];   //******************+QUESTO E' SBAGLIATO*****************
                //std::cout << i << " was reset in its t-1 state at time " << t << " as Int[i] was " << Int[i] << '\n';
            }

            Int[i] = 0;   //reset Interaction for the new step
            

        }*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//prova che avevo fatto per avere tutte le lucciole con 0 ed 1 ma non sapevo come farle continuare col loro periodo dopo che si erano stoppate e non farle flashare tutte a t = kT
    
/*if ( ( t!=0 ) && ( dmod(t , T , 10) == 0 ) )  {         //per t=T le lucciole cambiano stato...
            for (int i=0; i<N; ++i) {
                if ( Int[i] > -0.00000001 ) {
                    //std::cout << "changing x[" <<i<< "]" << " at time t=T= " <<t<< '\n';
                    //std::cout << "previous state was " << x[i] << '\n';
                    //changeState(x[i]);
                    x[i] += 0.2;
                    //std::cout << "new state is " << x[i] << '\n';
                }
                else if ( Int[i] < -0.00000001 ) {    // ...a meno che il termine di interazione sia negativo: in quel caso rimangono per un altro step
                    to_change[i] = 1;               //mi salvo l'indice della lucciola da cambiare allo step successivo
                    //std::cout << "interaction term for " <<i<< " was negative at t= " <<t<< " so it will stay in its state for another step" <<'\n';
                    t1 = t+dt;
                }
            }
        }

        if ( t == t1 ) {   //quindi vado a cambiargli stato a T+dt
            t1 = -15;
            for (int i=0; i<n; i++) {
                if ( (to_change[i] != 0) && (Int[i] > -0.000000001) ) {
                    //std::cout << "changing x[" <<i<< "] state at t= " <<t<< '\n';
                    //std::cout << "previous state was " << x[i] << '\n';
                    //changeState(x[i]);
                    x[i] += 0.2;
                    //std::cout << "new state is " << x[i] << '\n';
                    to_change[i] = 0;
                }
                else if ( (to_change[i] != 0) && (Int[i] < -0.000000001) ) {
                    //std::cout << "x[" <<i<< "] should have been changed but its interaction term was still negative at T+dt; it will wait another step" <<'\n';
                    t2 = t+dt;
                }
            }
        }

        if ( t == t2 ){
            for (int i = 0; i<n; ++i) {
                if ( to_change[i] != 0 ) {
                    //std::cout << "changing x[" <<i<< "] at t= " <<t<< "after 2 steps in its state" <<'\n';
                    //changeState(x[i]);
                    x[i] += 0.2;
                    to_change[i] = 0;
                }
                t2= -15;
            }
        } 
        for (int i=0; i<n; ++i) {
            Int[i] = 0;                 //resetto il termine di interazione
        } */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
                }
            }
            counter = 0;                   

        } */

