#include <iostream>
#include <fstream>
#include <string>

#include "sophia_interface.h"

const double pi = 3.1415926;

// speed of light value [cm/s]
const double c = 2.99792458e10; 

// value of m_{p} c^{2} [GeV]
const double mpc2 = 0.938272088;

// value of (m_{p} c^{2})^{2} [GeV^2]
const double mp2c4 = 0.880354511;


class sophia_histogram
{
    private:
        int eventCounter = 1, nbin, partID, N;
        double eta, theta, Ep, Emin, Emax, logEmin, logEmax, dlogE;
        std::vector <int> eventID;
        std::vector <int> PDG;
        std::vector <double> EGeV;
        std::vector <int> H;
    public:
        sophia_histogram(int npartID, double neta, double ntheta, double nEp, int nnbin, int nN);
        ~sophia_histogram(){};
        
        int find_index(double E);
        double E_i(int index);
        
        void add(const sophiaevent_output &seo);
        void print_data();
        void make_hist();
        void print_hist();

        double Nbin_Nall(double x);
};

sophia_histogram::sophia_histogram(int npartID, double neta, double ntheta, double nEp, int nnbin, int nN)
{
    eta = neta;
    Ep = nEp;
    theta = ntheta;

    nbin = nnbin;
    N = nN;
    partID = npartID;
    
    for(int i=0;i<nbin;i++)
    {
        H.push_back(0);
    }
}

void sophia_histogram::add(const sophiaevent_output &seo)
{
    for (int i = 0; i < seo.Nout; ++i) {
        eventID.push_back(eventCounter);
        PDG.push_back(ID_sophia_to_PDG(seo.outPartID[i]));
        EGeV.push_back(seo.outPartP[3][i]);    
    }
    eventCounter++;
}

void sophia_histogram::print_data()
{
    for(int i=0;i<eventID.size();i++)
    {
        std::cout << eventID[i] << " " << PDG[i]<< " " << EGeV[i] << "\n";
    }
}

void sophia_histogram::make_hist()
{
    bool check = false;
    
    // std::cout<<logEmin<<" "<<logEmax<<" "<<dlogE<<"\n";
    
    for(int i=0;i<eventID.size();i++)
    {

        if(PDG[i] == partID && check == false)
        {
            Emin = EGeV[0]; 
            Emax = EGeV[0];
            check = true;
        }
        if(PDG[i] == partID && check == true)
        {
            if(EGeV[i]<Emin){Emin = EGeV[i];}
            if(EGeV[i]>Emax){Emax = EGeV[i];}
        }
    }

    logEmin = std::log10(Emin);
    logEmax = std::log10(Emax);
    dlogE = (logEmax-logEmin)/nbin;

    for(int i=0;i<nbin;i++)
    {
        H[i] = 0;
    }

    int index; 

    for(int i=0;i<eventID.size();i++)
    {
        if(PDG[i] == partID)
        {
            index = find_index(EGeV[i]);

            H[index]++;
        }
    }

}

int sophia_histogram::find_index(double EGeV)
{
    int index = static_cast<int>((std::log10(EGeV)-logEmin)/dlogE);
        
    if(index == nbin){index--;}
    if(index < 0){index = 0;}

    return index;
}

double sophia_histogram::E_i(int index)
{
    return std::pow(10.,logEmin+double(index)*dlogE);
}

void sophia_histogram::print_hist()
{
    for(int i=0;i<nbin;i++)
    {
        std::cout<<E_i(i)/Ep<<" "<<H[i]<<"\n";
    }
}

double sophia_histogram::Nbin_Nall(double x)
{
    int index = find_index(x*Ep);

    return double(H[index])/double(N);
}


double eta(double Ep, double eps){
// ****** INPUT ******************************************
// Ep = energy of incident proton (in lab frame) [GeV]
// eps = energy of incident photon (in lab frame) [Gev]
// ****** OUTPUT *****************************************
// eta = normalized value of Ep*eps []
// *******************************************************
    return 4.*Ep*eps/mp2c4;
}

double eps_prime(double eta, double theta){
// ****** INPUT ******************************************
// eta = eta(Ep,eps) []
// theta = angle between incident proton and photon (in lab frame) [degrees]
// ****** OUTPUT *****************************************
// eps_prime = energy of photon in PRF (proton rest frame) [GeV]
// *******************************************************
    return 0.25*eta*mpc2*(1-std::cos(theta* pi / 180.));
}

double crossection(double eps_prime){
// ****** INPUT ******************************************
// eps_prime = eps_prime(eta, theta) energy of photon in PRF (proton rest frame) [GeV]
// ****** OUTPUT *****************************************
// total corssection of photomeson interaction [cm^2]
// *******************************************************
    static sophia_interface SI; 
                                          // micro barn to cm^2
    return SI.crossection(eps_prime,3,13)*1.e-30; 
}

// void generate_events(int partID, double N, double x, double eta, double theta, double Ep = 1.){

//     static sophia_interface SI;

//     sophia_histogram HIST(partID,eta,theta,Ep,10,N);

//     for(int i=0;i<N;i++){
//         sophiaevent_output seo = SI.sophiaevent_mod(Ep, (mp2c4*eta/(4.*Ep)), theta, false);
//         HIST.add(seo);
//     }


//     // HIST.print_data();

//     HIST.make_hist();
//     HIST.print_hist();
  
//     return;
// }

double compute_Nbin_Nall(int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1.){

    static sophia_interface SI;

    sophia_histogram HIST(partID,eta,theta,Ep,nbin,N);

    for(int i=0;i<N;i++){
        sophiaevent_output seo = SI.sophiaevent_mod(Ep, (mp2c4*eta/(4.*Ep)), theta, false);
        HIST.add(seo);
    }

    HIST.make_hist();
    // HIST.print_hist();


    return HIST.Nbin_Nall(x);
}
double integrand(int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1.){

    double eps_prim = eps_prime(eta,theta);
    return c*(1. - std::cos(theta* pi / 180.))*crossection(eps_prim)*compute_Nbin_Nall(partID,x,eta,theta,N,nbin,Ep);

}

int main() {
    std::cout.precision(10);

    // generate_events(22,10000,1.,1e1,90.,1e3);

    // double f = compute_Nbin_Nall(22,1.e-3,1.e1,90,10000,10,1e3);
    
    // double f = integrand(22,1.e-3,30.,90,10000,10,1e3);
    // std::cout<<"integrand = "<<f<<"\n";

    int k = 0;
    double f = 0.;
    for(int i=10;i<180;i=i+10)
    {
        f += integrand(22,1.e-3,30.,double(i),10000,10,1e3);
        // std::cout<<i<<" integrand = "<<f<<"\n";
        k++;
    }

    f = f/double(k);
    std::cout<<"phi = "<<f<<"\n";

    return 0;
}