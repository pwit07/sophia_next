#include <iostream>
#include <fstream>
#include <string>

#include "sophia_interface.h"

const double pi = 3.1415926;

// value of m_{p} c^{2} [GeV]
const double mpc2 = 0.938272088;

// value of (m_{p} c^{2})^{2} [GeV^2]
const double mp2c4 = 0.880354511;


class sophia_histogram
{
private:
    bool init_check = false;
    int eventCounter = 1, nbin, partID;
    double x, eta, theta, E0, Emin, Emax, logEmin, logEmax, dlogE;
    std::vector <int> eventID;
    std::vector <int> PDG;
    std::vector <double> EGeV;
    std::vector <int> H;
public:
    sophia_histogram(double nx, double neta, double ntheta, int nnbin, int npartID);
    ~sophia_histogram();
    
    double x_i(int index);
    void add(const sophiaevent_output &seo);
    void print_data();
    void make_hist();
    void print_hist();
};

sophia_histogram::sophia_histogram(double nx, double neta, double ntheta, int nnbin, int npartID)
{
    x = nx;
    eta = neta;
    theta = ntheta;

    nbin = nnbin;
    partID = npartID;
    
    for(int i=0;i<nbin;i++)
    {
        H.push_back(0);
    }
}

sophia_histogram::~sophia_histogram()
{
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
    
    std::cout<<logEmin<<" "<<logEmax<<" "<<dlogE<<"\n";
    
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
            index = static_cast<int>((std::log10(EGeV[i])-logEmin)/dlogE);
        
            if(index == nbin){index--;}
            if(index < 0){index = 0;}

            H[index]++;
        }

        // std::cout<<index<<"\n";
    }

}

double sophia_histogram::x_i(int index)
{
    return std::pow(10.,logEmin+double(index)*dlogE);
}

void sophia_histogram::print_hist()
{
    for(int i=0;i<nbin;i++)
    {
        std::cout<<x_i(i)<<" "<<H[i]<<"\n";
    }
}


double eta(double E0, double eps){
// ****** INPUT ******************************************
// E0 = energy of incident proton (in lab frame) [GeV]
// eps = energy of incident photon (in lab frame) [Gev]
// ****** OUTPUT *****************************************
// eta = normalized value of E0*eps []
// *******************************************************
    return 4.*E0*eps/mp2c4;
}

double eps_prime(double eta, double theta){
// ****** INPUT ******************************************
// eta = eta(E0,eps) []
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

void generate_events(double N, double x, double eta, double theta, double E0 = 1.0){

    static sophia_interface SI;

    sophia_histogram HIST(x,eta,theta,10,22);

    for(int i=0;i<N;i++){
        sophiaevent_output seo = SI.sophiaevent_mod(E0, (mp2c4*eta/(4.*E0)), theta, false);
        HIST.add(seo);
    }


    HIST.print_data();

    HIST.make_hist();
    HIST.print_hist();
  
    return;
}

int main() {
    std::cout.precision(10);

    bool onProton = true;
    // bool onProton = false;
    double Ein = 1e9;  // GeV
    double eps = 1e-5;  // GeV
    // bool declareChargedPionsStable = true;
    bool declareChargedPionsStable = false;

    std::ofstream outfile;
    outfile.open("outData.csv");

    outfile << "eventID\t" << "partID\t" << "Px\t" << "Py\t" << "Pz\t" << "EGeV\t" << "m\n";

    int nEvent = 10000;
    int Nout, eventID;
    for (int k = 0; k < nEvent; ++k) {

        sophia_interface SI;

        // std::cout << "Event # " << k + 1 << std::endl;

        sophiaevent_output seo = SI.sophiaevent(onProton, Ein, eps, declareChargedPionsStable);


        // std::cout << "\nevent result # " << (k + 1) << std::endl;
        Nout = seo.Nout;
        // std::cout << "N = " <<Nout << std::endl;
        for (int i = 0; i < Nout; ++i) {
        //     std::cout << "ID = " << ID_sophia_to_PDG(seo.outPartID[i]) << std::endl;
            eventID = k + 1;
                    outfile << eventID << "\t"
                            << ID_sophia_to_PDG(seo.outPartID[i]) << "\t"
                            << seo.outPartP[0][i] << "\t"
                            << seo.outPartP[1][i] << "\t"
                            << seo.outPartP[2][i] << "\t"
                            << seo.outPartP[3][i] << "\t"
                            << seo.outPartP[4][i] << "\n";
            // for (int j = 0; j < 5; ++j) {
            //     // if (j == 3) std::cout << seo.outPartP[j][i] << std::endl;
            //     if (j == 3) {
            //         int eventID = k + 1;
            //         outfile << eventID << "\t"
            //                 << ID_sophia_to_PDG(seo.outPartID[i]) << "\t"
            //                 << seo.outPartP[j][i] << "\n";
            //     }
            // }
        }
        // std::cout << "\n---------------------------------" << std::endl;
    }
    outfile.close();

    // std::cout<<"Done!  "<<crossection(1.)<<std::endl;
    // std::cout<<"Done!  "<<crossection(10.)<<std::endl;
    generate_events(1000,1.,10e1,0.);

    return 0;
}