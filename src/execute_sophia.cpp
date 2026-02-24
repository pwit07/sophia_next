#include <iostream>
#include <fstream>
#include <string>

#include "sophia_interface.h"

const double pi = 3.1415926;

// speed of light value [cm/s]
const double c = 2.99792458e10; 

// m_p c^2 value [GeV]
const double mpc2 = 0.938272088;

// (m_p c^2)^2 value [GeV^2]
const double mp2c4 = 0.880354511;

// class computing histogram and N_bin / N_all values
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

        double NbinOverNall(double x);
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

void sophia_histogram::add(const sophiaevent_output &seo){
// ****** OUTPUT *****************************************
// add sophiaevent_output seo data to this class
// *******************************************************
    for (int i = 0; i < seo.Nout; ++i) {
        eventID.push_back(eventCounter);
        PDG.push_back(ID_sophia_to_PDG(seo.outPartID[i]));
        EGeV.push_back(seo.outPartP[3][i]);    
    }
    eventCounter++;
}

void sophia_histogram::print_data(){
// ****** OUTPUT *****************************************
// printed data from the all sophiaevent_output seo passed to this class
// *******************************************************
    for(int i=0;i<eventID.size();i++)
    {
        std::cout << eventID[i] << " " << PDG[i]<< " " << EGeV[i] << "\n";
    }
}

void sophia_histogram::make_hist(){
// ****** OUTPUT *****************************************
// table H[i] filled with a number of particles corresponding to the energy E_i(i) 
// *******************************************************
    bool check = false;
     
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

int sophia_histogram::find_index(double EGeV){
// ****** INPUT ******************************************
// EGeV = energy of produced particle (in lab frame) [GeV]
// ****** OUTPUT *****************************************
// index = index corresponding to the energy EGeV of produced particle []
// *******************************************************
    int index = static_cast<int>((std::log10(EGeV)-logEmin)/dlogE);
        
    if(index == nbin){index--;}
    if(index < 0){index = 0;}

    return index;
}

double sophia_histogram::E_i(int index){
// ****** INPUT ******************************************
// index = index corresponding to the energy EGeV of produced particle []
// ****** OUTPUT *****************************************
// EGeV = energy of produced particle (in lab frame) [GeV]
// *******************************************************
    return std::pow(10.,logEmin+double(index)*dlogE);
}

void sophia_histogram::print_hist(){
// ****** OUTPUT *****************************************
// printed histogram in table
// *******************************************************
    for(int i=0;i<nbin;i++)
    {
        std::cout<<E_i(i)/Ep<<" "<<H[i]<<"\n";
    }
}

double sophia_histogram::NbinOverNall(double x){
// ****** INPUT ******************************************
// x = EGeV / Ep = energy EGeV of produced particle divided by energy of proton (both in lab frame) [GeV/GeV]
// ****** OUTPUT *****************************************
// N_bin / N_all value []
// *******************************************************
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

double compute_NbinOverNall(int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1.){
// ****** INPUT ******************************************
// partID = sophia ID of produced particles: 
//  1000010010 -> proton
// -1000010010 -> anti-proton
//  1000000010 -> neutron
// -1000000010 -> anti-neutron
//          22 -> photon
//          11 -> electron
//         -11 -> positron
//          12 -> nu_e
//         -12 -> anti-nu_e
//          14 -> nu_mu
//         -14 -> anti-nu_mu
//         111 -> pi0
//         211 -> pi+
//        -211 -> pi-
// *******************************************************
// x = EGeV / Ep = energy EGeV of produced particle divided by energy of proton (both in lab frame) [GeV/GeV]
// eta = 4.*eps*Ep / (m_p c^2)^2 []
// theta =angle between incident proton and photon (in lab frame) [degrees]
// N = number of simulated photomeson (proton-photon) interactions []
// nbin = number of energy bins for produced particles (of selected type by partID) []
// Ep = energy of incident proton (in lab frame) [GeV]
// ****** OUTPUT *****************************************
// N_bin / N_all value
// *******************************************************
    static sophia_interface SI;

    sophia_histogram HIST(partID,eta,theta,Ep,nbin,N);

    for(int i=0;i<N;i++){
        sophiaevent_output seo = SI.sophiaevent_mod(Ep, (mp2c4*eta/(4.*Ep)), theta, false);
        HIST.add(seo);
    }

    HIST.make_hist();
    // HIST.print_hist();


    return HIST.NbinOverNall(x);
}
double integrand(int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1.){
// ****** OUTPUT *****************************************
// integrand for computation Phi function (integral over the angles) [cm^3/s]
// *******************************************************

    double eps_prim = eps_prime(eta,theta);
    return c*(1. - std::cos(theta* pi / 180.))*crossection(eps_prim)*compute_NbinOverNall(partID,x,eta,theta,N,nbin,Ep);

}

double compute_Phi(int partID, double x, double eta, int N = 10000, int nbin = 10, double Ep = 1., int Nintegrate = 18){
// ****** OUTPUT *****************************************
// value of Phi function from K&A 2008 [cm^3/s]
// *******************************************************
    std::vector <double> X, Y;
    
    X.push_back(0.);
    Y.push_back(integrand(partID,x,eta,0.,N,nbin,Ep));

    double Xi = 0., dXi = 180./double(Nintegrate);

    for(int i=0;i<Nintegrate-1;i++)
    {
        Xi+=dXi;
        // std::cout<<xi<<"\n";
        X.push_back(Xi);
        Y.push_back(integrand(partID,x,eta,Xi,N,nbin,Ep));
    }

    X.push_back(180.);
    Y.push_back(integrand(partID,x,eta,180.,N,nbin,Ep));


    double F = 0.;
    
    for(int i=0;i<X.size()-1;i++)
    {
        // std::cout<<X[i]<<"  "<<Y[i]<<"\n";
        F+=(Y[i+1]+Y[i])*dXi*0.5;
    }
    
    return F*pi/180.;
}

void Phi2File(std::string path, int partID, double x_a, double x_b, double eta, int N)
{
    std::ofstream outfile;
    outfile.open(path.c_str());


    outfile << "#x\t" << "Phi [cm^3/s]\n";
    
    double logxmin = std::log10(x_a), logxmax = std::log10(x_b), dlogx = (logxmax - logxmin )/double(N);
    double logx = logxmin, x = std::pow(10.,logx);

    while (logx <= logxmax)
    {
        outfile << x << "\t" << compute_Phi(partID,x,eta,4000,50,1e4,18*3) << "\n";
        logx+=dlogx;
        x = std::pow(10.,logx);
    }
    

    outfile.close();
    return;
}

// If you want to see warnings please uncomment line 283 in sophia_interface.cpp file

int main() {
    std::cout.precision(10);

    // generate_events(22,10000,1.,1e1,90.,1e3);

    // double f = compute_Nbin_Nall(22,1.e-3,1.e1,90,10000,10,1e3);
    
    // double f = integrand(22,1.e-3,30.,90,10000,10,1e3);
    // std::cout<<"integrand = "<<f<<"\n";

    // int k = 0;
    // double f = 0.;
    // for(int i=0;i<180;i=i+10)
    // {
    //     f += integrand(22,1.e-3,30.,double(i),10000,10,1e4);
    //     // std::cout<<i<<" integrand = "<<f<<"\n";
    //     k++;
    // }

    // f = f/double(k);
    // std::cout<<"Phi = "<<f<<"\n";

    double eta_0 = 0.313;

    // double x = 1.e-1;
    // double eta = 1.5*eta_0;

    // std::cout<<"Phi = "<<x*compute_Phi(22,x,eta,10000,10,1e4,18)<<"\n";

    // Phi2File("./src/kelner_aharonian_2008/fig2_values/gamma_1.5eta0_sophia.txt",22,1.e-4,1.,1.5*eta_0,20);
    // Phi2File("./src/kelner_aharonian_2008/fig2_values/gamma_30eta0_sophia.txt",22,1.e-4,1.,30.*eta_0,20);

    // Phi2File("./src/kelner_aharonian_2008/fig3_values/positron_1.5eta0_sophia.txt",-11,1.e-4,1.,1.5*eta_0,20);
    Phi2File("./src/kelner_aharonian_2008/fig3_values/positron_30eta0_sophia.txt",-11,1.e-4,1.,30.*eta_0,20); ///!!!

    // Phi2File("./src/kelner_aharonian_2008/fig4_values/muon_antineutrino_1.5eta0_sophia.txt",-14,1.e-4,1.,1.5*eta_0,20);
    // Phi2File("./src/kelner_aharonian_2008/fig4_values/muon_antineutrino_30eta0_sophia.txt",-14,1.e-4,1.,30.*eta_0,20);

    // Phi2File("./src/kelner_aharonian_2008/fig5_values/muon_neutrino_1.5eta0_sophia.txt",14,1.e-4,1.,1.5*eta_0,20);
    // Phi2File("./src/kelner_aharonian_2008/fig5_values/muon_neutrino_30eta0_sophia.txt",14,1.e-4,1.,30.*eta_0,20); // !!!

    return 0;
}