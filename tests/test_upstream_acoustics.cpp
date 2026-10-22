#include "FlowModel.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

void check(bool value,const char* message) {
    if (!value) throw std::runtime_error(message);
}
void near(double a,double b) {
    check(std::isfinite(a+b) && std::abs(a-b)<=1e-14*std::max(std::abs(a),std::abs(b))+1e-25,
          "acoustic coefficient mismatch");
}
int main(int argc,char** argv) {
    try {
        if (argc!=1 && argc!=3) return 2;
        for (int which=1; which<=2; ++which) {
            SimulationParams p;std::string error;
            const bool legacy=which==2;
            if (argc==3) check(p.loadFromFile(argv[which],error)&&p.validate(error),error.c_str());
            else {
                p.rho=1.13;p.mu=1.83e-5;p.ps=1600.;p.dt=1e-5;
                p.L_vt=.175;p.N_vt=10;
                p.glottalMethods.configure("legacy_fluid_all",{});
                p.setUpstreamAcousticProfile(legacy?"legacy_main":"modern");
            }
            check(p.N_sub==(legacy?3:10),"wrong upstream partition count");
            check(p.N_vt==10 && p.L_vt==.175,"downstream circuit changed");
            Geometry g;g.nxsup=3;g.zmax=17.;g.points={{0.,0.,0.},{10.505,0.,17.}};
            FlowModel flow;flow.initialize(p,g,300);
            check(flow.upstreamFlow().size()==std::size_t(p.N_sub+1),"Uu size");
            check(flow.upstreamPressure().size()==std::size_t(p.N_sub+2),"Pu size");
            const double ai=M_PI*p.r_inlet*p.r_inlet,as=M_PI*p.r_sub*p.r_sub;
            const double av=M_PI*p.r_vt*p.r_vt,dx=p.L_sub/p.N_sub,dv=p.L_vt/p.N_vt;
            const auto c=flow.acousticCoefficients();
            near(c[0],p.rho*p.L_inlet/(2*ai));near(c[1],p.L_inlet*ai/(p.rho*p.c_sound*p.c_sound));
            near(c[2],p.rho*dx/((legacy?1.:2.)*as));near(c[3],dx*as/(p.rho*p.c_sound*p.c_sound));
            near(c[4],(-2.5e-5*p.ps+.185)/(as*as)*std::sqrt(p.rho*p.mu*p.c_sound));
            near(c[5],p.rho*dv/av);near(c[6],dv*av/(p.rho*p.c_sound*p.c_sound));
            near(c[7],p.rho*1.1*std::sqrt(av/M_PI)/av);
            near(c[8],(1.6e-3*p.ps+.6)*p.rho*p.c_sound/(9*M_PI*M_PI*av));
            ChannelSections sections;sections.x={0.,5.,10.};sections.area={30.,20.,30.};
            sections.gapCubed={100.,100.,100.};sections.valid={true,true,true};
            double old=0.;
            for (int n=0; n<300; ++n) {
                flow.advance(n*p.dt,p.dt,20.*1e-6,old,sections);old=flow.flowRate();
                check(std::isfinite(old+flow.glottalPressure()+flow.outletPressure()),"nonfinite circuit state");
                for (double x:flow.upstreamFlow())check(std::isfinite(x),"Uu finite");
                for (double x:flow.upstreamPressure())check(std::isfinite(x),"Pu finite");
                if (legacy)check(flow.newtonTrace().converged,"legacy NR convergence");
            }
        }
        std::cout << "PASS 3/10 upstream partitions, SI coefficients, 2Lu loop, same downstream, finite NR evolution\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n';return 1; }
}
