#include "KanayaDiagnostics.h"
#include "Simulation.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>

void expect(bool ok) { if(!ok)throw std::runtime_error("Kanaya unit expectation failed"); }
int main(int argc,char** argv) {
    if(argc==3 && std::string(argv[1])=="--smoke") {
        Simulation s;s.initialize(argv[2]);s.run();
        std::ofstream trace(s.runDir/"regression_state.bin",std::ios::binary);
        auto write=[&](const std::string& name,const std::vector<double>& values){
            char label[32]{};name.copy(label,std::min(name.size(),sizeof(label)-1));
            std::uint64_t count=values.size();trace.write(label,32);trace.write(reinterpret_cast<const char*>(&count),8);
            trace.write(reinterpret_cast<const char*>(values.data()),values.size()*8);
        };
        write("Ug_history",s.fCalc.Ug);write("section_area",s.fCalc.harea);write("pressure",s.fCalc.psurf);
        write("qL",s.stateL.q);write("qR",s.stateR.q);write("qdotL",s.stateL.qdot);write("qdotR",s.stateR.qdot);
        write("qddotL",s.stateL.qddot);write("qddotR",s.stateR.qddot);
        for(int side=0;side<2;++side){auto& st=side==0?s.stateL:s.stateR;std::vector<double> xyz;
            for(auto p:st.disp){xyz.push_back(p.ux);xyz.push_back(p.uy);xyz.push_back(p.uz);}write(side==0?"coordinatesL":"coordinatesR",xyz);
            std::vector<double> force;const auto& x=side==0?s.fCalc.prevFdisXL:s.fCalc.prevFdisXR;const auto& y=side==0?s.fCalc.prevFdisYL:s.fCalc.prevFdisYR;
            for(int i=0;i<int(x.size());++i)for(int j=0;j<int(x[i].size());++j){force.push_back(x[i][j]);force.push_back(y[i][j]);}write(side==0?"adopted_contactL":"adopted_contactR",force);}
        std::vector<double> pairs;for(int j=0;j<int(s.fCalc.contactPairCache.size());++j)for(auto p:s.fCalc.contactPairCache[j]){pairs.push_back(j);pairs.push_back(p.first);pairs.push_back(p.second);}write("contact_candidate_cache",pairs);
        if(!trace)return 1;return 0;
    }
#ifndef ENABLE_KANAYA_DIAGNOSTICS
    return 2;
#else
    if(argc!=2)return 2;
    for(double sign:{-1.0,1.0})for(bool left:{true,false}) {
        double outward=sign*(left?-1:1);
        expect(KanayaDiagnostics::openingForce(sign,left,2*outward)==2);
        expect(KanayaDiagnostics::openingForce(sign,left,-3*outward)==-3);
    }
    expect(KanayaDiagnostics::includedInterval(0,3));
    expect(KanayaDiagnostics::includedInterval(1,3));
    expect(!KanayaDiagnostics::includedInterval(2,3));
    expect(!KanayaDiagnostics::includedInterval(-1,3));
    expect(KanayaDiagnostics::scaledResidual(1e-8,1)>0.99);
    expect(KanayaDiagnostics::scaledResidual(2e-8,1)>1);
    expect(KanayaDiagnostics::scaledResidual(2e-12,0)>1);
    Simulation s;s.params.nstep=3;s.params.dt=1e-5;s.params.nmode=1;
    s.runDir=fs::path(argv[1])/std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
    for(int side=0;side<2;++side) {
        auto& g=side==0?s.geomL:s.geomR;auto& st=side==0?s.stateL:s.stateR;auto& md=side==0?s.mdataL:s.mdataR;
        g.nsurfl=g.nsurfz=g.nxsup=3;g.nPoints=9;g.surfp.assign(3,std::vector<int>(3));
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){g.surfp[i][j]=g.points.size();g.points.push_back({double(i),side==0?0.0:2.0,double(j)});}
        st.initialize(9,1,3,g);st.qf=st.qfdot=st.qfddot=std::vector<double>(1,0);
        md.nModes=1;md.nPoints=9;md.frequencies={100};md.modes.assign(1,std::vector<Displacement>(9));
        for(auto& u:md.modes[0])u.uy=side==0?1:-1;
    }
    s.omegaL=s.omegaR={2*M_PI*100};s.fCalc.initialize();s.fCalc.updateChannelSections();
    for(double a:s.fCalc.harea)expect(std::abs(a-4)<1e-12);
    s.fCalc.applyFluidLoads(0,0); // Establish the real separation/pressure metadata.
    expect(setenv("SIMULATION_KANAYA_DIAGNOSTICS","1",1)==0);
    KanayaDiagnostics diag(s);
    for(int n=0;n<3;++n) {
        diag.beginStep(n,n*s.params.dt);
        const double pressure=n==1?-100:100;
        for(int side=0;side<2;++side){bool left=side==0;double ay=left?-1e-6:1e-6;
            auto& fy=left?s.fCalc.fyL:s.fCalc.fyR;fy[1][1]=pressure*ay;
            KanayaDiagnostics::node(left,1,1,pressure,0,ay,0,0,fy[1][1],0);}
        diag.afterFluid();s.fCalc.projectLoadsToModes();diag.adoptedContact();
        for(int side=0;side<2;++side){auto& st=side==0?s.stateL:s.stateR;auto& g=side==0?s.geomL:s.geomR;auto& md=side==0?s.mdataL:s.mdataR;
            st.qf[0]=.001*(n+1);st.qfdot[0]=.1;
            st.mode2uf(g,md,n+1);}
        diag.endStep(0,2);s.stateL.uf2u();s.stateR.uf2u();
    }
    diag.finish();
    std::ifstream meta(s.runDir/"kanaya_diagnostic/run_metadata.json");
    std::string json((std::istreambuf_iterator<char>(meta)),{});
    expect(json.find("\"integrated_intervals\": 2")!=std::string::npos);
    expect(json.find("\"excluded_intervals\": 1")!=std::string::npos);
    std::cout<<"PASS constant geometry, mirrored +/- pressure, actual nodal/modal work, sign reversal, tolerance rejection, final-only exclusion\n";
#endif
}
