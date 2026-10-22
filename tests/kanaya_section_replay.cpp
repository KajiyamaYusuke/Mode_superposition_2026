// Offline geometry only: directly compile the unchanged production section
// implementation so its private intersection/segment rules are reused exactly.
#include "../src/fluid/channelSections.cpp"
#include "Simulation.h"
#include "SeparationSelection.h"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {
    try {
        if (argc!=5) return 2;
        Simulation sim;sim.initialize(argv[1]);
        if (sim.mdataL.nModes!=20 || sim.mdataR.nModes!=20) throw std::runtime_error("all 20 modes required");
        std::ifstream binary(argv[2],std::ios::binary),requests(argv[3]);
        const fs::path directory=argv[4];fs::create_directories(directory);
        std::ofstream sections(directory/"replay_sections.csv"),intervals(directory/"replay_intervals.csv"),
            checks(directory/"replay_checks.csv"),coords(directory/"control_coordinates.csv");
        sections<<"event,step,section,x_mm,valid,area_mm2,gap_cubed_mm4\n"<<std::setprecision(17);
        intervals<<"event,step,section,x_mm,piece,z0_mm,z1_mm,yL0_mm,yL1_mm,yR0_mm,yR1_mm,h0_mm,h1_mm\n"<<std::setprecision(17);
        checks<<"event,step,min_index,x_min_mm,Amin_mm2,probe_uxL_mm,probe_uyL_mm,probe_uxR_mm,probe_uyR_mm,reference_z_min_mm,reference_z_max_mm,probe_id\n"<<std::setprecision(17);
        coords<<"side,id,x_mm,y_mm,z_mm\n"<<std::setprecision(17);
        double distance=1e100,zMin=1e100,zMax=-1e100;int probe=-1;
        for (int i=0;i<sim.geomL.nxsup;++i) for (int j=0;j<sim.geomL.nsurfz;++j) {
            const int id=sim.geomL.surfp[i][j];if(id<0)continue;
            const auto& p=sim.geomL.points[id];
            zMin=std::min(zMin,p.z);zMax=std::max(zMax,p.z);
            const double d=(p.x-9.2)*(p.x-9.2)+(p.z-8.5)*(p.z-8.5);
            if(d<distance){distance=d;probe=id;}
        }
        int event,step;
        while (requests>>event>>step) {
            if(step<0 || step>=19999)throw std::runtime_error("excluded/native-missing step requested");
            binary.seekg(std::streamoff(step)*2*20*7*8);
            for (int side=0;side<2;++side) {
                auto& st=side==0?sim.stateL:sim.stateR;
                const auto& geom=side==0?sim.geomL:sim.geomR;
                const auto& md=side==0?sim.mdataL:sim.mdataR;
                st.qf.resize(20);st.qfdot.resize(20);
                for (int m=0;m<20;++m) {
                    std::array<double,7> row{};binary.read(reinterpret_cast<char*>(row.data()),sizeof(row));
                    if(!binary || !std::isfinite(row[4]+row[5]))throw std::runtime_error("incomplete native modal history");
                    st.qf[m]=row[4];st.qfdot[m]=row[5];
                }
                st.mode2uf(geom,md,step);st.disp=st.predictedDisp;
                if(event==99) for(int id=0;id<geom.nPoints;++id) {
                    const auto& p=st.disp[id];coords<<(side==0?'L':'R')<<','<<id<<','<<p.ux<<','<<p.uy<<','<<p.uz<<'\n';
                }
            }
            ChannelSections sec;sec.resize(sim.geomL,sim.geomR,50);
            ChannelSectionBuilder builder;builder.build(sim.geomL,sim.stateL,sim.geomR,sim.stateR,sec);
            const auto selection=selectSeparationByAreaRatio(sec.area,sec.valid,1.,.01);
            const auto& l=sim.stateL.disp[probe];const auto& r=sim.stateR.disp[probe];
            const auto& gl=sim.geomL.points[probe];const auto& gr=sim.geomR.points[probe];
            checks<<event<<','<<step<<','<<selection.minIndex<<','<<sec.x[selection.minIndex]<<','<<selection.minAreaMm2<<','
                <<l.ux-gl.x<<','<<l.uy-gl.y<<','<<r.ux-gr.x<<','<<r.uy-gr.y<<','<<zMin<<','<<zMax<<','<<probe<<'\n';
            for (int k=0;k<50;++k) {
                sections<<event<<','<<step<<','<<k<<','<<sec.x[k]<<','<<sec.valid[k]<<','<<sec.area[k]<<','<<sec.gapCubed[k]<<'\n';
                const auto left=sectionSegments(sim.geomL,sim.stateL,sec.x[k],sec.gapSign>0.,"left");
                const auto right=sectionSegments(sim.geomR,sim.stateR,sec.x[k],sec.gapSign<0.,"right");
                std::size_t il=0,ir=0;int piece=0;
                while(il<left.size() && ir<right.size()) {
                    const auto& a=left[il];const auto& b=right[ir];
                    const double low=std::max(a.a.z,b.a.z),high=std::min(a.b.z,b.b.z);
                    if(high>low) {
                        const double l0=interpolateYAtZ(a,low),l1=interpolateYAtZ(a,high);
                        const double r0=interpolateYAtZ(b,low),r1=interpolateYAtZ(b,high);
                        intervals<<event<<','<<step<<','<<k<<','<<sec.x[k]<<','<<piece++<<','<<low<<','<<high<<','
                            <<l0<<','<<l1<<','<<r0<<','<<r1<<','<<sec.gapSign*(r0-l0)<<','<<sec.gapSign*(r1-l1)<<'\n';
                    }
                    if(a.b.z<=b.b.z)++il;
                    if(b.b.z<=a.b.z)++ir;
                }
            }
        }
        if(!sections || !intervals || !checks || !coords)throw std::runtime_error("replay output error");
        std::cout<<"PASS offline all-mode geometry replay; no flow/structural time integration invoked\n";
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
