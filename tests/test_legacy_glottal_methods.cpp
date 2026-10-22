#include "LegacyGlottalMethods.h"
#include "Simulation.h"
#include "legacy_snapshot_fixture.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

void require(bool test,const char* text) { if(!test) throw std::runtime_error(text); }
void near(double a,double b,double abs=1e-12,double rel=1e-11) {
    if(!std::isfinite(a+b)||std::abs(a-b)>abs+rel*std::max(std::abs(a),std::abs(b)))
        throw std::runtime_error("legacy numeric assertion failed: "+std::to_string(a)+" versus "+std::to_string(b));
}

void mappingSnapshots(std::ostream& out) {
    out << std::setprecision(17)
        << "case,experiment,method,area_m2,centroid_x_mm,centroid_z_mm,Fx_N,Fy_N,Fz_N,Q0,Q1,W\n";
    for (int kind=0; kind<6; ++kind) {
        Geometry l,r; State ls,rs; ModeData md;
        legacyFixture(kind,l,ls,md); r=l; rs.disp=ls.disp;
        for (auto& p:r.points) p.y=-p.y;
        for (auto& p:rs.disp) p.uy=-p.uy;
        ChannelSections sec; sec.resize(l,r,7);
        LegacyGlottal::buildArea(l,ls,r,rs,sec);
        GlottalMethodOptions old,modernOnRows;
        old.configure("legacy_fluid_all",{});
        modernOnRows.configure("custom",{{"legacyarea",true}});
        const auto oldP=LegacyGlottal::wallPressure(sec,l,old,
            LegacyGlottal::select(sec,old,1.),1.13,1.83e-5,3e-4,1000.,1000.,20.);
        const auto newP=LegacyGlottal::wallPressure(sec,l,modernOnRows,
            LegacyGlottal::select(sec,modernOnRows,1.),1.13,1.83e-5,3e-4,1000.,1000.,20.);
        // Two same-pressure experiments separate the area-vector formula
        // from the historical material-row load support. The third compares
        // each pressure rule on the SAME material-row geometry and Ug/Pg;
        // it is intentionally not claimed to be the 50-plane modern model.
        for (int experiment=0; experiment<3; ++experiment) {
            for (int method=0; method<2; ++method) {
                const bool legacy=method==0;
                double area=0.,cx=0.,cz=0.,fx=0.,fy=0.,fz=0.;
                for (int i=1; i<6; ++i) for (int j=1; j<4; ++j) {
                    if (experiment>0 && legacy && i>=oldP.loadEndExclusive) continue;
                    const auto a=LegacyGlottal::areaVector(l,ls,true,i,j,legacy);
                    const auto pnt=ls.disp[l.surfp[i][j]];
                    const double magnitude=std::sqrt(a.ux*a.ux+a.uy*a.uy+a.uz*a.uz);
                    const double p=experiment<2?1000.:(legacy?oldP.pressure[i]:
                        LegacyGlottal::interpolatePressure(sec,newP.pressure,newP.sepIndex,pnt.ux,.5,20.,false));
                    area+=magnitude; cx+=magnitude*pnt.ux; cz+=magnitude*pnt.uz;
                    fx+=p*a.ux; fy+=p*a.uy; fz+=p*a.uz;
                }
                const double q0=fy,q1=.3*fx-.1*fy+.2*fz;
                out << kind << ',' << experiment << ',' << (legacy?"legacy":"modern_on_rows")
                    << ',' << area << ',' << (area>0?cx/area:0.) << ',' << (area>0?cz/area:0.)
                    << ',' << fx << ',' << fy << ',' << fz << ',' << q0 << ',' << q1
                    << ',' << q0*1e-6-q1*2e-6 << '\n';
            }
        }
    }
}

void snapshots(std::ostream& out) {
    out<<std::setprecision(17)<<"case,i,j,area,pressure,Fx,Fy,Fz,Q0,Q1,W\n";
    for(int kind=0;kind<6;++kind) {
        Geometry l,r; State ls,rs; ModeData md;
        legacyFixture(kind,l,ls,md); r=l; rs.disp=ls.disp;
        for(auto& p:r.points)p.y=-p.y;
        for(auto& p:rs.disp)p.uy=-p.uy;
        ChannelSections sec; sec.resize(l,r,7);
        LegacyGlottal::buildArea(l,ls,r,rs,sec);
        GlottalMethodOptions options; options.configure("legacy_fluid_all",{});
        const auto selection=LegacyGlottal::select(sec,options,1.);
        const auto pressure=LegacyGlottal::wallPressure(sec,l,options,selection,1.13,1.83e-5,3e-4,1000.,1000.,20.);
        double totalX=0,totalY=0,totalZ=0,q0=0,q1=0;
        for(int i=0;i<7;++i)for(int j=0;j<5;++j) {
            double fx=0,fy=0,fz=0;
            if(i>0 && i<pressure.loadEndExclusive && i<6 && j>0 && j<4) {
                const auto a=LegacyGlottal::areaVector(l,ls,true,i,j,true);
                const auto b=LegacyGlottal::areaVector(r,rs,false,i,j,true);
                near(a.ux,b.ux); near(a.uy,-b.uy); near(a.uz,b.uz);
                fx=pressure.pressure[i]*a.ux; fy=pressure.pressure[i]*a.uy; fz=pressure.pressure[i]*a.uz;
            }
            const double qi=fy, qj=.3*fx-.1*fy+.2*fz;
            out<<kind<<','<<i<<','<<j<<','<<sec.area[i]<<','<<pressure.pressure[i]<<','
               <<fx<<','<<fy<<','<<fz<<','<<qi<<','<<qj<<','<<qi*1e-6-qj*2e-6<<'\n';
            totalX+=fx;totalY+=fy;totalZ+=fz;q0+=qi;q1+=qj;
        }
        out<<kind<<",-1,-1,0,0,"<<totalX<<','<<totalY<<','<<totalZ<<','<<q0<<','<<q1<<','<<q0*1e-6-q1*2e-6<<'\n';
        if(kind==0) for(double a:sec.area)near(a,17.*.6);
        if(kind==5) for(double a:sec.area)near(a,0.);
        // Reject asymmetric geometry instead of silently applying a symmetric model.
        rs.disp[r.surfp[2][2]].ux+=.01;
        bool rejected=false;
        try{LegacyGlottal::buildArea(l,ls,r,rs,sec);}catch(const std::runtime_error&){rejected=true;}
        require(rejected,"asymmetric legacy area must be rejected");
    }
    TimeIntegrator integrator;
    for(int k=0;k<16;++k) {
        double q,v,a;
        integrator.newmarkStep(k%2?0.:.2,std::sin(k)*1e-5,std::cos(k)*1e-3,k*.001,
            1e-5,(45.+k*10.)*2*M_PI,.06,.275625,.55,q,v,a);
        out<<100+k<<",-2,-2,0,0,"<<q<<','<<v<<','<<a<<",0,0,0\n";
    }
}

int main(int argc,char** argv) {
    try {
        GlottalMethodOptions options;
        require(!options.any(),"modern must remain default");
        for(const std::string flag:{"legacycontact","legacyintegrator"}) {
            bool rejected=false;
            try{options.configure("custom",{{flag,true}});}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"unimplemented flag silently accepted");
        }
        for(const std::string profile:{"nonsense","modern","custom"}) {
            bool rejected=false;
            try{options.configure(profile,{{"legacysurfacepressureload",true}});}catch(const std::runtime_error&){rejected=true;}
            require(rejected,"invalid profile combination accepted");
        }
        options.configure("legacy_fluid_all",{});
        require(options.area&&options.flowCoefficients&&options.wallViscosity&&options.separation
                &&options.postSeparation&&options.surfacePressureLoad,"incomplete legacy profile");
        const auto a=LegacyGlottal::coefficients(1.13,1.83e-5,.3175,.012,.017,1e-5);
        const auto b=LegacyGlottal::coefficients(1.13,1.83e-5,.3175,.012,.017,2e-5);
        near(a.Lg,1.13*.5*.012/1e-5); near(a.Rv,12.*1.83e-5*.017*.017*.012/1e-15,1e-6);
        near(a.Lg/b.Lg,2.); near(a.Rv/b.Rv,8.); near(a.Rk/b.Rk,4.);
        Geometry g; State st; ModeData md; legacyFixture(0,g,st,md);
        auto normal=LegacyGlottal::areaVector(g,st,true,2,2,true);
        near(normal.ux,0.); near(normal.uy,-2.*4.25e-6); near(normal.uz,0.);
        const auto modernFlat=LegacyGlottal::areaVector(g,st,true,2,2,false);
        near(normal.ux,modernFlat.ux); near(normal.uy,modernFlat.uy); near(normal.uz,modernFlat.uz);
        // Known affine surface: test atan angles and projected area, not only mirror signs.
        for(auto& p:st.disp)p.uy=.2*p.ux+.1*p.uz;
        normal=LegacyGlottal::areaVector(g,st,true,2,2,true);
        const double ds=std::hypot(2.,.4), tz=std::atan(.1), tx=std::atan(.2);
        near(normal.ux,ds*4.25e-6*std::cos(tz)*std::sin(tx));
        near(normal.uy,-ds*4.25e-6*std::cos(tz)*std::cos(tx));
        near(normal.uz,ds*4.25e-6*std::sin(tz));
        const auto modernSlope=LegacyGlottal::areaVector(g,st,true,2,2,false);
        near(modernSlope.ux,.4*4.25e-6);
        near(modernSlope.uy,-2.*4.25e-6);
        near(modernSlope.uz,2.*.425e-6);
        // A known linear p(x) tests physical interpolation independently of
        // pressure formation or the area-vector formulas.
        ChannelSections line; line.x={0.,2.,4.};
        near(LegacyGlottal::interpolatePressure(line,{2.,8.,14.},2,1.,.5,20.,false),5.);
        const double next=std::nextafter(4.,std::numeric_limits<double>::infinity());
        near(LegacyGlottal::interpolatePressure(line,{2.,8.,14.},2,next,.5,20.,true),14.);
        near(LegacyGlottal::interpolatePressure(line,{2.,8.,14.},2,4.+1e-10,.5,20.,true),0.);
        // Existing endpoint behavior is detected but intentionally not changed in this task.
        st.nSteps=1; st.nPoints=35; st.predictedDisp.assign(35,{}); st.qf={1.,0.}; st.qfdot={0.,0.};
        st.mode2uf(g,md,1); near(st.predictedDisp[0].ux,0.);
        if(argc==3&&(std::string(argv[1])=="--snapshots"||std::string(argv[1])=="--mapping")) {
            if(std::filesystem::exists(argv[2]))throw std::runtime_error("snapshot already exists");
            std::ofstream f(argv[2]);
            if (std::string(argv[1])=="--mapping") mappingSnapshots(f); else snapshots(f);
        } else snapshots(std::cout);
        std::cout<<"PASS legacy profiles, area, coefficient scaling, pressure, mirrored loads, work, endpoint guard\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
