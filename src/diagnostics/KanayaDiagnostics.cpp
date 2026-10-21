#include "KanayaDiagnostics.h"
#include "Simulation.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>

namespace {
constexpr double kMissing = std::numeric_limits<double>::quiet_NaN();
struct Vec { double x=0, y=0, z=0; };
struct Node {
    int i, j, id;
    double p, area, x, z;
    Vec force;
};
struct Modal {
    double pos=0, neg=0, fluid=0, damping=0, contact=0, q2=0, v2=0, maxq=0;
    int count=0;
};
struct Contact {
    int pairs=0, excluded=0;
    double fx=0, fy=0, x=0, z=0, penetration=0;
};
struct Region {
    Vec force;
    double area=0, xmin=std::numeric_limits<double>::infinity(), xmax=-std::numeric_limits<double>::infinity();
    double pmin=std::numeric_limits<double>::infinity(), pmax=-std::numeric_limits<double>::infinity();
    int count=0, wrongNormal=0;
};
struct ForceStats {
    double sum=0,area=0,min=std::numeric_limits<double>::infinity(),max=-std::numeric_limits<double>::infinity(),minTime=0,maxTime=0;
    int count=0;
};
int sign(double p) { return p>0 ? 1 : p<0 ? -1 : 0; }
bool apEnd(int j, int nj) { return j<3 || j>=nj-3; }
double dot(const Vec& f, const Displacement& u) { return f.x*u.ux+f.y*u.uy+f.z*u.uz; }
void requireFinite(double x) { if(!std::isfinite(x)) throw std::runtime_error("Kanaya: nonfinite observation"); }
std::ofstream output(const fs::path& path, bool binary=false) {
    if(fs::exists(path)) throw std::runtime_error("Kanaya refuses to overwrite existing diagnostic: "+path.string());
    std::ofstream f(path, binary ? std::ios::binary : std::ios::out);
    if(!f) throw std::runtime_error("Kanaya cannot open "+path.string());
    f<<std::setprecision(17); return f;
}
}

struct KanayaDiagnostics::Data {
    Simulation& sim;
    fs::path dir;
    int n=-1, calls=0, integrated=0, excluded=0, iteration=0;
    double time=0, gapSign=1, lg=kMissing, rv=kMissing, rk=kMissing;
    bool coefficientsActive=false, contactDiagnosticOnly=false, fluidReady=false, contactReady=false;
    std::array<std::vector<Node>,2> nodes;
    std::array<std::vector<unsigned char>,2> seen;
    std::array<std::vector<Vec>,2> snapshot;
    std::array<std::vector<Vec>,2> velocitySnapshot;
    std::array<std::vector<double>,2> qpos,qneg,qfluid,qcontact;
    std::array<std::array<std::vector<Modal>,4>,2> sums;
    std::array<std::array<std::array<double,3>,4>,2> cumulativeNodal{};
    std::array<std::array<Contact,2>,2> contacts;
    std::array<std::array<std::array<ForceStats,2>,4>,2> forceStats;
    std::array<double,2> forceResidual{}, qResidual{}, forceScale{}, qScale{};
    std::ofstream regions,summary,consistency,contactFile,flow,raw;
    const Geometry& geom(int side) const { return side==0?sim.geomL:sim.geomR; }
    const ModeData& modes(int side) const { return side==0?sim.mdataL:sim.mdataR; }
    const State& state(int side) const { return side==0?sim.stateL:sim.stateR; }
    const std::vector<std::vector<double>>& fx(int s) const {return s==0?sim.fCalc.fxL:sim.fCalc.fxR;}
    const std::vector<std::vector<double>>& fy(int s) const {return s==0?sim.fCalc.fyL:sim.fCalc.fyR;}
    const std::vector<std::vector<double>>& fz(int s) const {return s==0?sim.fCalc.fzL:sim.fCalc.fzR;}
    explicit Data(Simulation& s):sim(s),dir(s.runDir/"kanaya_diagnostic") {
        if(sim.params.iforce!=0) throw std::runtime_error("Kanaya observes iforce=0 only");
        fs::create_directories(dir);
        gapSign=sim.fCalc.diagnosticSections().gapSign;
        regions=output(dir/"force_regions.csv"); summary=output(dir/"modal_energy_summary.csv");
        consistency=output(dir/"consistency.csv"); contactFile=output(dir/"contact_regions.csv");
        flow=output(dir/"flow_geometry.csv"); raw=output(dir/"modal_native.bin",true);
        regions<<"step,time_s,side,x_bin_mm,zone,ap_region,pressure_sign,Fx_N,Fy_N,Fz_N,opening_Fy_N,applied_node_count,area_m2,x_min_mm,x_max_mm,p_min_Pa,p_max_Pa,wrong_normal_count\n";
        summary<<"side,mode,source_mode,f_Hz,window,W_pos,W_neg,W_fluid,D_struct,W_contact_reference,q_rms,qdot_rms,q_max,interval_count,start_s,end_s\n";
        consistency<<"step,time_s,side,included,force_abs,force_scale,Q_abs,Q_scale,work_split_abs,nodal_work_abs,work_scale,contact_Q_abs,scaled_max,valid\n";
        contactFile<<"step,time_s,side,ap_region,adopted_iteration,pairs,excluded_pairs,Fx_N,Fy_N,opening_Fy_N,x_mean_mm,z_mean_mm,max_penetration_m,contact_stage\n";
        flow<<"step,time_s,interval_included,Amin_mm2,xMin_mm,xSep_mm,xBlendEnd_mm,Pg_Pa,Psep_Pa,Pd0_Pa,Ug_m3_s,Lg,Rv,Rk,coefficients_active,contact_iterations,stop_reason,contact_flag,force_residual,penetration_residual,monitor_x_mm,monitor_z_mm,uxL_mm,uyL_mm,uxR_mm,uyR_mm,monitor_gap_mm,mean_delta_gap_mm,max_delta_gap_mm,max_delta_gap_x_mm,sep_node_x_mm,sep_node_uyL_mm,sep_node_uyR_mm,sep_node_vyL_mm_s,sep_node_vyR_mm_s\n";
        for(int side=0;side<2;++side) {
            const int nm=modes(side).nModes;
            for(auto* a:{&qpos[side],&qneg[side],&qfluid[side],&qcontact[side]}) a->assign(nm,0);
            for(auto& window:sums[side])window.resize(nm);
            seen[side].resize(geom(side).nsurfl*geom(side).nsurfz);
            snapshot[side].resize(seen[side].size()); nodes[side].reserve(seen[side].size());
            velocitySnapshot[side].resize(seen[side].size());
        }
    }
    void check(double difference,double scale,const char* what) const {
        if(!std::isfinite(difference)||!std::isfinite(scale)||KanayaDiagnostics::scaledResidual(difference,scale)>1)
            throw std::runtime_error("Kanaya consistency failure step="+std::to_string(n)+" "+what+
                                     " difference="+std::to_string(difference));
    }
    void saveRegions(int side) {
        std::map<std::tuple<int,int,int,int>,Region> bins;
        for(const auto& node:nodes[side]) {
            int zone=node.x<=sim.fCalc.separationX()?0:node.x<sim.fCalc.separationBlendEndX()?1:2;
            auto& r=bins[{int(std::floor(node.x/.5)),zone,apEnd(node.j,geom(side).nsurfz)?0:1,sign(node.p)}];
            r.force.x+=node.force.x;r.force.y+=node.force.y;r.force.z+=node.force.z;
            r.area+=node.area;++r.count;r.xmin=std::min(r.xmin,node.x);r.xmax=std::max(r.xmax,node.x);
            r.pmin=std::min(r.pmin,node.p);r.pmax=std::max(r.pmax,node.p);
            if(node.p!=0 && KanayaDiagnostics::openingForce(gapSign,side==0,node.force.y)/node.p < -1e-18)++r.wrongNormal;
        }
        for(const auto& [key,r]:bins) {
            auto [bin,zone,ap,sgn]=key;
            regions<<n<<','<<time<<','<<(side==0?'L':'R')<<','<<bin*.5<<','
                   <<(zone==0?"upstream":zone==1?"blend":"downstream")<<','<<(ap==0?"end":"interior")<<','<<sgn<<','
                   <<r.force.x<<','<<r.force.y<<','<<r.force.z<<','<<KanayaDiagnostics::openingForce(gapSign,side==0,r.force.y)<<','
                   <<r.count<<','<<r.area<<','<<r.xmin<<','<<r.xmax<<','<<r.pmin<<','<<r.pmax<<','<<r.wrongNormal<<'\n';
        }
    }
    void saveFlow(int iterations,int stopReason) {
        const auto& l=sim.geomL;const auto& r=sim.geomR;
        double monitorDistance=1e100,sepDistance=1e100,mean=0,maxGap=-1e100,maxX=kMissing;
        int mi=-1,mj=-1,si=-1,sj=-1,count=0;
        for(int i=0;i<std::min(l.nxsup,r.nxsup);++i)for(int j=0;j<std::min(l.nsurfz,r.nsurfz);++j) {
            int il=l.surfp[i][j],ir=r.surfp[i][j];if(il<0||ir<0)continue;
            const auto& p=sim.stateL.disp[il];
            double d=(l.points[il].x-9.2)*(l.points[il].x-9.2)+(l.points[il].z-8.5)*(l.points[il].z-8.5);
            if(d<monitorDistance){monitorDistance=d;mi=i;mj=j;}
            d=(p.ux-sim.fCalc.separationX())*(p.ux-sim.fCalc.separationX())+(p.uz-8.5)*(p.uz-8.5);
            if(d<sepDistance){sepDistance=d;si=i;sj=j;}
            double dg=gapSign*((sim.stateR.disp[ir].uy-p.uy)-(r.points[ir].y-l.points[il].y));
            mean+=dg;++count;if(dg>maxGap){maxGap=dg;maxX=p.ux;}
        }
        if(mi<0||si<0)throw std::runtime_error("Kanaya missing paired probes");
        int il=l.surfp[mi][mj],ir=r.surfp[mi][mj],sl=l.surfp[si][sj],sr=r.surfp[si][sj];
        auto pl=sim.stateL.disp[il],pr=sim.stateR.disp[ir];
        const auto& f=sim.fCalc;const auto& sec=f.diagnosticSections();
        const auto select=selectSeparationByAreaRatio(sec.area,sec.valid,sim.params.flowSeparationAreaRatio,.01);
        flow<<n<<','<<time<<','<<KanayaDiagnostics::includedInterval(n,sim.params.nstep)<<','<<f.separationMinArea()<<','<<sec.x.at(select.minIndex)<<','
            <<f.separationX()<<','<<f.separationBlendEndX()<<','<<f.currentPg<<','<<f.separationPressure()<<','<<f.diagnosticFlow().downstreamPressure().front()<<','
            <<f.currentUg<<','<<(coefficientsActive?lg:kMissing)<<','<<(coefficientsActive?rv:kMissing)<<','<<(coefficientsActive?rk:kMissing)<<','<<coefficientsActive<<','
            <<iterations<<','<<stopReason<<','<<f.contactFlag<<','<<f.contact_force_residual<<','<<f.contact_penetration_residual<<','
            <<pl.ux<<','<<pl.uz<<','<<pl.ux-l.points[il].x<<','<<pl.uy-l.points[il].y<<','<<pr.ux-r.points[ir].x<<','<<pr.uy-r.points[ir].y<<','
            <<gapSign*(pr.uy-pl.uy)<<','<<mean/count<<','<<maxGap<<','<<maxX<<','<<sim.stateL.disp[sl].ux<<','
            <<sim.stateL.disp[sl].uy-l.points[sl].y<<','<<sim.stateR.disp[sr].uy-r.points[sr].y<<','<<velocitySnapshot[0][si*l.nsurfz+sj].y<<','<<velocitySnapshot[1][si*r.nsurfz+sj].y<<'\n';
    }
};

KanayaDiagnostics::Data* KanayaDiagnostics::active_=nullptr;
KanayaDiagnostics::KanayaDiagnostics(Simulation& s) {
    const char* flag=std::getenv("SIMULATION_KANAYA_DIAGNOSTICS");
    if(flag && std::string(flag)=="1") {
        if(active_)throw std::runtime_error("Kanaya duplicate observer");
        data_=std::make_unique<Data>(s);active_=data_.get();
    }
}
KanayaDiagnostics::~KanayaDiagnostics(){if(data_)active_=nullptr;}
double KanayaDiagnostics::openingForce(double gapSign,bool left,double fy){return gapSign*(left?-fy:fy);}
double KanayaDiagnostics::scaledResidual(double difference,double scale,double absoluteFloor){return std::abs(difference)/(absoluteFloor+1e-8*std::abs(scale));}
bool KanayaDiagnostics::includedInterval(int n,int count){return n>=0 && n<count-1;}

void KanayaDiagnostics::beginStep(int n,double t) {
    if(!data_)return;auto& d=*data_;
    if(n!=d.calls || t!=n*d.sim.params.dt)throw std::runtime_error("Kanaya nonsequential time index");
    ++d.calls;d.n=n;d.time=t;d.fluidReady=d.contactReady=false;
    d.contacts={};d.iteration=0;d.coefficientsActive=false;
    for(int side=0;side<2;++side){d.nodes[side].clear();std::fill(d.seen[side].begin(),d.seen[side].end(),0);
        const auto& g=d.geom(side);for(int i=0;i<g.nsurfl;++i)for(int j=0;j<g.nsurfz;++j){int id=g.surfp[i][j];if(id>=0){auto v=d.state(side).vel[id];d.velocitySnapshot[side][i*g.nsurfz+j]={v.ux,v.uy,v.uz};}}}
}
void KanayaDiagnostics::node(bool left,int i,int j,double p,double ax,double ay,double az,double fx,double fy,double fz) {
    if(!active_)return;auto& d=*active_;int s=left?0:1;const auto& g=d.geom(s);const auto& st=d.state(s);
    if(d.n<0||i<0||i>=g.nsurfl||j<0||j>=g.nsurfz)throw std::runtime_error("Kanaya invalid node index");
    const int k=i*g.nsurfz+j,id=g.surfp[i][j];if(id<0||d.seen[s][k]++)throw std::runtime_error("Kanaya duplicate/missing node observation");
    for(double a:{p,ax,ay,az,fx,fy,fz,st.disp[id].ux,st.disp[id].uz})requireFinite(a);
    d.check(fx-p*ax,std::abs(fx),"Fx != p*dA");d.check(fy-p*ay,std::abs(fy),"Fy != p*dA");d.check(fz-p*az,std::abs(fz),"Fz != p*dA");
    if(p==0 && (fx!=0||fy!=0||fz!=0))throw std::runtime_error("Kanaya zero pressure with nonzero load");
    d.nodes[s].push_back({i,j,id,p,std::sqrt(ax*ax+ay*ay+az*az),st.disp[id].ux,st.disp[id].uz,{fx,fy,fz}});
}
void KanayaDiagnostics::coefficients(double lg,double rv,double rk,bool on) {
    if(!active_)return;active_->lg=lg;active_->rv=rv;active_->rk=rk;active_->coefficientsActive=on;
    if(on){requireFinite(lg);requireFinite(rv);requireFinite(rk);}
}
void KanayaDiagnostics::beginContact(bool diagnosticOnly,int iteration) {
    if(!active_)return;auto& d=*active_;d.contactDiagnosticOnly=diagnosticOnly;
    if(!diagnosticOnly && includedInterval(d.n,d.sim.params.nstep)){d.contacts={};d.iteration=iteration;}
}
void KanayaDiagnostics::contact(int j,double x,double z,double pen,double fxl,double fyl,double fxr,double fyr,bool excluded) {
    if(!active_||active_->contactDiagnosticOnly||!includedInterval(active_->n,active_->sim.params.nstep))return;
    auto& d=*active_;for(double a:{x,z,pen,fxl,fyl,fxr,fyr})requireFinite(a);
    for(int s=0;s<2;++s){auto& r=d.contacts[s][apEnd(j,d.geom(s).nsurfz)?0:1];++r.pairs;r.excluded+=excluded;
        r.x+=x;r.z+=z;r.penetration=std::max(r.penetration,pen);if(!excluded){r.fx+=s==0?fxl:fxr;r.fy+=s==0?fyl:fyr;}}
}

void KanayaDiagnostics::afterFluid() {
    if(!data_)return;auto& d=*data_;if(d.fluidReady)throw std::runtime_error("Kanaya fluid double count");
    for(int s=0;s<2;++s) {
        const auto& g=d.geom(s);std::array<Vec,3> decomposed{};Vec total{};double scale=0;
        for(int i=0;i<g.nsurfl;++i)for(int j=0;j<g.nsurfz;++j){int k=i*g.nsurfz+j;
            Vec f{d.fx(s)[i][j],d.fy(s)[i][j],d.fz(s)[i][j]};d.snapshot[s][k]=f;
            if(!d.seen[s][k]&&(f.x!=0||f.y!=0||f.z!=0))throw std::runtime_error("Kanaya load without pressure observation");
            total.x+=f.x;total.y+=f.y;total.z+=f.z;scale+=std::abs(f.x)+std::abs(f.y)+std::abs(f.z);}
        for(const auto& p:d.nodes[s]){auto& f=decomposed[sign(p.p)+1];f.x+=p.force.x;f.y+=p.force.y;f.z+=p.force.z;}
        d.forceResidual[s]=std::max({std::abs(total.x-decomposed[0].x-decomposed[2].x),std::abs(total.y-decomposed[0].y-decomposed[2].y),std::abs(total.z-decomposed[0].z-decomposed[2].z)});
        d.forceScale[s]=scale;d.check(d.forceResidual[s],scale,"force split");
        if(includedInterval(d.n,d.sim.params.nstep))for(int sg=0;sg<2;++sg){
            double opening=openingForce(d.gapSign,s==0,decomposed[sg==0?0:2].y),area=0;
            for(const auto& p:d.nodes[s])if(sign(p.p)==(sg==0?-1:1))area+=p.area;
            for(int w:{0,d.time<.05?1:d.time<.1?2:3}){auto& stats=d.forceStats[s][w][sg];stats.sum+=opening;stats.area+=area;++stats.count;
                if(opening<stats.min){stats.min=opening;stats.minTime=d.time;}if(opening>stats.max){stats.max=opening;stats.maxTime=d.time;}}
        }
        if(d.n%20==0)d.saveRegions(s);
        if(!includedInterval(d.n,d.sim.params.nstep))continue;
        const auto& md=d.modes(s);
        #pragma omp parallel for schedule(static)
        for(int m=0;m<md.nModes;++m){double pos=0,neg=0,all=0;
            for(const auto& p:d.nodes[s]){double v=dot(p.force,md.modes[m][p.id]);if(p.p>0)pos+=v;else if(p.p<0)neg+=v;}
            for(int i=0;i<g.nsurfl;++i)for(int j=0;j<g.nsurfz;++j){int id=g.surfp[i][j];if(id<0)continue;all+=dot(d.snapshot[s][i*g.nsurfz+j],md.modes[m][id]);}
            d.qpos[s][m]=pos;d.qneg[s][m]=neg;d.qfluid[s][m]=all;}
        double res=0,qs=0;for(int m=0;m<md.nModes;++m){res=std::max(res,std::abs(d.qpos[s][m]+d.qneg[s][m]-d.qfluid[s][m]));qs=std::max(qs,std::abs(d.qpos[s][m])+std::abs(d.qneg[s][m]));}
        d.qResidual[s]=res;d.qScale[s]=qs;d.check(res,qs,"Q split");
    }
    d.fluidReady=true;
}

void KanayaDiagnostics::adoptedContact() {
    if(!data_)return;auto& d=*data_;if(!d.fluidReady||d.contactReady)throw std::runtime_error("Kanaya contact stage ordering");
    if(includedInterval(d.n,d.sim.params.nstep))for(int s=0;s<2;++s){const auto& g=d.geom(s);const auto& md=d.modes(s);
        const auto& cx=s==0?d.sim.fCalc.prevFdisXL:d.sim.fCalc.prevFdisXR;
        const auto& cy=s==0?d.sim.fCalc.prevFdisYL:d.sim.fCalc.prevFdisYR;
        #pragma omp parallel for schedule(static)
        for(int m=0;m<md.nModes;++m){double q=0;for(int i=0;i<int(cx.size());++i)for(int j=0;j<int(cx[i].size());++j){int id=g.surfp[i][j];if(id<0)continue;q+=cx[i][j]*md.modes[m][id].ux+cy[i][j]*md.modes[m][id].uy;}d.qcontact[s][m]=q;}
    }
    d.contactReady=true;
}

void KanayaDiagnostics::endStep(int iterations,int stopReason) {
    if(!data_)return;auto& d=*data_;if(!d.fluidReady||!d.contactReady)throw std::runtime_error("Kanaya incomplete observation stage");
    d.saveFlow(iterations,stopReason);
    if(!includedInterval(d.n,d.sim.params.nstep)) {
        ++d.excluded;for(int s=0;s<2;++s)d.consistency<<d.n<<','<<d.time<<','<<(s==0?'L':'R')<<",0,"<<d.forceResidual[s]<<','<<d.forceScale[s]<<",nan,nan,nan,nan,nan,nan,nan,0\n";
        return;
    }
    ++d.integrated;int window=d.time<.05?1:d.time<.1?2:3;
    for(int s=0;s<2;++s) {
        const auto& st=d.state(s);const auto& md=d.modes(s);const auto& fi=s==0?d.sim.fCalc.fiL:d.sim.fCalc.fiR;
        std::array<double,3> modalWork{},nodalWork{};double scale=0,contactRes=0,contactScale=0;
        for(int m=0;m<md.nModes;++m) {
            const double q=st.q[m],v=st.qdot[m],dq=st.qf[m]-q,qp=d.qpos[s][m],qn=d.qneg[s][m],qf=d.qfluid[s][m],qc=d.qcontact[s][m];
            for(double a:{q,v,dq,qp,qn,qf,qc})requireFinite(a);
            const double wp=qp*dq,wn=qn*dq,wf=qf*dq,wc=qc*dq;
            modalWork[0]+=wn;modalWork[1]+=wp;modalWork[2]+=wf;scale+=std::abs(wn)+std::abs(wp);
            contactRes=std::max(contactRes,std::abs(fi[m]-qf-qc));contactScale=std::max(contactScale,std::abs(fi[m])+std::abs(qf)+std::abs(qc));
            const double damping=2*(s==0?d.sim.params.zetaL:d.sim.params.zetaR)*(s==0?d.sim.omegaL[m]:d.sim.omegaR[m])*v*v*d.sim.params.dt;
            for(int w:{0,window}){auto& r=d.sums[s][w][m];r.pos+=wp;r.neg+=wn;r.fluid+=wf;r.damping+=damping;r.contact+=wc;r.q2+=q*q;r.v2+=v*v;r.maxq=std::max(r.maxq,std::abs(q));++r.count;}
            const std::array<double,7> row{qp,qn,qf,qc,q,v,dq};d.raw.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
        }
        for(const auto& p:d.nodes[s]){const auto& a=st.disp[p.id];const auto& b=st.predictedDisp[p.id];
            double w=(p.force.x*(b.ux-a.ux)+p.force.y*(b.uy-a.uy)+p.force.z*(b.uz-a.uz))*1e-3;
            requireFinite(w);nodalWork[sign(p.p)<0?0:1]+=w;nodalWork[2]+=w;}
        double ws=std::abs(modalWork[0]+modalWork[1]-modalWork[2]);
        double nw=std::max({std::abs(nodalWork[0]-modalWork[0]),std::abs(nodalWork[1]-modalWork[1]),std::abs(nodalWork[2]-modalWork[2])});
        d.check(ws,scale,"work split");d.check(nw,scale,"actual nodal versus modal work");d.check(contactRes,contactScale,"adopted fluid+contact Q");
        for(int w:{0,window})for(int k=0;k<3;++k)d.cumulativeNodal[s][w][k]+=nodalWork[k];
        double scaled=std::max({scaledResidual(d.forceResidual[s],d.forceScale[s]),scaledResidual(d.qResidual[s],d.qScale[s]),scaledResidual(ws,scale),scaledResidual(nw,scale),scaledResidual(contactRes,contactScale)});
        d.consistency<<d.n<<','<<d.time<<','<<(s==0?'L':'R')<<",1,"<<d.forceResidual[s]<<','<<d.forceScale[s]<<','<<d.qResidual[s]<<','<<d.qScale[s]<<','<<ws<<','<<nw<<','<<scale<<','<<contactRes<<','<<scaled<<",1\n";
        for(int ap=0;ap<2;++ap){const auto& r=d.contacts[s][ap];d.contactFile<<d.n<<','<<d.time<<','<<(s==0?'L':'R')<<','<<(ap==0?"end":"interior")<<','<<d.iteration<<','<<r.pairs<<','<<r.excluded<<','<<r.fx<<','<<r.fy<<','<<openingForce(d.gapSign,s==0,r.fy)<<','<<(r.pairs?r.x/r.pairs:kMissing)<<','<<(r.pairs?r.z/r.pairs:kMissing)<<','<<r.penetration<<",last_adopted_predicted_contact_iterate\n";}
    }
    if(d.n%1000==0){d.consistency.flush();d.flow.flush();std::cout<<"[Kanaya] validated interval "<<d.n<<" time_s="<<d.time<<std::endl;}
}

void KanayaDiagnostics::finish() {
    if(!data_)return;auto& d=*data_;if(d.calls!=d.sim.params.nstep||d.integrated!=d.sim.params.nstep-1||d.excluded!=1)throw std::runtime_error("Kanaya interval accounting failure");
    std::array<std::set<int>,2> selected;
    for(int s=0;s<2;++s){const auto& md=d.modes(s);std::vector<int> rank(md.nModes);for(int m=0;m<md.nModes;++m)rank[m]=m;
        std::sort(rank.begin(),rank.end(),[&](int a,int b){return std::abs(d.sums[s][0][a].fluid)>std::abs(d.sums[s][0][b].fluid);});
        for(int m:{0,3,4})if(m<md.nModes)selected[s].insert(m);for(int k=0;k<std::min(5,md.nModes);++k)selected[s].insert(rank[k]);
        for(int w=0;w<4;++w)for(int m=0;m<md.nModes;++m){const auto& r=d.sums[s][w][m];int source=md.sourceModeIndices.empty()?m:md.sourceModeIndices[m];
            double start=w==0||w==1?0:w==2?.05:.1,end=w==1?.05:w==2?.1:d.integrated*d.sim.params.dt;
            d.check(r.pos+r.neg-r.fluid,std::abs(r.pos)+std::abs(r.neg),"cumulative work split");
            d.summary<<(s==0?'L':'R')<<','<<m+1<<','<<source+1<<','<<md.frequencies[m]<<','<<(w==0?"all":w==1?"ramp":w==2?"transient":"evaluation")<<','
                     <<r.pos<<','<<r.neg<<','<<r.fluid<<','<<r.damping<<','<<r.contact<<','<<(r.count?std::sqrt(r.q2/r.count):kMissing)<<','<<(r.count?std::sqrt(r.v2/r.count):kMissing)<<','<<r.maxq<<','<<r.count<<','<<start<<','<<std::min(end,d.integrated*d.sim.params.dt)<<'\n';}
    }
    for(int s=0;s<2;++s)for(int w=0;w<4;++w){std::array<double,3> modal{};double scale=0;
        for(const auto& r:d.sums[s][w]){modal[0]+=r.neg;modal[1]+=r.pos;modal[2]+=r.fluid;scale+=std::abs(r.neg)+std::abs(r.pos);}
        for(int k=0;k<3;++k)d.check(d.cumulativeNodal[s][w][k]-modal[k],scale,"cumulative actual nodal/modal work");}
    auto meta=output(d.dir/"run_metadata.json");
    meta<<"{\n\"complete\": true,\n\"diagnostic_enabled\": true,\n\"nstep\": "<<d.calls<<",\n\"dt_s\": "<<d.sim.params.dt
        <<",\n\"integrated_intervals\": "<<d.integrated<<",\n\"excluded_intervals\": "<<d.excluded
        <<",\n\"excluded_step\": "<<d.sim.params.nstep-1<<",\n\"work_interval_start_s\": 0,\n\"work_interval_end_s\": "<<d.integrated*d.sim.params.dt
        <<",\n\"last_force_sample_time_s\": "<<(d.integrated-1)*d.sim.params.dt
        <<",\n\"exclusion_reason\": \"User-approved existing State endpoint guard bug: exclude last work interval and all final contact statistics; solver still executes all steps\",\n"
        <<"\"actual_ramp_time_s\": "<<(d.sim.params.pressureRampTimeSecExplicit?d.sim.params.pressureRampTimeSec:.05)
        <<",\n\"force_csv_stride\": 20,\n\"modal_binary_fields\": [\"Q_pos\",\"Q_neg\",\"Q_fluid\",\"Q_contact_adopted\",\"q_n\",\"qdot_n\",\"delta_q\"],\n"
        <<"\"modal_binary_order\": \"step, side L/R, active mode, 7 native little-endian float64\",\n"
        <<"\"ap_end_rows\": 3,\n\"x_bin_width_mm\": 0.5,\n\"relative_tolerance\": 1e-8,\n\"absolute_floor\": 1e-12,\n"
        <<"\"work_units\": \"algebraic SI force-displacement work (J); per-mode energy interpretation conditional on sampled mass Gram matrix\",\n"
        <<"\"damping_definition\": \"sum(2*zeta*omega*qdot_n^2*dt); normalized modal units, not proven independent physical dissipation\",\n"
        <<"\"contact_work_definition\": \"last adopted predicted contact iterate Q times q_(n+1)-q_n; reference work, not dissipation\",\n"
        <<"\"gapSign\": "<<d.gapSign<<",\n\"nmode\": ["<<d.sim.mdataL.nModes<<','<<d.sim.mdataR.nModes<<"],\n\"gram\": [";
    bool first=true;for(int s=0;s<2;++s){const auto& md=d.modes(s);double mass=d.sim.params.mass/d.geom(s).nPoints;
        for(int a:selected[s])for(int b:selected[s]){double value=0;for(int k=0;k<md.nPoints;++k){const auto& u=md.modes[a][k];const auto& v=md.modes[b][k];value+=u.ux*v.ux+u.uy*v.uy+u.uz*v.uz;}value*=mass;
            if(!first)meta<<',';first=false;meta<<"[\""<<(s==0?'L':'R')<<"\","<<a+1<<','<<b+1<<','<<value<<']';}}
    meta<<"],\n\"force_statistics\": [";first=true;
    for(int s=0;s<2;++s)for(int w=0;w<4;++w)for(int sg=0;sg<2;++sg){const auto& r=d.forceStats[s][w][sg];if(!r.count)continue;
        if(!first)meta<<',';first=false;meta<<"{\"side\":\""<<(s==0?'L':'R')<<"\",\"window\":"<<w<<",\"pressure_sign\":"<<(sg==0?-1:1)
            <<",\"count\":"<<r.count<<",\"mean_opening_N\":"<<r.sum/r.count<<",\"mean_area_m2\":"<<r.area/r.count
            <<",\"min_opening_N\":"<<r.min<<",\"max_opening_N\":"<<r.max<<",\"min_time_s\":"<<r.minTime<<",\"max_time_s\":"<<r.maxTime<<'}';}
    meta<<"],\n\"cumulative_work_checks\": [";first=true;
    for(int s=0;s<2;++s)for(int w=0;w<4;++w){std::array<double,3> modal{};
        for(const auto& r:d.sums[s][w]){modal[0]+=r.neg;modal[1]+=r.pos;modal[2]+=r.fluid;}
        for(int k=0;k<3;++k){if(!first)meta<<',';first=false;meta<<"[\""<<(s==0?'L':'R')<<"\","<<w<<','<<k<<','<<modal[k]<<','<<d.cumulativeNodal[s][w][k]<<','<<d.cumulativeNodal[s][w][k]-modal[k]<<']';}}
    meta<<"]\n}\n";
    for(auto* f:{&d.regions,&d.summary,&d.consistency,&d.contactFile,&d.flow,&d.raw}){f->flush();if(!*f)throw std::runtime_error("Kanaya write failure");}
}
