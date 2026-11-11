//-------------------------------------------------------------------
// $Id: m_gem2mtsep.cpp 968 2007-12-13 13:23:32Z gems $
//
// Implementation of TGEM2MT class, standalone functions
//
// Copyright (C) 2022-2026  S. Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization
//
// This file may be distributed under the GPL v.3 license
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------
//

#include <array>
#include <algorithm>
#include <cstring>
#include <limits>
#include <sstream>
#include <map>
#include <iomanip>
#include <cmath>
#include "m_gem2mt.h"
#include "GEMS3K/nodearray.h"
#include "GEMS3K/io_keyvalue.h"
#include "GEMS3K/io_simdjson.h"
#ifdef USE_NLOHMANNJSON
#include "GEMS3K/io_nlohmann.h"
#endif

TGEM2MT* TGEM2MT::pm;

TGEM2MT::TGEM2MT(char ps_mode, long n_nodes)
{
    mtp=&mt[0];
    set_def(0);
    na = 0;
    pa_mt = 0;
    mtp->PsMode = ps_mode;
    mtp->nC = n_nodes;

    math_transport_defaults();
}

TGEM2MT::~TGEM2MT()
{
    mem_kill(0);
    if( pa_mt )
        delete pa_mt;
}

//Calculate record
void TGEM2MT::RecCalc()
{
    try {

        // use particles
        if(mtp->PsMode == RMT_MODE_W) {
            na->SetGrid(mtp->sizeLc, mtp->grid);   // set up grid structure
            pa_mt = new TParticleArray(mtp->nPTypes, mtp->nProps,
                                       mtp->NPmean, mtp->ParTD, mtp->nPmin, mtp->nPmax, na.get());
            pa_mt->setUpCounters();
        }

        // put HydP
        if(mtp->PsMode != RMT_MODE_S && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B) {
            putHydP(na->pNodT0());
            putHydP(na->pNodT1());
        }

        if(mtp->PsVTK != S_OFF)  {
            u_create_directory(pathVTK+nameVTK+"/");
        }

        if(mtp->iStat != AS_RUN) {  // ???? Could we stop and restart calculations?
            mtp->gStat = GS_GOING;
            mt_reset();
            mtp->gStat = GS_DONE;
        }

        gem_failures = GemFailures();
        gem_failures.ran = true;

        // internal calc
        auto iret = internalCalc();
        if(!iret) {
            mtp->iStat = AS_DONE;
        }
        //else we have a stop point
    }
    catch(TError& xcpt)  {
        mtp->gStat = GS_ERR;
        mtp->iStat = AS_INDEF;
        throw xcpt;
    }
}

// read TGEM2MT structure from file
int TGEM2MT::ReadTask(const std::string& gem2mt_file, const std::string& vtk_dir)
{
    try {
        std::string gem2mt_in = gem2mt_file;
        std::fstream ff(gem2mt_in, std::ios::in );
        ErrorIf(!ff.good(), gem2mt_in, "Fileopen error");

        // clear old allocation
        mem_kill(0);

        if(gem2mt_in.rfind(".json") != std::string::npos) {
#ifdef USE_NLOHMANNJSON
            io_formats::NlohmannJsonRead in_format(ff, "", "gem2mt");
            from_text_file( in_format );
        }
#else
            io_formats::SimdJsonRead in_format(ff, "", "gem2mt");
            from_text_file( in_format );
        }
#endif
        else {
            io_formats::KeyValueRead in_format( ff );
            from_text_file( in_format );
        }

        read_task_from_file = true;
        alloc_na();

        pathVTK = vtk_dir;
        if(!pathVTK.empty()) {
            pathVTK += "/";
        }
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
    }
    return 1;
}

int TGEM2MT::ReadTaskString(const std::string json_string)
{
    if(json_string.empty()) {
        return 1;
    }

    try  {

        // clear old allocation
        mem_kill(0);

        std::stringstream ss;
        ss.str(json_string);
#ifdef USE_NLOHMANNJSON
        io_formats::NlohmannJsonRead in_format(ss, "", "gem2mt");
        from_text_file( in_format );
#else
        io_formats::SimdJsonRead in_format(ss, "", "gem2mt");
        from_text_file(in_format);
#endif
        read_task_from_file = true;
        alloc_na();
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() << std::endl;
    }
    return 1;
}

// Write TGEM2MT structure to file
int TGEM2MT::WriteTask(const std::string& gem2mt_file)
{
    try  {
        std::string gem2mt_out = gem2mt_file;
        std::fstream ff(gem2mt_out, std::ios::out);
        ErrorIf(!ff.good(), gem2mt_out, "Fileopen error");

        if(gem2mt_out.rfind(".json") != std::string::npos) {
#ifdef USE_NLOHMANNJSON
            io_formats::NlohmannJsonWrite out_format( ff, "");
            to_text_file(out_format, true, false );
        }
#else
            io_formats::SimdJsonWrite out_format(ff, "", true);
            to_text_file(out_format, true, false);
        }
#endif
        else  {
            io_formats::KeyValueWrite out_format(ff);
            to_text_file(out_format, true, false);
        }
        return 0;
    }
    catch(TError& err) {
        // ???? add logger
        std::fstream f_log("gem2mtlog.txt", std::ios::out|std::ios::app );
        f_log << err.title.c_str() << "  : " << err.mess.c_str() <<std:: endl;
    }
    return 1;
}

void TGEM2MT::default_VTK(const std::string& work_path)
{
    std::string folder, fname, ext;
    if(!work_path.empty()) {
        u_splitpath(work_path, folder, fname, ext);
    }
    auto pos = fname.rfind("-");
    if(pos != std::string::npos) {
        fname = fname.substr(0, pos);
    }
    if(pathVTK.empty()) {
        pathVTK = folder+"/VTK/";
    }
    if(fname.empty()) {
        fname = "vtk";
    }
    if(nameVTK.empty()) {
        nameVTK = fname;
    }
    prefixVTK = nameVTK;
}


// Set up math transport default values7sizes in constructor
void TGEM2MT::math_transport_defaults()
{
    // from set_def(int q), ask gem2mt users for better defaults
    memset( &mtp->Msysb1, 0, sizeof(double)*20 );
    memset( mtp->size[0], 0, sizeof(float)*8 );

    read_task_from_file = false;
    mtp->nVTKfld = 0;
    mtp->sizeLc[0] = 1.;
    mtp->Tau[START_] = 0.;
    mtp->Tau[STOP_] = 1000.;
    mtp->Tau[STEP_] = 1.;
    mtp->ntM =1000;
    mtp->cdv = 1e-9;
    mtp->cez = 1e-12;
    mtp->tf = 1.;
    mtp->eps_in = 1.;

    mtp->Msysb1 = 0.;
    mtp->Vsysb1 = 0.;
    mtp->Mwatb1 = 1.;
    mtp->Maqb1 = 1.;
    mtp->Vaqb1 = 1.;

    // Alloc important arrays
    mtp->DiCp = new long int[ mtp->nC][2];
    defaults_DiCp();

    alloc_na();

    if(mtp->PsMode == RMT_MODE_W) {
        mtp->nPTypes = 10;
    }

    mtp->PsMO =   S_ON;
    mtp->iStat =  AS_READY;

    if(mtp->PsMode == RMT_MODE_B || mtp->PsMode == RMT_MODE_S || mtp->PsMode == RMT_MODE_F) {
        setNumberPhaseGroups(1);
    }

}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
int TGEM2MT::gem3k_files_read(const std::string& ipm_lst_file, const std::string& dbr_lst_file)
{
    // define name of vtk file
    default_VTK(ipm_lst_file);

    // Prepare the array for initial conditions allocation.
    // With only one initial system (dbr file) all nodes start from it; DiCp[q][0] stays as a group label
    std::vector<long int> nodeType;
    for(int ii=0; ii<mtp->nC; ++ii) {
        nodeType.push_back(mtp->nIV == 1 ? 0 : mtp->DiCp[ii][0]);
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    // if mtp->iStat == AS_RUN we resume calculation ????
    if(na->GEM_init(ipm_lst_file.c_str(), dbr_lst_file.c_str(), nodeType.data(), mtp->iStat == AS_RUN)) {
        // ???? add logger
        return 1;  // error reading files
    }
    return 0;
}

// Set up NodeArray and ParticleArray classes after reading gems3k files
int TGEM2MT::restore_data_from_gems3k(const std::vector<std::string>& dbr_names)
{
    int ii;
    // Node identification: the node index (GEMS3K sets the index of the initial system), so that
    // GEMS3K messages ("Error Node:...") name the node
    for(long int q=0; q<mtp->nC; ++q) {
        if(na->pNodT0()[q]) na->pNodT0()[q]->NodeHandle = q;
        if(na->pNodT1()[q]) na->pNodT1()[q]->NodeHandle = q;
    }
    CalcIPM(NEED_GEM_AIA, 0, mtp->nC); //recalc all nodes ?

    // realloc gem2mt memory  (if not read exported gem2mt)
    if(!read_task_from_file) {

        // Restore sizes from gems3k export
        mtp->Lsf = na->pCSD()->nDCs;
        mtp->Nf = na->pCSD()->nIC;
        mtp->FIf = na->pCSD()->nPH;
        mtp->nTai = na->pCSD()->nTp;
        mtp->nPai = na->pCSD()->nPp;

        // allocate memory and setup default values
        mem_new(0);
        init_arrays(true);
        // initialization from scripts in GUI
        if(mtp->HydP) {
            defaults_HydP();
        }
        if(mtp->PvGrid != S_OFF && mtp->grid) {
            defaults_Grid();
        }

        // read names
        std::string name;
        for(ii=0; ii<mtp->nIV; ++ii) {
            if(ii<dbr_names.size()) {
                name = std::to_string(ii)+dbr_names[ii];
                strncpy(mtp->nam_i[ii], name.c_str(), MAXIDNAME);
            }
        }
    }

    for(ii=0; ii<mtp->nTai; ++ii) {
        mtp->Tval[ii] = na->pCSD()->TKval[ii]-C_to_K;
    }
    for(ii=0; ii<mtp->nPai; ++ii) {
        mtp->Pval[ii] = na->pCSD()->Pval[ii]/bar_to_Pa;
    }

    return 0;
}

// Here we read the MULTI structure, DATACH and DATABR strings prepared from GEMS
int TGEM2MT::gems3k_strings(const std::string& dch_json, const std::string& ipm_json,
                            const std::vector<std::string>& dbr_json)
{
    // define name of vtk file
    default_VTK("");

    // Prepare the array for initial conditions allocation (one dbr: all nodes start from it)
    std::vector<long int> nodeType;
    for(int ii=0; ii<mtp->nC; ++ii) {
        nodeType.push_back(dbr_json.size() == 1 ? 0 : mtp->DiCp[ii][0]);
    }

    // Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
    if(na->GEM_init(dch_json, ipm_json, dbr_json, nodeType.data())) {
        // ???? add logger
        return 1;  // error reading files
    }
    return 0;
}


int TGEM2MT::MassTransInit(const std::string &ipm_lst_file, const std::string &dbr_lst_file)
{    gem_failures = GemFailures();     // the calculations of the initialization are not part of a run
    // get dbr file names first: their number decides how the nodes are initialized
    GEMS3KGenerator generator(ipm_lst_file);
    mtp->nIV = generator.load_dbr_lst_file(dbr_lst_file);
    dbr_file_names = generator.dbr_names();

    if(gem3k_files_read(ipm_lst_file, dbr_lst_file)) {
        return 1;
    }

    restore_data_from_gems3k(generator.dbr_names());
    return 0;
}

// Here we read the MULTI structure, DATACH and DATABR files prepared from GEMS
// Set up NodeArray and ParticleArray classes
int TGEM2MT::MassTransStringInit(const std::string& dch_json, const std::string& ipm_json,
                                 const std::vector<std::string>& dbr_json)
{    gem_failures = GemFailures();
    if(gems3k_strings(dch_json, ipm_json, dbr_json)) {
        return 1;
    }

    mtp->nIV = dbr_json.size();
    dbr_file_names.clear();
    for(size_t i=0; i<dbr_json.size(); ++i)
        dbr_file_names.push_back("dbr string "+std::to_string(i));
    restore_data_from_gems3k({});
    return 0;
}

void TGEM2MT::setVTKfields(const std::vector<std::pair<int, int>> &vtk_fields)
{
    if(mtp->xVTKfld) {
        delete[]  mtp->xVTKfld;
        mtp->xVTKfld = nullptr;
    }
    mtp->nVTKfld = vtk_fields.size();
    if(mtp->nVTKfld>0) {
        mtp->xVTKfld = new long int[mtp->nVTKfld][2];
        for(int ii=0; ii<mtp->nVTKfld; ++ii) {
            mtp->xVTKfld[ii][0] =vtk_fields[ii].first;
            mtp->xVTKfld[ii][1] =vtk_fields[ii].second;
        }
        mtp->PvnVTK = S_ON;
        mtp->PsVTK = mtp->PvnVTK;
    }
    else {
        mtp->PvnVTK = S_OFF;
        mtp->PsVTK = mtp->PvnVTK;
    }
}

void TGEM2MT::setParticle(size_t pndx, long int pmean, long int pmin, long int pmax, const std::array<long int, 6> &pparam)
{
    if(mtp->PsMode == RMT_MODE_W && pndx<mtp->nPTypes) {
        mtp->NPmean[pndx] = pmean;
        mtp->nPmin[pndx] = pmin;
        mtp->nPmax[pndx] = pmax;
        std::copy(pparam.begin(), pparam.end(), mtp->ParTD[pndx]);
    }
}

std::array<double, SIZE_HYDP> TGEM2MT::hydraulicParameters(size_t pndx) const
{
    ErrorIf(!mtp->HydP, GetName(), "Hydraulic parameters are not allocated: call initNodes first");
    ErrorIf(pndx >= static_cast<size_t>(mtp->nC), GetName(),
            "Node index "+std::to_string(pndx)+" out of range (nodes: "+std::to_string(mtp->nC)+")");
    std::array<double, SIZE_HYDP> hp;
    for(size_t i=0; i<SIZE_HYDP; ++i)
        hp[i] = mtp->HydP[pndx][i];
    return hp;
}

double TGEM2MT::materialVolume(size_t q, const std::vector<double>& b, double T, double P, double limit_scale) const
{
    ErrorIf(!na || !na->pNodT1() || !na->pNodT1()[0], GetName(), "Nodes are not initialized: call initNodes first");
    ErrorIf(q >= static_cast<size_t>(mtp->nC), GetName(),
            "Node index "+std::to_string(q)+" out of range (nodes: "+std::to_string(mtp->nC)+")");
    const DATACH* ch = na->pCSD();
    ErrorIf(static_cast<long int>(b.size()) != ch->nIC, GetName(), "Wrong number of element amounts");
    std::vector<bool> kept(ch->nIC, false);
    std::vector<double> b_dbr(ch->nICb);
    for(long int j=0; j<ch->nICb; ++j) {
        // elements that are not in the material get the minimal amount (as in the node calculations)
        b_dbr[j] = std::max(b[ch->xic[j]], mtp->cez);
        kept[ch->xic[j]] = true;
    }
    for(long int i=0; i<ch->nIC; ++i)
        ErrorIf(!kept[i] && b[i] != 0., GetName(),
                "Element "+ch->ICNL[i]+" is not kept in the node data (DBR); its amount must be 0");

    // A copy of the work node, with the bounds of the species (scaled to the amounts of this calculation)
    // and the surface areas of node q
    const DATABR* nq = na->pNodT1()[q];
    std::vector<double> dul(nq->dul, nq->dul + ch->nDCb), dll(nq->dll, nq->dll + ch->nDCb);
    for(auto& x: dul) x *= limit_scale;
    for(auto& x: dll) x *= limit_scale;
    TNode work(*na->getCalcNode());
    work.GEM_from_MT(0, NEED_GEM_AIA, T, P, 0., 0., b_dbr.data(), dul.data(), dll.data(), nq->aPH);
    const long int st = work.GEM_run(false);
    ErrorIf(!(st == OK_GEM_AIA || st == OK_GEM_SIA), GetName(),
            "GEM calculation failed for the volume of a part of the node (status "+std::to_string(st)+"): "+
            work.ipmLogError());
    return work.cVs() * 1e3;     // m3 -> L
}

TGEM2MT::NodeInputs TGEM2MT::saveNodeInputs(size_t q) const
{
    ErrorIf(!na || !na->pNodT1() || !na->pNodT1()[0], GetName(), "Nodes are not initialized: call initNodes first");
    ErrorIf(q >= static_cast<size_t>(mtp->nC), GetName(),
            "Node index "+std::to_string(q)+" out of range (nodes: "+std::to_string(mtp->nC)+")");
    const DATABR* n = na->pNodT1()[q];
    NodeInputs in;
    in.bIC.assign(n->bIC, n->bIC + na->pCSD()->nICb);
    in.TK = n->TK;
    in.P = n->P;
    return in;
}

void TGEM2MT::restoreNodeInputs(size_t q, const NodeInputs& in)
{
    ErrorIf(!na || !na->pNodT0() || !na->pNodT1() || !na->pNodT0()[0], GetName(),
            "Nodes are not initialized: call initNodes first");
    ErrorIf(q >= static_cast<size_t>(mtp->nC) || static_cast<long int>(in.bIC.size()) != na->pCSD()->nICb, GetName(),
            "Cannot restore node "+std::to_string(q));
    DATABRPTR* C0 = na->pNodT0();
    DATABRPTR* C1 = na->pNodT1();
    std::copy(in.bIC.begin(), in.bIC.end(), C0[q]->bIC);
    std::copy(in.bIC.begin(), in.bIC.end(), C1[q]->bIC);
    C0[q]->TK = C1[q]->TK = in.TK;
    C0[q]->P = C1[q]->P = in.P;
    TestModeGEMParam mode(NEED_GEM_AIA, S_OFF, mtp->ct, mtp->cdv, mtp->cez);
    na->CalcIPM_One(mode, static_cast<long int>(q), nullptr);                 // the earlier inputs again
    na->CopyNodeFromTo(static_cast<long int>(q), mtp->nC, C1, C0);
}

void TGEM2MT::setNodeBulkComposition(size_t q, const std::vector<double>* b, double T, double P)
{
    ErrorIf(!na || !na->pNodT0() || !na->pNodT1() || !na->pNodT0()[0], GetName(),
            "Nodes are not initialized: call initNodes first");
    ErrorIf(q >= static_cast<size_t>(mtp->nC), GetName(),
            "Node index "+std::to_string(q)+" out of range (nodes: "+std::to_string(mtp->nC)+")");
    const DATACH* ch = na->pCSD();
    const long int nICb = ch->nICb;

    // New element amounts in DBR order: the DBR keeps nICb of the nIC elements (index list xic)
    std::vector<double> b_dbr;
    if(b) {
        ErrorIf(static_cast<long int>(b->size()) != ch->nIC, GetName(),
                "Bulk composition has "+std::to_string(b->size())+" values, the chemical system has "+
                std::to_string(ch->nIC)+" elements");
        std::vector<bool> kept(ch->nIC, false);
        b_dbr.resize(nICb);
        for(long int j=0; j<nICb; ++j) {
            b_dbr[j] = (*b)[ch->xic[j]];
            kept[ch->xic[j]] = true;
        }
        for(long int i=0; i<ch->nIC; ++i)
            ErrorIf(!kept[i] && (*b)[i] != 0., GetName(),
                    "Element "+ch->ICNL[i]+" is not kept in the node data (DBR); its amount must be 0");
    }

    DATABRPTR* C0 = na->pNodT0();
    DATABRPTR* C1 = na->pNodT1();
    // The previous inputs, to restore them if the calculation fails
    const NodeInputs previous = saveNodeInputs(q);

    if(b)
        for(long int j=0; j<nICb; ++j)
            C0[q]->bIC[j] = C1[q]->bIC[j] = b_dbr[j];
    if(std::isfinite(T))
        C0[q]->TK = C1[q]->TK = T;
    if(std::isfinite(P))
        C0[q]->P = C1[q]->P = P;

    // Equilibrate the node from scratch into C1, then make C0 the same
    TestModeGEMParam mode(NEED_GEM_AIA, S_OFF, mtp->ct, mtp->cdv, mtp->cez);
    bool ok = false;
    try {
        ok = na->CalcIPM_One(mode, static_cast<long int>(q), nullptr);
    }
    catch(...) {
        restoreNodeInputs(q, previous);
        throw;
    }
    if(!ok) {
        restoreNodeInputs(q, previous);
        Error(GetName(), "GEM calculation failed for the new composition of node "+std::to_string(q)+
                         "; the node keeps its previous composition");
    }
    na->CopyNodeFromTo(static_cast<long int>(q), mtp->nC, C1, C0);
}

void TGEM2MT::setHydraulicParameters(size_t pndx, double Vt, double vp, double eps, double Km, double al, double Dif, double nto)
{
    if(mtp->HydP && pndx<mtp->nC) {
        mtp->HydP[pndx][0] = Vt;
        mtp->HydP[pndx][1] = vp;
        mtp->HydP[pndx][2] = eps;
        mtp->HydP[pndx][3] = Km;
        mtp->HydP[pndx][4] = al;
        mtp->HydP[pndx][5] = Dif;
        mtp->HydP[pndx][6] = nto;
    }
}

void TGEM2MT::setPhaseGroupsID(size_t gndx, const std::string &ids)
{
    if(mtp->MGPid && gndx<mtp->nPG) {
        strncpy( mtp->MGPid[gndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setUnitsPhaseQuantities(size_t pndx, char units)
{
    if(mtp->UMGP && pndx<mtp->FIf) {
        mtp->UMGP[pndx] = units;
    }
}

void TGEM2MT::setPhaseGroupsQuantities(size_t gndx, size_t pndx, double quantity)
{
    if(mtp->MGPid && gndx<mtp->nPG && pndx<mtp->FIf) {
        mtp->PGT[gndx*mtp->FIf+pndx] = quantity;
    }
}

void TGEM2MT::setICsourceQuantities(size_t gndx, size_t indx, double quantity)
{
    if(mtp->MGPid && gndx<mtp->nSFD && indx<mtp->Nf) {
        mtp->BSF[gndx*mtp->Nf+indx] = quantity;
    }
}

void TGEM2MT::setFluxDefinition(size_t pndx, double order, double rate, double quantity, double val)
{
    if(mtp->FDLf && pndx<mtp->nFD) {
        mtp->FDLf[pndx][0] = order;
        mtp->FDLf[pndx][1] = rate;
        mtp->FDLf[pndx][2] = quantity;
        mtp->FDLf[pndx][3] = val;
    }
}

void TGEM2MT::setFluxMGPid(size_t pndx, const std::string &ids)
{
    if(mtp->FDLmp && pndx<mtp->nFD) {
        strncpy( mtp->FDLmp[pndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setFluxID(size_t pndx, const std::string &ids)
{
    if(mtp->FDLid && pndx<mtp->nFD) {
        strncpy( mtp->FDLid[pndx], ids.c_str(), MAXSYMB);
    }
}

void TGEM2MT::setGridPoint(size_t pndx, double x, double y, double z)
{
    if(mtp->grid && pndx<mtp->nC) {
        mtp->grid[pndx][0] = x;
        mtp->grid[pndx][1] = y;
        mtp->grid[pndx][2] = z;
    }
}

void TGEM2MT::setFluxSourceReceive(size_t pndx, long int source, long int receive)
{
    if(mtp->FDLi && pndx<mtp->nFD) {
        mtp->FDLi[pndx][0] = source;
        mtp->FDLi[pndx][1] = receive;
    }
}

//==========================================================================================

// free dynamic memory in objects and values
void TGEM2MT::mem_kill(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
            "E05GTrem: Attempt to access corrupted dynamic memory.");

    //- if( mtp->lNam) delete[] mtp->lNam;
    //- if( mtp->lNamE) delete[] mtp->lNamE;
    //- if( mtp->tExpr) delete[] mtp->tExpr;
    //- if( mtp->gExpr) delete[] mtp->gExpr;
    if(mtp->sdref) {
        delete[] mtp->sdref;
        mtp->sdref = nullptr;
    }
    if(mtp->sdval) {
        delete[] mtp->sdval;
        mtp->sdval = nullptr;
    }
    if(mtp->DiCp) {
        delete[] mtp->DiCp;
        mtp->DiCp = nullptr;
    }
    if(mtp->FDLi) {
        delete[] mtp->FDLi;
        mtp->FDLi = nullptr;
    }
    //- if( mtp->PTVm) delete[] mtp->PTVm;
    //- if( mtp->StaP) delete[] mtp->StaP;
    if(mtp->xVTKfld) {
        delete[] mtp->xVTKfld;
        mtp->xVTKfld = nullptr;
    }
    //- if( mtp->xEt) delete[] mtp->xEt;
    //- if( mtp->yEt) delete[] mtp->yEt;
    //- if( mtp->Bn) delete[] mtp->Bn;
    if(mtp->HydP) {
        delete[] mtp->HydP;
        mtp->HydP = nullptr;
    }
    //- if( mtp->qpi) delete[] mtp->qpi;
    //- if( mtp->qpc) delete[] mtp->qpc;
    //- if( mtp->xt) delete[] mtp->xt;
    //- if( mtp->yt) delete[] mtp->yt;
    //- if( mtp->CIb) delete[] mtp->CIb;
    //- if( mtp->CAb) delete[] mtp->CAb;
    if(mtp->FDLf) {
        delete[] mtp->FDLf;
        mtp->FDLf = nullptr;
    }
    if(mtp->PGT) {
        delete[] mtp->PGT;
        mtp->PGT = nullptr;
    }
    if(mtp->Tval) {
        delete[] mtp->Tval;
        mtp->Tval = nullptr;
    }
    if(mtp->Pval) {
        delete[] mtp->Pval;
        mtp->Pval = nullptr;
    }
    if(mtp->nam_i) {
        delete[] mtp->nam_i;
        mtp->nam_i = nullptr;
    }
    //- if( mtp->for_i) delete[] mtp->for_i;
    //- if( mtp->stld) delete[] mtp->stld;
    //- if( mtp->CIclb) delete[] mtp->CIclb;
    //- if( mtp->AUcln) delete[] mtp->AUcln;
    if(mtp->FDLid) {
        delete[] mtp->FDLid;
        mtp->FDLid = nullptr;
    }
    if(mtp->FDLop) {
        delete[] mtp->FDLop;
        mtp->FDLop = nullptr;
    }
    if(mtp->FDLmp) {
        delete[] mtp->FDLmp;
        mtp->FDLmp = nullptr;
    }
    if(mtp->MGPid) {
        delete[] mtp->MGPid;
        mtp->MGPid = nullptr;
    }
    if(mtp->UMGP) {
        delete[] mtp->UMGP;
        mtp->UMGP = nullptr;
    }
    //- if( mtp->SBM) delete[] mtp->SBM;
    if(mtp->BSF) {
        delete[] mtp->BSF;
        mtp->BSF = nullptr;
    }
    if(mtp->MB) {
        delete[] mtp->MB;
        mtp->MB = nullptr;
    }
    if(mtp->dMB) {
        delete[] mtp->dMB;
        mtp->dMB = nullptr;
    }
    if(mtp->DDc) {
        delete[] mtp->DDc;
        mtp->DDc = nullptr;
    }
    if(mtp->DIc) {
        delete[] mtp->DIc;
        mtp->DIc = nullptr;
    }
    if(mtp->DEl) {
        delete[] mtp->DEl;
        mtp->DEl = nullptr;
    }
    if(mtp->for_e) {
        delete[] mtp->for_e;
        mtp->for_e = nullptr;
    }
    //- if( mtp->xIC) delete[] mtp->xIC;
    //- if( mtp->xDC) delete[] mtp->xDC;
    //- if( mtp->xPH) delete[] mtp->xPH;
    if(mtp->grid) {
        delete[] mtp->grid;
        mtp->grid = nullptr;
    }
    if(mtp->NPmean) {
        delete[] mtp->NPmean;
        mtp->NPmean = nullptr;
    }
    if(mtp->nPmin) {
        delete[] mtp->nPmin;
        mtp->nPmin = nullptr;
    }
    if(mtp->nPmax) {
        delete[] mtp->nPmax;
        mtp->nPmax = nullptr;
    }
    if(mtp->ParTD) {
        delete[] mtp->ParTD;
        mtp->ParTD = nullptr;
    }
    if(mtp->BM) {
        delete[] mtp->BM;
        mtp->BM = nullptr;
    }
    if(mtp->BdM) {
        delete[] mtp->BdM;
        mtp->BdM = nullptr;
    }
    if(mtp->FmgpJ) {
        delete[] mtp->FmgpJ;
        mtp->FmgpJ = nullptr;
    }
    if(mtp->BmgpM) {
        delete[] mtp->BmgpM;
        mtp->BmgpM = nullptr;
    }
    // work
    //- if( mtp->An) delete[] mtp->An;
    //- if( mtp->Ae) delete[] mtp->Ae;
    if(mtp->gfc) {
        delete[] mtp->gfc;
        mtp->gfc = nullptr;
    }
    if(mtp->yfb) {
        delete[] mtp->yfb;
        mtp->yfb = nullptr;
    }
    if(mtp->tt) {
        delete[] mtp->tt;
        mtp->tt = nullptr;
    }
    //- if( mtp->etext) delete[] mtp->etext;
    //- if( mtp->tprn) delete[] mtp->tprn;
    //- FreeNa();
    //- freeNodeWork();
}

// realloc dynamic memory
void TGEM2MT::mem_new(int q)
{
    ErrorIf( mtp!=&mt[q], GetName(),
            "E04GTrem: Attempt to access corrupted dynamic memory.");

    //- mtp->xIC = new long int[mtp->nICb];
    //- mtp->xDC = new long int[mtp->nDCb];
    //- mtp->xPH = new long int[mtp->nPHb];

    if( mtp->PvGrid == S_OFF )
    { if(mtp->grid) delete[] mtp->grid;
        mtp->grid = 0;
    }
    else
        mtp->grid = new double[ mtp->nC][3];

    if(mtp->PsMode == RMT_MODE_W) {
        mtp->NPmean = new long int[ mtp->nPTypes];
        mtp->nPmin = new long int[ mtp->nPTypes];
        mtp->nPmax = new long int[ mtp->nPTypes];
        mtp->ParTD = new long int[mtp->nPTypes][6];
    }
    else {
        if(mtp->NPmean) delete[] mtp->NPmean;
        if(mtp->nPmin) delete[] mtp->nPmin;
        if(mtp->nPmax) delete[] mtp->nPmax;
        if(mtp->ParTD) delete[] mtp->ParTD;
        mtp->NPmean = nullptr;
        mtp->nPmin = nullptr;
        mtp->nPmax = nullptr;
        mtp->ParTD = nullptr;
    }

    mtp->nam_i= new char[mtp->nIV][ MAXIDNAME ];
    //- mtp->PTVm = new double[ mtp->nIV][5];
    if(!mtp->DiCp) { // could be allocated in constructor
        mtp->DiCp = new long int[ mtp->nC][2];
    }
    //- mtp->StaP = new double[ mtp->nC ][4];

    if(mtp->PvnVTK == S_OFF) {
        if(mtp->xVTKfld) {
            delete[] mtp->xVTKfld;
            mtp->xVTKfld = nullptr;
        }
    }
    else {
        if(!mtp->xVTKfld) { // could be allocated in function setVTKfields
            mtp->xVTKfld = new long int[mtp->nVTKfld][2];
        }
    }

    //- mtp->stld = new char[ mtp->nIV ][EQ_RKLEN];
    mtp->Tval  = new double[ mtp->nTai ];
    mtp->Pval  = new double[ mtp->nPai ];
    //- mtp->Bn = new double[ mtp->nIV][ mtp->Nb ];
    //- mtp->SBM = new char [ mtp->Nb][MAXICNAME+MAXSYMB];

    if(mtp->PsMode != RMT_MODE_S  && mtp->PsMode != RMT_MODE_F && mtp->PsMode != RMT_MODE_B) {
        mtp->HydP = new double[mtp->nC][SIZE_HYDP];
    }
    else {
        if(mtp->HydP) delete[] mtp->HydP;
        mtp->HydP = nullptr;
    }

    //-if( mtp->PvICi == S_OFF )
    //-   {
    //-    if(mtp->CIb) delete[] mtp->CIb;
    //-    if(mtp->CIclb) delete[] mtp->CIclb;
    //-    mtp->CIb = 0;
    //-    mtp->CIclb = 0;
    //-   }
    //-   else
    //-   {
    //-    mtp->CIb = new double[ mtp->nIV][mtp->Nb];
    //-    mtp->CIclb =  new char[ mtp->Nb ];
    //-   }

    //- if( mtp->PvAUi == S_OFF )
    //-    {
    //-     if(mtp->CAb) delete[] mtp->CAb;
    //-      if(mtp->for_i) delete[] mtp->for_i;
    //-      if(mtp->AUcln) delete[] mtp->AUcln;
    //-      if(mtp->An) delete[] mtp->An;
    //-      mtp->CAb = 0;
    //-      mtp->for_i = 0;
    //-      mtp->AUcln = 0;
    //-      mtp->An = 0;
    //-      mtp->Lbi = 0;
    //-    }
    //-    else
    //-    {
    //-      mtp->CAb = new double[ mtp->nIV][ mtp->Lbi ];
    //-      mtp->for_i = new char[ mtp->Lbi][ MAXFORMUNITDT ];
    //-      mtp->AUcln = new char[ mtp->Lbi ];
    //-      mtp->An = new double[ mtp->Lbi][ mtp->Nb ];
    //-   }

    if(mtp->PvFDL == S_OFF) {
        if(mtp->FDLi) delete[] mtp->FDLi;
        if(mtp->FDLf) delete[] mtp->FDLf;
        if(mtp->FDLid) delete[] mtp->FDLid;
        if(mtp->FDLop) delete[] mtp->FDLop;
        if(mtp->FDLmp) delete[] mtp->FDLmp;
        mtp->FDLi = nullptr;
        mtp->FDLf = nullptr;
        mtp->FDLid = nullptr;
        mtp->FDLop = nullptr;
        mtp->FDLmp = nullptr;
        mtp->nFD = 0;
    }
    else  {
        mtp->FDLi = new long int[ mtp->nFD][2];
        mtp->FDLf = new double[ mtp->nFD][4];
        mtp->FDLid= new char[ mtp->nFD][MAXSYMB];
        mtp->FDLop= new char[ mtp->nFD][MAXSYMB];
        mtp->FDLmp = new char[ mtp->nFD][MAXSYMB];
        for(long int ii=0; ii<mtp->nFD; ++ii) {
            fillValue(mtp->FDLid[ii], '\0', MAXSYMB);
            fillValue(mtp->FDLop[ii], '\0', MAXSYMB);
            fillValue(mtp->FDLmp[ii], '\0', MAXSYMB);
        }
    }

    if(mtp->PvPGD == S_OFF) {
        if(mtp->PGT) delete[] mtp->PGT;
        if(mtp->MGPid) delete[] mtp->MGPid;
        if(mtp->UMGP) delete[] mtp->UMGP;
        mtp->PGT = nullptr;
        mtp->MGPid = nullptr;
        mtp->UMGP = nullptr;
        mtp->nPG = 0;
    }
    else  {
        mtp->PGT  =  new double[ mtp->FIf*mtp->nPG ];
        mtp->MGPid = new char[ mtp->nPG ][MAXSYMB];
        mtp->UMGP = new char[ mtp->FIf ];
        for(long int ii=0; ii<mtp->nPG; ++ii) {
            fillValue(mtp->MGPid[ii], '\0', MAXSYMB);
        }
        for(long int ii=0; ii<mtp->FIf; ++ii) {
            mtp->UMGP[ii] = ' ';
        }
    }

    if( mtp->PvSFL == S_OFF )
    { if(mtp->BSF) delete[] mtp->BSF;
        mtp->BSF = 0;
    }
    else
        mtp->BSF = new double[ mtp->nSFD*mtp->Nf ];

    if(mtp->PvPGD != S_OFF && mtp->PvFDL != S_OFF) {
        mtp->MB =  new double[mtp->nC*mtp->Nf];
        mtp->dMB = new double[mtp->nC*mtp->Nf];
    }
    else {
        if(mtp->MB) delete[] mtp->MB;
        if(mtp->dMB) delete[] mtp->dMB;
        mtp->MB = nullptr;
        mtp->dMB = nullptr;
    }
    if( mtp->PvDDc == S_OFF )
    {
        if(mtp->DDc) delete[] mtp->DDc;
        mtp->DDc = 0;
    }
    else
        mtp->DDc = new double[mtp->Lsf];

    if( mtp->PvDIc == S_OFF )
    {
        if(mtp->DIc) delete[] mtp->DIc;
        mtp->DIc = 0;
    }
    else
        mtp->DIc = new double[ mtp->Nf ];

    if( mtp->nEl <= 0  )
    {
        if(mtp->DEl) delete[] mtp->DEl;
        if(mtp->for_e) delete[] mtp->for_e;
        //-   if(mtp->Ae) delete[] mtp->Ae;
        mtp->DEl = 0;
        mtp->for_e = 0;
        //-   mtp->Ae = 0;
        mtp->nEl = 0;
    }
    else
    {
        mtp->DEl = new double[ mtp->nEl ];
        mtp->for_e = new char[mtp->nEl][MAXFORMUNITDT];
        //-     mtp->Ae = new double[ mtp->nEl*mtp->Nb ];
    }

    //----------------------------------------------------------------
    //- if( mtp->Nqpt > 0  )
    //-  mtp->qpi   = new double[mtp->Nqpt];
    //- else
    //- { if(mtp->qpi) delete[] mtp->qpi; mtp->qpi = 0;}

    //- if( mtp->Nqpg > 0  )
    //-  mtp->qpc   = new double[mtp->Nqpg];
    //- else
    //- { if(mtp->qpc) delete[] mtp->qpc; mtp->qpc = 0;}

    //- if( mtp->PvMSt == S_OFF )
    //- { if(mtp->tExpr) delete[] mtp->tExpr; mtp->tExpr = 0;}
    //- else
    //-    mtp->tExpr = new char[4096];

    //-if( mtp->PvMSg == S_OFF )
    //-   {
    //-    if(mtp->lNam) delete[] mtp->lNam;
    //-    if(mtp->gExpr) delete[] mtp->gExpr;
    //-    if(mtp->xt) delete[] mtp->xt;
    //-    if(mtp->yt) delete[] mtp->yt;
    //-    mtp->lNam = 0;
    //-    mtp->gExpr = 0;
    //-    mtp->xt = 0;
    //-    mtp->yt = 0;
    //-   }
    //-   else
    //-   {
    //-        mtp->lNam = new char[ mtp->nYS][ MAXGRNAME];
    //-        mtp->gExpr = new char[2048];
    //-        mtp->xt   = new double[ mtp->nC];
    //-        mtp->yt   = new double[ mtp->nC*mtp->nYS];
    //-   }

    //-if( mtp->PvEF == S_OFF )
    //-   {
    //-    if(mtp->lNamE) delete[] mtp->lNamE;
    //-    if(mtp->xEt) delete[] mtp->xEt;
    //-    if(mtp->yEt) delete[] mtp->yEt;
    //-    mtp->lNamE = 0;
    //-    mtp->xEt = 0;
    //-    mtp->yEt = 0;
    //-   }
    //-   else
    //-   {
    //-     mtp->lNamE = new char[ mtp->nYE][ MAXGRNAME];
    //-     mtp->xEt   = new double[ mtp->nE ];
    //-     mtp->yEt   = new double[ mtp->nE*mtp->nYE ];
    //-   }

    if( mtp->Nsd > 0 )
    {
        mtp->sdref = new char[ mtp->Nsd ][ V_SD_RKLEN ];
        mtp->sdval = new char[ mtp->Nsd ][ V_SD_VALEN ];
    }
    else
    {
        if(mtp->sdref) delete[] mtp->sdref;
        if(mtp->sdval) delete[] mtp->sdval;
        mtp->sdref = 0;
        mtp->sdval = 0;
    }
    //- mtp->etext = new char[4096];
    //- mtp->tprn = new char[2048];
    //mtp->gfc = (double *)aObj[ o_mtgfc].Free();
    //mtp->yfb = (double *)aObj[ o_mtyfb].Free();
    //mtp->tt = (double *)aObj[ o_mttt].Free();
}

// Conversion of concentration units to moles
double TGEM2MT::Reduce_Conc( char UNITP, double Xe, double DCmw, double Vm,
                            double R1, double Msys, double Mwat, double Vaq, double Maq, double Vsys )
{
    double Xincr = 0.;
    switch( UNITP )
    {  // Quantities
    case QUAN_MKMOL: /*'Y'*/
        Xincr = Xe / 1e6;
        goto FINISH;
    case QUAN_MMOL:  /*'h'*/
        Xincr = Xe / 1e3;
        goto FINISH;
    case QUAN_MOL:   /*'M'*/
        Xincr = Xe;
        goto FINISH;
    }
    if( DCmw > 1e-12 )
        switch( UNITP )
        {
        case QUAN_MGRAM: /*'y'*/
            Xincr = Xe / DCmw / 1e3;
            goto FINISH;
        case QUAN_GRAM:  /*'g'*/
            Xincr = Xe / DCmw;
            goto FINISH;
        case QUAN_KILO:  /*'G'*/
            Xincr = Xe * 1e3 / DCmw;
            goto FINISH;
        }
    /* Concentrations */
    if( fabs( R1 ) > 1e-12 )
        switch( UNITP )
        { // mole fractions relative to total moles in the system
        case CON_MOLFR:  /*'n'*/
            Xincr = Xe * R1;
            goto FINISH;
        case CON_MOLPROC:/*'N'*/
            Xincr = Xe / 100. * R1;
            goto FINISH;
        case CON_pMOLFR: /*'f'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow(10., -Xe )* R1;
            goto FINISH;
        }
    if( fabs( Vsys ) > 1e-12 && Vm > 1e-12 )
        switch( UNITP )   /* Volumes */
        {
        case CON_VOLFR:  /*'v'*/
            Xincr = Xe * Vsys * 1e3 / Vm;
            goto FINISH;
        case CON_VOLPROC:/*'V'*/
            Xincr = Xe * Vsys * 10. / Vm;
            goto FINISH;
        case CON_pVOLFR: /*'u'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vsys / Vm;
            goto FINISH;
        }
    if( fabs( Msys ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP ) // Mass fractions relative to mass of the system
        {
        case CON_WTFR:   /*'w'*/
            Xincr = Xe * Msys * 1e3 / DCmw;
            goto FINISH;
        case CON_WTPROC: /*'%'*/
            Xincr = Xe * Msys *10. / DCmw;
            goto FINISH;
        case CON_PPM:    /*'P'*/
            Xincr = Xe * Msys / 1e3 / DCmw;
            goto FINISH;
        }
    if( fabs( Mwat ) > 1e-12 )
        switch( UNITP ) /* Molalities */
        {
        case CON_MOLAL:  /*'m'*/
            Xincr = Xe * Mwat;
            goto FINISH;
        case CON_MMOLAL: /*'i'*/
            Xincr = Xe / 1e3 * Mwat;
            goto FINISH;
        case CON_pMOLAL: /*'p'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Mwat;
            goto FINISH;
        }
    if( fabs( Vaq ) > 1e-12 )
        switch( UNITP )  /* Molarities */
        {
        case CON_MOLAR:  /*'L'*/
            Xincr = Xe * Vaq;
            goto FINISH;
        case CON_MMOLAR: /*'j'*/
            Xincr = Xe * Vaq / 1e3;
            goto FINISH;
        case CON_pMOLAR: /*'q'*/
            if( Xe > -1. && Xe < 15 )
                Xincr = pow( 10., -Xe ) * Vaq;
            goto FINISH;
            /* g/l, mg/l, mkg/l */
        case CON_AQGPL: /*'d'*/
            Xincr = Xe * Vaq / DCmw;
            goto FINISH;
        case CON_AQMGPL: /*'e'*/
            Xincr = Xe * Vaq / DCmw / 1e3;
            goto FINISH;
        case CON_AQMKGPL: /*'b'*/
            Xincr = Xe * Vaq / DCmw / 1e6;
            goto FINISH;
        }
    if( fabs( Maq ) > 1e-12 && DCmw > 1e-12 )
        switch( UNITP )     /* Weight concentrations */
        {
        case CON_AQWFR:  /*'C'*/
            Xincr = Xe * Maq * 1e3 / DCmw;
            goto FINISH;
        case CON_AQWPROC:/*'c'*/
            Xincr = Xe * 10. * Maq / DCmw;
            goto FINISH;
        case CON_AQPPM:  /*'a'*/
            Xincr = Xe * Maq / 1e3 / DCmw;
            goto FINISH;
        }
    /* Error */
FINISH:
    return Xincr;
}

//---------------------------------------------------------------------------


// ------------------------------------------------------------------------------------------------
// Text summary of the task setup

namespace {

std::string num(double v, int prec = 4)
{
    std::ostringstream s;
    s << std::setprecision(prec) << v;
    return s.str();
}

// Split text into lines of at most `width` characters, breaking after ", " or at a space
std::vector<std::string> wrapLines(const std::string& text, size_t width)
{
    width = std::max<size_t>(width, 20);
    std::vector<std::string> lines;
    std::string rest = text;
    while(rest.size() > width) {
        size_t pos = rest.rfind(", ", width-1);
        if(pos != std::string::npos && pos > width/3)
            pos += 1;                                   // keep the comma on this line
        else {
            pos = rest.rfind(' ', width);
            if(pos == std::string::npos || pos < width/3) {
                // no break before width: do not cut a word, break at the next space instead
                pos = rest.find(' ', width);
                if(pos == std::string::npos)
                    break;
            }
        }
        std::string line = rest.substr(0, pos);
        while(!line.empty() && line.back() == ' ') line.pop_back();
        lines.push_back(line);
        rest = rest.substr(pos);
        rest.erase(0, rest.find_first_not_of(' ') == std::string::npos ? rest.size() : rest.find_first_not_of(' '));
    }
    if(!rest.empty() || lines.empty())
        lines.push_back(rest);
    return lines;
}

// Text wrapped to `width`; the first line starts with `indent`, the others with `cont`
std::string wrap(const std::string& text, size_t width, const std::string& indent = "", const std::string& cont = "    ")
{
    std::ostringstream o;
    const auto lines = wrapLines(text, width - std::max(indent.size(), cont.size()));
    for(size_t i=0; i<lines.size(); ++i)
        o << (i ? cont : indent) << lines[i] << "\n";
    return o.str();
}

// Table with left-aligned columns. If a row is wider than `width` (0: no limit), its last
// column continues on the next lines, below the start of that column
std::string table(const std::vector<std::vector<std::string>>& rows, const std::string& indent = "  ",
                  size_t width = 0)
{
    size_t ncol = 0;
    for(const auto& r: rows) ncol = std::max(ncol, r.size());
    std::vector<size_t> w(ncol, 0);
    for(const auto& r: rows)
        if(r.size() == ncol)                            // a shorter row ("... more groups") does not set widths
            for(size_t c=0; c<r.size(); ++c)
                w[c] = std::max(w[c], r[c].size());
    size_t start = indent.size();                       // start of the last column
    for(size_t c=0; c+1<ncol; ++c) start += w[c]+2;

    std::ostringstream o;
    for(const auto& r: rows) {
        std::string line = indent;
        for(size_t c=0; c+1<r.size(); ++c)
            line += r[c] + std::string(w[c] > r[c].size() ? w[c]-r[c].size()+2 : 2, ' ');
        const std::string last = r.empty() ? "" : r.back();
        if(r.size() < ncol || width == 0 || line.size() + last.size() <= width) {
            line += last;
            while(!line.empty() && line.back() == ' ') line.pop_back();
            o << line << "\n";
        }
        else if(width >= start + 40)                    // wrap the last column in its place
        {
            const auto lines = wrapLines(last, width - start);
            for(size_t i=0; i<lines.size(); ++i)
                o << (i ? std::string(start, ' ') : line) << lines[i] << "\n";
        }
        else {                                          // the last column below the row
            while(!line.empty() && line.back() == ' ') line.pop_back();
            o << line << "\n" << wrap(last, width, indent+"    ", indent+"    ");
        }
    }
    return o.str();
}

// Header lines: label, value, and a hint in brackets on the same line if it fits, else below the value
std::string headerBlock(const std::vector<std::vector<std::string>>& head, size_t width)
{
    size_t lw = 0;
    for(const auto& h: head) lw = std::max(lw, h[0].size());
    lw += 2;
    const std::string pad(lw, ' ');
    std::ostringstream o;
    for(const auto& h: head) {
        auto lines = wrapLines(h[1], width > lw + 20 ? width - lw : 40);
        for(size_t i=0; i<lines.size(); ++i)
            o << (i ? pad : h[0] + std::string(lw - h[0].size(), ' ')) << lines[i]
              << (i+1 < lines.size() ? "\n" : "");
        const std::string& hint = h.size() > 2 ? h[2] : std::string();
        if(hint.empty())
            o << "\n";
        else if(lw + lines.back().size() + 2 + hint.size() <= width)
            o << "  " << hint << "\n";
        else
            o << "\n" << wrap(hint, width, pad, pad + "  ");
    }
    return o.str();
}

std::string nodeRange(long int a, long int b)
{
    return a == b ? std::to_string(a) : std::to_string(a)+"-"+std::to_string(b);
}

} // namespace

std::string TGEM2MT::summary(size_t width, bool phases) const
{
    if(width < 60) width = 60;
    static const std::map<char, std::string> mode_text = {
        {RMT_MODE_S, "chain of reactors with waves of moving fluid or solids"},
        {RMT_MODE_F, "chain of flow-through reactors (transport, then equilibrium)"},
        {RMT_MODE_B, "box-flux system (fluxes of mobile phase groups between boxes)"},
        {RMT_MODE_A, "1-D advection + dispersion, explicit finite differences"},
        {RMT_MODE_C, "1-D advection with Courant number 1 (one node per step), no dispersion"},
        {RMT_MODE_W, "1-D advection + diffusion with random-walk particles"} };
    const char mode = mtp->PsMode;
    const bool fixed_ends = (mode == RMT_MODE_A || mode == RMT_MODE_C);
    const bool one_d = fixed_ends || mode == RMT_MODE_W;
    const long int nC = mtp->nC;
    const bool initialized = na && na->pNodT0() && nC > 0 && na->pNodT0()[0];
    const double tau0 = mtp->Tau[START_], tau1 = mtp->Tau[STOP_], tau_step = mtp->Tau[STEP_];

    std::ostringstream o;
    std::string title = name();       // bounded: the name array may have no terminating NUL
    while(!title.empty() && (title.back() == ' ' || title.back() == '`')) title.pop_back();
    o << "MassTransport" << (title.empty() ? std::string() : " \"" + title + "\"") << "\n";
    auto mt = mode_text.find(mode);
    o << "Mode " << mode << ": " << (mt != mode_text.end() ? mt->second : "?") << "\n";

    std::vector<std::vector<std::string>> head;
    if(one_d) {
        const double L = mtp->sizeLc[0];
        head.push_back({"Column", num(L)+" m, "+std::to_string(nC)+" nodes, node length "+num(L/nC)+" m",
                        "[setSpatialDimensions(x), MassTransport(mode, nnodes)]"});
        head.push_back({"Velocity", num(mtp->fVel)+" m/s (pore velocity, all nodes)", "[setFluidAdvectionVelocity(velocity)]"});
    }
    else
        head.push_back({"Nodes", std::to_string(nC), "[MassTransport(mode, nnodes)]"});
    head.push_back({"Time", num(tau0, 8)+" .. "+num(tau1, 8)+", step "+num(tau_step, 8), "[setTau(start, stop, step)]"});

    std::string warning;
    if(fixed_ends) {
        const double dx = mtp->sizeLc[0]/nC;
        const double tf = (mode == RMT_MODE_A ? std::max(mtp->tf, 1.) : (mtp->tf > 0 ? mtp->tf : 1.));
        double dt = std::numeric_limits<double>::infinity(), dt_adv = dt, dt_disp = dt;
        std::string limit;
        if(mode == RMT_MODE_A) {
            if(mtp->fVel > 0.) dt_adv = 0.5*dx/mtp->fVel;
            double Dmax = 0.;
            for(long int q=0; q<nC; ++q) Dmax = std::max(Dmax, DispersionCoef(q));
            if(Dmax > 0.) dt_disp = 0.25*dx*dx/Dmax;
            dt = std::min(dt_adv, dt_disp);
            limit = (dt_disp < dt_adv ? "dispersion" : "advection");
            limit += "; advection limit " + (std::isfinite(dt_adv) ? num(dt_adv)+" s" : std::string("none")) +
                     ", dispersion limit " + (std::isfinite(dt_disp) ? num(dt_disp)+" s" : std::string("none"));
        }
        else if(mtp->fVel > 0.) {
            dt = dx/mtp->fVel;
            limit = "node length / velocity";
        }
        if(std::isfinite(dt)) {
            dt /= tf;
            const double span = tau1 - tau0;
            const long int steps = static_cast<long int>(std::ceil(span/dt - 1e-9));
            const long int max_steps = static_cast<long int>(span/tau_step) + 1;
            head.push_back({"Time step", num(dt)+" s ("+limit+"), reduction factor "+num(tf),
                            "[set by the program; divided by setTimeStepReductionFactor(factor)]"});
            head.push_back({"Steps", "about "+std::to_string(std::min(steps, max_steps))+", to t = "+
                            num(tau0+std::min(steps, max_steps)*dt, 8), "[stop time and setTau step]"});
            if(max_steps < steps)
                warning = "WARNING: the run stops after " + std::to_string(max_steps) + " steps (t = " +
                          num(tau0+max_steps*dt) + "), before the stop time. Make the setTau step smaller than the time step.\n";
        }
        else
            warning = "WARNING: the velocity and the dispersion are 0: no time step (runMassTransport will fail).\n";
    }

    // Initial systems
    std::string systems;
    if(initialized) {
        for(size_t i=0; i<dbr_file_names.size(); ++i)
            systems += (i ? ",  " : "") + std::to_string(i) + " = " + dbr_file_names[i];
        if(dbr_file_names.size() == 1)
            systems += "  (one file: all nodes start from it; the system numbers below are groups)";
    }
    else
        systems = "not read yet";
    head.push_back({"Initial systems", systems, "[initNodes(dat_lst) reads the dbr files in their order]"});
    o << headerBlock(head, width);
    if(!warning.empty()) {
        while(!warning.empty() && warning.back() == '\n') warning.pop_back();
        o << wrap(warning, width, "", "         ");
    }

    // Per node description; consecutive nodes with an equal description form one group
    auto role = [&](long int q) -> std::string {
        if(fixed_ends) {
            if(q <= 1) return "inflow, fixed";
            if(q == nC-1) return "outlet, fixed";
            return "column";
        }
        switch(mtp->DiCp[q][1]) {
        case 0: return one_d ? "column" : "box";
        case 1: case 2: case 3: return "source";
        case -1: case -2: case -3: return "sink";
        case 4: return "input function";
        default: return "?";
        }
    };
    auto transport = [&](long int q) {
        double eps = mtp->eps_in, al = mtp->al_in, Dif = mtp->Dif_in, nto = mtp->nto_in, vp = mtp->fVel;
        if(mtp->HydP) {
            vp = mtp->HydP[q][1]; eps = mtp->HydP[q][2]; al = mtp->HydP[q][4];
            Dif = mtp->HydP[q][5]; nto = mtp->HydP[q][6];
        }
        return std::vector<std::string>{num(eps), num(al), num(Dif), num(nto), num(vp)};
    };
    bool own_velocity = false;
    if(mtp->HydP)
        for(long int q=0; q<nC; ++q)
            if(mtp->HydP[q][1] != mtp->fVel) own_velocity = true;

    // Phases present in any node (DBR phase indices)
    std::vector<long int> present;
    if(initialized) {
        const DATACH* ch = na->pCSD();
        for(long int k=0; k<ch->nPHb; ++k)
            for(long int q=0; q<nC; ++q)
                if(na->pNodT1()[q]->xPH[k] > 1e-12) { present.push_back(k); break; }
    }

    auto describe = [&](long int q) {
        std::vector<std::string> r{std::to_string(mtp->DiCp[q][0]), role(q), std::to_string(mtp->DiCp[q][1])};
        if(one_d) {
            auto t = transport(q);
            r.insert(r.end(), t.begin(), t.begin()+4);
            if(own_velocity) r.push_back(t[4]);
        }
        if(initialized) {
            const DATABR* n = na->pNodT1()[q];   // equilibrated state (initNodes, setNodes, each step)
            r.push_back(num(n->TK, 5));
            r.push_back(num(n->P/1e5));
            r.push_back(num(n->Ms));
            r.push_back(num(n->Vs*1e3));
            if(phases)                           // amounts only when asked: nodes then group by them too
                for(auto k: present)
                    r.push_back(n->xPH[k] > 1e-12 ? num(n->xPH[k]) : "-");
        }
        return r;
    };

    std::vector<std::string> hdr{"system", "role", "type"};
    if(one_d) {
        for(auto h: {"porosity", "dispersivity m", "diffusivity m2/s", "tortuosity"}) hdr.push_back(h);
        if(own_velocity) hdr.push_back("velocity m/s");
    }
    const size_t n_setup = hdr.size();     // the columns before T, P, mass, volume
    size_t n_chem = n_setup;               // the columns before the phase amounts
    if(initialized) {
        for(auto h: {"T K", "P bar", "mass kg", "volume L"}) hdr.push_back(h);
        n_chem = hdr.size();
        if(phases)
            for(auto k: present) hdr.push_back(na->pCSD()->PHNL[na->pCSD()->xph[k]]);
    }
    std::vector<std::pair<std::pair<long int,long int>, std::vector<std::string>>> groups;
    for(long int q=0; q<nC; ++q) {
        auto d = describe(q);
        if(!groups.empty() && groups.back().second == d)
            groups.back().first.second = q;
        else
            groups.push_back({{q, q}, d});
    }
    // Two tables: the node setup, and the chemical state (if the nodes are read)
    const size_t max_rows = 30;
    auto one_table = [&](size_t c0, size_t c1) {
        std::vector<std::vector<std::string>> t{{"nodes"}};
        t[0].insert(t[0].end(), hdr.begin()+c0, hdr.begin()+c1);
        for(size_t g=0; g<groups.size(); ++g) {
            if(groups.size() > max_rows && g == max_rows-5) {
                t.push_back({"...", "("+std::to_string(groups.size()-max_rows+1)+" more groups)"});
                g = groups.size()-5;
            }
            std::vector<std::string> r{nodeRange(groups[g].first.first, groups[g].first.second)};
            r.insert(r.end(), groups[g].second.begin()+c0, groups[g].second.begin()+c1);
            t.push_back(r);
        }
        return t;
    };
    // Columns [c0, c1) as one or more tables that fit the width (the nodes column is repeated);
    // a last text column (phases) may wrap
    auto node_table = [&](size_t c0, size_t c1) {
        auto full = one_table(c0, c1);
        std::vector<size_t> w(c1-c0+1, 0);
        for(const auto& r: full)
            for(size_t c=0; c<r.size() && r.size() == w.size(); ++c)
                w[c] = std::max(w[c], r[c].size());
        std::string out;
        size_t c = c0;
        while(c < c1) {
            size_t used = 2 + w[0] + 2, e = c;
            while(e < c1 && (e == c || used + w[e-c0+1] <= width)) {
                used += w[e-c0+1] + 2;
                ++e;
            }
            if(!out.empty()) out += "\n";
            out += table(one_table(c, e), "  ", width);
            c = e;
        }
        return out;
    };

    // Simple picture of the column: the groups along the flow
    if(one_d && groups.size() <= 8) {
        std::vector<std::vector<std::string>> pic{{"flow ->"}, {"nodes"}, {"role"}, {"system"}};
        for(const auto& gr: groups) {
            pic[1].push_back("| "+nodeRange(gr.first.first, gr.first.second));
            pic[2].push_back("| "+gr.second[1]);
            pic[3].push_back("| "+gr.second[0]);
        }
        for(size_t i=1; i<pic.size(); ++i) pic[i].push_back("|");
        o << "\n" << table(pic, "  ", width);
    }
    o << "\n" << node_table(0, n_setup);
    if(initialized)
    {
        o << "\n" << node_table(n_setup, n_chem);
        if(phases && !present.empty())
            o << "\nPhases (mol), \"-\" = not present:\n" << node_table(n_chem, hdr.size());
        else if(!present.empty()) {
            std::string names;
            for(auto k: present)
                names += (names.empty() ? "" : ", ") + na->pCSD()->PHNL[na->pCSD()->xph[k]];
            if(names.size() > width/2)
                names = names.substr(0, names.rfind(", ", width/2)) + ", ...";
            o << wrap("Phases present: " + std::to_string(present.size()) + " (" + names +
                      ")  [amounts: summary(phases=True)]", width, "", "    ");
        }
    }

    if(!one_d) {
        auto id = [](const char* p) {
            std::string t(p, strnlen(p, MAXSYMB));
            while(!t.empty() && (t.back() == ' ' || t.back() == '`')) t.pop_back();
            return t;
        };
        const DATACH* ch = initialized ? na->pCSD() : nullptr;
        // Mobile phase groups: the phases (with quantities) that a flux moves
        o << "\nMobile phase groups: " << mtp->nPG << "\n"
          << wrap("[setNumberPhaseGroups(n), setPhaseGroupID(g, id), setPhaseGroupQuantities(g, phase, q), "
                  "setUnitsPhaseQuantities(phase, unit)]", width, "  ", "   ");
        if(mtp->nPG > 0 && mtp->MGPid && mtp->PGT) {
            std::vector<std::vector<std::string>> g_rows{{"group", "id", "phases (quantity unit)"}};
            for(long int g=0; g<mtp->nPG; ++g) {
                std::string ph;
                for(long int f=0; f<mtp->FIf; ++f) {
                    const double qv = mtp->PGT[f*mtp->nPG+g];
                    if(qv == 0.) continue;
                    std::string name = (ch && f < ch->nPHb) ? ch->PHNL[ch->xph[f]] : "phase "+std::to_string(f);
                    ph += (ph.empty() ? "" : ", ") + name + " " + num(qv) +
                          (mtp->UMGP ? " " + std::string(1, mtp->UMGP[f]) : "");
                }
                g_rows.push_back({std::to_string(g), id(mtp->MGPid[g]), ph.empty() ? "-" : ph});
            }
            o << table(g_rows, "  ", width);
        }
        // Flux definitions
        o << "Fluxes: " << mtp->nFD << "\n"
          << wrap("[setNumberMGPfluxes(n), setFluxSourceReceive(i, from, to), setFluxDefinition(i, order, rate, "
                  "quantity, value), setFluxMGPid(i, id), setFluxID(i, id)]", width, "  ", "   ");
        if(mtp->nFD > 0 && mtp->FDLi && mtp->FDLf) {
            // Consecutive fluxes along a chain (i -> i+1) with equal settings form one row
            struct Flux { long int i0, i1, from0, from1, to0, to1; std::vector<std::string> rest; };
            std::vector<Flux> fl;
            for(long int i=0; i<mtp->nFD; ++i) {
                const long int from = mtp->FDLi[i][0], to = mtp->FDLi[i][1];
                std::vector<std::string> rest{mtp->FDLid ? id(mtp->FDLid[i]) : "", mtp->FDLmp ? id(mtp->FDLmp[i]) : "",
                                              num(mtp->FDLf[i][0]), num(mtp->FDLf[i][1]), num(mtp->FDLf[i][2])};
                if(!fl.empty() && fl.back().rest == rest && fl.back().i1 == i-1 &&
                   from >= 0 && to >= 0 && fl.back().from1 >= 0 && fl.back().to1 >= 0 &&
                   fl.back().from1 == from-1 && fl.back().to1 == to-1) {
                    fl.back().i1 = i; fl.back().from1 = from; fl.back().to1 = to;
                }
                else
                    fl.push_back({i, i, from, from, to, to, rest});
            }
            auto node = [](long int a, long int b) {
                return a < 0 ? std::string("outside") : nodeRange(a, b);
            };
            std::vector<std::vector<std::string>> f_rows{{"fluxes", "id", "from node", "to node", "moves", "order", "rate", "quantity"}};
            for(const auto& f: fl)
                f_rows.push_back({nodeRange(f.i0, f.i1), f.rest[0], node(f.from0, f.from1), node(f.to0, f.to1),
                                  f.rest[1], f.rest[2], f.rest[3], f.rest[4]});
            o << table(f_rows, "  ", width);
            if(fl.size() < static_cast<size_t>(mtp->nFD))
                o << wrap("(a row \"1-49, from 1-49 to 2-50\" means: flux 1 goes from node 1 to 2, flux 2 from 2 to 3, ...)",
                          width, "  ", "   ");
        }
        if(mtp->nSFD > 0)
            o << "Elemental source fluxes: " << mtp->nSFD << "\n"
              << wrap("[setNumberICsourceFluxes(n), setICsourceQuantities(i, element, q)]", width, "  ", "   ");
    }

    if(gem_failures.ran) {
        const auto& f = gem_failures;
        o << "\n";
        if(f.count == 0)
            o << "Last run: all GEM calculations in the nodes converged.\n";
        else {
            std::string nodes;
            long int a = -2, b = -2;
            auto flush = [&]() { if(a >= 0) nodes += (nodes.empty() ? "" : ", ") + nodeRange(a, b); };
            for(auto q: f.nodes) {
                if(q == b+1) { b = q; continue; }
                flush(); a = b = q;
            }
            flush();
            o << wrap("WARNING: last run: " + std::to_string(f.count) + " GEM calculations failed or are not trustful, in nodes " +
                      nodes + ". First: step " + std::to_string(f.first_step) + ", t = " + num(f.first_time) + ", node " +
                      std::to_string(f.first_node) + ": " + f.first_message + "  [gemFailures(); see xGEM2MT.log]", width, "", "         ");
        }
    }

    std::vector<std::vector<std::string>> how{
        {"system of a node", "before initNodes", "setDistributing(node, system)"},
        {"composition, T, P", "after initNodes", "setSystem(system, composition=..., T=..., P=...), setNodes(nodes, ...)"}};
    if(one_d) {
        how.push_back({"porosity, dispersivity,", "before initNodes:", "setInitialNodeEffectivePorosity(v), setInitialDispersivity(v), ... (all nodes)"});
        how.push_back({"diffusivity, tortuosity", "after initNodes:", "setSystem/setNodes(..., transport=TransportProperties(porosity=..., ...))"});
    }
    how.push_back({"node type", fixed_ends ? "not used in this mode" : "", "setNodeType(node, type), setNodesType(nodes, type)"});
    o << "\nHow to change:\n" << table(how, "  ", width);
    if(fixed_ends)
        o << wrap(std::string("Mode ")+mode+": nodes 0 and 1 (inflow) and the last node (outlet) keep their composition during the run.",
                  width);
    if(!initialized)
        o << "Compositions are shown after initNodes.\n";
    return o.str();
}
