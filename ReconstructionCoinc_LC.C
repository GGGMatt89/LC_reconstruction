 // fichier qui reconstruit les coincidences trouvees par la fonction Coincidence.C
#include "TH1F.h"
#include "TF1.h"
#include "TMath.h"
#include "TFrame.h"
#include "TFile.h"
#include <iostream>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "Riostream.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include <TGraph.h>

#include "ReconstructionCC.h"
#include "ComptonFitHisto.h"
#include "Fit.h"


Double_t IsCloser(Double_t y1,Double_t y2,Double_t ye)
{

    //cout<<"iscloser ye "<<ye<< " y1 "<<y1<<" y2 "<<y2<<endl;
    Double_t difference_1=0;
    Double_t difference_2=0;


    difference_1= TMath::Abs(y1-ye);
    difference_2= TMath::Abs(y2-ye);
    if(TMath::Abs(difference_1)<TMath::Abs(difference_2))return y1;
    //if(TMath::Abs(y1)<TMath::Abs(y2)) return y1;
    else return y2;

}



//protons : 3,4,5,10,15,1,2
//protons deuxieme batch: 20,30,50,80,100,150,200,217

void ReconstructionCoinc(char* name, Int_t nbin, Int_t fileindex)
{

    string num_file;

    ostringstream obuffer1;
    obuffer1<<((Double_t)fileindex);
    num_file=obuffer1.str();

    string test = name;
    //string outpath = "/sps/hep/hadronth/mfontana/results_CCanalysis/";
    string outpath = "./results/";
    string  strname= outpath + test+"_2hit_";
    strname+=num_file;
    strname+=".root";

    /*string strname = name;
    strname+="_2hit.root";*/



    //-------------------oooooOOOOO00000OOOOOooooo---------------------#
    //                                                                 #
    //      	DEFINITION DU TYPE DE DONNEES 			   #
    //       	 			   			   #
    //                                                                 #
    //-------------------oooooOOOOO00000OOOOOooooo---------------------#

    //CUT TOF OR NOT

    Int_t cut =1;				// If cut = 0 -> NO cut TOF,   if cut = 1 -> Cut TOF
    Int_t save_data =1;		// If save_data = 0 -> no data fils create,   If save_data = 1 -> data fils create


    Int_t limit_high_cut_TOF =0;
    //Donnees fichiers root
    Double_t Nion = 100*pow(10.,6);
    TFile* file;
    TTree * tree;
    Int_t nentries;
    char* filename;
    string newstrname;
    string newstrname_MLEM;
    char* newname;
    Int_t num_coinc =0;

    Double_t energy_initiale =0;

    //Donnees coincidences
    Coinc mycoinc;
    Double_t edep1,edep2,x1,y1,z1,x2,y2,z2,time1,time2;
    bool samepscoinc = false;
    bool gammacoinc = false ;
    bool trueeventcoinc = false;
    Double_t ncoinc=0.;
    Int_t Coinc_number_per_event=0;
    Double_t nsamepscoinc=0.;
    Double_t ngammacoinc=0.;
    Double_t ntrueeventcoinc=0.;
    string gammatype,gammatype2, particletypesi,particletypelyso,creatorprocess;

    Int_t ngammatrue_ym =0;
    Int_t nvraie_ym =0;
    Int_t nrandom_ym =0;
    Double_t sum_edep=0;
    Double_t sum_edep_95=0;
    Double_t minimum_value =0;
    Int_t position_value[31];
    for (int m =0;m<31;m++) position_value[m]= -300+20*m;
    Double_t difference_position =0;
    Int_t up[2] ={0,0};
    Double_t weight[2]={0,0};

    Int_t test_Sum=0;
    Double_t sum_CC=0;
    Int_t Nion_Si = 0;
    Int_t Nion_LYSO = 0;
    Int_t n_fortuits =0;
    Int_t n_vrai =0;
    Int_t n_null=0;
    Int_t nombre_coinc_reconstruite =0;
    Int_t Sameparticule_number =0;
    Int_t Diffparticule_number =0;
    Int_t Sameparticule_number_histo =0;
    Int_t Diffparticule_number_histo =0;
    bool Sameparticule = false ;
    bool Diffparticule = false ;
    bool fortuitparticule = false;

    Double_t ne_moins=0;
    Double_t nemoins_sum=0;

    Double_t nprotons=0;
    Double_t nprotons_sum=0;

    Double_t ne_plus=0;
    Double_t neplus_sum=0;

    Double_t n_neutrons=0;
    Double_t n_neutrons_sum=0;

    string particule_init;
    string type_particule_init;


    //TFile *f = new TFile("test.root","RECREATE");

    Int_t count_gamma_out=0;
    bool gamma_out=false;
    Int_t count_sumAccurate = 0;

    //Donn�es geometrie de la camera
    Double_t d0,d1,d2,l1,e1,l2,e2;	//en mm
    d0=200;//100;
    d1=10;
    d2=400;
    l1=90;//80;
    e1=2;
    l2=380;//300;
    e2=30;//25;
    Double_t c=3*pow(10.,8)*pow(10.,3)/pow(10.,9);	// en mm.ns-1

    Int_t neventstore =0;
    Int_t verif_otherevents=0;

    //Resultat reconstruction
    Double_t* y ;
    Double_t ym;//ym[2]
    Double_t tofSi,tofLYSO;
    char ParticleTypeSi[64];
    char ParticleTypeLYSO[64];
    char gammaTypeSi[64];
    char gammaTypeLYSO[64];
    char ProcessSi[64];
    char ProcessLYSO[64];
    Int_t nsol;
    Double_t ym_background[2];
    Double_t ym_TrueGamma [2];
    Double_t ym_all [2];
    Double_t ym_sameparticle[2];
    Double_t ym_different[2];
    Double_t ym_fortuit[2];
    Double_t ym_gamma_prompt[2];
    Double_t ym_gamma_prompt_scatt[2];
    Double_t ym_gamma_other_value[2];

    Int_t num_gamma_init =0;
    Int_t num_gamma_scatt_init=0;
    Int_t num_gamma_prompt_ini =0;
    Int_t num_gamma_other =0;


    Int_t num_gamma_recons =0;
    Int_t num_gamma_scatt_recons=0;
    Int_t num_gamma_prompt_recons =0;
    Int_t num_gamma_other_recons =0;



    //Histo diff�rence solutions/au vertex
    Double_t Histo_difference_y0[100000];
    Double_t Histo_difference_y1[100000];
    Double_t Histo_difference_y0y1[100000];
    Double_t moyenne_Histo_difference_y0=0;
    Double_t moyenne_Histo_difference_y1=0;
    Double_t moyenne_Histo_difference_y0y1=0;


    //Hsto diff�rence entre solutions
    Double_t difference_y0y1=0;
    Double_t difference_y0y1_entertarget=0;
    Double_t difference_y0y1_outtarget=0;
    Double_t difference_y0y1_fallofftarget=0;
    Double_t difference_y0y1_cut0_1MeV=0;
    Double_t difference_y0y1_cut1_3MeV=0;
    Double_t difference_y0y1_cut3MeV_sup=0;
    Double_t difference_y0y1_cut1_3MeV_edep=0;
    Double_t difference_y0y1_cut3MeV_sup_edep=0;


    Double_t difference_weight =0;


    //efficacit� absolue de la CC

    Double_t value_efficacite_1MeV[31]= {9.98e-05, 0.0001083, 0.0001177, 0.0001343, 0.0001523, 0.0001749, 0.0001926, 0.00021, 0.0002345, 0.0002564, 0.0002777, 0.000305, 0.0003166, 0.0003225, 0.0003308, 0.0003362, 0.0003339, 0.000318, 0.000309, 0.0002983, 0.0002758, 0.0002542, 0.0002288, 0.0002169, 0.0001925, 0.0001701, 0.0001502, 0.0001388, 0.0001177, 0.0001037, 9.25e-05};


    Double_t value_efficacite_2MeV[31]={7.84e-05,8.71e-05, 0.0001038, 0.0001139, 0.0001298, 0.0001443, 0.0001706, 0.0001899, 0.0002212, 0.0002386, 0.0002715, 0.0002917, 0.0003223, 0.0003367, 0.0003481, 0.0003482, 0.0003492, 0.0003313, 0.0003144, 0.0002879, 0.0002688, 0.0002412, 0.0002164, 0.0001941, 0.0001646, 0.0001444, 0.0001351, 0.0001218, 0.0001007, 9.06e-05, 8.22e-05};

    Double_t value_efficacite_3MeV[31]={9.24e-05, 0.0001051, 0.000107, 0.0001333, 0.0001521, 0.0001665, 0.0001924, 0.0002134, 0.0002556, 0.0002796, 0.0003044, 0.0003349, 0.0003693, 0.0003914, 0.0004119, 0.0004231, 0.0004131, 0.0003918, 0.0003787, 0.0003295, 0.0003034, 0.0002811, 0.0002493, 0.0002161, 0.0001892, 0.0001623, 0.0001448, 0.0001256, 0.0001073, 9.77e-05, 9.25e-05};

    Double_t value_efficacite_4MeV[31]={0.0001046, 0.0001188, 0.0001235, 0.0001463, 0.0001726, 0.0001929, 0.000218, 0.0002467, 0.0002882, 0.0003335, 0.0003704, 0.0004142, 0.0004482, 0.0004842, 0.0005098, 0.0005283, 0.0005141, 0.0004835, 0.000462, 0.000423, 0.0003755, 0.0003334, 0.0002788, 0.0002496, 0.0002285, 0.0001938, 0.0001735, 0.0001484, 0.0001295, 0.0001076, 9.96e-05};


    Double_t value_efficacite_5MeV[31]={0.0001208, 0.0001391, 0.0001462, 0.0001698, 0.0001976, 0.0002335, 0.0002673, 0.000309, 0.0003515, 0.0004061, 0.0004552, 0.0005253, 0.0005657, 0.0006208, 0.0006305, 0.0006359, 0.0006376, 0.0006031, 0.0005477, 0.0005214, 0.0004533, 0.0004093, 0.0003422, 0.0002947, 0.0002628, 0.000241, 0.0001964, 0.0001735, 0.0001478, 0.00013, 0.0001158};

    Double_t value_efficacite_6MeV[31]={0.0001295, 0.0001445, 0.000172, 0.000198, 0.0002308, 0.0002626, 0.0003029, 0.0003599, 0.0004074, 0.000484, 0.0005392, 0.000621, 0.0006653, 0.0007635, 0.0007912, 0.0007962, 0.0007684, 0.0007374, 0.0006788, 0.0006085, 0.0005472, 0.0004795, 0.0004139, 0.0003693, 0.0003113, 0.0002637, 0.0002236, 0.0001972, 0.0001696, 0.0001509, 0.0001345};

    for( int r=0; r<100000;r++)
    {
        Histo_difference_y0[r]=0;
        Histo_difference_y1[r]=0;
        Histo_difference_y0y1[r]=0;
    }

    nbin=68;//68;
    Int_t limitym_histo_inf=-180;
    Int_t limitym_histo_sup=160;

    TH1* ym_histo = new TH1F("ym","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_prompt = new TH1F("ymprompt","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_true = new TH1F("ymtrue","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_true_same_gamma = new TH1F("ymtrue_same_gamma","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_false = new TH1F("ymfalse","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_true_same_particle_no_gamma = new TH1F("ymtruesameparticle","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_histo_true_different_particle = new TH1F("ymtruediffparticle","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_vertex = new TH1F("y vertex","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_gamma_scattered = new TH1F("ym_gamma_scattered","",nbin,limitym_histo_inf,limitym_histo_sup);
    TH1* ym_gamma_other = new TH1F("ym_gamma_other","",nbin,limitym_histo_inf,limitym_histo_sup);

    TH1* ym_histo_autre_solution = new TH1F("ym_histo_autre","",nbin,limitym_histo_inf,limitym_histo_sup);

    TH1* spectre_energy_protons = new TH1F("spectre_energy_protons","",500,0,100);
    TH1* spectre_energy_neutrons = new TH1F("spectre_energy_neutrons","",500,0,100);
    TH1* spectre_energy_electrons = new TH1F("spectre_energy_electrons","",500,0,100);

    nbin=1000;

    TH1* ye_histo = new TH1F("ye","",nbin,-500,500);
    TH1* ye_histo_true = new TH1F("yetrue","",nbin,-500,500);
    TH1* ye_histo_false = new TH1F("yefalse","",nbin,-500,500);

    TH1* xe_histo = new TH1F("xe","",nbin,-500,500);
    TH1* ze_histo = new TH1F("ze","",nbin,-500,500);

    nbin=1000;


    TH1* tof_histo = new TH1F("tof","",nbin,-5,20);
    TH1* tof_histo_true = new TH1F("toftrue","",nbin,-5,20);
    TH1* tof_histo_false = new TH1F("toffalse","",nbin,-5,20);
    TH1* tof_histo_samegamma = new TH1F("toffalse","",nbin,-5,20);
    TH1* tof_histo_sameparticle= new TH1F("toffalse","",nbin,-5,20);

    TH1* n_histo = new TH1F("n","",nbin,-500,0);
    TH1* n_histo_true = new TH1F("ntrue","",nbin,-500,0);
    TH1* n_histo_false = new TH1F("nfalse","",nbin,-500,0);

    TH2F* type_histo = new TH2F("type","",3,0,3,2,0,2);
    type_histo->SetBit(TH1::kCanRebin);
    type_histo->SetStats(0);
    TH2F* type_histo_true = new TH2F("type_true","",3,0,3,2,0,2);
    type_histo_true->SetBit(TH1::kCanRebin);
    type_histo_true->SetStats(0);
    TH2F* type_histo_false = new TH2F("type_false","",3,0,3,2,0,2);
    type_histo_false->SetBit(TH1::kCanRebin);
    type_histo_true->SetStats(0);

    TH2F* gammatype_histo = new TH2F("gammatype","",3,0,3,2,0,2);
    gammatype_histo->SetBit(TH1::kCanRebin);
    gammatype_histo->SetStats(0);
    TH2F* gammatype_histo_true = new TH2F("gammatype_true","",3,0,3,2,0,2);
    gammatype_histo_true->SetBit(TH1::kCanRebin);
    gammatype_histo_true->SetStats(0);
    TH2F* gammatype_histo_false = new TH2F("gammatype_false","",3,0,3,2,0,2);
    gammatype_histo_false->SetBit(TH1::kCanRebin);
    gammatype_histo_true->SetStats(0);

    TH2F* gammaproc_histo_false = new TH2F("gammaproc_false","",3,0,3,2,0,2);
    gammaproc_histo_false->SetBit(TH1::kCanRebin);
    gammaproc_histo_false->SetStats(0);

    TH2F* process_histo = new TH2F("process","",3,0,3,2,0,2);
    process_histo->SetBit(TH1::kCanRebin);
    process_histo->SetStats(0);
    TH2F* process_histo_true = new TH2F("process_true","",3,0,3,2,0,2);
    process_histo_true->SetBit(TH1::kCanRebin);
    process_histo_true->SetStats(0);
    TH2F* process_histo_false = new TH2F("process_false","",3,0,3,2,0,2);
    process_histo_false->SetBit(TH1::kCanRebin);
    process_histo_true->SetStats(0);


    TH2F* pp_histo = new TH2F("pp","",3,0,3,2,0,2);
    pp_histo->SetBit(TH1::kCanRebin);
    pp_histo->SetStats(0);

    TH1 * nion1_histo = new TH1F("nion1","",20000000,0.,20000000);
    TH2 * nion_histo = new TH2F("nion","",nbin,0.,20000000.,nbin,0.,20000000.);
    TH2 * nion_histo_true = new TH2F("niontrue","",nbin,0.,20000000.,nbin,0.,20000000.);
    TH2 * nion_histo_false = new TH2F("nionfalse","",nbin,0.,20000000.,nbin,0.,20000000.);

    nbin=200;
    TH2 * edep_histo = new TH2F("edep gamma","",nbin,0.,4.,nbin,0.,4.);
    TH2 * edep_histo_true = new TH2F("edep prompt no scatt","",nbin,0.,4.,nbin,0.,4.);
    //TH1 * edep_histo_true = new TH1F("edeptrue","",nbin,0.,100.);
    TH2 * edep_histo_false = new TH2F("edep prompt scatt","",nbin,0.,4.,nbin,0.,4.);
    TH2 * edep_histo_true_gamma = new TH2F("edep gamma other","",nbin,0.,4.,nbin,0.,4.);

     TH1F* spectre_energy_gamma = new TH1F("spectre_energy_gamma","",140,0,7);

    TH1F * histo_yVertex_y1y0 = new TH1F("vertex","",100,0.,500.);
    TH1F * histo_yVertex_y1 = new TH1F("vertex","",100,-300.,500.);
    TH1F * histo_yVertex_y0 = new TH1F("vertex","",100,-300.,300.);


    TH1F * histo_poisson = new TH1F("poisson","",200,0,20);

    //Distance entre 2 solutions
    TH1F * histo_distancey0y1 = new TH1F("distancey0y1","",100,0.,500.);
     TH1F * histo_distancey0y1_enterTarget = new TH1F("distancey0y1_enterTarget","",100,0.,500.);
    TH1F * histo_distancey0y1_outTarget = new TH1F("distancey0y1_outTarget","",100,0.,500.);
    TH1F * histo_distancey0y1_FallOffTarget = new TH1F("distancey0y1_FallOffTarget","",100,0.,500.);
    TH1F * histo_distancey0y1_cut0_1MeV = new TH1F("distancey0y1_0_1MeV","",100,0.,500.);
    TH1F * histo_distancey0y1_cut1MeV_3MeV = new TH1F("distancey0y1_1_3MeV","",100,0.,500.);
    TH1F * histo_distancey0y1_cut3MeV_sup = new TH1F("distancey0y1_3MeV_sup","",100,0.,500.);
    TH1F * histo_distancey0y1_cut1MeV_3MeV_edep = new TH1F("distancey0y1_1_3MeV_edep","",100,0.,500.);
    TH1F * histo_distancey0y1_cut3MeV_sup_edep = new TH1F("distancey0y1_3MeV_sup_edep","",100,0.,500.);

    //Difference weight ponderation
     TH1F * histo_difference_weight = new TH1F("histo_difference_weight","",100,0.,1.1);


    //-------------------oooooOOOOO00000OOOOOooooo---------------------#
    //                                                                 #
    //   		   DEBUT DU PROGRAMME D'ANALYSE			   #
    //                                                                 #
    //-------------------oooooOOOOO00000OOOOOooooo---------------------#

    //Enregistrement du resultat de la reconstruction dans un fichier root

    //if(cut ==0)newstrname="/sps/hep/hadronth/JLL/Macro_Root/Coincidences_reelles_fortuites/code/code2/Shift_camera_study/Results_Shift_plus5_cm/2015_03_24_Protons_160MeV_PMMACylinder_Reconstructed_profile_ALL_NoCutTOF_shift_camera_plus5cm_"+num_file;
    //if(cut ==1) newstrname="/sps/hep/hadronth/JLL/Macro_Root/Coincidences_reelles_fortuites/code/code_protons2/Shift_camera_study/Results_Shift_moins5_cm/2015_03_31_Protons_160MeV_PMMACylinder_Reconstructed_profile_ALL_CutTOF_amera_moins5cm_"+num_file;
    if(cut==0)newstrname =  outpath + name +"noCut";
    if(cut==1)newstrname =  outpath + name +"Cut";
    newstrname+=".root";
    //2014_12_15_Protons_160MeV_PMMACube_Reconstructed_profile_ALL_NoCutTof_
    //2014_12_15_Protons_160MeV_PMMACube_Reconstructed_profile_All_CutTof_
    //2014_12_15_Protons_160MeV_PMMACylinder_Reconstructed_profile_ALL_NoCutTof_
    //2014_12_15_Protons_160MeV_PMMACylinder_Reconstructed_profile_All_CutTof_
    //2014_12_27_Protons_160MeV_PMMACube20x20cm_Reconstructed_profile_ALL_NoCutTof_d020cmtest
    //2015_01_16_Protons_160MeV_PMMACylinder_Reconstructed_profile_ALL_Show_Gamma_NewFiles_CutTof_TESTT

    TFile* ReconstructionTree = new TFile(newstrname.c_str(),"RECREATE");
    TTree* reconstree=new TTree("T_Recons","A Root Tree");
    reconstree->Branch("B_ym_background",ym_background,"ym_background[2]/D");
    reconstree->Branch("B_ym_TrueGamma",ym_TrueGamma,"ym_TrueGamma[2]/D");
    reconstree->Branch("B_ym_sameparticle",ym_sameparticle,"ym_sameparticle[2]/D");
    reconstree->Branch("B_ym_different",ym_different,"ym_different[2]/D");
    reconstree->Branch("B_ym_fortuit",ym_fortuit,"ym_fortuit[2]/D");
    reconstree->Branch("B_ym_gamma_prompt",ym_gamma_prompt,"ym_gamma_prompt[2]/D");
    reconstree->Branch("B_ym_ym_gamma_prompt_scatt",ym_gamma_prompt_scatt,"ym_gamma_prompt_scatt[2]/D");
    reconstree->Branch("B_ym_gamma_other_value",ym_gamma_other_value,"ym_gamma_other_value[2]/D");
    reconstree->Branch("B_ym_weight",weight,"weight[2]/D");


    TTree* reconstree_all=new TTree("T_Recons_all","A Root Tree");
    reconstree_all->Branch("B_ym_all",ym_all,"ym_all[2]/D");
    // rectree->Branch("B_part",&ParticleTypeSi,"ParticleTypeSi/C");




    // -------------------------------------------------------------------------------
    // TEST JLL 22/01/2013

    Int_t n=100;
    string filename2;
    char * filenamechar2;
    string* numbern = new string[n];
    bool isfile;
    vector<Int_t> Nion_ini;
    vector<double_t> t0_ini;
    vector<Int_t> ID_ini;
    Int_t nentries2;
    Int_t nombre_nion =0;
    double_t nentries_total =0;
    double_t T0_diff;


    /*    for(Int_t k=0;k<2;k++)
     {
     cout<<k<<endl;
     ostringstream obuffer1;
     obuffer1<<(1+(Double_t)k);
 	   numbern[k]=obuffer1.str();

     double_t t0 = 0;
     Int_t Nion =0;
     Int_t interactionID =0;

     /*
     filename2="/sps/hep/hadronth/JLL/Hadrontherapy/2014_08_27_Etude_Coincidence_fortuites_Hadrontherapy_Carbon_ions_305MeV_avec_resol_spatial/2014_08_27_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_CarbonIons_305MeV_Real_Spatial_Resolution_"+numbern[k]+".root.root";

     ///sps/hep/hadronth/JLL/Hadrontherapy/2014_01_23_Etude_Coincidences_fortuites_ICTR_PHE_Protons_Geant4_9_6/2014_01_23_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_ICTR_PHE_Proton_160MeV_"+numbern[k]+".root.root"
     //   /sps/hep/hadronth/JLL/Hadrontherapy/2013_12_16_Etude_Coincidences_fortuites_ICTR_PHE_Carbon/2013_12_18_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d1_20cm_ICTR_PHE_Carbon_1_
     //2013_12_19_Etude_Coincidences_fortuites_ICTR_PHE_Carbon_200MeV/2013_12_19_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_ICTR_PHE_Carbon_200MeV_
     //	/sps/hep/hadronth/JLL/Hadrontherapy/2013_01_04_Etude_Coincidences_fortuites_ICTR_PHE_Carbon_400MeV/2014_01_04_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_ICTR_PHE_Carbon_400MeV_"+numbern[k]+".root.root";


     //	/sps/hep/hadronth/JLL/Hadrontherapy/2014_01_23_Etude_Coincidences_fortuites_ICTR_PHE_Protons_Geant4_9_6/2014_01_23_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_ICTR_PHE_Proton_160MeV_"+numbern[k]+".root.root";

     //	/sps/hep/hadronth/JLL/Hadrontherapy/2013_11_18_Etude_Coincidences_fortuites_ICTR_PHE/2013_11_18_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_ICTR_PHE_1_"


     //"/sps/hep/hadronth/JLL/Hadrontherapy/2014_02_05_Etude_Coincidences_fortuites_ICTR_PHE_CarbonIons_305MeV_Geant4_9_6_camera_Plus_7_5cm_Yaxis/2014_02_05_PlanPhase_HadronT_SansBox_BGO_0_17_Si9cm_ENC_200_d2_40cm_d0_20cm_ICTR_PHE_CarbonIons_305MeV_Shift_7_5cm_Yaxis_"+numbern[k]+".root.root"

     fstream fin;
     fin.open(filename2.c_str(),ios::in);
     if( !fin.is_open() )
     {
		   cerr<<"error: file "<<filename2.c_str()<<" does not exist"<<endl;
     isfile=false;
     }
     else isfile=true;
     fin.close();

     if(isfile)
		   {
		   cout<<filename2.c_str()<<endl;
		   file = new TFile(filename2.c_str());
		   tree = (TTree*)(file)->Get("T_hitinfo");
		   tree->SetBranchAddress("B_hitinfo_interactionID",&interactionID);
		   tree->SetBranchAddress("B_hitinfo_t0",&t0);
		   tree->SetBranchAddress("B_hitinfo_Nion",&Nion);
		   nentries2 = tree->GetEntries();
		   nentries_total = nentries_total +nentries2 ;


		   for(Int_t i=0;i<nentries2;i++)
		   {
			  tree->GetEntry(i);
			  t0_ini.push_back(t0) ;
			  Nion_ini.push_back(Nion+100000*k);
     ID_ini.push_back(interactionID);
		   }

		   }
     }
     cout<< " nentries 2 " <<nentries_total<<endl;
     */
    //---------------------------------------------------------------

    //Lecture du fichier conteant les coinc trouvees et reconstruction

    Int_t test_ncoinc =0;

    file = new TFile(strname.c_str());

    mycoinc.hitLYSO.posprim=new Double_t[3];
    mycoinc.hitLYSO.posh=new Double_t[3];
    mycoinc.hitSi.posprim=new Double_t[3];
    mycoinc.hitSi.posh=new Double_t[3];
    if((file)->Get("T_2hit"))
    {
        tree=(TTree*)(file)->Get("T_2hit");
        tree->SetBranchAddress("B_hit_LYSO_interactionID",&mycoinc.hitLYSO.interactionID);
        tree->SetBranchAddress("B_hit_LYSO_posprim",mycoinc.hitLYSO.posprim);
        tree->SetBranchAddress("B_hit_LYSO_posh",mycoinc.hitLYSO.posh);
        tree->SetBranchAddress("B_hit_LYSO_edep",&mycoinc.hitLYSO.edep);
        tree->SetBranchAddress("B_hit_LYSO_E0",&mycoinc.hitLYSO.E0);
        tree->SetBranchAddress("B_hit_LYSO_tmin",&mycoinc.hitLYSO.tmin);
        tree->SetBranchAddress("B_hit_LYSO_tmax",&mycoinc.hitLYSO.tmax);
        tree->SetBranchAddress("B_hit_LYSO_t0",&mycoinc.hitLYSO.t0);
        tree->SetBranchAddress("B_hit_LYSO_tof",&mycoinc.hitLYSO.tof);
        tree->SetBranchAddress("B_hit_LYSO_Nion",&mycoinc.hitLYSO.Nion);
        tree->SetBranchAddress("B_hit_LYSO_nameprim",&mycoinc.hitLYSO.nameprim);
        tree->SetBranchAddress("B_hit_LYSO_procprim",&mycoinc.hitLYSO.procprim);
        tree->SetBranchAddress("B_hit_LYSO_gammatype",&mycoinc.hitLYSO.gammatype);
        tree->SetBranchAddress("B_hit_Si_interactionID",&mycoinc.hitSi.interactionID);
        tree->SetBranchAddress("B_hit_Si_posprim",mycoinc.hitSi.posprim);
        tree->SetBranchAddress("B_hit_Si_posh",mycoinc.hitSi.posh);
        tree->SetBranchAddress("B_hit_Si_edep",&mycoinc.hitSi.edep);
        tree->SetBranchAddress("B_hit_Si_E0",&mycoinc.hitSi.E0);
        tree->SetBranchAddress("B_hit_Si_tmin",&mycoinc.hitSi.tmin);
        tree->SetBranchAddress("B_hit_Si_tmax",&mycoinc.hitSi.tmax);
        tree->SetBranchAddress("B_hit_Si_t0",&mycoinc.hitSi.t0);
        tree->SetBranchAddress("B_hit_Si_tof",&mycoinc.hitSi.tof);
        tree->SetBranchAddress("B_hit_Si_Nion",&mycoinc.hitSi.Nion);
        tree->SetBranchAddress("B_hit_Si_nameprim",&mycoinc.hitSi.nameprim);
        tree->SetBranchAddress("B_hit_Si_procprim",&mycoinc.hitSi.procprim);
        tree->SetBranchAddress("B_hit_Si_gammatype",&mycoinc.hitSi.gammatype);
        tree->SetBranchAddress("B_ncoinc",&mycoinc.ncoinc);
        nentries=tree->GetEntries();
        cout<<"reconstruction, "<<nentries<<" entrees"<<endl;
        for(Int_t i=0;i<nentries;i++)//nentries
        {
            T0_diff=0;
            sum_CC=0;
            sum_edep=0;
            sum_edep_95=0;

            // if (mycoinc.ncoinc!=0)cout<<" ncoinc "<<mycoinc.ncoinc<<endl;

            //  cout<<"mycoinc.hitLYSO.interactionID "<<mycoinc.hitLYSO.interactionID<<endl;

            /*ym_TrueGamma = {-1000,-1000};
            ym_all = {-1000,-1000};
            ym_sameparticle= {-1000,-1000};
            ym_different= {-1000,-1000};
            ym_fortuit= {-1000,-1000};
            ym_background= {-1000,-1000};
            ym_gamma_prompt= {-1000,-1000};
            ym_gamma_prompt_scatt= {-1000,-1000};
            ym_gamma_other_value= {-1000,-1000};*/
            for(Int_t vc = 0; vc<2; vc++){
              ym_TrueGamma[vc] = -1000;
              ym_all[vc] = -1000;
              ym_sameparticle[vc] = -1000;
              ym_different[vc] = -1000;
              ym_fortuit[vc] = -1000;
              ym_background[vc] = -1000;
              ym_gamma_prompt[vc] = -1000;
              ym_gamma_prompt_scatt[vc] = -1000;
              ym_gamma_other_value[vc] = -1000;
            }

            //if ( mycoinc.hitLYSO.tof< 1000 || mycoinc.hitSi.tof< 1000) cout<<"i "<<i<<"Tof si "<<mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0<<"Tof lyso "<<mycoinc.hitLYSO.tof<<endl;

            //cout<<"posprim BGO x "<<mycoinc.hitLYSO.posprim[0]<<" posprim BGO y "<<mycoinc.hitLYSO.posprim[1]<<" posprim BGO z "<<mycoinc.hitLYSO.posprim[2]<<endl;
            // cout<<"posprim Si x "<<mycoinc.hitSi.posprim[0]<<" posprim Si y "<<mycoinc.hitSi.posprim[1]<<" posprim Si z "<<mycoinc.hitSi.posprim[2]<<endl<<endl;

            samepscoinc = false;
            gammacoinc = false ;
            trueeventcoinc = false;
            Sameparticule = false ;
            Diffparticule = false ;
            if(i%10000==0)cout<<i<<endl;
            tree->GetEntry(i);


            particletypesi=mycoinc.hitSi.nameprim;
            particletypelyso=mycoinc.hitLYSO.nameprim;
            creatorprocess=mycoinc.hitSi.procprim;

            edep1=mycoinc.hitSi.edep;
            edep2=mycoinc.hitLYSO.edep;
            x1=mycoinc.hitSi.posh[0];
            y1=mycoinc.hitSi.posh[1];
            z1=mycoinc.hitSi.posh[2];
            x2=mycoinc.hitLYSO.posh[0];
            y2=mycoinc.hitLYSO.posh[1];
            z2=mycoinc.hitLYSO.posh[2];

            tofSi=mycoinc.hitSi.tmax-mycoinc.hitSi.t0;
            tofLYSO=mycoinc.hitLYSO.tmax-mycoinc.hitLYSO.t0;
            gammatype=mycoinc.hitSi.gammatype;
            gammatype2=mycoinc.hitLYSO.gammatype;

            energy_initiale =mycoinc.hitSi.E0;



            //if (1 && IsTrue(mycoinc)!=0)

            for(Int_t ii=0;ii<64;ii++)ParticleTypeSi[ii]=mycoinc.hitSi.nameprim[ii];
            for(Int_t ii=0;ii<64;ii++)ParticleTypeLYSO[ii]=mycoinc.hitLYSO.nameprim[ii];
            for(Int_t ii=0;ii<64;ii++)gammaTypeSi[ii]=mycoinc.hitSi.gammatype[ii];
            for(Int_t ii=0;ii<64;ii++)gammaTypeLYSO[ii]=mycoinc.hitLYSO.gammatype[ii];
            for(Int_t ii=0;ii<64;ii++)ProcessSi[ii]=mycoinc.hitSi.procprim[ii];
            for(Int_t ii=0;ii<64;ii++)ProcessLYSO[ii]=mycoinc.hitLYSO.procprim[ii];

            //if (tofSi<3 && tofSi>1.5 || tofLYSO <3 && tofLYSO <1.5) cout<<" part "<<mycoinc.hitLYSO.nameprim<<" gamma type "<<mycoinc.hitSi.gammatype<<endl;



            //////JLL 22/01/2014
            /*for (Int_t g=0; g<nentries_total;g++)
             {
             if (Nion_ini[g]==mycoinc.hitLYSO.Nion||Nion_ini[g]==mycoinc.hitSi.Nion)
             {
             nombre_nion++;
             T0_diff= mycoinc.hitLYSO.t0-t0_ini[g];
             // cout <<"t0 ini "<<t0_ini[g]<<" mycoinc.hitLYSO.t0 "<<mycoinc.hitLYSO.t0<< "Diff "<<T0_diff<<endl;
             // cout<<" tmax avant : "<<mycoinc.hitLYSO.tmin<<endl;
             mycoinc.hitLYSO.tmax= mycoinc.hitLYSO.tmax-T0_diff;
             mycoinc.hitLYSO.tmin= mycoinc.hitLYSO.tmin-T0_diff;
             // cout<<" tmax BGO apr�s : "<<mycoinc.hitLYSO.tmin<<endl;
             mycoinc.hitSi.tmax= mycoinc.hitSi.tmax-T0_diff;
             mycoinc.hitSi.tmin= mycoinc.hitSi.tmin-T0_diff;
             cout<<" tmax BGO "<<mycoinc.hitLYSO.tmax<<" tmin BGO "<<mycoinc.hitLYSO.tmax<<" tmax Si "<<mycoinc.hitSi.tmax<<" tmin Si "<<mycoinc.hitSi.tmin<<endl;
             //	cout<<" TOF ini "<<mycoinc.hitSi.tof<<endl;
             break;
             }
             //cout <<"nombre_nion " <<nombre_nion<<endl;

             }*/


            if(cut==0)limit_high_cut_TOF=40;
            if(cut==1)limit_high_cut_TOF=6;//6

            sum_edep=edep1+edep2;
            sum_edep_95=energy_initiale*0.90;





          // if (edep1 >=0.05 && edep2>=0.1 && energy_initiale==sum_edep && mycoinc.hitLYSO.tof<limit_high_cut_TOF &&mycoinc.hitLYSO.tof>0 && mycoinc.ncoinc==1)

           //  if (edep1 >=0.05 && edep2>=0.1 && sum_edep>1  && mycoinc.hitLYSO.tof<limit_high_cut_TOF &&mycoinc.hitLYSO.tof>0 && mycoinc.ncoinc==1)
                if (edep1 >=0.05 && edep2>=0.1 && sum_edep>=sum_edep_95 && mycoinc.hitLYSO.tof<limit_high_cut_TOF &&mycoinc.hitLYSO.tof>0 && mycoinc.ncoinc>0)


                    // if (edep1 >=0 && edep2>=0 && sum_edep>1  && mycoinc.hitLYSO.tof<limit_high_cut_TOF &&mycoinc.hitLYSO.tof>0 && mycoinc.ncoinc==1) //mycoinc.hitLYSO.tmin >0.5 &&mycoinc.hitLYSO.tmin <4) //edep1<=3 &&edep2 <=10  )//&&  edep2>= 0.45
                //mycoinc.hitLYSO.tmax >0 &&mycoinc.hitLYSO.tmax <3 &&
                //if (edep2 <0.1 || edep2 <0.05)cout<<" edep 1 :"<<edep1<<" edep2 "<<edep2<<endl;
                //mycoinc.hitLYSO.tmin+mycoinc.hitLYSO.t0<10 &&mycoinc.hitLYSO.tmin+mycoinc.hitLYSO.t0>0

            {
                num_coinc++;
                particule_init =mycoinc.hitSi.nameprim;
                type_particule_init=mycoinc.hitSi.gammatype;
                if(particule_init=="gamma")num_gamma_init ++;
                if(particule_init== "gamma" && type_particule_init == "PromptGamma") num_gamma_prompt_ini++;
                if(particule_init == "gamma" && type_particule_init == "ScattPrompt")num_gamma_scatt_init++;
                if(particule_init== "gamma" && type_particule_init == "OtherGamma")num_gamma_other++;

                if(particule_init=="gamma")edep_histo->Fill(edep1,edep2);
                if(particule_init== "gamma" && type_particule_init == "PromptGamma")edep_histo_true->Fill(edep1,edep2);
                if(particule_init == "gamma" && type_particule_init == "ScattPrompt")
                {
                    edep_histo_false ->Fill(edep1,edep2);
                    spectre_energy_gamma->Fill(mycoinc.hitSi.E0);

                }
                if(particule_init== "gamma" && type_particule_init == "OtherGamma")edep_histo_true_gamma->Fill(edep1,edep2);




                ncoinc++;
                //if(i>48000)cout<<"i avant "<<i<<endl;
                //on regarde quel type de coinc a eu lieu


                // cout<<" Si "<<mycoinc.hitSi.Nion<<" bgo "<<mycoinc.hitLYSO.Nion<<endl;
                if( IsTrue(mycoinc) || mycoinc.hitSi.Nion==mycoinc.hitLYSO.Nion)//
                {

                    //if(i >10000) cout<< "i "<<i<<endl;
                    /*cout<<"BGO : t0 "<<mycoinc.hitLYSO.t0<<" tmax "<<mycoinc.hitLYSO.tmax<<" tmin "<<mycoinc.hitLYSO.tmin<<endl;
                     cout<<"si : t0 "<<mycoinc.hitSi.t0<<" time "<<mycoinc.hitSi.tmin<<" "<<mycoinc.hitSi.tmax<<endl;
                     cout<< " diff BGO "<<mycoinc.hitLYSO.tmax-mycoinc.hitLYSO.t0<<" diff tmin "<<mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0<<endl;
                     cout<< " diff Si "<<mycoinc.hitSi.tmax-mycoinc.hitSi.t0<<" diff tmin "<<mycoinc.hitSi.tmin-mycoinc.hitSi.t0<<endl<<endl;
                     */
                    //cout<< "part "<<myhit.nameprim<<endl<<endl;*/


                    //neventstore++;
                    fortuitparticule = false;

                    n_vrai++;
                    samepscoinc=true;
                    nsamepscoinc++;
                    string typegammaSi =mycoinc.hitSi.gammatype;
                    string typegammaLYSO =mycoinc.hitLYSO.gammatype;


                    if(particletypesi=="gamma" && particletypelyso=="gamma" &&typegammaSi ==typegammaLYSO&&mycoinc.hitSi.E0==mycoinc.hitLYSO.E0)
                    {


                        //new 18/01/2015

                       // num_gamma_init++;

                        gammacoinc=true;
                        ngammacoinc++;

                        size_t found;
                        found=creatorprocess.find("nelastic");
                        if(found!=string::npos)creatorprocess="Inelastic";

                        neventstore++;



                        string gammatype_Si =mycoinc.hitSi.gammatype;
                        string gammatype_BGO =mycoinc.hitLYSO.gammatype;

                        string gammaprocess_Si =mycoinc.hitSi.procprim;
                        string gammaprocess_BGO =mycoinc.hitLYSO.procprim;
                        if(gammatype_Si=="PromptGamma" && gammatype_BGO=="PromptGamma")

                            //if(gammaprocess_Si!="annihil" && gammatype_Si!="PromptGamma")//gammaprocess_BGO=="annihil")
                            //if(mycoinc.hitSi.interactionID==0 && (mycoinc.hitLYSO.interactionID==0 ||mycoinc.hitLYSO.interactionID==1))
                        {
                            trueeventcoinc=true;
                            ntrueeventcoinc++;
                            //num_gamma_prompt_ini++;

                        }


                        if(gammatype=="OtherGamma")
                        {
                            //num_gamma_other++;
                        }

                        if(gammatype=="ScattPrompt")
                        {
                            //num_gamma_scatt_init++;
                        }


                        /*	ofstream
                         fichier("2014_02_25_Planphase_305MeVCarbonIons_BGO_d0_20cm_d2_40cm_Bunch_8Ions_Shift12cm_Gammas.txt",ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
                         if (fichier)
                         {
                         if (edep1 >0.0001 && edep2 >0.0001)
                         {
                         //if (abs(x)>0.00001 && abs(B_gauss_pos1[0])>0.00001 && abs(B_gauss_pos1[0])>0.00001)
                         //{
                         //if(abs(B_gauss_pos2[0])>0.00001 && abs(B_gauss_pos2[0])>0.00001 && abs(B_gauss_pos2[0])>0.00001)
                         //{
                         //cout <<"pos1[0] "<<x1<<endl;
                         fichier<<"---Nouvel evenement---"<<endl;
                         fichier<<"position Si "<<x1<<" "<<y1<<" "<<z1<<endl;
                         fichier<<"energie Si " <<edep1<<endl;
                         fichier<<"position Absorbeur "<<x2<<" "<<y2<<" "<<z2<<endl;
                         fichier<<"energie Absorbeur " <<edep2<<endl;
                         //}
                         //}
                         }
                         }
                         fichier.close();*/




                    }


                    if(particletypesi==particletypelyso &&particletypesi!="gamma" &&particletypelyso!="gamma" &&mycoinc.hitSi.E0==mycoinc.hitLYSO.E0)
                    {

                        Sameparticule_number ++; // meme particules autre que gamma
                        Sameparticule = true ;
                        tof_histo_sameparticle->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitSi.tmin-mycoinc.hitSi.t0);

                    }

                    if(particletypesi!=particletypelyso || particletypesi==particletypelyso && mycoinc.hitSi.E0!=mycoinc.hitLYSO.E0 || particletypesi=="gamma" &&particletypelyso=="gamma" && typegammaSi!=typegammaLYSO)
                    {

                        Diffparticule_number ++; // deux particules differentes
                        Diffparticule = true ;
                    }

                    if( Diffparticule==false && Sameparticule==false && gammacoinc==false && fortuitparticule==false)
                    {
                        cout<<"Si part type "<<particletypesi<<" Si nion "<<mycoinc.hitSi.Nion<<" E0 Si "<<mycoinc.hitSi.E0<<" BGO part type "<<particletypelyso<<" BGO nion "<<mycoinc.hitLYSO.Nion<<" Eo BGO " <<mycoinc.hitLYSO.E0<< endl;
                        //verif_otherevents ++;
                    }
                }
                else
                {
                    samepscoinc=false;
                    gammacoinc=false;
                    trueeventcoinc=false;
                    fortuitparticule = true;
                    //if(mycoinc.hitSi.Nion!=mycoinc.hitLYSO.Nion) cout<<"Nion differents "<<mycoinc.hitSi.Nion<<" "<<mycoinc.hitLYSO.Nion<<endl;
                    n_fortuits++;
                    tof_histo_false->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitSi.tmin-mycoinc.hitSi.t0);
                    //nrandom_ym++;
                }




                //cout<<"ntrue events :"<<ntrueeventcoinc<<endl;
                //on enregistre les infos relative aux differents types de coinc dans differents histos

                //cout <<" posprim : "<<mycoinc.hitSi.posprim[1]<<endl;

                xe_histo->Fill(mycoinc.hitSi.posprim[0]);
                ye_histo->Fill(mycoinc.hitSi.posprim[1]);
                ze_histo->Fill(mycoinc.hitSi.posprim[2]);
                if(samepscoinc)ye_histo_true->Fill(mycoinc.hitSi.posprim[1]);
                if(fortuitparticule)ye_histo_false->Fill(mycoinc.hitSi.posprim[1]);

                type_histo->Fill(ParticleTypeSi,ParticleTypeLYSO,1);
                if(gammacoinc)type_histo_true->Fill(ParticleTypeSi,ParticleTypeLYSO,1);
                if(fortuitparticule)
                {
                    type_histo_false->Fill(ParticleTypeSi,ParticleTypeLYSO,1);
                    if(particletypesi=="e-")
                    {
                        ne_moins=ne_moins +mycoinc.hitSi.E0;
                        nemoins_sum++;
                        spectre_energy_electrons->Fill(mycoinc.hitSi.E0);
                    }
                    if(particletypesi=="proton")
                    {
                        nprotons=nprotons +mycoinc.hitSi.E0;
                        spectre_energy_protons->Fill(mycoinc.hitSi.E0);
                        nprotons_sum++;
                    }
                    if(particletypesi=="e+")
                    {
                        ne_plus=ne_plus +mycoinc.hitSi.E0;
                        neplus_sum++;
                    }
                    if(particletypesi=="neutron")
                    {
                        n_neutrons=n_neutrons +mycoinc.hitSi.E0;
                        n_neutrons_sum++;
                        spectre_energy_neutrons->Fill(mycoinc.hitSi.E0);
                    }
                    //cout<<"edep neutron "<<mycoinc.hitSi.E0<<endl;
                    //if(particletypesi=="e-")cout<<"edep e- "<<mycoinc.hitSi.E0<<endl;
                }

                ;

                gammatype_histo->Fill(gammaTypeSi,gammaTypeLYSO,1);
                if(gammacoinc)gammatype_histo_true->Fill(gammaTypeSi,gammaTypeLYSO,1);
                //if(fortuitparticule)gammatype_histo_false->Fill(gammaTypeSi,gammaTypeLYSO,1);

                process_histo->Fill(ProcessSi,ProcessLYSO,1);
                if(gammacoinc)process_histo_true->Fill(ProcessSi,ProcessLYSO,1);
                if(fortuitparticule)process_histo_false->Fill(ProcessSi,ProcessLYSO,1);

                pp_histo->Fill(ParticleTypeSi,ProcessSi,1);


                nion1_histo->Fill(mycoinc.hitLYSO.Nion);
                if(mycoinc.hitLYSO.Nion==7801 || (mycoinc.hitLYSO.Nion>5331000 && mycoinc.hitLYSO.Nion<5332000 ))
                {
                    Print(mycoinc.hitSi);
                    Print(mycoinc.hitLYSO);
                }
                nion_histo->Fill(mycoinc.hitSi.Nion,mycoinc.hitLYSO.Nion);
                if(samepscoinc)nion_histo_true->Fill(mycoinc.hitSi.Nion,mycoinc.hitLYSO.Nion);
                if(fortuitparticule)nion_histo_false->Fill(mycoinc.hitSi.Nion,mycoinc.hitLYSO.Nion);

                //edep_histo->Fill(edep1,edep2);
                //edep_histo_true->Fill(edep1, edep2);//trueeventcoinc //samepscoinc && particletypesi=="neutron"
              //  if(samepscoinc)edep_histo_true->Fill(mycoinc.hitLYSO.tmin, edep2);
               // if(fortuitparticule)edep_histo_false->Fill(mycoinc.hitLYSO.tmin,edep2);
               // if(gammacoinc)edep_histo_true_gamma->Fill(mycoinc.hitLYSO.tmin,edep2);

                n_histo->Fill(mycoinc.hitSi.posh[2]);
                if(samepscoinc)n_histo_true->Fill(mycoinc.hitSi.posh[2]);
                if(fortuitparticule)n_histo_false->Fill(mycoinc.hitSi.posh[2]);

                //if(samepscoinc)nvraie_ym++;
                //if(!samepscoinc)nrandom_ym++;
                //if (gammacoinc)ngammatrue_ym++;


                //Reconstruction line cone


                y=ReconstructionCC(x1,y1,z1,x2,y2,z2,edep1,edep2,0,0);





                //if( Sameparticule ==true || gammacoinc== true) y=ReconstructionCC(x1,y1,z1,x2,y2,z2,edep1,edep2,mycoinc.hitSi.E0,0);


                gamma_out=false;

                //newstrname_MLEM="/sps/hep/hadronth/JLL/Macro_Root/Coincidences_reelles_fortuites/code/verification_article_oct_2016/2017_04_30_1Proton_160MeV_HIT_9_4nsBunch_2nsSizeBunch_40nswindowscoinc_ResolSi15ns_ResolBGO3nsHigh_Stat_1e10events_CutEnergies_NoCutSumEnergies_MultiCoinc_Timing_EdepEgalEinit90_Check_MLEM_Datas";//+num_file;
                newstrname_MLEM=outpath + "test_MLEM";
		            newstrname_MLEM+=".txt";

                //cout<<newstrname_MLEM.c_str()<<endl;
                ofstream fichier;
		            fichier.open(newstrname_MLEM.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
                if (fichier)
                {
                    if (edep1 >0.0001 && edep2 >0.0001)
                    {
                        //if (abs(x)>0.00001 && abs(B_gauss_pos1[0])>0.00001 && abs(B_gauss_pos1[0])>0.00001)
                        //{
                        //if(abs(B_gauss_pos2[0])>0.00001 && abs(B_gauss_pos2[0])>0.00001 && abs(B_gauss_pos2[0])>0.00001)
                        //{
                        //	cout <<"pos1[0] "<<B_gauss_pos1[0]<<endl;
                        fichier<<"---Nouvel evenement---"<<endl;
                        fichier<<"position Si "<<x1<<" "<<y1<<" "<<z1<<endl;
                        fichier<<"energie Si " <<edep1<<endl;
                        fichier<<"position Absorbeur "<<x2<<" "<<y2<<" "<<z2<<endl;
                        fichier<<"energie Absorbeur " <<edep2<<endl;
                        //}
                        //}
                    }
                }else{cout<<"Problem opening file FICHIER for MLEM"<<endl;}
                fichier.close();






                /*if (gammacoinc)
                 {
                 y=ReconstructionCC(x1,y1,z1,x2,y2,z2,edep1,edep2);
                 //cout <<"x1 "<<x1<<" x2 "<<x2<<" y1 "<<y1<<" y2 "<<y2<<" z1 "<<z1<<" z2 "<<z2<<"edep1 "<<edep1<<" edep2 "<<edep2<<endl;
                 count_gamma_out++;
                 if (y!=NULL)
                 {
                 gamma_out =true;

                 //ym=IsCloser(y[0],y[1],mycoinc.hitSi.posprim[1]);

                 //if(ym >15 )
                 //{
                 type_histo->Fill(ParticleTypeSi,ParticleTypeLYSO,1);
                 //if(gammacoinc)type_histo_true->Fill(ParticleTypeSi,ParticleTypeLYSO,1);
                 //if(!samepscoinc)type_histo_false->Fill(ParticleTypeSi,ParticleTypeLYSO,1);

                 gammatype_histo->Fill(gammaTypeSi,gammaTypeLYSO,1);
                 if(gammacoinc)gammatype_histo_true->Fill(gammaTypeSi,gammaTypeLYSO,1);
                 if(!samepscoinc)gammatype_histo_false->Fill(gammaTypeSi,gammaTypeLYSO,1);

                 process_histo->Fill(ProcessSi,ProcessLYSO,1);
                 if(gammacoinc)process_histo_true->Fill(ProcessSi,ProcessLYSO,1);
                 if(!samepscoinc)process_histo_false->Fill(ProcessSi,ProcessLYSO,1);

                 string gammatype_Si =mycoinc.hitSi.gammatype;
                 string gammatype_BGO =mycoinc.hitLYSO.gammatype;

                 if(gammatype_Si=="PromptGamma" && gammatype_BGO=="PromptGamma")
                 {
                 //cout<<"edep1 = "<<edep1<<" edep2 = "<<edep2 <<" somme = "<<edep1+edep2<<" E0BGo "<<mycoinc.hitLYSO.E0<<" E0Si "<<mycoinc.hitSi.E0<<endl;
                 if (edep1+edep2 >=(mycoinc.hitLYSO.E0-mycoinc.hitLYSO.E0*0.1) && edep1+edep2 <=(mycoinc.hitLYSO.E0+mycoinc.hitLYSO.E0*0.1)) count_sumAccurate ++;
                 //if (edep1+edep2 ==mycoinc.hitLYSO.E0) count_sumAccurate ++;
                 }
                 //}

                 }
                 }*/


                if(y!=NULL )//&& tofLYSO<40)//modif JLL 40 au lieu de

                    //10 //y!=NULL && tofSi<10 && tofLYSO<10
                {


                    if(y[0]==y[1])
                    {

                        /*ym=y[0];
                         ym_histo->Fill(ym);
                         if(samepscoinc)
                         {
                         ym_histo_true->Fill(ym);
                         nvraie_ym++;
                         }
                         //if(samepscoinc && !gammacoinc)ym_histo_true->Fill(ym);
                         if(!samepscoinc)
                         {
                         ym_histo_false->Fill(ym);
                         nrandom_ym++;
                         }
                         if (gammacoinc)
                         {
                         ym_histo_true_same_gamma->Fill(ym);
                         ngammatrue_ym++;
                         }
                         if(gammatype=="PromptGamma" && samepscoinc)ym_histo_prompt->Fill(ym);
                         if(TMath::Abs(ym)<300)rectree->Fill();*/


                    }
                    else
                    {
                        ym=IsCloser(y[0],y[1],mycoinc.hitLYSO.posprim[1]);

                        //cout<<"ym "<<ym<<endl;
                        //if(i>70000) cout<<"test i "<<i<<endl;
                        //cout<<"y0 "<<y[0]<<"y1 "<<y[1]<<endl;
                        //if(ym >-150 && ym <150 )

                        //Add le 04/02/2015
                      //  if (ym ==y[0])ym_histo_autre_solution->Fill(y[1]);
                       // if (ym ==y[1])ym_histo_autre_solution->Fill(y[0]);

                        //Inversion solutions : regarde reconstruction de la solution la plus �loign�e
                       // if (ym ==y[0])ym=y[1];
                       // if (ym ==y[1])ym=y[0];


                        //R�flexion sur la pond�ration avec l'efficacit�

                         sum_CC = edep1+edep2;

                        weight[0]=1;
                        weight[1]=1;
                        difference_weight=0;
                        up[0] =0;
                        up[1] =0;

                        for( int k=0;k<2;k++)
                        {
                            minimum_value =1000;


                         for (int p =0; p<30;p++)
                         {
                            difference_position = TMath::Abs((y[k]-50) - position_value[p]);
                            if( difference_position <= minimum_value)
                            {
                                minimum_value =difference_position;
                                up[k] = p;
                            }
                            /* difference_position = TMath::Abs((y[k]-50) - position_value[p]);
                             if( difference_position <= minimum_value)
                             {
                                 minimum_value =difference_position;
                                 up = p;
                             }*/
                         }
                        }
                           // cout<<"y[15] "<<position_value[15]<<endl;
                          //  cout<<" i " <<i<<" minimum value "<<minimum_value<< " difference position "<<difference_position<< " up "<<up<< " y["<<k<<"] "<<y[k]<<endl;
                         if (sum_CC >=1 && sum_CC <1.5)
                         {
                            // weight[k] =  value_efficacite_1MeV[up] /value_efficacite_1MeV[15];
                            //weight[k] =  1-(value_efficacite_1MeV[up] /value_efficacite_1MeV[15])+1;// ref
                            // weight[k] = 1/value_efficacite_1MeV[up];
                            /* if (value_efficacite_1MeV[up[0]]> value_efficacite_1MeV[up[1]])
                             {
                                 weight[0] =  1;
                                 weight[1] = value_efficacite_1MeV[up[1]]/value_efficacite_1MeV[up[0]] ;
                             }
                             if (value_efficacite_1MeV[up[1]]> value_efficacite_1MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] = value_efficacite_1MeV[up[0]]/value_efficacite_1MeV[up[1]] ;
                             }
                             if (value_efficacite_1MeV[up[1]]== value_efficacite_1MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] =  1;
                             }
                            // cout<<"1MeV  k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_1MeV "<<value_efficacite_1MeV[up]<< " value_efficacite_1MeV[15] "<<value_efficacite_1MeV[15]<<endl;
                             */
                             weight[0] =value_efficacite_1MeV[up[0]]/(value_efficacite_1MeV[up[0]]+value_efficacite_1MeV[up[1]]) ;
                             weight[1] = value_efficacite_1MeV[up[1]]/(value_efficacite_1MeV[up[0]]+value_efficacite_1MeV[up[1]]) ;

                         }

                         if (sum_CC >=1.5 && sum_CC <2.5)
                         {
                            // weight[k] =  value_efficacite_2MeV[up] /value_efficacite_2MeV[15];
                            // weight[k] =  1-(value_efficacite_2MeV[up] /value_efficacite_2MeV[15])+1;
                            // weight[k] = 1/value_efficacite_2MeV[up];
                             /*if (value_efficacite_2MeV[up[0]]> value_efficacite_2MeV[up[1]])
                             {
                                 weight[0] =  1;
                                 weight[1] = value_efficacite_2MeV[up[1]]/value_efficacite_2MeV[up[0]] ;
                             }
                             if (value_efficacite_2MeV[up[1]]> value_efficacite_2MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] = value_efficacite_2MeV[up[0]]/value_efficacite_2MeV[up[1]] ;
                             }
                             if (value_efficacite_2MeV[up[1]]== value_efficacite_2MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] =  1;
                             }*/
                            // cout<<" 2MeV k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_2MeV "<<value_efficacite_2MeV[up]<< " value_efficacite_2MeV[15] "<<value_efficacite_2MeV[15]<<endl;

                             weight[0] =value_efficacite_2MeV[up[0]]/(value_efficacite_2MeV[up[0]]+value_efficacite_2MeV[up[1]]) ;
                             weight[1] = value_efficacite_2MeV[up[1]]/(value_efficacite_2MeV[up[0]]+value_efficacite_2MeV[up[1]]) ;
                         }

                         if (sum_CC >=2.5 && sum_CC <3.5)
                         {
                             //weight[k] =  value_efficacite_3MeV[up] /value_efficacite_3MeV[15];
                            // weight[k] =  1-(value_efficacite_3MeV[up] /value_efficacite_3MeV[15])+1;
                            // weight[k] = 1/value_efficacite_3MeV[up];
                             /*if (value_efficacite_3MeV[up[0]]> value_efficacite_3MeV[up[1]])
                             {
                                 weight[0] =  1;
                                 weight[1] = value_efficacite_3MeV[up[1]]/value_efficacite_3MeV[up[0]] ;
                             }
                             if (value_efficacite_3MeV[up[1]]> value_efficacite_3MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] = value_efficacite_3MeV[up[0]]/value_efficacite_3MeV[up[1]] ;
                             }
                             if (value_efficacite_3MeV[up[1]]== value_efficacite_3MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] =  1;
                             }*/
                           // cout<<" 3MeV k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_3MeV "<<value_efficacite_3MeV[up]<< " value_efficacite_3MeV[15] "<<value_efficacite_3MeV[15]<<endl;

                             weight[0] =value_efficacite_3MeV[up[0]]/(value_efficacite_3MeV[up[0]]+value_efficacite_3MeV[up[1]]) ;
                             weight[1] = value_efficacite_3MeV[up[1]]/(value_efficacite_3MeV[up[0]]+value_efficacite_3MeV[up[1]]) ;
                         }
                         if (sum_CC >=3.5 && sum_CC <4.5)
                         {
                            // weight[k]=value_efficacite_4MeV[up] /value_efficacite_4MeV[15];
                            // weight[k] =  1-(value_efficacite_4MeV[up] /value_efficacite_4MeV[15])+1;
                             //weight[k] = 1/value_efficacite_4MeV[up];
                             /*if (value_efficacite_4MeV[up[0]]> value_efficacite_4MeV[up[1]])
                             {
                                 weight[0] =  1;
                                 weight[1] = value_efficacite_4MeV[up[1]]/value_efficacite_4MeV[up[0]] ;
                             }
                             if (value_efficacite_4MeV[up[1]]> value_efficacite_4MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] = value_efficacite_4MeV[up[0]]/value_efficacite_4MeV[up[1]] ;
                             }
                             if (value_efficacite_4MeV[up[1]]== value_efficacite_4MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] =  1;
                             }*/
                            // cout<<" 4MeV k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_4MeV "<<value_efficacite_4MeV[up]<< " value_efficacite_4MeV[15] "<<value_efficacite_4MeV[15]<<endl;

                             weight[0] =value_efficacite_4MeV[up[0]]/(value_efficacite_4MeV[up[0]]+value_efficacite_4MeV[up[1]]) ;
                             weight[1] = value_efficacite_4MeV[up[1]]/(value_efficacite_4MeV[up[0]]+value_efficacite_4MeV[up[1]]) ;
                         }
                        if (sum_CC >=4.5 && sum_CC <5.5)
                            {
                                //weight[k]=value_efficacite_5MeV[up] /value_efficacite_5MeV[15];
                                //weight[k] =  1-(value_efficacite_5MeV[up] /value_efficacite_5MeV[15])+1;
                               // weight[k] = 1/value_efficacite_5MeV[up];
                                /*if (value_efficacite_5MeV[up[0]]> value_efficacite_5MeV[up[1]])
                                {
                                    weight[0] =  1;
                                    weight[1] = value_efficacite_5MeV[up[1]]/value_efficacite_5MeV[up[0]] ;
                                }
                                if (value_efficacite_5MeV[up[1]]> value_efficacite_5MeV[up[0]])
                                {
                                    weight[1] =  1;
                                    weight[0] = value_efficacite_5MeV[up[0]]/value_efficacite_5MeV[up[1]] ;
                                }
                                if (value_efficacite_5MeV[up[1]]== value_efficacite_5MeV[up[0]])
                                {
                                    weight[1] =  1;
                                    weight[0] =  1;
                                }*/
                              //  cout<<" 5MeV k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_5MeV "<<value_efficacite_5MeV[up]<< " value_efficacite_5MeV[15] "<<value_efficacite_5MeV[15]<<endl;
                                weight[0] =value_efficacite_5MeV[up[0]]/(value_efficacite_5MeV[up[0]]+value_efficacite_5MeV[up[1]]) ;
                                weight[1] = value_efficacite_5MeV[up[1]]/(value_efficacite_5MeV[up[0]]+value_efficacite_5MeV[up[1]]) ;

                            }
                         if (sum_CC >=5.5 && sum_CC <20)
                         {
                            // weight [k] =  value_efficacite_6MeV[up] /value_efficacite_6MeV[15];
                            // weight[k] =  1-(value_efficacite_6MeV[up] /value_efficacite_6MeV[15])+1;
                            // weight[k] = 1/value_efficacite_6MeV[up];
                             /*if (value_efficacite_6MeV[up[0]]> value_efficacite_6MeV[up[1]])
                             {
                                 weight[0] =  1;
                                 weight[1] = value_efficacite_6MeV[up[1]]/value_efficacite_6MeV[up[0]] ;
                             }
                             if (value_efficacite_6MeV[up[1]]> value_efficacite_6MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] = value_efficacite_6MeV[up[0]]/value_efficacite_6MeV[up[1]] ;
                             }
                             if (value_efficacite_6MeV[up[1]]== value_efficacite_6MeV[up[0]])
                             {
                                 weight[1] =  1;
                                 weight[0] =  1;
                             }*/
                          //   cout<<" 6MeV k "<<k<<" sumCC " <<sum_CC<< " weight " <<weight[k]<< " value_efficacite_6MeV "<<value_efficacite_3MeV[up]<< " value_efficacite_6MeV[15] "<<value_efficacite_6MeV[15]<<endl;

                             weight[0] =value_efficacite_6MeV[up[0]]/(value_efficacite_6MeV[up[0]]+value_efficacite_6MeV[up[1]]) ;
                             weight[1] = value_efficacite_6MeV[up[1]]/(value_efficacite_6MeV[up[0]]+value_efficacite_6MeV[up[1]]) ;
                         }



                        //}
                         //-> rajouter weight � l'histogram : -> Fill (y[0],weight);

                        //cout<<" !!!!! weight 0 "<<weight[0]<< " weight 1 "<<weight[1]<<" diff yo-y1 "<<TMath::Abs(y[0]-y[1])<< endl;

                        //if(y[0] >0 && y[0] <60 || y[1] >0 && y[1] <60)
                        //{
                           difference_weight= TMath::Abs(weight[0]-weight[1]);
			 histo_difference_weight->Fill(difference_weight);
			// histo_difference_weight->Fill(weight[1]);
                        //}
                        //cout<<endl<<" !!!!! weight 0 "<<weight[0]<< " weight 1 "<<weight[1]<<" diff yo-y1 "<<TMath::Abs(y[0]-y[1])<< endl;

                        weight[0]=1;
                        weight[1]=1;



                        if (sum_CC <1) test_Sum++;

                        //cout <<"sum "<<sum_CC<<endl;
                        if(y[0] >-300 && y[0] <300 || y[1] >-300 && y[1] <300)
                        //if(ym >-300 && ym <300 )
                        {




                            ym_vertex->Fill(mycoinc.hitLYSO.posprim[1]);
                            nombre_coinc_reconstruite ++;

                            tof_histo->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0);
                            if(!gammacoinc)tof_histo_true->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0);
                            if(gammacoinc)	tof_histo_samegamma->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0);

                            //if(	particletypelyso!="gamma" && )cout <<"t diff " <<mycoinc.hitSi.tmin-mycoinc.hitSi.t0<<" Particle BGO : "<<ParticleTypeLYSO<< " Particule Si : "<<ParticleTypeSi<<endl;


                            //cout<<"y1 : " <<y[0]<<" y2: "<<y[1]<<endl;

                            //cout <<"mycoinc.hitSi.posprim[1] : "<<mycoinc.hitSi.posprim[1]<<endl;
                            //ym=y[0];
                            // ym_histo->Fill(ym-25);pour PTCOG
                            //ym_histo->Fill(y[0]);
                            //ym_histo->Fill(ym);
                            //if(trueeventcoinc)ym_histo_true->Fill(y);



                            ym_all[0] =y[0];
                            ym_all[1]= y[1];
                            /*Double_t ym_sameparticle=0;
                             Double_t ym_different=0;
                             Double_t =0;
                             ym_background
                             ym_TrueGamma*/

                            if(samepscoinc)
                            {

                                ym_histo_true->Fill(y[0],weight[0]);
                                ym_histo_true->Fill(y[1],weight[1]);
                               // ym_histo_true->Fill(ym);
                                nvraie_ym++;
                                //if(!gammacoinc && samepscoinc)ym_histo->Fill(y[0]-25);
                                //if(!gammacoinc && samepscoinc)ym_histo->Fill(y[1]-25);

                                //	nombre_coinc_reconstruite ++;
                                Histo_difference_y0[nombre_coinc_reconstruite] =  (y[0]-mycoinc.hitLYSO.posprim[1]);//TMath::Abs(y[0]-mycoinc.hitLYSO.posprim[1]);//
                                Histo_difference_y1[nombre_coinc_reconstruite] =(y[1]-mycoinc.hitLYSO.posprim[1]); // TMath::Abs(y[1]-mycoinc.hitLYSO.posprim[1]);//
                                Histo_difference_y0y1[nombre_coinc_reconstruite] = (y[1]-y[2]);//TMath::Abs(y[0]-y[1]);//

                                histo_yVertex_y1y0 ->Fill(Histo_difference_y0y1[nombre_coinc_reconstruite]);
                                histo_yVertex_y1 ->Fill(Histo_difference_y1[nombre_coinc_reconstruite]);
                                histo_yVertex_y0->Fill(Histo_difference_y0[nombre_coinc_reconstruite]);

                            }
                            if(!gammacoinc && samepscoinc)
                            {
                                //ym_histo->Fill(ym-25);
                                //ym_histo->Fill(y[0]-25);
                                //ym_histo->Fill(y[1]-25);
                            }

                            if(fortuitparticule)
                            {

                                ym_histo_false->Fill(y[0],weight[0]);
                                ym_histo_false->Fill(y[1],weight[1]);
                               // ym_histo_false->Fill(ym);
                                nrandom_ym++;
                                ym_histo->Fill(y[0],weight[0]);
                                ym_histo->Fill(y[1],weight[1]);
                                //ym_histo->Fill(ym);
                                tof_histo_true->Fill(mycoinc.hitLYSO.tof);//mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0);

                                ym_fortuit[0]=y[0];
                                ym_fortuit[1]=y[1];
                                ym_background[0]=y[0];
                                ym_background[1]=y[1];
                                reconstree->Fill();

                            }

                            /*if (samepscoinc)// ajout 30/05/2014 pour PTCOG
                             {
                             if(!gammacoinc)ym_histo_false->Fill(ym-25);
                             }*/

                            if (gammacoinc)
                            {
                                ym_histo_true_same_gamma->Fill(y[0],weight[0]);
                                ym_histo_true_same_gamma->Fill(y[1],weight[1]);
                                // ym_histo_true_same_gamma->Fill(ym);
                                ngammatrue_ym++;
                                num_gamma_recons++;


                                ym_TrueGamma[0] =y[0];
                                ym_TrueGamma[1] =y[1];


                                difference_y0y1 = TMath::Abs(y[0]-y[1]);
                               // cout <<" diff "<<difference_y0y1<<endl;
                                histo_distancey0y1->Fill(difference_y0y1);

                                if( mycoinc.hitLYSO.posprim[1]>-100 && mycoinc.hitLYSO.posprim[1]<0)
                                {
                                    difference_y0y1_entertarget =TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_enterTarget->Fill(difference_y0y1_entertarget);
                                }

                                if( mycoinc.hitLYSO.posprim[1]>0 && mycoinc.hitLYSO.posprim[1]<65)
                                {
                                    difference_y0y1_fallofftarget  =TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_FallOffTarget->Fill(difference_y0y1_fallofftarget);
                                }

                                if( mycoinc.hitLYSO.posprim[1]>65 && mycoinc.hitLYSO.posprim[1]<100)
                                {
                                    difference_y0y1_outtarget =TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_outTarget->Fill(difference_y0y1_outtarget);
                                }

                                if (mycoinc.hitSi.E0<0)
                                {
                                    difference_y0y1_cut0_1MeV= TMath::Abs(y[0]-y[1]);
                                    cout<<"test "<<difference_y0y1_cut0_1MeV<<endl;
                                    histo_distancey0y1_cut0_1MeV->Fill(difference_y0y1_cut0_1MeV);
                                }

                                if ( mycoinc.hitSi.E0>1 && mycoinc.hitSi.E0<=3)
                                {
                                    difference_y0y1_cut1_3MeV = TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_cut1MeV_3MeV->Fill(difference_y0y1_cut1_3MeV);
                                }


                                if ( mycoinc.hitSi.E0>3)
                                {
                                    difference_y0y1_cut3MeV_sup = TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_cut3MeV_sup->Fill(difference_y0y1_cut3MeV_sup);
                                }


                                if ( sum_CC>1 && sum_CC<=3)
                                {
                                    difference_y0y1_cut1_3MeV_edep = TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_cut1MeV_3MeV_edep->Fill(difference_y0y1_cut1_3MeV_edep);
                                }


                                if ( sum_CC>3)
                                {
                                    difference_y0y1_cut3MeV_sup_edep = TMath::Abs(y[0]-y[1]);
                                    histo_distancey0y1_cut3MeV_sup_edep->Fill(difference_y0y1_cut3MeV_sup_edep);
                                }

                                if(gammatype=="PromptGamma")
                                {
                                    ym_histo_prompt->Fill(y[0],weight[0]);
                                    ym_histo_prompt->Fill(y[1],weight[1]);
                                    //ym_histo_prompt->Fill(ym);
                                    ym_gamma_prompt[0]=y[0];
                                    ym_gamma_prompt[1]=y[1];
                                    num_gamma_prompt_recons++;

                                }
                                if(gammatype=="OtherGamma")
                                {
                                    ym_gamma_other->Fill(y[0],weight[0]);
                                    ym_gamma_other->Fill(y[1],weight[1]);
                                    //ym_gamma_other->Fill(ym);
                                    //verif_otherevents ++;
                                    ym_gamma_other_value[0]= y[0];
                                    ym_gamma_other_value[1]= y[1];
                                    num_gamma_other_recons++;
                                }

                                if(gammatype=="ScattPrompt")
                                {
                                    ym_gamma_scattered->Fill(y[0],weight[0]);
                                    ym_gamma_scattered->Fill(y[1],weight[1]);
                                    //ym_gamma_scattered->Fill(ym);
                                    ym_gamma_prompt_scatt[0] = y[0];
                                    ym_gamma_prompt_scatt[1] = y[1];
                                    num_gamma_scatt_recons++;

                                }
                                reconstree->Fill();

                            }

                            //if(TMath::Abs(ym)<300)rectree->Fill();

                            if (Sameparticule)
                            {
                                //ym_histo_false->Fill(ym-25);
                                ym_histo_true_same_particle_no_gamma->Fill(y[0],weight[0]);
                                ym_histo_true_same_particle_no_gamma->Fill(y[1],weight[1]);
                                //  ym_histo_true_same_particle_no_gamma->Fill(ym);

                                ym_histo->Fill(y[0],weight[0]);
                                ym_histo->Fill(y[1],weight[1]);
                               // ym_histo->Fill(ym);

                                // cout<<"ym_same "<<ym<<endl;
                                Sameparticule_number_histo ++;
                                //ym_histo_false->Fill(ym-25);
                                ym_sameparticle[0] =y[0];
                                ym_sameparticle[1] =y[1];
                                ym_background[0]=y[0];
                                ym_background[1]=y[1];

                                reconstree->Fill();
                            }

                            if (Diffparticule)
                            {
                                //ym_histo_true_same_gamma->Fill(ym-25);
                                // ym_histo_true_different_particle->Fill(y[0]-25);
                                // ym_histo_true_different_particle->Fill(y[1]-25);
                                ym_histo->Fill(y[0],weight[0]);
                                ym_histo->Fill(y[1],weight[1]);
                                //ym_histo->Fill(ym);

                                ym_histo_true_same_particle_no_gamma->Fill(y[0],weight[0]);
                                ym_histo_true_same_particle_no_gamma->Fill(y[1],weight[1]);
                               // ym_histo_true_different_particle->Fill(ym);
                                //cout<<"ym_diff "<<ym<<endl;
                                Diffparticule_number_histo ++;

                                ym_different[0] =y[0];
                                ym_different[1] =y[1];
                                ym_background[0]=y[0];
                                ym_background[1]=y[1];

                                reconstree->Fill();
                            }


                            //if(!samepscoinc)ym_histo->Fill(ym-25);


                            reconstree_all->Fill();
                        }



                        /*ym=y[1];
                         ym_histo->Fill(ym);
                         if(samepscoinc)ym_histo_true->Fill(ym);
                         if(!samepscoinc)ym_histo_false->Fill(ym);
                         if(gammatype=="PromptGamma" && samepscoinc)ym_histo_prompt->Fill(ym);
                         if(TMath::Abs(ym)<300)rectree->Fill();
                         */
                        //cout<<"resultat --- Nion "<<mycoinc.hitLYSO.Nion<<" "<<y[0]<<" "<<y[1]<<endl;
                    }
                }
                else n_null++;
                //cout<<"Tof Si : "<<tofSi<<" to : "<<mycoinc.hitLYSO.t0<<"t max "<<mycoinc.hitLYSO.tmax<<endl;

                //cout<<"tmin : "<<mycoinc.hitLYSO.tmin<<" tmax: "<<mycoinc.hitLYSO.tmax<<"diff:"<<mycoinc.hitLYSO.tmax-mycoinc.hitLYSO.tmin<<endl;

                /*tof_histo->Fill(tofSi);
                 if(gammacoinc)tof_histo_true->Fill(tofSi);
                 if(!samepscoinc)tof_histo_false->Fill(tofSi);*/

                //mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.tmax

                //if(ParticleTypeLYSO!="gamma")cout<<" diff BGO : "<<mycoinc.hitLYSO.tmin-mycoinc.hitLYSO.t0<<" Particle BGO : "<<ParticleTypeLYSO<< " Particule Si : "<<ParticleTypeSi<<endl<<endl;

                //tof_histo_true->Fill(mycoinc.hitSi.tmin-mycoinc.hitSi.t0);
                //if(samepscoinc)tof_histo_true->Fill(mycoinc.hitLYSO.tmax);


                /*tof_histo->Fill(tofLYSO);
                 if(gammacoinc)tof_histo_true->Fill(tofLYSO);
                 if(!samepscoinc)tof_histo_false->Fill(tofLYSO);*/
            }
        }
        //cout<<"reconstruction finie, fichier "<<newstrname.c_str()<<" ecrit"<<endl;
        //ComptonFitHisto(ym_histo,nbin);


        //---->>>>>>>>>>recfile->Write();
        ReconstructionTree->Write();

        //cout<<"Ecriture du fichier reconstruction finie"<<endl;
    }
    else cout<<" Pas de coincidences trouv�es "<<endl;


    cout<<" Nombre total de 2 hits "<<ncoinc<<" 2 hits gamma "<<ngammacoinc<<endl;
    cout<<" % de coinc avec la meme particule incidente "<<(100.*nsamepscoinc/ncoinc)<<endl;
    cout<<" % de coinc avec le meme gamma incident "<<(100.*ngammacoinc/ncoinc)<<endl;
    cout<<" % de coinc gamma - true event "<<(100.*ntrueeventcoinc/ncoinc)<<endl;
    cout<<" Nombre de 2 hit/ p incident "<<ncoinc/Nion<<endl;
    cout<<" Nombre de coinc avec la meme particule incidente / p incident "<<(nsamepscoinc/Nion)<<endl;
    cout<<" Nombre de coinc avec le meme gamma incident / p incident "<<(ngammacoinc/Nion)<<endl;
    cout<<" Nombre de coinc gamma - true event / p incident "<<(ntrueeventcoinc/Nion)<<endl;


    /* cout<<endl<<"--------data comme rapport MH-------"<<endl<<endl;
     cout<<" % de coinc avec la meme particule incidente "<<(100.*nsamepscoinc/ncoinc)-(100.*ngammacoinc/ncoinc)<<endl;
     cout<<" % de coinc avec le meme gamma incident "<<(100.*ngammacoinc/ncoinc)-(100.*ntrueeventcoinc/ncoinc)<<endl;
     cout<<" % de coinc gamma - true event "<<(100.*ntrueeventcoinc/ncoinc)<<endl;
     cout<<" % de fortuits "<<100*((ncoinc-nsamepscoinc)/ncoinc)<<endl;
     cout<<" Nombre de 2 hit/ p incident "<<ncoinc/Nion<<endl;
     cout<<" Nombre de coinc avec la meme particule incidente / p incident "<<(nsamepscoinc-ngammacoinc)/Nion<<endl;
     cout<<" Nombre de coinc avec le meme gamma incident / p incident "<<(ngammacoinc-ntrueeventcoinc)/Nion<<endl;
     cout<<" Nombre de coinc gamma - true event / p incident "<<(ntrueeventcoinc/Nion)<<endl;
     cout<<" Nombre de coinc fortuits / p incident "<<(ncoinc-nsamepscoinc)/Nion<<endl;
     cout<<" Nombre de coinc vraies / p incident "<<(nsamepscoinc)/Nion<<endl;
     //cout<<" nsamepscoinc "<<nsamepscoinc<<" Nion "<<Nion<<endl;*/

    cout << " nombre_nion : "<<nombre_nion<<endl;

    cout<<" nvrai ym : "<<nvraie_ym<<endl;
    cout<<" ngammatrue_ym : "<<ngammatrue_ym<<endl;
    cout<<" nrandom : "<<nrandom_ym<<endl;
    cout<<" n_null : "<<n_null<<endl;
    cout<<"nvrai " <<n_vrai<<endl;
    cout<<"nfortuit "<<n_fortuits<<endl;
    cout<< " n same particule " << Sameparticule_number<<endl;
    cout<< " n diff particule " << Diffparticule_number<<endl;
    cout<<" ngamma coin : "<<ngammacoinc<<endl;
    cout<<" ngamma prompt "<<ntrueeventcoinc<<endl;

    cout <<"ncoinc : "<<ncoinc<<" Nion : "<<Nion<<endl;

    cout<<" % de fortuits "<<100.*((n_fortuits)/ncoinc)<<endl;
    cout<<" % de vrais "<<100.*((n_vrai)/ncoinc)<<endl;
    cout<<" % de vrai  gamma "<<100.*(ngammacoinc/ncoinc)<<endl;
    cout<<" % de vrai autre "<<100.*((n_vrai-ngammacoinc)/ncoinc)<<endl;
    cout<<" % de vrai same particule "<<100.*((Sameparticule_number)/ncoinc)<<endl;

    cout<<" Nombre de coinc fortuits / p incident "<<(n_fortuits)/Nion<<endl;
    cout<<" Nombre de coinc vrai / p incident "<<n_vrai/Nion<<endl;
    cout<<" Nombre de vrai  gamma / p incident "<<(ngammacoinc)/Nion<<endl;
    cout<<" Nombre de vrai autre / p incident "<<(n_vrai-ngammacoinc)/Nion<<endl;

    cout<< " n same particule/ p incident " << Sameparticule_number/ncoinc<<endl;
    cout<< " n diff particule / p incident" << Diffparticule_number/ncoinc<<endl;
    cout<<" ngamma coin / p incident "<<ngammacoinc/ncoinc<<endl;

    cout<<"count_sumAccurate "<<count_sumAccurate<<endl;

    cout << " event MLEM " <<neventstore<<endl;

    for (int s =0;s< nombre_coinc_reconstruite;s++)
    {
        moyenne_Histo_difference_y0 = Histo_difference_y0[s] +moyenne_Histo_difference_y0;
        moyenne_Histo_difference_y1 = Histo_difference_y1[s] +moyenne_Histo_difference_y1;
        moyenne_Histo_difference_y0y1 = Histo_difference_y0y1[s] +moyenne_Histo_difference_y0y1;
    }

    cout<< " moyenne_Histo_difference_y0 av " <<moyenne_Histo_difference_y0 << " moyenne_Histo_difference_y1 av " << moyenne_Histo_difference_y1 <<" moyenne_Histo_difference_y0y1 av "<<moyenne_Histo_difference_y0y1 <<" nbr events " <<nombre_coinc_reconstruite<<endl;
    moyenne_Histo_difference_y0 = moyenne_Histo_difference_y0 / nombre_coinc_reconstruite;
    moyenne_Histo_difference_y1 = moyenne_Histo_difference_y1 / nombre_coinc_reconstruite;
    moyenne_Histo_difference_y0y1 = moyenne_Histo_difference_y0y1 / nombre_coinc_reconstruite;

    cout<< " moyenne_Histo_difference_y0 " <<moyenne_Histo_difference_y0 << " moyenne_Histo_difference_y1 " << moyenne_Histo_difference_y1 <<" moyenne_Histo_difference_y0y1 "<<moyenne_Histo_difference_y0y1 <<" nbr events " <<nombre_coinc_reconstruite<<endl;
    cout<<" test sum "<<test_Sum<<endl;
    cout<<"others "<<verif_otherevents<<endl;

    cout<<endl<<endl<<" gamma coinc "<<num_gamma_init<<" prompt coinc "<<num_gamma_prompt_ini<<" scatt coin "<<num_gamma_scatt_init<< " other init "<<num_gamma_other<<endl;

    cout<<" gamma recons "<<num_gamma_recons<<" prompt recons "<<num_gamma_prompt_recons<<" scatt recons"<<num_gamma_scatt_recons<< " other recons "<<num_gamma_other_recons<<endl;

    //Definition des noms des fichiers de sauvegarde pour les taux de coincidences absolus, relatifs et reconstruits (histo)

    string	  savefiletxtbasis;

    string	  savefile_total_relatif;
    string	  savefile_aleatoire_relatif;
    string	  savefile_vraiesGamma_relatif;
    string	  savefile_vraies_relatif;
    string	  savefile_vraiesAutres_relatif;
    string	  savefile_vraiesMemeParticule_relatif;

    string	  savefile_aleatoire_absolu;
    string	  savefile_vraies_absolu;
    string	  savefile_vraiesGamma_absolu;
    string	  savefile_vraiesAutres_absolu;
    string	  savefile_vraiesMemeParticule_absolu;

    string	  savefile_total_Histo;
    string	  savefile_vraies_Histo;
    string	  savefile_random_Histo;
    string	  savefile_Gamma_Histo;
    string	  savefile_vraiesMemeParticule_absolu_Histo;
    string	  savefile_DiffParticule_absolu_Histo;



    //Definition de la base des fichiers � sauvegarder

    if(cut == 0){
        //savefiletxtbasis="2015_01_16_1Proton_160MeV_9_42nsBunch_2nsSizeBunch_40nswindowscoinc_ResolSi15ns_ResolBGO3ns_CylinderPMMA_15cm20cm_High_stat_Gamma_NewFiles1e9_NOCutTOF";
          savefiletxtbasis= "/sps/hep/hadronth/mfontana/results_CCanalysis/saveFiletxtBasisNoCut";}
    if(cut == 1){
      savefiletxtbasis="/sps/hep/hadronth/mfontana/results_CCanalysis/saveFiletxtBasisCut";}
        //savefiletxtbasis="2015_01_16_1Proton_160MeV_9_42nsBunch_2nsSizeBunch_40nswindowscoinc_ResolSi15ns_ResolBGO3ns_CylinderPMMA_15cm20cm_High_stat_Gamma_NewFiles1e9_CutTOF_0_6ns";
    //2014_12_30_1Proton_160MeV_9_42nsBunch_2nsSizeBunch_40nswindowscoinc_ResolSi15ns_ResolBGO3ns_CylinderPMMA_15cm20cm_High_stat_NOCutTOF


    //relatif avec cut
    savefile_total_relatif = savefiletxtbasis+"_total_relatif.txt";
    savefile_aleatoire_relatif = savefiletxtbasis+"_fortuits_relatif.txt";
    savefile_vraies_relatif = savefiletxtbasis+"_vraies_relatif.txt";
    savefile_vraiesGamma_relatif = savefiletxtbasis+"_vraiesGamma_relatif.txt";
    savefile_vraiesAutres_relatif = savefiletxtbasis+"_vraiesAutres_relatif.txt";
    savefile_vraiesMemeParticule_relatif = savefiletxtbasis+"_vraiesMemeParticule_relatif.txt";

    //absolu avec cut
    savefile_aleatoire_absolu = savefiletxtbasis+"_aleatoire_absolu.txt";
    savefile_vraies_absolu = savefiletxtbasis+"_vraies_absolu.txt";
    savefile_vraiesGamma_absolu = savefiletxtbasis+"_vraiesGamma_absolu.txt";
    savefile_vraiesAutres_absolu = savefiletxtbasis+"_vraiesAutres_absolu.txt";
    savefile_vraiesMemeParticule_absolu = savefiletxtbasis+"_vraiesMemeParticule_absolu.txt";

    //histo avec cut
    savefile_total_Histo = savefiletxtbasis+"_total_Histo.txt";
    savefile_vraies_Histo = savefiletxtbasis+"_vraies_Histo.txt";
    savefile_random_Histo = savefiletxtbasis+"_random_Histo.txt";
    savefile_Gamma_Histo = savefiletxtbasis+"_Gamma_Histo.txt";
    savefile_vraiesMemeParticule_absolu_Histo = savefiletxtbasis+"_vraiesMemeParticule_absolu_Histo.txt";
    savefile_DiffParticule_absolu_Histo = savefiletxtbasis+"_DiffParticule_absolu_Histo.txt";

    if(save_data ==1)
    {

        if(cut==1)
        {
            //Relatif

            ofstream fichier26(savefile_total_relatif.c_str(),ios::out | ios::app);
            if (fichier26) fichier26<<ncoinc<<",";
            fichier26.close();

            ofstream fichier30(savefile_aleatoire_relatif.c_str(),ios::out | ios::app);
            if (fichier30) fichier30<<100.*((n_fortuits)/ncoinc)<<",";
            fichier30.close();

            ofstream fichier2(savefile_vraies_relatif.c_str(),ios::out | ios::app);
            if (fichier2) fichier2<<100.*((n_vrai)/ncoinc)<<",";
            fichier2.close();

            ofstream fichier3(savefile_vraiesGamma_relatif.c_str(),ios::out | ios::app);
            if (fichier3) fichier3<<100.*((ngammacoinc)/ncoinc)<<",";
            fichier3.close();

            ofstream fichier4(savefile_vraiesAutres_relatif.c_str(),ios::out | ios::app);
            if (fichier4) fichier4<<100.*(Diffparticule_number/ncoinc)<<",";
            fichier4.close();

            ofstream fichier5(savefile_vraiesMemeParticule_relatif.c_str(),ios::out | ios::app);
            if (fichier5) fichier5<<100.*((Sameparticule_number)/ncoinc)<<",";
            fichier5.close();


            //Absolu

            ofstream fichier6(savefile_aleatoire_absolu.c_str(),ios::out | ios::app);
            if (fichier6) fichier6<<((n_fortuits)/Nion)<<",";
            fichier6.close();

            ofstream fichier7(savefile_vraies_absolu.c_str(),ios::out | ios::app);
            if (fichier7) fichier7<<((n_vrai)/Nion)<<",";
            fichier7.close();

            ofstream fichier8(savefile_vraiesGamma_absolu.c_str(),ios::out | ios::app);
            if (fichier8) fichier8<<((ngammacoinc)/Nion)<<",";
            fichier8.close();

            ofstream fichier9(savefile_vraiesAutres_absolu.c_str(),ios::out | ios::app);
            if (fichier9) fichier9<<(Diffparticule_number/Nion)<<",";
            fichier9.close();

            ofstream fichier10(savefile_vraiesMemeParticule_absolu.c_str(),ios::out | ios::app);
            if (fichier10) fichier10<<((Sameparticule_number)/Nion)<<",";
            fichier10.close();


            //Events reconstruction

            ofstream fichier20(savefile_total_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier20) fichier20<<nombre_coinc_reconstruite<<",";
            fichier20.close();


            ofstream fichier21(savefile_vraies_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier21) fichier21<<nvraie_ym<<",";
            fichier21.close();


            ofstream fichier22(savefile_random_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier22) fichier22<<nrandom_ym<<",";
            fichier22.close();


            ofstream fichier23(savefile_Gamma_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier23) fichier23<<ngammatrue_ym<<",";
            fichier23.close();


            ofstream fichier24(savefile_vraiesMemeParticule_absolu_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier24) fichier24<<Sameparticule_number_histo<<",";
            fichier24.close();



            ofstream fichier25(savefile_DiffParticule_absolu_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier25) fichier25<<Diffparticule_number_histo<<",";
            fichier25.close();


        }



        //--------------------------------------
        //---------- Sans TOF cut --------------
        //--------------------------------------

        if(cut==0)
        {
            //Relatif

            ofstream fichier26(savefile_total_relatif.c_str(),ios::out | ios::app);
            if (fichier26) fichier26<<ncoinc<<",";
            fichier26.close();

            ofstream fichier30(savefile_aleatoire_relatif.c_str(),ios::out | ios::app);
            if (fichier30) fichier30<<100.*((n_fortuits)/ncoinc)<<",";
            fichier30.close();

            ofstream fichier2(savefile_vraies_relatif.c_str(),ios::out | ios::app);
            if (fichier2) fichier2<<100.*((n_vrai)/ncoinc)<<",";
            fichier2.close();

            ofstream fichier3(savefile_vraiesGamma_relatif.c_str(),ios::out | ios::app);
            if (fichier3) fichier3<<100.*((ngammacoinc)/ncoinc)<<",";
            fichier3.close();

            ofstream fichier4(savefile_vraiesAutres_relatif.c_str(),ios::out | ios::app);
            if (fichier4) fichier4<<100.*(Diffparticule_number/ncoinc)<<",";
            fichier4.close();

            ofstream fichier5(savefile_vraiesMemeParticule_relatif.c_str(),ios::out | ios::app);
            if (fichier5) fichier5<<100.*((Sameparticule_number)/ncoinc)<<",";
            fichier5.close();


            //Absolu

            ofstream fichier6(savefile_aleatoire_absolu.c_str(),ios::out | ios::app);
            if (fichier6) fichier6<<((n_fortuits)/Nion)<<",";
            fichier6.close();

            ofstream fichier7(savefile_vraies_absolu.c_str(),ios::out | ios::app);
            if (fichier7) fichier7<<((n_vrai)/Nion)<<",";
            fichier7.close();

            ofstream fichier8(savefile_vraiesGamma_absolu.c_str(),ios::out | ios::app);
            if (fichier8) fichier8<<((ngammacoinc)/Nion)<<",";
            fichier8.close();

            ofstream fichier9(savefile_vraiesAutres_absolu.c_str(),ios::out | ios::app);
            if (fichier9) fichier9<<(Diffparticule_number/Nion)<<",";
            fichier9.close();

            ofstream fichier10(savefile_vraiesMemeParticule_absolu.c_str(),ios::out | ios::app);
            if (fichier10) fichier10<<((Sameparticule_number)/Nion)<<",";
            fichier10.close();


            //Events reconstruction

            ofstream fichier20(savefile_total_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier20) fichier20<<nombre_coinc_reconstruite<<",";
            fichier20.close();


            ofstream fichier21(savefile_vraies_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier21) fichier21<<nvraie_ym<<",";
            fichier21.close();


            ofstream fichier22(savefile_random_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier22) fichier22<<nrandom_ym<<",";
            fichier22.close();


            ofstream fichier23(savefile_Gamma_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier23) fichier23<<ngammatrue_ym<<",";
            fichier23.close();


            ofstream fichier24(savefile_vraiesMemeParticule_absolu_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier24) fichier24<<Sameparticule_number_histo<<",";
            fichier24.close();



            ofstream fichier25(savefile_DiffParticule_absolu_Histo.c_str(),ios::out | ios::app); // pour enregsitrer les graines dans le fichier sauvegarde.txt
            if (fichier25) fichier25<<Diffparticule_number_histo<<",";
            fichier25.close();

        }
    }

    cout <<"test largeur TFH "<<nombre_coinc_reconstruite<<" "<<nvraie_ym<<" "<<nrandom_ym<<" "<<ngammatrue_ym<<" "<<Sameparticule_number_histo<<" "<<Diffparticule_number_histo<<endl;
    cout <<"test number total "<<ncoinc<<" "<<n_vrai<<" "<<n_fortuits<<" "<<ngammacoinc<<" "<<Sameparticule_number<<" "<<Diffparticule_number<<endl;
    cout<<"ngammatrue_ym histo "<<ngammatrue_ym<< " ngammacoinc normal "<<endl;
    cout<<"gamma out "<<count_gamma_out<<endl;

    cout<<"emoins "<<ne_moins/nemoins_sum<<endl;
    cout<<"proton "<<nprotons/nprotons_sum<<endl;
    cout<<"ne_plus "<<ne_plus/neplus_sum<<endl;
    cout<<"neutron "<<n_neutrons/n_neutrons_sum<<endl;

    cout<<"num coinc comparaison HIT " <<num_coinc<<endl;

    ReconstructionTree->Close();

    //-------------------oooooOOOOO00000OOOOOooooo---------------------#
    //                                                                 #
    //   		   	AFFICHAGE DES RESULTATS		           #
    //                                                                 #
    //-------------------oooooOOOOO00000OOOOOooooo---------------------#
    ym_histo->Scale(1/Nion);
    ym_histo_true->Scale(1/Nion);
    ym_histo_false->Scale(1/Nion);
    ym_histo_prompt->Scale(1/Nion);
    ym_histo_true_same_gamma->Scale(1/Nion);
    ym_histo_true_same_particle_no_gamma->Scale(1/Nion);
    ym_histo_true_different_particle->Scale(1/Nion);
    ym_vertex->Scale(1/Nion);
    ym_gamma_other->Scale(1/Nion);
    ym_gamma_scattered->Scale(1/Nion);
    edep_histo->Scale(1/Nion);
    edep_histo_true->Scale(1/Nion);
    edep_histo_false->Scale(1/Nion);
    edep_histo_true_gamma->Scale(1/Nion);
    histo_distancey0y1->Scale(1/Nion);
    histo_distancey0y1_cut1MeV_3MeV->Scale(1/Nion);
    histo_distancey0y1_cut0_1MeV->Scale(1/Nion);
    histo_distancey0y1_cut3MeV_sup->Scale(1/Nion);
    histo_distancey0y1_enterTarget->Scale(1/Nion);
    histo_distancey0y1_FallOffTarget->Scale(1/Nion);
    histo_distancey0y1_outTarget->Scale(1/Nion);
    histo_distancey0y1_cut1MeV_3MeV_edep->Scale(1/Nion);
    histo_distancey0y1_cut3MeV_sup_edep->Scale(1/Nion);
    ym_histo_autre_solution->Scale(1/Nion);


    cout<<"histo_distancey0y1 mean"<<histo_distancey0y1->GetMean(1)<<" rms "<<histo_distancey0y1->GetRMS(1)<<endl;
    cout<<"histo_distancey0y1_cut1MeV_3MeV mean"<<histo_distancey0y1_cut1MeV_3MeV->GetMean(1)<<" rms "<<histo_distancey0y1_cut1MeV_3MeV->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_cut3MeV_sup mean"<<histo_distancey0y1_cut3MeV_sup->GetMean(1)<<" rms "<<histo_distancey0y1_cut3MeV_sup->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_enterTarget mean"<<histo_distancey0y1_enterTarget->GetMean(1)<<" rms "<<histo_distancey0y1_enterTarget->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_FallOffTarget mean"<<histo_distancey0y1_FallOffTarget->GetMean(1)<<" rms "<<histo_distancey0y1_FallOffTarget->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_outTarget mean"<<histo_distancey0y1_outTarget->GetMean(1)<<" rms "<<histo_distancey0y1_outTarget->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_cut1MeV_3MeV_edep mean"<<histo_distancey0y1_cut1MeV_3MeV_edep->GetMean(1)<<" rms "<<histo_distancey0y1_cut1MeV_3MeV_edep->GetRMS(1)<<endl;
     cout<<"histo_distancey0y1_cut3MeV_sup_edep mean"<<histo_distancey0y1_cut3MeV_sup_edep->GetMean(1)<<" rms "<<histo_distancey0y1_cut3MeV_sup_edep->GetRMS(1)<<endl;

    cout<<"dadaddadad"<<endl;

    int binmax=0;
    int ValueHisto=0;
    //int =0; 01/12/2014 effacement de la variable initialis�... � retrouver ^^
    int maxHisto =0;
    int maxbincalibre =0;

     TH1 *frame = new TH1F("frame","1",200,-180,130);   //200,-180,130);
     frame->SetMinimum(0);
     frame->SetMaximum(0.0000007);//  0.007
     frame->SetDirectory(0);
     frame->SetStats(0);
     frame->SetTitle("");
     //	frame->	GetXaxis()->SetLabelSize(20);
     //	frame->	GetYaxis()->SetLabelSize(20);
     //frame->GetYaxis()->SetLabelSize(0.001);
     frame->GetXaxis()->SetTitleSize(0.05);
     frame->GetXaxis()->SetTitle("Beam axis [mm]");
     frame->GetXaxis()->SetTickLength(0.02);
     frame->GetXaxis()->SetLabelSize(0.05);
     frame->GetYaxis()->SetTitleSize(0.05);
     frame->GetYaxis()->SetTitle("#splitline{Number of events}{per incident proton}");
     frame->GetYaxis()->SetLabelSize(0.05);
     frame->GetYaxis()->SetTitleOffset(1.5);
     frame->Draw("");
     ym_histo->SetLineColor(kBlack);
     ym_histo->SetLineWidth(2);
     ym_histo->SetLineStyle(2);
     ym_histo->Draw("same");

     ym_histo_prompt->SetLineColor(kCyan);
     ym_histo_prompt->Draw("same");
     ym_gamma_other->SetLineColor(kRed);
     ym_gamma_other->Draw("same");
     ym_gamma_scattered->SetLineColor(kGreen);
     ym_gamma_scattered->Draw("same");
     //ym_histo_true->SetLineColor(kGreen);
     //ym_histo_true->Draw("same");
     // ym_histo_false->SetLineColor(kRed);//
     // ym_histo_false->Draw("same");
     ym_histo_true_same_gamma->SetLineColor(kBlue);
     ym_histo_true_same_gamma->SetLineWidth(2);
     ym_histo_true_same_gamma->Draw("same");

   // ym_histo_autre_solution->SetLineColor(kPink);
    //ym_histo_autre_solution->SetLineWidth(2);
    //ym_histo_autre_solution->Draw("same");

     TLegend * leg = new TLegend(0.3,0.75,0.85,0.89);
     leg->SetFillColor(10);
     leg->AddEntry(ym_histo, "Background","l");
     leg->AddEntry(ym_histo_true_same_gamma,"True gamma","l");
     leg->AddEntry(ym_histo_prompt,"Prompt gamma","l");
     leg->AddEntry(ym_gamma_scattered,"Prompt gamma scattered","l");
     leg->AddEntry(ym_gamma_other,"other gamma ","l");
    //leg->AddEntry(ym_histo_autre_solution,"other solutions all ","l");

     leg->Draw();


    TCanvas *c200 = new TCanvas("c200","Graph",0,0,700,600);
     c200->SetFillColor(0);
     c200->SetBorderMode(0);
     gPad->SetBorderMode(0);

    histo_distancey0y1->GetXaxis()->SetTitle("Solutions difference [mm]");
    histo_distancey0y1->GetYaxis()->SetTitle("#splitline{Number of events}{per incident protons}");
    histo_distancey0y1->SetLineColor(kBlack);
    histo_distancey0y1->SetLineWidth(2);
    histo_distancey0y1->SetLineStyle(1);
    histo_distancey0y1->Draw("same");
    /*histo_distancey0y1_cut0_1MeV->SetLineColor(kCyan);
    histo_distancey0y1_cut0_1MeV->SetLineWidth(2);
    histo_distancey0y1_cut0_1MeV->SetLineStyle(1);
    histo_distancey0y1_cut0_1MeV->Draw("same");*/
    histo_distancey0y1_cut1MeV_3MeV->SetLineColor(kRed);
    histo_distancey0y1_cut1MeV_3MeV->SetLineWidth(2);
    histo_distancey0y1_cut1MeV_3MeV->SetLineStyle(1);
    histo_distancey0y1_cut1MeV_3MeV->Draw("same");
    histo_distancey0y1_cut3MeV_sup->SetLineColor(kGreen);
    histo_distancey0y1_cut3MeV_sup->SetLineWidth(2);
    histo_distancey0y1_cut3MeV_sup->SetLineStyle(1);
    histo_distancey0y1_cut3MeV_sup->Draw("same");
    histo_distancey0y1_cut1MeV_3MeV_edep->SetLineColor(kRed);
    histo_distancey0y1_cut1MeV_3MeV_edep->SetLineWidth(2);
    histo_distancey0y1_cut1MeV_3MeV_edep->SetLineStyle(2);
    histo_distancey0y1_cut1MeV_3MeV_edep->Draw("same");
    histo_distancey0y1_cut3MeV_sup_edep->SetLineColor(kGreen);
    histo_distancey0y1_cut3MeV_sup_edep->SetLineWidth(2);
    histo_distancey0y1_cut3MeV_sup_edep->SetLineStyle(2);
    histo_distancey0y1_cut3MeV_sup_edep->Draw("same");

    TLegend * leg2 = new TLegend(0.3,0.75,0.85,0.89);
    leg2->SetFillColor(10);
    leg2->AddEntry(histo_distancey0y1, "All","l");
  //  leg2->AddEntry(histo_distancey0y1_cut0_1MeV,"Cut [0,1] MeV","l");
    leg2->AddEntry(histo_distancey0y1_cut1MeV_3MeV,"Cut ]1,3] MeV - at emission","l");
    leg2->AddEntry(histo_distancey0y1_cut3MeV_sup,"Cut > 3 MeV- at emission","l");
    leg2->AddEntry(histo_distancey0y1_cut1MeV_3MeV_edep,"Cut ]1,3] MeV - deposited","l");
    leg2->AddEntry(histo_distancey0y1_cut3MeV_sup_edep,"Cut > 3 MeV - deposited","l");

    /*leg2->AddEntry(histo_distancey0y1_enterTarget,"Enter target : -100 -> 0 mm","l");
    leg2->AddEntry(histo_distancey0y1_FallOffTarget,"Falloff target : 0 -> 65 mm","l");
    leg2->AddEntry(histo_distancey0y1_outTarget,"End target : 65 -> 100 mm","l");*/
    leg2->Draw();



    TCanvas *c300 = new TCanvas("c300","Graph",0,0,700,600);
    c300->SetFillColor(0);
    c300->SetBorderMode(0);
    gPad->SetBorderMode(0);

    histo_difference_weight->GetXaxis()->SetTitle("weight");
    histo_difference_weight->GetYaxis()->SetTitle("Number of events");
    histo_difference_weight->SetLineColor(kBlack);
    histo_difference_weight->SetLineWidth(2);
    histo_difference_weight->SetLineStyle(1);
    histo_difference_weight->Draw("");


    //Fit du profil reconstruit
    Double_t* parameters=new Double_t[10];

     TCanvas *c4 = new TCanvas("c4","Graph",0,0,700,600);
     c4->SetFillColor(0);
     c4->SetBorderMode(0);
     gPad->SetBorderMode(0);
     c4->Divide(2,2);
     c4->cd(1);
     edep_histo->Draw("col2z");
     edep_histo->GetXaxis()->SetTitle("E deposit Si [MeV]");
     edep_histo->GetYaxis()->SetTitle("E deposit BGO [MeV]");
     c4->cd(2);
     edep_histo_true->Draw("col2z");
     edep_histo_true->GetXaxis()->SetTitle("E deposit Si [MeV] ");
     edep_histo_true->GetYaxis()->SetTitle("E deposit BGO [MeV]");

     //edep_histo_true->GetXaxis()->SetTitle("E deposit Si [MeV]");
     //edep_histo_true->GetYaxis()->SetTitle("E deposit BGO [MeV]");
     c4->cd(3);
     edep_histo_false->Draw("col2z");
     edep_histo_false->GetXaxis()->SetTitle("E deposit Si [MeV]");
     edep_histo_false->GetYaxis()->SetTitle("E deposit BGO [MeV]");

     c4->cd(4);
     edep_histo_true_gamma->Draw("col2z");
     edep_histo_true_gamma->GetXaxis()->SetTitle("edep Si [MeV]");
     edep_histo_true_gamma->GetYaxis()->SetTitle("edep BGO [MeV]");



    TCanvas *c50 = new TCanvas("c3","Graph",0,0,700,600);
     c50->SetFillColor(0);
     c50->SetBorderMode(0);
     gPad->SetBorderMode(0);
     spectre_energy_gamma->GetXaxis()->SetTitle("Energy[MeV]");
     spectre_energy_gamma->GetYaxis()->SetTitle("Events");
     spectre_energy_gamma->SetLineColor(kBlack);
     spectre_energy_gamma->Draw("");


    //-------------------oooooOOOOO00000OOOOOooooo---------------------#
    //                                                                 #
    //   		  Enregistrement des histos		           #
    //                                                                 #
    //-------------------oooooOOOOO00000OOOOOooooo---------------------#

    string histo_reconstruction;

    if(cut==0){//histo_reconstruction="/sps/hep/hadronth/JLL/Macro_Root/Coincidences_reelles_fortuites/code/code2/Shift_camera_study/Results_Shift_plus5_cm/2015_03_24_Protons_160MeV_PMMACylinder_Reconstructed_profile_NoCutTOF_HighStat_amera_plus5cm_"+num_file+".root";
      histo_reconstruction =  outpath + name +"noCut_"+num_file+".root";
    }
    if(cut==1){//histo_reconstruction="/sps/hep/hadronth/JLL/Macro_Root/Coincidences_reelles_fortuites/code/code_protons2/Shift_camera_study/Results_Shift_moins5_cm/2015_03_31_Protons_160MeV_PMMACylinder_Reconstructed_profile_CutTOF_HighStat_camera_moins5cm_"+num_file+".root";
      histo_reconstruction =  outpath + name +"Cut_"+num_file+".root";
    }

    //2015_01_16_Protons_160MeV_PMMACylinder_Reconstructed_profile_Show_Gamma_NewFiles_CutTOF_TEST"+num_file+".root";

    //2014_12_01_Protons_160MeV_PMMACylinder_Reconstructed_profile_NoCutTOF_
    //2014_12_01_Protons_160MeV_PMMACylinder_Reconstructed_profile_CutTOF_
    //2014_12_27_Protons_160MeV_PMMACube20x20xm_d020cm_Reconstructed_profile_CutTOF_test
    //TFile *f = new TFile("/sps/hep/hadronth/mfontana/results_CCanalysis/reconstruction.root","RECREATE");
    TFile *f = new TFile("./results/reconstruction.root","RECREATE");

    ym_histo->Write("background");
    // ym_histo_true->Write("ym_histo_true");
    ym_histo_true_same_gamma->Write("True gamma");

    f->Close();

    //theApp->Run();

    return ;

}
