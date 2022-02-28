//TCanvas *c1 = new TCanvas();

TH1I* gHistPDG;
TH1D* gHistJpsiAll_pt;
TH1D* gHistJpsiAll_Y;
TH1D* gHistJpsiAll_eta;
TH1D* gHistJpsiPrompt_pt;
TH1D* gHistGammaGammaMass;



#include <iostream>
#include <string>

using std::cout; using std::cin;
using std::endl; using std::string;

void ReadKinematics()
{
	//TFile *input = new TFile("/home/paulo/testes/Grid/Starlight/run9/5/Kinematics.root","read");
	
/*    LOAD OVER ALL FOLDERS
	char folder[200];
	for(int i = 0 ; i<5;i++){
		sprintf(folder, "Event%d", i);
		
		TTree *tree = (TTree*)input->Get(Form("%s/TreeK", folder ));
	
		int entries = tree->GetEntries();} 
		
*/
		
	//TTree *tree = (TTree*)input->Get("Event0/TreeK");
	
	TFile* save = TFile::Open("Results_Kinematics_run12.root", "recreate");


//===================HISTOGRAMS=========================
			
	gHistPDG = new TH1I("HistPDG", "PDG codes", 457, -12., 445.);
	gHistJpsiAll_pt = new TH1D("HistJpsiAll_pt", "Inclusive J/#psi; p_{T} (GeV/c)", 1000, 0., 1.);
	gHistJpsiAll_Y = new TH1D("HistJpsiAll_Y", "Inclusive J/#psi; y", 100, 2., 6.);
	gHistJpsiAll_eta = new TH1D("HistJpsiAll_Eta", "Inclusive J/#psi; #eta", 100, 6., 10.);
	gHistJpsiPrompt_pt = new TH1D("HistJpsiPrompt_pt", "Prompt J/#psi; p_{T} (GeV/c)", 1000, 0., 1.);
	
	gHistGammaGammaMass = new TH1D("HistGammaGammaMass", "#gamma #gamma photo-produced continuum; mass (GeV/c^{2})", 100, 0., 6.);
	
//=========================================================	
	
	AliRunLoader *rl = AliRunLoader::Open("run12/1/galice.root");
	
	rl->LoadKinematics();
	
//======= LOOP OVER EVENTS =================================================================	
	
	for(int iev = 0; iev < rl->GetNumberOfEvents();iev++)
	{
		rl->GetEvent(iev);
		
		double p[2][3];
		int electronCount = 0;
		
		//========LOOP OVER PARTICLES===========
		AliStack *stack = rl->Stack();
		for(int ipart=0; ipart < stack->GetNtrack(); ipart++)
		{
			TParticle *particle = stack->Particle(ipart);
			//if(!particle) continue;
			
			int pdgcode = particle->GetPdgCode();
			
			gHistPDG->Fill(pdgcode);

			
			
			double mass;
			int motherlabel = particle->GetMother(0);
			bool isprompt = (motherlabel<=0);
			//cout << ipart << "\t" << pdgcode << "\t" << motherlabel << endl;
			
			//cout << pdgcode << "\t" << particle->Px() << "\t" << particle->Py() << "\t" << particle->Pz() << endl; 
			
			
			if(ipart==1) {
            	cout << "PDG1 = " << particle->GetPdgCode() << endl;
            	if(TMath::Abs(particle->GetPdgCode()) == 11) {
					p[electronCount][0] = particle->Px();
					p[electronCount][1] = particle->Py();
					p[electronCount][2] = particle->Pz();
					electronCount++;
					cout << "px,y,z = " << p[0][0] << " / " << p[0][1] << " / " << p[0][2] << endl;
				}
			}

			if(ipart==2) {
				cout << "PDG2 = " << particle->GetPdgCode() << endl;
				if(TMath::Abs(particle->GetPdgCode()) == 11) {
					p[electronCount][0] = particle->Px();
					p[electronCount][1] = particle->Py();
					p[electronCount][2] = particle->Pz();
					electronCount++;
					cout << "px,y,z = " << p[1][0] << " / " << p[1][1] << " / " << p[1][2] << endl;
				}
			}
			
			if(pdgcode == 443)
			{
				gHistJpsiAll_pt->Fill(particle->Pt());
				gHistJpsiAll_Y->Fill(particle->Y());
				gHistJpsiAll_eta->Fill(particle->Eta());
				//cout << particle->Eta()<<endl;
				
				if(isprompt) gHistJpsiPrompt_pt->Fill(particle->Pt());	
			}
		}//===END LOOP OVER PARTICLE=====================	
		
		if(electronCount==2)
			{
			double me = 0.511e-3;
				
			double pt1 = TMath::Sqrt(p[0][0]*p[0][0] + p[0][1]*p[0][1]);
			double p1 = TMath::Sqrt(pt1*pt1 + p[0][2]*p[0][2]);
			double pt2 = TMath::Sqrt(p[1][0]*p[1][0] + p[1][1]*p[1][1]);
			double p2 = TMath::Sqrt(pt2*pt2 + p[1][2]*p[1][2]);
			//double mass = 2.0*me*me + 2.0*TMath::Sqrt(me*me + p1*p1)*TMath::Sqrt(me*me + p2*p2) - p[0][0]*p[1][0] - p[0][1]*p[1][1] - p[0][2]*p[1][2];
			double pt = TMath::Sqrt((p[0][0]+p[1][0])*(p[0][0]+p[1][0]) + (p[0][1]+p[1][1])*(p[0][1]+p[1][1]));
			
			double mass= TMath::Sqrt((p1+p2)*(p1+p2) - (p[0][0]+p[1][0])*(p[0][0]+p[1][0]) - (p[0][1]+p[1][1])*(p[0][1]+p[1][1]) - (p[0][2]+p[1][2])*(p[0][2]+p[1][2]));
			
			gHistGammaGammaMass->Fill(mass);
			cout << "INVARIANT MASS = " << mass << endl;
	
			}
			//cout<<"Number of Tracks = " << stack->GetNtrack()<<endl;
		}//======END LOOP OVER EVENT==============
	
	//cout << "Number of events = " << rl->GetNumberOfEvents()<<endl;
	
//===================SAVE ON ROOT FILE =========================
	save->cd();
	
	gHistPDG->Write();
	gHistJpsiAll_pt->Write();
	gHistJpsiAll_Y->Write();
	gHistJpsiAll_eta->Write();
	gHistJpsiPrompt_pt->Write();
	gHistGammaGammaMass->Write();
	
	save->Close();
}



