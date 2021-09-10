///Script for assesing B-physics CMS models. Limited daughters + other FS particles are recorded. Only daughters that pass L1T and other FS particles in at least one daughter cone.
#include "Pythia8/Pythia.h"
//#include <vector>
#include <cmath>
using namespace Pythia8;
using namespace std;


bool is_in(vector<int> v, int iParticle){
	bool r_val = false;
	for(int i = 0; i < int(v.size()); i++){
		if(v[i] == iParticle) r_val = true;
	}
	return r_val;
}

double deltaR(Event event, int i, int j){
	double dPhi = event[i].phi() - event[j].phi();
	double dEta = event[i].eta() - event[j].eta();
	double dR = sqrt(dPhi*dPhi + dEta*dEta);
	return dR;
}

bool pass_trigger(Event event, int iParticle){
	//Get daughter ids
	int d1 = event[iParticle].daughter1();
	int d2 = event[iParticle].daughter2();
	
	//Compute trigger kinematic variables
	double dR = deltaR(event, d1, d2);
	double pT1 = event[d1].pT();
	double pT2 = event[d2].pT();
	double eta1 = event[d1].eta();
	double eta2 = event[d2].eta();

	//Returns true if muon pair passes the trigger requirements
	bool check1 = false;
	if(pT1 > 4. && pT2 > 4. && dR < 1.2) check1 = true;
	if(abs(eta1) < 1.4 && abs(eta2) < 1.4 && dR < 1.4) check1 = true;
	if(pT1 > 15. and pT2 > 7.) check1 = true;
	if(pT1 > 7. and pT2 > 15.) check1 = true;
	
	if(check1){
		if(pT1 > 3. && pT2 > 3. && abs(eta1) < 2.4 && abs(eta2) < 2.4) return true;
	}
	return false;
}

bool is_near_dter(Event event, vector<int> dters, int iTrack){
	double R = 0.3;
	for(int i = 0; i < int(dters.size()); i++){
		if(deltaR(event, dters[i], iTrack) < R) return true;
	}

	return false;
}

void write_jet(int iEvent, int iJet, SlowJet slowJet, ofstream &myFile){
	//Jets are written just like particles with 0. values for properties only relevant 
	//To particles. Pseudorapidity is replaced with rapidity
	myFile  <<  3 << " " << iEvent << " ";
	myFile  <<  iJet << " " << "0 ";
	myFile << scientific << setprecision(16);
	myFile << slowJet.p(iJet).px() << " " << slowJet.p(iJet).py() << " " << slowJet.p(iJet).pz() << " ";
	myFile << slowJet.pT(iJet) << " " << slowJet.y(iJet) << " " << slowJet.phi(iJet) << " ";
	myFile << slowJet.p(iJet).e() << " " << 0. << " " << 0. << " ";
	myFile << 0. << " " << 0. << " " << 0. << " " << 0. << " ";
	myFile << 0. << " " << 0. << " " << 0. << " " << 0. << " ";
	myFile << 0. << " " << 0. << " " << 0. << " " << 0. << " ";
	myFile << "\n";
}

void write_particle(int isAp, int iEvent, int iParticle, Event event, ofstream &myFile){
	//Record all standard info
	if(isAp == 0){
		myFile  << 0 << " " << iEvent << " ";
	}
	else if(isAp == 1){
		myFile << 1 << " " << iEvent << " ";
	}
	else{
		myFile << 2 << " " << iEvent << " ";
	}
	myFile  <<  iParticle << " " << event[iParticle].id() <<" ";
	myFile << scientific << setprecision(16);
	myFile << event[iParticle].px() << " " << event[iParticle].py() << " " << event[iParticle].pz() << " ";
	myFile << event[iParticle].pT() << " " << event[iParticle].eta() << " " << event[iParticle].phi() << " ";
	myFile << event[iParticle].e() << " " << event[iParticle].m() << " " << event[iParticle].tau() << " ";
	myFile << event[iParticle].tProd() << " " << event[iParticle].xProd() << " " << event[iParticle].yProd() << " " << event[iParticle].zProd() << " ";
	myFile << event[iParticle].tDec() << " " << event[iParticle].xDec() << " " << event[iParticle].yDec() << " " << event[iParticle].zDec() << " ";
	myFile << "\n";
}

int main(int argc, char* argv[]) {
  //Open file + get cardName 
  string cardName, outFile;
  cardName=argv[1];
  outFile="events_"+ cardName.substr(0,cardName.find(".",0)) +"_mm12.dat";
  ofstream myFile;
  myFile.open(outFile);
  
  std::cout << outFile <<"\n";
  
  //Instantiate pythia w/ cardName
  Pythia pythia;
  Event& event = pythia.event; // Generator; shorthand for event and particleData.
  pythia.readFile(cardName);
  //pythia.readString("Random:setSeed = on");
  //pythia.readString("Random:seed = 0");
  pythia.init();
  
  // Set up SlowJet jet finder, with anti-kT clustering(-1), R=0.4, pT=30, eta=3.0, observable 
  // FS particles analyzed, all particles given correct masses
  // SlowJet slowJet( -1, 0.4, 30., 3.0, 2, 2); This is correct - just using smaller pT to make sure I'm getting jets. 
  SlowJet slowJet( -1, 0.4, 30., 3.0, 2, 2);
  
  //int nEvents = pythia.mode("Main:numberOfEvents");
  int nEvents = 10;
  // Generate events
  int iEvent = 0;
  while(iEvent < nEvents){
    if (!pythia.next()) continue;
    bool foundOne = false;
    vector<int> dters; 
    for(int i = 0; i < event.size(); i++){
	//If found Ap
    	if(abs(event[i].id()) == 999999){
		if(pass_trigger(event, i)){
		  write_particle(0, iEvent, i, event, myFile);
		  write_particle(1, iEvent, event[i].daughter1(), event, myFile);
		  write_particle(1, iEvent, event[i].daughter2(), event, myFile);
		  dters.push_back(event[i].daughter1());
		  dters.push_back(event[i].daughter2());
		  foundOne = true;
		}
	}
    }

    if(foundOne){
      for(int i = 0; i < event.size(); i++){
	if(event[i].isFinal() && event[i].isCharged() && !is_in(dters, i)){
		if(is_near_dter(event, dters, i)){	
			write_particle(2, iEvent, i, event, myFile);
		}
	}
      }
    

    // Save all SlowJet jet data
      slowJet.analyze(event);

      for (int i = 0; i < slowJet.sizeJet(); i++) {
        write_jet(iEvent, i, slowJet, myFile);
      }

      iEvent = iEvent+1;
      cout << "iEvent: " << iEvent << endl;
    }
    
  }
  pythia.stat();
  myFile.close();
  // Done.
  return 0;
}
