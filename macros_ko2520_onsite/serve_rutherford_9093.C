#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "THttpServer.h"
#include "TKey.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TROOT.h"
#include "TError.h"

#include <fstream>
#include <cctype>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {
THttpServer* gRutherfordServer = nullptr;
std::set<std::string> gLoadedFiles;
std::vector<TObject*> gPublishedObjects;
std::vector<TFile*> gOpenFiles;

TString safeName(TString name)
{
  if (name.EndsWith(".root")) name.Remove(name.Length()-5);
  for (int i=0;i<name.Length();i++) {
    const unsigned char c=name[i];
    if (!std::isalnum(c) && c!='_' && c!='-') name[i]='_';
  }
  return name;
}

int publishDirectory(TDirectory* directory, const TString& webFolder)
{
  if (!directory) return 0;
  std::map<std::string,TKey*> newest;
  TIter nextKey(directory->GetListOfKeys());
  while (auto key=dynamic_cast<TKey*>(nextKey())) {
    const std::string name=key->GetName();
    auto found=newest.find(name);
    if (found==newest.end() || key->GetCycle()>found->second->GetCycle())
      newest[name]=key;
  }

  int published=0;
  for (const auto& entry:newest) {
    TKey* key=entry.second;
    TClass* klass=TClass::GetClass(key->GetClassName());
    if (!klass) continue;
    if (klass->InheritsFrom(TDirectory::Class())) {
      published += publishDirectory(
          directory->GetDirectory(key->GetName()),
          webFolder+"/"+safeName(key->GetName()));
      continue;
    }
    if (!klass->InheritsFrom(TH1::Class())
        && !klass->InheritsFrom(TCanvas::Class())) continue;
    TObject* object=key->ReadObj();
    if (!object) continue;
    if (auto histogram=dynamic_cast<TH1*>(object))
      histogram->SetDirectory(nullptr);
    gRutherfordServer->Register(webFolder,object);
    gPublishedObjects.push_back(object);
    published++;
  }
  return published;
}

void loadRecoFiles(const TString& directory)
{
  TSystemDirectory listing("rutherford_js_reco_files",directory);
  TList* entries=listing.GetListOfFiles();
  if (!entries) return;

  TIter next(entries);
  while (auto item=dynamic_cast<TSystemFile*>(next())) {
    if (item->IsDirectory()) continue;
    TString name=item->GetName();
    if (!name.BeginsWith("compass_0") || !name.EndsWith(".reco.root")) continue;
    TString full=directory+"/"+name;
    if (gLoadedFiles.count(full.Data())) continue;

    TFile* input=TFile::Open(full,"READ");
    if (!input || input->IsZombie()) continue;
    const TString folder="/Directory/"+safeName(name);
    const int count=publishDirectory(input,folder);
    // Keep legacy ROOT files open: closing them destroys nested canvas
    // primitives with inconsistent ownership flags and corrupts the server.
    gROOT->GetListOfFiles()->Remove(input);
    gOpenFiles.push_back(input);
    gLoadedFiles.insert(full.Data());
    Info("serve_rutherford_9093","loaded %s (%d drawable object(s))",
         full.Data(),count);
  }
  delete entries;
}

void loadFitFiles(const TString& directory)
{
  TSystemDirectory listing("rutherford_js_fit_files",directory);
  TList* entries=listing.GetListOfFiles();
  if (!entries) return;
  TIter next(entries);
  while (auto item=dynamic_cast<TSystemFile*>(next())) {
    if (item->IsDirectory()) continue;
    TString name=item->GetName();
    const bool compassFit=name.BeginsWith("compass_0") && name.EndsWith(".rutherford.root");
    const bool productionFit=name=="production_beam_rate.rutherford.root";
    if (!compassFit && !productionFit) continue;
    TString full=directory+"/"+name;
    if (gLoadedFiles.count(full.Data())) continue;

    TFile* input=TFile::Open(full,"READ");
    if (!input || input->IsZombie()) continue;
    TString folder;
    if (productionFit) {
      folder="/Directory/production_runs";
    } else {
      TString runFolder=name;
      runFolder.ReplaceAll(".rutherford.root",".reco.root");
      folder="/Directory/"+safeName(runFolder);
    }
    const int count=publishDirectory(input,folder);
    gROOT->GetListOfFiles()->Remove(input);
    gOpenFiles.push_back(input);
    gLoadedFiles.insert(full.Data());
    Info("serve_rutherford_9093","loaded %s (%d Rutherford canvas(es))",
         full.Data(),count);
  }
  delete entries;
}
}

void serve_rutherford_9093()
{
  gROOT->SetBatch(kTRUE);
  // The build contains a consolidated libLILAK.so although its ROOT maps
  // still name the old per-dictionary libraries.  Load the actual library
  // before any reconstructed file is opened.
  if (gSystem->Load("/root/lilak/build/libLILAK.so")<0) {
    Error("serve_rutherford_9093","cannot load libLILAK.so");
    return;
  }
  const char* address=gSystem->Getenv("RUTHERFORD_JS_ADDRESS");
  const char* directory=gSystem->Getenv("RUTHERFORD_JS_DIRECTORY");
  const char* fitDirectory=gSystem->Getenv("RUTHERFORD_JS_FIT_DIRECTORY");
  const char* pidFile=gSystem->Getenv("RUTHERFORD_JS_PIDFILE");
  if (!address) address="192.168.1.102:9093";
  if (!directory) directory="data_reco3";
  if (!fitDirectory) fitDirectory="data_reco3/rutherford_sidecars";

  TString engine=Form("http:%s?top=LILAK;noglobal",address);
  gRutherfordServer=new THttpServer(engine);
  if (!gRutherfordServer->IsAnyEngine()) {
    Error("serve_rutherford_9093","cannot listen at %s",address);
    return;
  }
  if (pidFile) {
    std::ofstream output(pidFile);
    output << gSystem->GetPid() << '\n';
  }

  Info("serve_rutherford_9093","serving %s at http://%s",directory,address);
  while (true) {
    loadRecoFiles(directory);
    loadFitFiles(fitDirectory);
    for (int i=0;i<100;i++) {
      gSystem->ProcessEvents();
      gSystem->Sleep(10);
    }
  }
}
