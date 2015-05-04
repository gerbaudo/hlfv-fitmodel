#include "map"
#include "TString.h"
class RooFitResult;
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class map<TString,RooFitResult*,less<TString>,allocator<pair<const TString,RooFitResult*> > >+;
#pragma link C++ class map<TString,RooFitResult*,less<TString>,allocator<pair<const TString,RooFitResult*> > >::*;
#pragma link C++ operators map<TString,RooFitResult*,less<TString>,allocator<pair<const TString,RooFitResult*> > >::iterator;
#pragma link C++ operators map<TString,RooFitResult*,less<TString>,allocator<pair<const TString,RooFitResult*> > >::const_iterator;
#pragma link C++ operators map<TString,RooFitResult*,less<TString>,allocator<pair<const TString,RooFitResult*> > >::reverse_iterator;
#endif
