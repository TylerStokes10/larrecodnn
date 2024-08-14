#include "classes.h"
#include "TClass.h"
#include "TROOT.h"
#include <iostream>

void checksum() {
    TClass *fvClass = TClass::GetClass("anab::FeatureVector<8>");
    TClass *mvaClass = TClass::GetClass("anab::MVADescription<8>");

    if (fvClass) {
        std::cout << "Checksum for anab::FeatureVector<8>: " << fvClass->GetCheckSum() << std::endl;
    } else {
        std::cerr << "Could not find class anab::FeatureVector<8>" << std::endl;
    }

    if (mvaClass) {
        std::cout << "Checksum for anab::MVADescription<8>: " << mvaClass->GetCheckSum() << std::endl;
    } else {
        std::cerr << "Could not find class anab::MVADescription<8>" << std::endl;
    }
}

