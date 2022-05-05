#include "pengrobinson.h"

#include "Math/poly34.h"

#include "pugixml.hpp"

PengRobinson::PengRobinson() = default;

void PengRobinson::init(const std::vector<std::string>& speciesNames, const std::vector<double> binCoeff) {
  this->gas_constant_universal = 8314.4621; // [J kmol^-1 K^-1]
  this->N_Av = 6.022e26; // kmol^-1
  this->Boltzman = 1.380649e-23;  // J/K

  this->n_species = speciesNames.size();

  // array initialization
  this->species.resize(this->n_species);

  this->X.resize(this->n_species);
  this->Y.resize(this->n_species);
  this->MW.resize(this->n_species);

  this->Tcrit.resize(this->n_species);
  this->Pcrit.resize(this->n_species);
  this->rhocrit.resize(this->n_species);
  this->Vcrit.resize(this->n_species);
  this->Zcrit.resize(this->n_species);
  this->omega.resize(this->n_species);

  this->NasaCoef = 7;
  this->nasa_poly_coeff.resize(this->NasaCoef * this->n_species * 2);
  this->nasa_poly_bounds.resize(this->n_species * 2);
  this->hspecies_ig.resize(this->n_species);
  this->cpspecies_ig.resize(this->n_species);
  this->sspecies_ig.resize(this->n_species);

  this->cst_b.resize(this->n_species);
  this->cst_a.resize(this->n_species);
  this->cst_c.resize(this->n_species);

  this->A_I.resize(this->n_species);
  this->dA_IdT.resize(this->n_species);
  this->d2A_IdT2.resize(this->n_species);
  this->A_IJ.resize(this->n_species * this->n_species);
  this->dAmdN.resize(this->n_species);
  this->d2AmdTdN.resize(this->n_species);
  this->dK1dN.resize(this->n_species);
  this->dPdN.resize(this->n_species);
  this->dVdN.resize(this->n_species);

  this->hspecies.resize(this->n_species);
  this->sspecies.resize(this->n_species);
  this->muspecies.resize(this->n_species);
  this->edspecies.resize(this->n_species);

  this->binCoeff = binCoeff;

  // init constants
  LOOP_l_N(this->n_species)
    this->species[l] = speciesNames[l];

  this->ReadNasaPolynomials();
  this->ReadCriticalProperties();
  this->SetRealFluidConstants();
}

//----------------------------------------------------------------------------

PengRobinson::~PengRobinson() = default;

//----------------------------------------------------------------------------

void PengRobinson::ReadNasaPolynomials() {
  std::string nasa_poly_dbpath = NASAPOLYPATH;

  std::vector<double> readFromXMLFile_low;
  std::vector<double> readFromXMLFile_high;

  pugi::xml_document xmlDoc;
  pugi::xml_parse_result result = xmlDoc.load_file(nasa_poly_dbpath.c_str());

  LOOP_l_N(this->n_species) {
    // search for species
    pugi::xml_node
        species = xmlDoc.child("ctml").child("speciesData").first_child();
    while (species.attribute("name").value() != this->species[l]) {
      species = species.next_sibling();
    }

    pugi::xml_node lowTemp = species.child("thermo").first_child();
    pugi::xml_node highTemp = lowTemp.next_sibling();

    // read temperature ranges
    this->nasa_poly_bounds[0 + 2 * l].first =
        lowTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[0 + 2 * l].second =
        highTemp.attribute("Tmin").as_double();
    this->nasa_poly_bounds[1 + 2 * l].first =
        lowTemp.attribute("Tmax").as_double();
    this->nasa_poly_bounds[1 + 2 * l].second =
        highTemp.attribute("Tmax").as_double();

    readFromXMLFile_low =
        Tokenize(lowTemp.child("floatArray").child_value(), " ,\n");
    readFromXMLFile_high =
        Tokenize(highTemp.child("floatArray").child_value(), " ,\n");

    LOOP_k_N(this->NasaCoef) {
      this->nasa_poly_coeff[k + this->NasaCoef * l] = readFromXMLFile_low[k];
      this->nasa_poly_coeff[k + this->NasaCoef * l
          + this->NasaCoef * this->n_species] = readFromXMLFile_high[k];
    }
  }
}

//----------------------------------------------------------------------------

void PengRobinson::ReadCriticalProperties() {
  //! \brief Reads the XML data-base and sets the coefficients
  LOOP_l_N(this->n_species) {
    if (this->species[l] == "O2") {
      this->MW[l] = 31.9988; //O2
      this->Tcrit[l] = 154.5800;
      this->Pcrit[l] = 5.0430e+6;
      this->rhocrit[l] = 436.140;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.0222;
    } else if (this->species[l] == "H2") {
      this->MW[l] = 2.01588; //H2
      this->Tcrit[l] = 33.1450;
      this->Pcrit[l] = 1.2964e+6;
      this->rhocrit[l] = 31.262;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = -0.219;
    } else if (this->species[l] == "H2O") {
      this->MW[l] = 18.01528;
      this->Tcrit[l] = 647.096;
      this->Pcrit[l] = 22.0640e+6;
      this->rhocrit[l] = 322.0;
      this->Vcrit[l] = MW[l] / rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.3443;
    } else if (this->species[l] == "OH") {
      this->MW[l] = 17.0073; //OH kg/kmol
      this->Tcrit[l] = 0.0;
      this->Pcrit[l] = 0.0;
      this->rhocrit[l] = 0.0;
      this->Vcrit[l] = 0.0;
      this->Zcrit[l] = 0.0;
      this->omega[l] = 0.0;
    } else if (this->species[l] == "N2") {
      this->MW[l] = 28.0134; //N2 kg/kmol
      this->Tcrit[l] = 126.1900;
      this->Pcrit[l] = 3.3958e+6;
      this->rhocrit[l] = 313.3;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.03720;
    } else if (this->species[l] == "NH3") {
      this->MW[l] = 17.031; //NH3 kg/kmol
      this->Tcrit[l] = 405.40;
      this->Pcrit[l] = 11.3330e+6;
      this->rhocrit[l] = 225.000;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.25601;
    } else if (this->species[l] == "C8H18,isooctane") {
      this->MW[l] = 114.23;
      this->Tcrit[l] = 543.9;
      this->Pcrit[l] = 25.7e5;
      this->rhocrit[l] = 244.4522;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.394;
    } else if (this->species[l] == "NC12H26") {
      this->MW[l] = 170.33484;
      this->Tcrit[l] = 658.10;
      this->Pcrit[l] = 1.817e+6;
      this->rhocrit[l] = 226.5453372;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.574;
    } else if (this->species[l] == "CH4") {
      this->MW[l] = 16.04;
      this->Tcrit[l] = 190.6;
      this->Pcrit[l] = 46.1e5;
      this->rhocrit[l] = 162.0;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.011;
    } else if (this->species[l] == "C2H6") {
      this->MW[l] = 30.0690;
      this->Tcrit[l] = 305.3;
      this->Pcrit[l] = 49e5;
      this->rhocrit[l] = 6.9*30.0690;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.099;
    } else if (this->species[l] == "C3H8") {
      this->MW[l] = 44.0956;
      this->Tcrit[l] = 369.9;
      this->Pcrit[l] = 42.5e5;
      this->rhocrit[l] = 5.1*44.0956;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.153;
    } else if (this->species[l] == "C4H10,n-butane") {
      this->MW[l] = 58.1222;
      this->Tcrit[l] = 425.;
      this->Pcrit[l] = 38.0e5;
      this->rhocrit[l] = 3.92*58.1222;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.199;
    } else if (this->species[l] == "C5H12,n-pentane") {
      this->MW[l] = 72.1488;
      this->Tcrit[l] = 469.8;
      this->Pcrit[l] = 33.6e5;
      this->rhocrit[l] = 3.22*72.1488;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.251;
    } else if (this->species[l] == "C6H14,n-hexane") {
      this->MW[l] = 86.1754;
      this->Tcrit[l] = 507.6;
      this->Pcrit[l] = 30.2e5;
      this->rhocrit[l] = 2.71*86.1754;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.299;
    } else if (this->species[l] == "C7H16,n-heptane") {
      this->MW[l] = 100.2019;
      this->Tcrit[l] = 540;
      this->Pcrit[l] = 27.4e5;
      this->rhocrit[l] = 2.35*100.2019;
      this->Vcrit[l] = this->MW[l] / this->rhocrit[l];
      this->Zcrit[l] =
          (Pcrit[l] * Vcrit[l]) / (this->gas_constant_universal * Tcrit[l]);
      this->omega[l] = 0.349;
    } else {
      std::cout << " WARNING -> Unknown species :[" << this->species[l]
           << "]. No critical properties found." << std::endl;
    }
  }

  return;
}

//----------------------------------------------------------------------------

void PengRobinson::SetRealFluidConstants() {
  LOOP_k_N(this->n_species) {
    this->cst_b[k] = 0.077796 * this->gas_constant_universal * this->Tcrit[k] / this->Pcrit[k];

    this->cst_a[k] = 0.457236 * std::pow(this->gas_constant_universal * this->Tcrit[k], 2)  / this->Pcrit[k];

    if (this->omega[k] <= 0.49)
      this->cst_c[k] = 0.37464 + 1.54226 * this->omega[k] - 0.26992 * std::pow(this->omega[k], 2);
    else
      this->cst_c[k] = 0.379642 + 1.485030*this->omega[k] - 0.164423*std::pow(this->omega[k], 2) + 0.016666*std::pow(this->omega[k], 3);
  }
}

//----------------------------------------------------------------------------

double PengRobinson::GetRho_SetMixture_TPY(const double T, const double P, const double *Y) {
  this->T = T;
  this->P = P;
  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncRZFromPressureTemperature();
  return this->rho;
}

double PengRobinson::GetRho_SetMixture_TPX(const double T, const double P, const double *X) {
  this->T = T;
  this->P = P;
  this->SetMolarFractionFromX(X);

  this->SetMolecularWeightMixtureFromX();
  this->SetMassFractionFromX();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncRZFromPressureTemperature();
  return this->rho;
}

//----------------------------------------------------------------------------

void PengRobinson::GetTdew_PXvapor(const double P, const double *Xvapor, const double TGuess, const double *XliqGuess, double &Tdew_out, double *Xliq_out) {
  assert(this->n_species == 2);
  // First iter
  double T = TGuess;
  double Xliq [2] = {XliqGuess[0], XliqGuess[1]};

  double phil[2], phiv[2];
  this->getPhifromTandPandX(T, P, Xliq, phil);
  this->getPhifromTandPandX(T, P, Xvapor, phiv);
  double K[2] = {phil[0] / phiv[0], phil[1] / phiv[1]};
  double xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];

  double xT_oldold = xT;
  double T_oldold = T;

  // Second iter
  T = T * xT;
  Xliq[0] = Xvapor[0] / K[0] / xT; Xliq[1] = 1. - Xliq[0];

  this->getPhifromTandPandX(T, P, Xliq, phil);
  this->getPhifromTandPandX(T, P, Xvapor, phiv);
  K[0] = phil[0] / phiv[0]; K[1] = phil[1] / phiv[1];
  xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];

  double xT_old = xT;
  double T_old = T;

  // If we need more iter
  while (std::fabs(xT-1) > 1e-12) {
    T = T_oldold + (1-xT_oldold)/(xT_old-xT_oldold)*(T_old-T_oldold);
    Xliq[0] = Xvapor[0] / K[0] / xT; Xliq[1] = 1. - Xliq[0];

    this->getPhifromTandPandX(T, P, Xliq, phil);
    this->getPhifromTandPandX(T, P, Xvapor, phiv);
    K[0] = phil[0] / phiv[0]; K[1] = phil[1] / phiv[1];
    xT = Xvapor[0] / K[0] + (1.-Xvapor[0]) / K[1];

    xT_oldold = xT_old;
    T_oldold = T_old;
    xT_old = xT;
    T_old = T;
  }
  Tdew_out = T;
  Xliq_out[0] = Xliq[0]; Xliq_out[1] = Xliq[1];
}

void PengRobinson::getPhifromTandPandX(const double T_in, const double P_in, const double *X_in, double *phiOut) {
  double B = 0.;
  LOOP_k_N(this->n_species)
      B += this->cst_b[k] * X_in[k];

  double A_I [this->n_species];
  LOOP_k_N(this->n_species) {
    A_I[k] = this->cst_a[k] * std::pow(1.0 + this->cst_c[k] * (1.0 - std::sqrt(T_in / this->Tcrit[k])), 2);
  }

  double A = 0.;
  double A_IJ [this->n_species * this->n_species];
  double dAdN [this->n_species];
  LOOP_k_N(this->n_species) {
    dAdN[k] = 0.;
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = X_in[l] * X_in[k];
      A_IJ[apos] = (1.-this->binCoeff[apos])*std::sqrt(A_I[k]*A_I[l]);
      A += X_X * A_IJ[apos];
      dAdN[k] += X_in[l] * A_IJ[apos];
    }
    dAdN[k] *= 2.;
  }

  double Amix = A * P_in / std::pow(this->gas_constant_universal, 2) / std::pow(T_in, 2);
  double Bmix = B * P_in / this->gas_constant_universal / T_in;

  double a0 = -(Amix * Bmix - Bmix * Bmix - Bmix * Bmix * Bmix);
  double a1 = Amix - 3 * Bmix * Bmix - 2 * Bmix;
  double a2 = -(1. - Bmix);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > Bmix)
      Z.push_back(xZ[i]);
  double Zout;
  if (Z.size() == 1) {
    Zout = Z[0];
  } else {
    std::vector<double> lnPhi(Z.size());

    int iMin = 0;
    for (int i = 0; i < Z.size(); i++) {
      lnPhi[i] = -std::log(Z[i] - Bmix)
          - Amix / Bmix / std::sqrt(8) * std::log((Z[i] + (1 + std::sqrt(2)) * Bmix) / (Z[i] + (1 - std::sqrt(2)) * Bmix))
          + Z[i] - 1;
      if (lnPhi[i] < lnPhi[iMin]) iMin = i;
    }
    Zout = Z[iMin];
  }

  LOOP_k_N(this->n_species) {
    double Bi = this->cst_b[k] * P_in / this->gas_constant_universal / T_in;
    phiOut[k] = std::exp(Bi / Bmix * (Zout - 1.) - std::log(Zout - Bmix)
                         - Amix / Bmix / 2.8284 * std::log((Zout+2.4142*Bmix)/(Zout-.4142*Bmix)) * (dAdN[k] / A - Bi/Bmix));
  }

}


//----------------------------------------------------------------------------

void PengRobinson::SetMixture_TRY(const double T, const double R, const double *Y) {
  this->T = T;
  this->rho = R;
  this->SetMassFractionFromY(Y);

  this->SetMolecularWeightMixtureFromY();
  this->SetMolarFractionFromY();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
}

void PengRobinson::SetMixture_TRX(const double T, const double R, const double *X) {
  this->T = T;
  this->rho = R;
  this->SetMolarFractionFromX(X);

  this->SetMolecularWeightMixtureFromX();
  this->SetMassFractionFromX();

  this->SyncRealFluidThermodynamicsFromTemperatureDensity();
  this->SyncPZFromTemperatureDensity();

  this->SyncIdealFluidThermodynamicsFromTemperature();
  this->SyncEFromTemperatureDensity();

  this->SyncPartialPropertiesFromTemperatureDensity();
}

//----------------------------------------------------------------------------

void PengRobinson::SyncRealFluidThermodynamicsFromTemperatureDensity() {
  double v = this->MW_M / this->rho;

  this->Am = 0.;
  this->Bm = 0.;
  this->dAmdT = 0.;
  this->d2AmdT2 = 0.;


  LOOP_k_N(this->n_species) {
    this->A_I[k] = this->cst_a[k] * std::pow(1.0 + this->cst_c[k] * (1.0 - std::sqrt(this->T / this->Tcrit[k])), 2);

    this->dA_IdT[k] = -this->cst_a[k] * (1.0 + this->cst_c[k] * (1.0 - std::sqrt(this->T / this->Tcrit[k])))
        * this->cst_c[k] / std::sqrt(this->T * this->Tcrit[k]);

    this->d2A_IdT2[k] = this->cst_a[k] * 0.5 * this->cst_c[k] *
        (1+this->cst_c[k]) / std::sqrt(this->Tcrit[k] * std::pow(this->T,3));
  }

  LOOP_k_N(this->n_species) {
    this->dAmdN[k] = 0.;
    this->d2AmdTdN[k] = 0.;

    this->Bm += this->X[k] * this->cst_b[k];

    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->X[l] * this->X[k];
      this->A_IJ[apos] = (1.-this->binCoeff[apos])*std::sqrt(this->A_I[k]*this->A_I[l]);

      this->Am += X_X * this->A_IJ[apos];
      this->dAmdT += X_X * 0.5 / std::sqrt(this->A_I[k]*this->A_I[l])
          * (this->A_I[k]*this->dA_IdT[l] + this->A_I[l]*this->dA_IdT[k]);
      this->d2AmdT2 += X_X *
          (0.5 / std::sqrt(this->A_I[k]*this->A_I[l]) * (2*this->dA_IdT[k]*this->dA_IdT[l] + this->A_I[k]*this->d2A_IdT2[l] + this->A_I[l]*this->d2A_IdT2[k])
           - 0.25 * (this->A_I[k]*this->dA_IdT[l] + this->A_I[l]*this->dA_IdT[k]) * (this->A_I[k]*this->dA_IdT[l] + this->A_I[l]*this->dA_IdT[k]) * std::pow(this->A_I[k]*this->A_I[l], -3./2.));

      this->dAmdN[k] += this->X[l] * this->A_IJ[apos];
      this->d2AmdTdN[k] += this->X[l] * 0.5 / std::sqrt(this->A_I[k]*this->A_I[l])
          * (this->A_I[k]*this->dA_IdT[l] + this->A_I[l]*this->dA_IdT[k]);
    }
    this->dAmdN[k] *= 2.0;
    this->d2AmdTdN[k] *= 2.0;
  }

  this->dPdT = this->gas_constant_universal / (v - this->Bm)
      - this->dAmdT / (std::pow(v, 2) + 2.0 * v * this->Bm - std::pow(this->Bm, 2));
  double arg = this->gas_constant_universal * this->T * (v + this->Bm)
      * std::pow((v / (v - this->Bm) + this->Bm / (v + this->Bm)), 2);
  this->dPdV = -this->gas_constant_universal * this->T / std::pow((v - this->Bm), 2)
      * (1.0 - 2.0 * this->Am / arg);
  this->expansivity = -this->dPdT / (v * this->dPdV); //ideal gas: equal to 3.34E-3 (1/K)
  this->K1 = 1.0 / (std::sqrt(8.0) * this->Bm)
      * std::log((v + (1 - std::sqrt(2.0)) * this->Bm) / (v + (1 + std::sqrt(2.0)) * this->Bm));

  double temp = v * v + 2.0 * this->Bm * v - this->Bm * this->Bm;
  LOOP_k_N(this->n_species) {
    this->dPdN[k] = this->gas_constant_universal * this->T / (v - this->Bm) +
        this->gas_constant_universal * this->T * this->cst_b[k]
            / std::pow((v - this->Bm), 2) - this->dAmdN[k] / temp
        + 2.0 * this->Am * this->cst_b[k] * (v - this->Bm) / std::pow(temp, 2);
    this->dVdN[k] = -this->dPdN[k] / this->dPdV;
    this->dK1dN[k] = 1.0 / temp * this->dVdN[k] - this->cst_b[k] / this->Bm *
        (this->K1 + v / temp);
  }
}

//----------------------------------------------------------------------------
void PengRobinson::SyncIdealFluidThermodynamicsFromTemperature() {
  double T1 = this->T;
  double T2 = this->T * this->T;
  double T3 = T2 * this->T;
  double T4 = T3 * this->T;

  this->h_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1 / 2.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 3.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 4.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 5.0;
      this->hspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 5] / T1;
    } else {
      this->hspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species];
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1 / 2.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2 / 3.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3 / 4.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4 / 5.0;
      this->hspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 5
          + this->NasaCoef * this->n_species] / T1;
    }
    this->h_ig += this->X[l] * (hspecies_ig[l]);
  }
  this->h_ig *= (this->T * this->gas_constant_universal) / this->MW_M;

  this->cp_ig = 0.;
  this->cv_ig = 0.;
  LOOP_l_N(this->n_species) {
    if (this->T <= 1000.0) {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0];
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3;
      this->cpspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4;
    } else {
      this->cpspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species];
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3;
      this->cpspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4;
    }
    double tmp = cpspecies_ig[l] * this->gas_constant_universal;
    this->cp_ig += this->X[l] * tmp; //cp_mol
    this->cv_ig += this->X[l] * (tmp - this->gas_constant_universal);
  }
  this->cp_ig /= this->MW_M;
  this->cv_ig /= this->MW_M;

  LOOP_l_N(this->n_species) {
    if (this->T < 1000.0) {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0] * std::log(T1);
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 1] * T1;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 2] * T2 / 2.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 3] * T3 / 3.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 4] * T4 / 4.0;
      this->sspecies_ig[l] +=
          this->nasa_poly_coeff[l * this->NasaCoef + 6];
    } else {
      this->sspecies_ig[l] = this->nasa_poly_coeff[l * this->NasaCoef + 0
          + this->NasaCoef * this->n_species] * std::log(T);
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 1
          + this->NasaCoef * this->n_species] * T1;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 2
          + this->NasaCoef * this->n_species] * T2 / 2.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 3
          + this->NasaCoef * this->n_species] * T3 / 3.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 4
          + this->NasaCoef * this->n_species] * T4 / 4.0;
      this->sspecies_ig[l] += this->nasa_poly_coeff[l * this->NasaCoef + 6
          + this->NasaCoef * this->n_species];
    }
    this->sspecies_ig[l] -= std::log(this->X[l]*this->P/1e5);
    if (this->X[l] < 1e-10)
      sspecies_ig[l] = 0.;
    if (std::isnan(sspecies_ig[l])) {
      std::cout << "Something wrong Ideal. T = " << this->T << ", rho = " << this->rho
                << ", YF = " << this->Y[0] << " , P = " << this->P << std::endl;
      throw -1;
    }
  }
}

//----------------------------------------------------------------------------

void PengRobinson::SyncRZFromPressureTemperature() {
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  if (Z.size() == 1) {
    this->Z = Z[0];
    this->rho = this->P / this->gas_constant / this->T / this->Z;
    return;
  } else {
    std::vector<double> lnPhi(Z.size());

    int iMin = 0;
    for (int i = 0; i < Z.size(); i++) {
      lnPhi[i] = -std::log(Z[i] - B)
          - A / B / std::sqrt(8) * std::log((Z[i] + (1 + std::sqrt(2)) * B) / (Z[i] + (1 - std::sqrt(2)) * B))
          + Z[i] - 1;
      if (lnPhi[i] < lnPhi[iMin]) iMin = i;
    }
    this->Z =  Z[iMin];
    this->rho = this->P / this->gas_constant / this->T / this->Z;

    return;
  }
}

//----------------------------------------------------------------------------

void PengRobinson::SyncPZFromTemperatureDensity() {
  double v = this->MW_M / this->rho; // m^3 / kmol

  this->P = (this->gas_constant_universal * this->T) / (v - this->Bm)
      - this->Am / (std::pow(v, 2) + 2.0 * v * this->Bm - pow(this->Bm, 2));

  this->Z = this->P * v / this->gas_constant_universal / this->T;

  double Tcmix, Pcmix, wcmix, Zcmix;
  this->ComputeMixtureCrit(Tcmix, Pcmix, wcmix, Zcmix);
  // Test if we might be in vapor dome
  if (this->P < 0 || this->CheckThreeRoots()) {
    double Tcmix, pcmix, wcmix, Zcmix;
    this->ComputeMixtureCrit(Tcmix, pcmix, wcmix, Zcmix);
    double psat = pcmix * pow(10.0, (7 / 3 * (1 + wcmix)) * (1 - Tcmix / this->T));

    double a = this->Am;
    double b = this->Bm;

    double A = a * psat / pow(this->gas_constant_universal, 2) / pow(this->T, 2);
    double B = b * psat / this->gas_constant_universal / this->T;

    double a0 = -(A * B - B * B - B * B * B);
    double a1 = A - 3 * B * B - 2 * B;
    double a2 = -(1 - B);

    std::vector<double> Z(3);
    int n = SolveP3(&Z[0], a2, a1, a0);
    if (n == 1)
      this->P = psat; // short-cut estimate wrong with PR
    else {
      double rhoL = psat / this->gas_constant / this->T / Z[0];
      double rhoV = psat / this->gas_constant / this->T / Z[2];

      if (this->P < 0 || (this->rho < rhoL && this->rho > rhoV))
        this->P = psat;
    }
  }

  if (this->P < 0 || std::isnan(this->P)) {
    std::cout << "Invalid pressure" << std::endl;
    throw -1;
  }

}

//----------------------------------------------------------------------------

void PengRobinson::SyncEFromTemperatureDensity() {
  double dep = (this->Am - this->T * this->dAmdT) * this->K1 / this->MW_M;
  this->E = this->h_ig - this->gas_constant * this->T + dep;

  double departureCv = (-this->T * this->d2AmdT2 * this->K1) / this->MW_M;
  double departureCp = departureCv
      + (-this->T * std::pow(this->dPdT, 2) / this->dPdV - this->gas_constant_universal)
          / this->MW_M;
  this->cv = this->cv_ig + departureCv;
  this->cp = this->cp_ig + departureCp;
  this->gamma = this->cp / this->cv;
  double v = this->MW_M / this->rho;
  double isocompressibility = -1. / (v * this->dPdV);
  this->sos = std::sqrt(this->gamma / (this->rho * std::fabs(isocompressibility)));
}

//----------------------------------------------------------------------------

void PengRobinson::SyncPartialPropertiesFromTemperatureDensity() {
  // hspecies_ig is molar, without rt
  double temp = this->Am - this->T * this->dAmdT;
  double rt = this->gas_constant_universal * this->T;

  LOOP_k_N(this->n_species) {
    this->edspecies[k] = this->dK1dN[k] * temp +
        this->K1 * (this->dAmdN[k] - this->T * this->d2AmdTdN[k]);
    this->hspecies[k] = this->hspecies_ig[k] * rt - rt + this->edspecies[k] +
        this->P * dVdN[k];
    this->hspecies[k] /= this->MW[k];
    this->edspecies[k] /= this->MW[k];
  }

  // muspecies, sspecies
  double temp1;
  LOOP_k_N(this->n_species) {
    double logphik = this->cst_b[k] / this->Bm * (this->Z - 1.) - std::log(this->Z - this->Bm*this->P/rt)
                              - this->Am / this->Bm / rt / 2.8284 * std::log((this->Z+2.4142*this->Bm*this->P/rt)/(this->Z-.4142*this->Bm*this->P/rt)) * (this->dAmdN[k]/this->Am - this->cst_b[k]/this->Bm);

    this->muspecies[k] = this->hspecies_ig[k]*rt - this->sspecies_ig[k]*rt + rt * logphik;
    this->muspecies[k] /= this->MW[k];
    this->sspecies[k] = (this->hspecies[k] - this->muspecies[k]) / this->T;
  }
}

//----------------------------------------------------------------------------

bool PengRobinson::CheckThreeRoots() {
  double a = this->Am;
  double b = this->Bm;

  double A = a * this->P / std::pow(this->gas_constant_universal, 2) / std::pow(this->T, 2);
  double B = b * this->P / this->gas_constant_universal / this->T;

  double a0 = -(A * B - B * B - B * B * B);
  double a1 = A - 3 * B * B - 2 * B;
  double a2 = -(1 - B);

  std::vector<double> xZ(3);
  double n = SolveP3(&xZ[0], a2, a1, a0);
  std::vector<double> Z;
  for (int i = 0; i < n; i++)
    if (xZ[i] > B)
      Z.push_back(xZ[i]);

  return Z.size() > 1;
}

void PengRobinson::ComputeMixtureCrit(double &Tcmix, double &Pcmix, double &wcmix, double &Zcmix) {
  Tcmix = 0.;
  Pcmix = 0.;
  wcmix = 0.;
  Zcmix = 0.;
  double Vcmix = 0.;

  LOOP_k_N(this->n_species) {
    Tcmix += this->X[k] * this->Tcrit[k];
    Zcmix += this->X[k] * this->Zcrit[k];
    Vcmix += this->X[k] * this->Vcrit[k];
    wcmix += this->X[k] * this->omega[k];
  }
  Pcmix = Zcmix * this->gas_constant_universal * Tcmix / Vcmix;
}

//----------------------------------------------------------------------------

void PengRobinson::getCP(double *X, double* kout, double* lout, double kappaGuess_in, double lambdaGuess_in) {
  this->SetMolarFractionFromX(X);

  this->SetMolecularWeightMixtureFromX();
  this->SetMassFractionFromX();

  // get Bm
  this->Bm = 0;
  LOOP_k_N(this->n_species)
    this->Bm += this->X[k] * this->cst_b[k];

  // Run iteration loop here
  double Tcmix = 0;
  LOOP_k_N(this->n_species) {
    Tcmix += this->X[k] * this->Tcrit[k];
  }

  double kappaGuess = 3.5;
  if (kappaGuess_in > 0)
    kappaGuess = kappaGuess_in;
  double firstlGuess = std::sqrt(1.3*Tcmix);
  if (lambdaGuess_in > 0)
    firstlGuess = lambdaGuess_in;

  double k, l;
  this->solveCP(kappaGuess, firstlGuess, k, l);

  double D,C;
  this->objectiveFun(l, k, D, C);

  double v = k * this->Bm;
  double rho = this->MW_M / v;

  this->SetMixture_TRX(l*l, rho, X);

  if (kout!=nullptr)
    *kout = k;

  if (lout!=nullptr)
    *lout = l;

}


void PengRobinson::solveCP(double kappaGuess, double firstlGuess, double& kout, double& lout) {
  double k, nextk, error, nextError, dummy;
  k = kappaGuess;
  nextk = kappaGuess*1.00001;

  double l = firstlGuess, nextl=firstlGuess;

  do {
    l = this->lambdaFromKappa(k, l);
    this->objectiveFun(l, k, dummy, error);

    nextl = this->lambdaFromKappa(nextk, nextl);
    this->objectiveFun(nextl, nextk, dummy, nextError);

    double temp = k;
    k = nextk;
    nextk = temp + (nextk - temp) * (0. - error) / (nextError - error);
  } while (std::fabs((k-nextk)/k) > this->eps);

  kout = k;
  lout = this->lambdaFromKappa(kout, l);


}


double PengRobinson::lambdaFromKappa(double kappa_in, double lambdaGuess) {
  double l, nextl, error, nextError, dummy;
  l = lambdaGuess;
  nextl = lambdaGuess*1.00001;

  do {
    this->objectiveFun(l, kappa_in, error, dummy);
    this->objectiveFun(nextl, kappa_in, nextError, dummy);
    double temp = l;
    l = nextl;
    nextl = temp + (nextl - temp) * (0. - error) / (nextError - error);
  } while (std::fabs((l-nextl)/l) > this->eps);

  return l;
}


void PengRobinson::objectiveFun(double lambda_in, double kappa_in, double& Dout, double& Cout) {
  this->T = lambda_in*lambda_in;
  double k = kappa_in;

  // get Am
  this->Am = 0.;

  LOOP_k_N(this->n_species) {
    this->A_I[k] = this->cst_a[k] * std::pow(1.0 + this->cst_c[k] * (1.0 - std::sqrt(this->T / this->Tcrit[k])), 2);
  }

  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      double X_X = this->X[l] * this->X[k];
      this->A_IJ[apos] = (1.-this->binCoeff[apos])*std::sqrt(this->A_I[k]*this->A_I[l]);

      this->Am += X_X * this->A_IJ[apos];
    }
  }

  // get F1-8 (eq A, M&H 1981)
  double F1, F2, F3, F4, F5, F6, F7, F8;
  double d1 = 1.+std::sqrt(2.), d2 = 1.-std::sqrt(2.);
//  double k = v / this->Bm;

  double g1 = d1/(k+d1), g2 = d2/(k+d2);

  F1 = 1./(k-1.);
  F2 = 2.*(g1 - g2)/(d1-d2);
  F3 = (g1*g1 - g2*g2)/(d1-d2);
  F4 = (g1*g1*g1 - g2*g2*g2)/(d1-d2);
  F5 = 2.*std::log((k+d1)/(k+d2))/(d1-d2);
  F6 = F2-F5;
  F7 = -F2/(1.+F1);
  F8 = F3/(1.+F1);

  // get alpha, beta (eq 15, M&H 1981)
  Eigen::VectorXd alpha = Eigen::VectorXd::Zero(this->n_species);
  Eigen::VectorXd beta = Eigen::VectorXd::Zero(this->n_species);
  LOOP_k_N(this->n_species) {
    beta[k] = this->cst_b[k]/this->Bm;
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      alpha[k] += this->X[l] * this->A_IJ[apos];
    }
    alpha[k] /= this->Am;
  }

  // get gamma_a, gamma_b (eq 15, B&L 1986)
  Eigen::VectorXd gamma_a = Eigen::VectorXd::Zero(this->n_species);
  Eigen::VectorXd gamma_b = Eigen::VectorXd::Zero(this->n_species);
  Eigen::VectorXd u = Eigen::VectorXd::Ones(this->n_species);

  gamma_a = (this->Am/F5) * (F1*F5-F6)/(1.+F1) * beta;
  gamma_b = (this->gas_constant_universal*this->T*this->Bm*F1/F5)*u
      + (this->Am/F5)*(F3/(1.+F1)+F6)*beta - (this->Am*F6/F5)*alpha;

  // get Matrix U and its inverse (eq 13, B&L 1986)
  Eigen::MatrixXd U = Eigen::MatrixXd::Zero(this->n_species, this->n_species);
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k * this->n_species + l;
      if (k!=l)
        U(k,l) = this->A_IJ[apos];
      else // k==l
        U(k,l) = this->A_IJ[apos]
            - (this->gas_constant_universal*this->T/this->X[k])*(this->Bm/F5);
    }
  }

  Eigen::MatrixXd iU = U.inverse();

  // get D1-4 (eq 18, B&L 1986)
  double D1, D2, D3, D4;

  D1 = (alpha.transpose()*iU*gamma_b)(0,0);
  D2 = (alpha.transpose()*iU*gamma_a)(0,0) - 1.;
  D3 = (beta.transpose()*iU*gamma_b)(0,0) - 1.;
  D4 = (beta.transpose()*iU*gamma_a)(0,0);

  // get alphaBar, betaBar (eq 20 21, B&L 1986)
  double alphaBar, betaBar;
  alphaBar = D1;
  betaBar = -D2;

  // get dN (eq 11, B&L 1986)
  Eigen::VectorXd dN = iU*gamma_b*betaBar + iU*gamma_a*alphaBar;

  // get nBar (eq A3, B&L 1986)
  double nBar = -F1*betaBar
      - (this->Am/this->Bm)*(betaBar*F3 - (F5+F6)*alphaBar)
        /(this->gas_constant_universal*this->T*(1.+F1));

  // get aBar
  double aBar = 0.;
  LOOP_k_N(this->n_species) {
    LOOP_l_N(this->n_species) {
      int apos = k*this->n_species + l;
      aBar += dN[k]*dN[l]*this->A_IJ[apos];
    }
  }
  aBar /= this->Am;

  // Output (eq 7, 19, B&L 1986)
  double D, C;
  D = D1*D4 - D2*D3;

  double sum = 0;
  LOOP_k_N(this->n_species) {
    sum += dN[k]*dN[k]*dN[k]/this->X[k]/this->X[k];
  }
  C = this->gas_constant_universal*this->T*(-sum + 3.*nBar*(betaBar*F1)*(betaBar*F1) + 2.*(betaBar*F1)*(betaBar*F1)*(betaBar*F1))
      + this->Am/this->Bm * (3.*betaBar*betaBar*(2*alphaBar-betaBar)*(F3+F6) - 2.*betaBar*betaBar*betaBar*F4 - 3*betaBar*aBar*F6);

  Dout = D;
  Cout = C;
}


//----------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& os, const PengRobinson& thermal) {
  os << "--------------------------------------------------------------------------------" << std::endl;
  os << "Number of species: " << thermal.n_species << " ";
  LOOP_k_N(thermal.n_species)
      os << thermal.species[k] << " ";
  os << std::endl;

  os  << "P = " << thermal.P << std::endl
      << "T = " << thermal.T << std::endl
      << "rho = " << thermal.rho << std::endl
      << "Z = " << thermal.Z << std::endl
      << "X = ";
  LOOP_k_N(thermal.n_species)
      os << thermal.X[k] << " ";
  os << std::endl;
  os << "Y = ";
  LOOP_k_N(thermal.n_species)
      os << thermal.Y[k] << " ";
  os << std::endl;
  os  << "MWm = " << thermal.MW_M << std::endl
      << "E = " << thermal.E << std::endl
      << "Partial H = ";
  LOOP_k_N(thermal.n_species)
      os << thermal.hspecies[k] << " ";
  os << std::endl;
  os << "Partial mu = ";
  LOOP_k_N(thermal.n_species)
      os << thermal.muspecies[k] << " ";
  os << std::endl;
  os << "Partial s = ";
  LOOP_k_N(thermal.n_species)
      os << thermal.sspecies[k] << " ";
  os << std::endl;
  os << "sos = " << thermal.sos << std::endl;
  os << "cp = " << thermal.cp << std::endl;
  os << "--------------------------------------------------------------------------------" << std::endl;
  return os;
}

