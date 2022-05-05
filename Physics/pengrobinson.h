#ifndef PENGROBINSON_H
#define PENGROBINSON_H

#include "Common/commonMacros.h"

class PengRobinson {
 public:
  PengRobinson();
  ~PengRobinson();

  void init(const std::vector<std::string>& speciesNames, const std::vector<double> binCoeff);

  // Friend insertion to ostream
  friend std::ostream& operator<<(std::ostream& os, const PengRobinson& thermal);

  // These functions will NOT result in a consistent class state for all member variables
  double GetRho_SetMixture_TPY(const double T, const double P, const double* Y);
  double GetRho_SetMixture_TPX(const double T, const double P, const double* X);

  void GetTdew_PXvapor(const double P, const double* Xvapor, const double TGuess, const double* XliqGuess, double& Tdew, double* Xliq);

  // These functions will result in a consistent class state for all member variables
  void SetMixture_TRY(const double T, const double R, const double* Y);
  void SetMixture_TRX(const double T, const double R, const double* X);

  // Get methods
  inline double GetUniversalGasCosntant() const { return this->gas_constant_universal; }

  inline int GetNSpecies() const { return this->n_species; }
  inline const std::vector<double>& GetMW() const { return this->MW; }

  inline double GetT() const { return this->T; }
  inline double GetP() const { return this->P; }
  inline double GetRho() const { return this->rho; }
  inline double GetZ() const { return this->Z; }
  inline double GetBv() const { return this->expansivity; }
  inline const std::vector<double> GetY() const { return this->Y; }
  inline const std::vector<double> GetX() const { return this->X; }

  inline double GetMWm() const { return this->MW_M; }
  inline double GetInternalEnergy() const { return this->E; }

  inline const std::vector<double> GetPartialEnthalpies() const { return this->hspecies; }
  inline const std::vector<double> GetPartialChemEnergies() const { return this->muspecies; }
  inline const std::vector<double> GetPartialEntropies() const { return this->sspecies; }

  inline double getSos() const { return this->sos; }


 private:
  void ReadNasaPolynomials();
  void ReadCriticalProperties();
  void SetRealFluidConstants();

  inline void SetMassFractionFromY(const double* Y) {
    this->Y[this->n_species - 1] = 1.;
    LOOP_l_N(this->n_species - 1) {
      this->Y[l] = Y[l];
      this->Y[this->n_species - 1] -= this->Y[l];
    }
  }

  inline void SetMolecularWeightMixtureFromY() {
    this->MW_M = 0.;
    LOOP_l_N(this->n_species)
        this->MW_M += this->Y[l] / MW[l];
    this->MW_M = 1. / this->MW_M;

    this->gas_constant = this->gas_constant_universal / this->MW_M;
  }

  inline void SetMolarFractionFromY() {
    LOOP_l_N(this->n_species)
        this->X[l] = this->Y[l] * this->MW_M / this->MW[l];
  }

  inline void SetMolarFractionFromX(const double* X) {
    this->X[this->n_species - 1] = 1.;
    LOOP_l_N(this->n_species - 1) {
      this->X[l] = X[l];
      this->X[this->n_species - 1] -= this->X[l];
    }
  }

  inline void SetMolecularWeightMixtureFromX() {
    this->MW_M = 0.;
    LOOP_l_N(this->n_species)
        this->MW_M += this->X[l] * this->MW[l];

    this->gas_constant = this->gas_constant_universal / this->MW_M;
  }

  inline void SetMassFractionFromX() {
    LOOP_l_N(this->n_species)
        this->Y[l] = this->X[l] * this->MW[l] / this->MW_M;
  }

  void SyncRealFluidThermodynamicsFromTemperatureDensity();
  void SyncIdealFluidThermodynamicsFromTemperature();

  void SyncRZFromPressureTemperature();
  void SyncPZFromTemperatureDensity();
  void SyncEFromTemperatureDensity();
  void SyncPartialPropertiesFromTemperatureDensity();

  // Utilities
  // Test if current T and P can give more than one solution
  bool CheckThreeRoots();
  void ComputeMixtureCrit(double& Tcmix, double& Pcmix, double& wcmix, double& Zcmix);

  // For dew point calulcation only
  void getPhifromTandPandX(const double T_in, const double P_in, const double* X_in, double* phiOut);

  // For CP
public:
  // This will create a consistent thermo state at CP
  void getCP(double* X, double* kout=nullptr, double* lout=nullptr, double kappaGuess=-1., double lambdaGuess=-1.);
private:
  double eps = 1e-10;
  void solveCP(double kappaGuess, double firstlGuess, double& kout, double& lout);
  double lambdaFromKappa(double kappa_in, double lambdaGuess);
  void objectiveFun(double lambda_in, double kappa_in, double& Dout, double& Cout);

 private:
  double gas_constant_universal;
  double N_Av;
  double Boltzman;

  int n_species;
  std::vector<std::string> species;

  double T;
  double P;
  double MW_M;
  double gas_constant;

  double rho;
  double Z;

  double E;

  double cp, cv, gamma;
  double sos;

  std::vector<double> A_I;
  std::vector<double> dA_IdT;
  std::vector<double> d2A_IdT2;
  std::vector<double> A_IJ;
  double Am;
  double Bm;
  double dAmdT;
  double d2AmdT2;
  std::vector<double> dAmdN;
  std::vector<double> d2AmdTdN;
  double K1;
  std::vector<double> dK1dN;

  double BT;
  double betaT;
  double expansivity;
  double dPdT;
  double dPdV;
  std::vector<double> dPdN;
  std::vector<double> dVdN;

  std::vector<double> hspecies;
  std::vector<double> sspecies;
  std::vector<double> muspecies;
  std::vector<double> edspecies;

  std::vector<double> X;
  std::vector<double> Y;
  std::vector<double> MW;

  std::vector<double> Tcrit;
  std::vector<double> Pcrit;
  std::vector<double> rhocrit;
  std::vector<double> Vcrit;
  std::vector<double> Zcrit;
  std::vector<double> omega;

  int NasaCoef;
  std::vector<double> nasa_poly_coeff;
  std::vector<std::pair<double, double>> nasa_poly_bounds;
  std::vector<double> cpspecies_ig;
  std::vector<double> hspecies_ig;
  std::vector<double> sspecies_ig;
  double h_ig;
  double cp_ig;
  double cv_ig;

  std::vector<double> cst_a;
  std::vector<double> cst_b;
  std::vector<double> cst_c; // k

  std::vector<double> binCoeff;
};



#endif // PENGROBINSON_H
