#include "Common/commonMacros.h"

#include "Physics/pengrobinson.h"

int main(int argc, char* argv[]) {

  {
  PengRobinson thermo;
  std::vector<double> binCoeff = {0, 0.1561, 0.1561, 0};
  thermo.init({"NC12H26", "N2"}, binCoeff);

  double X[2] = {0.37, 0.63};
  thermo.getCP(X);
  std::cout << thermo;
  }

  { // Peng & Robinson 1977 Mixture 27
  PengRobinson thermo;
  std::string
      C2 = "C2H6",
      C3 = "C3H8",
      nC4 = "C4H10,n-butane",
      nC5 = "C5H12,n-pentane",
      nC6 = "C6H14,n-hexane"
      ;

  thermo.init({C2, C3, nC4, nC5, nC6}, std::vector<double> (5*5, 0.));
  double X[5] = {0.3977, 0.2926, 0.1997, 0.0713, 0.0387};
  thermo.getCP(X);
  std::cout << thermo;
  }


}
