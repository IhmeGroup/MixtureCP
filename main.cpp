#include "Common/commonMacros.h"

#include "Physics/pengrobinson.h"

int main(int argc, char* argv[]) {

  // nC12, N2 mixture (Cordova et al 2011)
  {
  PengRobinson thermo;
  std::vector<double> binCoeff = {0, 0.1561, 0.1561, 0};
  thermo.init({"NC12H26", "N2"}, binCoeff);

  Eigen::VectorXd XC12Vec = Eigen::VectorXd::LinSpaced(500, 0.999, 0.3);
  Eigen::VectorXd PCVec = Eigen::VectorXd::Zero(500);
  Eigen::VectorXd TCVec = Eigen::VectorXd::Zero(500);

  double kout, lout, kGuess=3.5, lGuess=std::sqrt(700.);
  LOOP_k_N(500) {
    double X[2] = {XC12Vec[k], 1.-XC12Vec[k]};
    thermo.getCP(X, &kout, &lout, kGuess,lGuess);
    kGuess = kout;
    lGuess = lout;

    PCVec[k] = thermo.GetP();
    TCVec[k] = thermo.GetT();
  }

  std::cout << "nC12/N2 Validation:" << std::endl;

  std::cout << "XC12 = ";
  LOOP_k_N(500)
    std::cout << XC12Vec[k] << " ";
  std::cout << std::endl;

  std::cout << "PC = ";
  LOOP_k_N(500)
    std::cout << PCVec[k] << " ";
  std::cout << std::endl;

  std::cout << "TC = ";
  LOOP_k_N(500)
    std::cout << TCVec[k] << " ";
  std::cout << std::endl;
  std::cout << "==============================================================" << std::endl;
  }

  // Ethane - n-Heptane CP line
  {
  PengRobinson thermo;
  std::string
      C2 = "C2H6",
      nC7 = "C7H16,n-heptane"
      ;

  thermo.init({C2, nC7}, std::vector<double> (2*2, 0.));

  Eigen::VectorXd XC2Vec = Eigen::VectorXd::LinSpaced(500, 0.001, 0.999);
  Eigen::VectorXd PCVec = Eigen::VectorXd::Zero(500);
  Eigen::VectorXd TCVec = Eigen::VectorXd::Zero(500);

  double kout, lout, kGuess=3.5, lGuess=std::sqrt(700.);
  LOOP_k_N(500) {
    double X[2] = {XC2Vec[k], 1.-XC2Vec[k]};
    thermo.getCP(X, &kout, &lout, kGuess,lGuess);
    kGuess = kout;
    lGuess = lout;

    PCVec[k] = thermo.GetP();
    TCVec[k] = thermo.GetT();
  }

  std::cout << "C2/C7 Validation:" << std::endl;

  std::cout << "XC2 = ";
  LOOP_k_N(500)
    std::cout << XC2Vec[k] << " ";
  std::cout << std::endl;

  std::cout << "PC = ";
  LOOP_k_N(500)
    std::cout << PCVec[k] << " ";
  std::cout << std::endl;

  std::cout << "TC = ";
  LOOP_k_N(500)
    std::cout << TCVec[k] << " ";
  std::cout << std::endl;
  std::cout << "==============================================================" << std::endl;
  }

  // Peng & Robinson 1977 Mixture 27
  {
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

  std::cout << "Peng & Robinson 1977, Mixture 27" << std::endl;
  std::cout << thermo;
  }
}
