import 'dart:math';
import 'package:askreatif_app/group_params.dart';
import 'package:flutter/material.dart';
import 'package:fl_chart/fl_chart.dart';

import 'amn_selected_groups.dart';

void main() => runApp(MyApp());

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      debugShowCheckedModeBanner: false,
      title: 'POS Perfume Calculator',
      theme: ThemeData(
        primarySwatch: Colors.blue,
        visualDensity: VisualDensity.adaptivePlatformDensity,
      ),
      home: SolverPage(),
    );
  }
}

class Compound {
  final String name;
  final double MW;
  final double Psat;
  final double Thr;
  final Map<String, int> groups;

  Compound({
    required this.name,
    required this.MW,
    required this.Psat,
    required this.Thr,
    required this.groups,
  });
}

final List<Compound> compounds = [
  Compound(
    name: "Geraniol",
    MW: 154.25,
    Psat: 2.67, // Using values closer to MATLAB
    Thr: 0.0000248,
    groups: {"CH3": 3, "CH2": 3, "CH=C": 2, "OH": 1},
  ),
  Compound(
    name: "Romandolide",
    MW: 270.36,
    Psat: 0.1,
    Thr: 0.00004,
    groups: {
      "CH3": 3,
      "CH2": 4,
      "CH": 1,
      "C": 1,
      "CH3CO": 1,
      "CH2CO": 1,
      "CH3COO": 1,
    },
  ),
  Compound(
    name: "Cashmeran",
    MW: 206.32,
    Psat: 1.24,
    Thr: 0.00001,
    groups: {"CH3": 5, "CH2": 2, "CH": 1, "C": 2, "C=C": 1, "CH2CO": 1},
  ),
  Compound(
    name: "Javanol",
    MW: 222.37,
    Psat: 0.03,
    Thr: 0.00000002,
    groups: {"CH3": 4, "CH2": 5, "CH": 3, "C": 3, "OH": 1},
  ),
  Compound(
    name: "Benzyl Acetate",
    MW: 150.18,
    Psat: 21.86,
    Thr: 0.000332,
    groups: {"ACH": 5, "ACCH2": 1, "CH3COO": 1},
  ),
  Compound(
    name: "Vanillin",
    MW: 152.2,
    Psat: 0.016,
    Thr: 0.000000187,
    groups: {"ACH": 3, "AC": 2, "ACOH": 1, "CHO": 1, "CH3O": 1},
  ),
  Compound(
    name: "Iso E Super",
    MW: 234.38,
    Psat: 0.231,
    Thr: 0.00000093,
    groups: {"CH3": 5, "CH2": 2, "CH": 1, "C": 2, "C=C": 1, "CH3CO": 1},
  ),
  Compound(
    name: "Isoamyl Acetate",
    MW: 130.19,
    Psat: 5.33,
    Thr: 0.002,
    groups: {"CH3": 2, "CH2": 2, "CH": 1, "CH3COO": 1},
  ),
  Compound(
    name: "Macrolide",
    MW: 240.39,
    Psat: 0.08,
    Thr: 0.001,
    groups: {"CH3": 3, "CH": 5, "CH=C": 2, "OH": 2, "CH2COO": 1, "COO": 2},
  ),
  Compound(
    name: "Benzaldehyde",
    MW: 164.2,
    Psat: 0.8,
    Thr: 0.042,
    groups: {"ACH": 5, "AC": 2, "CHO": 1},
  ),
  Compound(
    name: "Eugenol",
    MW: 158.24,
    Psat: 0.00985,
    Thr: 0.011,
    groups: {"ACH": 3, "AC": 2, "ACCH2": 1, "ACOH": 1, "OCH3": 1, "CH=CH2": 1},
  ),
  Compound(
    name: "Oenanthic Ether",
    MW: 142.2,
    Psat: 0.0427,
    Thr: 0.002,
    groups: {"CH3": 2, "CH2": 5, "CH2COO": 1},
  ),
  Compound(
    name: "Cis-3-Hexenyl Acetate",
    MW: 184.32,
    Psat: 2.14,
    Thr: 0.000015,
    groups: {"CH3": 1, "CH2": 3, "CH=CH": 1, "CH3COO": 1},
  ),
  Compound(
    name: "Dodecanal",
    MW: 194,
    Psat: 1.99,
    Thr: 0.0005,
    groups: {"CH3": 1, "CH2": 3, "CH": 1},
  ),
  Compound(
    name: "Sandalmysore Core",
    MW: 192.3,
    Psat: 0.128,
    Thr: 0.051,
    groups: {"CH3": 4, "CH2": 3, "CH": 1, "C": 1, "CH=C": 2, "OH": 1},
  ),
  Compound(
    name: "a- ionone",
    MW: 136.2,
    Psat: 1.8,
    Thr: 0.00003,
    groups: {
      "CH3": 3,
      "CH2": 2,
      "CH": 1,
      "C": 1,
      "CH=C": 2,
      "CH=CH": 1,
      "CH3CO": 1,
    },
  ),
  Compound(
    name: "Limonene",
    MW: 174.2,
    Psat: 2.05,
    Thr: 0.00245,
    groups: {"CH3": 2, "CH2": 3, "CH": 1, "C": 1, "CH2=C": 1, "CH=C": 1},
  ),
  Compound(
    name: "Frukton",
    MW: 236.40,
    Psat: 0.181,
    Thr: 0.00000406,
    groups: {"CH3": 2, "C": 1, "CH2CO": 1, "CH2O": 1},
  ),
  Compound(
    name: "Ambroxan",
    MW: 284.44,
    Psat: 1.25,
    Thr: 0.0000029,
    groups: {"CH3": 4, "CH2": 6, "CH": 2, "C": 3, "CH2O": 1},
  ),
  Compound(
    name: "Helvetolida",
    MW: 154.25,
    Psat: 1.33,
    Thr: 0.0000017,
    groups: {"CH3": 6, "CH2": 5, "CH": 1, "C": 2, "CHO": 1, "CH2COO": 1},
  ),
  Compound(
    name: "Heksil Acetate",
    MW: 226.31,
    Psat: 0.0185,
    Thr: 0.000115,
    groups: {"CH3": 1, "CH2": 5, "CH3COO": 1},
  ),
  Compound(
    name: "Hedione",
    MW: 154.25,
    Psat: 0.00071,
    Thr: 0.015,
    groups: {"CH3": 2, "CH2": 5, "CH": 2, "CH2CO": 1, "CH2COO": 1},
  ),
  Compound(
    name: "Nerol",
    MW: 156.26,
    Psat: 2.39,
    Thr: 0.0000006,
    groups: {"CH3": 3, "CH2": 3, "CH=C": 2, "OH": 1},
  ),
  Compound(
    name: "Citronellol",
    MW: 172.27,
    Psat: 8.6,
    Thr: 0.018,
    groups: {"CH3": 3, "CH2": 4, "CH=C": 1, "CH": 1, "OH": 1},
  ),
  Compound(
    name: "Octyl Acetate",
    MW: 172,
    Psat: 25.86,
    Thr: 0.00019,
    groups: {"CH3": 1, "CH2": 7, "CH3COO": 1},
  ),
  Compound(
    name: "Cinnamaldehyde",
    MW: 196.20,
    Psat: 3.85,
    Thr: 0.05,
    groups: {"ACH": 5, "AC": 1, "CH=CH": 1, "CHO": 1},
  ),
  Compound(
    name: "b-caryophyllenol",
    MW: 222.37,
    Psat: 0.2173,
    Thr: 0.000845,
    groups: {"CH3": 3, "CH2": 7, "CH": 2, "C": 3, "OH": 1},
  ),
  Compound(
    name: "Coumarin",
    MW: 132.16,
    Psat: 1.33,
    Thr: 0.034,
    groups: {"ACH": 4, "AC": 2, "CH=CH": 1, "COO": 1},
  ),
  Compound(
    name: "Evernyl",
    MW: 208.34,
    Psat: 0.002399,
    Thr: 0.000000247,
    groups: {"CH3": 3, "ACH": 1, "AC": 5, "COO": 1, "OH": 2},
  ),
  Compound(
    name: "Bacdanol",
    MW: 156.26,
    Psat: 0.142,
    Thr: 0.00085,
    groups: {"CH3": 4, "CH2": 4, "CH": 1, "C": 1, "CH=C": 2, "OH": 1},
  ),
  Compound(
    name: "Menthol",
    MW: 156.27,
    Psat: 7.99,
    Thr: 0.0021,
    groups: {"CH3": 3, "CH2": 3, "CH": 4, "OH": 1},
  ),
  Compound(
    name: "Heptyl Acetate",
    MW: 158.24,
    Psat: 0.1599,
    Thr: 0.00042,
    groups: {"CH3": 1, "CH2": 6, "CHCOO": 1},
  ),
  Compound(
    name: "Caryophyllene Oxide",
    MW: 220.35,
    Psat: 1,
    Thr: 0.41,
    groups: {"CH3": 3, "CH2": 5, "C": 3, "CH": 2, "CHO": 1, "CH2=C": 1},
  ),
  Compound(name: "Water", MW: 18.15, Psat: 0.317, Thr: 10, groups: {"H2O": 1}),
  Compound(
    name: "Ethanol",
    MW: 46.07,
    Psat: 0.727,
    Thr: 0.0533,
    groups: {"CH3": 1, "CH2": 1, "OH": 1},
  ),
  Compound(
    name: "DPG",
    MW: 76.09,
    Psat: 0.013,
    Thr: 0.00143,
    groups: {"CH3": 2, "CH2": 1, "CH": 2, "OH": 2, "CH2O": 1},
  ),
];

enum SolventOption { waterEthanol, water, ethanol, dpg }

Map<SolventOption, String> solventNames = {
  SolventOption.waterEthanol: "Water + Ethanol",
  SolventOption.water: "Water",
  SolventOption.ethanol: "Ethanol",
  SolventOption.dpg: "DPG",
};

class SolverPage extends StatefulWidget {
  @override
  _SolverPageState createState() => _SolverPageState();
}

class _SolverPageState extends State<SolverPage> {
  bool isLoading = false;
  double progressValue = 0.0;

  List<String> selectedCompounds = [
    "Geraniol",
    "Romandolide",
    "Cashmeran",
    "Javanol",
    "Benzyl Acetate",
    "Vanillin",
  ];
  SolventOption selectedSolvent = SolventOption.waterEthanol;

  double totalSolventFraction = 0.3;
  double ethanolRatio = 0.4;
  double waterRatio = 0.6;

  List<double> fractions = [];
  String result = "";
  List<FlSpot> ternaryPoints = [];

  @override
  void initState() {
    super.initState();
    generateRandomFractions();
  }

  void generateRandomFractions() {
    Random random = Random();
    List<double> tempFractions = List.generate(
      6,
      (_) => random.nextDouble(),
    ); // 6 untuk fragrance compounds
    double tempSum = tempFractions.reduce((a, b) => a + b);
    tempFractions =
        tempFractions
            .map((f) => f / tempSum * (1 - totalSolventFraction))
            .toList();
    fractions = tempFractions;
  }

  // Calculate UNIFAC activity coefficients
  List<double> calculateImprovedUnifacCoefficients(
    List<Compound> comps,
    List<double> moleFractions,
  ) {
    const double T = 298.15; // Temperature in K
    const double R = 8.314; // Gas constant
    final int nComps = comps.length;

    // Initialize activity coefficients
    List<double> gammas = List.filled(nComps, 1.0);

    // Step 1: Calculate molecular parameters (r and q)
    List<double> r = List.filled(nComps, 0.0);
    List<double> q = List.filled(nComps, 0.0);

    for (int i = 0; i < nComps; i++) {
      Compound comp = comps[i];
      comp.groups.forEach((group, count) {
        if (groupParams.containsKey(group)) {
          r[i] += count * groupParams[group]![0]; // Volume parameter
          q[i] += count * groupParams[group]![1]; // Surface area parameter
        }
      });
    }

    // Step 2: Calculate mixture properties
    double sumRX = 0.0;
    double sumQX = 0.0;

    for (int i = 0; i < nComps; i++) {
      sumRX += r[i] * moleFractions[i];
      sumQX += q[i] * moleFractions[i];
    }

    // Prevent division by zero
    if (sumRX <= 0) sumRX = 1e-10;
    if (sumQX <= 0) sumQX = 1e-10;

    // Calculate volume and surface area fractions
    List<double> phi = List.filled(nComps, 0.0); // Volume fractions
    List<double> theta = List.filled(nComps, 0.0); // Surface area fractions

    for (int i = 0; i < nComps; i++) {
      phi[i] = (r[i] * moleFractions[i]) / sumRX;
      theta[i] = (q[i] * moleFractions[i]) / sumQX;
    }

    // Step 3: Calculate combinatorial part (sesuai Persamaan 2.12)
    List<double> lnGammaC = List.filled(nComps, 0.0);

    // Calculate l_i parameter for each component
    List<double> l = List.filled(nComps, 0.0);
    for (int i = 0; i < nComps; i++) {
      l[i] = 5.0 * (r[i] - q[i]) - (r[i] - 1.0);
    }

    // Calculate sum of x_j * l_j
    double sumXL = 0.0;
    for (int j = 0; j < nComps; j++) {
      sumXL += moleFractions[j] * l[j];
    }

    for (int i = 0; i < nComps; i++) {
      if (phi[i] > 0 && theta[i] > 0 && moleFractions[i] > 0) {
        double term1 = log(phi[i] / moleFractions[i]);
        double term2 = 5.0 * q[i] * log(theta[i] / phi[i]);
        double term3 = l[i];
        double term4 = (phi[i] / moleFractions[i]) * sumXL;

        lnGammaC[i] = term1 + term2 + term3 - term4;
      } else {
        lnGammaC[i] = 0.0;
      }
    }

    // Step 4: Calculate residual part - Enhanced for PQ2D
    List<double> lnGammaR = List.filled(nComps, 0.0);

    // Get all unique groups across all components
    Set<String> allGroups = {};
    for (var comp in comps) {
      allGroups.addAll(comp.groups.keys);
    }

    // Calculate group mole fractions in the mixture
    Map<String, double> X_m = {}; // Group mole fractions in mixture
    double totalGroupMoles = 0.0;

    for (String group in allGroups) {
      double groupMoles = 0.0;
      for (int i = 0; i < nComps; i++) {
        if (comps[i].groups.containsKey(group)) {
          groupMoles += comps[i].groups[group]! * moleFractions[i];
        }
      }
      X_m[group] = groupMoles;
      totalGroupMoles += groupMoles;
    }

    // Normalize group mole fractions
    if (totalGroupMoles > 0) {
      X_m.forEach((group, moles) {
        X_m[group] = moles / totalGroupMoles;
      });
    }

    // Calculate group activity coefficients in mixture (Persamaan 2.14)
    Map<String, double> lnGamma_k = {};

    for (String k in allGroups) {
      if (!groupParams.containsKey(k)) continue;

      double Q_k = groupParams[k]![1];

      // Calculate theta_m for all groups in mixture
      Map<String, double> theta_m = {};
      double sumQX = 0.0;

      for (String m in allGroups) {
        if (groupParams.containsKey(m)) {
          double Q_m = groupParams[m]![1];
          double X_m_val = X_m[m] ?? 0.0;
          sumQX += Q_m * X_m_val;
        }
      }

      // Calculate theta_m fractions
      for (String m in allGroups) {
        if (groupParams.containsKey(m)) {
          double Q_m = groupParams[m]![1];
          double X_m_val = X_m[m] ?? 0.0;
          if (sumQX > 0) {
            theta_m[m] = (Q_m * X_m_val) / sumQX;
          } else {
            theta_m[m] = 0.0;
          }
        }
      }

      // Calculate sum in ln term
      double sumThetaPsi = 0.0;
      for (String m in allGroups) {
        double theta_m_val = theta_m[m] ?? 0.0;
        double a_mk = 0.0;
        if (amn.containsKey(m) && amn[m]!.containsKey(k)) {
          a_mk = amn[m]![k]!;
        }
        double psi_mk = exp(-a_mk / (R * T));
        sumThetaPsi += theta_m_val * psi_mk;
      }

      // Calculate second sum term
      double secondSum = 0.0;
      for (String m in allGroups) {
        double theta_m_val = theta_m[m] ?? 0.0;
        double a_km = 0.0;
        if (amn.containsKey(k) && amn[k]!.containsKey(m)) {
          a_km = amn[k]![m]!;
        }
        double psi_km = exp(-a_km / (R * T));

        // Calculate denominator for this term
        double denominator = 0.0;
        for (String n in allGroups) {
          double theta_n_val = theta_m[n] ?? 0.0;
          double a_nm = 0.0;
          if (amn.containsKey(n) && amn[n]!.containsKey(m)) {
            a_nm = amn[n]![m]!;
          }
          double psi_nm = exp(-a_nm / (R * T));
          denominator += theta_n_val * psi_nm;
        }

        if (denominator > 0) {
          secondSum += (theta_m_val * psi_km) / denominator;
        }
      }

      // Apply Persamaan 2.14
      double lnGamma_k_val = 0.0;
      if (sumThetaPsi > 0) {
        lnGamma_k_val = Q_k * (1.0 - log(sumThetaPsi) - secondSum);
      }

      lnGamma_k[k] = lnGamma_k_val;
    }

    // Calculate group activity coefficients in pure components (Persamaan 2.14)
    Map<String, Map<String, double>> lnGamma_k_pure = {};

    for (int i = 0; i < nComps; i++) {
      lnGamma_k_pure[i.toString()] = {};

      // Get groups for this component
      Set<String> compGroups = comps[i].groups.keys.toSet();

      // Calculate group mole fractions in pure component
      Map<String, double> X_pure = {};
      double totalPureGroups = 0.0;

      for (String group in compGroups) {
        double count = comps[i].groups[group]?.toDouble() ?? 0.0;
        X_pure[group] = count;
        totalPureGroups += count;
      }

      // Normalize
      if (totalPureGroups > 0) {
        X_pure.forEach((group, count) {
          X_pure[group] = count / totalPureGroups;
        });
      }

      // Calculate group activity coefficients in pure component using Persamaan 2.14
      for (String k in compGroups) {
        if (!groupParams.containsKey(k)) continue;

        double Q_k = groupParams[k]![1];

        // Calculate theta_m for pure component
        Map<String, double> theta_m_pure = {};
        double sumQX_pure = 0.0;

        for (String m in compGroups) {
          if (groupParams.containsKey(m)) {
            double Q_m = groupParams[m]![1];
            double X_m_val = X_pure[m] ?? 0.0;
            sumQX_pure += Q_m * X_m_val;
          }
        }

        // Calculate theta_m fractions for pure component
        for (String m in compGroups) {
          if (groupParams.containsKey(m)) {
            double Q_m = groupParams[m]![1];
            double X_m_val = X_pure[m] ?? 0.0;
            if (sumQX_pure > 0) {
              theta_m_pure[m] = (Q_m * X_m_val) / sumQX_pure;
            } else {
              theta_m_pure[m] = 0.0;
            }
          }
        }

        // Calculate sum in ln term for pure component
        double sumThetaPsi_pure = 0.0;
        for (String m in compGroups) {
          double theta_m_val = theta_m_pure[m] ?? 0.0;
          double a_mk = 0.0;
          if (amn.containsKey(m) && amn[m]!.containsKey(k)) {
            a_mk = amn[m]![k]!;
          }
          double psi_mk = exp(-a_mk / (R * T));
          sumThetaPsi_pure += theta_m_val * psi_mk;
        }

        // Calculate second sum term for pure component
        double secondSum_pure = 0.0;
        for (String m in compGroups) {
          double theta_m_val = theta_m_pure[m] ?? 0.0;
          double a_km = 0.0;
          if (amn.containsKey(k) && amn[k]!.containsKey(m)) {
            a_km = amn[k]![m]!;
          }
          double psi_km = exp(-a_km / (R * T));

          // Calculate denominator for this term
          double denominator_pure = 0.0;
          for (String n in compGroups) {
            double theta_n_val = theta_m_pure[n] ?? 0.0;
            double a_nm = 0.0;
            if (amn.containsKey(n) && amn[n]!.containsKey(m)) {
              a_nm = amn[n]![m]!;
            }
            double psi_nm = exp(-a_nm / (R * T));
            denominator_pure += theta_n_val * psi_nm;
          }

          if (denominator_pure > 0) {
            secondSum_pure += (theta_m_val * psi_km) / denominator_pure;
          }
        }

        // Apply Persamaan 2.14 for pure component
        double lnGamma_k_pure_val = 0.0;
        if (sumThetaPsi_pure > 0) {
          lnGamma_k_pure_val =
              Q_k * (1.0 - log(sumThetaPsi_pure) - secondSum_pure);
        }

        lnGamma_k_pure[i.toString()]![k] = lnGamma_k_pure_val;
      }
    }

    // Calculate residual part for each component (Persamaan 2.13)
    for (int i = 0; i < nComps; i++) {
      double residualSum = 0.0;

      // Apply Persamaan 2.13: ln γᵢᴿ = Σₖ νₖ⁽ⁱ⁾[ln Γₖ - ln Γₖ⁽ⁱ⁾]
      comps[i].groups.forEach((group, count) {
        double lnGamma_mix = lnGamma_k[group] ?? 0.0;
        double lnGamma_pure = lnGamma_k_pure[i.toString()]?[group] ?? 0.0;

        // νₖ⁽ⁱ⁾ is the count of group k in component i
        residualSum += count * (lnGamma_mix - lnGamma_pure);
      });

      lnGammaR[i] = residualSum;
    }

    // Step 5: Combine parts and apply PQ2D specific adjustments
    for (int i = 0; i < nComps; i++) {
      // Check for numerical issues
      if (!lnGammaC[i].isFinite) lnGammaC[i] = 0.0;
      if (!lnGammaR[i].isFinite) lnGammaR[i] = 0.0;

      // Apply PQ2D specific correction factor
      double pq2dCorrection = calculatePQ2DCorrection(
        comps[i],
        moleFractions[i],
      );

      // Total activity coefficient
      double lnGammaTotal = lnGammaC[i] + lnGammaR[i] + pq2dCorrection;

      gammas[i] = exp(lnGammaTotal);

      // Safety bounds for numerical stability
      if (!gammas[i].isFinite || gammas[i] <= 0) {
        gammas[i] = 1.0;
      }

      // Apply reasonable bounds for perfumery applications
      if (gammas[i] > 100.0) gammas[i] = 100.0;
      if (gammas[i] < 0.01) gammas[i] = 0.01;
    }

    return gammas;
  }

  // PQ2D specific correction factor for perfumery compounds
  double calculatePQ2DCorrection(Compound comp, double moleFraction) {
    // For standard UNIFAC calculation without octant factors
    // This can be extended later if needed for specific PQ2D adjustments
    return 0.0;
  }

  // Helper function to check if a number is finite
  bool isFinite(double value) {
    return !value.isNaN && !value.isInfinite;
  }

  // Modified Newton-Raphson method that accepts an initial guess and target ratios
  List<double> solveForCustomOdorValuesWithInitialGuess(
    List<Compound> selectedComps,
    double solventFraction,
    List<double> initialGuess,
    List<double> targetRatios,
  ) {
    const int maxIterations = 100;
    const double tolerance = 1e-6;

    int numComps = selectedComps.length;
    List<double> x = List.from(initialGuess);

    try {
      for (int iter = 0; iter < maxIterations; iter++) {
        // Calculate current activity coefficients
        List<double> gamma = calculateImprovedUnifacCoefficients(
          selectedComps,
          x,
        );

        // Calculate current odor values
        List<double> odorValues = calculateOdorValues(selectedComps, x, gamma);

        // Check for NaN values
        if (odorValues.any((value) => value.isNaN || value.isInfinite)) {
          throw Exception("Non-finite odor values detected");
        }

        // Check if odor values match target ratios (within tolerance)
        bool ratiosMatch = true;
        for (int i = 1; i < odorValues.length; i++) {
          double targetRatio = targetRatios[i - 1];
          double actualRatio = odorValues[i] / odorValues[0];
          if ((actualRatio - targetRatio).abs() > tolerance) {
            ratiosMatch = false;
            break;
          }
        }

        if (ratiosMatch) {
          return x; // Solution found
        }

        // Set up Jacobian matrix untuk Newton-Raphson method
        List<List<double>> jacobian = List.generate(
          numComps - 1,
          (_) => List.filled(numComps - 1, 0.0),
        );

        List<double> f = List.filled(numComps - 1, 0.0);

        // Equations: OV_i = targetRatio_i * OV_1 for i = 2...n
        for (int i = 0; i < numComps - 1; i++) {
          f[i] = odorValues[i + 1] - targetRatios[i] * odorValues[0];

          // Numerical approximation of Jacobian
          for (int j = 0; j < numComps - 1; j++) {
            double h = max(1e-6, x[j + 1] * 1e-4);
            List<double> xPerturbed = List.from(x);
            xPerturbed[j + 1] += h;

            // Ensure sum remains 1 - solventFraction
            double sum = xPerturbed.sublist(1).reduce((a, b) => a + b);
            double excessOrDeficit = sum - (1.0 - solventFraction - x[0]);
            xPerturbed[0] = max(0.0, x[0] - excessOrDeficit);

            // Re-normalize if needed
            double newSum =
                xPerturbed[0] + xPerturbed.sublist(1).reduce((a, b) => a + b);
            if (newSum != 1.0 - solventFraction) {
              double scaleFactor = (1.0 - solventFraction) / newSum;
              for (int k = 0; k < xPerturbed.length; k++) {
                xPerturbed[k] *= scaleFactor;
              }
            }

            List<double> gammaPerturbed = calculateImprovedUnifacCoefficients(
              selectedComps,
              xPerturbed,
            );
            List<double> odorValuesPerturbed = calculateOdorValues(
              selectedComps,
              xPerturbed,
              gammaPerturbed,
            );

            // Check for non-finite values
            if (odorValuesPerturbed.any((val) => val.isNaN || val.isInfinite)) {
              jacobian[i][j] = 0.0;
            } else {
              double fPerturbed =
                  odorValuesPerturbed[i + 1] -
                  targetRatios[i] * odorValuesPerturbed[0];
              jacobian[i][j] = (fPerturbed - f[i]) / h;
            }
          }
        }

        // Check if Jacobian has non-finite values
        bool hasNonFiniteValues = false;
        for (var row in jacobian) {
          for (var elem in row) {
            if (elem.isNaN || elem.isInfinite) {
              hasNonFiniteValues = true;
              break;
            }
          }
          if (hasNonFiniteValues) break;
        }

        if (hasNonFiniteValues) {
          throw Exception("Non-finite values in Jacobian matrix");
        }

        // Solve J * dx = -f
        List<double> dx;
        try {
          dx = solveLinearSystem(jacobian, f.map((v) => -v).toList());
        } catch (e) {
          throw Exception("Failed to solve linear system: $e");
        }

        // Check for non-finite values in dx
        if (dx.any((value) => value.isNaN || value.isInfinite)) {
          throw Exception("Non-finite step values");
        }

        // Update x dengan damping factor
        double dampingFactor = 0.5;
        for (int i = 0; i < numComps - 1; i++) {
          x[i + 1] += dampingFactor * dx[i];
          x[i + 1] = max(0.0, min(1.0 - solventFraction, x[i + 1]));
        }

        // Adjust x[0] to maintain sum = 1 - solventFraction
        double sum = x.sublist(1).reduce((a, b) => a + b);
        x[0] = 1.0 - solventFraction - sum;

        // Ensure x[0] is also non-negative
        if (x[0] < 0) {
          double scale = (1.0 - solventFraction) / sum;
          for (int i = 1; i < numComps; i++) {
            x[i] *= scale;
          }
          x[0] = 0.0;
        }
      }

      return x;
    } catch (e) {
      print("Error in Newton-Raphson solver: $e");
      throw Exception("Newton-Raphson failed to converge: $e");
    }
  }

  // Enhanced solveLinearSystem to handle ill-conditioned matrices
  List<double> solveLinearSystem(List<List<double>> a, List<double> b) {
    int n = b.length;

    // Check for empty system
    if (n == 0) return [];

    // Create copies to avoid modifying the original arrays
    List<List<double>> matrix = List.generate(
      n,
      (i) => List<double>.from(a[i]),
    );
    List<double> rhs = List<double>.from(b);

    // Simple Gaussian elimination with partial pivoting
    for (int i = 0; i < n; i++) {
      // Find pivot
      int maxRow = i;
      double maxVal = matrix[i][i].abs();

      for (int k = i + 1; k < n; k++) {
        if (matrix[k][i].abs() > maxVal) {
          maxVal = matrix[k][i].abs();
          maxRow = k;
        }
      }

      // Check for singular matrix
      if (maxVal < 1e-10) {
        throw Exception("Matrix is nearly singular");
      }

      // Swap rows
      if (maxRow != i) {
        List<double> temp = matrix[i];
        matrix[i] = matrix[maxRow];
        matrix[maxRow] = temp;

        double t = rhs[i];
        rhs[i] = rhs[maxRow];
        rhs[maxRow] = t;
      }

      // Elimination with scaling
      for (int k = i + 1; k < n; k++) {
        // Prevent division by zero
        if (matrix[i][i].abs() < 1e-10) continue;

        double factor = matrix[k][i] / matrix[i][i];
        rhs[k] -= factor * rhs[i];

        for (int j = i; j < n; j++) {
          matrix[k][j] -= factor * matrix[i][j];
        }
      }
    }

    // Back substitution
    List<double> x = List.filled(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
      double sum = 0.0;
      for (int j = i + 1; j < n; j++) {
        sum += matrix[i][j] * x[j];
      }

      // Prevent division by zero
      if (matrix[i][i].abs() < 1e-10) {
        x[i] = 0.0; // Default value for singular solutions
      } else {
        x[i] = (rhs[i] - sum) / matrix[i][i];
      }

      // Check for non-finite value
      if (x[i].isNaN || x[i].isInfinite) {
        x[i] = 0.0;
      }
    }

    return x;
  }

  // Calculate odor values based on UNIFAC model
  List<double> calculateOdorValues(
    List<Compound> selectedComps,
    List<double> moleFractions,
    List<double> activityCoeffs,
  ) {
    const double R = 8.314462618; // Gas constant (J/mol·K)
    const double T = 298.15; // Temperature (K)

    // Step 1: Calculate partial pressure of each compound
    List<double> partialPressures = [];
    for (int i = 0; i < selectedComps.length; i++) {
      double Pi = activityCoeffs[i] * moleFractions[i] * selectedComps[i].Psat;

      // Make sure Pi is finite and positive
      Pi = (Pi.isNaN || Pi.isInfinite || Pi < 0) ? 1e-6 : Pi;
      partialPressures.add(Pi);
    }

    // Step 2: Calculate total pressure
    double Ptotal = partialPressures.reduce((a, b) => a + b);

    // Ensure Ptotal is not zero or too small
    if (Ptotal < 1e-10) Ptotal = 1e-10;

    // Step 3: Calculate vapor fraction and OV
    List<double> odorValues = [];
    for (int i = 0; i < selectedComps.length; i++) {
      double yi = partialPressures[i] / Ptotal;

      double threshold = selectedComps[i].Thr;

      // Hanya handle kasus threshold yang benar-benar tidak valid
      if (threshold <= 0 || threshold.isNaN || threshold.isInfinite) {
        threshold = 1e-6; // Gunakan nilai kecil untuk threshold tidak valid
      }
      // TIDAK mengubah threshold yang sudah valid

      double ov = (yi * selectedComps[i].MW * Ptotal) / (R * T * threshold);

      // Ensure OV is finite
      ov = (ov.isNaN || ov.isInfinite) ? 1e-6 : ov;
      odorValues.add(ov);
    }

    return odorValues;
  }

  void checkOdorValueRequirement(
    List<Compound> allCompounds,
    List<double> allOdorValues,
  ) {
    print("\n=== CHECKING OV REQUIREMENT ===");

    // ✅ PERBAIKAN: Identifikasi pelarut berdasarkan nama yang dikenal
    List<String> knownSolvents = ["Ethanol", "Water", "DPG"];
    List<int> solventIndices = [];
    List<int> fragranceIndices = [];

    for (int i = 0; i < allCompounds.length; i++) {
      if (knownSolvents.contains(allCompounds[i].name)) {
        solventIndices.add(i);
      } else {
        fragranceIndices.add(i);
      }
    }

    print("Solvent OVs:");
    for (int idx in solventIndices) {
      print(
        "  ${allCompounds[idx].name}: ${allOdorValues[idx].toStringAsFixed(6)}",
      );
    }

    print("Fragrance OVs:");
    for (int idx in fragranceIndices) {
      print(
        "  ${allCompounds[idx].name}: ${allOdorValues[idx].toStringAsFixed(6)}",
      );
    }

    // Check requirement: OV pelarut < OV senyawa fragrance
    bool requirementMet = true;
    for (int solventIdx in solventIndices) {
      for (int fragranceIdx in fragranceIndices) {
        if (allOdorValues[solventIdx] >= allOdorValues[fragranceIdx]) {
          print(
            "❌ REQUIREMENT VIOLATED: ${allCompounds[solventIdx].name} OV (${allOdorValues[solventIdx].toStringAsFixed(6)}) >= ${allCompounds[fragranceIdx].name} OV (${allOdorValues[fragranceIdx].toStringAsFixed(6)})",
          );
          requirementMet = false;
        }
      }
    }

    if (requirementMet) {
      print("✅ REQUIREMENT MET: All solvent OVs < fragrance OVs");
    }
  }

  void solveEquations() async {
    setState(() {
      isLoading = true;
      progressValue = 0.0;
    });

    // ✅ PERBAIKAN: Pastikan hanya mengambil senyawa fragrance, bukan DPG
    List<Compound> selectedList =
        selectedCompounds
            .map((name) => compounds.firstWhere((c) => c.name == name))
            .toList();

    // ✅ PERBAIKAN: DPG diperlakukan sebagai pelarut
    List<Compound> solvents = [];
    List<double> solventRatios = [];

    switch (selectedSolvent) {
      case SolventOption.waterEthanol:
        solvents = [
          compounds.firstWhere((c) => c.name == "Ethanol"),
          compounds.firstWhere((c) => c.name == "Water"),
        ];
        solventRatios = [ethanolRatio, waterRatio];
        break;
      case SolventOption.water:
        solvents = [compounds.firstWhere((c) => c.name == "Water")];
        solventRatios = [1.0];
        break;
      case SolventOption.ethanol:
        solvents = [compounds.firstWhere((c) => c.name == "Ethanol")];
        solventRatios = [1.0];
        break;
      case SolventOption.dpg:
        solvents = [compounds.firstWhere((c) => c.name == "DPG")];
        solventRatios = [1.0];
        break;
    }

    // Validasi bahwa DPG tidak ada dalam selectedList
    bool dpgInFragrance = selectedList.any((comp) => comp.name == "DPG");
    if (dpgInFragrance) {
      setState(() {
        result =
            "ERROR: DPG ditemukan dalam daftar senyawa fragrance. DPG harus hanya sebagai pelarut.";
        isLoading = false;
      });
      return;
    }

    // Generate random initial guess dengan time-based seed
    Random random = Random(DateTime.now().millisecondsSinceEpoch);
    List<double> initialGuess = List.generate(
      selectedList.length,
      (_) => random.nextDouble(),
    );
    double sum = initialGuess.reduce((a, b) => a + b);
    initialGuess =
        initialGuess
            .map((f) => f / sum * (1.0 - totalSolventFraction))
            .toList();

    // Generate random target ratios
    List<double> targetRatios = List.generate(
      selectedList.length - 1,
      (_) => 0.8 + 0.4 * random.nextDouble(),
    );

    // Try Newton-Raphson dengan multiple initial guesses
    List<double> optimizedFractions = [];
    bool success = false;
    String errorMessage = "";

    const int maxAttempts = 5;

    for (int attempt = 0; attempt < maxAttempts; attempt++) {
      try {
        if (attempt > 0) {
          initialGuess = List.generate(
            selectedList.length,
            (_) => random.nextDouble(),
          );
          sum = initialGuess.reduce((a, b) => a + b);
          initialGuess =
              initialGuess
                  .map((f) => f / sum * (1.0 - totalSolventFraction))
                  .toList();

          if (attempt >= 3) {
            targetRatios = List.generate(
              selectedList.length - 1,
              (_) => 0.7 + 0.6 * random.nextDouble(),
            );
          }

          setState(() {
            progressValue = attempt / maxAttempts;
          });
        }

        optimizedFractions = solveForCustomOdorValuesWithInitialGuess(
          selectedList,
          totalSolventFraction,
          initialGuess,
          targetRatios,
        );

        success = true;
        break;
      } catch (e) {
        errorMessage = e.toString();
        print("Attempt $attempt failed: $errorMessage");
      }
    }

    if (!success) {
      setState(() {
        result =
            "Failed to find solution after $maxAttempts attempts.\nLast error: $errorMessage";
        isLoading = false;
      });
      return;
    }

    // ✅ PERBAIKAN: Recalculate dengan pembagian yang jelas antara fragrance dan solvent
    List<double> totalFractions = [
      ...optimizedFractions, // Fragrance compounds
      ...solventRatios.map(
        (ratio) => ratio * totalSolventFraction,
      ), // Solvent compounds
    ];

    List<Compound> totalCompounds = [
      ...selectedList,
      ...solvents,
    ]; // Fragrance + Solvent

    List<double> gammas = calculateImprovedUnifacCoefficients(
      totalCompounds,
      totalFractions,
    );

    // ✅ PERBAIKAN: Hitung OV hanya untuk fragrance compounds
    List<double> ov = calculateOdorValues(
      selectedList, // Hanya fragrance compounds
      optimizedFractions, // Hanya fragrance fractions
      gammas.sublist(0, selectedList.length), // Hanya fragrance gammas
    );

    // ✅ PERBAIKAN: Update result dengan pembagian yang jelas
    String newResult = "PERFUME SOLUTION DETAILS:\n";
    newResult += "\n=== MOLE FRACTIONS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult +=
          "${selectedList[i].name}: ${optimizedFractions[i].toStringAsFixed(4)}\n";
    }

    newResult += "\n=== SOLVENT COMPOUNDS - MOLE FRACTIONS ===\n";
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${(solventRatios[i] * totalSolventFraction).toStringAsFixed(4)}\n";
    }

    newResult += "\n=== COMPOUNDS - ACTIVITY COEFFICIENTS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${gammas[i].toStringAsFixed(4)}\n";
    }

    newResult += "\n=== SOLVENT COMPOUNDS - ACTIVITY COEFFICIENTS ===\n";
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${gammas[selectedList.length + i].toStringAsFixed(4)}\n";
    }

    newResult += "\n=== COMPOUNDS - ODOR VALUES ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${ov[i].toStringAsFixed(4)}\n";
    }

    // ✅ PERBAIKAN: Hitung dan tampilkan OV untuk semua pelarut
    newResult += "\n=== SOLVENT COMPOUNDS - ODOR VALUES ===\n";
    for (int i = 0; i < solvents.length; i++) {
      double solventOV = calculateSolventOV(
        solvents[i],
        solventRatios[i] * totalSolventFraction,
        selectedList, // Fragrance compounds
        optimizedFractions, // Fragrance fractions
      );
      newResult += "${solvents[i].name}: ${solventOV.toStringAsFixed(6)}\n";
    }

    // ✅ PERBAIKAN: Validasi requirement dengan pembagian yang jelas
    checkOdorValueRequirement(totalCompounds, [
      ...ov,
      ...List.generate(
        solvents.length,
        (i) => calculateSolventOV(
          solvents[i],
          solventRatios[i] * totalSolventFraction,
          selectedList,
          optimizedFractions,
        ),
      ),
    ]);

    setState(() {
      result = newResult;
      isLoading = false;
    });
  }

  double calculateSolventOV(
    Compound solvent,
    double fraction,
    List<Compound> fragranceCompounds,
    List<double> fragranceFractions,
  ) {
    const double R = 8.314462618;
    const double T = 298.15;

    // Gabungkan semua komponen (fragrance + solvent)
    List<Compound> allCompounds = [...fragranceCompounds, solvent];
    List<double> allFractions = [...fragranceFractions, fraction];

    // Hitung koefisien aktivitas untuk SEMUA komponen sekaligus
    List<double> gammas = calculateImprovedUnifacCoefficients(
      allCompounds,
      allFractions,
    );

    // Ambil gamma untuk pelarut (index terakhir)
    double gammaSolvent = gammas.last;

    // Hitung tekanan parsial untuk SEMUA komponen
    List<double> partialPressures = [];
    for (int i = 0; i < allCompounds.length; i++) {
      double Pi = gammas[i] * allFractions[i] * allCompounds[i].Psat;
      Pi = (Pi.isNaN || Pi.isInfinite || Pi < 0) ? 1e-10 : Pi;
      partialPressures.add(Pi);
    }

    // Hitung tekanan total dari SEMUA komponen
    double Ptotal = partialPressures.reduce((a, b) => a + b);
    Ptotal = max(Ptotal, 1e-10);

    // Hitung fraksi uap pelarut
    double yi = partialPressures.last / Ptotal;

    // ✅ PERBAIKAN: Gunakan threshold asli tanpa modifikasi
    double threshold = solvent.Thr;
    if (threshold <= 0 || threshold.isNaN || threshold.isInfinite) {
      threshold = 1e-6; // Gunakan nilai kecil untuk threshold tidak valid
    }
    // TIDAK mengubah threshold yang sudah valid

    // Hitung OV tanpa pembatasan artificial
    double ov = (yi * solvent.MW * Ptotal) / (R * T * threshold);

    // Hanya pastikan nilai finite, TIDAK batasi dengan nilai maksimum
    if (!isFinite(ov) || ov < 0) {
      ov = 1e-6; // Nilai minimal untuk stabilitas numerik
    }

    return ov;
  }

  void debugOVCalculation(
    List<Compound> allCompounds,
    List<double> allFractions,
    List<double> odorValues,
  ) {
    print("\n=== DEBUG OV CALCULATION ===");
    print("Total compounds: ${allCompounds.length}");

    for (int i = 0; i < allCompounds.length; i++) {
      print("${allCompounds[i].name}:");
      print("  Fraction: ${allFractions[i].toStringAsFixed(6)}");
      print("  Psat: ${allCompounds[i].Psat}");
      print("  Threshold: ${allCompounds[i].Thr}");
      if (i < odorValues.length) {
        print("  OV: ${odorValues[i].toStringAsFixed(6)}");
      }
    }
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(
        title: Text(
          "Perfumery Octonary System",
          style: TextStyle(color: Colors.white),
        ),
        backgroundColor: Colors.purple[800],
      ),
      body: Padding(
        padding: const EdgeInsets.all(16.0),
        child: SingleChildScrollView(
          child: Column(
            crossAxisAlignment: CrossAxisAlignment.start,
            children: [
              Card(
                elevation: 4,
                margin: EdgeInsets.only(bottom: 16),
                child: Padding(
                  padding: const EdgeInsets.all(16.0),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        "Perfumery Octonary System (6 Senyawa + 2 Pelarut)",
                        style: TextStyle(
                          fontWeight: FontWeight.bold,
                          fontSize: 16,
                        ),
                      ),
                      Divider(),
                      Text(
                        "Perhitungan ini menggunakan metode POS (Perfumery Octonary System) dengan 6 senyawa pewangi dan rasio tetap 60:40 air:ethanol sebagai pelarut.",
                        style: TextStyle(fontSize: 14),
                      ),
                    ],
                  ),
                ),
              ),

              Text(
                "Pilih 6 Senyawa:",
                style: TextStyle(fontWeight: FontWeight.bold, fontSize: 16),
              ),
              SizedBox(height: 10),

              // Compound selection dropdowns
              for (int i = 0; i < 6; i++)
                Padding(
                  padding: const EdgeInsets.only(bottom: 8.0),
                  child: DropdownButton<String>(
                    isExpanded: true,
                    value: selectedCompounds[i],
                    onChanged: (newValue) {
                      setState(() {
                        selectedCompounds[i] = newValue!;
                        generateRandomFractions();
                      });
                    },
                    items:
                        compounds
                            .where(
                              (c) =>
                                  !["Water", "Ethanol", "DPG"].contains(c.name),
                            )
                            .map(
                              (c) => DropdownMenuItem(
                                value: c.name,
                                child: Text(c.name),
                              ),
                            )
                            .toList(),
                  ),
                ),

              SizedBox(height: 16),

              Text(
                "Pilih Pelarut:",
                style: TextStyle(fontWeight: FontWeight.bold, fontSize: 16),
              ),
              SizedBox(height: 10),
              DropdownButton<SolventOption>(
                isExpanded: true,
                value: selectedSolvent,
                onChanged: (newValue) {
                  setState(() {
                    selectedSolvent = newValue!;
                    // Reset water/ethanol ratio to default if switching back to dual solvent
                    if (newValue == SolventOption.waterEthanol) {
                      waterRatio = 0.6;
                      ethanolRatio = 0.4;
                    }
                  });
                },
                items:
                    SolventOption.values
                        .map(
                          (option) => DropdownMenuItem(
                            value: option,
                            child: Text(solventNames[option]!),
                          ),
                        )
                        .toList(),
              ),
              SizedBox(height: 16),

              // Solvent information
              // Replace the existing solvent information Card with this:
              Card(
                color: Colors.blue[50],
                child: Padding(
                  padding: const EdgeInsets.all(16.0),
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        "Solvent Information",
                        style: TextStyle(fontWeight: FontWeight.bold),
                      ),
                      SizedBox(height: 8),

                      // Show ratio sliders only for water+ethanol option
                      if (selectedSolvent == SolventOption.waterEthanol) ...[
                        Text("Water : Ethanol Ratio"),
                        Row(
                          children: [
                            Expanded(
                              flex: (waterRatio * 100).round(),
                              child: Container(
                                height: 20,
                                color: Colors.blue[200],
                                child: Center(
                                  child: Text(
                                    "${(waterRatio * 100).round()}%",
                                    style: TextStyle(fontSize: 12),
                                  ),
                                ),
                              ),
                            ),
                            Expanded(
                              flex: (ethanolRatio * 100).round(),
                              child: Container(
                                height: 20,
                                color: Colors.amber[200],
                                child: Center(
                                  child: Text(
                                    "${(ethanolRatio * 100).round()}%",
                                    style: TextStyle(fontSize: 12),
                                  ),
                                ),
                              ),
                            ),
                          ],
                        ),
                        Slider(
                          value: waterRatio,
                          min: 0.1,
                          max: 0.9,
                          divisions: 8,
                          label: waterRatio.toStringAsFixed(1),
                          onChanged: (value) {
                            setState(() {
                              waterRatio = value;
                              ethanolRatio = 1.0 - value;
                            });
                          },
                        ),
                        Text(
                          "Water: ${(waterRatio * 100).round()}% | Ethanol: ${(ethanolRatio * 100).round()}%",
                          style: TextStyle(fontStyle: FontStyle.italic),
                        ),
                        SizedBox(height: 12),
                      ] else ...[
                        Text(
                          "Selected Solvent: ${solventNames[selectedSolvent]}",
                        ),
                        SizedBox(height: 12),
                      ],

                      Text("Total Solvent Fraction:"),
                      Slider(
                        value: totalSolventFraction,
                        min: 0.1,
                        max: 0.5,
                        divisions: 40,
                        label: totalSolventFraction.toStringAsFixed(2),
                        onChanged: (value) {
                          setState(() {
                            totalSolventFraction = value;
                            generateRandomFractions();
                          });
                        },
                      ),
                      if (selectedSolvent == SolventOption.waterEthanol)
                        Text(
                          "Ethanol: ${(ethanolRatio * totalSolventFraction).toStringAsFixed(3)} | Water: ${(waterRatio * totalSolventFraction).toStringAsFixed(3)}",
                          style: TextStyle(fontStyle: FontStyle.italic),
                        )
                      else
                        Text(
                          "${solventNames[selectedSolvent]}: ${totalSolventFraction.toStringAsFixed(3)}",
                          style: TextStyle(fontStyle: FontStyle.italic),
                        ),
                    ],
                  ),
                ),
              ),

              SizedBox(height: 24),

              // Calculate button
              Center(
                child:
                    isLoading
                        ? Column(
                          children: [
                            LinearProgressIndicator(value: progressValue),
                            SizedBox(height: 8),
                            Text(
                              "Calculating... ${(progressValue * 100).toStringAsFixed(1)}%",
                            ),
                          ],
                        )
                        : ElevatedButton(
                          onPressed: solveEquations,
                          style: ElevatedButton.styleFrom(
                            padding: EdgeInsets.symmetric(
                              horizontal: 32,
                              vertical: 12,
                            ),
                            backgroundColor: Colors.purple[700],
                            textStyle: TextStyle(fontSize: 16),
                          ),
                          child: Text(
                            "Calculate Optimal Odor Values",
                            style: TextStyle(color: Colors.white),
                          ),
                        ),
              ),

              SizedBox(height: 24),

              // Results display
              result.isNotEmpty
                  ? Container(
                    padding: EdgeInsets.all(16),
                    decoration: BoxDecoration(
                      color: Colors.grey.shade100,
                      borderRadius: BorderRadius.circular(8),
                      border: Border.all(color: Colors.grey.shade300),
                    ),
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        Text(
                          "POS Optimization Results",
                          style: TextStyle(
                            fontWeight: FontWeight.bold,
                            fontSize: 18,
                          ),
                        ),
                        Divider(),
                        Text(
                          result,
                          style: TextStyle(fontFamily: 'Courier', fontSize: 14),
                        ),
                      ],
                    ),
                  )
                  : Container(),

              SizedBox(height: 24),
            ],
          ),
        ),
      ),
    );
  }
}
