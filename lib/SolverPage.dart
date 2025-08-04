import 'dart:math';

import 'package:askreatif_app/auth_service.dart';
import 'package:fl_chart/fl_chart.dart';
import 'package:flutter/material.dart';
import 'package:flutter/services.dart';

import 'amn_selected_groups.dart';
import 'group_params.dart';

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

class CompositionResult {
  final List<double> moleFractions;
  final List<double> massCompositions;
  final List<double> grams;
  final double totalMass;
  final List<String> compoundNames;

  CompositionResult({
    required this.moleFractions,
    required this.massCompositions,
    required this.grams,
    required this.totalMass,
    required this.compoundNames,
  });
}

class SolverPage extends StatefulWidget {
  @override
  _SolverPageState createState() => _SolverPageState();
}

class _SolverPageState extends State<SolverPage> {
  bool isLoading = false;
  double progressValue = 0.0;
  TextEditingController volumeController = TextEditingController(text: '30');
  double selectedVolume = 30.0; // mL - akan diupdate dari text input
  double perfumeConcentration = 15.0; // jP dari MATLAB (15% default)

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

    // TAMBAHKAN LISTENER UNTUK TEXT CONTROLLER
    volumeController.addListener(() {
      double? newVolume = double.tryParse(volumeController.text);
      if (newVolume != null && newVolume > 0) {
        setState(() {
          selectedVolume = newVolume;
        });
      }
    });
  }

  @override
  void dispose() {
    volumeController.dispose();
    super.dispose();
  }

  void generateRandomFractions() {
    Random random = Random();
    List<double> tempFractions = List.generate(6, (_) => random.nextDouble());
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
    final int nComps = comps.length;
    List<double> gammas = List.filled(nComps, 1.0);

    // Step 1: Calculate combinatorial part
    List<double> r = List.filled(nComps, 0.0);
    List<double> q = List.filled(nComps, 0.0);

    // Calculate r and q for each compound from its groups
    for (int i = 0; i < nComps; i++) {
      Compound comp = comps[i];
      comp.groups.forEach((group, count) {
        if (groupParams.containsKey(group)) {
          r[i] += count * groupParams[group]![0];
          q[i] += count * groupParams[group]![1];
        }
      });
    }

    // Calculate volume and surface area fractions
    List<double> volumeFractions = List.filled(nComps, 0.0);
    List<double> surfaceAreaFractions = List.filled(nComps, 0.0);

    double sumRX = 0.0;
    double sumQX = 0.0;

    for (int i = 0; i < nComps; i++) {
      sumRX += r[i] * moleFractions[i];
      sumQX += q[i] * moleFractions[i];
    }

    for (int i = 0; i < nComps; i++) {
      volumeFractions[i] = r[i] * moleFractions[i] / sumRX;
      surfaceAreaFractions[i] = q[i] * moleFractions[i] / sumQX;
    }

    // Calculate combinatorial part of activity coefficient
    List<double> lnGammaC = List.filled(nComps, 0.0);
    for (int i = 0; i < nComps; i++) {
      lnGammaC[i] =
          1 -
          volumeFractions[i] +
          log(volumeFractions[i]) -
          5 *
              q[i] *
              (1 -
                  volumeFractions[i] / surfaceAreaFractions[i] +
                  log(volumeFractions[i] / surfaceAreaFractions[i]));
    }

    // Step 2: Calculate residual part (simplified for implementation ease)
    List<double> lnGammaR = List.filled(nComps, 0.0);

    // Calculate group fractions across the mixture
    Map<String, double> groupFractions = {};
    Map<String, double> totalGroups = {}; // Changed from int to double

    // First, count total groups in the mixture
    for (int i = 0; i < nComps; i++) {
      comps[i].groups.forEach((group, count) {
        // FIX: Use the mole fraction directly without converting to int
        totalGroups[group] =
            (totalGroups[group] ?? 0.0) + count * moleFractions[i];
      });
    }

    // Calculate group fractions
    double totalGroupCount = totalGroups.values.reduce((a, b) => a + b);

    // Check for zero division
    if (totalGroupCount > 0) {
      totalGroups.forEach((group, count) {
        groupFractions[group] = count / totalGroupCount;
      });
    } else {
      // Default equal distribution if we have zero total
      totalGroups.keys.forEach((group) {
        groupFractions[group] = 1.0 / totalGroups.length;
      });
    }

    // Simple residual calculation (this is highly simplified)
    for (int i = 0; i < nComps; i++) {
      double residualSum = 0.0;

      comps[i].groups.forEach((group1, count1) {
        groupFractions.forEach((group2, fraction) {
          if (amn.containsKey(group1) && amn[group1]!.containsKey(group2)) {
            double a = amn[group1]![group2]!;
            double tau = exp(-a / (8.314 * T));
            // Protect against extreme values or zero fractions
            if (fraction > 0 && isFinite(tau) && isFinite(fraction)) {
              residualSum +=
                  count1 *
                  fraction *
                  (log(tau) - log(1 + (tau - 1) * fraction));
            }
          }
        });
      });

      lnGammaR[i] = residualSum;
    }

    // Combine combinatorial and residual parts
    for (int i = 0; i < nComps; i++) {
      // Make sure we don't have any NaN or infinity values
      if (!isFinite(lnGammaC[i])) lnGammaC[i] = 0;
      if (!isFinite(lnGammaR[i])) lnGammaR[i] = 0;

      gammas[i] = exp(lnGammaC[i] + lnGammaR[i]);

      // Safety check for extreme values
      if (!isFinite(gammas[i])) gammas[i] = 1.0;
    }

    return gammas;
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
    List<double> targetRatios, // Add target odor value ratios
  ) {
    const int maxIterations = 100;
    const double tolerance = 1e-6;

    int numComps = selectedComps.length;
    List<double> x = List.from(initialGuess); // Use the provided initial guess

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

        // Set up Jacobian matrix for Newton-Raphson method
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
              jacobian[i][j] = 0.0; // Use a safe default
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

        // Update x with damping factor to improve stability
        double dampingFactor = 0.5;
        for (int i = 0; i < numComps - 1; i++) {
          x[i + 1] += dampingFactor * dx[i];
          // Ensure no negative values
          x[i + 1] = max(0.0, min(1.0 - solventFraction, x[i + 1]));
        }

        // Adjust x[0] to maintain sum = 1 - solventFraction
        double sum = x.sublist(1).reduce((a, b) => a + b);
        x[0] = 1.0 - solventFraction - sum;

        // Ensure x[0] is also non-negative
        if (x[0] < 0) {
          // If x[0] becomes negative, rescale all other values
          double scale = (1.0 - solventFraction) / sum;
          for (int i = 1; i < numComps; i++) {
            x[i] *= scale;
          }
          x[0] = 0.0;
        }
      }

      // If we get here, Newton-Raphson didn't converge
      // Return current best estimate
      return x;
    } catch (e) {
      print("Error in Newton-Raphson solver: $e");
      // Rethrow to allow for retry with different initial guess
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

      // ✅ PERBAIKAN: Jangan ubah threshold yang sudah valid!
      double threshold = selectedComps[i].Thr;

      // Hanya ubah jika threshold benar-benar tidak valid (0, NaN, atau infinity)
      if (threshold <= 0 || threshold.isNaN || threshold.isInfinite) {
        threshold = 1.0; // Set ke nilai tinggi untuk threshold tidak valid
      }
      // JANGAN ubah threshold yang sudah valid meski kecil!

      double ov = (yi * selectedComps[i].MW * Ptotal) / (R * T * threshold);

      // Ensure OV is finite
      ov = (ov.isNaN || ov.isInfinite) ? 1.0 : ov;
      odorValues.add(ov);
    }

    return odorValues;
  }

  void checkOdorValueRequirement(
    List<Compound> selectedComps,
    List<double> odorValues,
  ) {
    print("\n=== CHECKING OV REQUIREMENT ===");

    // Cari index pelarut (Ethanol dan Water)
    List<int> solventIndices = [];
    List<int> fragranceIndices = [];

    for (int i = 0; i < selectedComps.length; i++) {
      if (selectedComps[i].name == "Ethanol" ||
          selectedComps[i].name == "Water") {
        solventIndices.add(i);
      } else {
        fragranceIndices.add(i);
      }
    }

    print("Solvent OVs:");
    for (int idx in solventIndices) {
      print(
        "  ${selectedComps[idx].name}: ${odorValues[idx].toStringAsFixed(6)}",
      );
    }

    print("Fragrance OVs:");
    for (int idx in fragranceIndices) {
      print(
        "  ${selectedComps[idx].name}: ${odorValues[idx].toStringAsFixed(6)}",
      );
    }

    // Check requirement: OV pelarut < OV senyawa
    bool requirementMet = true;
    for (int solventIdx in solventIndices) {
      for (int fragranceIdx in fragranceIndices) {
        if (odorValues[solventIdx] >= odorValues[fragranceIdx]) {
          print(
            "❌ REQUIREMENT VIOLATED: ${selectedComps[solventIdx].name} OV (${odorValues[solventIdx].toStringAsFixed(6)}) >= ${selectedComps[fragranceIdx].name} OV (${odorValues[fragranceIdx].toStringAsFixed(6)})",
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

    List<Compound> selectedList =
        selectedCompounds
            .map((name) => compounds.firstWhere((c) => c.name == name))
            .toList();

    // Get solvent compound(s) based on selection
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

    // Generate random initial guess with time-based seed
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

    // Try Newton-Raphson with multiple initial guesses if needed
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

    // Update the fractions variable for composition calculation
    fractions = optimizedFractions;

    // Now recalculate all the values based on the optimized fractions
    List<double> totalFractions = [
      ...optimizedFractions,
      ...solventRatios.map((ratio) => ratio * totalSolventFraction),
    ];

    List<Compound> totalCompounds = [...selectedList, ...solvents];

    List<double> gammas = calculateImprovedUnifacCoefficients(
      totalCompounds,
      totalFractions,
    );
    List<double> ov = calculateOdorValues(
      selectedList,
      optimizedFractions,
      gammas.sublist(0, selectedList.length),
    );

    // ============ TAMBAHKAN PERHITUNGAN KOMPOSISI MASSA DI SINI ============

    // Calculate composition and mass using MATLAB formulas
    CompositionResult compositionResult = calculateMassComposition();

    // Calculate standard deviation to check how well we met our target ratios
    double firstOV = ov[0];
    List<double> actualRatios = ov.sublist(1).map((v) => v / firstOV).toList();
    List<double> ratioDiffs = List.generate(
      targetRatios.length,
      (i) => (actualRatios[i] - targetRatios[i]).abs(),
    );
    double avgRatioDiff =
        ratioDiffs.reduce((a, b) => a + b) / ratioDiffs.length;

    // ============ UPDATE RESULT STRING DENGAN KEDUA HASIL ============

    String newResult = "PERFUME CALCULATION RESULTS (Complete Analysis):\n";
    newResult += "Volume: ${selectedVolume.toStringAsFixed(0)} mL\n";
    newResult +=
        "Concentration (jP): ${perfumeConcentration.toStringAsFixed(1)}%\n\n";

    newResult += "=== OPTIMIZED MOLE FRACTIONS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult +=
          "${selectedList[i].name}: ${optimizedFractions[i].toStringAsFixed(4)}\n";
    }

    // Show solvent fractions
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${(solventRatios[i] * totalSolventFraction).toStringAsFixed(4)}\n";
    }

    newResult += "\n=== ACTIVITY COEFFICIENTS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${gammas[i].toStringAsFixed(4)}\n";
    }

    // Show solvent activity coefficients
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${gammas[selectedList.length + i].toStringAsFixed(4)}\n";
    }

    newResult += "\n=== ODOR VALUES (OV) ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${ov[i].toStringAsFixed(4)}\n";
    }

    // Calculate and show solvent odor values
    newResult += "\n=== SOLVENT ODOR VALUES ===\n";
    List<double> solventOVs = [];

    for (int i = 0; i < solvents.length; i++) {
      double solventOV = calculateSolventOV(
        solvents[i],
        solventRatios[i] * totalSolventFraction,
        selectedList,
        optimizedFractions,
      );
      solventOVs.add(solventOV);
      newResult += "${solvents[i].name}: ${solventOV.toStringAsFixed(6)}\n";
    }

    // ============ TAMBAHKAN HASIL KOMPOSISI MASSA ============

    newResult += "\n" + "=" * 50 + "\n";
    newResult += "MASS COMPOSITION CALCULATIONS\n";
    newResult += "=" * 50 + "\n";

    newResult += "\n=== MASS COMPOSITIONS (x * MW) ===\n";
    for (int i = 0; i < compositionResult.compoundNames.length; i++) {
      newResult +=
          "${compositionResult.compoundNames[i]}: ${compositionResult.massCompositions[i].toStringAsFixed(4)}\n";
    }

    double totalMassComposition = compositionResult.massCompositions.reduce(
      (a, b) => a + b,
    );
    double jP = perfumeConcentration / 100.0;
    double za = totalMassComposition / (jP * selectedVolume);

    newResult +=
        "\nTotal Mass Composition: ${totalMassComposition.toStringAsFixed(4)}\n";
    newResult +=
        "jP (Perfume Concentration): ${perfumeConcentration.toStringAsFixed(1)}%\n";
    newResult += "Volume: ${selectedVolume.toStringAsFixed(0)} mL\n";
    newResult += "za Factor: ${za.toStringAsFixed(4)}\n";

    newResult +=
        "\n=== GRAMS NEEDED FOR ${selectedVolume.toStringAsFixed(0)}mL PERFUME ===\n";
    for (int i = 0; i < compositionResult.compoundNames.length; i++) {
      newResult +=
          "${compositionResult.compoundNames[i]}: ${compositionResult.grams[i].toStringAsFixed(4)} g\n";
    }

    newResult +=
        "\nTotal Mass: ${compositionResult.totalMass.toStringAsFixed(4)} g\n";
    newResult +=
        "Absolute Mass (for ${perfumeConcentration.toStringAsFixed(1)}%): ${(compositionResult.totalMass * jP).toStringAsFixed(4)} g\n";

    setState(() {
      result = newResult;
      isLoading = false;
    });
    return;
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

    // Gunakan threshold yang benar
    double threshold = solvent.Thr;
    if (threshold <= 0 || threshold.isNaN || threshold.isInfinite) {
      threshold = 1.0;
    }

    // Hitung OV
    double ov = (yi * solvent.MW * Ptotal) / (R * T * threshold);
    return isFinite(ov) ? ov : 0.0;
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

  CompositionResult calculateMassComposition() {
    List<Compound> selectedList =
        selectedCompounds
            .map((name) => compounds.firstWhere((c) => c.name == name))
            .toList();

    // Get solvent compound(s) based on selection
    List<Compound> solvents = [];
    List<double> solventRatios = [];
    List<String> solventNames = [];

    switch (selectedSolvent) {
      case SolventOption.waterEthanol:
        solvents = [
          compounds.firstWhere((c) => c.name == "Ethanol"),
          compounds.firstWhere((c) => c.name == "Water"),
        ];
        solventRatios = [ethanolRatio, waterRatio];
        solventNames = ["Ethanol", "Water"];
        break;
      case SolventOption.water:
        solvents = [compounds.firstWhere((c) => c.name == "Water")];
        solventRatios = [1.0];
        solventNames = ["Water"];
        break;
      case SolventOption.ethanol:
        solvents = [compounds.firstWhere((c) => c.name == "Ethanol")];
        solventRatios = [1.0];
        solventNames = ["Ethanol"];
        break;
      case SolventOption.dpg:
        solvents = [compounds.firstWhere((c) => c.name == "DPG")];
        solventRatios = [1.0];
        solventNames = ["DPG"];
        break;
    }

    // Combine all compounds
    List<Compound> allCompounds = [...selectedList, ...solvents];

    // Combine fractions (fragrance compounds + solvent)
    List<double> allFractions = [
      ...fractions, // fractions yang sudah ada untuk fragrance compounds
      ...solventRatios.map((ratio) => ratio * totalSolventFraction),
    ];

    // Get MW from each compound (seperti variabel M di MATLAB)
    List<double> MW = allCompounds.map((compound) => compound.MW).toList();

    // Hitung komposisi massa: komposisi = x * xT * MW
    List<double> massCompositions = [];
    double xT = 1.0; // dari MATLAB

    for (int i = 0; i < allFractions.length; i++) {
      double composition = allFractions[i] * xT * MW[i];
      massCompositions.add(composition);
    }

    double totalMassComposition = massCompositions.reduce((a, b) => a + b);

    // Rumus dari MATLAB: za = totalKomposisi / (jP * mL)
    double jP = perfumeConcentration / 100.0; // Convert % to decimal
    double za = totalMassComposition / (jP * selectedVolume);

    // Hitung gram: gram = komposisi / za
    List<double> grams = massCompositions.map((comp) => comp / za).toList();

    // Combine compound names
    List<String> allCompoundNames = [...selectedCompounds, ...solventNames];

    return CompositionResult(
      moleFractions: allFractions,
      massCompositions: massCompositions,
      grams: grams,
      totalMass: grams.reduce((a, b) => a + b),
      compoundNames: allCompoundNames,
    );
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
        automaticallyImplyLeading: false,
        actions: [
          // User info
          Padding(
            padding: EdgeInsets.symmetric(horizontal: 8.0),
            child: Row(
              children: [
                Icon(Icons.account_circle, color: Colors.white),
                SizedBox(width: 4),
                Text(
                  'askreatifperfume',
                  style: TextStyle(color: Colors.white70, fontSize: 12),
                ),
              ],
            ),
          ),
          // Logout button
          IconButton(
            icon: Icon(Icons.logout, color: Colors.white),
            tooltip: 'Logout',
            onPressed: () {
              _showLogoutDialog(context);
            },
          ),
        ],
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

              buildVolumeAndConcentrationControls(),

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
                            "Calculate Complete Analysis (OV + Composition)",
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
                          "Complete Perfume Analysis Results",
                          style: TextStyle(
                            fontWeight: FontWeight.bold,
                            fontSize: 18,
                            color: Colors.purple[800],
                          ),
                        ),
                        Divider(),
                        Container(
                          width: MediaQuery.of(context).size.width,
                          child: Text(
                            result,
                            style: TextStyle(
                              fontFamily: 'Courier',
                              fontSize: 12,
                            ),
                          ),
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

  Widget buildVolumeAndConcentrationControls() {
    return Card(
      elevation: 4,
      margin: EdgeInsets.only(bottom: 16),
      child: Padding(
        padding: const EdgeInsets.all(16.0),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text(
              "Volume & Concentration Settings",
              style: TextStyle(fontWeight: FontWeight.bold, fontSize: 16),
            ),
            Divider(),

            // Volume Input (Text Field)
            Row(
              children: [
                Expanded(
                  flex: 2,
                  child: Column(
                    crossAxisAlignment: CrossAxisAlignment.start,
                    children: [
                      Text(
                        "Volume Parfum (mL):",
                        style: TextStyle(fontWeight: FontWeight.w500),
                      ),
                      SizedBox(height: 8),
                      TextField(
                        controller: volumeController,
                        keyboardType: TextInputType.numberWithOptions(
                          decimal: true,
                        ),
                        inputFormatters: [
                          FilteringTextInputFormatter.allow(
                            RegExp(r'^\d*\.?\d*'),
                          ),
                        ],
                        decoration: InputDecoration(
                          border: OutlineInputBorder(),
                          hintText: '30',
                          suffixText: 'mL',
                          contentPadding: EdgeInsets.symmetric(
                            horizontal: 12,
                            vertical: 8,
                          ),
                        ),
                        onChanged: (value) {
                          double? newVolume = double.tryParse(value);
                          if (newVolume != null && newVolume > 0) {
                            setState(() {
                              selectedVolume = newVolume;
                            });
                          }
                        },
                      ),
                    ],
                  ),
                ),
                SizedBox(width: 16),
                Expanded(
                  flex: 1,
                  child: Container(
                    padding: EdgeInsets.all(12),
                    decoration: BoxDecoration(
                      color: Colors.blue[50],
                      borderRadius: BorderRadius.circular(8),
                      border: Border.all(color: Colors.blue[200]!),
                    ),
                    child: Column(
                      children: [
                        Icon(Icons.local_drink, color: Colors.blue[600]),
                        SizedBox(height: 4),
                        Text(
                          "${selectedVolume.toStringAsFixed(0)} mL",
                          style: TextStyle(
                            fontWeight: FontWeight.bold,
                            color: Colors.blue[800],
                          ),
                        ),
                      ],
                    ),
                  ),
                ),
              ],
            ),

            SizedBox(height: 20),

            // Concentration Slider (jP from MATLAB)
            Text(
              "Konsentrasi Parfum (jP): ${perfumeConcentration.toStringAsFixed(1)}%",
              style: TextStyle(fontWeight: FontWeight.w500),
            ),
            SizedBox(height: 8),

            // Slider dengan indikator visual
            Column(
              children: [
                Slider(
                  value: perfumeConcentration,
                  min: 5.0,
                  max: 30.0,
                  divisions: 25,
                  label: "${perfumeConcentration.toStringAsFixed(1)}%",
                  activeColor: Colors.purple[600],
                  inactiveColor: Colors.purple[100],
                  onChanged: (value) {
                    setState(() {
                      perfumeConcentration = value;
                    });
                  },
                ),

                SizedBox(height: 8),

                // Kategori parfum berdasarkan konsentrasi
                Container(
                  padding: EdgeInsets.all(8),
                  decoration: BoxDecoration(
                    color: _getPerfumeTypeColor().withOpacity(0.1),
                    borderRadius: BorderRadius.circular(6),
                    border: Border.all(
                      color: _getPerfumeTypeColor().withOpacity(0.3),
                    ),
                  ),
                  child: Row(
                    mainAxisAlignment: MainAxisAlignment.center,
                    children: [
                      Icon(
                        _getPerfumeTypeIcon(),
                        size: 16,
                        color: _getPerfumeTypeColor(),
                      ),
                      SizedBox(width: 6),
                      Text(
                        _getPerfumeTypeName(),
                        style: TextStyle(
                          fontWeight: FontWeight.w500,
                          color: _getPerfumeTypeColor(),
                          fontSize: 13,
                        ),
                      ),
                    ],
                  ),
                ),
              ],
            ),

            SizedBox(height: 12),

            // Info tambahan
            Container(
              padding: EdgeInsets.all(12),
              decoration: BoxDecoration(
                color: Colors.grey[50],
                borderRadius: BorderRadius.circular(6),
                border: Border.all(color: Colors.grey[300]!),
              ),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    "Informasi:",
                    style: TextStyle(fontWeight: FontWeight.w500, fontSize: 12),
                  ),
                  SizedBox(height: 4),
                  Text(
                    "• jP menentukan konsentrasi parfum dalam produk akhir\n"
                    "• Volume mempengaruhi total massa yang dibutuhkan\n"
                    "• Rumus: za = totalKomposisi / (jP × mL)",
                    style: TextStyle(fontSize: 11, color: Colors.grey[700]),
                  ),
                ],
              ),
            ),
          ],
        ),
      ),
    );
  }

  Color _getPerfumeTypeColor() {
    if (perfumeConcentration <= 8) {
      return Colors.lightBlue;
    } else if (perfumeConcentration <= 15) {
      return Colors.green;
    } else if (perfumeConcentration <= 20) {
      return Colors.orange;
    } else {
      return Colors.red;
    }
  }

  IconData _getPerfumeTypeIcon() {
    if (perfumeConcentration <= 8) {
      return Icons.water_drop_outlined;
    } else if (perfumeConcentration <= 15) {
      return Icons.local_florist;
    } else if (perfumeConcentration <= 20) {
      return Icons.star;
    } else {
      return Icons.diamond;
    }
  }

  String _getPerfumeTypeName() {
    if (perfumeConcentration <= 8) {
      return "Eau de Cologne (5-8%)";
    } else if (perfumeConcentration <= 15) {
      return "Eau de Parfum (8-15%)";
    } else if (perfumeConcentration <= 20) {
      return "Parfum de Toilette (15-20%)";
    } else {
      return "Pure Parfum (20-30%)";
    }
  }

  // ============ TAMBAHKAN VALIDATION UNTUK VOLUME INPUT ============

  bool _isVolumeValid() {
    return selectedVolume > 0 && selectedVolume <= 1000; // Max 1 liter
  }

  Widget _buildVolumeValidationMessage() {
    if (!_isVolumeValid()) {
      return Container(
        margin: EdgeInsets.only(top: 8),
        padding: EdgeInsets.all(8),
        decoration: BoxDecoration(
          color: Colors.red[50],
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: Colors.red[300]!),
        ),
        child: Row(
          children: [
            Icon(Icons.warning, size: 16, color: Colors.red[600]),
            SizedBox(width: 6),
            Text(
              "Volume harus antara 1-1000 mL",
              style: TextStyle(fontSize: 12, color: Colors.red[600]),
            ),
          ],
        ),
      );
    }
    return SizedBox.shrink();
  }

  void _showLogoutDialog(BuildContext context) {
    showDialog(
      context: context,
      builder: (BuildContext context) {
        return AlertDialog(
          title: Row(
            children: [
              Icon(Icons.logout, color: Colors.orange[600]),
              SizedBox(width: 8),
              Text('Konfirmasi Logout'),
            ],
          ),
          content: Text('Apakah Anda yakin ingin keluar dari sistem?'),
          actions: [
            TextButton(
              onPressed: () {
                Navigator.of(context).pop(); // Tutup dialog
              },
              child: Text('Batal'),
            ),
            ElevatedButton(
              onPressed: () {
                Navigator.of(context).pop(); // Tutup dialog
                _logout(context);
              },
              style: ElevatedButton.styleFrom(backgroundColor: Colors.red[600]),
              child: Text('Logout', style: TextStyle(color: Colors.white)),
            ),
          ],
        );
      },
    );
  }

  void _logout(BuildContext context) {
    AuthService.logout(); // Set status logout

    // Kembali ke login page
    Navigator.of(
      context,
    ).pushNamedAndRemoveUntil('/login', (Route<dynamic> route) => false);
  }
}
