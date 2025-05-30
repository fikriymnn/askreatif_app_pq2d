import 'dart:math';
import 'package:flutter/material.dart';
import 'package:fl_chart/fl_chart.dart';

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
  final double R;
  final double Q;
  final Map<String, int> groups;

  Compound({
    required this.name,
    required this.MW,
    required this.Psat,
    required this.Thr,
    required this.R,
    required this.Q,
    required this.groups,
  });
}

// Expanded UNIFAC interaction parameters (amn) - using values from MATLAB
final Map<String, Map<String, double>> amn = {
  "CH3": {
    "CH3": 0.0,
    "CH2": 0.0,
    "OH": 476.4,
    "C=O": 986.5,
    "C6H5": 1318.0,
    "NH2": 1500.0,
    "H2O": 3000.0,
  },
  "CH2": {
    "CH3": 0.0,
    "CH2": 0.0,
    "OH": 1318.0,
    "C=O": 986.5,
    "C6H5": 1318.0,
    "NH2": 1500.0,
    "H2O": 3000.0,
  },
  "OH": {
    "CH3": 476.4,
    "CH2": 1318.0,
    "OH": 0.0,
    "C=O": -229.1,
    "C6H5": -229.1,
    "NH2": 0.0,
    "H2O": 0.0,
  },
  "C=O": {
    "CH3": 986.5,
    "CH2": 986.5,
    "OH": -229.1,
    "C=O": 0.0,
    "C6H5": 1200.0,
    "NH2": 1200.0,
    "H2O": 0.0,
  },
  "C6H5": {
    "CH3": 1318.0,
    "CH2": 1318.0,
    "OH": -229.1,
    "C=O": 1200.0,
    "C6H5": 0.0,
    "NH2": 0.0,
    "H2O": 0.0,
  },
  "NH2": {
    "CH3": 1500.0,
    "CH2": 1500.0,
    "OH": 0.0,
    "C=O": 1200.0,
    "C6H5": 0.0,
    "NH2": 0.0,
    "H2O": 0.0,
  },
  "H2O": {
    "CH3": 3000.0,
    "CH2": 3000.0,
    "OH": 0.0,
    "C=O": 0.0,
    "C6H5": 0.0,
    "NH2": 0.0,
    "H2O": 0.0,
  },
};

// UNIFAC group parameters - similar to MATLAB R and Q values
final Map<String, List<double>> groupParams = {
  "CH3": [0.9011, 0.848],
  "CH2": [0.6744, 0.540],
  "OH": [1.0000, 1.200],
  "C=O": [1.6724, 1.488],
  "C6H5": [3.1878, 2.400],
  "NH2": [0.6948, 0.696],
  "H2O": [0.9200, 1.400],
};

final List<Compound> compounds = [
  Compound(
    name: "Geraniol",
    MW: 154.25,
    Psat: 8.60, // Using values closer to MATLAB
    Thr: 0.018,
    R: 7.5037,
    Q: 6.716,
    groups: {"CH3": 2, "CH2": 4, "OH": 1},
  ),
  Compound(
    name: "Romandolide",
    MW: 156.26,
    Psat: 7.999343,
    Thr: 0.0021,
    R: 11.0885,
    Q: 9.328,
    groups: {"CH3": 3, "CH2": 6, "C=O": 1},
  ),
  Compound(
    name: "Cashmeran",
    MW: 192.30,
    Psat: 2.5331254,
    Thr: 0.0004,
    R: 8.8464,
    Q: 7.213,
    groups: {"CH3": 3, "C6H5": 1},
  ),
  Compound(
    name: "Javanol",
    MW: 222.37,
    Psat: 1.333224,
    Thr: 0.000027,
    R: 9.9756,
    Q: 7.976,
    groups: {"CH3": 2, "OH": 1, "C6H5": 1},
  ),
  Compound(
    name: "Benzyl Acetate",
    MW: 284.44,
    Psat: 1.33,
    Thr: 0.0000017,
    R: 5.5992,
    Q: 4.388,
    groups: {"CH3": 1, "C=O": 1},
  ),
  Compound(
    name: "Vanillin",
    MW: 284.44,
    Psat: 1.33,
    Thr: 0.0000017,
    R: 5.3625,
    Q: 4.156,
    groups: {"CH3": 1, "OH": 1, "C=O": 1},
  ),
  Compound(
    name: "Water",
    MW: 18.015,
    Psat: 3170,
    Thr: double.infinity,
    R: 0.92,
    Q: 1.4,
    groups: {"H2O": 1},
  ),
  Compound(
    name: "Ethanol",
    MW: 46.07,
    Psat: 7270,
    Thr: 0.0553,
    R: 2.5755,
    Q: 2.588,
    groups: {"CH3": 1, "OH": 1},
  ),
];

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
    List<double> tempFractions = List.generate(6, (_) => random.nextDouble());
    double tempSum = tempFractions.reduce((a, b) => a + b);
    tempFractions = tempFractions
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
      lnGammaC[i] = 1 -
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
              residualSum += count1 *
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
              double fPerturbed = odorValuesPerturbed[i + 1] -
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

      // Ensure threshold value is not zero
      double threshold = selectedComps[i].Thr;
      if (threshold < 1e-10) threshold = 1e-10;

      double ov = (yi * selectedComps[i].MW * Ptotal) / (R * T * threshold);

      // Ensure OV is finite
      ov = (ov.isNaN || ov.isInfinite) ? 1.0 : ov;
      odorValues.add(ov);
    }

    return odorValues;
  }

  void solveEquations() async {
    setState(() {
      isLoading = true;
      progressValue = 0.0;
    });

    List<Compound> selectedList = selectedCompounds
        .map((name) => compounds.firstWhere((c) => c.name == name))
        .toList();

    Compound ethanol = compounds.firstWhere((c) => c.name == "Ethanol");
    Compound water = compounds.firstWhere((c) => c.name == "Water");

    // Generate random initial guess with time-based seed
    Random random = Random(DateTime.now().millisecondsSinceEpoch);
    List<double> initialGuess = List.generate(
      selectedList.length,
      (_) => random.nextDouble(),
    );
    double sum = initialGuess.reduce((a, b) => a + b);
    initialGuess = initialGuess
        .map((f) => f / sum * (1.0 - totalSolventFraction))
        .toList();

    // THE KEY CHANGE: Generate random target ratios instead of equal odor values
    // These ratios determine the relative odor strengths between compounds
    List<double> targetRatios = List.generate(
      selectedList.length - 1,
      (_) => 0.8 + 0.4 * random.nextDouble(), // Ratios between 0.8 and 1.2
    );

    // Try Newton-Raphson with multiple initial guesses if needed
    List<double> optimizedFractions = [];
    bool success = false;
    String errorMessage = "";

    // Maximum number of attempts with different initial guesses
    const int maxAttempts = 5;

    for (int attempt = 0; attempt < maxAttempts; attempt++) {
      try {
        if (attempt > 0) {
          // Generate a new random initial guess for this attempt
          initialGuess = List.generate(
            selectedList.length,
            (_) => random.nextDouble(),
          );
          sum = initialGuess.reduce((a, b) => a + b);
          initialGuess = initialGuess
              .map((f) => f / sum * (1.0 - totalSolventFraction))
              .toList();

          // For later attempts, also vary the target ratios
          if (attempt >= 3) {
            targetRatios = List.generate(
              selectedList.length - 1,
              (_) => 0.7 + 0.6 * random.nextDouble(), // Wider range: 0.7-1.3
            );
          }

          // Update progress indicator
          setState(() {
            progressValue = attempt / maxAttempts;
          });
        }

        // Attempt Newton-Raphson with this initial guess and target ratios
        optimizedFractions = solveForCustomOdorValuesWithInitialGuess(
          selectedList,
          totalSolventFraction,
          initialGuess,
          targetRatios,
        );

        // If we get here without an exception, it worked
        success = true;
        break;
      } catch (e) {
        errorMessage = e.toString();
        print("Attempt $attempt failed: $errorMessage");
        // Continue to next attempt with a new initial guess
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

    // Now recalculate all the values based on the optimized fractions
    List<double> totalFractions = [
      ...optimizedFractions,
      ethanolRatio * totalSolventFraction,
      waterRatio * totalSolventFraction,
    ];

    List<Compound> totalCompounds = [...selectedList, ethanol, water];

    List<double> gammas = calculateImprovedUnifacCoefficients(
      totalCompounds,
      totalFractions,
    );
    List<double> ov = calculateOdorValues(
      selectedList,
      optimizedFractions,
      gammas.sublist(0, selectedList.length),
    );

    // Calculate standard deviation to check how well we met our target ratios
    double firstOV = ov[0];
    List<double> actualRatios = ov.sublist(1).map((v) => v / firstOV).toList();
    List<double> ratioDiffs = List.generate(
      targetRatios.length,
      (i) => (actualRatios[i] - targetRatios[i]).abs(),
    );
    double avgRatioDiff =
        ratioDiffs.reduce((a, b) => a + b) / ratioDiffs.length;

    // Calculate solvent OVs
    double ethanolOV = calculateSolventOV(
      ethanol,
      ethanolRatio * totalSolventFraction,
      selectedList,
      optimizedFractions,
    );
    double waterOV = calculateSolventOV(
      water,
      waterRatio * totalSolventFraction,
      selectedList,
      optimizedFractions,
    );

    // Update result string with comprehensive information
    String newResult =
        "PERFUME SOLUTION DETAILS (Modified Newton-Raphson Method):\n";
    newResult += "\n=== MOLE FRACTIONS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult +=
          "${selectedList[i].name}: ${optimizedFractions[i].toStringAsFixed(4)}\n";
    }
    newResult +=
        "Ethanol: ${(ethanolRatio * totalSolventFraction).toStringAsFixed(4)}\n";
    newResult +=
        "Water: ${(waterRatio * totalSolventFraction).toStringAsFixed(4)}\n";

    newResult += "\n=== ACTIVITY COEFFICIENTS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${gammas[i].toStringAsFixed(4)}\n";
    }
    newResult += "Ethanol: ${gammas[selectedList.length].toStringAsFixed(4)}\n";
    newResult +=
        "Water: ${gammas[selectedList.length + 1].toStringAsFixed(4)}\n";

    newResult += "\n=== ODOR VALUES ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${ov[i].toStringAsFixed(4)}\n";
    }

    // newResult += "\n=== TARGET RATIOS ===\n";
    // newResult += "Base Compound: ${selectedList[0].name}\n";
    // for (int i = 0; i < targetRatios.length; i++) {
    //   newResult +=
    //       "${selectedList[i + 1].name}: ${targetRatios[i].toStringAsFixed(4)}\n";
    // }
    // newResult +=
    //     "Average Ratio Difference: ${avgRatioDiff.toStringAsFixed(6)}\n";

    newResult += "\n=== SOLVENT ODOR VALUES ===\n";
    newResult += "Ethanol: ${ethanolOV.toStringAsFixed(4)}\n";
    newResult += "Water: ${waterOV.toStringAsFixed(4)}\n";

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

    // Combine fragrance compounds + solvent
    List<Compound> allCompounds = [...fragranceCompounds, solvent];
    List<double> allFractions = [...fragranceFractions, fraction];

    // Calculate gamma for all
    List<double> gammas = calculateImprovedUnifacCoefficients(
      allCompounds,
      allFractions,
    );
    double gammaSolvent = gammas.last;

    // Calculate partial pressures for all components
    List<double> partialPressures = [];
    for (int i = 0; i < allCompounds.length; i++) {
      double Pi = gammas[i] * allFractions[i] * allCompounds[i].Psat;
      // Handle potential non-finite values
      Pi = isFinite(Pi) ? Pi : 1e-10;
      partialPressures.add(Pi);
    }

    // Calculate total pressure
    double Ptotal = partialPressures.reduce((a, b) => a + b);
    // Prevent division by zero
    Ptotal = max(Ptotal, 1e-10);

    // Calculate partial pressure and vapor fraction of solvent
    double vaporPressure = gammaSolvent * fraction * solvent.Psat;
    double yi = vaporPressure / Ptotal;

    // Calculate and return OV, handling potential non-finite values
    double ov = (yi * solvent.MW * Ptotal) / (R * T * max(solvent.Thr, 1e-10));
    return isFinite(ov) ? ov : 0.0;
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
                    items: compounds
                        .where(
                          (c) => !["Water", "Ethanol"].contains(c.name),
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

              // Solvent information
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
                      Text("Water : Ethanol Ratio = 60 : 40"),
                      SizedBox(height: 12),
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
                      Text(
                        "Ethanol: ${(ethanolRatio * totalSolventFraction).toStringAsFixed(3)} | Water: ${(waterRatio * totalSolventFraction).toStringAsFixed(3)}",
                        style: TextStyle(fontStyle: FontStyle.italic),
                      ),
                    ],
                  ),
                ),
              ),

              SizedBox(height: 24),

              // Calculate button
              Center(
                child: isLoading
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
                            style:
                                TextStyle(fontFamily: 'Courier', fontSize: 14),
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
