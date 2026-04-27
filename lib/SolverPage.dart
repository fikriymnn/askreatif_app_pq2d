import 'dart:math';

import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/amn_selected_groups.dart';
import 'package:askreatif_app/classifyNote.dart';
import 'package:askreatif_app/colorTheme.dart';
import 'package:askreatif_app/group_params.dart';
import 'package:drop_down_search_field/drop_down_search_field.dart';
import 'package:flutter/material.dart';

class SolverPage extends StatefulWidget {
  @override
  _SolverPageState createState() => _SolverPageState();
}

class _SolverPageState extends State<SolverPage> with TickerProviderStateMixin {
  bool isLoading = false;
  double progressValue = 0.0;
  late AnimationController _sidebarController;

  SolventOption selectedSolvent = SolventOption.waterEthanol;
  List<String> selectedCompounds = List.filled(6, '');

  double totalSolventFraction = 0.3;
  double ethanolRatio = 0.4;
  double waterRatio = 0.6;

  List<double> fractions = [];
  String result = "";

  List<Compound> _lastSelectedList = [];
  List<double> _lastFractions = [];
  List<double> _lastOV = [];
  List<Compound> _lastSolvents = [];
  List<double> _lastSolventRatios = [];
  List<double> _lastSolventOVs = [];

  late List<TextEditingController> _compoundControllers;

  @override
  void initState() {
    super.initState();
    _compoundControllers = List.generate(
      6,
      (i) => TextEditingController(text: selectedCompounds[i]),
    );
    _sidebarController = AnimationController(
      vsync: this,
      duration: Duration(milliseconds: 300),
      value: 1.0,
    );
  }

  @override
  void dispose() {
    for (var c in _compoundControllers) c.dispose();
    _sidebarController.dispose();
    super.dispose();
  }

  static const double _minFraction = 0.01;

  void generateRandomFractions() {
    final int filled = selectedCompounds.where((c) => c.isNotEmpty).length;
    if (filled == 0) return;
    Random random = Random();
    final int n = filled;
    final double available = 1.0 - totalSolventFraction;
    final double reserved = _minFraction * n;
    List<double> temp = List.generate(n, (_) => random.nextDouble());
    double sum = temp.reduce((a, b) => a + b);
    fractions =
        temp
            .map((f) => _minFraction + (f / sum) * (available - reserved))
            .toList();
  }

  List<double> calculateImprovedUnifacCoefficients(
    List<Compound> comps,
    List<double> moleFractions,
  ) {
    const double T = 298.15;
    final int nComps = comps.length;
    List<double> gammas = List.filled(nComps, 1.0);

    List<double> r = List.filled(nComps, 0.0);
    List<double> q = List.filled(nComps, 0.0);

    for (int i = 0; i < nComps; i++) {
      Compound comp = comps[i];
      comp.groups.forEach((group, count) {
        if (groupParams.containsKey(group)) {
          r[i] += count * groupParams[group]![0];
          q[i] += count * groupParams[group]![1];
        }
      });
    }

    List<double> volumeFractions = List.filled(nComps, 0.0);
    List<double> surfaceAreaFractions = List.filled(nComps, 0.0);
    double sumRX = 0.0, sumQX = 0.0;
    for (int i = 0; i < nComps; i++) {
      sumRX += r[i] * moleFractions[i];
      sumQX += q[i] * moleFractions[i];
    }
    for (int i = 0; i < nComps; i++) {
      volumeFractions[i] = r[i] * moleFractions[i] / sumRX;
      surfaceAreaFractions[i] = q[i] * moleFractions[i] / sumQX;
    }

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

    List<double> lnGammaR = List.filled(nComps, 0.0);
    Map<String, double> groupFractions = {};
    Map<String, double> totalGroups = {};
    for (int i = 0; i < nComps; i++) {
      comps[i].groups.forEach((group, count) {
        totalGroups[group] =
            (totalGroups[group] ?? 0.0) + count * moleFractions[i];
      });
    }
    double totalGroupCount = totalGroups.values.reduce((a, b) => a + b);
    if (totalGroupCount > 0) {
      totalGroups.forEach((group, count) {
        groupFractions[group] = count / totalGroupCount;
      });
    } else {
      totalGroups.keys.forEach((group) {
        groupFractions[group] = 1.0 / totalGroups.length;
      });
    }

    for (int i = 0; i < nComps; i++) {
      double residualSum = 0.0;
      comps[i].groups.forEach((group1, count1) {
        groupFractions.forEach((group2, fraction) {
          if (amn.containsKey(group1) && amn[group1]!.containsKey(group2)) {
            double a = amn[group1]![group2]!;
            double tau = exp(-a / (8.314 * T));
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

    for (int i = 0; i < nComps; i++) {
      if (!isFinite(lnGammaC[i])) lnGammaC[i] = 0;
      if (!isFinite(lnGammaR[i])) lnGammaR[i] = 0;
      gammas[i] = exp(lnGammaC[i] + lnGammaR[i]);
      if (!isFinite(gammas[i])) gammas[i] = 1.0;
    }

    return gammas;
  }

  bool isFinite(double value) {
    return !value.isNaN && !value.isInfinite;
  }

  List<double> solveLinearSystem(List<List<double>> a, List<double> b) {
    int n = b.length;
    if (n == 0) return [];
    List<List<double>> matrix = List.generate(
      n,
      (i) => List<double>.from(a[i]),
    );
    List<double> rhs = List<double>.from(b);

    for (int i = 0; i < n; i++) {
      int maxRow = i;
      double maxVal = matrix[i][i].abs();
      for (int k = i + 1; k < n; k++) {
        if (matrix[k][i].abs() > maxVal) {
          maxVal = matrix[k][i].abs();
          maxRow = k;
        }
      }
      if (maxVal < 1e-10) throw Exception("Matrix is nearly singular");
      if (maxRow != i) {
        List<double> temp = matrix[i];
        matrix[i] = matrix[maxRow];
        matrix[maxRow] = temp;
        double t = rhs[i];
        rhs[i] = rhs[maxRow];
        rhs[maxRow] = t;
      }
      for (int k = i + 1; k < n; k++) {
        if (matrix[i][i].abs() < 1e-10) continue;
        double factor = matrix[k][i] / matrix[i][i];
        rhs[k] -= factor * rhs[i];
        for (int j = i; j < n; j++) {
          matrix[k][j] -= factor * matrix[i][j];
        }
      }
    }

    List<double> x = List.filled(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
      double sum = 0.0;
      for (int j = i + 1; j < n; j++) {
        sum += matrix[i][j] * x[j];
      }
      if (matrix[i][i].abs() < 1e-10) {
        x[i] = 0.0;
      } else {
        x[i] = (rhs[i] - sum) / matrix[i][i];
      }
      if (x[i].isNaN || x[i].isInfinite) x[i] = 0.0;
    }
    return x;
  }

  List<double> calculateOdorValues(
    List<Compound> selectedComps,
    List<double> moleFractions,
    List<double> activityCoeffs,
  ) {
    const double R = 8.314462618;
    const double T = 298.15;

    List<double> partialPressures = [];
    for (int i = 0; i < selectedComps.length; i++) {
      double Pi = activityCoeffs[i] * moleFractions[i] * selectedComps[i].Psat;
      Pi = (Pi.isNaN || Pi.isInfinite || Pi < 0) ? 1e-6 : Pi;
      partialPressures.add(Pi);
    }

    double Ptotal = partialPressures.reduce((a, b) => a + b);
    if (Ptotal < 1e-10) Ptotal = 1e-10;

    List<double> odorValues = [];
    for (int i = 0; i < selectedComps.length; i++) {
      double yi = partialPressures[i] / Ptotal;
      double threshold = selectedComps[i].Thr;
      if (threshold <= 0 || threshold.isNaN || threshold.isInfinite) {
        threshold = 1.0;
      }
      double ov = (yi * selectedComps[i].MW * Ptotal) / (R * T * threshold);
      ov = (ov.isNaN || ov.isInfinite) ? 1.0 : ov;
      odorValues.add(ov);
    }

    return odorValues;
  }

  double calculateSolventOV(
    Compound solvent,
    double fraction,
    List<Compound> fragranceCompounds,
    List<double> fragranceFractions,
  ) {
    const double R = 8.314462618;
    const double T = 298.15;

    List<Compound> allCompounds = [...fragranceCompounds, solvent];
    List<double> allFractions = [...fragranceFractions, fraction];

    List<double> gammas = calculateImprovedUnifacCoefficients(
      allCompounds,
      allFractions,
    );

    List<double> partialPressures = [];
    for (int i = 0; i < allCompounds.length; i++) {
      double Pi = gammas[i] * allFractions[i] * allCompounds[i].Psat;
      Pi = (Pi.isNaN || Pi.isInfinite || Pi < 0) ? 1e-10 : Pi;
      partialPressures.add(Pi);
    }

    double Ptotal = partialPressures.reduce((a, b) => a + b);
    Ptotal = max(Ptotal, 1e-10);

    double yi = partialPressures.last / Ptotal;
    double threshold = solvent.Thr;
    if (threshold <= 0 || threshold.isNaN || threshold.isInfinite)
      threshold = 1.0;

    double ov = (yi * solvent.MW * Ptotal) / (R * T * threshold);
    return isFinite(ov) ? ov : 0.0;
  }

  void solveEquations() async {
    final filled = selectedCompounds.where((c) => c.isNotEmpty).toList();
    if (filled.length < 2) {
      ScaffoldMessenger.of(context).showSnackBar(
        SnackBar(
          content: Text('Pilih minimal 2 senyawa terlebih dahulu'),
          backgroundColor: AppTheme.green800,
        ),
      );
      return;
    }

    setState(() {
      isLoading = true;
      progressValue = 0.0;
    });

    await Future.delayed(Duration(milliseconds: 50));

    List<Compound> selectedList =
        selectedCompounds
            .where((name) => name.isNotEmpty)
            .map((name) => compounds.firstWhere((c) => c.name == name))
            .toList();

    List<Compound> solvents = [];
    List<double> solventRatios = [];

    if (selectedSolvent case SolventOption.waterEthanol) {
      solvents = [
        compounds.firstWhere((c) => c.name == "Ethanol"),
        compounds.firstWhere((c) => c.name == "Water"),
      ];
      solventRatios = [ethanolRatio, waterRatio];
    } else if (selectedSolvent case SolventOption.water) {
      solvents = [compounds.firstWhere((c) => c.name == "Water")];
      solventRatios = [1.0];
    } else if (selectedSolvent case SolventOption.ethanol) {
      solvents = [compounds.firstWhere((c) => c.name == "Ethanol")];
      solventRatios = [1.0];
    } else if (selectedSolvent case SolventOption.dpg) {
      solvents = [compounds.firstWhere((c) => c.name == "DPG")];
      solventRatios = [1.0];
    }

    final int n = selectedList.length;
    final double fragrantFraction = 1.0 - totalSolventFraction;

    List<double> intrinsicOV =
        selectedList.map((c) => c.Psat / (c.Thr > 0 ? c.Thr : 1.0)).toList();
    List<double> logIntrinsic =
        intrinsicOV.map((v) => log(v > 0 ? v : 1e-10)).toList();
    double logMeanIntrinsic = logIntrinsic.reduce((a, b) => a + b) / n;

    List<double> xInit = List.filled(n, 0.0);
    for (int i = 0; i < n; i++) {
      double correction = -(logIntrinsic[i] - logMeanIntrinsic);
      correction = correction.clamp(-10.0, 10.0);
      xInit[i] = exp(correction);
    }

    double sumInit = xInit.reduce((a, b) => a + b);
    xInit = xInit.map((v) => v / sumInit * fragrantFraction).toList();
    for (int i = 0; i < n; i++) {
      if (xInit[i] < _minFraction) xInit[i] = _minFraction;
    }
    xInit = _renormalize(xInit, fragrantFraction);

    List<double> x = xInit;
    const int maxIter = 800;
    const double alpha = 0.4;
    const double convergenceTol = 0.05;

    List<double> bestX = List.from(x);
    double bestScore = double.infinity;

    for (int iter = 0; iter < maxIter; iter++) {
      List<double> totalFrac = [
        ...x,
        ...solventRatios.map((r) => r * totalSolventFraction),
      ];
      List<Compound> totalComps = [...selectedList, ...solvents];
      List<double> gammas = calculateImprovedUnifacCoefficients(
        totalComps,
        totalFrac,
      );

      List<double> ov = calculateOdorValues(
        selectedList,
        x,
        gammas.sublist(0, n),
      );

      List<double> ovSafe =
          ov
              .map((v) => (v <= 0 || v.isNaN || v.isInfinite) ? 1e-30 : v)
              .toList();
      List<double> logOV = ovSafe.map((v) => log(v)).toList();
      double logMean = logOV.reduce((a, b) => a + b) / n;
      double score = logOV.map((lv) => (lv - logMean).abs()).reduce(max);

      if (score < bestScore) {
        bestScore = score;
        bestX = List.from(x);
      }

      if (iter % 20 == 0) {
        setState(() {
          progressValue = iter / maxIter;
        });
        await Future.delayed(Duration(milliseconds: 1));
      }

      if (score < convergenceTol) break;

      double adaptiveAlpha = score > 5.0 ? 0.7 : (score > 2.0 ? 0.5 : alpha);

      List<double> xNew = List.filled(n, 0.0);
      for (int i = 0; i < n; i++) {
        double correction = -adaptiveAlpha * (logOV[i] - logMean);
        correction = correction.clamp(-3.0, 3.0);
        xNew[i] = x[i] * exp(correction);
      }

      xNew = _renormalize(xNew, fragrantFraction);
      x = xNew;
    }

    List<double> optimizedFractions = bestX;

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

    String newResult = "PERFUME SOLUTION DETAILS:\n";
    newResult += "\n=== MOLE FRACTIONS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult +=
          "${selectedList[i].name}: ${optimizedFractions[i].toStringAsFixed(4)}\n";
    }
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${(solventRatios[i] * totalSolventFraction).toStringAsFixed(4)}\n";
    }
    newResult += "\n=== ACTIVITY COEFFICIENTS ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${gammas[i].toStringAsFixed(4)}\n";
    }
    for (int i = 0; i < solvents.length; i++) {
      newResult +=
          "${solvents[i].name}: ${gammas[selectedList.length + i].toStringAsFixed(4)}\n";
    }
    newResult += "\n=== ODOR VALUES ===\n";
    for (int i = 0; i < selectedList.length; i++) {
      newResult += "${selectedList[i].name}: ${ov[i].toStringAsFixed(4)}\n";
    }
    newResult += "\n=== SOLVENT ODOR VALUES ===\n";
    for (int i = 0; i < solvents.length; i++) {
      double solventOV = calculateSolventOV(
        solvents[i],
        solventRatios[i] * totalSolventFraction,
        selectedList,
        optimizedFractions,
      );
      newResult += "${solvents[i].name}: ${solventOV.toStringAsFixed(6)}\n";
    }

    List<double> ovSafe =
        ov.map((v) => (v <= 0 || v.isNaN || v.isInfinite) ? 1e-30 : v).toList();
    List<double> logOVFinal = ovSafe.map((v) => log(v)).toList();
    double logMeanFinal = logOVFinal.reduce((a, b) => a + b) / ov.length;
    double finalScore = logOVFinal
        .map((lv) => (lv - logMeanFinal).abs())
        .reduce(max);

    if (finalScore > 2.0) {
      newResult +=
          "\n⚠️ Beberapa senyawa memiliki OV intrinsik yang sangat berbeda\n"
          "sehingga tidak bisa diseimbangkan sempurna.\n";
    }

    List<double> solventOVsList = [];
    for (int i = 0; i < solvents.length; i++) {
      double sOV = calculateSolventOV(
        solvents[i],
        solventRatios[i] * totalSolventFraction,
        selectedList,
        optimizedFractions,
      );
      solventOVsList.add(sOV);
    }

    _lastSelectedList = selectedList;
    _lastFractions = optimizedFractions;
    _lastOV = ov;
    _lastSolvents = solvents;
    _lastSolventRatios = solventRatios;
    _lastSolventOVs = solventOVsList;

    setState(() {
      result = newResult;
      isLoading = false;
      progressValue = 1.0;
    });
  }

  List<double> _renormalize(List<double> x, double targetSum) {
    final int n = x.length;
    List<double> result = List.from(x);
    for (int i = 0; i < n; i++) {
      if (result[i] < _minFraction) result[i] = _minFraction;
    }
    double clampedSum = 0, freeSum = 0;
    for (int i = 0; i < n; i++) {
      if (result[i] <= _minFraction)
        clampedSum += _minFraction;
      else
        freeSum += result[i];
    }
    double remainForFree = targetSum - clampedSum;
    if (remainForFree <= 0) return List.filled(n, targetSum / n);
    if (freeSum > 0) {
      for (int i = 0; i < n; i++) {
        if (result[i] > _minFraction) {
          result[i] = result[i] / freeSum * remainForFree;
        }
      }
    }
    return result;
  }

  List<String> _getSuggestions(String query) {
    final allFragrances =
        compounds
            .where((c) => !["Water", "Ethanol", "DPG"].contains(c.name))
            .map((c) => c.name)
            .toList();
    if (query.isEmpty) return allFragrances;
    return allFragrances
        .where((name) => name.toLowerCase().contains(query.toLowerCase()))
        .toList();
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      backgroundColor: AppTheme.bg,
      body: Column(
        children: [
          _buildTopBar(),
          Expanded(
            child: Row(
              children: [
                // LEFT SIDEBAR
                _buildSidebar(),
                // MAIN CONTENT
                Expanded(child: _buildMainContent()),
              ],
            ),
          ),
        ],
      ),
    );
  }

  Widget _buildTopBar() {
    return Container(
      height: 56,
      decoration: BoxDecoration(
        color: AppTheme.bgSurface,
        border: Border(bottom: BorderSide(color: AppTheme.border, width: 1)),
      ),
      child: Padding(
        padding: EdgeInsets.symmetric(horizontal: 24),
        child: Row(
          children: [
            // Logo
            Container(
              width: 32,
              height: 32,
              decoration: BoxDecoration(
                shape: BoxShape.circle,
                gradient: AppTheme.primaryGradient,
                border: Border.all(color: AppTheme.green400.withOpacity(0.5)),
              ),
              child: Center(
                child: Text(
                  'ask',
                  style: TextStyle(
                    fontFamily: AppTheme.fontSerif,
                    color: Colors.white,
                    fontSize: 14,
                    fontStyle: FontStyle.italic,
                  ),
                ),
              ),
            ),
            SizedBox(width: 12),
            Text(
              'ASKREATIF Essential Oil',
              style: TextStyle(
                fontFamily: AppTheme.fontSerif,
                color: AppTheme.textPrimary,
                fontSize: 15,
                letterSpacing: 1,
              ),
            ),
            Container(
              margin: EdgeInsets.symmetric(horizontal: 16),
              width: 1,
              height: 20,
              color: AppTheme.border,
            ),
            Text(
              'POS Calculator',
              style: TextStyle(color: AppTheme.textSecondary, fontSize: 13),
            ),
            Spacer(),
            if (isLoading) ...[
              SizedBox(
                width: 200,
                child: Column(
                  mainAxisAlignment: MainAxisAlignment.center,
                  children: [
                    ClipRRect(
                      borderRadius: BorderRadius.circular(2),
                      child: LinearProgressIndicator(
                        value: progressValue,
                        backgroundColor: AppTheme.border,
                        color: AppTheme.green400,
                        minHeight: 3,
                      ),
                    ),
                    SizedBox(height: 3),
                    Text(
                      'Menghitung... ${(progressValue * 100).toStringAsFixed(0)}%',
                      style: TextStyle(
                        color: AppTheme.textSecondary,
                        fontSize: 10,
                      ),
                    ),
                  ],
                ),
              ),
              SizedBox(width: 16),
            ],
            _StatusDot(label: 'UNIFAC Ready', active: !isLoading),
          ],
        ),
      ),
    );
  }

  Widget _buildSidebar() {
    return Container(
      width: 300,
      decoration: BoxDecoration(
        color: AppTheme.bgSurface,
        border: Border(right: BorderSide(color: AppTheme.border, width: 1)),
      ),
      child: Column(
        children: [
          // Sidebar header
          Container(
            padding: EdgeInsets.all(20),
            child: Row(
              children: [
                Icon(Icons.tune, color: AppTheme.green400, size: 16),
                SizedBox(width: 8),
                Text(
                  'PARAMETER',
                  style: TextStyle(
                    color: AppTheme.green400,
                    fontSize: 11,
                    letterSpacing: 2,
                    fontWeight: FontWeight.w600,
                  ),
                ),
              ],
            ),
          ),
          Divider(color: AppTheme.border, height: 1),

          Expanded(
            child: SingleChildScrollView(
              padding: EdgeInsets.all(16),
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  // Compound selection
                  _SectionLabel(
                    label: 'SENYAWA PEWANGI',
                    icon: Icons.science_outlined,
                  ),
                  SizedBox(height: 10),
                  for (int i = 0; i < 6; i++) ...[
                    _buildCompoundSearchField(i),
                    SizedBox(height: 8),
                  ],

                  SizedBox(height: 20),
                  _SectionLabel(
                    label: 'PELARUT',
                    icon: Icons.water_drop_outlined,
                  ),
                  SizedBox(height: 10),

                  // Solvent dropdown
                  Container(
                    padding: EdgeInsets.symmetric(horizontal: 12, vertical: 4),
                    decoration: BoxDecoration(
                      color: AppTheme.bgCard,
                      borderRadius: BorderRadius.circular(4),
                      border: Border.all(color: AppTheme.border),
                    ),
                    child: DropdownButton<SolventOption>(
                      isExpanded: true,
                      value: selectedSolvent,
                      dropdownColor: AppTheme.bgCard,
                      underline: SizedBox(),
                      style: TextStyle(
                        color: AppTheme.textPrimary,
                        fontSize: 13,
                      ),
                      iconEnabledColor: AppTheme.green400,
                      onChanged: (newValue) {
                        setState(() {
                          selectedSolvent = newValue!;
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
                  ),

                  SizedBox(height: 16),

                  // Solvent config card
                  Container(
                    padding: EdgeInsets.all(14),
                    decoration: BoxDecoration(
                      color: AppTheme.bgCard,
                      borderRadius: BorderRadius.circular(4),
                      border: Border.all(color: AppTheme.border),
                    ),
                    child: Column(
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        if (selectedSolvent == SolventOption.waterEthanol) ...[
                          Text(
                            'Rasio Water : Ethanol',
                            style: TextStyle(
                              color: AppTheme.textSecondary,
                              fontSize: 11,
                              letterSpacing: 0.5,
                            ),
                          ),
                          SizedBox(height: 8),
                          _WaterEthanolBar(
                            waterRatio: waterRatio,
                            ethanolRatio: ethanolRatio,
                          ),
                          SliderTheme(
                            data: SliderThemeData(
                              trackHeight: 2,
                              activeTrackColor: AppTheme.green400,
                              inactiveTrackColor: AppTheme.border,
                              thumbColor: AppTheme.green400,
                              overlayColor: AppTheme.green400.withOpacity(0.1),
                              thumbShape: RoundSliderThumbShape(
                                enabledThumbRadius: 6,
                              ),
                            ),
                            child: Slider(
                              value: waterRatio,
                              min: 0.1,
                              max: 0.9,
                              divisions: 8,
                              onChanged: (value) {
                                setState(() {
                                  waterRatio = value;
                                  ethanolRatio = 1.0 - value;
                                });
                              },
                            ),
                          ),
                          SizedBox(height: 8),
                        ],

                        Text(
                          'Total Solvent Fraction',
                          style: TextStyle(
                            color: AppTheme.textSecondary,
                            fontSize: 11,
                            letterSpacing: 0.5,
                          ),
                        ),
                        SizedBox(height: 4),
                        SliderTheme(
                          data: SliderThemeData(
                            trackHeight: 2,
                            activeTrackColor: AppTheme.green600,
                            inactiveTrackColor: AppTheme.border,
                            thumbColor: AppTheme.green600,
                            overlayColor: AppTheme.green600.withOpacity(0.1),
                            thumbShape: RoundSliderThumbShape(
                              enabledThumbRadius: 6,
                            ),
                          ),
                          child: Slider(
                            value: totalSolventFraction,
                            min: 0.1,
                            max: 0.5,
                            divisions: 40,
                            onChanged: (value) {
                              setState(() {
                                totalSolventFraction = value;
                                generateRandomFractions();
                              });
                            },
                          ),
                        ),
                        Row(
                          mainAxisAlignment: MainAxisAlignment.spaceBetween,
                          children: [
                            Text(
                              'Fraksi Pelarut',
                              style: TextStyle(
                                color: AppTheme.textMuted,
                                fontSize: 10,
                              ),
                            ),
                            Container(
                              padding: EdgeInsets.symmetric(
                                horizontal: 8,
                                vertical: 2,
                              ),
                              decoration: BoxDecoration(
                                color: AppTheme.green900.withOpacity(0.5),
                                borderRadius: BorderRadius.circular(3),
                                border: Border.all(color: AppTheme.green800),
                              ),
                              child: Text(
                                totalSolventFraction.toStringAsFixed(2),
                                style: TextStyle(
                                  color: AppTheme.green200,
                                  fontSize: 11,
                                  fontFamily: 'Courier',
                                ),
                              ),
                            ),
                          ],
                        ),
                      ],
                    ),
                  ),

                  SizedBox(height: 24),

                  // Calculate button
                  _CalcButton(isLoading: isLoading, onPressed: solveEquations),
                ],
              ),
            ),
          ),
        ],
      ),
    );
  }

  Widget _buildMainContent() {
    return Container(
      color: AppTheme.bg,
      child:
          result.isEmpty && !isLoading
              ? _buildEmptyState()
              : SingleChildScrollView(
                padding: EdgeInsets.all(28),
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    if (result.isNotEmpty) ...[
                      _buildResultsSection(),
                      SizedBox(height: 24),
                      if (_lastSelectedList.isNotEmpty)
                        buildNotesAndTernary(
                          _lastSelectedList,
                          _lastFractions,
                          _lastOV,
                        ),
                      SizedBox(height: 24),
                      buildOctonaryDiagram(
                        _lastSelectedList,
                        _lastFractions,
                        _lastOV,
                        _lastSolvents,
                        _lastSolventRatios,
                        _lastSolventOVs,
                      ),
                      SizedBox(height: 40),
                    ],
                  ],
                ),
              ),
    );
  }

  Widget _buildEmptyState() {
    return Center(
      child: Column(
        mainAxisAlignment: MainAxisAlignment.center,
        children: [
          Container(
            width: 80,
            height: 80,
            decoration: BoxDecoration(
              shape: BoxShape.circle,
              color: AppTheme.bgCard,
              border: Border.all(color: AppTheme.border),
            ),
            child: Icon(
              Icons.science_outlined,
              color: AppTheme.textMuted,
              size: 32,
            ),
          ),
          SizedBox(height: 20),
          Text(
            'Siap untuk Kalkulasi',
            style: TextStyle(
              fontFamily: AppTheme.fontSerif,
              color: AppTheme.textPrimary,
              fontSize: 22,
            ),
          ),
          SizedBox(height: 8),
          Text(
            'Pilih senyawa dan konfigurasi pelarut di panel kiri,\nlalu klik "Hitung" untuk memulai optimasi.',
            textAlign: TextAlign.center,
            style: TextStyle(
              color: AppTheme.textMuted,
              fontSize: 13,
              height: 1.6,
            ),
          ),
          SizedBox(height: 32),
          // Mini feature pills
          Wrap(
            spacing: 8,
            runSpacing: 8,
            alignment: WrapAlignment.center,
            children: [
              _Pill(label: 'UNIFAC Model'),
              _Pill(label: 'Log-space Balancing'),
              _Pill(label: 'Ternary Diagram'),
              _Pill(label: 'Octonary System'),
            ],
          ),
        ],
      ),
    );
  }

  Widget _buildResultsSection() {
    final lines = result.split('\n').where((l) => l.trim().isNotEmpty).toList();

    // Parse sections
    Map<String, List<String>> sections = {};
    String currentSection = '';
    for (String line in lines) {
      if (line.startsWith('===')) {
        currentSection = line.replaceAll('=', '').trim();
        sections[currentSection] = [];
      } else if (currentSection.isNotEmpty && !line.startsWith('PERFUME')) {
        sections[currentSection]!.add(line);
      }
    }

    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Row(
          children: [
            Icon(
              Icons.check_circle_outline,
              color: AppTheme.green400,
              size: 18,
            ),
            SizedBox(width: 8),
            Text(
              'Hasil Optimasi',
              style: TextStyle(
                fontFamily: AppTheme.fontSerif,
                color: AppTheme.textPrimary,
                fontSize: 20,
              ),
            ),
          ],
        ),
        SizedBox(height: 16),
        // Result cards in a row
        Row(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            for (var entry in sections.entries)
              Expanded(
                child: Padding(
                  padding: EdgeInsets.only(right: 12),
                  child: _ResultCard(title: entry.key, lines: entry.value),
                ),
              ),
          ],
        ),
      ],
    );
  }

  Widget _buildCompoundSearchField(int index) {
    return DropDownSearchFormField<String>(
      textFieldConfiguration: TextFieldConfiguration(
        controller: _compoundControllers[index],
        style: TextStyle(color: AppTheme.textPrimary, fontSize: 13),
        decoration: InputDecoration(
          labelText: 'Senyawa ${index + 1}',
          labelStyle: TextStyle(color: AppTheme.textMuted, fontSize: 12),
          hintText: 'Cari senyawa...',
          hintStyle: TextStyle(color: AppTheme.textMuted, fontSize: 12),
          filled: true,
          fillColor: AppTheme.bgCard,
          contentPadding: EdgeInsets.symmetric(horizontal: 12, vertical: 10),
          border: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.border),
          ),
          enabledBorder: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.border),
          ),
          focusedBorder: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.green400),
          ),
          prefixIcon: Icon(Icons.search, color: AppTheme.textMuted, size: 16),
          suffixIcon:
              selectedCompounds[index].isNotEmpty
                  ? Icon(Icons.check_circle, color: AppTheme.green400, size: 16)
                  : Icon(
                    Icons.arrow_drop_down,
                    color: AppTheme.textMuted,
                    size: 20,
                  ),
        ),
      ),
      suggestionsCallback: (pattern) => _getSuggestions(pattern),
      itemBuilder: (context, String suggestion) {
        final isSelected = suggestion == selectedCompounds[index];
        return Container(
          color:
              isSelected ? AppTheme.green900.withOpacity(0.5) : AppTheme.bgCard,
          child: ListTile(
            dense: true,
            title: Text(
              suggestion,
              style: TextStyle(
                color: isSelected ? AppTheme.green200 : AppTheme.textPrimary,
                fontSize: 13,
              ),
            ),
            trailing:
                isSelected
                    ? Icon(Icons.check, color: AppTheme.green400, size: 14)
                    : null,
          ),
        );
      },
      transitionBuilder:
          (context, suggestionsBox, controller) => suggestionsBox,
      onSuggestionSelected: (String suggestion) {
        setState(() {
          selectedCompounds[index] = suggestion;
          _compoundControllers[index].text = suggestion;
          generateRandomFractions();
        });
      },
      hideOnEmpty: true,
      suggestionsBoxDecoration: SuggestionsBoxDecoration(
        elevation: 8,
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        constraints: BoxConstraints(maxHeight: 220),
        //border: Border.all(color: AppTheme.border),
      ),
    );
  }
}

// =========================================================
// SMALL COMPONENTS
// =========================================================

class _StatusDot extends StatelessWidget {
  final String label;
  final bool active;
  const _StatusDot({required this.label, required this.active});

  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        Container(
          width: 7,
          height: 7,
          decoration: BoxDecoration(
            shape: BoxShape.circle,
            color: active ? AppTheme.green400 : Colors.orange,
            boxShadow: [
              BoxShadow(
                color: (active ? AppTheme.green400 : Colors.orange).withOpacity(
                  0.5,
                ),
                blurRadius: 4,
              ),
            ],
          ),
        ),
        SizedBox(width: 6),
        Text(
          label,
          style: TextStyle(color: AppTheme.textSecondary, fontSize: 11),
        ),
      ],
    );
  }
}

class _SectionLabel extends StatelessWidget {
  final String label;
  final IconData icon;
  const _SectionLabel({required this.label, required this.icon});

  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        Icon(icon, color: AppTheme.green400, size: 13),
        SizedBox(width: 6),
        Text(
          label,
          style: TextStyle(
            color: AppTheme.green400,
            fontSize: 10,
            letterSpacing: 2,
            fontWeight: FontWeight.w600,
          ),
        ),
      ],
    );
  }
}

class _WaterEthanolBar extends StatelessWidget {
  final double waterRatio;
  final double ethanolRatio;
  const _WaterEthanolBar({
    required this.waterRatio,
    required this.ethanolRatio,
  });

  @override
  Widget build(BuildContext context) {
    return ClipRRect(
      borderRadius: BorderRadius.circular(2),
      child: Row(
        children: [
          Expanded(
            flex: (waterRatio * 100).round(),
            child: Container(
              height: 18,
              color: Color(0xFF1565C0).withOpacity(0.7),
              child: Center(
                child: Text(
                  'W ${(waterRatio * 100).round()}%',
                  style: TextStyle(fontSize: 10, color: Colors.white70),
                ),
              ),
            ),
          ),
          Expanded(
            flex: (ethanolRatio * 100).round(),
            child: Container(
              height: 18,
              color: AppTheme.green800.withOpacity(0.8),
              child: Center(
                child: Text(
                  'E ${(ethanolRatio * 100).round()}%',
                  style: TextStyle(fontSize: 10, color: Colors.white70),
                ),
              ),
            ),
          ),
        ],
      ),
    );
  }
}

class _CalcButton extends StatefulWidget {
  final bool isLoading;
  final VoidCallback onPressed;
  const _CalcButton({required this.isLoading, required this.onPressed});

  @override
  __CalcButtonState createState() => __CalcButtonState();
}

class __CalcButtonState extends State<_CalcButton> {
  bool _hovered = false;

  @override
  Widget build(BuildContext context) {
    return MouseRegion(
      onEnter: (_) => setState(() => _hovered = true),
      onExit: (_) => setState(() => _hovered = false),
      child: GestureDetector(
        onTap: widget.isLoading ? null : widget.onPressed,
        child: AnimatedContainer(
          duration: Duration(milliseconds: 200),
          width: double.infinity,
          height: 44,
          decoration: BoxDecoration(
            gradient: LinearGradient(
              colors:
                  widget.isLoading
                      ? [AppTheme.border, AppTheme.border]
                      : _hovered
                      ? [AppTheme.green600, AppTheme.green800]
                      : [AppTheme.green800, AppTheme.green900],
            ),
            borderRadius: BorderRadius.circular(4),
            boxShadow:
                _hovered && !widget.isLoading
                    ? [
                      BoxShadow(
                        color: AppTheme.green800.withOpacity(0.4),
                        blurRadius: 12,
                      ),
                    ]
                    : [],
            border: Border.all(
              color:
                  widget.isLoading
                      ? AppTheme.border
                      : AppTheme.green600.withOpacity(_hovered ? 0.8 : 0.4),
            ),
          ),
          child: Row(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              Icon(
                widget.isLoading ? Icons.hourglass_top : Icons.play_arrow,
                color: widget.isLoading ? AppTheme.textMuted : Colors.white,
                size: 16,
              ),
              SizedBox(width: 8),
              Text(
                widget.isLoading ? 'Menghitung...' : 'HITUNG OPTIMAL',
                style: TextStyle(
                  color: widget.isLoading ? AppTheme.textMuted : Colors.white,
                  fontSize: 12,
                  fontWeight: FontWeight.w600,
                  letterSpacing: 2,
                ),
              ),
            ],
          ),
        ),
      ),
    );
  }
}

class _ResultCard extends StatelessWidget {
  final String title;
  final List<String> lines;
  const _ResultCard({required this.title, required this.lines});

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: AppTheme.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            title,
            style: TextStyle(
              color: AppTheme.green400,
              fontSize: 10,
              letterSpacing: 1.5,
              fontWeight: FontWeight.w600,
            ),
          ),
          Divider(color: AppTheme.border, height: 16),
          ...lines.map((line) {
            final parts = line.split(':');
            if (parts.length >= 2) {
              return Padding(
                padding: EdgeInsets.only(bottom: 6),
                child: Row(
                  children: [
                    Expanded(
                      flex: 3,
                      child: Text(
                        parts[0].trim(),
                        style: TextStyle(
                          color: AppTheme.textSecondary,
                          fontSize: 11,
                        ),
                      ),
                    ),
                    Text(
                      parts.sublist(1).join(':').trim(),
                      style: TextStyle(
                        color: AppTheme.textPrimary,
                        fontSize: 11,
                        fontFamily: 'Courier',
                      ),
                    ),
                  ],
                ),
              );
            }
            return Text(
              line,
              style: TextStyle(
                color:
                    line.startsWith('⚠')
                        ? Colors.amber[300]
                        : AppTheme.textMuted,
                fontSize: 11,
              ),
            );
          }),
        ],
      ),
    );
  }
}

class _Pill extends StatelessWidget {
  final String label;
  const _Pill({required this.label});

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 12, vertical: 5),
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(20),
        border: Border.all(color: AppTheme.border),
      ),
      child: Text(
        label,
        style: TextStyle(color: AppTheme.textSecondary, fontSize: 11),
      ),
    );
  }
}
