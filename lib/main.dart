import 'dart:math';
import 'dart:math' as math;
import 'package:flutter/material.dart';
import 'package:fl_chart/fl_chart.dart';

void main() {
  runApp(MyApp());
}

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(debugShowCheckedModeBanner: false, home: SolverPage());
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

final Map<String, Map<String, double>> anm = {
  "CH3": {"OH": 476.4, "C=O": 986.5, "CH2": 0.0},
  "CH2": {"OH": 1318.0, "CH3": 0.0},
  "C=O": {"OH": -229.1, "CH3": 986.5},
  "C6H5": {"CH3": 1318.0, "OH": -229.1},
  "NH2": {"CH3": 1500.0, "C=O": 1200.0},
  "H2O": {"CH3": 3000.0, "OH": 0.0},
};

final List<Compound> compounds = [
  Compound(
    name: "Limonene",
    MW: 136.23,
    Psat: 2100,
    Thr: 0.0093,
    R: 2.2,
    Q: 2.5,
    groups: {"CH3": 2, "C=C": 1},
  ),
  Compound(
    name: "Geraniol",
    MW: 154.25,
    Psat: 26.7,
    Thr: 2.48e-5,
    R: 2.5,
    Q: 2.8,
    groups: {"CH3": 2, "CH2": 4, "OH": 1},
  ),
  Compound(
    name: "Vanillin",
    MW: 152.15,
    Psat: 0.23,
    Thr: 0.0002,
    R: 2.0,
    Q: 2.4,
    groups: {"CH3": 1, "C=O": 1, "OH": 1},
  ),
  Compound(
    name: "Galaxolide",
    MW: 258.4,
    Psat: 0.12,
    Thr: 0.0003,
    R: 3.2,
    Q: 3.6,
    groups: {"CH3": 3, "C6H5": 1},
  ),
  Compound(
    name: "Eugenol",
    MW: 164.2,
    Psat: 20.0,
    Thr: 0.00015,
    R: 2.7,
    Q: 3.0,
    groups: {"CH3": 1, "OH": 1, "C6H5": 1},
  ),
  Compound(
    name: "Cinnamaldehyde",
    MW: 132.16,
    Psat: 14.5,
    Thr: 0.0004,
    R: 2.3,
    Q: 2.6,
    groups: {"CH3": 1, "C=O": 1},
  ),
  Compound(
    name: "Methyl Anthranilate",
    MW: 151.16,
    Psat: 5.2,
    Thr: 0.0001,
    R: 2.1,
    Q: 2.3,
    groups: {"CH3": 1, "C=O": 1, "NH2": 1},
  ),
  Compound(
    name: "Ethanol",
    MW: 46.07,
    Psat: 7270,
    Thr: 5.53e-2,
    R: 1.2,
    Q: 1.6,
    groups: {"CH3": 1, "OH": 1},
  ),
  Compound(
    name: "Water",
    MW: 18.015,
    Psat: 3170,
    Thr: 1.0e-2,
    R: 1.4,
    Q: 1.9,
    groups: {"H2O": 1},
  ),
];

class SolverPage extends StatefulWidget {
  @override
  _SolverPageState createState() => _SolverPageState();
}

class _SolverPageState extends State<SolverPage> {
  bool isLoading = false;
  double progressValue = 0.0; // used for LinearProgressIndicator
  List<String> selectedCompounds = List.generate(
    7,
    (index) => compounds[index].name,
  );
  String selectedSolvent = "Ethanol";
  double solventFraction = 0.3;
  List<double> fractions = [];
  String result = "";
  List<FlSpot> ternaryPoints = [];

  @override
  void initState() {
    super.initState();
    generateRandomFractions();
  }

  void generateRandomFractions() {
    List<double> tempFractions;
    do {
      tempFractions = List.generate(7, (_) => getRandomFraction());
      double sum = tempFractions.reduce((a, b) => a + b) + solventFraction;
      if ((sum - 1).abs() < 0.0001) {
        fractions = tempFractions;
        break;
      }
    } while (true);
  }

  bool isOVUniform(List<double> ovList, {double tolerance = 0.05}) {
    double mean = ovList.reduce((a, b) => a + b) / ovList.length;
    return ovList.every((v) => (v - mean).abs() / mean <= tolerance);
  }

  void solveEquations() async {
    setState(() {
      isLoading = true;
      progressValue = 0.0;
    });

    await Future.delayed(Duration(milliseconds: 50));

    int maxTries = 300;
    int attempt = 0;

    List<Compound> bestSelected = [];
    List<double> bestFractions = [];
    List<double> bestGammas = [];
    List<double> bestOV = [];
    double lowestOVStdDev = double.infinity;

    while (attempt < maxTries) {
      await Future.delayed(Duration(milliseconds: 1));

      generateRandomFractions();

      List<Compound> selected =
          selectedCompounds
              .map((name) => compounds.firstWhere((c) => c.name == name))
              .toList();
      List<double> gammas =
          selected
              .map((comp) => calculateActivityCoefficient(comp, selected))
              .toList();
      List<double> OV = List.generate(
        selected.length,
        (i) => calculateOdorValue(selected[i], fractions[i], gammas[i]),
      );

      // Hitung standard deviation OV
      double mean = OV.reduce((a, b) => a + b) / OV.length;
      double stdDev = sqrt(
        OV.map((v) => pow(v - mean, 2)).reduce((a, b) => a + b) / OV.length,
      );

      if (stdDev < lowestOVStdDev) {
        lowestOVStdDev = stdDev;
        bestSelected = selected;
        bestFractions = List.from(fractions);
        bestGammas = gammas;
        bestOV = OV;
      }

      attempt++;
      setState(() {
        progressValue = attempt / maxTries;
      });
    }

    // Final result
    Compound air = compounds.firstWhere((c) => c.name == "Water");
    Compound etanol = compounds.firstWhere((c) => c.name == "Ethanol");

    Compound solvent;
    if (selectedSolvent == "Ethanol+Water") {
      solvent = Compound(
        name: "Ethanol+Water",
        MW: 0,
        Psat: 0,
        Thr: 1,
        R: 0,
        Q: 0,
        groups: {},
      );
    } else {
      solvent = compounds.firstWhere((c) => c.name == selectedSolvent);
    }

    double solventGamma = calculateActivityCoefficient(solvent, bestSelected);
    double solventOV = calculateOdorValue(
      solvent,
      solventFraction,
      solventGamma,
    );

    setState(() {
      result = "Fraksi Mol:\n";
      for (int i = 0; i < bestSelected.length; i++) {
        result +=
            "${bestSelected[i].name} = ${bestFractions[i].toStringAsPrecision(4)}\n";
      }
      result +=
          "Pelarut (${solvent.name}) = ${solventFraction.toStringAsPrecision(4)}\n\n";

      result += "\nKoefisien Aktivitas:\n";
      for (int i = 0; i < bestSelected.length; i++) {
        result +=
            "γ(${bestSelected[i].name}) = ${bestGammas[i].toStringAsPrecision(4)}\n";
      }

      result += "\nOdor Value:\n";
      for (int i = 0; i < bestSelected.length; i++) {
        result +=
            "OV(${bestSelected[i].name}) = ${bestOV[i].toStringAsPrecision(4)}\n";
      }

      result +=
          "\nKoefisien Aktivitas Pelarut (${solvent.name}): ${solventGamma.toStringAsPrecision(4)}";
      result +=
          "\nOdor Value Pelarut (${solvent.name}): ${solventOV.toStringAsPrecision(4)}";

      result += "\n\nKomparasi OV Pelarut:\n";
      if (selectedSolvent == "Water") {
        double ovAir = hitungOVPelarut(air, solventFraction, bestSelected);
        result += "OV Air (Water)         : ${ovAir.toStringAsPrecision(4)}\n";
      } else if (selectedSolvent == "Ethanol") {
        double ovEtanol = hitungOVPelarut(
          etanol,
          solventFraction,
          bestSelected,
        );
        result +=
            "OV Etanol              : ${ovEtanol.toStringAsPrecision(4)}\n";
      } else if (selectedSolvent == "Ethanol+Water") {
        double frEtanol = 0.6 * solventFraction;
        double frAir = 0.4 * solventFraction;
        double ovCampuran =
            hitungOVPelarut(etanol, frEtanol, bestSelected) +
            hitungOVPelarut(air, frAir, bestSelected);
        result +=
            "OV Campuran (60:40)    : ${ovCampuran.toStringAsPrecision(4)}\n";
      }

      fractions = bestFractions;
      ternaryPoints = generateTernaryPoints();
      isLoading = false;
    });
  }

  double hitungOVPelarut(
    Compound solvent,
    double fraction,
    List<Compound> selectedCompounds,
  ) {
    double gamma = calculateActivityCoefficient(solvent, selectedCompounds);
    return calculateOdorValue(solvent, fraction, gamma);
  }

  double getRandomFraction() {
    return 0.001 + Random().nextDouble() * (0.7 - 0.001);
  }

  double calculateActivityCoefficient(Compound comp, List<Compound> others) {
    double tauSum = 0.0;
    comp.groups.forEach((group1, count1) {
      others.forEach((other) {
        other.groups.forEach((group2, count2) {
          if (anm.containsKey(group1) && anm[group1]!.containsKey(group2)) {
            double tau = exp(-anm[group1]![group2]! / (8.314 * 298.15));
            tauSum += tau * count1 * count2;
          }
        });
      });
    });
    return exp(tauSum);
  }

  double calculateOdorValue(Compound comp, double x, double gamma) {
    double yi = gamma * x * comp.Psat / 101325;
    return (yi * comp.MW * 101325) / (8.314 * 298.15) / comp.Thr;
  }

  List<FlSpot> generateTernaryPoints() {
    List<FlSpot> points = [];
    for (int i = 0; i < fractions.length; i++) {
      double xTernary = fractions[i] * 100;
      double yTernary = (fractions[i] + solventFraction) * 100;
      points.add(FlSpot(xTernary, yTernary));
    }
    return points;
  }

  @override
  Widget build(BuildContext context) {
    return Scaffold(
      appBar: AppBar(title: Text("PQ2D Application")),
      body: Padding(
        padding: EdgeInsets.all(16.0),
        child: SingleChildScrollView(
          child: Column(
            children: [
              Text("Pilih 7 Senyawa:"),
              for (int i = 0; i < 7; i++)
                DropdownButton<String>(
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
                                c.name != selectedSolvent && c.name != "Water",
                          )
                          .map(
                            (c) => DropdownMenuItem(
                              value: c.name,
                              child: Text(c.name),
                            ),
                          )
                          .toList(),
                ),
              Text("\nPilih Pelarut:"),
              DropdownButton<String>(
                value: selectedSolvent,
                onChanged: (newValue) {
                  setState(() {
                    selectedSolvent = newValue!;
                  });
                },
                items: [
                  DropdownMenuItem(
                    value: "Ethanol",
                    child: Text("Ethanol (100%)"),
                  ),
                  DropdownMenuItem(
                    value: "Ethanol+Water",
                    child: Text("Ethanol + Water (60:40)"),
                  ),
                  DropdownMenuItem(value: "Water", child: Text("Air (100%)")),
                ],
              ),
              TextField(
                decoration: InputDecoration(labelText: "Fraksi Mol Pelarut"),
                keyboardType: TextInputType.number,
                onChanged: (value) {
                  double? parsed = double.tryParse(value);
                  if (parsed != null && parsed >= 0.001 && parsed <= 0.7) {
                    setState(() {
                      solventFraction = parsed;
                    });
                  }
                },
              ),
              SizedBox(height: 16.0),
              isLoading
                  ? Column(
                    children: [
                      SizedBox(height: 20),
                      LinearProgressIndicator(value: progressValue),
                      SizedBox(height: 12),
                      Text(
                        "Menghitung OV seragam... \${(progressValue * 100).toStringAsFixed(1)}%",
                      ),
                    ],
                  )
                  : Column(
                    children: [
                      ElevatedButton(
                        onPressed: solveEquations,
                        child: Text("Solve"),
                      ),
                      SizedBox(height: 16.0),
                      Text(result),
                    ],
                  ),
              Container(
                height: 300,
                child:
                    ternaryPoints.isEmpty
                        ? Center(child: Text("Generate Ternary Plot"))
                        : LineChart(
                          LineChartData(
                            gridData: FlGridData(show: true),
                            borderData: FlBorderData(show: true),
                            lineBarsData: [
                              LineChartBarData(
                                spots: ternaryPoints,
                                isCurved: true,
                                belowBarData: BarAreaData(show: true),
                              ),
                            ],
                          ),
                        ),
              ),
            ],
          ),
        ),
      ),
    );
  }
}
