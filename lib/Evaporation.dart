// ============================================================
// ASKREATIF Perfumery Engine — Phase 1 + 2
// Dynamic Evaporation Model & Headspace Concentration
// ============================================================
// dN_i/dt = -k_i * N_i  →  N_i(t) = N_i(0) * exp(-k_i * t)
// C_i(t) = (y_i(t) * P_total * MW_i) / (R * T)
// ============================================================
import 'dart:math';
import 'package:askreatif_app/Compound.dart';

const double _R = 8.314462618; // J/(mol·K)
const double _T = 298.15; // K (25°C)

// Evaporation rate constant k_i
// Derived from vapour pressure + molecular weight
//   k_i = α * Psat_i / (MW_i^0.5)
// α is a scaling factor that maps to real-world evaporation
// timescales (hours). Tune α to match experimental data.
const double _alpha = 0.0012; // h⁻¹ normaliser

double evaporationConstant(Compound c) {
  if (c.Psat <= 0) return 1e-6;
  return _alpha * c.Psat / sqrt(c.MW);
}

// ── Evaporation time points (hours) ───────────────────────
// Fine-grained for first 2 h, then broader intervals up to 24 h
List<double> defaultTimePoints() {
  final List<double> t = [];
  for (double h = 0; h <= 2.0; h += 0.1) t.add(h);
  for (double h = 2.5; h <= 6.0; h += 0.5) t.add(h);
  for (double h = 7.0; h <= 24.0; h += 1.0) t.add(h);
  return t;
}

// ── Single time-step concentration decay ──────────────────
/// Returns mole-fraction vector at time t (hours).
/// x_i(t) = x_i(0) * exp(-k_i * t)   then renormalised so Σ = 1
List<double> moleFractionsAtTime({
  required List<Compound> compounds,
  required List<double> x0,
  required double t,
}) {
  final int n = compounds.length;
  List<double> xt = List.filled(n, 0.0);
  for (int i = 0; i < n; i++) {
    final double k = evaporationConstant(compounds[i]);
    xt[i] = x0[i] * exp(-k * t);
  }
  double sum = xt.reduce((a, b) => a + b);
  if (sum < 1e-30) return List.filled(n, 1.0 / n);
  for (int i = 0; i < n; i++) xt[i] /= sum;
  return xt;
}

// ── Headspace concentration (mg/m³) ───────────────────────
/// C_i(t) = (y_i(t) × P_total(t) × MW_i) / (R × T)
/// y_i(t) = γ_i × x_i(t) × Psat_i  /  P_total(t)
/// Returns list of C_i in mg/m³
List<double> headspaceConcentration({
  required List<Compound> compounds,
  required List<double> xt, // mole fractions at time t
  required List<double> gammas, // UNIFAC activity coefficients
}) {
  final int n = compounds.length;

  // Partial pressures (kPa)
  List<double> pp = List.filled(n, 0.0);
  for (int i = 0; i < n; i++) {
    double p = gammas[i] * xt[i] * compounds[i].Psat;
    pp[i] = (p.isNaN || p.isInfinite || p < 0) ? 1e-12 : p;
  }
  double Ptotal = pp.reduce((a, b) => a + b);
  if (Ptotal < 1e-12) Ptotal = 1e-12;

  // C_i [mol/m³] = (y_i × P_total [Pa]) / (R × T)
  // P in kPa → ×1000 for Pa; MW in g/mol → ×1000 for mg/mol → mg/m³
  List<double> C = List.filled(n, 0.0);
  for (int i = 0; i < n; i++) {
    double yi = pp[i] / Ptotal;
    double Ci = (yi * Ptotal * 1000.0 * compounds[i].MW) / (_R * _T);
    C[i] = (Ci.isNaN || Ci.isInfinite) ? 0.0 : Ci; // mg/m³
  }
  return C;
}

// ── IMPROVED: Activity-coefficient-aware evaporation ─────
// k_i sekarang dimodulasi oleh γ_i(t) karena saat komposisi
// berubah, activity coefficient berubah → driving force berubah
//
// Effective evaporation rate:
//   k_eff_i(t) = α × γ_i(t) × Psat_i / √MW_i
//
// Diffusion resistance correction (stagnant film model):
//   k_diff_i = k_eff_i / (1 + δ × MW_i^0.33)
// di mana δ adalah diffusion damping factor

double evaporationConstantCorrected({
  required Compound c,
  required double activityCoefficient,
  double alpha = 0.0012,
  double diffusionDamping = 0.002,
}) {
  if (c.Psat <= 0) return 1e-6;
  // Activity-coefficient correction
  final double kEff = alpha * activityCoefficient * c.Psat / sqrt(c.MW);
  // Diffusion resistance: heavier molecules diffuse slower
  final double diffFactor = 1.0 / (1.0 + diffusionDamping * pow(c.MW, 0.33));
  final double k = kEff * diffFactor;
  return (k.isFinite && k > 0) ? k : alpha * c.Psat / sqrt(c.MW);
}

// Solvent retention effect:
// Ethanol memperlambat evaporasi senyawa polar melalui
// kompetisi untuk vapor phase
// Efek ini dimodelkan sebagai faktor multiplicative pada k_i
double solventRetentionFactor({
  required Compound c,
  required double ethanolFraction,
  required double waterFraction,
}) {
  // Polar compounds (OH, ACOH, COO) lebih tahan dalam
  // campuran ethanol-water
  final double polarity = _groupPolarityScore(c);
  // Retention lebih kuat di kondisi basah (water)
  final double retentionStrength =
      ethanolFraction * 0.3 * polarity + waterFraction * 0.5 * polarity;
  return (1.0 - retentionStrength.clamp(0.0, 0.7));
}

double _groupPolarityScore(Compound c) {
  const Map<String, double> polarity = {
    'OH': 1.0,
    'ACOH': 0.9,
    'COOH': 1.0,
    'CH3CO': 0.4,
    'CHO': 0.5,
    'COO': 0.35,
    'CH3COO': 0.3,
    'CH2O': 0.4,
    'H2O': 1.0,
  };
  double score = 0.0, total = 0.0;
  c.groups.forEach((g, n) {
    score += (polarity[g] ?? 0.05) * n;
    total += n;
  });
  return total > 0 ? (score / total).clamp(0.0, 1.0) : 0.05;
}

// ── Full evaporation trajectory ────────────────────────────
/// Returns one EvaporationPoint per time-step.
/// gammaFn: a callback that computes UNIFAC gammas for a given
///          mole-fraction vector (reuse your existing UNIFAC logic).
class EvaporationPoint {
  final double time; // hours
  final List<double> moleFractions; // normalised, per compound
  final List<double> headspaceConc; // mg/m³, per compound
  final List<double> odorValues; // dimensionless OV
  final List<String> noteLabels; // dynamic note classification

  EvaporationPoint({
    required this.time,
    required this.moleFractions,
    required this.headspaceConc,
    required this.odorValues,
    required this.noteLabels,
  });
}

List<EvaporationPoint> computeEvaporationTrajectory({
  required List<Compound> compounds,
  required List<double> x0,
  required List<double> Function(List<double>) gammaFn,
  List<double>? timePoints,
}) {
  final List<double> ts = timePoints ?? defaultTimePoints();
  final int n = compounds.length;

  return ts.map((t) {
    final List<double> xt = moleFractionsAtTime(
      compounds: compounds,
      x0: x0,
      t: t,
    );
    final List<double> gammas = gammaFn(xt);
    final List<double> C = headspaceConcentration(
      compounds: compounds,
      xt: xt,
      gammas: gammas,
    );

    // Odor value from headspace: OV_i = C_i / (Thr_i × 1e6)
    // Thr in μg/L = mg/m³  → direct ratio
    List<double> ov = List.filled(n, 0.0);
    for (int i = 0; i < n; i++) {
      double thr = compounds[i].Thr;
      if (thr <= 0 || thr.isNaN || thr.isInfinite) thr = 1.0;
      ov[i] = C[i] / thr;
      if (ov[i].isNaN || ov[i].isInfinite) ov[i] = 0.0;
    }

    // Dynamic note classification via ROI (Phase 3)
    final List<String> notes =
        compounds.map((c) {
          return _dynamicNoteLabel(c, t);
        }).toList();

    return EvaporationPoint(
      time: t,
      moleFractions: xt,
      headspaceConc: C,
      odorValues: ov,
      noteLabels: notes,
    );
  }).toList();
}

// ── Dynamic note label (Phase 3 ROI preview) ──────────────
// Full ROI model lives in roi_model.dart; this is a lightweight
// time-aware classification used inside the trajectory.
String _dynamicNoteLabel(Compound c, double t) {
  final double k = evaporationConstant(c);
  // Residual fraction at time t
  final double residual = exp(-k * t);

  if (residual < 0.05) return 'Dry-down / Dry';
  if (c.Psat >= 5.0) return 'Top';
  if (c.Psat >= 0.5) return t < 1.0 ? 'Top-Mid' : 'Heart';
  if (c.Psat >= 0.05) return t < 3.0 ? 'Heart' : 'Heart-Base';
  return 'Base';
}
