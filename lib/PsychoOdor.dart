// ============================================================
// ASKREATIF Perfumery Engine — Phase 4 + 5
// Psychophysical Odor Model (Stevens' Power Law) +
// Odor Interaction Matrix (Synergy / Masking / Inhibition)
// ============================================================
// I = k × C^n        (Stevens, 1957)
// I_mix ≠ Σ I_i      (Cain & Drexler, 1974)
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';

// ════════════════════════════════════════════════════════════
// PHASE 4 — STEVENS' POWER LAW
// ════════════════════════════════════════════════════════════

/// Empirical Stevens' exponent for smell.
/// Literature range: 0.4 – 0.7 depending on odorant class.
/// Default: 0.55 (geometric mean of reported values).
const double _defaultStevensN = 0.55;
const double _defaultStevensK = 1.0;  // relative scale

/// Perceptual intensity for a single compound.
/// C is the supra-threshold concentration ratio: C_i / Threshold_i
double stevensIntensity({
  required double concentration,  // C_i (mg/m³ headspace)
  required double threshold,      // Thr_i (mg/m³)
  double k = _defaultStevensK,
  double n = _defaultStevensN,
}) {
  if (threshold <= 0 || concentration <= 0) return 0.0;
  final double C = concentration / threshold;
  if (C < 1.0) return 0.0; // sub-threshold
  final double I = k * pow(C, n).toDouble();
  return (I.isNaN || I.isInfinite) ? 0.0 : I;
}

/// Perceptual intensity vector for the full mixture.
List<double> intensityVector({
  required List<Compound> compounds,
  required List<double> headspaceConc, // mg/m³
  double k = _defaultStevensK,
  double n = _defaultStevensN,
}) {
  return [
    for (int i = 0; i < compounds.length; i++)
      stevensIntensity(
        concentration: headspaceConc[i],
        threshold: compounds[i].Thr,
        k: k,
        n: n,
      )
  ];
}

// ════════════════════════════════════════════════════════════
// PHASE 5 — ODOR INTERACTION MATRIX
// ════════════════════════════════════════════════════════════

/// Interaction type between two odorants.
enum OdorInteraction { synergy, neutral, masking, inhibition }

/// Pairwise interaction coefficient matrix.
/// Values > 1  → synergy   (mutual enhancement)
/// Values = 1  → neutral   (additive)
/// 0 < val < 1 → masking   (partial suppression)
/// Values ≤ 0  → inhibition (near-complete suppression)
///
/// Currently populated with evidence-based heuristics.
/// Replace individual entries with measured data as available.
class OdorInteractionMatrix {
  final Map<String, Map<String, double>> _matrix = {};

  OdorInteractionMatrix() {
    _populateDefaults();
  }

  /// Set or override a pairwise coefficient.
  void set(String a, String b, double coefficient) {
    _matrix[a] ??= {};
    _matrix[b] ??= {};
    _matrix[a]![b] = coefficient;
    _matrix[b]![a] = coefficient; // symmetric
  }

  /// Get coefficient; returns 1.0 (neutral) if undefined.
  double get(String a, String b) {
    return _matrix[a]?[b] ?? 1.0;
  }

  OdorInteraction classify(double coeff) {
    if (coeff > 1.05) return OdorInteraction.synergy;
    if (coeff < 0.0)  return OdorInteraction.inhibition;
    if (coeff < 0.85) return OdorInteraction.masking;
    return OdorInteraction.neutral;
  }

  // ── Evidence-based default interactions ────────────────
  void _populateDefaults() {
    // Floral × Musk → synergy (classic accord)
    set('Geraniol',     'Ambroxan',       1.25);
    set('Linalool',     'Ambroxan',       1.30);
    set('Nerol',        'Galaxolide',     1.20);
    set('Citronellol',  'Ambroxan',       1.15);

    // Citrus top × fresh heart → slight synergy
    set('Limonene',     'Linalool',       1.10);
    set('Limonene',     'Geraniol',       1.08);

    // Woody × amber → synergy
    set('Cashmeran',    'Romandolide',    1.18);
    set('Iso E Super',  'Ambroxan',       1.35);
    set('Javanol',      'Ambroxan',       1.22);

    // Aldehyde × floral → slight masking (aldehydes can overpower)
    set('CIS-6-NONEAL', 'Geraniol',       0.85);
    set('Benzaldehyde', 'Vanillin',       0.90);

    // Phenolic × musk → inhibition
    set('Eugenol',      'Galaxolide',     0.60);
    set('Eugenol',      'Ambroxan',       0.70);

    // Coumarin × musk → slight synergy
    set('Coumarin',     'Ambroxan',       1.12);
    set('Coumarin',     'Galaxolide',     1.10);

    // Vanilla × amber → synergy
    set('Vanillin',     'Ambroxan',       1.28);
    set('Vanillin',     'Romandolide',    1.15);

    // Fresh green × citrus → neutral/mild synergy
    set('Cis-3-Hexenyl Acetate', 'Limonene', 1.05);

    // Heavy orientals can mask lighter florals
    set('Patchouli Alcohol', 'Nerol',     0.75);
    set('Patchouli Alcohol', 'Linalool',  0.78);
  }
}

// Shared singleton instance
final OdorInteractionMatrix globalInteractionMatrix = OdorInteractionMatrix();

// ── Mixture perceived intensity (non-additive model) ───────
/// Computes the blended perceived intensity for each compound
/// accounting for pairwise interaction coefficients.
///
/// Model: I_mix_i = I_i × Π_j≠i  f(interaction_ij, I_j)
/// where f = 1 + (coeff - 1) × (I_j / ΣI)   (weighted blend)
List<double> blendedIntensity({
  required List<Compound> compounds,
  required List<double> rawIntensities, // from stevensIntensity
  OdorInteractionMatrix? matrix,
}) {
  final OdorInteractionMatrix m = matrix ?? globalInteractionMatrix;
  final int n = compounds.length;
  if (n == 0) return [];

  double totalI = rawIntensities.fold(0.0, (a, b) => a + b);
  if (totalI < 1e-30) return List.from(rawIntensities);

  List<double> blended = List.from(rawIntensities);

  for (int i = 0; i < n; i++) {
    double modFactor = 1.0;
    for (int j = 0; j < n; j++) {
      if (i == j) continue;
      final double coeff = m.get(compounds[i].name, compounds[j].name);
      // Weighted influence: compounds with higher intensity
      // exert stronger interaction effects on compound i.
      final double wj = rawIntensities[j] / totalI;
      modFactor *= (1.0 + (coeff - 1.0) * wj);
    }
    blended[i] = rawIntensities[i] * modFactor;
    if (blended[i].isNaN || blended[i].isInfinite || blended[i] < 0) {
      blended[i] = rawIntensities[i];
    }
  }
  return blended;
}

/// Total perceived intensity of the blend (sum of blended).
double totalBlendedIntensity({
  required List<Compound> compounds,
  required List<double> rawIntensities,
  OdorInteractionMatrix? matrix,
}) {
  final List<double> b = blendedIntensity(
    compounds: compounds,
    rawIntensities: rawIntensities,
    matrix: matrix,
  );
  return b.fold(0.0, (a, v) => a + v);
}