// ============================================================
// ASKREATIF Perfumery Engine — Phase 3
// Advanced ROI / Odor Impact Model & Dynamic Note Classification
// ============================================================
// ROI_i ∝ (x_i × Psat_i) / (MW_i × Thr_i)
// Notes are no longer static labels — they are continuous,
// time-dependent perceptual states.
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';

// ── ROI Calculation ───────────────────────────────────────
/// Relative Odor Impact for compound i in a mixture.
/// Higher ROI = greater perceptual dominance.
double roiValue({required Compound c, required double moleFraction}) {
  if (c.Thr <= 0 || c.MW <= 0) return 0.0;
  final double roi = (moleFraction * c.Psat) / (c.MW * c.Thr);
  return (roi.isNaN || roi.isInfinite) ? 0.0 : roi;
}

/// ROI vector for all compounds in a mixture.
List<double> roiVector({
  required List<Compound> compounds,
  required List<double> moleFractions,
}) {
  return [
    for (int i = 0; i < compounds.length; i++)
      roiValue(c: compounds[i], moleFraction: moleFractions[i]),
  ];
}

/// Normalised ROI (sums to 1.0) — odor dominance weight.
List<double> normalisedRoi({
  required List<Compound> compounds,
  required List<double> moleFractions,
}) {
  final List<double> roi = roiVector(
    compounds: compounds,
    moleFractions: moleFractions,
  );
  final double total = roi.reduce((a, b) => a + b);
  if (total < 1e-30) return List.filled(compounds.length, 0.0);
  return roi.map((v) => v / total).toList();
}

// ── CORRECTED ROI — thermodynamically consistent ──────────
// Perbaikan: ROI menggunakan Pi = γ_i × x_i × Psat_i
// bukan x_i × Psat_i langsung (yang mengabaikan γ_i)
//
// ROI_corrected_i = Pi / (MW_i × Thr_i)
//                 = (γ_i × x_i × Psat_i) / (MW_i × Thr_i)
//
// Ini konsisten dengan OV karena Pi adalah driving force
// yang sama untuk VLE dan headspace release.

double roiValueCorrected({
  required Compound c,
  required double moleFraction,
  required double activityCoefficient, // γ_i dari UNIFAC
}) {
  if (c.Thr <= 0 || c.MW <= 0) return 0.0;
  final double Pi = activityCoefficient * moleFraction * c.Psat;
  final double roi = Pi / (c.MW * c.Thr);
  return (roi.isNaN || roi.isInfinite) ? 0.0 : roi;
}

List<double> roiVectorCorrected({
  required List<Compound> compounds,
  required List<double> moleFractions,
  required List<double> activityCoefficients,
}) {
  return [
    for (int i = 0; i < compounds.length; i++)
      roiValueCorrected(
        c: compounds[i],
        moleFraction: moleFractions[i],
        activityCoefficient: activityCoefficients[i],
      ),
  ];
}

List<double> normalisedRoiCorrected({
  required List<Compound> compounds,
  required List<double> moleFractions,
  required List<double> activityCoefficients,
}) {
  final List<double> roi = roiVectorCorrected(
    compounds: compounds,
    moleFractions: moleFractions,
    activityCoefficients: activityCoefficients,
  );
  final double total = roi.fold(0.0, (a, b) => a + b);
  if (total < 1e-30) return List.filled(compounds.length, 0.0);
  return roi.map((v) => v / total).toList();
}

// ── Dynamic Note Classification ────────────────────────────
/// Note state encodes both classical (Top/Heart/Base) and
/// perceptual (Dominant/Supporting/Trace) dimensions.
class NoteState {
  final String classicNote; // Top | Heart | Base
  final String perceptualRole; // Dominant | Supporting | Trace
  final double roiWeight; // 0–1 normalised dominance
  final double evaporatedFraction; // 0–1, how much has evaporated

  NoteState({
    required this.classicNote,
    required this.perceptualRole,
    required this.roiWeight,
    required this.evaporatedFraction,
  });

  /// A short readable label combining classic + perceptual.
  String get fullLabel => '$classicNote · $perceptualRole';

  /// CSS/UI colour hint based on note class.
  String get colorHex {
    switch (classicNote) {
      case 'Top':
        return '#4FC3F7'; // cyan
      case 'Heart':
        return '#A5D6A7'; // green
      case 'Base':
        return '#CE93D8'; // purple
      default:
        return '#78909C'; // grey
    }
  }
}

/// Classify a single compound's dynamic note state.
NoteState classifyDynamicNote({
  required Compound c,
  required double moleFraction,
  required double normRoi, // 0–1
  required double t, // elapsed time (hours)
}) {
  // ── Classic note by Psat ──────────────────────────────
  String classic;
  if (c.Psat >= 5.0) {
    classic = 'Top';
  } else if (c.Psat >= 0.3) {
    classic = 'Heart';
  } else {
    classic = 'Base';
  }

  // ── Evaporation fraction ──────────────────────────────
  final double k = _kConst(c);
  final double evap = 1.0 - exp(-k * t); // fraction evaporated

  // ── Perceptual role by normalised ROI ─────────────────
  String role;
  if (normRoi >= 0.25) {
    role = 'Dominant';
  } else if (normRoi >= 0.08) {
    role = 'Supporting';
  } else if (normRoi >= 0.01) {
    role = 'Trace';
  } else {
    role = 'Subliminal';
  }

  return NoteState(
    classicNote: classic,
    perceptualRole: role,
    roiWeight: normRoi,
    evaporatedFraction: evap,
  );
}

/// Classify all compounds in a mixture at time t.
List<NoteState> classifyMixture({
  required List<Compound> compounds,
  required List<double> moleFractions,
  required double t,
}) {
  final List<double> nRoi = normalisedRoi(
    compounds: compounds,
    moleFractions: moleFractions,
  );
  return [
    for (int i = 0; i < compounds.length; i++)
      classifyDynamicNote(
        c: compounds[i],
        moleFraction: moleFractions[i],
        normRoi: nRoi[i],
        t: t,
      ),
  ];
}

double _kConst(Compound c) {
  const double alpha = 0.0012;
  if (c.Psat <= 0) return 1e-6;
  return alpha * c.Psat / sqrt(c.MW);
}

// ── Accord Suggestion based on ROI profile ─────────────────
/// Returns the top 3 accord archetypes that best match the
/// current ROI-weighted compound profile.
List<String> suggestAccords(List<Compound> compounds, List<double> normRoi) {
  // Simple group-based heuristic. Each accord maps to UNIFAC groups.
  const Map<String, List<String>> accordGroups = {
    'Woody': ['CH3', 'CH2', 'CH', 'C', 'C=C'],
    'Floral': ['OH', 'ACH', 'ACOH', 'CH=C'],
    'Amber': ['CH3CO', 'CH2CO', 'COO', 'CH2COO'],
    'Musk': ['CH2O', 'CH3COO', 'CH2COO'],
    'Fresh': ['CH=CH', 'CH2=CH', 'CHO', 'CH=C'],
    'Gourmand': ['CH3CO', 'CH2CO', 'CH3COO'],
    'Animalic': ['ACOH', 'OH'],
    'Aquatic': ['CH2', 'CH3', 'CH=CH'],
  };

  Map<String, double> accordScore = {};
  for (final accord in accordGroups.keys) {
    double score = 0.0;
    for (int i = 0; i < compounds.length; i++) {
      for (final group in accordGroups[accord]!) {
        if (compounds[i].groups.containsKey(group)) {
          score += normRoi[i] * compounds[i].groups[group]!;
        }
      }
    }
    accordScore[accord] = score;
  }

  final sorted =
      accordScore.entries.toList()..sort((a, b) => b.value.compareTo(a.value));
  return sorted.take(3).map((e) => e.key).toList();
}
