// ============================================================
// ASKREATIF Perfumery Engine — Phase 6 (Improved)
// Accord Intelligence System
// ============================================================
// Changes from v1:
//   1. ROI menggunakan γ_i-corrected (thermodynamically consistent)
//   2. Group scoring menggunakan IDF-style weighting — group yang
//      ada di semua senyawa (CH3, CH2) diberi bobot lebih rendah
//   3. Accord conflict detection
//   4. Accord harmony scoring
//   5. Temporal accord evolution
//   6. Accord hierarchy modeling (dominant → supporting → trace)
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/Roi.dart';

// ── Accord Archetype Definitions (unchanged) ──────────────
class AccordProfile {
  final String name;
  final String icon;
  final String description;
  final List<String> keyGroups;
  final List<String> keyCompounds;

  // NEW: accord family untuk conflict/harmony detection
  final String family; // 'fresh', 'floral', 'woody', 'oriental', 'aromatic'

  // NEW: typical note position
  final String noteAffinity; // 'top', 'heart', 'base', 'heart-base'

  const AccordProfile({
    required this.name,
    required this.icon,
    required this.description,
    required this.keyGroups,
    required this.keyCompounds,
    required this.family,
    required this.noteAffinity,
  });
}

const List<AccordProfile> accordArchetypes = [
  AccordProfile(
    name: 'Floral',
    icon: '🌸',
    description: 'Delicate, petal-like, rosy, jasmine, ylang character.',
    keyGroups: ['OH', 'ACH', 'ACOH', 'CH=C', 'ACCH2'],
    keyCompounds: ['Geraniol', 'Nerol', 'Citronellol', 'Linalool'],
    family: 'floral',
    noteAffinity: 'heart',
  ),
  AccordProfile(
    name: 'Woody',
    icon: '🌲',
    description: 'Dry, cedar, sandalwood, vetiver, earthy warmth.',
    keyGroups: ['CH3', 'CH2', 'CH', 'C', 'C=C'],
    keyCompounds: ['Javanol', 'Sandalmysore Core', 'Cashmeran', 'Iso E Super'],
    family: 'woody',
    noteAffinity: 'base',
  ),
  AccordProfile(
    name: 'Amber / Oriental',
    icon: '🍯',
    description: 'Warm, resinous, sweet-balsamic, animalic depth.',
    keyGroups: ['COO', 'CH2COO', 'CH3COO', 'CH3CO', 'CH2CO'],
    keyCompounds: ['Ambroxan', 'Vanillin', 'Coumarin', 'Romandolide'],
    family: 'oriental',
    noteAffinity: 'base',
  ),
  AccordProfile(
    name: 'Musk',
    icon: '🌫️',
    description: 'Soft, skin-like, powdery, clean, diffusive.',
    keyGroups: ['CH2O', 'CH3COO', 'CH2COO', 'COO'],
    keyCompounds: ['Galaxolide', 'Ambroxan', 'Romandolide'],
    family: 'musk',
    noteAffinity: 'base',
  ),
  AccordProfile(
    name: 'Fresh / Citrus',
    icon: '🍃',
    description: 'Bright, airy, green, citrus, ozonic, uplifting.',
    keyGroups: ['CH=CH', 'CH2=CH', 'CH2=C', 'CH=C'],
    keyCompounds: [
      'Limonene',
      'Linalool',
      'Cis-3-Hexenyl Acetate',
      'CIS-6-NONEAL',
    ],
    family: 'fresh',
    noteAffinity: 'top',
  ),
  AccordProfile(
    name: 'Gourmand',
    icon: '🍮',
    description: 'Sweet, edible, caramel, vanilla, lactonic.',
    keyGroups: ['CH3CO', 'CH2CO', 'CH3COO'],
    keyCompounds: ['Vanillin', 'Isoamyl Acetate', 'Benzyl Acetate', 'Coumarin'],
    family: 'oriental',
    noteAffinity: 'heart-base',
  ),
  AccordProfile(
    name: 'Aromatic / Spicy',
    icon: '🌿',
    description: 'Herbal, medicinal, warm spice, clove, cinnamon.',
    keyGroups: ['ACH', 'AC', 'ACOH', 'CHO'],
    keyCompounds: ['Eugenol', 'Cinnamaldehyde', 'Vanillin', 'Benzaldehyde'],
    family: 'aromatic',
    noteAffinity: 'heart',
  ),
  AccordProfile(
    name: 'Fougère',
    icon: '☘️',
    description: 'Mossy, lavender-like, hay, coumarinic green.',
    keyGroups: ['COO', 'CH=CH', 'ACH', 'OH'],
    keyCompounds: ['Coumarin', 'Geraniol', 'Evernyl', 'Linalool'],
    family: 'aromatic',
    noteAffinity: 'heart',
  ),
];

// ── IDF-style Group Weights ───────────────────────────────
// Masalah v1: CH3 dan CH2 ada di hampir semua senyawa
// sehingga semua campuran cenderung terklasifikasi Woody.
//
// Fix: Bobot grup = 1 / (1 + ln(df_group))
// di mana df_group = jumlah senyawa dalam database yang mengandung grup ini.
// Grup yang lebih spesifik (ACOH, COO) → bobot lebih tinggi.
// Grup yang sangat umum (CH3, CH2) → bobot lebih rendah.

final Map<String, double> _groupSpecificity = _computeGroupSpecificity();

Map<String, double> _computeGroupSpecificity() {
  // Hitung document frequency setiap grup dari database
  final Map<String, int> df = {};
  for (final c in compounds) {
    for (final group in c.groups.keys) {
      df[group] = (df[group] ?? 0) + 1;
    }
  }
  final int totalDocs = compounds.length;
  final Map<String, double> weights = {};
  df.forEach((group, freq) {
    // IDF = log(N / df) — grup spesifik → IDF tinggi
    weights[group] = log(totalDocs / freq.toDouble()).clamp(0.1, 3.0);
  });
  return weights;
}

double _groupWeight(String group) {
  return _groupSpecificity[group] ?? 1.0;
}

// ── Accord Scoring (Improved) ─────────────────────────────
class AccordScore {
  final AccordProfile accord;
  final double score; // 0–1 normalised
  final double roiContrib; // raw weighted score
  final AccordRole role; // dominant / supporting / trace

  const AccordScore({
    required this.accord,
    required this.score,
    required this.roiContrib,
    required this.role,
  });
}

enum AccordRole { dominant, supporting, trace, absent }

/// Scores each accord using:
/// 1. γ_i-corrected ROI weights (thermodynamically consistent)
/// 2. IDF-weighted group specificity (no more CH3 bias)
/// 3. Key compound bonus
List<AccordScore> scoreAccords({
  required List<Compound> compounds,
  required List<double> moleFractions,
  List<double>? activityCoefficients, // NEW: optional γ_i
}) {
  // Gunakan ROI corrected jika γ_i tersedia, fallback ke lama
  final List<double> nRoi =
      activityCoefficients != null
          ? normalisedRoiCorrected(
            compounds: compounds,
            moleFractions: moleFractions,
            activityCoefficients: activityCoefficients,
          )
          : normalisedRoi(compounds: compounds, moleFractions: moleFractions);

  List<AccordScore> scores = [];

  for (final accord in accordArchetypes) {
    double score = 0.0;

    for (int i = 0; i < compounds.length; i++) {
      // IDF-weighted group scoring
      double groupScore = 0.0;
      double totalGroupWeight = 0.0;

      for (final group in accord.keyGroups) {
        final double idf = _groupWeight(group);
        totalGroupWeight += idf;
        if (compounds[i].groups.containsKey(group)) {
          final int count = compounds[i].groups[group]!;
          // Score = IDF × tanh(count/3) — tanh prevents dominance
          // dari senyawa dengan banyak grup yang sama
          groupScore += idf * tan(count / 3.0);
        }
      }

      // Normalize group score to [0,1]
      final double normalizedGroup =
          totalGroupWeight > 0
              ? (groupScore / totalGroupWeight).clamp(0.0, 1.0)
              : 0.0;

      // Key compound bonus (spesifik, tidak bisa digeneralisasi)
      final double compoundBonus =
          accord.keyCompounds.contains(compounds[i].name) ? 0.25 : 0.0;

      // ROI-weighted contribution
      score += nRoi[i] * (normalizedGroup * 0.75 + compoundBonus);
    }

    scores.add(
      AccordScore(
        accord: accord,
        score: score.clamp(0.0, 2.0),
        roiContrib: score,
        role: AccordRole.absent, // assigned after normalisation
      ),
    );
  }

  // Normalise relative to highest scorer
  final double maxScore = scores
      .map((s) => s.roiContrib)
      .fold(0.0, (a, b) => a > b ? a : b);

  if (maxScore > 0) {
    scores =
        scores.map((s) {
          final double normalized = s.roiContrib / maxScore;
          // Assign perceptual role based on relative score
          final AccordRole role = _assignRole(normalized);
          return AccordScore(
            accord: s.accord,
            score: normalized,
            roiContrib: s.roiContrib,
            role: role,
          );
        }).toList();
  }

  scores.sort((a, b) => b.score.compareTo(a.score));
  return scores;
}

AccordRole _assignRole(double normalizedScore) {
  if (normalizedScore >= 0.65) return AccordRole.dominant;
  if (normalizedScore >= 0.30) return AccordRole.supporting;
  if (normalizedScore >= 0.10) return AccordRole.trace;
  return AccordRole.absent;
}

// ── Accord Harmony Scoring ────────────────────────────────
// Harmony matrix: seberapa baik dua accord keluarga berpadu.
// Berdasarkan classical perfumery knowledge (fragrance wheel).
// Nilai 1.0 = perfect harmony, 0.0 = clash.

const Map<String, Map<String, double>> _harmonyMatrix = {
  'floral': {
    'floral': 1.00,
    'fresh': 0.85,
    'woody': 0.80,
    'oriental': 0.75,
    'musk': 0.90,
    'aromatic': 0.65,
  },
  'fresh': {
    'floral': 0.85,
    'fresh': 1.00,
    'woody': 0.60,
    'oriental': 0.40,
    'musk': 0.70,
    'aromatic': 0.80,
  },
  'woody': {
    'floral': 0.80,
    'fresh': 0.60,
    'woody': 1.00,
    'oriental': 0.90,
    'musk': 0.85,
    'aromatic': 0.75,
  },
  'oriental': {
    'floral': 0.75,
    'fresh': 0.40,
    'woody': 0.90,
    'oriental': 1.00,
    'musk': 0.85,
    'aromatic': 0.70,
  },
  'musk': {
    'floral': 0.90,
    'fresh': 0.70,
    'woody': 0.85,
    'oriental': 0.85,
    'musk': 1.00,
    'aromatic': 0.65,
  },
  'aromatic': {
    'floral': 0.65,
    'fresh': 0.80,
    'woody': 0.75,
    'oriental': 0.70,
    'musk': 0.65,
    'aromatic': 1.00,
  },
};

/// Harmony score antara dua accord (0 = clash, 1 = perfect).
double accordHarmony(AccordProfile a, AccordProfile b) {
  return _harmonyMatrix[a.family]?[b.family] ?? 0.5;
}

// ── Accord Intelligence Result ────────────────────────────
class AccordIntelligenceResult {
  final List<AccordScore> scores;
  final double harmonyScore; // 0–1 overall blend harmony
  final List<String> conflicts; // accord pairs that clash
  final List<String> synergies; // accord pairs that enhance
  final String dominantProfile; // primary accord description
  final String blendCharacter; // perfumer's summary
  final Map<String, double> balance; // accord → share %

  const AccordIntelligenceResult({
    required this.scores,
    required this.harmonyScore,
    required this.conflicts,
    required this.synergies,
    required this.dominantProfile,
    required this.blendCharacter,
    required this.balance,
  });
}

/// Full accord intelligence analysis.
AccordIntelligenceResult analyzeAccords({
  required List<Compound> compounds,
  required List<double> moleFractions,
  List<double>? activityCoefficients,
}) {
  final List<AccordScore> scores = scoreAccords(
    compounds: compounds,
    moleFractions: moleFractions,
    activityCoefficients: activityCoefficients,
  );

  // Only consider non-absent accords for analysis
  final List<AccordScore> active =
      scores.where((s) => s.role != AccordRole.absent).toList();

  // ── Harmony Score ─────────────────────────────────────
  double harmonySum = 0.0;
  int harmonyCount = 0;
  final List<String> conflicts = [];
  final List<String> synergies = [];

  for (int i = 0; i < active.length; i++) {
    for (int j = i + 1; j < active.length; j++) {
      final double h = accordHarmony(active[i].accord, active[j].accord);
      // Weight by product of scores (dominant accords matter more)
      harmonySum += h * active[i].score * active[j].score;
      harmonyCount++;

      final String pair = '${active[i].accord.name} × ${active[j].accord.name}';
      if (h < 0.50) conflicts.add(pair);
      if (h >= 0.85) synergies.add(pair);
    }
  }

  final double harmonyScore =
      harmonyCount > 0 ? (harmonySum / harmonyCount).clamp(0.0, 1.0) : 1.0;

  // ── Accord Balance Map ────────────────────────────────
  final double totalScore = active.fold(0.0, (a, s) => a + s.score);
  final Map<String, double> balance = {
    for (final s in active)
      s.accord.name: totalScore > 0 ? s.score / totalScore : 0.0,
  };

  // ── Dominant Profile ──────────────────────────────────
  final String dominantProfile =
      scores.isNotEmpty
          ? '${scores.first.accord.icon} ${scores.first.accord.name}'
          : 'Undefined';

  // ── Blend Character (perfumer summary) ───────────────
  final String blendCharacter = _generateBlendCharacter(scores, harmonyScore);

  return AccordIntelligenceResult(
    scores: scores,
    harmonyScore: harmonyScore,
    conflicts: conflicts,
    synergies: synergies,
    dominantProfile: dominantProfile,
    blendCharacter: blendCharacter,
    balance: balance,
  );
}

String _generateBlendCharacter(List<AccordScore> scores, double harmony) {
  final List<AccordScore> dominant =
      scores.where((s) => s.role == AccordRole.dominant).toList();
  final List<AccordScore> supporting =
      scores.where((s) => s.role == AccordRole.supporting).toList();

  if (dominant.isEmpty) return 'Diffuse, undefined character.';

  final String domName = dominant.first.accord.name;
  final String harmonyDesc =
      harmony >= 0.80
          ? 'well-harmonised'
          : harmony >= 0.60
          ? 'moderately balanced'
          : 'contrasting';

  String character = 'Primarily $domName';
  if (supporting.isNotEmpty) {
    final String suppNames = supporting.map((s) => s.accord.name).join(' and ');
    character += ', supported by $suppNames';
  }
  character += '. A $harmonyDesc blend.';
  return character;
}

// ── Temporal Accord Evolution ─────────────────────────────
/// Accord profile at time t, using evaporation-adjusted
/// mole fractions. Models how the scent character evolves.
class TemporalAccordSnapshot {
  final double time;
  final AccordIntelligenceResult accordResult;
  final List<double> moleFractionsAtT;

  const TemporalAccordSnapshot({
    required this.time,
    required this.accordResult,
    required this.moleFractionsAtT,
  });
}

List<TemporalAccordSnapshot> computeAccordEvolution({
  required List<Compound> compounds,
  required List<double> x0,
  required List<double> Function(List<double>) gammaFn,
  List<double> timePoints = const [0, 0.5, 1, 2, 4, 8, 24],
}) {
  return timePoints.map((t) {
    // Evaporated mole fractions
    final List<double> xt = _moleFractionsAtTime(
      compounds: compounds,
      x0: x0,
      t: t,
    );
    final List<double> gammas = gammaFn(xt);

    final AccordIntelligenceResult result = analyzeAccords(
      compounds: compounds,
      moleFractions: xt,
      activityCoefficients: gammas,
    );

    return TemporalAccordSnapshot(
      time: t,
      accordResult: result,
      moleFractionsAtT: xt,
    );
  }).toList();
}

// Internal helper (avoids circular import with Evaporation.dart)
List<double> _moleFractionsAtTime({
  required List<Compound> compounds,
  required List<double> x0,
  required double t,
}) {
  const double alpha = 0.0012;
  final int n = compounds.length;
  List<double> xt = List.filled(n, 0.0);
  for (int i = 0; i < n; i++) {
    final double k =
        compounds[i].Psat > 0
            ? alpha * compounds[i].Psat / sqrt(compounds[i].MW)
            : 1e-6;
    xt[i] = x0[i] * exp(-k * t);
  }
  double sum = xt.reduce((a, b) => a + b);
  if (sum < 1e-30) return List.filled(n, 1.0 / n);
  for (int i = 0; i < n; i++) xt[i] /= sum;
  return xt;
}

// ── Legacy API (backward compatible) ─────────────────────
/// Kept for backward compatibility with existing SolverPage calls.
List<String> topAccordNames(
  List<Compound> compounds,
  List<double> moleFractions, {
  int n = 3,
}) {
  final List<AccordScore> scores = scoreAccords(
    compounds: compounds,
    moleFractions: moleFractions,
  );
  return scores.take(n).map((s) => s.accord.name).toList();
}

Map<String, double> accordBalance(
  List<Compound> compounds,
  List<double> moleFractions,
) {
  final result = analyzeAccords(
    compounds: compounds,
    moleFractions: moleFractions,
  );
  return result.balance;
}
