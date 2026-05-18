// ============================================================
// ASKREATIF Perfumery Engine — Phase 6
// Accord Intelligence System
// ============================================================
// Classifies, recommends, and scores fragrance accords.
// Accords are emergent compound groupings with collective
// perceptual identity — not individual compound labels.
// ============================================================

import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/Roi.dart';

// ── Accord Archetype Definitions ──────────────────────────
class AccordProfile {
  final String name;
  final String icon;           // unicode emoji / icon hint
  final String description;
  final List<String> keyGroups;  // UNIFAC groups that contribute
  final List<String> keyCompounds; // benchmark compounds

  const AccordProfile({
    required this.name,
    required this.icon,
    required this.description,
    required this.keyGroups,
    required this.keyCompounds,
  });
}

const List<AccordProfile> accordArchetypes = [
  AccordProfile(
    name: 'Floral',
    icon: '🌸',
    description: 'Delicate, petal-like, rosy, jasmine, ylang character.',
    keyGroups: ['OH', 'ACH', 'ACOH', 'CH=C', 'ACCH2'],
    keyCompounds: ['Geraniol', 'Nerol', 'Citronellol', 'Linalool'],
  ),
  AccordProfile(
    name: 'Woody',
    icon: '🌲',
    description: 'Dry, cedar, sandalwood, vetiver, earthy warmth.',
    keyGroups: ['CH3', 'CH2', 'CH', 'C', 'C=C'],
    keyCompounds: ['Javanol', 'Sandalmysore Core', 'Cashmeran', 'Iso E Super'],
  ),
  AccordProfile(
    name: 'Amber / Oriental',
    icon: '🍯',
    description: 'Warm, resinous, sweet-balsamic, animalic depth.',
    keyGroups: ['COO', 'CH2COO', 'CH3COO', 'CH3CO', 'CH2CO'],
    keyCompounds: ['Ambroxan', 'Vanillin', 'Coumarin', 'Romandolide'],
  ),
  AccordProfile(
    name: 'Musk',
    icon: '🌫',
    description: 'Soft, skin-like, powdery, clean, diffusive.',
    keyGroups: ['CH2O', 'CH3COO', 'CH2COO', 'COO'],
    keyCompounds: ['Galaxolide', 'Ambroxan', 'Habanolide', 'Romandolide'],
  ),
  AccordProfile(
    name: 'Fresh / Citrus',
    icon: '🍃',
    description: 'Bright, airy, green, citrus, ozonic, uplifting.',
    keyGroups: ['CH=CH', 'CH2=CH', 'CH2=C', 'CH=C'],
    keyCompounds: ['Limonene', 'Linalool', 'Cis-3-Hexenyl Acetate', 'CIS-6-NONEAL'],
  ),
  AccordProfile(
    name: 'Gourmand',
    icon: '🍮',
    description: 'Sweet, edible, caramel, vanilla, lactonic.',
    keyGroups: ['CH3CO', 'CH2CO', 'CH3COO'],
    keyCompounds: ['Vanillin', 'Isoamyl Acetate', 'Benzyl Acetate', 'Coumarin'],
  ),
  AccordProfile(
    name: 'Aromatic / Spicy',
    icon: '🌿',
    description: 'Herbal, medicinal, warm spice, clove, cinnamon.',
    keyGroups: ['ACH', 'AC', 'ACOH', 'CHO'],
    keyCompounds: ['Eugenol', 'Cinnamaldehyde', 'Vanillin', 'Benzaldehyde'],
  ),
  AccordProfile(
    name: 'Fougère',
    icon: '☘️',
    description: 'Mossy, lavender-like, hay, coumarinic green.',
    keyGroups: ['COO', 'CH=CH', 'ACH', 'OH'],
    keyCompounds: ['Coumarin', 'Geraniol', 'Evernyl', 'Linalool'],
  ),
];

// ── Accord Scoring ─────────────────────────────────────────
class AccordScore {
  final AccordProfile accord;
  final double score;         // 0–1 normalised
  final double roiContrib;   // sum of ROI weights matching accord

  const AccordScore({
    required this.accord,
    required this.score,
    required this.roiContrib,
  });
}

/// Scores each accord archetype against the current mixture.
/// Uses ROI-weighted group counting.
List<AccordScore> scoreAccords({
  required List<Compound> compounds,
  required List<double> moleFractions,
}) {
  final List<double> nRoi = normalisedRoi(
    compounds: compounds,
    moleFractions: moleFractions,
  );

  List<AccordScore> scores = [];

  for (final accord in accordArchetypes) {
    double score = 0.0;

    for (int i = 0; i < compounds.length; i++) {
      for (final group in accord.keyGroups) {
        if (compounds[i].groups.containsKey(group)) {
          // ROI-weighted group contribution
          score += nRoi[i] * (compounds[i].groups[group]! / 10.0);
        }
      }
      // Bonus if compound is a known key compound for this accord
      if (accord.keyCompounds.contains(compounds[i].name)) {
        score += nRoi[i] * 0.3;
      }
    }

    scores.add(AccordScore(
      accord: accord,
      score: score.clamp(0.0, 1.0),
      roiContrib: score,
    ));
  }

  // Normalise to [0, 1] relative to highest scorer
  final double maxScore = scores
      .map((s) => s.roiContrib)
      .fold(0.0, (a, b) => a > b ? a : b);
  if (maxScore > 0) {
    scores = scores.map((s) => AccordScore(
      accord: s.accord,
      score: s.roiContrib / maxScore,
      roiContrib: s.roiContrib,
    )).toList();
  }

  scores.sort((a, b) => b.score.compareTo(a.score));
  return scores;
}

/// Top N accord names from the scoring.
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

/// Accord balance check: returns a map of dominant accordname → share.
Map<String, double> accordBalance(
  List<Compound> compounds,
  List<double> moleFractions,
) {
  final List<AccordScore> scores = scoreAccords(
    compounds: compounds,
    moleFractions: moleFractions,
  );
  final double totalScore = scores.fold(0.0, (a, s) => a + s.score);
  if (totalScore <= 0) return {};
  return {
    for (final s in scores.where((s) => s.score > 0.05))
      s.accord.name: s.score / totalScore
  };
}