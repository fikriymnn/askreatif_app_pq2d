// ============================================================
// ASKREATIF Perfumery Engine
// Unified Perceptual State Framework
// Replaces the 4 conflicting note classification systems with
// one coherent, continuous, time-aware perceptual model.
// ============================================================
// Model basis:
//   P_i(t) = I_i(t) × σ_adapt(t) × σ_fatigue(t)
// where:
//   I_i(t)        = Stevens intensity from headspace at time t
//   σ_adapt(t)    = sensory adaptation decay [Weber-Fechner]
//   σ_fatigue(t)  = olfactory fatigue accumulation
//
// Note classification is CONTINUOUS, not discrete:
//   noteScore_i(t) ∈ [0,1] for Top, Heart, Base simultaneously
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/Evaporation.dart';
import 'package:askreatif_app/PsychoOdor.dart';
import 'package:askreatif_app/Roi.dart';

// ── Continuous Note Score ─────────────────────────────────
/// A compound can simultaneously contribute to multiple notes.
/// Scores sum to ≤ 1.0 per compound (not across compounds).
class ContinuousNoteScore {
  final double top; // 0–1
  final double heart; // 0–1
  final double base; // 0–1

  const ContinuousNoteScore({
    required this.top,
    required this.heart,
    required this.base,
  });

  /// Dominant note label (highest score)
  String get dominant {
    if (top >= heart && top >= base) return 'Top';
    if (heart >= base) return 'Heart';
    return 'Base';
  }

  /// Secondary note (second highest)
  String get secondary {
    final List<MapEntry<String, double>> entries = [
      MapEntry('Top', top),
      MapEntry('Heart', heart),
      MapEntry('Base', base),
    ]..sort((a, b) => b.value.compareTo(a.value));
    return entries[1].key;
  }
}

// ── Sensory Adaptation Model ──────────────────────────────
/// Weber-Fechner adaptation: prolonged exposure reduces
/// perceived intensity. Recovery follows exponential kinetics.
///
/// σ_adapt(t) = σ_min + (1 - σ_min) × exp(-λ_adapt × t)
///
/// λ_adapt depends on compound class (aldehydes adapt faster
/// than musks, for instance).
double sensoryadaptation({
  required Compound c,
  required double t, // hours of exposure
  double sigmaMin = 0.15,
}) {
  // Adaptation rate: faster for high-Psat volatiles
  // Literature: top notes adapt in ~5-15 min, base notes ~2-4h
  final double kAdapt = _adaptationRate(c);
  final double sigma = sigmaMin + (1.0 - sigmaMin) * exp(-kAdapt * t);
  return sigma.clamp(sigmaMin, 1.0);
}

double _adaptationRate(Compound c) {
  // λ_adapt [h⁻¹]: empirically derived from psychophysics literature
  // High Psat → fast adaptation (top notes)
  // Low Psat → slow adaptation (base notes)
  if (c.Psat >= 10.0) return 3.0; // top: full adapt ~20 min
  if (c.Psat >= 1.0) return 1.2; // top-mid: ~50 min
  if (c.Psat >= 0.1) return 0.5; // heart: ~2h
  if (c.Psat >= 0.01) return 0.15; // heart-base: ~7h
  return 0.05; // base: ~20h
}

// ── Olfactory Fatigue ─────────────────────────────────────
/// Cumulative fatigue from total mixture intensity over time.
/// Reduces perceived intensity of all compounds after
/// sustained high-intensity exposure.
///
/// σ_fatigue(t) = exp(-κ × ∫₀ᵗ I_total(τ) dτ)
///
/// Approximated as: σ_fatigue(t) ≈ exp(-κ × I_mean × t)
double olfactoryFatigue({
  required double cumulativeIntensity, // integral of I_total over time
  double kappa = 0.008, // fatigue sensitivity
}) {
  final double sigma = exp(-kappa * cumulativeIntensity);
  return sigma.clamp(0.05, 1.0);
}

// ── Unified Perceptual State ──────────────────────────────
class UnifiedPerceptualState {
  final double time; // hours
  final List<double> perceivedIntensity; // per compound, post-adaptation
  final List<ContinuousNoteScore> noteScores;
  final List<double> perceptualDominance; // normalised perceived contribution
  final double totalPerceivedIntensity;
  final double adaptationLevel; // 0=fresh, 1=fully adapted
  final double fatigueLevel; // 0=fresh, 1=fatigued

  const UnifiedPerceptualState({
    required this.time,
    required this.perceivedIntensity,
    required this.noteScores,
    required this.perceptualDominance,
    required this.totalPerceivedIntensity,
    required this.adaptationLevel,
    required this.fatigueLevel,
  });
}

// ── Continuous Note Score from Physics ───────────────────
/// Maps compound properties to continuous note scores using
/// smooth sigmoid functions rather than hard cutoffs.
///
/// Top score:   logistic centered at log(Psat) = 1.0 (Psat~10)
/// Heart score: bell curve centered at log(Psat) = 0.0 (Psat~1)
/// Base score:  logistic centered at log(Psat) = -1.5 (Psat~0.03)
ContinuousNoteScore continuousNoteScore(Compound c) {
  final double lp = log(c.Psat.clamp(1e-6, 1e6)) / ln10;

  // Sigmoid functions
  double _sig(double x, double center, double width) {
    return 1.0 / (1.0 + exp(-(x - center) / width));
  }

  // Top: high Psat compounds
  final double topScore = _sig(lp, 0.7, 0.6);

  // Base: low Psat compounds
  final double baseScore = 1.0 - _sig(lp, -1.2, 0.6);

  // Heart: bell curve (neither extreme)
  final double heartScore = (1.0 - topScore) * (1.0 - baseScore) * 2.5;

  // Normalise so max = 1 per compound
  final num maxS = [topScore, heartScore.clamp(0, 1), baseScore].reduce(max);
  if (maxS <= 0) {
    return const ContinuousNoteScore(top: 0, heart: 1, base: 0);
  }

  return ContinuousNoteScore(
    top: topScore / maxS,
    heart: heartScore.clamp(0, 1) / maxS,
    base: baseScore / maxS,
  );
}

// ── Full Unified Perceptual State Computation ─────────────
/// Computes the full unified perceptual state at time t.
/// This is the SINGLE authoritative note/perception source.
UnifiedPerceptualState computeUnifiedPerceptualState({
  required List<Compound> compounds,
  required List<double> x0, // initial mole fractions
  required List<double> Function(List<double>) gammaFn,
  required double t, // hours
  double cumulativeIntensity = 0.0, // ∫ I_total dt (pass from trajectory)
}) {
  final int n = compounds.length;

  // 1. Physics: mole fractions at time t
  final List<double> xt = moleFractionsAtTime(
    compounds: compounds,
    x0: x0,
    t: t,
  );

  // 2. Activity coefficients
  final List<double> gammas = gammaFn(xt);

  // 3. Headspace concentrations [mg/m³]
  final List<double> C = headspaceConcentration(
    compounds: compounds,
    xt: xt,
    gammas: gammas,
  );

  // 4. Raw Stevens intensity per compound
  final List<double> rawI = intensityVector(
    compounds: compounds,
    headspaceConc: C,
  );

  // 5. Sensory adaptation per compound
  final List<double> adaptFactors =
      compounds.map((c) => sensoryadaptation(c: c, t: t)).toList();
  final double avgAdapt = adaptFactors.fold(0.0, (a, b) => a + b) / n;

  // 6. Olfactory fatigue (mixture-level)
  final double fatigueF = olfactoryFatigue(
    cumulativeIntensity: cumulativeIntensity,
  );

  // 7. Perceived intensity = raw × adaptation × fatigue
  final List<double> perceivedI = [
    for (int i = 0; i < n; i++)
      (rawI[i] * adaptFactors[i] * fatigueF).clamp(0.0, double.infinity),
  ];

  // 8. Total perceived intensity
  final double totalI = perceivedI.fold(0.0, (a, b) => a + b);

  // 9. Perceptual dominance (normalised)
  final List<double> dominance =
      totalI > 1e-30
          ? perceivedI.map((v) => v / totalI).toList()
          : List.filled(n, 1.0 / n);

  // 10. Continuous note scores
  final List<ContinuousNoteScore> noteScores =
      compounds.map(continuousNoteScore).toList();

  // 11. Fatigue level (0=fresh, 1=fatigued)
  final double fatigueLevel = 1.0 - fatigueF;

  return UnifiedPerceptualState(
    time: t,
    perceivedIntensity: perceivedI,
    noteScores: noteScores,
    perceptualDominance: dominance,
    totalPerceivedIntensity: totalI,
    adaptationLevel: 1.0 - avgAdapt,
    fatigueLevel: fatigueLevel,
  );
}

// ── Perceptual Trajectory ─────────────────────────────────
/// Full timeline of unified perceptual states.
List<UnifiedPerceptualState> computePerceptualTrajectory({
  required List<Compound> compounds,
  required List<double> x0,
  required List<double> Function(List<double>) gammaFn,
  List<double>? timePoints,
}) {
  final List<double> ts = timePoints ?? defaultTimePoints();
  final List<UnifiedPerceptualState> trajectory = [];
  double cumulativeI = 0.0;
  double prevI = 0.0;
  double prevT = 0.0;

  for (final t in ts) {
    final UnifiedPerceptualState state = computeUnifiedPerceptualState(
      compounds: compounds,
      x0: x0,
      gammaFn: gammaFn,
      t: t,
      cumulativeIntensity: cumulativeI,
    );
    // Trapezoidal integration of I_total
    final double dt = t - prevT;
    cumulativeI += 0.5 * (prevI + state.totalPerceivedIntensity) * dt;
    prevI = state.totalPerceivedIntensity;
    prevT = t;
    trajectory.add(state);
  }
  return trajectory;
}
