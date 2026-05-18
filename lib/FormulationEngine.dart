// ============================================================
// ASKREATIF Perfumery Engine — Phase 8
// Real-World Formulation Engine
// ============================================================
// Converts mole fractions → mass fractions → gram weights.
// Supports batch scaling, density correction, mL output,
// concentration percentage, and lab-ready formulation sheets.
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';

// ── Default densities (g/mL) at 25 °C ────────────────────
const Map<String, double> _defaultDensities = {
  'Water': 1.000,
  'Ethanol': 0.789,
  'DPG': 1.023,
  'Geraniol': 0.880,
  'Linalool': 0.868,
  'Citronellol': 0.857,
  'Limonene': 0.842,
  'Nerol': 0.877,
  'Benzyl Acetate': 1.059,
  'Eugenol': 1.066,
  'Cinnamaldehyde': 1.050,
  'Vanillin': 1.056,
  'Coumarin': 0.935,
  'Menthol': 0.890,
  'Ambroxan': 1.000,
  'Galaxolide': 0.990,
  'Cashmeran': 0.960,
  'Iso E Super': 0.965,
  'Patchouli Alcohol': 0.983,
  'Javanol': 0.952,
  'Romandolide': 0.960,
  'Hedione': 0.980,
  'Isoamyl Acetate': 0.876,
  'Hexyl Acetate': 0.878,
  'Heptyl Acetate': 0.875,
  'Octyl Acetate': 0.872,
  'Benzaldehyde': 1.044,
  'a- ionone': 0.933,
  'Cedryl Acetate': 1.000,
  'Cis-3-Hexenyl Acetate': 0.898,
  'CIS-6-NONEAL': 0.850,
  'Sandalmysore Core': 0.960,
  'Bacdanol': 0.940,
  'Frukton': 0.970,
  'Helvetolida': 0.965,
  'Evernyl': 1.050,
  'b-caryophyllenol': 0.970,
  'Caryophyllene Oxide': 1.000,
  'Macrolide': 0.965,
  'Oenanthic Ether': 0.875,
  'Dodecanal': 0.830,
};

double densityForCompound(Compound c) {
  return _defaultDensities[c.name] ?? 1.00;
}

// ── Formulation Entry ─────────────────────────────────────
class FormulationEntry {
  final String name;
  final bool isSolvent;
  final double moleFraction;
  final double massFraction; // 0–1
  final double massGram; // grams in batch
  final double volumeMl; // mL in batch
  final double concentrationPct; // w/w percent (0–100)
  final double density; // g/mL
  final double MW;

  const FormulationEntry({
    required this.name,
    required this.isSolvent,
    required this.moleFraction,
    required this.massFraction,
    required this.massGram,
    required this.volumeMl,
    required this.concentrationPct,
    required this.density,
    required this.MW,
  });
}

// ── Batch Formulation Builder ─────────────────────────────
/// Generates a full lab-ready formulation for a given batch mass.
/// Returns one [FormulationEntry] per compound (fragrance + solvent).
List<FormulationEntry> buildFormulation({
  required List<Compound> fragranceComps,
  required List<double> fragranceFractions,
  required List<Compound> solventComps,
  required List<double> solventRatios,
  required double totalSolventFraction,
  required double batchGram,
}) {
  if (fragranceComps.isEmpty) return [];

  // Combine all components
  final List<Compound> allComps = [...fragranceComps, ...solventComps];
  final List<double> allFractions = [
    ...fragranceFractions,
    ...solventRatios.map((r) => r * totalSolventFraction),
  ];

  final Set<String> solventNames = solventComps.map((c) => c.name).toSet();

  // MW-weighted mass fractions: w_i ∝ x_i × MW_i
  final List<double> wMW = [
    for (int i = 0; i < allComps.length; i++) allFractions[i] * allComps[i].MW,
  ];

  double totalWMW = wMW.fold(0.0, (a, b) => a + b);
  if (totalWMW <= 0) totalWMW = 1.0;

  final List<FormulationEntry> entries = [];

  for (int i = 0; i < allComps.length; i++) {
    final double massFrac = wMW[i] / totalWMW;
    final double massG = massFrac * batchGram;
    final double density = densityForCompound(allComps[i]);
    final double volMl = density > 0 ? massG / density : 0.0;

    entries.add(
      FormulationEntry(
        name: allComps[i].name,
        isSolvent: solventNames.contains(allComps[i].name),
        moleFraction: allFractions[i],
        massFraction: massFrac,
        massGram: massG,
        volumeMl: volMl,
        concentrationPct: massFrac * 100.0,
        density: density,
        MW: allComps[i].MW,
      ),
    );
  }

  return entries;
}

// ── Plain-Text Lab Formulation Sheet ─────────────────────
/// Generates a human-readable lab formulation sheet string.
String generateFormulationSheet({
  required List<FormulationEntry> entries,
  required double batchGram,
  required String solventSystem,
  String? formulaName,
}) {
  final StringBuffer sb = StringBuffer();
  final String sep = '─' * 72;

  sb.writeln(
    '╔══════════════════════════════════════════════════════════════════╗',
  );
  sb.writeln(
    '║  ASKREATIF Perfumery Engine — Lab Formulation Sheet              ║',
  );
  sb.writeln(
    '╚══════════════════════════════════════════════════════════════════╝',
  );
  sb.writeln('Formula : ${formulaName ?? "Unnamed Formula"}');
  sb.writeln('Batch   : ${batchGram.toStringAsFixed(2)} g');
  sb.writeln('Solvent : $solventSystem');
  sb.writeln('Generated: ${DateTime.now().toLocal()}');
  sb.writeln(sep);

  sb.writeln(
    '${'Component'.padRight(26)}'
    '${'Mol%'.padLeft(6)}  '
    '${'w/w%'.padLeft(6)}  '
    '${'Mass (g)'.padLeft(10)}  '
    '${'Vol (mL)'.padLeft(9)}  '
    'Type',
  );
  sb.writeln(sep);

  // Fragrance rows
  sb.writeln('[FRAGRANCE COMPOUNDS]');
  for (final e in entries.where((e) => !e.isSolvent)) {
    sb.writeln(
      '${e.name.padRight(26)}'
      '${(e.moleFraction * 100).toStringAsFixed(2).padLeft(6)}  '
      '${e.concentrationPct.toStringAsFixed(2).padLeft(6)}  '
      '${e.massGram.toStringAsFixed(4).padLeft(10)}  '
      '${e.volumeMl.toStringAsFixed(4).padLeft(9)}  '
      'Fragrance',
    );
  }

  sb.writeln('');
  sb.writeln('[SOLVENT / CARRIER]');
  for (final e in entries.where((e) => e.isSolvent)) {
    sb.writeln(
      '${e.name.padRight(26)}'
      '${(e.moleFraction * 100).toStringAsFixed(2).padLeft(6)}  '
      '${e.concentrationPct.toStringAsFixed(2).padLeft(6)}  '
      '${e.massGram.toStringAsFixed(4).padLeft(10)}  '
      '${e.volumeMl.toStringAsFixed(4).padLeft(9)}  '
      'Solvent',
    );
  }

  sb.writeln(sep);

  final double totalMass = entries.fold(0.0, (a, e) => a + e.massGram);
  final double totalVol = entries.fold(0.0, (a, e) => a + e.volumeMl);

  sb.writeln(
    '${'TOTAL'.padRight(26)}'
    '${'100.00'.padLeft(6)}  '
    '${'100.00'.padLeft(6)}  '
    '${totalMass.toStringAsFixed(4).padLeft(10)}  '
    '${totalVol.toStringAsFixed(4).padLeft(9)}',
  );

  sb.writeln(sep);
  sb.writeln('Notes: Mole fractions from UNIFAC-optimised OV balancing.');
  sb.writeln('       Gram weights via MW-weighted conversion.');
  sb.writeln('       Volumes assume literature densities at 25 °C.');

  return sb.toString();
}

// ── Projected Longevity Estimate ─────────────────────────
/// Estimates hours until ~90 % of total odor value has decayed.
/// Uses OV-weighted effective first-order decay constant:
///   k_eff = Σ(w_i × k_i)   where w_i = OV_i / ΣOV
///   t_90% = -ln(0.10) / k_eff
double projectedLongevityHours({
  required List<Compound> compounds,
  required List<double> odorValues,
}) {
  if (compounds.isEmpty || odorValues.isEmpty) return 0.0;

  final double totalOV = odorValues.fold(0.0, (a, b) => a + b);
  if (totalOV <= 0) return 0.0;

  double kEff = 0.0;
  for (int i = 0; i < compounds.length; i++) {
    final double psat = compounds[i].Psat;
    final double mw = compounds[i].MW;
    if (psat <= 0 || mw <= 0) continue;

    // First-order evaporation constant from Graham-law model
    final double k = 0.0012 * psat / sqrt(mw);
    final double w = odorValues[i] / totalOV;
    kEff += w * k;
  }

  if (kEff <= 0) return 99.0;

  // t at which 90 % has decayed: N(t) = 0.10 × N(0)
  return -log(0.10) / kEff;
}
