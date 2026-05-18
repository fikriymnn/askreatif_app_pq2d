// ============================================================
// ASKREATIF Perfumery Engine
// PW-Integrated ROI Model + Compound Profile Widget
// ============================================================
// Replaces the heuristic ROI model with PerfumersWorld-
// calibrated Relative Impact (RI) values and Odour Life data.
// ============================================================

import 'dart:math';
import 'package:flutter/material.dart';
import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/colorTheme.dart';
import 'pw_roi_database.dart';

// ════════════════════════════════════════════════════════════
// ENHANCED ROI CALCULATION (PW-calibrated)
// ════════════════════════════════════════════════════════════

/// PW-corrected ROI vector for all compounds in a mixture.
/// Integrates thermodynamic mole fraction × Psat with
/// perceptual RI calibration from PerfumersWorld.
List<double> pwRoiVector({
  required List<Compound> compounds,
  required List<double> moleFractions,
}) {
  // Guard: use the shorter length to prevent RangeError
  final int n =
      compounds.length < moleFractions.length
          ? compounds.length
          : moleFractions.length;

  return [
    for (int i = 0; i < n; i++)
      pwCorrectedRoi(
        name: compounds[i].name,
        moleFraction: moleFractions[i],
        psat: compounds[i].Psat,
        mw: compounds[i].MW,
        threshold: compounds[i].Thr,
      ),
  ];
}

/// Normalised PW-ROI (sums to 1.0) — perceptual dominance share.
List<double> normalisedPwRoi({
  required List<Compound> compounds,
  required List<double> moleFractions,
}) {
  if (compounds.isEmpty || moleFractions.isEmpty) {
    return List.filled(compounds.length, 0.0);
  }

  final int n =
      compounds.length < moleFractions.length
          ? compounds.length
          : moleFractions.length;

  final List<Compound> safeComps = compounds.sublist(0, n);
  final List<double> safeFracs = moleFractions.sublist(0, n);

  final List<double> roi = pwRoiVector(
    compounds: safeComps,
    moleFractions: safeFracs,
  );
  final double total = roi.fold(0.0, (a, b) => a + b);
  if (total <= 0) return List.filled(n, 0.0);
  return roi.map((v) => v / total).toList();
}

// ════════════════════════════════════════════════════════════
// DYNAMIC NOTE CLASSIFICATION — PW-calibrated
// ════════════════════════════════════════════════════════════

class PwNoteState {
  final String compoundName;
  final String abcClass; // A–Z PerfumersWorld class
  final String abcLabel; // Full ABC description
  final double relativeImpact; // PW RI value
  final double odourLifeHours; // PW odour life
  final double pwRoiWeight; // 0–1 normalised PW-ROI
  final String classicNote; // Top / Heart / Base
  final String perceptualRole; // Dominant / Supporting / Trace
  final String colorHex;

  const PwNoteState({
    required this.compoundName,
    required this.abcClass,
    required this.abcLabel,
    required this.relativeImpact,
    required this.odourLifeHours,
    required this.pwRoiWeight,
    required this.classicNote,
    required this.perceptualRole,
    required this.colorHex,
  });

  String get fullLabel => '$classicNote · $perceptualRole  [$abcClass]';
}

List<PwNoteState> classifyWithPwData({
  required List<Compound> compounds,
  required List<double> moleFractions,
  double t = 0.0,
}) {
  // Clamp to shorter list
  final int n =
      compounds.length < moleFractions.length
          ? compounds.length
          : moleFractions.length;

  final List<Compound> safeComps = compounds.sublist(0, n);
  final List<double> safeFracs = moleFractions.sublist(0, n);

  final List<double> nRoi = normalisedPwRoi(
    compounds: safeComps,
    moleFractions: safeFracs,
  );

  return [
    for (int i = 0; i < n; i++)
      _buildPwNoteState(c: safeComps[i], normRoi: nRoi[i], t: t),
  ];
}

PwNoteState _buildPwNoteState({
  required Compound c,
  required double normRoi,
  required double t,
}) {
  final PwRoiEntry? entry = pwRoiFor(c.name);
  final double ri = entry?.relativeImpact ?? 100.0;
  final double ol = entry?.odourLifeHours ?? 12.0;
  final String abc = entry?.abcClass ?? '?';
  final String abcLbl = entry?.abcLabel ?? c.name;

  // Classic note classification — odour life based (more accurate
  // than pure Psat when PW data is available)
  String classic;
  if (ol <= 8.0) {
    classic = 'Top';
  } else if (ol <= 48.0) {
    classic = 'Heart';
  } else {
    classic = 'Base';
  }

  // Perceptual role by PW-ROI dominance
  String role;
  if (normRoi >= 0.25)
    role = 'Dominant';
  else if (normRoi >= 0.08)
    role = 'Supporting';
  else if (normRoi >= 0.01)
    role = 'Accent';
  else
    role = 'Subliminal';

  // Colour by classic note
  String hex;
  switch (classic) {
    case 'Top':
      hex = '#4FC3F7';
      break;
    case 'Heart':
      hex = '#A5D6A7';
      break;
    default:
      hex = '#CE93D8';
      break;
  }

  return PwNoteState(
    compoundName: c.name,
    abcClass: abc,
    abcLabel: abcLbl,
    relativeImpact: ri,
    odourLifeHours: ol,
    pwRoiWeight: normRoi,
    classicNote: classic,
    perceptualRole: role,
    colorHex: hex,
  );
}

// ════════════════════════════════════════════════════════════
// COMPOUND PROFILE WIDGET — shows ABC, RI, OL, ROI bar
// ════════════════════════════════════════════════════════════

class CompoundProfileCard extends StatelessWidget {
  final Compound compound;
  final double moleFraction;
  final double pwRoiNorm; // 0–1 normalised

  const CompoundProfileCard({
    Key? key,
    required this.compound,
    required this.moleFraction,
    required this.pwRoiNorm,
  }) : super(key: key);

  @override
  Widget build(BuildContext context) {
    final PwRoiEntry? entry = pwRoiFor(compound.name);
    final double ri = entry?.relativeImpact ?? 100.0;
    final double ol = entry?.odourLifeHours ?? 12.0;
    final String abc = entry?.abcClass ?? '?';
    final String lbl = entry?.abcLabel ?? 'Unknown class';
    final String src = entry?.source ?? 'EST';

    // Impact strength label (inverse of RI)
    String impactLabel;
    if (ri <= 10)
      impactLabel = 'ULTRA HIGH';
    else if (ri <= 50)
      impactLabel = 'HIGH';
    else if (ri <= 150)
      impactLabel = 'MEDIUM';
    else if (ri <= 500)
      impactLabel = 'LOW';
    else
      impactLabel = 'VERY LOW';

    // Note tier by odour life
    String noteTier;
    if (ol <= 8)
      noteTier = 'TOP';
    else if (ol <= 48)
      noteTier = 'HEART';
    else
      noteTier = 'BASE';

    // Tier colour
    Color tierColor;
    switch (noteTier) {
      case 'TOP':
        tierColor = const Color(0xFF4FC3F7);
        break;
      case 'HEART':
        tierColor = const Color(0xFFA5D6A7);
        break;
      default:
        tierColor = const Color(0xFFCE93D8);
        break;
    }

    // RI normalised bar (log scale, Linalool=100 as midpoint)
    final double riBar = (1.0 - (log(ri.clamp(1, 10000)) / log(10000))).clamp(
      0.0,
      1.0,
    );

    return Container(
      margin: const EdgeInsets.only(bottom: 10),
      padding: const EdgeInsets.all(12),
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: tierColor.withOpacity(0.3)),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // ── Header row ───────────────────────────────────
          Row(
            children: [
              // ABC badge
              Container(
                width: 26,
                height: 26,
                alignment: Alignment.center,
                decoration: BoxDecoration(
                  color: tierColor.withOpacity(0.15),
                  borderRadius: BorderRadius.circular(4),
                  border: Border.all(color: tierColor.withOpacity(0.5)),
                ),
                child: Text(
                  abc,
                  style: TextStyle(
                    color: tierColor,
                    fontSize: 12,
                    fontWeight: FontWeight.bold,
                    fontFamily: 'Courier',
                  ),
                ),
              ),
              const SizedBox(width: 8),
              Expanded(
                child: Column(
                  crossAxisAlignment: CrossAxisAlignment.start,
                  children: [
                    Text(
                      compound.name,
                      style: TextStyle(
                        color: AppTheme.textPrimary,
                        fontSize: 12,
                        fontWeight: FontWeight.w600,
                      ),
                    ),
                    Text(
                      lbl,
                      style: TextStyle(color: AppTheme.textMuted, fontSize: 9),
                    ),
                  ],
                ),
              ),
              // Note tier chip
              Container(
                padding: const EdgeInsets.symmetric(horizontal: 7, vertical: 3),
                decoration: BoxDecoration(
                  color: tierColor.withOpacity(0.1),
                  borderRadius: BorderRadius.circular(3),
                  border: Border.all(color: tierColor.withOpacity(0.4)),
                ),
                child: Text(
                  noteTier,
                  style: TextStyle(
                    color: tierColor,
                    fontSize: 9,
                    fontWeight: FontWeight.bold,
                    letterSpacing: 0.8,
                  ),
                ),
              ),
            ],
          ),

          const SizedBox(height: 10),

          // ── Data grid ────────────────────────────────────
          Row(
            children: [
              _dataCell('RI', ri.toStringAsFixed(0), sub: 'vs Linalool=100'),
              _dataCell('OL', '${ol.toStringAsFixed(0)} h', sub: 'Odour life'),
              _dataCell('IMPACT', impactLabel, sub: 'Strength tier'),
              _dataCell('SOURCE', src, sub: 'Data source'),
            ],
          ),

          const SizedBox(height: 8),

          // ── Relative Impact bar ───────────────────────────
          _labeledBar(
            label: 'Relative Impact',
            value: riBar,
            color: tierColor,
            hint: 'Higher bar = stronger material',
          ),

          const SizedBox(height: 4),

          // ── PW-ROI dominance bar ──────────────────────────
          _labeledBar(
            label: 'Perceptual Dominance',
            value: pwRoiNorm.clamp(0.0, 1.0),
            color: AppTheme.green400,
            hint: '${(pwRoiNorm * 100).toStringAsFixed(1)}% of blend impact',
          ),
        ],
      ),
    );
  }

  Widget _dataCell(String label, String value, {String? sub}) {
    return Expanded(
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            label,
            style: TextStyle(
              color: AppTheme.textMuted,
              fontSize: 8,
              letterSpacing: 0.8,
            ),
          ),
          Text(
            value,
            style: TextStyle(
              color: AppTheme.textPrimary,
              fontSize: 11,
              fontFamily: 'Courier',
              fontWeight: FontWeight.w600,
            ),
          ),
          if (sub != null)
            Text(sub, style: TextStyle(color: AppTheme.textMuted, fontSize: 8)),
        ],
      ),
    );
  }

  Widget _labeledBar({
    required String label,
    required double value,
    required Color color,
    String? hint,
  }) {
    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Row(
          mainAxisAlignment: MainAxisAlignment.spaceBetween,
          children: [
            Text(
              label,
              style: TextStyle(color: AppTheme.textMuted, fontSize: 9),
            ),
            if (hint != null)
              Text(
                hint,
                style: TextStyle(color: AppTheme.textMuted, fontSize: 9),
              ),
          ],
        ),
        const SizedBox(height: 3),
        Stack(
          children: [
            Container(
              height: 5,
              decoration: BoxDecoration(
                color: AppTheme.bg,
                borderRadius: BorderRadius.circular(2),
              ),
            ),
            FractionallySizedBox(
              widthFactor: value,
              child: Container(
                height: 5,
                decoration: BoxDecoration(
                  color: color.withOpacity(0.7),
                  borderRadius: BorderRadius.circular(2),
                ),
              ),
            ),
          ],
        ),
      ],
    );
  }
}

// ════════════════════════════════════════════════════════════
// FULL BLEND PROFILE PANEL
// ════════════════════════════════════════════════════════════
/// Shows all compound profile cards with PW data for a mixture.
class BlendProfilePanel extends StatelessWidget {
  final List<Compound> compounds;
  final List<double> moleFractions;

  const BlendProfilePanel({
    Key? key,
    required this.compounds,
    required this.moleFractions,
  }) : super(key: key);

  @override
  Widget build(BuildContext context) {
    final List<double> nRoi = normalisedPwRoi(
      compounds: compounds,
      moleFractions: moleFractions,
    );

    final int n = nRoi.length;

    // Sort only over the safe range
    final List<int> order = List.generate(n, (i) => i)
      ..sort((a, b) => nRoi[b].compareTo(nRoi[a]));

    return Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Padding(
          padding: const EdgeInsets.fromLTRB(0, 0, 0, 12),
          child: Row(
            children: [
              Container(width: 3, height: 18, color: AppTheme.green400),
              const SizedBox(width: 10),
              Text(
                'COMPOUND PROFILES  ·  PW CALIBRATED',
                style: TextStyle(
                  color: AppTheme.textSecondary,
                  fontSize: 10,
                  letterSpacing: 1.5,
                  fontWeight: FontWeight.w600,
                ),
              ),
            ],
          ),
        ),
        for (final i in order)
          CompoundProfileCard(
            compound: compounds[i],
            moleFraction: moleFractions[i],
            pwRoiNorm: nRoi[i],
          ),
        Padding(
          padding: const EdgeInsets.only(top: 6),
          child: Text(
            'RI = Relative Impact vs Linalool (100).  '
            'OL = Odour Life (hrs).  '
            'Sources: PW = PerfumersWorld, LIT = Literature, EST = Estimated.',
            style: TextStyle(
              color: AppTheme.textMuted,
              fontSize: 8.5,
              height: 1.5,
            ),
          ),
        ),
      ],
    );
  }
}

// ════════════════════════════════════════════════════════════
// ABC WHEEL OVERVIEW WIDGET
// ════════════════════════════════════════════════════════════
/// Quick reference view showing the ABC class distribution
/// of the current blend with ROI weight per class.
class AbcWheelWidget extends StatelessWidget {
  final List<Compound> compounds;
  final List<double> moleFractions;

  const AbcWheelWidget({
    Key? key,
    required this.compounds,
    required this.moleFractions,
  }) : super(key: key);

  @override
  Widget build(BuildContext context) {
    // Guard: nRoi length may be shorter than compounds if lists were mismatched
    final List<double> nRoi = normalisedPwRoi(
      compounds: compounds,
      moleFractions: moleFractions,
    );

    // Use nRoi.length as the loop bound — it is always the safe minimum
    final int n = nRoi.length;

    final Map<String, double> classShare = {};
    for (int i = 0; i < n; i++) {
      // ← n, not compounds.length
      final String cls = abcClassOf(compounds[i].name);
      classShare[cls] = (classShare[cls] ?? 0.0) + nRoi[i];
    }

    // Sort descending
    final sorted =
        classShare.entries.toList()..sort((a, b) => b.value.compareTo(a.value));

    const Map<String, String> abcNames = {
      'A': 'Aliphatic',
      'B': 'Iceberg',
      'C': 'Citrus',
      'D': 'Dairy',
      'E': 'Edible',
      'F': 'Fruit',
      'G': 'Green',
      'H': 'Herb',
      'I': 'Iris',
      'J': 'Jasmine',
      'K': 'Konifer',
      'L': 'Linalool',
      'M': 'Muguet',
      'N': 'Narcotic',
      'O': 'Orchid',
      'P': 'Phenol',
      'Q': 'Oriental',
      'R': 'Rose',
      'S': 'Spice',
      'T': 'Tar',
      'U': 'Animal',
      'V': 'Vanilla',
      'W': 'Wood',
      'X': 'Musk',
      'Y': 'Earthy',
      'Z': 'Solvent',
      '?': 'Unknown',
    };

    return Container(
      padding: const EdgeInsets.all(14),
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: AppTheme.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Row(
            children: [
              const Icon(
                Icons.auto_awesome_mosaic_outlined,
                color: Color(0xFF4FC3F7),
                size: 13,
              ),
              const SizedBox(width: 6),
              Text(
                'ABC CLASS DISTRIBUTION',
                style: TextStyle(
                  color: const Color(0xFF4FC3F7),
                  fontSize: 10,
                  letterSpacing: 1.5,
                  fontWeight: FontWeight.w600,
                ),
              ),
            ],
          ),
          const SizedBox(height: 12),
          for (final e in sorted)
            Padding(
              padding: const EdgeInsets.only(bottom: 6),
              child: Row(
                children: [
                  // Class badge
                  Container(
                    width: 22,
                    height: 22,
                    alignment: Alignment.center,
                    decoration: BoxDecoration(
                      color: AppTheme.bg,
                      borderRadius: BorderRadius.circular(3),
                      border: Border.all(color: AppTheme.border),
                    ),
                    child: Text(
                      e.key,
                      style: TextStyle(
                        color: AppTheme.textPrimary,
                        fontSize: 11,
                        fontFamily: 'Courier',
                        fontWeight: FontWeight.bold,
                      ),
                    ),
                  ),
                  const SizedBox(width: 8),
                  SizedBox(
                    width: 70,
                    child: Text(
                      abcNames[e.key] ?? e.key,
                      style: TextStyle(
                        color: AppTheme.textSecondary,
                        fontSize: 10,
                      ),
                    ),
                  ),
                  Expanded(
                    child: Stack(
                      children: [
                        Container(
                          height: 6,
                          decoration: BoxDecoration(
                            color: AppTheme.bg,
                            borderRadius: BorderRadius.circular(2),
                          ),
                        ),
                        FractionallySizedBox(
                          widthFactor: e.value.clamp(0.0, 1.0),
                          child: Container(
                            height: 6,
                            decoration: BoxDecoration(
                              color: AppTheme.green400.withOpacity(0.7),
                              borderRadius: BorderRadius.circular(2),
                            ),
                          ),
                        ),
                      ],
                    ),
                  ),
                  const SizedBox(width: 8),
                  SizedBox(
                    width: 38,
                    child: Text(
                      '${(e.value * 100).toStringAsFixed(1)}%',
                      style: TextStyle(
                        color: AppTheme.green400,
                        fontSize: 10,
                        fontFamily: 'Courier',
                      ),
                      textAlign: TextAlign.right,
                    ),
                  ),
                ],
              ),
            ),
        ],
      ),
    );
  }
}
