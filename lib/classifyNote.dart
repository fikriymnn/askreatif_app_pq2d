import 'dart:math';

import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/OctonaryPainter.dart';
import 'package:askreatif_app/TernaryPainter.dart';
import 'package:askreatif_app/AccordIntelligence.dart';
import 'package:askreatif_app/ClimateModel.dart';
import 'package:flutter/material.dart';

// =========================================================
// THEME (copy of AppTheme for use in classifyNote)
// =========================================================
class _T {
  static const Color bg = Color(0xFF0D2818);
  static const Color bgSurface = Color(0xFF122B1D);
  static const Color bgCard = Color(0xFF173322);
  static const Color green900 = Color(0xFF1B5E20);
  static const Color green800 = Color(0xFF2E7D32);
  static const Color green600 = Color(0xFF43A047);
  static const Color green400 = Color(0xFF66BB6A);
  static const Color green200 = Color(0xFF81C784);
  static const Color green100 = Color(0xFFA5D6A7);
  static const Color textPrimary = Color(0xFFF0F7F0);
  static const Color textSecondary = Color(0xFF8BAE8B);
  static const Color textMuted = Color(0xFF4D7A54);
  static const Color border = Color(0xFF234D2E);
}

Map<String, List<int>> groupNoteIndices(List<Compound> comps) {
  // === TAHAP 1: Klasifikasi berdasarkan Psat standar industri ===
  Map<String, List<int>> groups = {"Top": [], "Middle": [], "Base": []};

  for (int i = 0; i < comps.length; i++) {
    double psat = comps[i].Psat;
    if (psat > 1.0) {
      groups["Top"]!.add(i);
    } else if (psat > 0.01) {
      groups["Middle"]!.add(i);
    } else {
      groups["Base"]!.add(i);
    }
  }

  // === TAHAP 2: Fallback jika ada group yang kosong ===
  // Urutkan semua index berdasarkan Psat descending
  List<int> sortedByPsat = List.generate(comps.length, (i) => i);
  sortedByPsat.sort((a, b) => comps[b].Psat.compareTo(comps[a].Psat));

  // Jika Top kosong → ambil senyawa dengan Psat tertinggi
  if (groups["Top"]!.isEmpty) {
    int candidate = sortedByPsat.firstWhere(
      (i) => !groups["Top"]!.contains(i),
      orElse: () => sortedByPsat.first,
    );
    // Pindahkan dari group lain ke Top
    groups["Middle"]!.remove(candidate);
    groups["Base"]!.remove(candidate);
    groups["Top"]!.add(candidate);
  }

  // Jika Base kosong → ambil senyawa dengan Psat terendah
  if (groups["Base"]!.isEmpty) {
    int candidate = sortedByPsat.lastWhere(
      (i) => !groups["Base"]!.contains(i),
      orElse: () => sortedByPsat.last,
    );
    groups["Middle"]!.remove(candidate);
    groups["Top"]!.remove(candidate);
    groups["Base"]!.add(candidate);
  }

  // Jika Middle kosong → ambil dari group terbesar (yang punya > 1 senyawa)
  if (groups["Middle"]!.isEmpty) {
    // Cari group dengan senyawa terbanyak
    String donorGroup =
        groups["Top"]!.length >= groups["Base"]!.length ? "Top" : "Base";

    if (groups[donorGroup]!.length > 1) {
      // Ambil senyawa paling "tengah" dari donor group
      List<int> donorSorted = List.from(groups[donorGroup]!);
      donorSorted.sort((a, b) => comps[b].Psat.compareTo(comps[a].Psat));

      // Dari Top ambil yang Psat-nya paling rendah (paling mendekati Middle)
      // Dari Base ambil yang Psat-nya paling tinggi (paling mendekati Middle)
      int candidate =
          donorGroup == "Top" ? donorSorted.last : donorSorted.first;

      groups[donorGroup]!.remove(candidate);
      groups["Middle"]!.add(candidate);
    }
  }

  return groups;
}

Widget buildNotesAndTernary(
  List<Compound> selectedList,
  List<double> optimizedFractions,
  List<double> ovValues,
) {
  final noteGroups = groupNoteIndices(selectedList);

  double topFrac = 0, middleFrac = 0, baseFrac = 0;
  for (int i in noteGroups["Top"]!) topFrac += optimizedFractions[i];
  for (int i in noteGroups["Middle"]!) middleFrac += optimizedFractions[i];
  for (int i in noteGroups["Base"]!) baseFrac += optimizedFractions[i];

  double totalFrac = topFrac + middleFrac + baseFrac;
  if (totalFrac <= 0) totalFrac = 1.0;

  double tTop = topFrac / totalFrac;
  double tMid = middleFrac / totalFrac;
  double tBase = baseFrac / totalFrac;

  return Column(
    crossAxisAlignment: CrossAxisAlignment.start,
    children: [
      SizedBox(height: 8),

      // Section header
      _SectionHeader(
        icon: Icons.local_florist_outlined,
        label: 'Notes Classification',
      ),
      SizedBox(height: 16),

      // Notes cards - horizontal layout
      Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Expanded(
            child: _buildNotesCard(
              '🌿 Top Notes',
              noteGroups["Top"]!,
              selectedList,
              optimizedFractions,
              ovValues,
              Color(0xFF1B3A1B),
              Color(0xFF2E7D32),
            ),
          ),
          SizedBox(width: 12),
          Expanded(
            child: _buildNotesCard(
              '🌸 Middle Notes',
              noteGroups["Middle"]!,
              selectedList,
              optimizedFractions,
              ovValues,
              Color(0xFF3A1B2E),
              Color(0xFF7B2D5E),
            ),
          ),
          SizedBox(width: 12),
          Expanded(
            child: _buildNotesCard(
              '🌳 Base Notes',
              noteGroups["Base"]!,
              selectedList,
              optimizedFractions,
              ovValues,
              Color(0xFF2E1F0A),
              Color(0xFF795548),
            ),
          ),
        ],
      ),

      SizedBox(height: 28),

      // Ternary diagram
      _SectionHeader(
        icon: Icons.change_history_outlined,
        label: 'Ternary Diagram (Top / Middle / Base)',
      ),
      SizedBox(height: 16),

      Container(
        padding: EdgeInsets.all(20),
        decoration: BoxDecoration(
          color: _T.bgCard,
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: _T.border),
        ),
        child: Column(
          children: [
            // Fraction chips row
            Row(
              mainAxisAlignment: MainAxisAlignment.center,
              children: [
                _FractionChip(
                  label: 'Top',
                  value: tTop,
                  color: Color(0xFF2E7D32),
                ),
                SizedBox(width: 12),
                _FractionChip(
                  label: 'Middle',
                  value: tMid,
                  color: Color(0xFF7B2D5E),
                ),
                SizedBox(width: 12),
                _FractionChip(
                  label: 'Base',
                  value: tBase,
                  color: Color(0xFF795548),
                ),
              ],
            ),
            SizedBox(height: 20),
            Center(
              child: CustomPaint(
                size: Size(320, 290),
                painter: TernaryPainter(top: tTop, middle: tMid, base: tBase),
              ),
            ),
            SizedBox(height: 12),
            Row(
              mainAxisAlignment: MainAxisAlignment.center,
              children: [
                _legendDot(Colors.green[600]!, "Top"),
                SizedBox(width: 16),
                _legendDot(Colors.pink[400]!, "Middle"),
                SizedBox(width: 16),
                _legendDot(Colors.brown[400]!, "Base"),
                SizedBox(width: 16),
                _legendDot(Colors.purple[400]!, "Composition"),
              ],
            ),
          ],
        ),
      ),
    ],
  );
}

Widget _buildNotesCard(
  String title,
  List<int> indices,
  List<Compound> comps,
  List<double> fractions,
  List<double> ovValues,
  Color bgColor,
  Color accentColor,
) {
  if (indices.isEmpty) {
    return Container(
      padding: EdgeInsets.all(14),
      decoration: BoxDecoration(
        color: bgColor.withOpacity(0.4),
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: accentColor.withOpacity(0.2)),
      ),
      child: Text(
        "$title\n—",
        style: TextStyle(color: _T.textMuted, fontSize: 12),
      ),
    );
  }

  return Container(
    padding: EdgeInsets.all(14),
    decoration: BoxDecoration(
      color: bgColor.withOpacity(0.5),
      borderRadius: BorderRadius.circular(4),
      border: Border.all(color: accentColor.withOpacity(0.3)),
    ),
    child: Column(
      crossAxisAlignment: CrossAxisAlignment.start,
      children: [
        Row(
          children: [
            Container(
              width: 8,
              height: 8,
              decoration: BoxDecoration(
                shape: BoxShape.circle,
                color: accentColor,
                boxShadow: [
                  BoxShadow(color: accentColor.withOpacity(0.5), blurRadius: 4),
                ],
              ),
            ),
            SizedBox(width: 8),
            Text(
              title,
              style: TextStyle(
                color: _T.textPrimary,
                fontSize: 12,
                fontWeight: FontWeight.w600,
                letterSpacing: 0.5,
              ),
            ),
          ],
        ),
        SizedBox(height: 10),
        ...indices.map(
          (i) => Container(
            margin: EdgeInsets.only(bottom: 8),
            padding: EdgeInsets.all(10),
            decoration: BoxDecoration(
              color: Colors.black.withOpacity(0.15),
              borderRadius: BorderRadius.circular(3),
            ),
            child: Column(
              crossAxisAlignment: CrossAxisAlignment.start,
              children: [
                Text(
                  comps[i].name,
                  style: TextStyle(
                    color: _T.textPrimary,
                    fontSize: 12,
                    fontWeight: FontWeight.w500,
                  ),
                ),
                SizedBox(height: 4),
                Row(
                  children: [
                    _DataBadge(
                      label: 'x',
                      value: fractions[i].toStringAsExponential(2),
                    ),
                    SizedBox(width: 6),
                    _DataBadge(
                      label: 'OV',
                      value: ovValues[i].toStringAsExponential(2),
                    ),
                    SizedBox(width: 6),
                    _DataBadge(label: 'Psat', value: '${comps[i].Psat}'),
                  ],
                ),
              ],
            ),
          ),
        ),
      ],
    ),
  );
}

Widget _legendDot(Color color, String label) {
  return Row(
    children: [
      Container(
        width: 10,
        height: 10,
        decoration: BoxDecoration(color: color, shape: BoxShape.circle),
      ),
      SizedBox(width: 5),
      Text(label, style: TextStyle(fontSize: 11, color: _T.textSecondary)),
    ],
  );
}

Widget buildOctonaryDiagram(
  List<Compound> selectedList,
  List<double> optimizedFractions,
  List<double> ovValues,
  List<Compound> solvents,
  List<double> solventRatios,
  List<double> solventOVs,
) {
  final noteGroups = groupNoteIndices(selectedList);

  Color _colorForCompound(int idx) {
    if (noteGroups["Top"]!.contains(idx)) return Color(0xFF43A047);
    if (noteGroups["Middle"]!.contains(idx)) return Color(0xFFAB47BC);
    return Color(0xFF8D6E63);
  }

  List<String> labels = [];
  List<double> rawOVs = [];
  List<Color> colors = [];

  for (String note in ["Top", "Middle", "Base"]) {
    for (int idx in noteGroups[note]!) {
      String name = selectedList[idx].name;
      if (name.length > 12) name = name.substring(0, 10) + "..";
      labels.add(name);
      rawOVs.add(ovValues[idx]);
      colors.add(_colorForCompound(idx));
    }
  }

  for (int i = 0; i < solvents.length; i++) {
    labels.add(solvents[i].name);
    rawOVs.add(solventOVs[i]);
    colors.add(Color(0xFF29B6F6));
  }

  while (labels.length < 8) {
    labels.add("-");
    rawOVs.add(0);
    colors.add(_T.border);
  }

  List<double> logOVs =
      rawOVs.map((v) {
        double val = (v <= 0 || v.isNaN || v.isInfinite) ? 1e-10 : v;
        return log(val);
      }).toList();

  double logMin = logOVs.reduce(min);
  double logMax = logOVs.reduce(max);
  double logRange = logMax - logMin;

  List<double> normalizedOVs =
      logOVs.map((lv) {
        if (logRange <= 0) return 0.5;
        return 0.05 + 0.95 * (lv - logMin) / logRange;
      }).toList();

  return Column(
    crossAxisAlignment: CrossAxisAlignment.start,
    children: [
      _SectionHeader(icon: Icons.radar, label: 'Octonary Diagram'),
      SizedBox(height: 16),

      Container(
        padding: EdgeInsets.all(20),
        decoration: BoxDecoration(
          color: _T.bgCard,
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: _T.border),
        ),
        child: Column(
          children: [
            // Legend
            Row(
              mainAxisAlignment: MainAxisAlignment.center,
              children: [
                _legendDot(Color(0xFF43A047), "Top"),
                SizedBox(width: 12),
                _legendDot(Color(0xFFAB47BC), "Middle"),
                SizedBox(width: 12),
                _legendDot(Color(0xFF8D6E63), "Base"),
                SizedBox(width: 12),
                _legendDot(Color(0xFF29B6F6), "Solvent"),
              ],
            ),
            SizedBox(height: 6),
            Text(
              'Nilai OV dinormalisasi dalam log-scale. Area lebih besar = kontribusi aroma lebih kuat.',
              textAlign: TextAlign.center,
              style: TextStyle(fontSize: 11, color: _T.textMuted),
            ),
            SizedBox(height: 16),
            Center(
              child: CustomPaint(
                size: Size(380, 380),
                painter: OctonaryPainter(
                  labels: labels,
                  values: normalizedOVs,
                  colors: colors,
                  rawValues: rawOVs,
                ),
              ),
            ),
          ],
        ),
      ),

      SizedBox(height: 12),

      // Summary table
      Container(
        padding: EdgeInsets.all(16),
        decoration: BoxDecoration(
          color: _T.bgCard,
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: _T.border),
        ),
        child: Column(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Text(
              'RINGKASAN OV',
              style: TextStyle(
                color: _T.green400,
                fontSize: 10,
                letterSpacing: 2,
                fontWeight: FontWeight.w600,
              ),
            ),
            SizedBox(height: 12),
            ...List.generate(labels.length, (i) {
              if (labels[i] == "-") return SizedBox.shrink();
              return Padding(
                padding: EdgeInsets.only(bottom: 8),
                child: Row(
                  children: [
                    Container(
                      width: 8,
                      height: 8,
                      decoration: BoxDecoration(
                        color: colors[i],
                        shape: BoxShape.circle,
                      ),
                    ),
                    SizedBox(width: 10),
                    Expanded(
                      flex: 3,
                      child: Text(
                        labels[i],
                        style: TextStyle(color: _T.textSecondary, fontSize: 12),
                      ),
                    ),
                    Container(
                      padding: EdgeInsets.symmetric(horizontal: 8, vertical: 3),
                      decoration: BoxDecoration(
                        color: _T.bg,
                        borderRadius: BorderRadius.circular(3),
                        border: Border.all(color: _T.border),
                      ),
                      child: Text(
                        'OV: ${rawOVs[i].toStringAsExponential(2)}',
                        style: TextStyle(
                          color: _T.textPrimary,
                          fontSize: 11,
                          fontFamily: 'Courier',
                        ),
                      ),
                    ),
                  ],
                ),
              );
            }),
          ],
        ),
      ),
    ],
  );
}
// ============================================================
// UI: Accord Intelligence + Climate Performance
// Tambahkan di bawah buildOctonaryDiagram()
// ============================================================

// ── Accord Intelligence Widget ────────────────────────────
Widget buildAccordIntelligence(
  List<Compound> selectedList,
  List<double> optimizedFractions,
  List<double> activityCoefficients,
) {
  final AccordIntelligenceResult result = analyzeAccords(
    compounds: selectedList,
    moleFractions: optimizedFractions,
    activityCoefficients: activityCoefficients,
  );

  return Column(
    crossAxisAlignment: CrossAxisAlignment.start,
    children: [
      SizedBox(height: 24),
      _SectionHeader(icon: Icons.auto_awesome, label: 'Accord Intelligence'),
      SizedBox(height: 16),

      // ── Top row: dominant profile + harmony score ──────
      Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Dominant profile card
          Expanded(flex: 3, child: _AccordDominantCard(result: result)),
          SizedBox(width: 12),
          // Harmony score card
          Expanded(
            flex: 2,
            child: _HarmonyScoreCard(
              score: result.harmonyScore,
              conflicts: result.conflicts,
              synergies: result.synergies,
            ),
          ),
        ],
      ),

      SizedBox(height: 12),

      // ── Accord scores bar chart ────────────────────────
      _AccordScoreBars(scores: result.scores),

      SizedBox(height: 12),

      // ── Conflicts & Synergies ──────────────────────────
      if (result.conflicts.isNotEmpty || result.synergies.isNotEmpty)
        _AccordInteractionBadges(
          conflicts: result.conflicts,
          synergies: result.synergies,
        ),
    ],
  );
}

// ── Accord Dominant Card ──────────────────────────────────
class _AccordDominantCard extends StatelessWidget {
  final AccordIntelligenceResult result;
  const _AccordDominantCard({required this.result});

  @override
  Widget build(BuildContext context) {
    final dominant = result.scores.isNotEmpty ? result.scores.first : null;

    return Container(
      padding: EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: _T.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: _T.green800.withOpacity(0.5)),
        gradient: LinearGradient(
          begin: Alignment.topLeft,
          end: Alignment.bottomRight,
          colors: [_T.bgCard, _T.green900.withOpacity(0.3)],
        ),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Row(
            children: [
              Icon(Icons.auto_awesome, color: _T.green400, size: 13),
              SizedBox(width: 6),
              Text(
                'DOMINANT PROFILE',
                style: TextStyle(
                  color: _T.green400,
                  fontSize: 10,
                  letterSpacing: 2,
                  fontWeight: FontWeight.w600,
                ),
              ),
            ],
          ),
          SizedBox(height: 12),
          if (dominant != null) ...[
            Text(dominant.accord.icon, style: TextStyle(fontSize: 32)),
            SizedBox(height: 6),
            Text(
              dominant.accord.name,
              style: TextStyle(
                color: _T.textPrimary,
                fontSize: 18,
                fontFamily: 'Georgia',
              ),
            ),
            SizedBox(height: 4),
            Text(
              dominant.accord.description,
              style: TextStyle(
                color: _T.textSecondary,
                fontSize: 11,
                height: 1.5,
              ),
            ),
          ],
          SizedBox(height: 12),
          Container(
            padding: EdgeInsets.all(10),
            decoration: BoxDecoration(
              color: Colors.black.withOpacity(0.2),
              borderRadius: BorderRadius.circular(3),
            ),
            child: Text(
              result.blendCharacter,
              style: TextStyle(
                color: _T.textSecondary,
                fontSize: 11,
                fontStyle: FontStyle.italic,
                height: 1.5,
              ),
            ),
          ),

          // Accord role hierarchy
          SizedBox(height: 12),
          Wrap(
            spacing: 6,
            runSpacing: 6,
            children:
                result.scores
                    .where((s) => s.role != AccordRole.absent)
                    .map((s) => _AccordRoleBadge(score: s))
                    .toList(),
          ),
        ],
      ),
    );
  }
}

// ── Harmony Score Card ────────────────────────────────────
class _HarmonyScoreCard extends StatelessWidget {
  final double score;
  final List<String> conflicts;
  final List<String> synergies;

  const _HarmonyScoreCard({
    required this.score,
    required this.conflicts,
    required this.synergies,
  });

  Color get _harmonyColor {
    if (score >= 0.80) return Color(0xFF43A047);
    if (score >= 0.60) return Color(0xFFFFA726);
    return Color(0xFFEF5350);
  }

  String get _harmonyLabel {
    if (score >= 0.80) return 'Harmonious';
    if (score >= 0.60) return 'Balanced';
    return 'Contrasting';
  }

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: _T.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: _T.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'HARMONY',
            style: TextStyle(
              color: _T.green400,
              fontSize: 10,
              letterSpacing: 2,
              fontWeight: FontWeight.w600,
            ),
          ),
          SizedBox(height: 16),

          // Circular harmony indicator
          Center(
            child: SizedBox(
              width: 90,
              height: 90,
              child: CustomPaint(
                painter: _HarmonyRingPainter(
                  score: score,
                  color: _harmonyColor,
                ),
                child: Center(
                  child: Column(
                    mainAxisAlignment: MainAxisAlignment.center,
                    children: [
                      Text(
                        '${(score * 100).toStringAsFixed(0)}',
                        style: TextStyle(
                          color: _harmonyColor,
                          fontSize: 22,
                          fontWeight: FontWeight.bold,
                          fontFamily: 'Courier',
                        ),
                      ),
                      Text(
                        '%',
                        style: TextStyle(color: _T.textMuted, fontSize: 10),
                      ),
                    ],
                  ),
                ),
              ),
            ),
          ),

          SizedBox(height: 10),
          Center(
            child: Text(
              _harmonyLabel,
              style: TextStyle(
                color: _harmonyColor,
                fontSize: 12,
                fontWeight: FontWeight.w600,
                letterSpacing: 1,
              ),
            ),
          ),

          SizedBox(height: 16),
          if (synergies.isNotEmpty) ...[
            _MiniLabel(
              icon: Icons.add_circle_outline,
              label:
                  '${synergies.length} synerg${synergies.length > 1 ? "ies" : "y"}',
              color: Color(0xFF43A047),
            ),
            SizedBox(height: 4),
          ],
          if (conflicts.isNotEmpty) ...[
            _MiniLabel(
              icon: Icons.warning_amber_outlined,
              label:
                  '${conflicts.length} conflict${conflicts.length > 1 ? "s" : ""}',
              color: Color(0xFFFFA726),
            ),
          ],
        ],
      ),
    );
  }
}

// ── Harmony Ring Painter ──────────────────────────────────
class _HarmonyRingPainter extends CustomPainter {
  final double score;
  final Color color;
  const _HarmonyRingPainter({required this.score, required this.color});

  @override
  void paint(Canvas canvas, Size size) {
    final center = Offset(size.width / 2, size.height / 2);
    final radius = size.width / 2 - 6;

    // Background ring
    canvas.drawCircle(
      center,
      radius,
      Paint()
        ..color = _T.border
        ..strokeWidth = 6
        ..style = PaintingStyle.stroke,
    );

    // Score arc
    final rect = Rect.fromCircle(center: center, radius: radius);
    canvas.drawArc(
      rect,
      -pi / 2,
      2 * pi * score,
      false,
      Paint()
        ..color = color
        ..strokeWidth = 6
        ..style = PaintingStyle.stroke
        ..strokeCap = StrokeCap.round,
    );
  }

  @override
  bool shouldRepaint(_HarmonyRingPainter old) => old.score != score;
}

// ── Accord Score Bars ─────────────────────────────────────
class _AccordScoreBars extends StatelessWidget {
  final List<AccordScore> scores;
  const _AccordScoreBars({required this.scores});

  @override
  Widget build(BuildContext context) {
    final visible = scores.where((s) => s.score > 0.05).toList();

    return Container(
      padding: EdgeInsets.all(16),
      decoration: BoxDecoration(
        color: _T.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: _T.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Text(
            'ACCORD PROFILE',
            style: TextStyle(
              color: _T.green400,
              fontSize: 10,
              letterSpacing: 2,
              fontWeight: FontWeight.w600,
            ),
          ),
          SizedBox(height: 12),
          ...visible.map(
            (s) => Padding(
              padding: EdgeInsets.only(bottom: 10),
              child: _AccordBar(score: s),
            ),
          ),
        ],
      ),
    );
  }
}

class _AccordBar extends StatelessWidget {
  final AccordScore score;
  const _AccordBar({required this.score});

  Color get _roleColor {
    switch (score.role) {
      case AccordRole.dominant:
        return Color(0xFF43A047);
      case AccordRole.supporting:
        return Color(0xFF29B6F6);
      case AccordRole.trace:
        return Color(0xFF8D6E63);
      default:
        return _T.border;
    }
  }

  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        // Icon + name
        SizedBox(
          width: 130,
          child: Row(
            children: [
              Text(score.accord.icon, style: TextStyle(fontSize: 14)),
              SizedBox(width: 6),
              Expanded(
                child: Text(
                  score.accord.name,
                  style: TextStyle(color: _T.textSecondary, fontSize: 11),
                  overflow: TextOverflow.ellipsis,
                ),
              ),
            ],
          ),
        ),
        SizedBox(width: 8),
        // Bar
        Expanded(
          child: Stack(
            children: [
              Container(
                height: 6,
                decoration: BoxDecoration(
                  color: _T.border,
                  borderRadius: BorderRadius.circular(3),
                ),
              ),
              FractionallySizedBox(
                widthFactor: score.score.clamp(0.0, 1.0),
                child: Container(
                  height: 6,
                  decoration: BoxDecoration(
                    color: _roleColor,
                    borderRadius: BorderRadius.circular(3),
                    boxShadow: [
                      BoxShadow(
                        color: _roleColor.withOpacity(0.4),
                        blurRadius: 4,
                      ),
                    ],
                  ),
                ),
              ),
            ],
          ),
        ),
        SizedBox(width: 8),
        // Score value + role badge
        SizedBox(
          width: 70,
          child: Row(
            mainAxisAlignment: MainAxisAlignment.end,
            children: [
              Text(
                '${(score.score * 100).toStringAsFixed(0)}%',
                style: TextStyle(
                  color: _T.textPrimary,
                  fontSize: 11,
                  fontFamily: 'Courier',
                ),
              ),
              SizedBox(width: 4),
              Container(
                padding: EdgeInsets.symmetric(horizontal: 4, vertical: 1),
                decoration: BoxDecoration(
                  color: _roleColor.withOpacity(0.15),
                  borderRadius: BorderRadius.circular(2),
                  border: Border.all(color: _roleColor.withOpacity(0.4)),
                ),
                child: Text(
                  _roleLabel(score.role),
                  style: TextStyle(
                    color: _roleColor,
                    fontSize: 8,
                    letterSpacing: 0.5,
                  ),
                ),
              ),
            ],
          ),
        ),
      ],
    );
  }

  String _roleLabel(AccordRole role) {
    switch (role) {
      case AccordRole.dominant:
        return 'DOM';
      case AccordRole.supporting:
        return 'SUP';
      case AccordRole.trace:
        return 'TRC';
      default:
        return '';
    }
  }
}

// ── Accord Interaction Badges ─────────────────────────────
class _AccordInteractionBadges extends StatelessWidget {
  final List<String> conflicts;
  final List<String> synergies;
  const _AccordInteractionBadges({
    required this.conflicts,
    required this.synergies,
  });

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.all(14),
      decoration: BoxDecoration(
        color: _T.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: _T.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          if (synergies.isNotEmpty) ...[
            _MiniLabel(
              icon: Icons.add_circle_outline,
              label: 'SYNERGIES',
              color: Color(0xFF43A047),
            ),
            SizedBox(height: 6),
            Wrap(
              spacing: 6,
              runSpacing: 6,
              children:
                  synergies
                      .map(
                        (s) => _InteractionChip(
                          label: s,
                          color: Color(0xFF43A047),
                        ),
                      )
                      .toList(),
            ),
            SizedBox(height: 12),
          ],
          if (conflicts.isNotEmpty) ...[
            _MiniLabel(
              icon: Icons.warning_amber_outlined,
              label: 'CONFLICTS',
              color: Color(0xFFFFA726),
            ),
            SizedBox(height: 6),
            Wrap(
              spacing: 6,
              runSpacing: 6,
              children:
                  conflicts
                      .map(
                        (s) => _InteractionChip(
                          label: s,
                          color: Color(0xFFFFA726),
                        ),
                      )
                      .toList(),
            ),
          ],
        ],
      ),
    );
  }
}

class _InteractionChip extends StatelessWidget {
  final String label;
  final Color color;
  const _InteractionChip({required this.label, required this.color});

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 8, vertical: 4),
      decoration: BoxDecoration(
        color: color.withOpacity(0.1),
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: color.withOpacity(0.3)),
      ),
      child: Text(label, style: TextStyle(color: color, fontSize: 10)),
    );
  }
}

// ── Accord Evolution Timeline ─────────────────────────────
Widget buildAccordEvolution(
  List<Compound> selectedList,
  List<double> optimizedFractions,
  List<double> Function(List<double>) gammaFn,
) {
  final List<TemporalAccordSnapshot> timeline = computeAccordEvolution(
    compounds: selectedList,
    x0: optimizedFractions,
    gammaFn: gammaFn,
    timePoints: [0, 0.25, 0.5, 1, 2, 4, 8],
  );

  return Column(
    crossAxisAlignment: CrossAxisAlignment.start,
    children: [
      SizedBox(height: 24),
      _SectionHeader(icon: Icons.timeline, label: 'Accord Evolution Timeline'),
      SizedBox(height: 4),
      Text(
        'Bagaimana karakter aroma berubah dari top note hingga drydown.',
        style: TextStyle(color: _T.textMuted, fontSize: 11),
      ),
      SizedBox(height: 16),

      Container(
        padding: EdgeInsets.all(16),
        decoration: BoxDecoration(
          color: _T.bgCard,
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: _T.border),
        ),
        child: Column(
          children: [
            // Timeline header labels
            Row(
              children: [
                SizedBox(width: 90),
                ...timeline.map(
                  (snap) => Expanded(
                    child: Text(
                      _timeLabel(snap.time),
                      textAlign: TextAlign.center,
                      style: TextStyle(
                        color: _T.textMuted,
                        fontSize: 9,
                        letterSpacing: 0.5,
                      ),
                    ),
                  ),
                ),
              ],
            ),
            SizedBox(height: 10),
            // Accord rows
            ...accordArchetypes
                .where((a) => _hasPresence(a, timeline))
                .map(
                  (accord) =>
                      _AccordTimelineRow(accord: accord, timeline: timeline),
                ),
          ],
        ),
      ),

      SizedBox(height: 12),

      // Time markers legend
      _TimelineLegend(),
    ],
  );
}

bool _hasPresence(AccordProfile accord, List<TemporalAccordSnapshot> timeline) {
  return timeline.any((snap) {
    final score =
        snap.accordResult.scores
            .firstWhere(
              (s) => s.accord.name == accord.name,
              orElse:
                  () => AccordScore(
                    accord: accord,
                    score: 0,
                    roiContrib: 0,
                    role: AccordRole.absent,
                  ),
            )
            .score;
    return score > 0.1;
  });
}

String _timeLabel(double t) {
  if (t < 1) return '${(t * 60).toStringAsFixed(0)}m';
  return '${t.toStringAsFixed(0)}h';
}

class _AccordTimelineRow extends StatelessWidget {
  final AccordProfile accord;
  final List<TemporalAccordSnapshot> timeline;

  const _AccordTimelineRow({required this.accord, required this.timeline});

  @override
  Widget build(BuildContext context) {
    return Padding(
      padding: EdgeInsets.only(bottom: 8),
      child: Row(
        children: [
          // Accord label
          SizedBox(
            width: 90,
            child: Row(
              children: [
                Text(accord.icon, style: TextStyle(fontSize: 12)),
                SizedBox(width: 4),
                Expanded(
                  child: Text(
                    accord.name,
                    style: TextStyle(color: _T.textSecondary, fontSize: 10),
                    overflow: TextOverflow.ellipsis,
                  ),
                ),
              ],
            ),
          ),
          // Score cells
          ...timeline.map((snap) {
            final AccordScore s = snap.accordResult.scores.firstWhere(
              (s) => s.accord.name == accord.name,
              orElse:
                  () => AccordScore(
                    accord: accord,
                    score: 0,
                    roiContrib: 0,
                    role: AccordRole.absent,
                  ),
            );
            return Expanded(child: _TimelineCell(score: s.score, role: s.role));
          }),
        ],
      ),
    );
  }
}

class _TimelineCell extends StatelessWidget {
  final double score;
  final AccordRole role;
  const _TimelineCell({required this.score, required this.role});

  Color get _cellColor {
    if (score < 0.1) return Colors.transparent;
    switch (role) {
      case AccordRole.dominant:
        return Color(0xFF43A047);
      case AccordRole.supporting:
        return Color(0xFF29B6F6);
      case AccordRole.trace:
        return Color(0xFF8D6E63);
      default:
        return Colors.transparent;
    }
  }

  @override
  Widget build(BuildContext context) {
    return Padding(
      padding: EdgeInsets.symmetric(horizontal: 2),
      child: Container(
        height: 22,
        decoration: BoxDecoration(
          color: _cellColor.withOpacity(score.clamp(0.0, 1.0) * 0.7),
          borderRadius: BorderRadius.circular(2),
          border:
              score > 0.1
                  ? Border.all(color: _cellColor.withOpacity(0.4))
                  : null,
        ),
        child:
            score > 0.1
                ? Center(
                  child: Text(
                    '${(score * 100).toStringAsFixed(0)}',
                    style: TextStyle(
                      color: Colors.white.withOpacity(0.8),
                      fontSize: 8,
                      fontFamily: 'Courier',
                    ),
                  ),
                )
                : null,
      ),
    );
  }
}

class _TimelineLegend extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        _legendDot(Color(0xFF43A047), 'Dominant'),
        SizedBox(width: 12),
        _legendDot(Color(0xFF29B6F6), 'Supporting'),
        SizedBox(width: 12),
        _legendDot(Color(0xFF8D6E63), 'Trace'),
      ],
    );
  }
}

// ── Climate Performance Widget ────────────────────────────
Widget buildClimatePerformance(List<ClimatePerformanceReport> reports) {
  if (reports.isEmpty) return SizedBox.shrink();

  return Column(
    crossAxisAlignment: CrossAxisAlignment.start,
    children: [
      SizedBox(height: 24),
      _SectionHeader(
        icon: Icons.thermostat_outlined,
        label: 'Climate Performance',
      ),
      SizedBox(height: 4),
      Text(
        'Proyeksi performa parfum di berbagai kondisi iklim.',
        style: TextStyle(color: _T.textMuted, fontSize: 11),
      ),
      SizedBox(height: 16),

      // Climate cards
      ...reports.map(
        (r) => Padding(
          padding: EdgeInsets.only(bottom: 10),
          child: _ClimateCard(report: r),
        ),
      ),
    ],
  );
}

class _ClimateCard extends StatelessWidget {
  final ClimatePerformanceReport report;
  const _ClimateCard({required this.report});

  @override
  Widget build(BuildContext context) {
    final c = report.climate;

    return Container(
      padding: EdgeInsets.all(14),
      decoration: BoxDecoration(
        color: _T.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: _T.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // Header
          Row(
            children: [
              Icon(_climateIcon(c.region), color: _T.green400, size: 14),
              SizedBox(width: 6),
              Text(
                c.name,
                style: TextStyle(
                  color: _T.textPrimary,
                  fontSize: 13,
                  fontWeight: FontWeight.w500,
                ),
              ),
              Spacer(),
              _ClimateBadge(
                label: '${c.ambientTempC.toStringAsFixed(0)}°C',
                icon: Icons.thermostat,
              ),
              SizedBox(width: 6),
              _ClimateBadge(
                label: '${(c.relativeHumidity * 100).toStringAsFixed(0)}% RH',
                icon: Icons.water_drop_outlined,
              ),
            ],
          ),

          SizedBox(height: 12),

          // Metrics row
          Row(
            children: [
              Expanded(
                child: _ClimateMetric(
                  label: 'Longevity',
                  value: _formatLongevity(report.projectedLongevityHours),
                  icon: Icons.access_time,
                  color: _longevityColor(report.projectedLongevityHours),
                ),
              ),
              Expanded(
                child: _ClimateMetric(
                  label: 'Projection',
                  value:
                      '${(report.projectionStrength * 100).toStringAsFixed(0)}%',
                  icon: Icons.spatial_audio_off,
                  color: _projectionColor(report.projectionStrength),
                ),
              ),
              Expanded(
                child: _ClimateMetric(
                  label: 'Drydown',
                  value:
                      report.drydownBalance < 1.0
                          ? 'Balanced'
                          : report.drydownBalance < 2.0
                          ? 'Moderate'
                          : 'Uneven',
                  icon: Icons.leaderboard_outlined,
                  color:
                      report.drydownBalance < 1.0
                          ? Color(0xFF43A047)
                          : report.drydownBalance < 2.0
                          ? Color(0xFFFFA726)
                          : Color(0xFFEF5350),
                ),
              ),
            ],
          ),

          SizedBox(height: 10),

          // Summary text
          Text(
            report.performanceSummary,
            style: TextStyle(
              color: _T.textMuted,
              fontSize: 10,
              height: 1.5,
              fontStyle: FontStyle.italic,
            ),
          ),

          // Longevity bar
          SizedBox(height: 8),
          _LongevityBar(hours: report.projectedLongevityHours),
        ],
      ),
    );
  }

  IconData _climateIcon(String region) {
    switch (region) {
      case 'Tropical':
        return Icons.wb_sunny;
      case 'Arid':
        return Icons.landscape;
      case 'Cold':
        return Icons.ac_unit;
      default:
        return Icons.cloud;
    }
  }

  String _formatLongevity(double h) {
    if (h >= 99) return '>24h';
    if (h >= 1) return '${h.toStringAsFixed(1)}h';
    return '${(h * 60).toStringAsFixed(0)}m';
  }

  Color _longevityColor(double h) {
    if (h >= 8) return Color(0xFF43A047);
    if (h >= 4) return Color(0xFF29B6F6);
    if (h >= 2) return Color(0xFFFFA726);
    return Color(0xFFEF5350);
  }

  Color _projectionColor(double p) {
    if (p >= 0.6) return Color(0xFF43A047);
    if (p >= 0.3) return Color(0xFF29B6F6);
    return Color(0xFF8D6E63);
  }
}

class _ClimateMetric extends StatelessWidget {
  final String label;
  final String value;
  final IconData icon;
  final Color color;
  const _ClimateMetric({
    required this.label,
    required this.value,
    required this.icon,
    required this.color,
  });

  @override
  Widget build(BuildContext context) {
    return Container(
      margin: EdgeInsets.symmetric(horizontal: 3),
      padding: EdgeInsets.symmetric(vertical: 8, horizontal: 6),
      decoration: BoxDecoration(
        color: color.withOpacity(0.08),
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: color.withOpacity(0.2)),
      ),
      child: Column(
        children: [
          Icon(icon, color: color, size: 14),
          SizedBox(height: 4),
          Text(
            value,
            style: TextStyle(
              color: color,
              fontSize: 12,
              fontWeight: FontWeight.w600,
              fontFamily: 'Courier',
            ),
          ),
          Text(label, style: TextStyle(color: _T.textMuted, fontSize: 9)),
        ],
      ),
    );
  }
}

class _LongevityBar extends StatelessWidget {
  final double hours;
  const _LongevityBar({required this.hours});

  @override
  Widget build(BuildContext context) {
    final double maxH = 12.0;
    final double fraction = (hours / maxH).clamp(0.0, 1.0);
    final Color barColor =
        hours >= 8
            ? Color(0xFF43A047)
            : hours >= 4
            ? Color(0xFF29B6F6)
            : Color(0xFFFFA726);

    return Row(
      children: [
        Text('0h', style: TextStyle(color: _T.textMuted, fontSize: 9)),
        SizedBox(width: 6),
        Expanded(
          child: Stack(
            children: [
              Container(
                height: 4,
                decoration: BoxDecoration(
                  color: _T.border,
                  borderRadius: BorderRadius.circular(2),
                ),
              ),
              FractionallySizedBox(
                widthFactor: fraction,
                child: Container(
                  height: 4,
                  decoration: BoxDecoration(
                    gradient: LinearGradient(
                      colors: [barColor.withOpacity(0.6), barColor],
                    ),
                    borderRadius: BorderRadius.circular(2),
                    boxShadow: [
                      BoxShadow(
                        color: barColor.withOpacity(0.4),
                        blurRadius: 4,
                      ),
                    ],
                  ),
                ),
              ),
            ],
          ),
        ),
        SizedBox(width: 6),
        Text('12h+', style: TextStyle(color: _T.textMuted, fontSize: 9)),
      ],
    );
  }
}

class _ClimateBadge extends StatelessWidget {
  final String label;
  final IconData icon;
  const _ClimateBadge({required this.label, required this.icon});

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 6, vertical: 2),
      decoration: BoxDecoration(
        color: _T.bg,
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: _T.border),
      ),
      child: Row(
        mainAxisSize: MainAxisSize.min,
        children: [
          Icon(icon, size: 10, color: _T.textMuted),
          SizedBox(width: 3),
          Text(label, style: TextStyle(color: _T.textSecondary, fontSize: 10)),
        ],
      ),
    );
  }
}

// ── Shared mini widgets ───────────────────────────────────
class _MiniLabel extends StatelessWidget {
  final IconData icon;
  final String label;
  final Color color;
  const _MiniLabel({
    required this.icon,
    required this.label,
    required this.color,
  });

  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        Icon(icon, size: 11, color: color),
        SizedBox(width: 4),
        Text(
          label,
          style: TextStyle(
            color: color,
            fontSize: 10,
            letterSpacing: 1,
            fontWeight: FontWeight.w600,
          ),
        ),
      ],
    );
  }
}

class _AccordRoleBadge extends StatelessWidget {
  final AccordScore score;
  const _AccordRoleBadge({required this.score});

  Color get _color {
    switch (score.role) {
      case AccordRole.dominant:
        return Color(0xFF43A047);
      case AccordRole.supporting:
        return Color(0xFF29B6F6);
      case AccordRole.trace:
        return Color(0xFF8D6E63);
      default:
        return _T.border;
    }
  }

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 7, vertical: 3),
      decoration: BoxDecoration(
        color: _color.withOpacity(0.12),
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: _color.withOpacity(0.35)),
      ),
      child: Text(
        '${score.accord.icon} ${score.accord.name}',
        style: TextStyle(color: _color, fontSize: 10),
      ),
    );
  }
}

// Mini widgets used in classifyNote
class _SectionHeader extends StatelessWidget {
  final IconData icon;
  final String label;
  const _SectionHeader({required this.icon, required this.label});

  @override
  Widget build(BuildContext context) {
    return Row(
      children: [
        Icon(icon, color: _T.green400, size: 16),
        SizedBox(width: 8),
        Text(
          label,
          style: TextStyle(
            fontFamily: 'Georgia',
            color: _T.textPrimary,
            fontSize: 18,
          ),
        ),
      ],
    );
  }
}

class _FractionChip extends StatelessWidget {
  final String label;
  final double value;
  final Color color;
  const _FractionChip({
    required this.label,
    required this.value,
    required this.color,
  });

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 14, vertical: 7),
      decoration: BoxDecoration(
        color: color.withOpacity(0.15),
        borderRadius: BorderRadius.circular(20),
        border: Border.all(color: color.withOpacity(0.4)),
      ),
      child: Row(
        children: [
          Container(
            width: 7,
            height: 7,
            decoration: BoxDecoration(shape: BoxShape.circle, color: color),
          ),
          SizedBox(width: 6),
          Text(
            '$label  ${(value * 100).toStringAsFixed(1)}%',
            style: TextStyle(color: _T.textPrimary, fontSize: 12),
          ),
        ],
      ),
    );
  }
}

class _DataBadge extends StatelessWidget {
  final String label;
  final String value;
  const _DataBadge({required this.label, required this.value});

  @override
  Widget build(BuildContext context) {
    return Container(
      padding: EdgeInsets.symmetric(horizontal: 6, vertical: 2),
      decoration: BoxDecoration(
        color: Colors.black.withOpacity(0.2),
        borderRadius: BorderRadius.circular(2),
      ),
      child: RichText(
        text: TextSpan(
          children: [
            TextSpan(
              text: '$label: ',
              style: TextStyle(color: _T.textMuted, fontSize: 10),
            ),
            TextSpan(
              text: value,
              style: TextStyle(
                color: _T.textPrimary,
                fontSize: 10,
                fontFamily: 'Courier',
              ),
            ),
          ],
        ),
      ),
    );
  }
}
