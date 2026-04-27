import 'dart:math';

import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/OctonaryPainter.dart';
import 'package:askreatif_app/TernaryPainter.dart';
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
