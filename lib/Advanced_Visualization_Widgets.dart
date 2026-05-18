// ============================================================
// ASKREATIF Perfumery Engine — Phase 10
// Advanced Visualization Widgets
// ============================================================
// • EvolutionChartWidget    — fragrance evolution over time
// • DynamicNotesBadges      — time-aware note labels
// • AccordRadarWidget        — accord balance radar chart
// • SafetyAlertsPanel        — IFRA/EU alerts display
// • HeadspaceBarWidget       — headspace concentration bars
// ============================================================

import 'dart:math';
import 'package:askreatif_app/AccordIntelligence.dart';
import 'package:askreatif_app/Evaporation.dart';
import 'package:askreatif_app/Roi.dart';
import 'package:flutter/material.dart';
import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/colorTheme.dart';


import 'ifra_safety.dart';

// ══════════════════════════════════════════════════════════
// 1. FRAGRANCE EVOLUTION CHART
// ══════════════════════════════════════════════════════════
/// Interactive time-slider chart showing OV decay curves.
class EvolutionChartWidget extends StatefulWidget {
  final List<EvaporationPoint> trajectory;
  final List<Compound> compounds;

  const EvolutionChartWidget({
    Key? key,
    required this.trajectory,
    required this.compounds,
  }) : super(key: key);

  @override
  State<EvolutionChartWidget> createState() => _EvolutionChartWidgetState();
}

class _EvolutionChartWidgetState extends State<EvolutionChartWidget> {
  double _sliderT = 0.0; // hours

  // Colour palette per compound index
  static const List<Color> _palette = [
    Color(0xFF4FC3F7), // cyan
    Color(0xFFA5D6A7), // green
    Color(0xFFCE93D8), // purple
    Color(0xFFFFCC80), // amber
    Color(0xFFEF9A9A), // red
    Color(0xFF80DEEA), // teal
  ];

  @override
  Widget build(BuildContext context) {
    final maxT = widget.trajectory.last.time;

    // Current snapshot at slider position
    final EvaporationPoint snap = _snapshotAt(_sliderT);

    return Container(
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: AppTheme.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          // ── Header ─────────────────────────────────────
          _cardHeader(
            Icons.timeline,
            'FRAGRANCE EVOLUTION  —  t = ${_sliderT.toStringAsFixed(1)} h',
          ),

          // ── Chart canvas ───────────────────────────────
          SizedBox(
            height: 200,
            child: CustomPaint(
              painter: _EvolutionCurvePainter(
                trajectory: widget.trajectory,
                compounds: widget.compounds,
                palette: _palette,
                currentT: _sliderT,
              ),
            ),
          ),

          // ── Time slider ────────────────────────────────
          Padding(
            padding: const EdgeInsets.symmetric(horizontal: 16),
            child: Column(
              children: [
                SliderTheme(
                  data: SliderThemeData(
                    trackHeight: 2,
                    activeTrackColor: AppTheme.green400,
                    inactiveTrackColor: AppTheme.border,
                    thumbColor: AppTheme.green400,
                    overlayColor: AppTheme.green400.withOpacity(0.12),
                    thumbShape:
                        const RoundSliderThumbShape(enabledThumbRadius: 6),
                  ),
                  child: Slider(
                    value: _sliderT,
                    min: 0,
                    max: maxT,
                    onChanged: (v) => setState(() => _sliderT = v),
                  ),
                ),
                Row(
                  mainAxisAlignment: MainAxisAlignment.spaceBetween,
                  children: [
                    Text('0 h',
                        style: TextStyle(
                            color: AppTheme.textMuted, fontSize: 10)),
                    Text('${maxT.toStringAsFixed(0)} h',
                        style: TextStyle(
                            color: AppTheme.textMuted, fontSize: 10)),
                  ],
                ),
              ],
            ),
          ),

          const SizedBox(height: 8),

          // ── Snapshot bar: OV at current t ──────────────
          Padding(
            padding:
                const EdgeInsets.symmetric(horizontal: 16, vertical: 8),
            child: _OvSnapshotBars(
              snap: snap,
              compounds: widget.compounds,
              palette: _palette,
            ),
          ),

          // ── Legend ─────────────────────────────────────
          Padding(
            padding:
                const EdgeInsets.fromLTRB(16, 0, 16, 16),
            child: Wrap(
              spacing: 12,
              runSpacing: 4,
              children: [
                for (int i = 0;
                    i < widget.compounds.length;
                    i++)
                  _LegendDot(
                    color: _palette[i % _palette.length],
                    label: widget.compounds[i].name,
                  ),
              ],
            ),
          ),
        ],
      ),
    );
  }

  EvaporationPoint _snapshotAt(double t) {
    if (widget.trajectory.isEmpty) {
      return EvaporationPoint(
        time: t,
        moleFractions: [],
        headspaceConc: [],
        odorValues: [],
        noteLabels: [],
      );
    }
    EvaporationPoint best = widget.trajectory.first;
    for (final p in widget.trajectory) {
      if ((p.time - t).abs() < (best.time - t).abs()) best = p;
    }
    return best;
  }
}

// ── Custom painter for OV curves ──────────────────────────
class _EvolutionCurvePainter extends CustomPainter {
  final List<EvaporationPoint> trajectory;
  final List<Compound> compounds;
  final List<Color> palette;
  final double currentT;

  _EvolutionCurvePainter({
    required this.trajectory,
    required this.compounds,
    required this.palette,
    required this.currentT,
  });

  @override
  void paint(Canvas canvas, Size size) {
    if (trajectory.isEmpty || compounds.isEmpty) return;

    final double maxT = trajectory.last.time;
    final int n = compounds.length;

    // Compute max OV across all time/compounds (for y scaling)
    double maxOV = 1e-10;
    for (final p in trajectory) {
      for (final ov in p.odorValues) {
        if (ov > maxOV) maxOV = ov;
      }
    }

    double tx(double t) => t / maxT * size.width;
    double ty(double ov) =>
        size.height - (log(max(ov, 1e-20)) / log(max(maxOV, 1e-10))) * size.height * 0.9;

    // Grid lines
    final Paint gridPaint = Paint()
      ..color = AppTheme.border.withOpacity(0.4)
      ..strokeWidth = 0.5;
    for (int h in [2, 4, 8, 12, 24]) {
      final double x = tx(h.toDouble());
      if (x > size.width) break;
      canvas.drawLine(Offset(x, 0), Offset(x, size.height), gridPaint);
    }

    // Curves
    for (int i = 0; i < n; i++) {
      final Color c = palette[i % palette.length];
      final Paint paint = Paint()
        ..color = c
        ..strokeWidth = 1.5
        ..style = PaintingStyle.stroke;

      final Path path = Path();
      bool started = false;
      for (final p in trajectory) {
        if (i >= p.odorValues.length) continue;
        final double x = tx(p.time);
        final double y = ty(p.odorValues[i]);
        if (!started) {
          path.moveTo(x, y);
          started = true;
        } else {
          path.lineTo(x, y);
        }
      }
      canvas.drawPath(path, paint);
    }

    // Current-time vertical line
    final Paint tLine = Paint()
      ..color = AppTheme.green400.withOpacity(0.7)
      ..strokeWidth = 1.0
      ..style = PaintingStyle.stroke;
    final double cx = tx(currentT);
    canvas.drawLine(Offset(cx, 0), Offset(cx, size.height), tLine);
  }

  @override
  bool shouldRepaint(covariant _EvolutionCurvePainter old) =>
      old.currentT != currentT || old.trajectory != trajectory;
}

// ── OV snapshot horizontal bars ───────────────────────────
class _OvSnapshotBars extends StatelessWidget {
  final EvaporationPoint snap;
  final List<Compound> compounds;
  final List<Color> palette;

  const _OvSnapshotBars({
    required this.snap,
    required this.compounds,
    required this.palette,
  });

  @override
  Widget build(BuildContext context) {
    if (snap.odorValues.isEmpty) return const SizedBox.shrink();
    final double maxOV =
        snap.odorValues.fold(0.0, (a, b) => a > b ? a : b);
    if (maxOV <= 0) return const SizedBox.shrink();

    return Column(
      children: [
        for (int i = 0; i < compounds.length; i++)
          Padding(
            padding: const EdgeInsets.symmetric(vertical: 2),
            child: Row(
              children: [
                SizedBox(
                  width: 110,
                  child: Text(
                    compounds[i].name,
                    style: TextStyle(
                        color: AppTheme.textSecondary, fontSize: 10),
                    overflow: TextOverflow.ellipsis,
                  ),
                ),
                Expanded(
                  child: Stack(
                    children: [
                      Container(
                        height: 8,
                        decoration: BoxDecoration(
                          color: AppTheme.bg,
                          borderRadius: BorderRadius.circular(2),
                        ),
                      ),
                      FractionallySizedBox(
                        widthFactor: (i < snap.odorValues.length
                                ? snap.odorValues[i] / maxOV
                                : 0.0)
                            .clamp(0.0, 1.0),
                        child: Container(
                          height: 8,
                          decoration: BoxDecoration(
                            color: palette[i % palette.length]
                                .withOpacity(0.75),
                            borderRadius: BorderRadius.circular(2),
                          ),
                        ),
                      ),
                    ],
                  ),
                ),
                SizedBox(
                  width: 54,
                  child: Text(
                    i < snap.odorValues.length
                        ? snap.odorValues[i].toStringAsFixed(2)
                        : '-',
                    style: TextStyle(
                        color: palette[i % palette.length],
                        fontSize: 10,
                        fontFamily: 'Courier'),
                    textAlign: TextAlign.right,
                  ),
                ),
              ],
            ),
          ),
      ],
    );
  }
}

// ══════════════════════════════════════════════════════════
// 2. ACCORD RADAR WIDGET
// ══════════════════════════════════════════════════════════
class AccordRadarWidget extends StatelessWidget {
  final List<AccordScore> scores;

  const AccordRadarWidget({Key? key, required this.scores}) : super(key: key);

  @override
  Widget build(BuildContext context) {
    final top = scores.take(6).toList();

    return Container(
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: AppTheme.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          _cardHeader(Icons.radar, 'ACCORD PROFILE'),
          Padding(
            padding: const EdgeInsets.fromLTRB(16, 0, 16, 16),
            child: Column(
              children: [
                SizedBox(
                  height: 180,
                  child: CustomPaint(
                    painter: _RadarPainter(scores: top),
                    size: const Size(double.infinity, 180),
                  ),
                ),
                const SizedBox(height: 12),
                Wrap(
                  spacing: 10,
                  runSpacing: 6,
                  children: top.map((s) {
                    return Row(
                      mainAxisSize: MainAxisSize.min,
                      children: [
                        Text(s.accord.icon, style: const TextStyle(fontSize: 14)),
                        const SizedBox(width: 4),
                        Text(s.accord.name,
                            style: TextStyle(
                                color: AppTheme.textSecondary,
                                fontSize: 11)),
                        const SizedBox(width: 4),
                        Text(
                            '${(s.score * 100).toStringAsFixed(0)}%',
                            style: TextStyle(
                                color: AppTheme.green400,
                                fontSize: 11,
                                fontFamily: 'Courier')),
                      ],
                    );
                  }).toList(),
                ),
              ],
            ),
          ),
        ],
      ),
    );
  }
}

class _RadarPainter extends CustomPainter {
  final List<AccordScore> scores;
  _RadarPainter({required this.scores});

  @override
  void paint(Canvas canvas, Size size) {
    if (scores.isEmpty) return;
    final int n = scores.length;
    final Offset center = Offset(size.width / 2, size.height / 2);
    final double radius = min(size.width, size.height) / 2 - 24;
    final double angleStep = 2 * pi / n;

    // Grid rings
    final Paint gridPaint = Paint()
      ..color = AppTheme.border
      ..style = PaintingStyle.stroke
      ..strokeWidth = 0.5;
    for (int r = 1; r <= 4; r++) {
      final Path ring = Path();
      for (int k = 0; k < n; k++) {
        final double angle = -pi / 2 + k * angleStep;
        final double rr = radius * r / 4;
        final Offset pt = center + Offset(rr * cos(angle), rr * sin(angle));
        k == 0 ? ring.moveTo(pt.dx, pt.dy) : ring.lineTo(pt.dx, pt.dy);
      }
      ring.close();
      canvas.drawPath(ring, gridPaint);
    }

    // Axes
    for (int k = 0; k < n; k++) {
      final double angle = -pi / 2 + k * angleStep;
      canvas.drawLine(
        center,
        center + Offset(radius * cos(angle), radius * sin(angle)),
        gridPaint,
      );
    }

    // Data polygon
    final Paint fillPaint = Paint()
      ..color = AppTheme.green400.withOpacity(0.15)
      ..style = PaintingStyle.fill;
    final Paint strokePaint = Paint()
      ..color = AppTheme.green400
      ..style = PaintingStyle.stroke
      ..strokeWidth = 1.5;

    final Path dataPath = Path();
    for (int k = 0; k < n; k++) {
      final double angle = -pi / 2 + k * angleStep;
      final double rr = radius * scores[k].score;
      final Offset pt = center + Offset(rr * cos(angle), rr * sin(angle));
      k == 0 ? dataPath.moveTo(pt.dx, pt.dy) : dataPath.lineTo(pt.dx, pt.dy);
    }
    dataPath.close();
    canvas.drawPath(dataPath, fillPaint);
    canvas.drawPath(dataPath, strokePaint);

    // Labels
    for (int k = 0; k < n; k++) {
      final double angle = -pi / 2 + k * angleStep;
      final Offset labelPos = center +
          Offset((radius + 18) * cos(angle), (radius + 18) * sin(angle));
      final TextPainter tp = TextPainter(
        text: TextSpan(
          text: scores[k].accord.name.split(' ').first,
          style: TextStyle(
              color: AppTheme.textSecondary, fontSize: 9),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      tp.paint(canvas, labelPos - Offset(tp.width / 2, tp.height / 2));
    }
  }

  @override
  bool shouldRepaint(covariant _RadarPainter old) => old.scores != scores;
}

// ══════════════════════════════════════════════════════════
// 3. SAFETY ALERTS PANEL
// ══════════════════════════════════════════════════════════
class SafetyAlertsPanel extends StatelessWidget {
  final List<SafetyAlert> alerts;

  const SafetyAlertsPanel({Key? key, required this.alerts}) : super(key: key);

  @override
  Widget build(BuildContext context) {
    if (alerts.isEmpty) {
      return Container(
        padding: const EdgeInsets.all(14),
        decoration: BoxDecoration(
          color: const Color(0xFF1B2B1E),
          borderRadius: BorderRadius.circular(4),
          border: Border.all(color: AppTheme.green800),
        ),
        child: Row(
          children: [
            const Icon(Icons.check_circle_outline,
                color: Color(0xFF81C784), size: 16),
            const SizedBox(width: 8),
            Text(
              'No IFRA / EU safety issues detected for this formula.',
              style: TextStyle(
                  color: AppTheme.green400,
                  fontSize: 12),
            ),
          ],
        ),
      );
    }

    return Container(
      decoration: BoxDecoration(
        color: AppTheme.bgCard,
        borderRadius: BorderRadius.circular(4),
        border: Border.all(color: AppTheme.border),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          _cardHeader(Icons.shield_outlined, 'SAFETY & REGULATORY ALERTS'),
          for (final alert in alerts) _alertTile(alert),
          const SizedBox(height: 4),
          Padding(
            padding: const EdgeInsets.fromLTRB(16, 0, 16, 12),
            child: Text(
              '⚠️  Data is for research only. Verify against IFRA 51st Amendment '
              'and consult a regulatory toxicologist before commercial use.',
              style: TextStyle(
                  color: AppTheme.textMuted, fontSize: 9.5, height: 1.5),
            ),
          ),
        ],
      ),
    );
  }

  Widget _alertTile(SafetyAlert alert) {
    Color bg, border, icon;
    IconData ico;
    switch (alert.severity) {
      case AlertSeverity.critical:
        bg = const Color(0xFF3E1C1C);
        border = const Color(0xFFEF9A9A);
        icon = const Color(0xFFEF9A9A);
        ico = Icons.error_outline;
        break;
      case AlertSeverity.warning:
        bg = const Color(0xFF3E2E0A);
        border = const Color(0xFFFFCC80);
        icon = const Color(0xFFFFCC80);
        ico = Icons.warning_amber_outlined;
        break;
      case AlertSeverity.info:
        bg = const Color(0xFF0D2B3A);
        border = const Color(0xFF4FC3F7);
        icon = const Color(0xFF4FC3F7);
        ico = Icons.info_outline;
        break;
    }

    return Container(
      margin: const EdgeInsets.fromLTRB(12, 6, 12, 0),
      padding: const EdgeInsets.all(10),
      decoration: BoxDecoration(
        color: bg,
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: border.withOpacity(0.4)),
      ),
      child: Row(
        crossAxisAlignment: CrossAxisAlignment.start,
        children: [
          Icon(ico, color: icon, size: 14),
          const SizedBox(width: 8),
          Expanded(
            child: Text(
              alert.message,
              style: TextStyle(color: icon, fontSize: 11, height: 1.4),
            ),
          ),
        ],
      ),
    );
  }
}

// ══════════════════════════════════════════════════════════
// 4. DYNAMIC NOTE BADGES
// ══════════════════════════════════════════════════════════
class DynamicNoteBadgesWidget extends StatelessWidget {
  final List<NoteState> notes;
  final List<Compound> compounds;

  const DynamicNoteBadgesWidget({
    Key? key,
    required this.notes,
    required this.compounds,
  }) : super(key: key);

  @override
  Widget build(BuildContext context) {
    return Wrap(
      spacing: 8,
      runSpacing: 6,
      children: [
        for (int i = 0; i < compounds.length; i++)
          _NoteBadge(
            compoundName: compounds[i].name,
            note: i < notes.length ? notes[i] : null,
          ),
      ],
    );
  }
}

class _NoteBadge extends StatelessWidget {
  final String compoundName;
  final NoteState? note;

  const _NoteBadge({required this.compoundName, this.note});

  @override
  Widget build(BuildContext context) {
    final Color col = note == null
        ? AppTheme.textMuted
        : Color(int.parse(
            note!.colorHex.replaceFirst('#', 'FF'), radix: 16));
    return Container(
      padding: const EdgeInsets.symmetric(horizontal: 8, vertical: 4),
      decoration: BoxDecoration(
        color: col.withOpacity(0.1),
        borderRadius: BorderRadius.circular(3),
        border: Border.all(color: col.withOpacity(0.4)),
      ),
      child: Column(
        crossAxisAlignment: CrossAxisAlignment.start,
        mainAxisSize: MainAxisSize.min,
        children: [
          Text(compoundName,
              style: TextStyle(
                  color: AppTheme.textPrimary, fontSize: 10)),
          if (note != null)
            Text(note!.fullLabel,
                style: TextStyle(color: col, fontSize: 9)),
        ],
      ),
    );
  }
}

// ── Shared helpers ─────────────────────────────────────────
Widget _cardHeader(IconData icon, String title) {
  return Padding(
    padding: const EdgeInsets.fromLTRB(16, 14, 16, 10),
    child: Row(
      children: [
        Icon(icon, color: AppTheme.green400, size: 14),
        const SizedBox(width: 6),
        Text(
          title,
          style: TextStyle(
            color: AppTheme.green400,
            fontSize: 10,
            letterSpacing: 1.5,
            fontWeight: FontWeight.w600,
          ),
        ),
      ],
    ),
  );
}

class _LegendDot extends StatelessWidget {
  final Color color;
  final String label;
  const _LegendDot({required this.color, required this.label});

  @override
  Widget build(BuildContext context) {
    return Row(
      mainAxisSize: MainAxisSize.min,
      children: [
        Container(
          width: 8,
          height: 8,
          decoration: BoxDecoration(shape: BoxShape.circle, color: color),
        ),
        const SizedBox(width: 4),
        Text(label,
            style: TextStyle(
                color: AppTheme.textSecondary, fontSize: 10)),
      ],
    );
  }
}