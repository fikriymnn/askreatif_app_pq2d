import 'dart:math';

import 'package:flutter/material.dart';

class OctonaryPainter extends CustomPainter {
  final List<String> labels;
  final List<double> values; // normalized 0..1
  final List<Color> colors;
  final List<double> rawValues; // nilai asli untuk label

  OctonaryPainter({
    required this.labels,
    required this.values,
    required this.colors,
    required this.rawValues,
  });

  @override
  void paint(Canvas canvas, Size size) {
    final int n = labels.length; // 8
    final Offset center = Offset(size.width / 2, size.height / 2);
    final double maxRadius = min(size.width, size.height) / 2 - 50;

    // === Gambar grid polygon (20%, 40%, 60%, 80%, 100%) ===
    final gridPaint =
        Paint()
          ..color = Colors.grey[300]!
          ..strokeWidth = 0.8
          ..style = PaintingStyle.stroke;

    for (int ring = 1; ring <= 5; ring++) {
      double r = maxRadius * ring / 5;
      Path gridPath = Path();
      for (int i = 0; i < n; i++) {
        double angle = 2 * pi * i / n - pi / 2;
        double x = center.dx + r * cos(angle);
        double y = center.dy + r * sin(angle);
        if (i == 0)
          gridPath.moveTo(x, y);
        else
          gridPath.lineTo(x, y);
      }
      gridPath.close();
      canvas.drawPath(gridPath, gridPaint);

      // Label % di ring
      double labelAngle = -pi / 2;
      double lx = center.dx + r * cos(labelAngle) + 4;
      double ly = center.dy + r * sin(labelAngle) - 8;
      final pct = TextPainter(
        text: TextSpan(
          text: "${ring * 20}%",
          style: TextStyle(fontSize: 8, color: Colors.grey[500]),
        ),
        textDirection: TextDirection.ltr,
      )..layout();
      pct.paint(canvas, Offset(lx, ly));
    }

    // === Gambar sumbu dari center ke tiap sudut ===
    final axisPaint =
        Paint()
          ..color = Colors.grey[400]!
          ..strokeWidth = 0.8;

    for (int i = 0; i < n; i++) {
      double angle = 2 * pi * i / n - pi / 2;
      double x = center.dx + maxRadius * cos(angle);
      double y = center.dy + maxRadius * sin(angle);
      canvas.drawLine(center, Offset(x, y), axisPaint);
    }

    // === Gambar area per senyawa (tiap segmen diberi warna berbeda) ===
    for (int i = 0; i < n; i++) {
      int next = (i + 1) % n;
      double angle1 = 2 * pi * i / n - pi / 2;
      double angle2 = 2 * pi * next / n - pi / 2;

      double r1 = maxRadius * values[i];
      double r2 = maxRadius * values[next];

      Path segPath = Path();
      segPath.moveTo(center.dx, center.dy);
      segPath.lineTo(
        center.dx + r1 * cos(angle1),
        center.dy + r1 * sin(angle1),
      );

      // Arc kecil antara dua titik
      int arcSteps = 10;
      for (int s = 0; s <= arcSteps; s++) {
        double t = s / arcSteps;
        double angle = angle1 + t * (angle2 - angle1);
        double r = r1 + t * (r2 - r1);
        segPath.lineTo(center.dx + r * cos(angle), center.dy + r * sin(angle));
      }
      segPath.close();

      // Fill dengan warna note
      final fillPaint =
          Paint()
            ..color = colors[i].withOpacity(0.35)
            ..style = PaintingStyle.fill;
      canvas.drawPath(segPath, fillPaint);

      // Border segmen
      final strokePaint =
          Paint()
            ..color = colors[i].withOpacity(0.8)
            ..strokeWidth = 1.5
            ..style = PaintingStyle.stroke;
      canvas.drawPath(segPath, strokePaint);
    }

    // === Gambar titik nilai di tiap sumbu ===
    for (int i = 0; i < n; i++) {
      double angle = 2 * pi * i / n - pi / 2;
      double r = maxRadius * values[i];
      double px = center.dx + r * cos(angle);
      double py = center.dy + r * sin(angle);

      final dotPaint = Paint()..color = colors[i];
      canvas.drawCircle(Offset(px, py), 5, dotPaint);
      final dotBorder =
          Paint()
            ..color = Colors.white
            ..strokeWidth = 1.5
            ..style = PaintingStyle.stroke;
      canvas.drawCircle(Offset(px, py), 5, dotBorder);
    }

    // === Label nama senyawa + nilai OV di tiap sudut ===
    for (int i = 0; i < n; i++) {
      double angle = 2 * pi * i / n - pi / 2;
      double labelR = maxRadius + 18;
      double lx = center.dx + labelR * cos(angle);
      double ly = center.dy + labelR * sin(angle);

      // Format OV: jika sangat besar gunakan scientific
      String ovText;
      if (rawValues[i] >= 1000) {
        ovText = "${rawValues[i].toStringAsExponential(2)}";
      } else {
        ovText = rawValues[i].toStringAsExponential(2);
      }

      // Nama senyawa
      final namePainter = TextPainter(
        text: TextSpan(
          text: labels[i],
          style: TextStyle(
            fontSize: 10,
            fontWeight: FontWeight.bold,
            color: colors[i].withOpacity(0.9),
          ),
        ),
        textDirection: TextDirection.ltr,
        textAlign: TextAlign.center,
      )..layout(maxWidth: 80);

      // Nilai OV
      final ovPainter = TextPainter(
        text: TextSpan(
          text: "OV: $ovText",
          style: TextStyle(fontSize: 9, color: Colors.grey[600]),
        ),
        textDirection: TextDirection.ltr,
        textAlign: TextAlign.center,
      )..layout(maxWidth: 80);

      // Posisi label berdasarkan sudut
      double offsetX = -namePainter.width / 2;
      double offsetY = -namePainter.height / 2;

      // Geser label agar tidak tumpang tindih dengan diagram
      if (cos(angle) > 0.3)
        offsetX = 4;
      else if (cos(angle) < -0.3)
        offsetX = -namePainter.width - 4;

      if (sin(angle) > 0.3)
        offsetY = 4;
      else if (sin(angle) < -0.3)
        offsetY = -namePainter.height - ovPainter.height - 2;

      namePainter.paint(canvas, Offset(lx + offsetX, ly + offsetY));
      ovPainter.paint(
        canvas,
        Offset(lx + offsetX, ly + offsetY + namePainter.height + 1),
      );
    }
  }

  @override
  bool shouldRepaint(OctonaryPainter old) => true;
}
