import 'package:flutter/material.dart';

class TernaryPainter extends CustomPainter {
  final double top;
  final double middle;
  final double base;

  TernaryPainter({required this.top, required this.middle, required this.base});

  // Konversi koordinat ternary → Cartesian
  // Sumbu: Top = puncak, Middle = kiri bawah, Base = kanan bawah
  Offset ternaryToCart(double t, double m, double b, Size size) {
    final double margin = 40.0;
    final double w = size.width - 2 * margin;
    final double h = size.height - 2 * margin;

    // Koordinat puncak segitiga sama sisi
    final Offset apex = Offset(margin + w / 2, margin);
    final Offset left = Offset(margin, margin + h);
    final Offset right = Offset(margin + w, margin + h);

    return Offset(
      t * apex.dx + m * left.dx + b * right.dx,
      t * apex.dy + m * left.dy + b * right.dy,
    );
  }

  @override
  void paint(Canvas canvas, Size size) {
    final margin = 40.0;
    final w = size.width - 2 * margin;
    final h = size.height - 2 * margin;

    final Offset apex = Offset(margin + w / 2, margin);
    final Offset left = Offset(margin, margin + h);
    final Offset right = Offset(margin + w, margin + h);

    // === Gambar segitiga ===
    final borderPaint =
        Paint()
          ..color = Colors.grey[400]!
          ..strokeWidth = 1.5
          ..style = PaintingStyle.stroke;
    final path =
        Path()
          ..moveTo(apex.dx, apex.dy)
          ..lineTo(left.dx, left.dy)
          ..lineTo(right.dx, right.dy)
          ..close();
    canvas.drawPath(path, borderPaint);

    // === Grid lines (20% intervals) ===
    final gridPaint =
        Paint()
          ..color = Colors.grey[300]!
          ..strokeWidth = 0.5
          ..style = PaintingStyle.stroke;

    for (int i = 1; i < 5; i++) {
      double f = i / 5.0;

      // Grid paralel tiap sumbu
      Offset a1 = ternaryToCart(f, 1 - f, 0, size);
      Offset a2 = ternaryToCart(f, 0, 1 - f, size);
      canvas.drawLine(a1, a2, gridPaint);

      Offset b1 = ternaryToCart(0, f, 1 - f, size);
      Offset b2 = ternaryToCart(1 - f, f, 0, size);
      canvas.drawLine(b1, b2, gridPaint);

      Offset c1 = ternaryToCart(1 - f, 0, f, size);
      Offset c2 = ternaryToCart(0, 1 - f, f, size);
      canvas.drawLine(c1, c2, gridPaint);
    }

    // === Label sudut ===
    final labelStyle = TextStyle(fontSize: 13, fontWeight: FontWeight.bold);

    void drawLabel(String text, Offset pos, {double dx = 0, double dy = 0}) {
      final tp = TextPainter(
        text: TextSpan(text: text, style: labelStyle),
        textDirection: TextDirection.ltr,
      )..layout();
      tp.paint(canvas, Offset(pos.dx - tp.width / 2 + dx, pos.dy + dy));
    }

    drawLabel("Top", apex, dy: -18);
    drawLabel("Middle", left, dx: -10, dy: 4);
    drawLabel("Base", right, dx: 10, dy: 4);

    // === Tick labels 20% ===
    final tickStyle = TextStyle(fontSize: 9, color: Colors.grey[600]);
    for (int i = 1; i < 5; i++) {
      double f = i / 5.0;
      String pct = "${(f * 100).round()}%";

      // Label kiri (Top axis)
      Offset tPos = ternaryToCart(f, 1 - f, 0, size);
      final tp1 = TextPainter(
        text: TextSpan(text: pct, style: tickStyle),
        textDirection: TextDirection.ltr,
      )..layout();
      tp1.paint(canvas, Offset(tPos.dx - tp1.width - 4, tPos.dy - 6));

      // Label bawah (Base axis)
      Offset bPos = ternaryToCart(0, f, 1 - f, size);
      final tp2 = TextPainter(
        text: TextSpan(text: pct, style: tickStyle),
        textDirection: TextDirection.ltr,
      )..layout();
      tp2.paint(canvas, Offset(bPos.dx - tp2.width / 2, bPos.dy + 4));
    }

    // === Sudut warna ===
    final topPaint = Paint()..color = Colors.green[300]!;
    final midPaint = Paint()..color = Colors.pink[300]!;
    final basePaint = Paint()..color = Colors.brown[300]!;

    canvas.drawCircle(apex, 8, topPaint);
    canvas.drawCircle(left, 8, midPaint);
    canvas.drawCircle(right, 8, basePaint);

    // === Plot titik komposisi ===
    final Offset point = ternaryToCart(top, middle, base, size);

    // Shadow
    final shadowPaint =
        Paint()
          ..color = Colors.black26
          ..maskFilter = MaskFilter.blur(BlurStyle.normal, 4);
    canvas.drawCircle(point + Offset(2, 2), 10, shadowPaint);

    // Titik utama
    final pointPaint = Paint()..color = Colors.purple[700]!;
    canvas.drawCircle(point, 10, pointPaint);

    // Border putih
    final borderWhite =
        Paint()
          ..color = Colors.white
          ..strokeWidth = 2
          ..style = PaintingStyle.stroke;
    canvas.drawCircle(point, 10, borderWhite);

    // Label koordinat di titik
    final coordStyle = TextStyle(
      fontSize: 10,
      color: Colors.purple[900],
      fontWeight: FontWeight.bold,
    );
    final coordText =
        "T:${(top * 100).toStringAsExponential(2)}% "
        "M:${(middle * 100).toStringAsExponential(2)}% "
        "B:${(base * 100).toStringAsExponential(2)}%";
    final tp = TextPainter(
      text: TextSpan(text: coordText, style: coordStyle),
      textDirection: TextDirection.ltr,
    )..layout();

    // Posisi label agar tidak keluar canvas
    double labelX = (point.dx + 14).clamp(0, size.width - tp.width);
    double labelY = (point.dy - 8).clamp(0, size.height - tp.height);
    tp.paint(canvas, Offset(labelX, labelY));
  }

  @override
  bool shouldRepaint(TernaryPainter old) =>
      old.top != top || old.middle != middle || old.base != base;
}
