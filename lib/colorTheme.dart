import 'package:flutter/material.dart';

class AppTheme {
  // Core colors
  static const Color bg = Color(0xFF0D2818);
  static const Color bgSurface = Color(0xFF122B1D);
  static const Color bgCard = Color(0xFF173322);
  static const Color bgCardHover = Color(0xFF1C3D29);

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
  static const Color borderLight = Color(0xFF2E6038);

  // Typography
  static const String fontSerif = 'Georgia';

  // Decorative gradient
  static LinearGradient get primaryGradient => LinearGradient(
    colors: [green900, green800],
    begin: Alignment.topLeft,
    end: Alignment.bottomRight,
  );

  static LinearGradient get bgGradient => LinearGradient(
    colors: [Color(0xFF0D2818), Color(0xFF0F3020)],
    begin: Alignment.topLeft,
    end: Alignment.bottomRight,
  );
}
