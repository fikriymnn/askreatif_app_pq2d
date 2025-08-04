import 'dart:math';
import 'package:askreatif_app/LoginPage.dart';
import 'package:askreatif_app/SolverPage.dart';
import 'package:askreatif_app/auth_guard.dart';
import 'package:askreatif_app/auth_service.dart';
import 'package:askreatif_app/group_params.dart';
import 'package:flutter/material.dart';
import 'package:fl_chart/fl_chart.dart';
import 'package:flutter/services.dart';

import 'amn_selected_groups.dart';

void main() => runApp(MyApp());

class MyApp extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return MaterialApp(
      debugShowCheckedModeBanner: false,
      title: 'POS Perfume Calculator',
      theme: ThemeData(
        primarySwatch: Colors.blue,
        visualDensity: VisualDensity.adaptivePlatformDensity,
      ),
      initialRoute: AuthService.isLoggedIn ? '/calculator' : '/login',
      onGenerateRoute: (settings) {
        switch (settings.name) {
          case '/login':
            return MaterialPageRoute(
              builder: (context) => LoginPage(),
              settings: settings,
            );
          case '/calculator':
            return MaterialPageRoute(
              builder: (context) => AuthGuard(child: SolverPage()),
              settings: settings,
            );
          default:
            // Route tidak ditemukan, redirect ke login
            return MaterialPageRoute(
              builder: (context) => LoginPage(),
              settings: settings,
            );
        }
      },
    );
  }
}
