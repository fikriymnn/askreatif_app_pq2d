import 'package:askreatif_app/auth_service.dart';
import 'package:flutter/material.dart';

class AuthGuard extends StatelessWidget {
  final Widget child;

  const AuthGuard({Key? key, required this.child}) : super(key: key);

  @override
  Widget build(BuildContext context) {
    // Cek apakah user sudah login
    if (AuthService.isLoggedIn) {
      return child;
    } else {
      // Jika belum login, redirect ke login page
      WidgetsBinding.instance.addPostFrameCallback((_) {
        Navigator.of(
          context,
        ).pushNamedAndRemoveUntil('/login', (Route<dynamic> route) => false);
      });

      // Tampilkan loading sementara redirect
      return Scaffold(body: Center(child: CircularProgressIndicator()));
    }
  }
}
