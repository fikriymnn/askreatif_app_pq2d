class AuthService {
  static bool _isLoggedIn = false;

  static bool get isLoggedIn => _isLoggedIn;

  static void login() {
    _isLoggedIn = true;
  }

  static void logout() {
    _isLoggedIn = false;
  }
}
