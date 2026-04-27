import 'package:flutter/material.dart';
import 'dart:math' as math;
import 'package:askreatif_app/SolverPage.dart';

class WelcomePage extends StatefulWidget {
  @override
  _WelcomePageState createState() => _WelcomePageState();
}

class _WelcomePageState extends State<WelcomePage>
    with TickerProviderStateMixin {
  late AnimationController _fadeController;
  late AnimationController _floatController;
  late AnimationController _shimmerController;
  late Animation<double> _fadeAnim;
  late Animation<double> _floatAnim;
  late Animation<double> _shimmerAnim;

  @override
  void initState() {
    super.initState();

    _fadeController = AnimationController(
      vsync: this,
      duration: Duration(milliseconds: 1800),
    );
    _floatController = AnimationController(
      vsync: this,
      duration: Duration(seconds: 4),
    )..repeat(reverse: true);
    _shimmerController = AnimationController(
      vsync: this,
      duration: Duration(seconds: 3),
    )..repeat();

    _fadeAnim = CurvedAnimation(parent: _fadeController, curve: Curves.easeOut);
    _floatAnim = Tween<double>(begin: -8, end: 8).animate(
      CurvedAnimation(parent: _floatController, curve: Curves.easeInOut),
    );
    _shimmerAnim = _shimmerController;

    _fadeController.forward();
  }

  @override
  void dispose() {
    _fadeController.dispose();
    _floatController.dispose();
    _shimmerController.dispose();
    super.dispose();
  }

  @override
  Widget build(BuildContext context) {
    final size = MediaQuery.of(context).size;

    return Scaffold(
      body: Stack(
        children: [
          // Background gradient
          Container(
            decoration: BoxDecoration(
              gradient: LinearGradient(
                begin: Alignment.topLeft,
                end: Alignment.bottomRight,
                colors: [
                  Color(0xFF0D2818),
                  Color(0xFF1A3A2A),
                  Color(0xFF0F3020),
                  Color(0xFF082010),
                ],
                stops: [0.0, 0.35, 0.65, 1.0],
              ),
            ),
          ),

          // Decorative botanical circles
          Positioned(
            top: -80,
            right: -80,
            child: _BotanicalCircle(size: 320, opacity: 0.07),
          ),
          Positioned(
            bottom: -60,
            left: -60,
            child: _BotanicalCircle(size: 280, opacity: 0.06),
          ),
          Positioned(
            top: size.height * 0.4,
            right: -40,
            child: _BotanicalCircle(size: 180, opacity: 0.05),
          ),

          // Leaf pattern overlay
          CustomPaint(size: size, painter: _LeafPatternPainter()),

          // Main content
          FadeTransition(
            opacity: _fadeAnim,
            child: Row(
              children: [
                // LEFT PANEL - Branding
                Expanded(
                  flex: 5,
                  child: Container(
                    decoration: BoxDecoration(
                      border: Border(
                        right: BorderSide(
                          color: Color(0xFF4CAF50).withOpacity(0.2),
                          width: 1,
                        ),
                      ),
                    ),
                    child: Padding(
                      padding: EdgeInsets.symmetric(
                        horizontal: 60,
                        vertical: 40,
                      ),
                      child: Column(
                        mainAxisAlignment: MainAxisAlignment.center,
                        crossAxisAlignment: CrossAxisAlignment.start,
                        children: [
                          // Logo area
                          AnimatedBuilder(
                            animation: _floatAnim,
                            builder: (context, child) {
                              return Transform.translate(
                                offset: Offset(0, _floatAnim.value),
                                child: child,
                              );
                            },
                            child: _LogoWidget(),
                          ),

                          SizedBox(height: 48),

                          // Brand name
                          Text(
                            'ASKREATIF',
                            style: TextStyle(
                              fontFamily: 'Georgia',
                              fontSize: 52,
                              fontWeight: FontWeight.w300,
                              color: Colors.white,
                              letterSpacing: 24,
                              height: 1.0,
                            ),
                          ),
                          Text(
                            'Essential Oil',
                            style: TextStyle(
                              fontFamily: 'Georgia',
                              fontSize: 22,
                              fontWeight: FontWeight.w300,
                              color: Color(0xFF81C784),
                              letterSpacing: 8,
                            ),
                          ),
                          SizedBox(height: 4),
                          Text(
                            'Health and Beauty',
                            style: TextStyle(
                              fontSize: 13,
                              color: Colors.white.withOpacity(0.5),
                              letterSpacing: 4,
                            ),
                          ),

                          SizedBox(height: 48),

                          // Divider line
                          Container(
                            width: 80,
                            height: 1,
                            color: Color(0xFF4CAF50).withOpacity(0.5),
                          ),

                          SizedBox(height: 48),

                          // Description
                          Text(
                            'Perfumery Octonary\nSystem',
                            style: TextStyle(
                              fontFamily: 'Georgia',
                              fontSize: 32,
                              fontWeight: FontWeight.w400,
                              color: Colors.white.withOpacity(0.9),
                              height: 1.3,
                            ),
                          ),
                          SizedBox(height: 16),
                          Text(
                            'Formulasikan parfum sempurna dengan\nkalkulator UNIFAC & Perfumery Octonary System(POS) canggih.',
                            style: TextStyle(
                              fontSize: 15,
                              color: Colors.white.withOpacity(0.55),
                              height: 1.7,
                              letterSpacing: 0.3,
                            ),
                          ),
                        ],
                      ),
                    ),
                  ),
                ),

                // RIGHT PANEL - Features & CTA
                Expanded(
                  flex: 4,
                  child: Padding(
                    padding: EdgeInsets.symmetric(horizontal: 56, vertical: 40),
                    child: Column(
                      mainAxisAlignment: MainAxisAlignment.center,
                      crossAxisAlignment: CrossAxisAlignment.start,
                      children: [
                        Text(
                          'Selamat Datang',
                          style: TextStyle(
                            fontSize: 13,
                            color: Color(0xFF81C784),
                            letterSpacing: 4,
                          ),
                        ),
                        SizedBox(height: 16),
                        Text(
                          'Formulasi Parfum\nDigital',
                          style: TextStyle(
                            fontFamily: 'Georgia',
                            fontSize: 40,
                            fontWeight: FontWeight.w300,
                            color: Colors.white,
                            height: 1.2,
                          ),
                        ),

                        SizedBox(height: 40),

                        // Feature cards
                        _FeatureCard(
                          icon: Icons.science_outlined,
                          title: 'UNIFAC Activity Coefficient',
                          desc:
                              'Kalkulasi koefisien aktivitas dengan model termodinamika presisi tinggi',
                          delay: 200,
                        ),
                        SizedBox(height: 12),
                        _FeatureCard(
                          icon: Icons.bubble_chart_outlined,
                          title: 'Octonary System (POS)',
                          desc:
                              'Optimasi 6 senyawa pewangi + 2 pelarut dalam satu sistem',
                          delay: 350,
                        ),
                        SizedBox(height: 12),
                        _FeatureCard(
                          icon: Icons.auto_graph,
                          title: 'Diagram Ternary & Octonary',
                          desc:
                              'Visualisasi komposisi dan nilai OV secara interaktif',
                          delay: 500,
                        ),

                        SizedBox(height: 48),

                        // CTA Button
                        AnimatedBuilder(
                          animation: _shimmerAnim,
                          builder: (context, child) {
                            return _CTAButton(
                              onPressed: () {
                                Navigator.of(context).pushReplacement(
                                  PageRouteBuilder(
                                    pageBuilder: (_, anim, __) => SolverPage(),
                                    transitionsBuilder: (_, anim, __, child) {
                                      return FadeTransition(
                                        opacity: anim,
                                        child: child,
                                      );
                                    },
                                    transitionDuration: Duration(
                                      milliseconds: 600,
                                    ),
                                  ),
                                );
                              },
                            );
                          },
                        ),

                        SizedBox(height: 24),

                        // Version info
                        Text(
                          'POS Calculator v2.0  ·  © ASKREATIF Essential Oil',
                          style: TextStyle(
                            fontSize: 11,
                            color: Colors.white.withOpacity(0.25),
                            letterSpacing: 1,
                          ),
                        ),
                      ],
                    ),
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

class _CTAButton extends StatefulWidget {
  final VoidCallback onPressed;
  const _CTAButton({required this.onPressed});

  @override
  __CTAButtonState createState() => __CTAButtonState();
}

class __CTAButtonState extends State<_CTAButton> {
  bool _hovered = false;

  @override
  Widget build(BuildContext context) {
    return MouseRegion(
      onEnter: (_) => setState(() => _hovered = true),
      onExit: (_) => setState(() => _hovered = false),
      child: GestureDetector(
        onTap: widget.onPressed,
        child: AnimatedContainer(
          duration: Duration(milliseconds: 250),
          width: double.infinity,
          height: 56,
          decoration: BoxDecoration(
            gradient: LinearGradient(
              colors:
                  _hovered
                      ? [Color(0xFF2E7D32), Color(0xFF43A047)]
                      : [Color(0xFF1B5E20), Color(0xFF2E7D32)],
            ),
            borderRadius: BorderRadius.circular(4),
            boxShadow:
                _hovered
                    ? [
                      BoxShadow(
                        color: Color(0xFF2E7D32).withOpacity(0.5),
                        blurRadius: 20,
                        offset: Offset(0, 6),
                      ),
                    ]
                    : [],
            border: Border.all(
              color: Color(0xFF4CAF50).withOpacity(_hovered ? 0.8 : 0.4),
              width: 1,
            ),
          ),
          child: Row(
            mainAxisAlignment: MainAxisAlignment.center,
            children: [
              Text(
                'MULAI KALKULASI',
                style: TextStyle(
                  color: Colors.white,
                  fontSize: 13,
                  fontWeight: FontWeight.w600,
                  letterSpacing: 3,
                ),
              ),
              SizedBox(width: 12),
              AnimatedContainer(
                duration: Duration(milliseconds: 250),
                transform: Matrix4.translationValues(_hovered ? 4 : 0, 0, 0),
                child: Icon(Icons.arrow_forward, color: Colors.white, size: 18),
              ),
            ],
          ),
        ),
      ),
    );
  }
}

class _FeatureCard extends StatefulWidget {
  final IconData icon;
  final String title;
  final String desc;
  final int delay;

  const _FeatureCard({
    required this.icon,
    required this.title,
    required this.desc,
    required this.delay,
  });

  @override
  __FeatureCardState createState() => __FeatureCardState();
}

class __FeatureCardState extends State<_FeatureCard> {
  bool _hovered = false;

  @override
  Widget build(BuildContext context) {
    return MouseRegion(
      onEnter: (_) => setState(() => _hovered = true),
      onExit: (_) => setState(() => _hovered = false),
      child: AnimatedContainer(
        duration: Duration(milliseconds: 200),
        padding: EdgeInsets.all(16),
        decoration: BoxDecoration(
          color:
              _hovered
                  ? Colors.white.withOpacity(0.08)
                  : Colors.white.withOpacity(0.04),
          borderRadius: BorderRadius.circular(4),
          border: Border.all(
            color:
                _hovered
                    ? Color(0xFF4CAF50).withOpacity(0.4)
                    : Colors.white.withOpacity(0.08),
            width: 1,
          ),
        ),
        child: Row(
          crossAxisAlignment: CrossAxisAlignment.start,
          children: [
            Container(
              padding: EdgeInsets.all(8),
              decoration: BoxDecoration(
                color: Color(0xFF2E7D32).withOpacity(0.4),
                borderRadius: BorderRadius.circular(3),
              ),
              child: Icon(widget.icon, color: Color(0xFF81C784), size: 18),
            ),
            SizedBox(width: 14),
            Expanded(
              child: Column(
                crossAxisAlignment: CrossAxisAlignment.start,
                children: [
                  Text(
                    widget.title,
                    style: TextStyle(
                      color: Colors.white.withOpacity(0.9),
                      fontSize: 12,
                      fontWeight: FontWeight.w600,
                      letterSpacing: 0.5,
                    ),
                  ),
                  SizedBox(height: 3),
                  Text(
                    widget.desc,
                    style: TextStyle(
                      color: Colors.white.withOpacity(0.45),
                      fontSize: 11,
                      height: 1.5,
                    ),
                  ),
                ],
              ),
            ),
          ],
        ),
      ),
    );
  }
}

class _LogoWidget extends StatelessWidget {
  @override
  Widget build(BuildContext context) {
    return Container(
      width: 90,
      height: 90,
      decoration: BoxDecoration(
        shape: BoxShape.circle,
        gradient: RadialGradient(
          colors: [
            Color(0xFF2E7D32).withOpacity(0.6),
            Color(0xFF1B5E20).withOpacity(0.3),
          ],
        ),
        border: Border.all(
          color: Color(0xFF4CAF50).withOpacity(0.5),
          width: 1.5,
        ),
        boxShadow: [
          BoxShadow(
            color: Color(0xFF2E7D32).withOpacity(0.4),
            blurRadius: 24,
            spreadRadius: 4,
          ),
        ],
      ),
      child: Center(
        child: Text(
          'ask',
          style: TextStyle(
            fontFamily: 'Georgia',
            fontSize: 28,
            fontWeight: FontWeight.w400,
            color: Colors.white.withOpacity(0.95),
            fontStyle: FontStyle.italic,
          ),
        ),
      ),
    );
  }
}

class _BotanicalCircle extends StatelessWidget {
  final double size;
  final double opacity;
  const _BotanicalCircle({required this.size, required this.opacity});

  @override
  Widget build(BuildContext context) {
    return Opacity(
      opacity: opacity,
      child: Container(
        width: size,
        height: size,
        decoration: BoxDecoration(
          shape: BoxShape.circle,
          border: Border.all(color: Color(0xFF4CAF50), width: 1),
        ),
        child: Center(
          child: Container(
            width: size * 0.7,
            height: size * 0.7,
            decoration: BoxDecoration(
              shape: BoxShape.circle,
              border: Border.all(color: Color(0xFF4CAF50), width: 1),
            ),
            child: Center(
              child: Container(
                width: size * 0.4,
                height: size * 0.4,
                decoration: BoxDecoration(
                  shape: BoxShape.circle,
                  border: Border.all(color: Color(0xFF4CAF50), width: 1),
                ),
              ),
            ),
          ),
        ),
      ),
    );
  }
}

class _LeafPatternPainter extends CustomPainter {
  @override
  void paint(Canvas canvas, Size size) {
    final paint =
        Paint()
          ..color = Color(0xFF4CAF50).withOpacity(0.04)
          ..strokeWidth = 1
          ..style = PaintingStyle.stroke;

    // Draw subtle diagonal lines
    for (double i = -size.height; i < size.width + size.height; i += 60) {
      canvas.drawLine(
        Offset(i, 0),
        Offset(i + size.height, size.height),
        paint,
      );
    }
  }

  @override
  bool shouldRepaint(_LeafPatternPainter old) => false;
}
