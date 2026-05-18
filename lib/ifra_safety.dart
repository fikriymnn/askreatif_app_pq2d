// ============================================================
// ASKREATIF Perfumery Engine — Phase 9
// IFRA Safety & Regulatory System
// ============================================================
// IFRA 51st Amendment limits (2023) — leave-on skin product
// category as primary reference. Extend with other categories.
// Allergen flags per EU Regulation (EC) No 1223/2009.
// Phototoxicity warnings based on published data.
// ============================================================
// ⚠️  Data is provided for RESEARCH purposes only.
//     Always verify against the official IFRA standard and
//     consult a qualified regulatory toxicologist before
//     commercial use.
// ============================================================

// ── IFRA Product Categories ────────────────────────────────
enum IfraCategory {
  cat1,  // Hydroalcoholic products applied to non-protected skin (EDT etc.)
  cat2,  // Deodorant, body spray
  cat4,  // Fine fragrance (Parfum, EDP)
  cat11, // Rinse-off: shampoo, bath gel
}

// ── IFRA Limit Entry ──────────────────────────────────────
class IfraLimit {
  final String compound;
  final Map<IfraCategory, double> limits; // % w/w in finished product
  final bool isRestricted;
  final bool isBanned;
  final String? note;

  const IfraLimit({
    required this.compound,
    required this.limits,
    this.isRestricted = false,
    this.isBanned = false,
    this.note,
  });
}

// ── IFRA Database (51st Amendment, selected compounds) ─────
const List<IfraLimit> ifraDatabase = [
  IfraLimit(
    compound: 'Geraniol',
    limits: {
      IfraCategory.cat1: 6.4,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 5.3,
    },
    isRestricted: true,
    note: 'Skin sensitiser at high concentrations.',
  ),
  IfraLimit(
    compound: 'Linalool',
    limits: {
      IfraCategory.cat1: 100.0,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 100.0,
    },
  ),
  IfraLimit(
    compound: 'Citronellol',
    limits: {
      IfraCategory.cat1: 13.0,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 5.0,
    },
    isRestricted: true,
  ),
  IfraLimit(
    compound: 'Eugenol',
    limits: {
      IfraCategory.cat1: 0.5,
      IfraCategory.cat4: 0.5,
      IfraCategory.cat11: 0.5,
    },
    isRestricted: true,
    note: 'Strong allergen; EU-listed sensitiser; Declaration required ≥0.001%.',
  ),
  IfraLimit(
    compound: 'Cinnamaldehyde',
    limits: {
      IfraCategory.cat1: 0.05,
      IfraCategory.cat4: 0.05,
      IfraCategory.cat11: 0.2,
    },
    isRestricted: true,
    note: 'Potent sensitiser. Likely prohibited in leave-on above 0.05%.',
  ),
  IfraLimit(
    compound: 'Isoamyl Acetate',
    limits: {
      IfraCategory.cat1: 100.0,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 100.0,
    },
  ),
  IfraLimit(
    compound: 'Coumarin',
    limits: {
      IfraCategory.cat1: 0.3,
      IfraCategory.cat4: 0.3,
      IfraCategory.cat11: 0.5,
    },
    isRestricted: true,
    note: 'EU-listed allergen; Declaration required ≥0.001% rinse-off, ≥0.01% leave-on.',
  ),
  IfraLimit(
    compound: 'Benzaldehyde',
    limits: {
      IfraCategory.cat1: 0.1,
      IfraCategory.cat4: 0.1,
      IfraCategory.cat11: 3.0,
    },
    isRestricted: true,
  ),
  IfraLimit(
    compound: 'Benzyl Acetate',
    limits: {
      IfraCategory.cat1: 4.9,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 100.0,
    },
  ),
  IfraLimit(
    compound: 'Limonene',
    limits: {
      IfraCategory.cat1: 1.0,
      IfraCategory.cat4: 10.0,
      IfraCategory.cat11: 2.0,
    },
    isRestricted: true,
    note: 'Oxidises to allergen limonene hydroperoxide; use antioxidant.',
  ),
  IfraLimit(
    compound: 'Vanillin',
    limits: {
      IfraCategory.cat1: 1.5,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 100.0,
    },
  ),
  IfraLimit(
    compound: 'Galaxolide',
    limits: {
      IfraCategory.cat1: 7.0,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 3.9,
    },
  ),
  IfraLimit(
    compound: 'Ambroxan',
    limits: {
      IfraCategory.cat1: 100.0,
      IfraCategory.cat4: 100.0,
      IfraCategory.cat11: 100.0,
    },
  ),
];

// ── EU Allergen Declaration List ──────────────────────────
const Set<String> euAllergens = {
  'Geraniol',
  'Linalool',
  'Citronellol',
  'Eugenol',
  'Cinnamaldehyde',
  'Coumarin',
  'Limonene',
  'Benzaldehyde',
};

// ── Phototoxic Compounds ──────────────────────────────────
const Set<String> phototoxicCompounds = {
  // Furanocoumarins, bergapten-containing materials
  // These are fragrance raw materials, not explicitly in our
  // compound list, but included as a safeguard pattern.
  'Benzyl Acetate', // mild phototoxicity reports
};

// ── Safety Check Result ───────────────────────────────────
class SafetyAlert {
  final String compound;
  final AlertSeverity severity;
  final String message;
  final double? limitPct;      // IFRA limit for category
  final double? actualPct;     // Actual concentration in formula

  const SafetyAlert({
    required this.compound,
    required this.severity,
    required this.message,
    this.limitPct,
    this.actualPct,
  });
}

enum AlertSeverity { info, warning, critical }

// ── Safety Checker ────────────────────────────────────────
class FormulaSafetyChecker {
  final IfraCategory category;

  const FormulaSafetyChecker({
    this.category = IfraCategory.cat1,
  });

  /// Check a formula's safety. Returns list of alerts (may be empty).
  List<SafetyAlert> check({
    required List<String> compoundNames,
    required List<double> massFractionsPct, // 0–100 w/w%
  }) {
    final List<SafetyAlert> alerts = [];

    for (int i = 0; i < compoundNames.length; i++) {
      final String name = compoundNames[i];
      final double actualPct = massFractionsPct[i];

      // ── IFRA limit check ──────────────────────────────
      final IfraLimit? entry = _findEntry(name);
      if (entry != null) {
        if (entry.isBanned) {
          alerts.add(SafetyAlert(
            compound: name,
            severity: AlertSeverity.critical,
            message: '$name is BANNED by IFRA. Remove from formula.',
          ));
          continue;
        }
        final double? limit = entry.limits[category];
        if (limit != null && actualPct > limit) {
          alerts.add(SafetyAlert(
            compound: name,
            severity: AlertSeverity.critical,
            message: '$name exceeds IFRA limit for category '
                '${_catName(category)}: '
                '${actualPct.toStringAsFixed(2)}% > ${limit.toStringAsFixed(2)}%.',
            limitPct: limit,
            actualPct: actualPct,
          ));
        } else if (entry.isRestricted) {
          alerts.add(SafetyAlert(
            compound: name,
            severity: AlertSeverity.warning,
            message: '$name is restricted. '
                '${entry.note ?? "Check IFRA 51st Amendment."}',
            limitPct: limit,
            actualPct: actualPct,
          ));
        }
      }

      // ── EU Allergen declaration ───────────────────────
      if (euAllergens.contains(name) && actualPct >= 0.01) {
        alerts.add(SafetyAlert(
          compound: name,
          severity: AlertSeverity.info,
          message: '$name is an EU-listed allergen. '
              'Declaration required on label (≥0.01% leave-on).',
          actualPct: actualPct,
        ));
      }

      // ── Phototoxicity ─────────────────────────────────
      if (phototoxicCompounds.contains(name)) {
        alerts.add(SafetyAlert(
          compound: name,
          severity: AlertSeverity.warning,
          message: '$name has reported phototoxicity. '
              'Avoid high sun-exposure application.',
          actualPct: actualPct,
        ));
      }
    }

    // Sort: critical → warning → info
    alerts.sort((a, b) =>
        b.severity.index.compareTo(a.severity.index));
    return alerts;
  }

  IfraLimit? _findEntry(String name) {
    try {
      return ifraDatabase.firstWhere((e) => e.compound == name);
    } catch (_) {
      return null;
    }
  }

  String _catName(IfraCategory c) {
    switch (c) {
      case IfraCategory.cat1:  return 'Cat 1 (EDT/Spray)';
      case IfraCategory.cat2:  return 'Cat 2 (Deo Spray)';
      case IfraCategory.cat4:  return 'Cat 4 (Parfum)';
      case IfraCategory.cat11: return 'Cat 11 (Rinse-off)';
    }
  }
}