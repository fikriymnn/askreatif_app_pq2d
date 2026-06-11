// ============================================================
// ASKREATIF Perfumery Engine
// PerfumersWorld-Calibrated ROI + Odour Life Database
// ============================================================
//
// SCIENTIFIC BASIS
// ─────────────────────────────────────────────────────────────
// PerfumersWorld defines "Relative Impact" (RI) as the
// amount of a compound (in µL) needed to match the perceived
// impact of a fixed dose of Linalool Synthetic on a smelling
// strip.  Linalool = 100 (midpoint reference).
//
// Scale interpretation (exponential):
//   RI < 10       → Ultra-high impact (use trace amounts)
//   RI 10–50      → High impact
//   RI 50–150     → Medium impact (Linalool zone)
//   RI 150–500    → Low impact
//   RI > 500      → Very low impact (fixative/blender range)
//
// "Odour Life" (OL) = hours on chromatography paper strip
// until the material becomes weak and uncharacteristic.
//
// DATA SOURCES
// ─────────────────────────────────────────────────────────────
// 1. PerfumersWorld ABCs system (Dowthwaite, 1999)
// 2. Arctander "Perfume and Flavor Chemicals" (1969)
// 3. Leffingwell & Associates Flavour-Base / Odour Units
// 4. BACIS / Givaudan Impact data (published literature)
// 5. Estimated from UNIFAC Psat / Threshold correlation
//    for compounds not directly listed in published tables.
//
// NOTE: Where exact PW values require login, values are
// cross-referenced with Arctander + Leffingwell and annotated.
// ─────────────────────────────────────────────────────────────

class PwRoiEntry {
  final String name;

  /// PerfumersWorld Relative Impact (RI)
  /// Linalool = 100 (reference midpoint)
  final double relativeImpact;

  /// Odour Life on smelling strip (hours)
  final double odourLifeHours;

  /// ABCs classification letter (PerfumersWorld system)
  final String abcClass;

  /// Full ABC label
  final String abcLabel;

  /// Source confidence: 'PW' = direct PW, 'LIT' = literature, 'EST' = estimated
  final String source;

  const PwRoiEntry({
    required this.name,
    required this.relativeImpact,
    required this.odourLifeHours,
    required this.abcClass,
    required this.abcLabel,
    this.source = 'LIT',
  });
}

// ── MASTER DATABASE ────────────────────────────────────────
// Ordered alphabetically. All 40 fragrance compounds from
// Compound.dart plus the 3 solvent compounds.
const List<PwRoiEntry> pwRoiDatabase = [
  // ── W — WOOD ─────────────────────────────────────────────
  PwRoiEntry(
    name: 'Ambroxan',
    relativeImpact: 400, // Extremely high RI — superpotent
    odourLifeHours: 750, // ~7 days; legendary fixative
    abcClass: 'U', // Animal/Amber in PW (Uw)
    abcLabel: 'U-Animal — Ambergris facet',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Bacdanol',
    relativeImpact: 80.0, // High impact sandalwood
    odourLifeHours: 300.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sandalwood facet',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Cashmeran',
    relativeImpact: 80.0, // Near-Linalool range
    odourLifeHours: 200.0,
    abcClass: 'X',
    abcLabel: 'X-Musk Sexy, Musk, Sensual Skin-like',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Cedryl Acetate',
    relativeImpact: 60.0, // Low-medium; soft cedar
    odourLifeHours: 40.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Cedarwood acetate',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Galaxolide',
    relativeImpact: 60.0, // High RI musk — strong impact
    odourLifeHours: 200.0,
    abcClass: 'X',
    abcLabel: 'X-Musk — Polycyclic musk',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Iso E Super',
    relativeImpact: 40.0, // High RI; diffusive woody amber
    odourLifeHours: 260.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Woody cedarwood amber',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Javanol',
    relativeImpact: 250.0, // Very high impact sandalwood
    odourLifeHours: 400.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sandalwood (Givaudan)',
    source: 'LIT',
  ),
  PwRoiEntry(
    name: 'Patchouli Alcohol',
    relativeImpact: 100.0, // = Linalool reference
    odourLifeHours: 240.0, // Very long-lasting base
    abcClass: 'W',
    abcLabel: 'W-Wood — Patchouli',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Sandalmysore Core',
    relativeImpact: 98.0,
    odourLifeHours: 42.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sandalwood core',
    source: 'PW',
  ),

  // ── R — ROSE ─────────────────────────────────────────────
  PwRoiEntry(
    name: 'Citronellol',
    relativeImpact: 100.0,
    odourLifeHours: 20.0,
    abcClass: 'R',
    abcLabel: 'R-Rose — Citronellol rose',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Geraniol',
    relativeImpact: 100.0, // Very close to Linalool reference
    odourLifeHours: 18.0,
    abcClass: 'R',
    abcLabel: 'R-Rose — Geraniol / Rose',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Nerol',
    relativeImpact: 100.0,
    odourLifeHours: 18.0,
    abcClass: 'C',
    abcLabel: 'C-Citrus Sour Sharp Citrus',
    source: 'PW',
  ),

  // ── L — LINALOOL ─────────────────────────────────────────
  PwRoiEntry(
    name: 'Linalool',
    relativeImpact: 100.0, // THE reference material
    odourLifeHours: 10.0,
    abcClass: 'L',
    abcLabel: 'L-Linalool — Fresh floral reference',
    source: 'PW',
  ),

  // ── J — JASMINE ──────────────────────────────────────────
  PwRoiEntry(
    name: 'Benzyl Acetate',
    relativeImpact: 120.0, // Moderate; jasmine solvent-like
    odourLifeHours: 15.0,
    abcClass: 'J',
    abcLabel: 'J-Jasmine — Jasmine/fruity',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Hedione',
    relativeImpact: 65.0, // Low impact; very diffusive
    odourLifeHours: 160.0,
    abcClass: 'J',
    abcLabel: 'J-Jasmine — Methyl dihydrojasmonate',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Alpha Hexyl Cinnamaldehyde',
    relativeImpact: 100.0,
    odourLifeHours: 460.0,
    abcClass: 'J',
    abcLabel: 'J-Jasmine — Methyl dihydrojasmonate',
    source: 'PW',
  ),

  // ── V — VANILLA ──────────────────────────────────────────
  PwRoiEntry(
    name: 'Coumarin',
    relativeImpact: 120.0,
    odourLifeHours: 85.0,
    abcClass: 'V',
    abcLabel: 'V-Vanilla — Coumarin / hay',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Vanillin',
    relativeImpact: 60.0,
    odourLifeHours: 400.0,
    abcClass: 'V',
    abcLabel: 'V-Vanilla — Vanillin',
    source: 'PW',
  ),

  // ── I — IRIS / IONONE ────────────────────────────────────
  PwRoiEntry(
    name: 'a- ionone',
    relativeImpact: 100.0, // High impact violet/woody
    odourLifeHours: 100.0,
    abcClass: 'I',
    abcLabel: 'I-Iris — Alpha-Ionone / violet',
    source: 'PW',
  ),

  // ── S — SPICE ────────────────────────────────────────────
  PwRoiEntry(
    name: 'Cinnamaldehyde',
    relativeImpact: 130.0, // High impact spice
    odourLifeHours: 75.0,
    abcClass: 'S',
    abcLabel: 'S-Spice — Cinnamon aldehyde',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Eugenol',
    relativeImpact: 175.0,
    odourLifeHours: 48.0,
    abcClass: 'S',
    abcLabel: 'S-Spice — Eugenol / clove',
    source: 'PW',
  ),

  // ── P — PHENOL ───────────────────────────────────────────
  PwRoiEntry(
    name: 'Evernyl',
    relativeImpact: 60.0, // Oakmoss phenolic — high impact
    odourLifeHours: 240.0,
    abcClass: 'Y',
    abcLabel: 'Y-Yeast Mossy Fungal Marine Ozone',
    source: 'PW',
  ),

  // ── G — GREEN ────────────────────────────────────────────
  PwRoiEntry(
    name: 'CIS-6-NONEAL',
    relativeImpact: 700, // Ultra-high impact green aldehyde
    odourLifeHours: 1.5, // Very short-lived top note
    abcClass: 'F',
    abcLabel: 'F-Fruit Sweet & Sour Berry to Banana',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Cis-3-Hexenyl Acetate',
    relativeImpact: 800.0,
    odourLifeHours: 2.0,
    abcClass: 'G',
    abcLabel: 'G-Green — Cis-3-hexenyl acetate / cut grass',
    source: 'PW',
  ),

  // ── C — CITRUS ───────────────────────────────────────────
  PwRoiEntry(
    name: 'Limonene',
    relativeImpact: 70.0, // Low impact; high usage dose
    odourLifeHours: 1.0, // Evaporates very fast
    abcClass: 'C',
    abcLabel: 'C-Citrus — Limonene / orange peel',
    source: 'PW',
  ),

  // ── A — ALIPHATIC ────────────────────────────────────────
  PwRoiEntry(
    name: 'Dodecanal',
    relativeImpact: 600, // High impact fatty aldehyde
    odourLifeHours: 80,
    abcClass: 'A',
    abcLabel: 'A-Aliphatic — Dodecanal / fatty aldehydic',
    source: 'PW',
  ),

  // ── F — FRUIT / ESTER ────────────────────────────────────
  PwRoiEntry(
    name: 'Isoamyl Acetate',
    relativeImpact: 800.0, // Medium; banana/pear
    odourLifeHours: 0.1,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Isoamyl acetate / banana',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Hexyl Acetate',
    relativeImpact: 220.0,
    odourLifeHours: 0.2,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Hexyl acetate / fruity',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Heptyl Acetate',
    relativeImpact: 200.0,
    odourLifeHours: 6.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Heptyl acetate / fruity herbal',
    source: 'EST',
  ),
  PwRoiEntry(
    name: 'Octyl Acetate',
    relativeImpact: 180.0,
    odourLifeHours: 8.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Octyl acetate / orange-fruity',
    source: 'EST',
  ),
  PwRoiEntry(
    name: 'Oenanthic Ether',
    relativeImpact: 160.0,
    odourLifeHours: 6.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Ethyl heptanoate / cognac',
    source: 'LIT',
  ),

  // ── X — MUSK ─────────────────────────────────────────────
  PwRoiEntry(
    name: 'Frukton',
    relativeImpact: 80.0,
    odourLifeHours: 2.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Fruit Sweet & Sour Berry to Banana',
    source: 'EST',
  ),

  // ── Q — ORIENTAL / AMBER ─────────────────────────────────
  PwRoiEntry(
    name: 'Romandolide',
    relativeImpact: 100.0,
    odourLifeHours: 400.0,
    abcClass: 'X',
    abcLabel: 'X-Musk Sexy, Musk, Sensual Skin-like',
    source: 'LIT',
  ),
  PwRoiEntry(
    name: 'Macrolide',
    relativeImpact: 250.0,
    odourLifeHours: 400.0,
    abcClass: 'X',
    abcLabel: 'X-Musk Sexy, Musk, Sensual Skin-like',
    source: 'EST',
  ),

  // ── M — MUGUET ─────────────────────────────────────────
  PwRoiEntry(
    name: 'Lily Aldehyde',
    relativeImpact: 80.0, // Moderate impact; floral green
    odourLifeHours: 72.0,
    abcClass: 'M',
    abcLabel: 'M-Muguet — Lily / floral green',
    source: 'PW',
  ),

  // ── N — NARCOTIC ─────────────────────────────────────────
  PwRoiEntry(
    name: 'Benzaldehyde',
    relativeImpact: 500.0, // High impact; sharp almond
    odourLifeHours: 1.5,
    abcClass: 'F', // PW classifies under narcotic/almond
    abcLabel: 'F-Fruit — Benzaldehyde / almond cherry',
    source: 'PW',
  ),

  // ── Compound entries that span multiple classes ───────────
  PwRoiEntry(
    name: 'Menthol',
    relativeImpact: 120.0,
    odourLifeHours: 18.0,
    abcClass: 'B',
    abcLabel: 'B-iceBerg — Menthol / mint cooling',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'b-caryophyllenol',
    relativeImpact: 50.0,
    odourLifeHours: 28.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sesquiterpene alcohol / woody spicy',
    source: 'EST',
  ),
  PwRoiEntry(
    name: 'Caryophyllene Oxide',
    relativeImpact: 200.0,
    odourLifeHours: 48.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Caryophyllene oxide / spicy woody',
    source: 'EST',
  ),
  PwRoiEntry(
    name: 'Helvetolida',
    relativeImpact: 120.0,
    odourLifeHours: 80.0,
    abcClass: 'X',
    abcLabel: 'X-Musk — Helvetolide macrocyclic',
    source: 'EST',
  ),

  // ── SOLVENTS / CARRIERS ──────────────────────────────────
  PwRoiEntry(
    name: 'Ethanol',
    relativeImpact: 0, // Very low odour impact
    odourLifeHours: 0.05,
    abcClass: 'Z',
    abcLabel: 'Z-Solvent — Ethanol / carrier',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Water',
    relativeImpact: 0, // Near-odourless
    odourLifeHours: 0.01,
    abcClass: 'Z',
    abcLabel: 'Z-Solvent — Water',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'DPG',
    relativeImpact: 3.0, // Odourless carrier
    odourLifeHours: 7,
    abcClass: 'Z',
    abcLabel: 'Z-Solvent — Dipropylene Glycol / carrier',
    source: 'PW',
  ),
];

// ── Lookup helpers ─────────────────────────────────────────

/// Get ROI entry by compound name (null if not found).
PwRoiEntry? pwRoiFor(String name) {
  try {
    return pwRoiDatabase.firstWhere(
      (e) => e.name.toLowerCase() == name.toLowerCase(),
    );
  } catch (_) {
    return null;
  }
}

/// Relative Impact value; returns 100.0 (neutral/Linalool)
/// if compound is unknown.
double relativeImpactOf(String name) => pwRoiFor(name)?.relativeImpact ?? 100.0;

/// Odour life in hours; returns 12.0 (Linalool default)
/// if compound is unknown.
double odourLifeOf(String name) => pwRoiFor(name)?.odourLifeHours ?? 12.0;

/// PW ABC class letter (A–Z).
String abcClassOf(String name) => pwRoiFor(name)?.abcClass ?? '?';

// ── PW-Corrected ROI for mixture ──────────────────────────
//
// The PerfumersWorld "Impact" is inversely related to the
// Relative Impact number:
//   PerceptualStrength_i ∝ 1 / RI_i
//
// This is because RI measures "how much you need" — a lower
// RI means you need less (stronger material).
//
// Integration with the thermodynamic ROI formula:
//   ROI_PW_i = (x_i × Psat_i) / (MW_i × Thr_i × RI_i / 100)
//
// The RI/100 normalises relative to Linalool.

double pwCorrectedRoi({
  required String name,
  required double moleFraction,
  required double psat,
  required double mw,
  required double threshold,
}) {
  final double ri = relativeImpactOf(name);
  if (ri <= 0 || mw <= 0 || threshold <= 0) return 0.0;
  final double roi = (moleFraction * psat) / (mw * threshold * (ri / 100.0));
  return roi.isNaN || roi.isInfinite ? 0.0 : roi;
}
