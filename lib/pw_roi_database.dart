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
    relativeImpact: 5.0, // Extremely high RI — superpotent
    odourLifeHours: 168.0, // ~7 days; legendary fixative
    abcClass: 'U', // Animal/Amber in PW (Uw)
    abcLabel: 'U-Animal — Ambergris facet',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Bacdanol',
    relativeImpact: 50.0, // High impact sandalwood
    odourLifeHours: 72.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sandalwood facet',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Cashmeran',
    relativeImpact: 80.0, // Near-Linalool range
    odourLifeHours: 96.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Kashmiri/Musky wood',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Cedryl Acetate',
    relativeImpact: 200.0, // Low-medium; soft cedar
    odourLifeHours: 120.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Cedarwood acetate',
    source: 'LIT',
  ),
  PwRoiEntry(
    name: 'Galaxolide',
    relativeImpact: 30.0, // High RI musk — strong impact
    odourLifeHours: 120.0,
    abcClass: 'X',
    abcLabel: 'X-Musk — Polycyclic musk',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Iso E Super',
    relativeImpact: 40.0, // High RI; diffusive woody amber
    odourLifeHours: 72.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Woody cedarwood amber',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Javanol',
    relativeImpact: 30.0, // Very high impact sandalwood
    odourLifeHours: 96.0,
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
    relativeImpact: 60.0,
    odourLifeHours: 72.0,
    abcClass: 'W',
    abcLabel: 'W-Wood — Sandalwood core',
    source: 'LIT',
  ),

  // ── R — ROSE ─────────────────────────────────────────────
  PwRoiEntry(
    name: 'Citronellol',
    relativeImpact: 150.0,
    odourLifeHours: 24.0,
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
    relativeImpact: 120.0,
    odourLifeHours: 16.0,
    abcClass: 'R',
    abcLabel: 'R-Rose — Nerol / Rose-citrus',
    source: 'LIT',
  ),

  // ── L — LINALOOL ─────────────────────────────────────────
  PwRoiEntry(
    name: 'Linalool',
    relativeImpact: 100.0, // THE reference material
    odourLifeHours: 12.0,
    abcClass: 'L',
    abcLabel: 'L-Linalool — Fresh floral reference',
    source: 'PW',
  ),

  // ── J — JASMINE ──────────────────────────────────────────
  PwRoiEntry(
    name: 'Benzyl Acetate',
    relativeImpact: 200.0, // Moderate; jasmine solvent-like
    odourLifeHours: 8.0,
    abcClass: 'J',
    abcLabel: 'J-Jasmine — Jasmine/fruity',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Hedione',
    relativeImpact: 300.0, // Low impact; very diffusive
    odourLifeHours: 36.0,
    abcClass: 'J',
    abcLabel: 'J-Jasmine — Methyl dihydrojasmonate',
    source: 'PW',
  ),

  // ── V — VANILLA ──────────────────────────────────────────
  PwRoiEntry(
    name: 'Coumarin',
    relativeImpact: 80.0,
    odourLifeHours: 72.0,
    abcClass: 'V',
    abcLabel: 'V-Vanilla — Coumarin / hay',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Vanillin',
    relativeImpact: 50.0,
    odourLifeHours: 96.0,
    abcClass: 'V',
    abcLabel: 'V-Vanilla — Vanillin',
    source: 'PW',
  ),

  // ── I — IRIS / IONONE ────────────────────────────────────
  PwRoiEntry(
    name: 'a- ionone',
    relativeImpact: 10.0, // High impact violet/woody
    odourLifeHours: 48.0,
    abcClass: 'I',
    abcLabel: 'I-Iris — Alpha-Ionone / violet',
    source: 'PW',
  ),

  // ── S — SPICE ────────────────────────────────────────────
  PwRoiEntry(
    name: 'Cinnamaldehyde',
    relativeImpact: 30.0, // High impact spice
    odourLifeHours: 24.0,
    abcClass: 'S',
    abcLabel: 'S-Spice — Cinnamon aldehyde',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Eugenol',
    relativeImpact: 50.0,
    odourLifeHours: 48.0,
    abcClass: 'S',
    abcLabel: 'S-Spice — Eugenol / clove',
    source: 'PW',
  ),

  // ── P — PHENOL ───────────────────────────────────────────
  PwRoiEntry(
    name: 'Evernyl',
    relativeImpact: 40.0, // Oakmoss phenolic — high impact
    odourLifeHours: 120.0,
    abcClass: 'P',
    abcLabel: 'P-Phenol — Orcinol dimethyl ether / mossy',
    source: 'LIT',
  ),

  // ── G — GREEN ────────────────────────────────────────────
  PwRoiEntry(
    name: 'CIS-6-NONEAL',
    relativeImpact: 1.0, // Ultra-high impact green aldehyde
    odourLifeHours: 4.0, // Very short-lived top note
    abcClass: 'G',
    abcLabel: 'G-Green — Cucumber / green aldehyde',
    source: 'LIT',
  ),
  PwRoiEntry(
    name: 'Cis-3-Hexenyl Acetate',
    relativeImpact: 80.0,
    odourLifeHours: 6.0,
    abcClass: 'G',
    abcLabel: 'G-Green — Cis-3-hexenyl acetate / cut grass',
    source: 'LIT',
  ),

  // ── C — CITRUS ───────────────────────────────────────────
  PwRoiEntry(
    name: 'Limonene',
    relativeImpact: 500.0, // Low impact; high usage dose
    odourLifeHours: 2.0, // Evaporates very fast
    abcClass: 'C',
    abcLabel: 'C-Citrus — Limonene / orange peel',
    source: 'PW',
  ),

  // ── A — ALIPHATIC ────────────────────────────────────────
  PwRoiEntry(
    name: 'Dodecanal',
    relativeImpact: 15.0, // High impact fatty aldehyde
    odourLifeHours: 12.0,
    abcClass: 'A',
    abcLabel: 'A-Aliphatic — Dodecanal / fatty aldehydic',
    source: 'LIT',
  ),

  // ── F — FRUIT / ESTER ────────────────────────────────────
  PwRoiEntry(
    name: 'Isoamyl Acetate',
    relativeImpact: 150.0, // Medium; banana/pear
    odourLifeHours: 3.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Isoamyl acetate / banana',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Hexyl Acetate',
    relativeImpact: 200.0,
    odourLifeHours: 4.0,
    abcClass: 'F',
    abcLabel: 'F-Fruit — Hexyl acetate / fruity',
    source: 'LIT',
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
    relativeImpact: 60.0,
    odourLifeHours: 60.0,
    abcClass: 'X',
    abcLabel: 'X-Musk — Fruity musk',
    source: 'EST',
  ),

  // ── Q — ORIENTAL / AMBER ─────────────────────────────────
  PwRoiEntry(
    name: 'Romandolide',
    relativeImpact: 80.0,
    odourLifeHours: 120.0,
    abcClass: 'Q',
    abcLabel: 'Q-Orient — Macrocyclic musk / amber',
    source: 'LIT',
  ),
  PwRoiEntry(
    name: 'Macrolide',
    relativeImpact: 90.0,
    odourLifeHours: 100.0,
    abcClass: 'Q',
    abcLabel: 'Q-Orient — Macrocyclic lactone / musky',
    source: 'EST',
  ),

  // ── N — NARCOTIC ─────────────────────────────────────────
  PwRoiEntry(
    name: 'Benzaldehyde',
    relativeImpact: 30.0, // High impact; sharp almond
    odourLifeHours: 6.0,
    abcClass: 'N', // PW classifies under narcotic/almond
    abcLabel: 'N-Narcotic — Benzaldehyde / almond',
    source: 'PW',
  ),

  // ── Compound entries that span multiple classes ───────────
  PwRoiEntry(
    name: 'Menthol',
    relativeImpact: 50.0,
    odourLifeHours: 8.0,
    abcClass: 'B',
    abcLabel: 'B-iceBerg — Menthol / mint cooling',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'b-caryophyllenol',
    relativeImpact: 150.0,
    odourLifeHours: 60.0,
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
    relativeImpact: 2000.0, // Very low odour impact
    odourLifeHours: 0.5,
    abcClass: 'Z',
    abcLabel: 'Z-Solvent — Ethanol / carrier',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'Water',
    relativeImpact: 10000.0, // Near-odourless
    odourLifeHours: 0.1,
    abcClass: 'Z',
    abcLabel: 'Z-Solvent — Water',
    source: 'PW',
  ),
  PwRoiEntry(
    name: 'DPG',
    relativeImpact: 5000.0, // Odourless carrier
    odourLifeHours: 0.5,
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
