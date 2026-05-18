// ============================================================
// ASKREATIF Perfumery Engine — Phase 7
// Machine Learning Infrastructure Layer
// ============================================================
// Clean data-capture layer for future ML/AI integration.
// Stores formula history, OV results, perceptual ratings,
// and accord profiling for downstream model training.
// ============================================================
// ⚠️  No fake AI here. This is a proper data substrate for
//     future recommendation / generative perfumery AI.
// ============================================================

import 'dart:convert';
import 'package:askreatif_app/Compound.dart';

// ── Formula Record ─────────────────────────────────────────
class FormulaRecord {
  final String id;                 // UUID or timestamp-based key
  final DateTime createdAt;
  final String solventSystem;      // e.g. "Water + Ethanol"
  final double totalSolventFraction;

  // Input
  final List<String> compoundNames;
  final List<double> moleFractions;

  // Thermodynamic outputs
  final List<double> activityCoefficients;
  final List<double> odorValues;
  final List<double> roiWeights;
  final List<double> headspaceConc; // mg/m³ at t=0

  // Accord profile
  final List<String> topAccords;   // e.g. ['Floral', 'Woody', 'Amber']
  final Map<String, double> accordBalance;

  // User ratings (nullable until user provides them)
  final int? longevityRating;      // 1–10
  final int? freshRating;          // 1–10
  final int? luxuryRating;         // 1–10
  final int? overallRating;        // 1–10
  final String? notes;             // free-text perfumer notes

  // Projected performance (physics-derived)
  final double projectedLongevityHours; // hours until OV < threshold
  final double ovBalance;              // std dev of log(OV) — lower = better

  const FormulaRecord({
    required this.id,
    required this.createdAt,
    required this.solventSystem,
    required this.totalSolventFraction,
    required this.compoundNames,
    required this.moleFractions,
    required this.activityCoefficients,
    required this.odorValues,
    required this.roiWeights,
    required this.headspaceConc,
    required this.topAccords,
    required this.accordBalance,
    this.longevityRating,
    this.freshRating,
    this.luxuryRating,
    this.overallRating,
    this.notes,
    required this.projectedLongevityHours,
    required this.ovBalance,
  });

  /// Deep copy with user-provided rating update.
  FormulaRecord withRatings({
    int? longevity,
    int? fresh,
    int? luxury,
    int? overall,
    String? notes,
  }) {
    return FormulaRecord(
      id: id,
      createdAt: createdAt,
      solventSystem: solventSystem,
      totalSolventFraction: totalSolventFraction,
      compoundNames: compoundNames,
      moleFractions: moleFractions,
      activityCoefficients: activityCoefficients,
      odorValues: odorValues,
      roiWeights: roiWeights,
      headspaceConc: headspaceConc,
      topAccords: topAccords,
      accordBalance: accordBalance,
      longevityRating: longevity ?? longevityRating,
      freshRating: fresh ?? freshRating,
      luxuryRating: luxury ?? luxuryRating,
      overallRating: overall ?? overallRating,
      notes: notes ?? this.notes,
      projectedLongevityHours: projectedLongevityHours,
      ovBalance: ovBalance,
    );
  }

  Map<String, dynamic> toJson() => {
    'id': id,
    'createdAt': createdAt.toIso8601String(),
    'solventSystem': solventSystem,
    'totalSolventFraction': totalSolventFraction,
    'compoundNames': compoundNames,
    'moleFractions': moleFractions,
    'activityCoefficients': activityCoefficients,
    'odorValues': odorValues,
    'roiWeights': roiWeights,
    'headspaceConc': headspaceConc,
    'topAccords': topAccords,
    'accordBalance': accordBalance,
    'longevityRating': longevityRating,
    'freshRating': freshRating,
    'luxuryRating': luxuryRating,
    'overallRating': overallRating,
    'notes': notes,
    'projectedLongevityHours': projectedLongevityHours,
    'ovBalance': ovBalance,
  };

  static FormulaRecord fromJson(Map<String, dynamic> json) => FormulaRecord(
    id: json['id'] as String,
    createdAt: DateTime.parse(json['createdAt'] as String),
    solventSystem: json['solventSystem'] as String,
    totalSolventFraction: (json['totalSolventFraction'] as num).toDouble(),
    compoundNames: List<String>.from(json['compoundNames'] as List),
    moleFractions: List<double>.from(
      (json['moleFractions'] as List).map((e) => (e as num).toDouble())),
    activityCoefficients: List<double>.from(
      (json['activityCoefficients'] as List).map((e) => (e as num).toDouble())),
    odorValues: List<double>.from(
      (json['odorValues'] as List).map((e) => (e as num).toDouble())),
    roiWeights: List<double>.from(
      (json['roiWeights'] as List).map((e) => (e as num).toDouble())),
    headspaceConc: List<double>.from(
      (json['headspaceConc'] as List).map((e) => (e as num).toDouble())),
    topAccords: List<String>.from(json['topAccords'] as List),
    accordBalance: Map<String, double>.from(
      (json['accordBalance'] as Map).map(
        (k, v) => MapEntry(k as String, (v as num).toDouble()))),
    longevityRating: json['longevityRating'] as int?,
    freshRating: json['freshRating'] as int?,
    luxuryRating: json['luxuryRating'] as int?,
    overallRating: json['overallRating'] as int?,
    notes: json['notes'] as String?,
    projectedLongevityHours:
        (json['projectedLongevityHours'] as num).toDouble(),
    ovBalance: (json['ovBalance'] as num).toDouble(),
  );
}

// ── In-Memory Formula History Store ───────────────────────
/// Thread-safe in-memory store.
/// Swap the backing store to SQLite / Hive / SharedPreferences
/// when persistent storage is needed.
class FormulaHistoryStore {
  FormulaHistoryStore._internal();
  static final FormulaHistoryStore instance = FormulaHistoryStore._internal();

  final List<FormulaRecord> _records = [];

  List<FormulaRecord> get all => List.unmodifiable(_records);

  void add(FormulaRecord record) {
    _records.insert(0, record); // newest first
    if (_records.length > 200) _records.removeLast(); // keep ≤200
  }

  FormulaRecord? findById(String id) {
    try {
      return _records.firstWhere((r) => r.id == id);
    } catch (_) {
      return null;
    }
  }

  void updateRatings(String id, {
    int? longevity,
    int? fresh,
    int? luxury,
    int? overall,
    String? notes,
  }) {
    final idx = _records.indexWhere((r) => r.id == id);
    if (idx < 0) return;
    _records[idx] = _records[idx].withRatings(
      longevity: longevity,
      fresh: fresh,
      luxury: luxury,
      overall: overall,
      notes: notes,
    );
  }

  void clear() => _records.clear();

  // ── Feature extraction for future ML ──────────────────
  /// Returns a flat feature vector per record, suitable for
  /// feeding into a classifier or regressor.
  /// Features are normalised where sensible.
  List<Map<String, double>> extractFeatures() {
    return _records.map((r) {
      final Map<String, double> features = {};
      features['n_compounds'] = r.compoundNames.length.toDouble();
      features['solvent_fraction'] = r.totalSolventFraction;
      features['ov_balance'] = r.ovBalance;
      features['projected_longevity'] = r.projectedLongevityHours;
      // Accord indicators (binary)
      for (final accord in [
        'Floral', 'Woody', 'Amber / Oriental', 'Musk',
        'Fresh / Citrus', 'Gourmand', 'Aromatic / Spicy', 'Fougère'
      ]) {
        features['accord_${accord.toLowerCase().replaceAll(' ', '_')}'] =
            r.topAccords.contains(accord) ? 1.0 : 0.0;
      }
      // ROI distribution stats
      if (r.roiWeights.isNotEmpty) {
        features['roi_max'] = r.roiWeights.reduce(
          (a, b) => a > b ? a : b);
        features['roi_min'] = r.roiWeights.reduce(
          (a, b) => a < b ? a : b);
        features['roi_std'] = _std(r.roiWeights);
      }
      return features;
    }).toList();
  }

  double _std(List<double> v) {
    if (v.isEmpty) return 0.0;
    final double mean = v.fold(0.0, (a, b) => a + b) / v.length;
    return (v.map((x) => (x - mean) * (x - mean))
        .fold(0.0, (a, b) => a + b) / v.length);
  }
}