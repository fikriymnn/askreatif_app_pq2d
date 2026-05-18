// ============================================================
// _CompoundSearchField — Native Flutter autocomplete
// Replaces drop_down_search_field to fix AXTree errors and
// broken selection on Windows desktop.
// ============================================================
// Drop-in replacement. Add this class to SolverPage.dart
// (or a separate file and import it).
// ============================================================

import 'package:flutter/material.dart';
import 'package:askreatif_app/Compound.dart';
import 'package:askreatif_app/colorTheme.dart';

class CompoundSearchField extends StatefulWidget {
  final int index;
  final String selectedValue;
  final TextEditingController controller;
  final ValueChanged<String> onSelected;

  const CompoundSearchField({
    Key? key,
    required this.index,
    required this.selectedValue,
    required this.controller,
    required this.onSelected,
  }) : super(key: key);

  @override
  State<CompoundSearchField> createState() => _CompoundSearchFieldState();
}

class _CompoundSearchFieldState extends State<CompoundSearchField> {
  final FocusNode _focusNode = FocusNode();
  OverlayEntry? _overlayEntry;
  final LayerLink _layerLink = LayerLink();
  List<String> _filtered = [];

  // All fragrance names (excluding solvents)
  static final List<String> _allFragrances =
      compounds
          .where((c) => !['Water', 'Ethanol', 'DPG'].contains(c.name))
          .map((c) => c.name)
          .toList();

  @override
  void initState() {
    super.initState();
    _filtered = _allFragrances;
    _focusNode.addListener(_onFocusChange);
  }

  @override
  void dispose() {
    _focusNode.removeListener(_onFocusChange);
    _focusNode.dispose();
    _removeOverlay();
    super.dispose();
  }

  void _onFocusChange() {
    if (!_focusNode.hasFocus) {
      // Delay so tap on item registers before overlay is removed
      Future.delayed(const Duration(milliseconds: 200), () {
        _removeOverlay();
      });
    }
  }

  void _removeOverlay() {
    _overlayEntry?.remove();
    _overlayEntry = null;
  }

  void _showOverlay() {
    _removeOverlay();
    if (_filtered.isEmpty) return;

    final RenderBox? renderBox = context.findRenderObject() as RenderBox?;
    if (renderBox == null) return;
    final double fieldWidth = renderBox.size.width;

    _overlayEntry = OverlayEntry(
      builder: (ctx) {
        return GestureDetector(
          // Transparent barrier — tap outside closes dropdown
          behavior: HitTestBehavior.translucent,
          onTap: () {
            _removeOverlay();
            _focusNode.unfocus();
          },
          child: Stack(
            children: [
              // Invisible full-screen barrier
              const SizedBox.expand(),
              // The actual dropdown
              CompositedTransformFollower(
                link: _layerLink,
                showWhenUnlinked: false,
                offset: Offset(0, renderBox.size.height + 4),
                child: Align(
                  alignment: Alignment.topLeft,
                  child: Material(
                    color: Colors.transparent,
                    child: Container(
                      width: fieldWidth,
                      constraints: const BoxConstraints(maxHeight: 240),
                      decoration: BoxDecoration(
                        color: AppTheme.bgCard,
                        borderRadius: BorderRadius.circular(4),
                        border: Border.all(color: AppTheme.border),
                        boxShadow: [
                          BoxShadow(
                            color: Colors.black.withOpacity(0.5),
                            blurRadius: 12,
                            offset: const Offset(0, 4),
                          ),
                        ],
                      ),
                      child: ClipRRect(
                        borderRadius: BorderRadius.circular(4),
                        child: ListView.builder(
                          padding: EdgeInsets.zero,
                          shrinkWrap: true,
                          itemCount: _filtered.length,
                          itemBuilder: (context, i) {
                            final String name = _filtered[i];
                            final bool isSelected =
                                name == widget.selectedValue;
                            return _DropdownItem(
                              name: name,
                              isSelected: isSelected,
                              onTap: () {
                                // Remove overlay FIRST, then call onSelected
                                _removeOverlay();
                                _focusNode.unfocus();
                                widget.controller.text = name;
                                widget.onSelected(name);
                              },
                            );
                          },
                        ),
                      ),
                    ),
                  ),
                ),
              ),
            ],
          ),
        );
      },
    );

    Overlay.of(context, rootOverlay: true).insert(_overlayEntry!);
  }

  void _onTextChanged(String query) {
    setState(() {
      _filtered =
          query.isEmpty
              ? _allFragrances
              : _allFragrances
                  .where((n) => n.toLowerCase().contains(query.toLowerCase()))
                  .toList();
    });
    // Rebuild overlay with new filtered list
    if (_overlayEntry != null) {
      _overlayEntry!.markNeedsBuild();
    } else if (_focusNode.hasFocus) {
      _showOverlay();
    }
  }

  @override
  Widget build(BuildContext context) {
    return CompositedTransformTarget(
      link: _layerLink,
      child: TextField(
        controller: widget.controller,
        focusNode: _focusNode,
        style: TextStyle(color: AppTheme.textPrimary, fontSize: 13),
        decoration: InputDecoration(
          labelText: 'Senyawa ${widget.index + 1}',
          labelStyle: TextStyle(color: AppTheme.textMuted, fontSize: 12),
          hintText: 'Cari senyawa...',
          hintStyle: TextStyle(color: AppTheme.textMuted, fontSize: 12),
          filled: true,
          fillColor: AppTheme.bgCard,
          contentPadding: const EdgeInsets.symmetric(
            horizontal: 12,
            vertical: 10,
          ),
          border: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.border),
          ),
          enabledBorder: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.border),
          ),
          focusedBorder: OutlineInputBorder(
            borderRadius: BorderRadius.circular(4),
            borderSide: BorderSide(color: AppTheme.green400),
          ),
          prefixIcon: Icon(Icons.search, color: AppTheme.textMuted, size: 16),
          suffixIcon:
              widget.selectedValue.isNotEmpty
                  ? Icon(Icons.check_circle, color: AppTheme.green400, size: 16)
                  : Icon(
                    Icons.arrow_drop_down,
                    color: AppTheme.textMuted,
                    size: 20,
                  ),
        ),
        onChanged: _onTextChanged,
        onTap: () {
          // Filter by current text and show dropdown
          _onTextChanged(widget.controller.text);
          _showOverlay();
        },
      ),
    );
  }
}

// ── Single dropdown item ───────────────────────────────────
class _DropdownItem extends StatefulWidget {
  final String name;
  final bool isSelected;
  final VoidCallback onTap;

  const _DropdownItem({
    required this.name,
    required this.isSelected,
    required this.onTap,
  });

  @override
  State<_DropdownItem> createState() => _DropdownItemState();
}

class _DropdownItemState extends State<_DropdownItem> {
  bool _hovered = false;

  @override
  Widget build(BuildContext context) {
    return MouseRegion(
      onEnter: (_) => setState(() => _hovered = true),
      onExit: (_) => setState(() => _hovered = false),
      child: GestureDetector(
        onTap: widget.onTap,
        child: Container(
          padding: const EdgeInsets.symmetric(horizontal: 12, vertical: 10),
          color:
              widget.isSelected
                  ? AppTheme.green900.withOpacity(0.6)
                  : _hovered
                  ? AppTheme.green900.withOpacity(0.3)
                  : Colors.transparent,
          child: Row(
            children: [
              Expanded(
                child: Text(
                  widget.name,
                  style: TextStyle(
                    color:
                        widget.isSelected
                            ? AppTheme.green200
                            : AppTheme.textPrimary,
                    fontSize: 13,
                  ),
                ),
              ),
              if (widget.isSelected)
                Icon(Icons.check, color: AppTheme.green400, size: 14),
            ],
          ),
        ),
      ),
    );
  }
}
