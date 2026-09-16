#!/usr/bin/env python3
"""Prune O2Physics DPL configuration JSONs down to the values that are not defaults.

Reads the defaults straight out of the task's .cxx (Configurable, ConfigurableAxis,
ConfigurableGroup and PROCESS_SWITCH), compares them against every .json in a folder,
and removes the keys that only restate a default. Keys with no counterpart in the
source are removed as well: DPL only ever queries the options a device registered, so
a stale key is inert, and a renamed one silently falls back to the compiled default.

Only the device keys owned by the analysed source are touched. Everything else in the
file (internal-dpl-*, sibling devices of the same workflow) is copied over verbatim.

Dry run unless --apply is given.
"""

import argparse
import datetime
import json
import math
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

# The DPL JSON writer serializes with max_digits10 of the declared type, so it never
# round-trips exactly: a double loses about one ulp (relative ~2e-16) and a float is cut
# to 9 significant digits (relative ~5e-9). The silent-equality band therefore follows
# the type. Between it and EPS_WARN the key is still pruned, but loudly -- see README.
EPS_EXACT_DOUBLE = 1e-12
EPS_EXACT_FLOAT = 1e-8
EPS_WARN = 1e-6

BACKUP_DIRNAME = "prePruningJsons"
PROBE_TAG = "JSONCLEANER"

VERDICT_KEEP = "KEEP"
VERDICT_PRUNE = "PRUNE"
VERDICT_NEAR = "NEAR"
VERDICT_ORPHAN = "ORPHAN"


def fatal(message):
    """The one way this tool is allowed to stop early."""
    sys.stderr.write("FATAL: %s\n" % message)
    sys.exit(1)


def warn(message):
    sys.stderr.write("WARNING: %s\n" % message)


# ---------------------------------------------------------------------------
# C++ scanning
# ---------------------------------------------------------------------------

def skip_literal(src, i):
    """Return the index just past the string or char literal starting at src[i]."""
    quote = src[i]
    i += 1
    while i < len(src):
        if src[i] == "\\":
            i += 2
            continue
        if src[i] == quote:
            return i + 1
        i += 1
    fatal("unterminated literal in source")


def strip_comments(src):
    """Drop // and /* */ comments, keeping string literals and line numbering intact.

    Has to be literal-aware: ccdbUrl's default is "http://alice-ccdb.cern.ch".
    """
    out = []
    i = 0
    while i < len(src):
        c = src[i]
        if c in "\"'":
            j = skip_literal(src, i)
            out.append(src[i:j])
            i = j
            continue
        if c == "/" and i + 1 < len(src):
            if src[i + 1] == "/":
                j = src.find("\n", i)
                i = len(src) if j < 0 else j
                continue
            if src[i + 1] == "*":
                j = src.find("*/", i + 2)
                j = len(src) if j < 0 else j + 2
                out.append("\n" * src.count("\n", i, j))
                i = j
                continue
        out.append(c)
        i += 1
    return "".join(out)


def match_bracket(src, i):
    """src[i] is an opening bracket; return the index just past its match."""
    pairs = {"{": "}", "(": ")", "[": "]", "<": ">"}
    opener = src[i]
    closer = pairs[opener]
    depth = 0
    while i < len(src):
        c = src[i]
        if c in "\"'":
            i = skip_literal(src, i)
            continue
        if c == opener:
            depth += 1
        elif c == closer:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    fatal("unbalanced '%s' in source" % opener)


def split_arguments(text):
    """Split an argument list on its top-level commas."""
    args = []
    depth = 0
    start = 0
    i = 0
    while i < len(text):
        c = text[i]
        if c in "\"'":
            i = skip_literal(text, i)
            continue
        if c in "{([":
            depth += 1
        elif c in "})]":
            depth -= 1
        elif c == "," and depth == 0:
            args.append(text[start:i].strip())
            start = i + 1
        i += 1
    args.append(text[start:].strip())
    return [a for a in args if a != ""]


def extract_enums(src):
    """Lift whole enum blocks verbatim, so the ROOT probe numbers them, not us."""
    blocks = []
    for m in re.finditer(r"^enum\s+(?:class\s+)?\w+[^{;]*\{", src, re.M):
        end = match_bracket(src, m.end() - 1)
        semi = src.find(";", end)
        blocks.append(src[m.start():(end if semi < 0 else semi + 1)])
    return blocks


def device_name_from_struct(struct_name):
    """O2's task-name rule: lowercase, with a dash before each inner capital."""
    out = []
    for i, ch in enumerate(struct_name):
        if ch.isupper() and i > 0:
            out.append("-")
        out.append(ch.lower())
    return "".join(out)


def find_task_name_overrides(src):
    """Pick up adaptAnalysisTask<X>(cfgc, TaskName{"y"}) so we do not guess wrongly."""
    overrides = {}
    for m in re.finditer(r"adaptAnalysisTask\s*<", src):
        angle_end = match_bracket(src, m.end() - 1)
        struct_name = src[m.end():angle_end - 1].strip().split("::")[-1]
        paren = src.find("(", angle_end)
        if paren < 0:
            continue
        args = src[paren:match_bracket(src, paren)]
        name = re.search(r'TaskName\s*\{\s*"([^"]*)"', args)
        if name:
            overrides[struct_name] = name.group(1)
    return overrides


def parse_template_argument(src, i):
    """At src[i] == '<', return (argument text, index just past '>')."""
    end = match_bracket(src, i)
    return src[i + 1:end - 1].strip(), end


def scan_configurables(body, offset=0):
    """Every Configurable / ConfigurableAxis declaration in a struct body."""
    found = []
    for m in re.finditer(r"\b(ConfigurableAxis|Configurable)\b", body):
        kind = m.group(1)
        i = m.end()
        while i < len(body) and body[i].isspace():
            i += 1
        template = None
        if kind == "Configurable":
            if i >= len(body) or body[i] != "<":
                continue  # not a declaration (e.g. the ConfigurableGroup base name)
            template, i = parse_template_argument(body, i)
        decl = re.match(r"\s*([A-Za-z_]\w*)\s*\{", body[i:])
        if not decl:
            continue
        brace = i + decl.end() - 1
        args = split_arguments(body[brace + 1:match_bracket(body, brace) - 1])
        if len(args) < 2 or not args[0].startswith('"'):
            continue
        found.append({
            "kind": kind,
            "template": template,
            "key": json.loads(args[0]),
            "default_expr": args[1],
            "position": offset + m.start(),
        })
    return found


def scan_process_switches(body, offset=0):
    """PROCESS_SWITCH(struct, name, help, default) and its _FULL variant."""
    found = []
    for m in re.finditer(r"\bPROCESS_SWITCH(_FULL)?\s*\(", body):
        paren = m.end() - 1
        args = split_arguments(body[paren + 1:match_bracket(body, paren) - 1])
        if m.group(1):  # _FULL(struct, method, name, help, default)
            if len(args) != 5:
                continue
            key, default_expr = args[2], args[4]
        else:
            if len(args) != 4:
                continue
            key, default_expr = args[1], args[3]
        found.append({
            "kind": "Configurable",
            "template": "bool",
            "key": key.strip(),
            "default_expr": default_expr,
            "position": offset + m.start(),
        })
    return found


def scan_source(path):
    """Return {device name: {json key path: declaration}} plus the source's enums."""
    try:
        with open(path, "r") as handle:
            raw = handle.read()
    except OSError as exc:
        fatal("cannot read source %s: %s" % (path, exc))

    src = strip_comments(raw)
    overrides = find_task_name_overrides(src)
    devices = {}

    for m in re.finditer(r"^struct\s+([A-Za-z_]\w*)\s*(?::[^{;]*)?\{", src, re.M):
        struct_name = m.group(1)
        brace = m.end() - 1
        body = src[brace + 1:match_bracket(src, brace) - 1]

        # ConfigurableGroups first: they claim the declarations inside their span.
        groups = []
        for g in re.finditer(r"struct\s*:\s*(?:public\s+)?ConfigurableGroup\s*\{", body):
            gbrace = g.end() - 1
            gend = match_bracket(body, gbrace)
            gbody = body[gbrace + 1:gend - 1]
            prefix = re.search(r'\bprefix\s*=\s*"([^"]*)"', gbody)
            if not prefix:
                fatal("ConfigurableGroup in %s has no prefix member" % struct_name)
            groups.append((gbrace, gend, prefix.group(1), gbody, gbrace + 1))

        declarations = {}
        claimed = []
        for gbrace, gend, prefix, gbody, goffset in groups:
            claimed.append((gbrace, gend))
            for decl in scan_configurables(gbody, goffset):
                decl["group"] = prefix
                declarations["%s.%s" % (prefix, decl["key"])] = decl

        for decl in scan_configurables(body) + scan_process_switches(body):
            if any(lo <= decl["position"] < hi for lo, hi in claimed):
                continue
            decl["group"] = None
            declarations[decl["key"]] = decl

        if declarations:
            name = overrides.get(struct_name, device_name_from_struct(struct_name))
            devices[name] = declarations

    if not devices:
        fatal("no task struct with configurables found in %s" % path)
    return devices, extract_enums(raw)


# ---------------------------------------------------------------------------
# Default values: literals in Python, anything else in ROOT
# ---------------------------------------------------------------------------

NUMBER_RE = re.compile(
    r"^[+-]?(?:0[xX][0-9a-fA-F]+|(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?)[fFlLuU]*$")


def float32(value):
    return struct.unpack("f", struct.pack("f", value))[0]


def parse_number_literal(text):
    """Return (value, is_float_literal) applying the literal's own type, or None."""
    text = text.strip()
    if not NUMBER_RE.match(text):
        return None
    suffix = re.search(r"[fFlLuU]*$", text).group(0)
    body = text[:len(text) - len(suffix)] if suffix else text
    if body.lower().startswith(("0x", "+0x", "-0x")):
        return int(body, 16), False
    if "f" in suffix or "F" in suffix:
        return float32(float(body)), True
    if any(c in body for c in ".eE"):
        return float(body), True
    return int(body), False


def literal_value(expr):
    """Evaluate a plain C++ literal, or return None if ROOT has to do it."""
    expr = expr.strip()
    if expr == "true":
        return True
    if expr == "false":
        return False
    if expr.startswith('"'):
        return json.loads(expr)
    number = parse_number_literal(expr)
    return None if number is None else number[0]


def axis_elements(expr):
    """The brace-enclosed element list of a ConfigurableAxis default."""
    expr = expr.strip()
    if not expr.startswith("{"):
        fatal("ConfigurableAxis default is not a brace list: %s" % expr)
    return split_arguments(expr[1:-1])


def collect_unresolved(devices):
    """Every default expression that is not a plain literal, deduplicated."""
    pending = set()
    for declarations in devices.values():
        for decl in declarations.values():
            exprs = (axis_elements(decl["default_expr"])
                     if decl["kind"] == "ConfigurableAxis" else [decl["default_expr"]])
            for expr in exprs:
                if expr == "VARIABLE_WIDTH":
                    continue
                if literal_value(expr) is None:
                    pending.add(expr)
    return sorted(pending)


def resolve_with_root(expressions, enum_blocks, cache_path):
    """Ask ROOT for the value of each non-literal default, caching the answers.

    Whole expressions go over, not just symbols: constants::math::TwoPI / 18 is a
    default too. The cast to double happens only at the print boundary, so the
    expression's own type (and any float32 narrowing it implies) is preserved.
    """
    cache = {}
    if os.path.exists(cache_path):
        try:
            with open(cache_path, "r") as handle:
                cache = json.load(handle)
        except (OSError, ValueError) as exc:
            warn("ignoring unreadable symbol cache %s: %s" % (cache_path, exc))

    missing = [e for e in expressions if e not in cache]
    if not missing:
        return cache

    if shutil.which("root") is None:
        fatal("these defaults are not literals and need ROOT to be evaluated:\n  %s\n"
              "Enter the O2Physics environment first (alienv enter O2Physics/latest), "
              "or point --symbol-cache at a cache that already has them."
              % "\n  ".join(missing))

    lines = ['#include <CommonConstants/MathConstants.h>',
             '#include <CommonConstants/PhysicsConstants.h>',
             '#include <cstdio>',
             'using namespace o2;',
             ""]
    lines.extend(enum_blocks)
    lines.append("")
    lines.append("void jsonCleanerProbe()")
    lines.append("{")
    for index, expr in enumerate(missing):
        lines.append('  std::printf("%s\\t%d\\t%%.17g\\n", static_cast<double>(%s));'
                     % (PROBE_TAG, index, expr))
    lines.append("}")

    workdir = tempfile.mkdtemp(prefix="jsonCleaner-")
    macro = os.path.join(workdir, "jsonCleanerProbe.C")
    try:
        with open(macro, "w") as handle:
            handle.write("\n".join(lines) + "\n")
        result = subprocess.run(["root", "-l", "-b", "-q", macro],
                                capture_output=True, text=True)
        harvested = {}
        for line in result.stdout.splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 3 and parts[0] == PROBE_TAG:
                harvested[missing[int(parts[1])]] = float(parts[2])
        unresolved = [e for e in missing if e not in harvested]
        if unresolved:
            fatal("ROOT could not evaluate:\n  %s\nROOT said:\n%s"
                  % ("\n  ".join(unresolved), result.stderr.strip() or "(nothing)"))
        cache.update(harvested)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    try:
        with open(cache_path, "w") as handle:
            json.dump(cache, handle, indent=2, sort_keys=True)
            handle.write("\n")
    except OSError as exc:
        warn("could not write symbol cache %s: %s" % (cache_path, exc))
    return cache


def type_category(template):
    """Map a Configurable's template argument onto how its value has to be compared."""
    name = re.sub(r"\b(o2::framework::|std::)", "", template or "").strip()
    if name == "bool":
        return "bool"
    if name == "string":
        return "string"
    if name == "float":
        return "float"
    if name in ("double", "long double"):
        return "double"
    if re.match(r"^(unsigned\s+)?(int|short|long|long\s+long|char)$", name) or \
       re.match(r"^u?int(8|16|32|64)_t$", name) or name in ("size_t", "unsigned"):
        return "int"
    return None


def evaluate_default(decl, symbols):
    """The compiled default, canonicalized the way C++ would store it.

    Two casts, because narrowing happens in both directions in real code:
    Configurable<double> radiusJet{..., 0.4f} widens a float literal, while
    Configurable<float> v0cospa{..., 0.995} narrows a double one.
    """
    if decl["kind"] == "ConfigurableAxis":
        values = []
        for expr in axis_elements(decl["default_expr"]):
            if expr == "VARIABLE_WIDTH":
                values.append(0.0)  # the sentinel is a plain 0 in the stored vector
                continue
            value = literal_value(expr)
            if value is None:
                value = symbols[expr]
            values.append(float(value))
        return "axis", values

    category = type_category(decl["template"])
    if category is None:
        fatal("unsupported Configurable type '%s' for key '%s'"
              % (decl["template"], decl["key"]))

    value = literal_value(decl["default_expr"])
    if value is None:
        value = symbols[decl["default_expr"]]  # already (double)(expr) from ROOT

    if category == "bool":
        return category, bool(value)
    if category == "string":
        if not isinstance(value, str):
            fatal("default of '%s' is not a string literal" % decl["key"])
        return category, value
    if category == "int":
        return category, int(value)
    if category == "float":
        return category, float32(float(value))
    return category, float(value)


# ---------------------------------------------------------------------------
# JSON side and comparison
# ---------------------------------------------------------------------------

def parse_json_value(raw, category, key):
    """DPL writes everything as strings; bring it back to the compared type."""
    try:
        if category == "axis":
            if not isinstance(raw, dict) or "values" not in raw:
                return None
            return [float(v) for v in raw["values"]]
        if isinstance(raw, dict) or isinstance(raw, list):
            return None
        if category == "bool":
            if isinstance(raw, bool):
                return raw
            text = str(raw).strip().lower()
            if text in ("true", "1"):
                return True
            if text in ("false", "0"):
                return False
            return None
        if category == "string":
            return str(raw)
        if category == "int":
            return int(str(raw).strip(), 0)
        return float(str(raw).strip())
    except (TypeError, ValueError):
        warn("could not read '%s' as %s from the JSON: %r" % (key, category, raw))
        return None


def compare_float(reference, candidate, category):
    """A float default is only written to 9 significant digits, a double to 17."""
    if reference == candidate:
        return "exact"
    if reference == 0.0 or not math.isfinite(reference):
        return "differ"  # a relative tolerance means nothing against zero
    relative = abs(candidate - reference) / abs(reference)
    if relative <= (EPS_EXACT_FLOAT if category == "float" else EPS_EXACT_DOUBLE):
        return "exact"
    if relative <= EPS_WARN:
        return "near"
    return "differ"


def compare_values(category, reference, candidate):
    """Return 'exact', 'near' or 'differ', plus the worst relative difference seen."""
    if category in ("bool", "string", "int"):
        return ("exact" if reference == candidate else "differ"), 0.0
    if category == "axis":
        if len(reference) != len(candidate):
            return "differ", 0.0
        outcome, worst = "exact", 0.0
        for ref, got in zip(reference, candidate):
            element = compare_float(ref, got, "double")  # axes are vectors of double
            if element != "exact" and ref != 0.0:
                worst = max(worst, abs(got - ref) / abs(ref))
            if element == "differ":
                outcome = "differ"
            elif element == "near" and outcome == "exact":
                outcome = "near"
        return outcome, worst
    outcome = compare_float(reference, candidate, category)
    worst = 0.0 if reference == 0.0 or outcome == "exact" \
        else abs(candidate - reference) / abs(reference)
    return outcome, worst


def describe(category, value):
    if category == "axis":
        head = ", ".join("%g" % v for v in value[:4])
        return "{%s%s} (%d values)" % (head, ", ..." if len(value) > 4 else "", len(value))
    if category == "bool":
        return "true" if value else "false"
    if category == "string":
        return '"%s"' % value
    if category in ("float", "double"):
        return "%.9g" % value
    return str(value)


# ---------------------------------------------------------------------------
# Per-file decisions
# ---------------------------------------------------------------------------

def flatten_device(block):
    """{key: value} and {key: (group, subkey)} for a device's JSON block."""
    flat, origin = {}, {}
    for key, value in block.items():
        if isinstance(value, dict) and "values" not in value:
            for subkey, subvalue in value.items():
                flat["%s.%s" % (key, subkey)] = subvalue
                origin["%s.%s" % (key, subkey)] = (key, subkey)
        else:
            flat[key] = value
            origin[key] = (None, key)
    return flat, origin


def decide_device(block, declarations, defaults):
    """Classify every key of one device block."""
    flat, origin = flatten_device(block)
    decisions = []
    for key in sorted(flat):
        if key not in declarations:
            decisions.append({"key": key, "verdict": VERDICT_ORPHAN, "detail":
                              "no such configurable in the source"})
            continue
        category, reference = defaults[key]
        candidate = parse_json_value(flat[key], category, key)
        if candidate is None:
            decisions.append({"key": key, "verdict": VERDICT_KEEP,
                              "detail": "unreadable value, kept untouched"})
            continue
        outcome, worst = compare_values(category, reference, candidate)
        if outcome == "exact":
            decisions.append({"key": key, "verdict": VERDICT_PRUNE,
                              "detail": "= default %s" % describe(category, reference)})
        elif outcome == "near":
            decisions.append({"key": key, "verdict": VERDICT_NEAR, "relative": worst,
                              "detail": "within %.2g of default %s"
                              % (worst, describe(category, reference))})
        else:
            decisions.append({"key": key, "verdict": VERDICT_KEEP,
                              "detail": "%s --> %s" % (describe(category, reference),
                                                       describe(category, candidate))})
    return decisions, flat, origin


def build_pruned_block(block, decisions, origin):
    """Rebuild a device block keeping only the KEEP verdicts, groups included."""
    keep = {d["key"] for d in decisions if d["verdict"] == VERDICT_KEEP}
    pruned = {}
    for key in keep:
        group, subkey = origin[key]
        if group is None:
            pruned[subkey] = block[subkey]
        else:
            pruned.setdefault(group, {})[subkey] = block[group][subkey]
    return pruned


def verify_invariant(name, flat, pruned, declarations, defaults, decisions):
    """merge(defaults, pruned) must reproduce merge(defaults, original), key by key.

    Runs in both modes, before anything is written. NEAR keys are allowed to move
    within the warning band -- that is exactly what pruning them means.
    """
    tolerance = {d["key"]: d["verdict"] for d in decisions}
    pruned_flat, _ = flatten_device(pruned)

    for key, decl in declarations.items():
        category, reference = defaults[key]
        before = parse_json_value(flat[key], category, key) if key in flat else reference
        if before is None:
            before = reference
        after = parse_json_value(pruned_flat[key], category, key) \
            if key in pruned_flat else reference
        outcome, worst = compare_values(category, before, after)
        allowed = "near" if tolerance.get(key) == VERDICT_NEAR else "exact"
        if outcome == "differ" or (allowed == "exact" and outcome != "exact"):
            fatal("invariant broken on %s, key '%s': effective value would change "
                  "from %s to %s (relative %.3g). Nothing was written."
                  % (name, key, describe(category, before), describe(category, after),
                     worst))

    for decision in decisions:
        if decision["verdict"] == VERDICT_ORPHAN and decision["key"] in declarations:
            fatal("invariant broken on %s: key '%s' was called an orphan but exists "
                  "in the source. Nothing was written." % (name, decision["key"]))


# ---------------------------------------------------------------------------
# Driving
# ---------------------------------------------------------------------------

def report_device(file_name, name, decisions):
    counts = {v: 0 for v in (VERDICT_KEEP, VERDICT_PRUNE, VERDICT_NEAR, VERDICT_ORPHAN)}
    print("  device: %s" % name)
    for decision in sorted(decisions, key=lambda d: (d["verdict"], d["key"])):
        counts[decision["verdict"]] += 1
        print("    %-7s %-52s %s"
              % (decision["verdict"], decision["key"], decision["detail"]))
        if decision["verdict"] == VERDICT_NEAR:
            warn("%s: '%s' is %s -- pruned as a default. If that difference was "
                 "deliberate, restore it from the backup."
                 % (file_name, decision["key"], decision["detail"]))
    print("    -- kept %d, pruned %d, pruned-with-warning %d, orphans removed %d"
          % (counts[VERDICT_KEEP], counts[VERDICT_PRUNE],
             counts[VERDICT_NEAR], counts[VERDICT_ORPHAN]))
    return counts


def process_file(path, devices, defaults, apply_changes, backup_dir):
    """Classify one config, verify it, and rewrite it if --apply was given."""
    name = os.path.basename(path)
    try:
        with open(path, "r") as handle:
            config = json.load(handle)
    except (OSError, ValueError) as exc:
        fatal("cannot read config %s: %s" % (path, exc))

    present = [d for d in devices if d in config and isinstance(config[d], dict)]
    print("\n=== %s ===" % name)
    if not present:
        warn("%s has no device key from this source (%s); skipped."
             % (name, ", ".join(sorted(devices))))
        return None

    pruned_config = dict(config)
    totals = {v: 0 for v in (VERDICT_KEEP, VERDICT_PRUNE, VERDICT_NEAR, VERDICT_ORPHAN)}
    changed = False

    for device in present:
        decisions, flat, origin = decide_device(
            config[device], devices[device], defaults[device])
        pruned_block = build_pruned_block(config[device], decisions, origin)
        verify_invariant(name, flat, pruned_block,
                         devices[device], defaults[device], decisions)
        counts = report_device(name, device, decisions)
        for verdict, count in counts.items():
            totals[verdict] += count
        if pruned_block != config[device]:
            changed = True
        pruned_config[device] = pruned_block

    untouched = [k for k in config if k not in present]
    if untouched:
        print("  untouched: %s" % ", ".join(untouched))

    if not apply_changes:
        print("  (dry run -- nothing written)")
        return totals
    if not changed:
        print("  already minimal -- left alone")
        return totals

    os.makedirs(backup_dir, exist_ok=True)
    shutil.copy2(path, os.path.join(backup_dir, name))
    with open(path, "w") as handle:
        json.dump(pruned_config, handle, indent=4)
        handle.write("\n")
    print("  written (backup in %s)" % backup_dir)
    return totals


def serialize_default(category, value):
    """Write a default back in the DPL string style, for the reference dumps."""
    if category == "axis":
        return {"values": ["%.17g" % v for v in value]}
    if category == "bool":
        return "true" if value else "false"
    if category == "string":
        return value
    if category == "int":
        return str(value)
    return "%.17g" % value


def dump_defaults(devices, defaults, directory, label):
    os.makedirs(directory, exist_ok=True)
    written = []
    for device in sorted(devices):
        block = {}
        for key in sorted(devices[device]):
            category, value = defaults[device][key]
            group = devices[device][key]["group"]
            plain = key.split(".", 1)[1] if group else key
            if group:
                block.setdefault(group, {})[plain] = serialize_default(category, value)
            else:
                block[plain] = serialize_default(category, value)
        stem = label if label else device
        path = os.path.join(directory, "dpl-config-%s-Defaults.json" % stem)
        with open(path, "w") as handle:
            json.dump({device: block}, handle, indent=4)
            handle.write("\n")
        written.append(path)
    return written


def main():
    global EPS_WARN

    parser = argparse.ArgumentParser(
        description="Prune O2Physics DPL configuration JSONs down to non-defaults.")
    parser.add_argument("--source", required=True,
                        help="the task .cxx whose defaults are the reference")
    parser.add_argument("--configs",
                        help="a .json file or a folder of them to prune")
    parser.add_argument("--apply", action="store_true",
                        help="actually rewrite the files (default is a dry run)")
    parser.add_argument("--dump-defaults", metavar="DIR",
                        help="also write the full compiled defaults as a reference JSON")
    parser.add_argument("--defaults-label", metavar="LABEL",
                        help="name the dump dpl-config-LABEL-Defaults.json")
    parser.add_argument("--symbol-cache", metavar="PATH",
                        default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "symbolCache.json"),
                        help="where ROOT-resolved defaults are remembered")
    parser.add_argument("--eps-warn", type=float, default=EPS_WARN,
                        help="relative band below which a difference is pruned but "
                             "warned about (default %g)" % EPS_WARN)
    args = parser.parse_args()
    EPS_WARN = args.eps_warn

    if not args.configs and not args.dump_defaults:
        fatal("nothing to do: give --configs, --dump-defaults, or both")

    devices, enum_blocks = scan_source(args.source)
    symbols = resolve_with_root(collect_unresolved(devices), enum_blocks,
                                args.symbol_cache)

    defaults = {}
    for device, declarations in devices.items():
        defaults[device] = {k: evaluate_default(d, symbols)
                            for k, d in declarations.items()}

    print("Source: %s" % args.source)
    for device in sorted(devices):
        print("  %s: %d configurables" % (device, len(devices[device])))

    if args.dump_defaults:
        for path in dump_defaults(devices, defaults, args.dump_defaults,
                                  args.defaults_label):
            print("  defaults written to %s" % path)

    if not args.configs:
        return

    if os.path.isdir(args.configs):
        targets = sorted(os.path.join(args.configs, f)
                         for f in os.listdir(args.configs) if f.endswith(".json"))
        root = args.configs
    else:
        targets = [args.configs]
        root = os.path.dirname(os.path.abspath(args.configs))
    if not targets:
        fatal("no .json files found in %s" % args.configs)

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    backup_dir = os.path.join(root, BACKUP_DIRNAME, stamp)

    grand = {v: 0 for v in (VERDICT_KEEP, VERDICT_PRUNE, VERDICT_NEAR, VERDICT_ORPHAN)}
    for target in targets:
        if os.path.join(BACKUP_DIRNAME, "") in target:
            continue  # never prune our own backups
        totals = process_file(target, devices, defaults, args.apply, backup_dir)
        if totals:
            for verdict, count in totals.items():
                grand[verdict] += count

    print("\nTotal: kept %d, pruned %d, pruned-with-warning %d, orphans removed %d"
          % (grand[VERDICT_KEEP], grand[VERDICT_PRUNE],
             grand[VERDICT_NEAR], grand[VERDICT_ORPHAN]))
    if not args.apply:
        print("Dry run. Re-run with --apply to write the pruned files.")


if __name__ == "__main__":
    main()