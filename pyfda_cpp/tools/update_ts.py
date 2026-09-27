#!/usr/bin/env python3
"""Update the translation files translations/pyfda_cpp_<lang>.ts from the GUI sources.

Collects the string literals of all tr("...") and QT_TR_NOOP("...") calls in src/gui,
keeps the existing translations, adds new strings as unfinished and drops strings that
are no longer used. The .ts files have Qt Linguist's format, but the application reads
them directly (src/gui/translator.cpp), so neither lupdate nor lrelease are needed.

  python3 tools/update_ts.py           update all .ts files
  python3 tools/update_ts.py --check   only check: exit 1 if a string is new, obsolete,
                                       untranslated or its %1 ... placeholders differ
"""
import glob
import os
import re
import sys
import xml.etree.ElementTree as ET

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT = re.compile(r'\s*"((?:[^"\\]|\\.)*)"')
CALL = re.compile(r'\b(?:tr|QT_TR_NOOP)\(')
ESC = {'n': 10, 't': 9, '"': 34, '\\': 92, "'": 39}


def unescape(s):
    """C string literal (UTF-8 source) -> str"""
    raw, out, i = s.encode('utf-8'), bytearray(), 0
    while i < len(raw):
        if raw[i] == 0x5c:
            n = chr(raw[i + 1])
            if n == 'x':
                j, h = i + 2, ''
                while j < len(raw) and len(h) < 2 and chr(raw[j]) in '0123456789abcdefABCDEF':
                    h += chr(raw[j])
                    j += 1
                out.append(int(h, 16))
                i = j
                continue
            out.append(ESC[n])
            i += 2
            continue
        out.append(raw[i])
        i += 1
    return out.decode('utf-8')


def extract():
    """source strings in the order of their first appearance"""
    strings = {}
    for f in sorted(glob.glob(os.path.join(ROOT, 'src/gui/*.[ch]pp'))):
        text = open(f, encoding='utf-8').read()
        for m in CALL.finditer(text):
            p, parts = m.end(), []
            while (lm := LIT.match(text, p)):
                parts.append(lm.group(1))
                p = lm.end()
            if parts:
                strings.setdefault(unescape(''.join(parts)), os.path.basename(f))
    return strings


def read_ts(path):
    tr = {}
    if os.path.exists(path):
        for msg in ET.parse(path).getroot().iter('message'):
            t = msg.find('translation')
            if t is not None and t.get('type') not in ('unfinished', 'obsolete') and t.text:
                tr[msg.findtext('source')] = t.text
    return tr


def placeholders(s):
    return sorted(re.findall(r'%\d', s))


def write_ts(path, lang, strings, tr):
    root = ET.Element('TS', version='2.1', language=lang)
    ctx = ET.SubElement(root, 'context')
    ET.SubElement(ctx, 'name').text = 'pyfda'
    for src, fname in strings.items():
        msg = ET.SubElement(ctx, 'message')
        ET.SubElement(msg, 'location', filename='../src/gui/' + fname)
        ET.SubElement(msg, 'source').text = src
        t = ET.SubElement(msg, 'translation')
        if src in tr:
            t.text = tr[src]
        else:
            t.set('type', 'unfinished')
    ET.indent(root)
    with open(path, 'w', encoding='utf-8', newline='\n') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n<!DOCTYPE TS>\n')
        f.write(ET.tostring(root, encoding='unicode') + '\n')


def main():
    check = '--check' in sys.argv
    strings = extract()
    ok = True
    for path in sorted(glob.glob(os.path.join(ROOT, 'translations/pyfda_cpp_*.ts'))):
        lang = re.search(r'pyfda_cpp_(\w+)\.ts$', path).group(1)
        tr = read_ts(path)
        missing = [s for s in strings if s not in tr]
        obsolete = [s for s in tr if s not in strings]
        bad = [s for s in strings if s in tr and placeholders(s) != placeholders(tr[s])]
        name = os.path.relpath(path, ROOT)
        print(f'{name}: {len(strings)} strings, {len(missing)} untranslated, '
              f'{len(obsolete)} obsolete, {len(bad)} with wrong placeholders')
        for s in missing:
            print('  untranslated:', repr(s))
        for s in bad:
            print('  placeholders:', repr(s), '->', repr(tr[s]))
        if check:
            ok = ok and not (missing or obsolete or bad)
        else:
            write_ts(path, lang, strings, {s: t for s, t in tr.items() if s in strings})
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
