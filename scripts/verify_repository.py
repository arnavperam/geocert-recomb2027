#!/usr/bin/env python3
"""Verify the published checkout, paper identities and local documentation links."""
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]

def main():
    provenance = json.loads((ROOT / 'PROVENANCE.json').read_text())
    project = provenance['repository'].removeprefix('https://github.com/')
    owner = project.split('/')[0]
    project_pattern = re.compile(re.escape(owner) + r'/[A-Za-z0-9_.-]+')
    files = []
    for line in (ROOT / 'MANIFEST.sha256').read_text().splitlines():
        expected, relative = line.split('  ', 1)
        path = ROOT / relative
        if not path.resolve().is_relative_to(ROOT):
            raise ValueError(f'unsafe manifest path: {relative}')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f'checksum mismatch: {relative}')
        files.append(path)
        if path.suffix in {'.md', '.tex', '.json', '.py', '.yml', '.yaml', '.txt', '.cff'}:
            for identifier in project_pattern.findall(unquote(path.read_text())):
                if identifier.rstrip('.') != project:
                    raise ValueError(f'Unexpected project reference in {relative}: {identifier}')
    for path in files:
        if path.suffix == '.md':
            for raw in re.findall(r'\[[^\]]*\]\(([^)]+)\)', path.read_text()):
                url = raw.split(' "', 1)[0].strip('<>')
                if re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', url) or url.startswith('#'):
                    continue
                target = unquote(url.split('#', 1)[0])
                if not (path.parent / target).exists():
                    raise ValueError(f'broken local link in {path.relative_to(ROOT)}: {url}')
        if path.suffix == '.tex':
            for target in re.findall(r'\\(?:input|includegraphics)(?:\[[^\]]*\])?\{([^}]+)\}', path.read_text()):
                if not (path.parent / target).exists():
                    raise ValueError(f'missing LaTeX dependency: {target}')
    atlas = ROOT / 'recomb2027-atlas'
    record = json.loads((atlas / 'release/PUBLICATION.json').read_text())
    archive = json.loads((atlas / 'artifact_parts/manifest.json').read_text())
    assert archive['distribution_version'] == provenance['distribution_version'] == (atlas / 'VERSION').read_text().strip()
    assert archive['sha256'] == record['base_archive_sha256'] == provenance['base_archive_sha256']
    assert record['scientific_source_url'] == provenance['repository']
    assert record['document_sha256'] == provenance['document_sha256']
    for name, expected in record['document_sha256'].items():
        assert hashlib.sha256((atlas / 'manuscript' / name).read_bytes()).hexdigest() == expected
    for name, expected in record['extension_records_sha256'].items():
        assert hashlib.sha256((atlas / name).read_bytes()).hexdigest() == expected
    print(json.dumps({'status': 'passed', 'manifest_entries': len(files), 'paper_hashes': 'matched', 'local_links': 'resolved', 'latex_dependencies': 'present'}))

if __name__ == '__main__':
    main()
