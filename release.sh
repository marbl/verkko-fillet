#!/usr/bin/env bash

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: ./release.sh <version>"
    echo "Example: ./release.sh 0.2.0"
    exit 1
fi

VERSION="$1"
TAG="v${VERSION}"
TODAY=$(date +"%Y-%m-%d")

# Validate version
if ! echo "$VERSION" | grep -Eq '^[0-9]+\.[0-9]+\.[0-9]+$'; then
    echo "Error: version must look like 0.2.0"
    exit 1
fi

# Required files
for FILE in version.txt CHANGELOG.md CITATION.cff pyproject.toml; do
    if [ ! -f "$FILE" ]; then
        echo "Error: required file not found: $FILE"
        exit 1
    fi
done

# Require clean working tree
if [ -n "$(git status --porcelain)" ]; then
    echo "Error: working tree is not clean."
    echo "Commit or stash your changes first."
    exit 1
fi

# Make sure tag does not already exist
if git rev-parse "$TAG" >/dev/null 2>&1; then
    echo "Error: tag $TAG already exists."
    exit 1
fi

CURRENT_VERSION=$(cat version.txt)

echo "Preparing verkko-fillet release"
echo "  Current version: $CURRENT_VERSION"
echo "  New version:     $VERSION"
echo

# Find previous tag
PREVIOUS_TAG=$(git describe --tags --abbrev=0 2>/dev/null || true)

if [ -n "$PREVIOUS_TAG" ]; then
    RANGE="${PREVIOUS_TAG}..HEAD"
    echo "Generating changelog from $PREVIOUS_TAG to HEAD"
else
    RANGE="HEAD"
    echo "No previous tag found. Using full git history."
fi

# Collect commits
FEATURES=$(git log "$RANGE" --pretty=format:"%s" | \
    grep '^feat:' | sed 's/^feat:[[:space:]]*/- /' || true)

FIXES=$(git log "$RANGE" --pretty=format:"%s" | \
    grep '^fix:' | sed 's/^fix:[[:space:]]*/- /' || true)

DOCS=$(git log "$RANGE" --pretty=format:"%s" | \
    grep '^docs:' | sed 's/^docs:[[:space:]]*/- /' || true)

CHANGES=$(git log "$RANGE" --pretty=format:"%s" | \
    grep '^refactor:' | sed 's/^refactor:[[:space:]]*/- /' || true)

PERFORMANCE=$(git log "$RANGE" --pretty=format:"%s" | \
    grep '^perf:' | sed 's/^perf:[[:space:]]*/- /' || true)

# Build changelog section
CHANGELOG_SECTION="## [$VERSION] - $TODAY
"

if [ -n "$FEATURES" ]; then
    CHANGELOG_SECTION+="
### Added
$FEATURES
"
fi

if [ -n "$FIXES" ]; then
    CHANGELOG_SECTION+="
### Fixed
$FIXES
"
fi

if [ -n "$CHANGES" ]; then
    CHANGELOG_SECTION+="
### Changed
$CHANGES
"
fi

if [ -n "$PERFORMANCE" ]; then
    CHANGELOG_SECTION+="
### Performance
$PERFORMANCE
"
fi

if [ -n "$DOCS" ]; then
    CHANGELOG_SECTION+="
### Documentation
$DOCS
"
fi

# If nothing matched
if [ -z "$FEATURES$FIXES$DOCS$CHANGES$PERFORMANCE" ]; then
    echo "Error: no conventional commits found since previous release."
    echo
    echo "Expected commit messages such as:"
    echo "  feat: add gap-filling for chr21"
    echo "  fix: correct chromosome assignment"
    echo "  docs: update README"
    exit 1
fi

# Update version.txt
printf '%s\n' "$VERSION" > version.txt

# Update pyproject.toml
python3 - "$VERSION" <<'PY'
import re
import sys
from pathlib import Path

version = sys.argv[1]
path = Path("pyproject.toml")
text = path.read_text()

new_text, count = re.subn(
    r'(?m)^version = "[^"]*"',
    f'version = "{version}"',
    text,
    count=1,
)

if count == 0:
    raise SystemExit("Error: version field not found in pyproject.toml")

path.write_text(new_text)
PY

# Update CITATION.cff
python3 - "$VERSION" <<'PY'
import sys
from pathlib import Path

version = sys.argv[1]
path = Path("CITATION.cff")

lines = path.read_text().splitlines()

found = False
output = []

for line in lines:
    if line.startswith("version:"):
        output.append(f'version: "{version}"')
        found = True
    else:
        output.append(line)

if not found:
    raise SystemExit("Error: version field not found in CITATION.cff")

path.write_text("\n".join(output) + "\n")
PY

# Insert new changelog after title/introduction
python3 - "$VERSION" "$TODAY" "$CHANGELOG_SECTION" <<'PY'
import sys
from pathlib import Path

version = sys.argv[1]
section = sys.argv[3]

path = Path("CHANGELOG.md")
text = path.read_text()

if f"## [{version}]" in text:
    raise SystemExit(f"Error: version {version} already exists in CHANGELOG.md")

lines = text.splitlines()

insert_at = None

for i, line in enumerate(lines):
    if line.startswith("## "):
        insert_at = i
        break

if insert_at is None:
    new_text = text.rstrip() + "\n\n" + section.rstrip() + "\n"
else:
    new_lines = (
        lines[:insert_at]
        + section.rstrip().splitlines()
        + [""]
        + lines[insert_at:]
    )
    new_text = "\n".join(new_lines) + "\n"

path.write_text(new_text)
PY

echo
echo "Generated CHANGELOG entry:"
echo
echo "$CHANGELOG_SECTION"

git add version.txt CHANGELOG.md CITATION.cff pyproject.toml

git commit -m "chore: prepare release $TAG"

git tag -a "$TAG" -m "verkko-fillet $TAG"

echo
echo "Release prepared successfully."
echo
echo "Created:"
echo "  commit: chore: prepare release $TAG"
echo "  tag:    $TAG"
echo
echo "Review:"
echo "  git show $TAG"
echo
echo "Then push:"
echo "  git push origin main --tags"
